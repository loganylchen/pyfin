"""Pipeline orchestrator: interval -> candidates -> scoring -> EM -> quantification."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pysam

from fin.analysis.assignments import em_with_coherence
from fin.analysis.quantification import (
    QuantResult,
    _exons_from_candidate,
    aggregate_across_intervals,
    compute_fulllen_frac,
    fulllen_fraction_drops,
    isoform_fraction_drops,
    mono_exon_drops,
    soft_mass_ratio_drops,
    polya5p_drops,
    quantify_transcripts,
)
from fin.candidates.canonical import chain_all_canonical, parse_motifs
from fin.candidates.dataclasses import CandidateSet
from fin.candidates.discovery import discover_candidates, merge_fusion_candidates
from fin.io.interval_manager import GenomicInterval, generate_isolated_intervals
from fin.pipeline.config import PipelineConfig
from fin.scoring.em_inputs import build_em_matrices
from fin.scoring.krill_tiebreak import krill_tiebreak
from fin.scoring.mappy_distance import compute_mappy_distance

logger = logging.getLogger(__name__)

# Pure tie-break junction-NLL M2-EM constants (SIRV-validated experiment values).
# K (transcript-bp each side of the wobbling junction) comes from the configurable
# m2_tiebreak_junction_k; FLANK is the diff-window padding and PAD the penalty for
# an unscorable tie cell (a tied candidate the eventalign window could not score).
_M2_EM_FLANK = 2
_M2_EM_PAD = 1.0


class PipelineRunner:
    """Orchestrates the full pyfin pipeline."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self._gtf_reader = None
        self._genome_fasta = None
        self._signal_reader = None
        # Genome-wide mappy aligner for fusion soft-clip arm re-alignment; built
        # lazily on first fusion interval and cached (indexing is expensive).
        self._fusion_genome_aligner = None
        self._fusion_aligner_built = False

    def setup(self):
        """Open file handles and load references.

        All signal scoring is in-memory krill (no f5c CLI / external-tool
        validation / f5c index build required).
        """
        # Load GTF
        if self.config.gtf_path:
            from fin.io.io_gtf import GTFReader

            self._gtf_reader = GTFReader(self.config.gtf_path)
            self._gtf_reader.open()
            self._gtf_reader.parse()

        # Load genome FASTA (full sequence per chrom)
        if self.config.genome_fasta_path:
            self._genome_fasta = self._load_genome_fasta(self.config.genome_fasta_path)

        # Open signal reader
        if self.config.signal_path:
            self._signal_reader = self._open_signal_reader()

        # Create work directory
        Path(self.config.work_dir).mkdir(parents=True, exist_ok=True)

        logger.info("Pipeline setup complete")

    def run(self) -> Dict[str, QuantResult]:
        """Run the full pipeline across all intervals.

        Returns:
            Aggregated quantification results.
        """
        # Generate intervals
        result = generate_isolated_intervals(
            self.config.bam_path,
            gtf_path=self.config.gtf_path,
            max_gap=self.config.max_gap,
            max_reads=self.config.max_reads,
        )
        intervals = result["intervals"]
        logger.info("Generated %d intervals", len(intervals))

        # Process each interval (serial by default; opt-in process parallelism).
        all_quant_results: List[List[QuantResult]] = []
        if self.config.threads <= 1:
            for i, interval in enumerate(intervals):
                logger.info("Processing interval %d/%d: %s", i + 1, len(intervals), interval.region_string)
                quant = self.process_interval(interval)
                if quant:
                    all_quant_results.append(quant)
        else:
            from fin.pipeline.parallel import run_parallel

            log_level = logging.getLevelName(logging.getLogger("fin").getEffectiveLevel())
            all_quant_results = run_parallel(
                self.config,
                intervals,
                self.config.threads,
                self.config.gpu_workers,
                log_level,
            )

        # Aggregate across intervals
        aggregated = aggregate_across_intervals(all_quant_results)
        return self._finalize_and_write(
            aggregated,
            output_gtf=self.config.output_gtf,
            output_tsv=self.config.output_tsv,
        )

    def _finalize_and_write(
        self,
        aggregated: Dict[str, QuantResult],
        output_gtf: Optional[str],
        output_tsv: Optional[str],
    ) -> Dict[str, QuantResult]:
        """Apply post-EM filters, resolve gene_ids, optionally write GTF/TSV."""
        # R5 ablation: enable_score_filter=False disables all post-EM filters so
        # the ablation row sees unfiltered EM output (AC7).
        _score_filter_on = getattr(self.config, "enable_score_filter", True)

        # Diagnostic: write unfiltered scores TSV BEFORE the post-EM filter
        # cascade so downstream FN-root-cause analysis can distinguish
        # "candidate never reached EM" from "candidate dropped by a filter".
        if getattr(self.config, "write_unfiltered_scores", False) and output_tsv:
            from fin.io.io_tsv import write_scoring_tsv

            unfiltered_path = str(Path(output_tsv).with_suffix(".unfiltered.tsv"))
            transcript_lengths_uf = {
                cid: sum(end - start for start, end in qr.exons) if qr.exons else 0
                for cid, qr in aggregated.items()
            }
            # Resolve gene_ids on the unfiltered snapshot too so downstream
            # consumers don't see empty gene_id columns.
            for cid, qr in aggregated.items():
                if qr.source == "gtf" and self._gtf_reader:
                    tx = self._gtf_reader.get_transcript(cid)
                    if tx and not qr.gene_id:
                        qr.gene_id = tx.gene_id
                if not qr.gene_id:
                    qr.gene_id = qr.candidate_id
            write_scoring_tsv(aggregated, transcript_lengths_uf, unfiltered_path)
            logger.info(
                "Wrote unfiltered scoring TSV (pre-filter): %s (n=%d)",
                unfiltered_path,
                len(aggregated),
            )

        # A4: post-EM abundance floor with per-source thresholds (on soft EM
        # abundance). NOVEL transcripts face min_abundance (the `fin` CLI defaults
        # it to 3). GTF transcripts face their own lighter floor min_gtf_abundance
        # (CLI default 1) so a GTF candidate whose EM soft-mass is below 1 read is
        # dropped — pyfin is not a copy-annotation tool — while genuine annotated
        # transcripts with real EM support are kept. When floor_gtf_abundance is
        # set (--floor-gtf-abundance) GTF is raised to the NOVEL floor — precisely
        # max(min_gtf_abundance, min_abundance) so the toggle never LOWERS an
        # explicit GTF floor — reproducing the SIRV gffcompare T>=3 with-GTF
        # operating point. Fusion is ALWAYS exempt so opt-in fusion calls are
        # never silently dropped. SIRV-tuned: these floors are overfit and gut
        # genuine low-abundance isoform recall on real data; disable via
        # --min-abundance 0 / --min-gtf-abundance 0.
        # NOTE: the separate ablation harness (fin/ablation/runner.py) keeps GTF
        # exempt regardless — that divergence is intentional (no-GTF de novo).
        floor_gtf = self.config.floor_gtf_abundance
        gtf_floor = (
            max(self.config.min_gtf_abundance, self.config.min_abundance)
            if floor_gtf
            else self.config.min_gtf_abundance
        )
        if _score_filter_on and (self.config.min_abundance > 0.0 or gtf_floor > 0.0):
            before = len(aggregated)
            aggregated = {
                cid: qr
                for cid, qr in aggregated.items()
                if qr.source == "fusion"
                or (
                    qr.abundance >= gtf_floor
                    if qr.source == "gtf"
                    else qr.abundance >= self.config.min_abundance
                )
            }
            dropped = before - len(aggregated)
            if dropped:
                logger.info(
                    "Dropped %d transcripts with abundance < floor "
                    "(novel < %.3f, gtf < %.3f)",
                    dropped,
                    self.config.min_abundance,
                    gtf_floor,
                )

        # Minimum isoform fraction (locus-relative abundance) filter. Drops
        # NOVEL multi-exon transcripts whose abundance is below
        # `min_isoform_fraction` of the dominant overlapping novel isoform at
        # the same locus — the standard Cufflinks --min-isoform-fraction /
        # StringTie -f heuristic (the low-fraction tail is enriched for
        # incompletely-spliced precursors and assembly artifacts).
        # gtf/fusion/mono candidates are EXEMPT. SIRV WARNING: the conservative
        # literature default (0.01) is used; SIRV's F1-optimal ~0.4 is overfit.
        if _score_filter_on and self.config.min_isoform_fraction > 0.0:
            drop_ids = isoform_fraction_drops(
                aggregated, self.config.min_isoform_fraction
            )
            if drop_ids:
                aggregated = {
                    cid: qr
                    for cid, qr in aggregated.items()
                    if cid not in drop_ids
                }
                logger.info(
                    "Dropped %d novel transcripts with isoform fraction < %.3f",
                    len(drop_ids),
                    self.config.min_isoform_fraction,
                )

        # Soft-mass / hard-read ratio filter. Drops a NOVEL multi-exon candidate
        # whose EM soft abundance is inflated far above its honest hard-read
        # (argmax) count — the signature of a wobble shadow that steals little
        # real read support but accumulates fractional soft crumbs from a
        # high-abundance structural near-copy (anchor). Catches the HIGH-relative
        # -abundance shadows that the isoform_fraction / cluster-recheck levers
        # cannot (a shadow at 32% of its anchor passes a fraction gate). Pure EM
        # evidence; gtf/fusion/mono exempt. max_soft_mass_ratio<=0 disables.
        if _score_filter_on and self.config.max_soft_mass_ratio > 0.0:
            drop_ids = soft_mass_ratio_drops(
                aggregated, self.config.max_soft_mass_ratio
            )
            if drop_ids:
                aggregated = {
                    cid: qr
                    for cid, qr in aggregated.items()
                    if cid not in drop_ids
                }
                logger.info(
                    "Dropped %d novel transcripts with soft/hard read ratio >= %.2f",
                    len(drop_ids),
                    self.config.max_soft_mass_ratio,
                )

        # Lever 3 — mono-exon (single-exon) read-support gate. Drop a NOVEL
        # single-exon candidate whose hard read count < min_mono_exon_reads OR
        # genomic length < min_mono_exon_length. Suppresses single-exon de novo
        # noise (IsoQuant drops novel unspliced by default for ONT) WITHOUT a
        # blanket drop — a high-support/long real intronless gene survives.
        # gtf/fusion/multi-exon exempt. Needs the master switch AND >=1 threshold
        # > 0; otherwise no drops (byte-identical).
        if (
            _score_filter_on
            and getattr(self.config, "drop_mono_exon_novel", False)
            and (
                self.config.min_mono_exon_reads > 0
                or self.config.min_mono_exon_length > 0
            )
        ):
            drop_ids = mono_exon_drops(
                aggregated,
                min_reads=self.config.min_mono_exon_reads,
                min_len=self.config.min_mono_exon_length,
            )
            if drop_ids:
                aggregated = {
                    cid: qr
                    for cid, qr in aggregated.items()
                    if cid not in drop_ids
                }
                logger.info(
                    "Dropped %d novel single-exon transcripts (mono gate: "
                    "min_reads=%d min_len=%d)",
                    len(drop_ids),
                    self.config.min_mono_exon_reads,
                    self.config.min_mono_exon_length,
                )

        # Full-length end-coherence filter (FLAIR/TALON-style full-length read
        # support). Drops NOVEL multi-exon (>=2 intron) transcripts whose
        # fraction of full-length assigned reads (read genomic 5' AND 3' both
        # within fulllen_window_bp of the candidate's ends) is below
        # min_fulllen_fraction. The fulllen_frac METRIC itself is signal-free
        # (BAM primary-alignment spans only), but it is computed per-interval
        # over the quant assignment population (the m2_em/m1_em EM hard
        # assignments or, under argmax, the M2-tiebreak picks); candidates with
        # fulllen_frac < 0 (unreachable or never scored) are EXEMPT and never
        # dropped. gtf/fusion/mono exempt. ORTHOGONAL to the isoform-fraction
        # filter above (the two stack on the competitor metric). SIRV WARNING:
        # the default 0.1 is SIRV-tuned; re-tune or disable on real data.
        if _score_filter_on and self.config.min_fulllen_fraction > 0.0:
            drop_ids = fulllen_fraction_drops(
                aggregated, self.config.min_fulllen_fraction
            )
            if drop_ids:
                aggregated = {
                    cid: qr
                    for cid, qr in aggregated.items()
                    if cid not in drop_ids
                }
                logger.info(
                    "Dropped %d novel transcripts with full-length fraction < %.3f",
                    len(drop_ids),
                    self.config.min_fulllen_fraction,
                )

        # polyA + 5'-proximity candidate-retention filter. Drop a candidate
        # unless >= min_polya5p_reads of its assigned reads BOTH have a krill
        # whole-read polyA tail (qc PASS & length > min_polya_length) AND map
        # with their genomic 5' end within polya5p_window_bp of the candidate's
        # 5' end. GTF candidates are exempt by default (polya5p_exempt_gtf=True,
        # like fusion); only novel candidates are gated. Requires signal; no-ops
        # when absent.
        if (
            _score_filter_on
            and getattr(self.config, "min_polya5p_reads", 0) > 0
            and self.config.signal_path
            and Path(self.config.signal_path).exists()
        ):
            aggregated = self._apply_polya5p_filter(aggregated)

        # Resolve gene_ids from GTF annotation
        for cid, qr in aggregated.items():
            if qr.source == "gtf" and self._gtf_reader:
                tx = self._gtf_reader.get_transcript(cid)
                if tx:
                    qr.gene_id = tx.gene_id
            if not qr.gene_id:
                qr.gene_id = qr.candidate_id

        # Write GTF output
        if output_gtf:
            from fin.io.io_gtf import write_gtf

            write_gtf(aggregated, output_gtf)
            logger.info("Wrote GTF output: %s", output_gtf)

        # US-013: Additional output writers
        if output_tsv:
            from fin.io.io_tsv import write_scoring_tsv

            transcript_lengths = {
                cid: sum(end - start for start, end in qr.exons) if qr.exons else 0
                for cid, qr in aggregated.items()
            }
            write_scoring_tsv(aggregated, transcript_lengths, output_tsv)
            logger.info("Wrote TSV output: %s", output_tsv)

        if self.config.fusion_enabled and self.config.output_bedpe:
            from fin.io.io_bedpe import write_fusion_bedpe

            write_fusion_bedpe(aggregated, self.config.output_bedpe)
            logger.info("Wrote BEDPE output: %s", self.config.output_bedpe)

        logger.info(
            "Pipeline complete: %d transcripts quantified", len(aggregated)
        )
        return aggregated

    def process_interval(
        self, interval: GenomicInterval
    ) -> Optional[List[QuantResult]]:
        """Discover candidates and quantify a single interval via ``quant_mode``.

        Returns:
            - List[QuantResult] from the dispatched quant_mode engine.
            - None when the interval has no candidates or no reads.
        """
        work_dir = Path(self.config.work_dir) / interval.region_string.replace(":", "_").replace("-", "_")
        work_dir.mkdir(parents=True, exist_ok=True)

        # Get genome sequence for this chromosome
        chrom_seq = ""
        if self._genome_fasta and interval.chrom in self._genome_fasta:
            chrom_seq = self._genome_fasta[interval.chrom]

        # --- Phase 1: Candidate discovery ---
        candidate_set = discover_candidates(
            interval=interval,
            bam_path=self.config.bam_path,
            gtf_reader=self._gtf_reader,
            genome_fasta=chrom_seq,
            threshold=self.config.three_prime_threshold,
            min_novel_reads=self.config.min_novel_reads,
            chain_cluster=getattr(self.config, "chain_cluster_discovery", False),
            chain_cluster_wobble_bp=getattr(self.config, "chain_cluster_wobble_bp", 6),
            chain_cluster_cassette_max_exon_bp=getattr(
                self.config, "chain_cluster_cassette_max_exon_bp", 70),
            chain_cluster_fold_monoexon=(
                getattr(self.config, "chain_cluster_fold_monoexon", False)
                # defer the mono fold to post-EM when mono_resolve_post_em is on
                and not getattr(self.config, "mono_resolve_post_em", False)),
            chain_cluster_fold_span_guard=getattr(
                self.config, "chain_cluster_fold_span_guard", False),
            canonical_search_bp=self.config.canonical_search_bp,
            max_chains_per_read=self.config.max_chains_per_read,
            canonical_motifs=self.config.canonical_motifs,
            denovo_wobble_tol=getattr(self.config, "denovo_wobble_tol", 0),
            denovo_wobble_shadow_ratio=getattr(
                self.config, "denovo_wobble_shadow_ratio", 1.0),
            denovo_graph=getattr(self.config, "denovo_graph", False),
            denovo_graph_tol=getattr(self.config, "denovo_graph_tol", 6),
            denovo_graph_min_edge_reads=getattr(
                self.config, "denovo_graph_min_edge_reads", 2),
            denovo_graph_tss_brake=getattr(
                self.config, "denovo_graph_tss_brake", True),
            denovo_graph_tss_tol=getattr(self.config, "denovo_graph_tss_tol", 20),
            denovo_graph_tss_min_reads=getattr(
                self.config, "denovo_graph_tss_min_reads", 3),
            denovo_graph_tss_frac=getattr(
                self.config, "denovo_graph_tss_frac", 0.4),
        )

        # --- Phase 1.5: Fusion candidate augmentation (optional) ---
        if self.config.fusion_enabled:
            candidate_set = self._augment_with_fusion_candidates(
                candidate_set, interval
            )

        # --- Phase 1.6: Canonical-motif gate (Stage B) ---
        # Drop NOVEL multi-exon candidates whose intron chain isn't all-canonical.
        # GTF-passthrough and fusion candidates are exempt (annotated junctions
        # are trusted); mono candidates trivially pass.
        if getattr(self.config, "canonical_gate", False):
            self._apply_canonical_gate(candidate_set, chrom_seq)

        # --- Phase 1.7: Junction-dominance gate (PRE-EM, junction-first) ---
        # Drop NOVEL multi-exon candidates with a junction that is weak
        # (< min_reads observed) or non-dominant (a different observed junction
        # within window bp carries strictly more reads). Removes multi-read wobble
        # shadows before EM so they never compete for reads. gtf/fusion/mono
        # exempt. OFF (default) -> no drops (byte-identical).
        if getattr(self.config, "junction_dominance_filter", False):
            self._apply_junction_dominance_gate(candidate_set, interval)

        if candidate_set.num_candidates == 0:
            logger.info("No candidates for interval %s", interval.region_string)
            return None

        if not candidate_set.read_ids:
            logger.info("No reads for interval %s", interval.region_string)
            return None

        # Sort once; canonical read axis for every distance matrix.
        read_ids = sorted(candidate_set.read_ids)
        candidate_ids = candidate_set.candidate_ids()

        # Quantification engine dispatch (quant_mode). All three modes are
        # krill-only (no f5c CLI).
        quant_mode = getattr(self.config, "quant_mode", "m2_em")
        if quant_mode == "argmax":
            results = self._quant_argmax_keep(candidate_set, read_ids, interval)
        elif quant_mode == "m1_em":
            results = self._quant_m1_em(candidate_set, read_ids, interval)
        elif quant_mode == "m2_em":
            results = self._quant_m2_em(candidate_set, read_ids, interval)
        elif quant_mode == "cluster":
            results = self._quant_cluster(candidate_set, read_ids, interval)
        else:
            raise ValueError(f"unknown quant_mode: {quant_mode!r}")
        # Full-length end-coherence: compute fulllen_frac per candidate over its
        # assigned reads (the non-circular population; the fulllen METRIC itself
        # uses BAM spans only — no signal). Gated by the same switches as the
        # drop in _finalize_and_write so a filter-disabled run pays no BAM-fetch
        # cost and leaves the -1.0 sentinel untouched.
        _filter_on = getattr(self.config, "enable_score_filter", True)
        if (
            results
            and _filter_on
            and getattr(self.config, "min_fulllen_fraction", 0.0) > 0.0
        ):
            self._annotate_fulllen_frac(results, interval)
        return results

    def _quant_argmax_first(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> List[QuantResult]:
        """M1-first hard argmin assignment (production default; no signal, no EM).

        Each read is aligned against every candidate with the configured mappy
        preset (``MAPPY_PRESET``, default map-ont); the reconstructed map-ont AS
        (``score_hit``, single-indel > cap → rejected) is taken as the cell
        score. The read is assigned to its single best-AS candidate; AS TIES are
        broken by lowest candidate index — which is the implicit GTF prior, since
        ``discover_candidates`` returns ``gtf_passthrough + collapsed_novel +
        fusion`` (GTF first). Abundance = integer read count.

        This mirrors the ablation harness ``_build_m1`` + ``_quant_tie(mode=
        "first")`` so production and the harness produce the same assignment for
        the same candidate set.
        """
        import mappy

        from fin.scoring.mappy_preset import get_m1_preset
        from fin.scoring.mappy_score import score_hit

        cand_list = list(candidate_set.candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        n_c = len(cand_list)

        preset = get_m1_preset()
        aligners = [
            mappy.Aligner(seq=c.sequence, preset=preset) if c.sequence else None
            for c in cand_list
        ]

        counts = [0] * n_c
        assigned: List[List[str]] = [[] for _ in range(n_c)]
        for rid in read_ids:
            seq = read_sequences.get(rid)
            if not seq:
                continue
            best_as = None      # highest reconstructed AS over all candidates
            best_j = -1
            for j, aln in enumerate(aligners):
                if aln is None:
                    continue
                cell = None
                for h in aln.map(seq):
                    v = score_hit(h)
                    if v is None:
                        continue
                    if cell is None or v > cell:
                        cell = v
                if cell is None:
                    continue
                # strict '>' keeps the LOWEST-index candidate on an AS tie.
                if best_as is None or cell > best_as:
                    best_as = cell
                    best_j = j
            if best_j < 0:
                continue
            counts[best_j] += 1
            assigned[best_j].append(rid)

        quant_results: List[QuantResult] = []
        for j, cand in enumerate(cand_list):
            has = counts[j] > 0
            qr = QuantResult(
                candidate_id=cand.candidate_id,
                abundance=float(counts[j]),
                confidence=1.0 if has else 0.0,
                num_assigned_reads=len(assigned[j]),
                source=cand.source,
                chrom=cand.chrom,
                strand=cand.strand,
                start=cand.start,
                end=cand.end,
                exons=_exons_from_candidate(cand),
                assigned_read_ids=tuple(assigned[j]),
                breakpoint_left=cand.breakpoint_left,
                breakpoint_right=cand.breakpoint_right,
                fusion_junction=cand.fusion_junction,
            )
            qr.max_R = 1.0 if has else 0.0
            quant_results.append(qr)

        logger.info(
            "R1 M1-first interval %s: %d reads -> %d candidates",
            interval.region_string, len(read_ids), len(cand_list),
        )
        return quant_results

    def _quant_argmax_keep(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> List[QuantResult]:
        """M1-keep SPLIT assignment (production default; no signal, no EM).

        Each read is aligned against every candidate with the configured mappy
        preset; the reconstructed AS (``score_hit``) is the cell score. The read
        is assigned to ALL of its simultaneously-best-AS candidates (within a
        1e-9 tolerance). Abundance is split 1/K across the K tied winners so the
        per-read mass is conserved (== 1.0 per read total).

        Mirrors the ablation harness ``_build_m1`` + ``_quant_tie(mode="split")``
        (config "M1-split"): production and harness produce the same split
        assignment for the same candidate set.
        """
        import mappy

        from fin.scoring.mappy_preset import get_m1_preset
        from fin.scoring.mappy_score import score_hit

        cand_list = list(candidate_set.candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        n_c = len(cand_list)

        preset = get_m1_preset()
        aligners = [
            mappy.Aligner(seq=c.sequence, preset=preset) if c.sequence else None
            for c in cand_list
        ]

        # Optional M2 tie resolution (Stage M2-1): when a read is simultaneously
        # best-AS across >=2 candidates, the validated junction-window mean-NLL
        # metric picks the single true wobble sibling and (if confident) takes the
        # read's FULL mass instead of the 1/K split. Lazily build a shared non-HMM
        # krill aligner; auto-skip the whole leg if signal is absent.
        m2_krill = None
        m2_gpu = False
        m2_on = bool(getattr(self.config, "m2_tiebreak", False)) and bool(
            self.config.signal_path
        )
        m2_threads = 0
        if m2_on:
            try:
                import krill

                from fin.scoring.krill_aligner import (
                    krill_thread_count,
                    make_krill_aligner,
                )

                m2_threads = krill_thread_count()
                m2_krill, m2_gpu = make_krill_aligner(
                    krill, self.config.krill_pore, self.config.use_gpu,
                    hmm_confidence=False, num_thread=m2_threads,
                )
            except Exception as exc:  # krill not importable -> keep 1/K split
                logger.warning("M2 tiebreak disabled (krill import failed): %s", exc)
                m2_krill = None
            if m2_krill is None:  # signal stack unavailable -> keep 1/K split
                logger.warning("M2 tiebreak disabled (krill init failed)")
                m2_on = False

        counts = [0.0] * n_c
        assigned: List[List[str]] = [[] for _ in range(n_c)]
        # Track the max single-read weight any candidate receives (== max_R).
        max_weight = [0.0] * n_c
        n_m2_override = 0
        for rid in read_ids:
            seq = read_sequences.get(rid)
            if not seq:
                continue
            cells: List[Optional[float]] = [None] * n_c
            best_as = None
            for j, aln in enumerate(aligners):
                if aln is None:
                    continue
                cell = None
                for h in aln.map(seq):
                    v = score_hit(h)
                    if v is None:
                        continue
                    if cell is None or v > cell:
                        cell = v
                if cell is None:
                    continue
                cells[j] = cell
                if best_as is None or cell > best_as:
                    best_as = cell
            if best_as is None:
                continue
            tied = [
                j for j in range(n_c)
                if cells[j] is not None and cells[j] >= best_as - 1e-9
            ]
            if not tied:
                continue
            # M2 single-winner override on a genuine (>=2) tie set.
            if m2_on and len(tied) >= 2:
                from fin.scoring.m2_junction_nll import m2_resolve_tie

                tied_cands = [cand_list[j] for j in tied]
                tied_aligners = [aligners[j] for j in tied]
                best_local, margin, scored_local = m2_resolve_tie(
                    rid, seq, tied_cands, self.config.signal_path,
                    pore=self.config.krill_pore,
                    junction_k=self.config.m2_tiebreak_junction_k,
                    krill_aligner=m2_krill, mappy_aligners=tied_aligners,
                    return_scored=True, use_gpu=m2_gpu, num_thread=m2_threads,
                )
                if best_local is not None and margin >= self.config.m2_tiebreak_margin:
                    j = tied[best_local]
                    counts[j] += 1.0
                    assigned[j].append(rid)
                    if 1.0 > max_weight[j]:
                        max_weight[j] = 1.0
                    n_m2_override += 1
                    continue
                # Score-gated fallback split (mk1): restrict the 1/K split to the
                # candidates eventalign could score. If none scored, keep the full
                # tied set (== M1). Recall-safe: the read is never dropped.
                if (
                    getattr(self.config, "m2_tie_scoregate_split", False)
                    and scored_local
                ):
                    tied = [tied[k] for k in scored_local]
            w = 1.0 / len(tied)
            for j in tied:
                counts[j] += w
                assigned[j].append(rid)
                if w > max_weight[j]:
                    max_weight[j] = w

        quant_results: List[QuantResult] = []
        for j, cand in enumerate(cand_list):
            has = counts[j] > 0.0
            qr = QuantResult(
                candidate_id=cand.candidate_id,
                abundance=float(counts[j]),
                confidence=1.0 if has else 0.0,
                num_assigned_reads=len(assigned[j]),
                source=cand.source,
                chrom=cand.chrom,
                strand=cand.strand,
                start=cand.start,
                end=cand.end,
                exons=_exons_from_candidate(cand),
                assigned_read_ids=tuple(assigned[j]),
                breakpoint_left=cand.breakpoint_left,
                breakpoint_right=cand.breakpoint_right,
                fusion_junction=cand.fusion_junction,
            )
            qr.max_R = max_weight[j]
            quant_results.append(qr)

        logger.info(
            "R1 M1-keep(split) interval %s: %d reads -> %d candidates "
            "(M2 overrides=%d)",
            interval.region_string, len(read_ids), len(cand_list), n_m2_override,
        )
        return quant_results

    def _quant_cluster(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> List[QuantResult]:
        """Within-cluster read assignment (quant_mode="cluster").

        Each generation cluster (``CandidateSet.clusters``) is scored in isolation:
        a read is aligned ONLY against that cluster's members (M1), straddling ties
        are resolved by the summed-LLR junction signal (M2), and containment/
        truncation ties fall to the cluster's main peak. Per-cluster assignment is
        delegated to :func:`fin.pipeline.cluster_quant.assign_cluster_reads` (pure/
        signal-agnostic); the M2 tiebreak is injected via ``m2_resolve_tie``. Cluster
        results are folded back into per-candidate abundance/survival and emitted as
        QuantResults (mirrors ``_quant_argmax_keep``).
        """
        import mappy

        from fin.pipeline.cluster_quant import assign_cluster_reads
        from fin.scoring.m2_junction_nll import m2_resolve_tie
        from fin.scoring.mappy_preset import get_m1_preset
        from fin.scoring.mappy_score import score_hit

        cand_list = list(candidate_set.candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        n_c = len(cand_list)

        # candidate_set.clusters holds candidate_ids per generation cluster. Translate
        # to CURRENT list indices (ids dropped by post-discovery gates are skipped).
        idx_of = {c.candidate_id: j for j, c in enumerate(cand_list)}
        if candidate_set.clusters is None:
            # Not chain-cluster mode: run unscoped over the whole candidate list.
            logger.warning(
                "quant_mode=cluster: interval %s has no cluster grouping "
                "(clusters=None); treating all %d candidates as one cluster",
                interval.region_string, n_c,
            )
            clusters = [list(range(n_c))]
        else:
            clusters = [
                [idx_of[cid] for cid in ids if cid in idx_of]
                for ids in candidate_set.clusters
            ]
            # Candidates added AFTER discovery (e.g. fusion augmentation) are in
            # cand_list but not in any generation cluster; group them into ONE extra
            # cluster so they are scored (never silently zero) AND a read shared by
            # several of them competes once (one singleton each would let the same
            # read take weight 1.0 in every singleton -> inflated abundance).
            referenced = {j for ml in clusters for j in ml}
            extra = [j for j in range(n_c) if j not in referenced]
            if extra:
                clusters.append(extra)

        # Read genomic spans over the interval (primary mapped alignments only;
        # first occurrence per read id wins). Mirrors _annotate_fulllen_frac.
        spans: Dict[str, Tuple[int, int]] = {}
        with pysam.AlignmentFile(self.config.bam_path, "rb") as bam:
            for r in bam.fetch(interval.chrom, interval.start, interval.end):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                rid = r.query_name
                if rid is None or rid in spans:
                    continue
                rs = r.reference_start
                re = r.reference_end
                if rs is None or re is None:
                    continue
                spans[rid] = (int(rs), int(re))

        # One shared non-HMM krill aligner for the whole interval; when signal is
        # absent m2_resolve stays None (straddling ties defer, which is fine).
        krill_aligner = None
        eff_gpu = False
        krill_threads = 0
        if self.config.signal_path and getattr(self.config, "cluster_use_m2", True):
            try:
                import krill

                from fin.scoring.krill_aligner import (
                    krill_thread_count,
                    make_krill_aligner,
                )

                krill_threads = krill_thread_count()
                krill_aligner, eff_gpu = make_krill_aligner(
                    krill, self.config.krill_pore, self.config.use_gpu,
                    hmm_confidence=False, num_thread=krill_threads,
                )
            except Exception as exc:  # krill unavailable -> defer straddling ties
                logger.warning(
                    "quant_mode=cluster: krill unavailable, M2 tiebreak off: %s", exc)
                krill_aligner = None

        preset = get_m1_preset()
        read_id_set = set(read_ids)

        # Accumulators over ALL candidates (global index j).
        counts = [0.0] * n_c
        assigned: List[List[str]] = [[] for _ in range(n_c)]
        max_weight = [0.0] * n_c
        survivor = [False] * n_c

        # Genome chromosome sequence for 5'-TSS short-isoform recovery (empty ->
        # recovery skipped for this interval).
        genome_seq = (
            self._genome_fasta[interval.chrom]
            if self._genome_fasta and interval.chrom in self._genome_fasta
            else ""
        )
        shadows_map = candidate_set.shadows or {}

        tot_unique = tot_containment = tot_m2_confident = tot_deferred = 0
        recovered_total = 0
        for members_idx in clusters:
            if not members_idx:
                continue
            member_candidates = [cand_list[i] for i in members_idx]
            local_aligners = [
                mappy.Aligner(seq=c.sequence, preset=preset) if c.sequence else None
                for c in member_candidates
            ]
            member_chains = [
                tuple(cand_list[i].intron_chain.introns) for i in members_idx
            ]

            # Cluster reads = union of members' supporting reads, restricted to this
            # interval's read set and reads that carry a sequence.
            cluster_read_ids: set = set()
            for i in members_idx:
                cluster_read_ids.update(cand_list[i].supporting_read_ids)
            cluster_read_ids &= read_id_set
            cluster_read_ids = {
                rid for rid in cluster_read_ids if read_sequences.get(rid)
            }

            # M1 within-cluster: best AS per (read, member), local member index.
            as_rows: Dict[str, Dict[int, float]] = {}
            for rid in cluster_read_ids:
                seq = read_sequences.get(rid)
                if not seq:
                    continue
                row: Dict[int, float] = {}
                for local_j, aln in enumerate(local_aligners):
                    if aln is None:
                        continue
                    cell = None
                    for h in aln.map(seq):
                        v = score_hit(h)
                        if v is None:
                            continue
                        if cell is None or v > cell:
                            cell = v
                    if cell is not None:
                        row[local_j] = cell
                if row:
                    as_rows[rid] = row

            # M2 resolver for straddling ties (summed-LLR); None when no krill.
            m2_resolve = None
            if krill_aligner is not None:
                def m2_resolve(rid, tie_local, _mc=member_candidates,
                               _la=local_aligners):
                    tied_cands = [_mc[j] for j in tie_local]
                    tied_aligners = [_la[j] for j in tie_local]
                    out = m2_resolve_tie(
                        rid, read_sequences[rid], tied_cands,
                        self.config.signal_path, pore=self.config.krill_pore,
                        krill_aligner=krill_aligner, mappy_aligners=tied_aligners,
                        use_gpu=eff_gpu, num_thread=krill_threads,
                        metric="summed_llr",
                    )
                    best, margin = out
                    if best is None:
                        return None
                    return (tie_local[best], margin)

            a = assign_cluster_reads(
                member_chains, as_rows,
                {rid: spans[rid] for rid in as_rows if rid in spans},
                m2_resolve=m2_resolve,
                min_support=self.config.cluster_min_support,
                m1_tie_margin=self.config.cluster_m1_tie_margin,
            )

            # Fold local results back to GLOBAL candidate indices.
            for rid, row in a.weights.items():
                for local_j, w in row.items():
                    gi = members_idx[local_j]
                    counts[gi] += w
                    assigned[gi].append(rid)
                    if w > max_weight[gi]:
                        max_weight[gi] = w
            for local_j in a.survivors:
                survivor[members_idx[local_j]] = True

            tot_unique += a.n_unique
            tot_containment += a.n_containment
            tot_m2_confident += a.n_m2_confident
            tot_deferred += a.n_deferred

            # 5'-TSS short-isoform recovery: each EM-confirmed survivor member may
            # have folded SHADOW sub-chains (exact 5'-truncations pooled into it at
            # generation). Re-test each shadow's truncation boundary against the
            # member's read-5'-end pileup; a SHARP peak (real TSS, not a degradation
            # ramp) re-activates the short isoform as a novel candidate. Mass is
            # conserved: the recovered peak-excess is subtracted from the parent.
            if genome_seq:
                for lj in sorted(a.survivors):
                    gj = members_idx[lj]
                    cand = cand_list[gj]
                    cid = cand.candidate_id
                    strand = cand.strand
                    shadows = shadows_map.get(cid)
                    if not shadows:
                        continue

                    # 5' end of a read (strand-aware): ref_start for +, ref_end for -.
                    def _p5(rid, _s=strand):
                        rs, re = spans[rid]
                        return rs if _s == "+" else re

                    # Member's EM-weighted assigned reads (this member's column).
                    member_reads = {
                        rid: row[lj]
                        for rid, row in a.weights.items()
                        if lj in row
                    }
                    ends = [
                        (_p5(rid), w)
                        for rid, w in member_reads.items()
                        if rid in spans
                    ]
                    if len(ends) < 8:
                        continue

                    # Build TSS proposals (median 5' end of each shadow's reads).
                    proposals: List[int] = []
                    prop_shadow: Dict[int, Tuple[tuple, tuple]] = {}
                    for schain, sreads in shadows:
                        s_ends = sorted(_p5(rid) for rid in sreads if rid in spans)
                        if not s_ends:
                            continue
                        tss = s_ends[len(s_ends) // 2]
                        prev = prop_shadow.get(tss)
                        if prev is None or len(sreads) > len(prev[1]):
                            prop_shadow[tss] = (schain, sreads)
                            if prev is None:
                                proposals.append(tss)

                    if not proposals:
                        continue

                    from fin.candidates.isoform_recovery import recover_5p_peaks

                    peaks = recover_5p_peaks(ends, list(proposals))
                    for peak in peaks:
                        # Nearest proposal within 20 bp.
                        near = min(
                            prop_shadow.keys(),
                            key=lambda p: abs(p - peak.pos),
                        )
                        if abs(near - peak.pos) > 20:
                            continue
                        schain, sreads = prop_shadow[near]

                        # Short-isoform genomic span: 5' end -> TSS, 3' end = parent's.
                        if strand == "+":
                            start, end = peak.pos, cand.end
                        else:
                            start, end = cand.start, peak.pos
                        if end <= start:
                            continue

                        from fin.candidates.discovery import (
                            _build_spliced_sequence,
                            _generate_novel_id,
                        )
                        from fin.candidates.dataclasses import (
                            IntronChain,
                            TranscriptCandidate,
                        )

                        seq = _build_spliced_sequence(
                            genome_seq, start, end,
                            IntronChain(introns=schain), strand,
                        )
                        if not seq:
                            continue

                        new_id = _generate_novel_id()
                        new_cand = TranscriptCandidate(
                            candidate_id=new_id,
                            intron_chain=IntronChain(introns=schain),
                            three_prime_pos=(end if strand == "+" else start),
                            sequence=seq,
                            source="novel",
                            supporting_read_ids=set(sreads),
                            chrom=cand.chrom,
                            strand=strand,
                            start=start,
                            end=end,
                        )
                        # Append in lockstep with the parallel accumulators so the
                        # emission loop (which enumerates cand_list) picks it up.
                        cand_list.append(new_cand)
                        counts.append(peak.excess)
                        assigned.append(list(sreads))
                        max_weight.append(1.0)
                        survivor.append(True)
                        # Conserve mass: remove the recovered excess from the parent.
                        counts[gj] = max(0.0, counts[gj] - peak.excess)
                        recovered_total += 1

        if recovered_total:
            logger.info(
                "R1 cluster interval %s: recovered %d 5'-TSS short isoforms",
                interval.region_string, recovered_total,
            )

        quant_results: List[QuantResult] = []
        for j, cand in enumerate(cand_list):
            # Honor the per-cluster survival gate: a member whose total assigned
            # weight is below cluster_min_support is NOT a survivor -> zero it so the
            # threshold is a real gate (not overridden here to confidence 1.0).
            kept = survivor[j]
            qr = QuantResult(
                candidate_id=cand.candidate_id,
                abundance=float(counts[j]) if kept else 0.0,
                confidence=1.0 if kept else 0.0,
                num_assigned_reads=len(assigned[j]) if kept else 0,
                source=cand.source,
                chrom=cand.chrom,
                strand=cand.strand,
                start=cand.start,
                end=cand.end,
                exons=_exons_from_candidate(cand),
                assigned_read_ids=tuple(assigned[j]) if kept else (),
                breakpoint_left=cand.breakpoint_left,
                breakpoint_right=cand.breakpoint_right,
                fusion_junction=cand.fusion_junction,
            )
            qr.max_R = max_weight[j]
            quant_results.append(qr)

        logger.info(
            "R1 cluster interval %s: %d clusters, %d reads -> %d candidates "
            "(unique=%d containment=%d m2_confident=%d deferred=%d)",
            interval.region_string, len(clusters), len(read_ids), len(cand_list),
            tot_unique, tot_containment, tot_m2_confident, tot_deferred,
        )
        return quant_results

    def _annotate_fulllen_frac(
        self,
        results: List[QuantResult],
        interval: GenomicInterval,
    ) -> None:
        """Store fulllen_frac on NOVEL multi-exon QuantResults (in place).

        full-length read support: the fraction of a candidate's argmax-assigned
        reads whose primary genomic 5' AND 3' alignment ends both fall within
        ``fulllen_window_bp`` of the candidate's genomic 5'/3' ends. Signal-free
        (BAM primary spans only). A single fetch over the interval builds the
        read-end map shared by every candidate. Candidates with fewer than
        ``fulllen_min_reads`` reads carrying a span keep the -1.0 sentinel
        (unreachable -> never dropped). gtf/fusion/mono candidates are skipped.
        """
        novel_multi = [
            qr for qr in results
            if qr.source == "novel" and len(qr.exons) >= 3
        ]
        if not novel_multi:
            return
        window = getattr(self.config, "fulllen_window_bp", 25)
        min_reads = getattr(self.config, "fulllen_min_reads", 4)
        # Single BAM fetch over the interval -> genomic read-end map (primary
        # mapped alignments only; first occurrence per read id wins).
        read_ends: Dict[str, Tuple[int, int]] = {}
        with pysam.AlignmentFile(self.config.bam_path, "rb") as bam:
            for r in bam.fetch(interval.chrom, interval.start, interval.end):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                rid = r.query_name
                if rid is None or rid in read_ends:
                    continue
                rs = r.reference_start
                re = r.reference_end
                if rs is None or re is None:
                    continue
                read_ends[rid] = (int(rs), int(re))
        for qr in novel_multi:
            qr.fulllen_frac = compute_fulllen_frac(
                qr, read_ends, window, min_reads
            )

    def _fetch_read_seqs_and_ends(
        self,
    ) -> Tuple[Dict[str, str], Dict[str, Tuple[int, int]]]:
        """Single BAM pass -> ({rid: query_sequence}, {rid: (ref_start, ref_end)}).

        Primary mapped alignments only; first occurrence per read id wins. The
        sequences feed krill whole-read polyA; the genomic spans feed 5'
        proximity (mirrors ``_annotate_fulllen_frac``'s read-end map).
        """
        read_seqs: Dict[str, str] = {}
        read_ends: Dict[str, Tuple[int, int]] = {}
        with pysam.AlignmentFile(self.config.bam_path, "rb") as bam:
            for r in bam.fetch():
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                rid = r.query_name
                if rid is None or rid in read_ends:
                    continue
                rs, re = r.reference_start, r.reference_end
                if rs is None or re is None:
                    continue
                read_ends[rid] = (int(rs), int(re))
                if r.query_sequence:
                    read_seqs[rid] = r.query_sequence
        return read_seqs, read_ends

    def _apply_polya5p_filter(
        self, aggregated: Dict[str, QuantResult]
    ) -> Dict[str, QuantResult]:
        """Drop candidates failing the polyA + 5'-proximity gate.

        Runs a krill whole-read polyA pass over all reads, then removes any
        gated candidate with fewer than ``min_polya5p_reads`` reads that both
        have a confident polyA tail and map 5'-flush to the candidate. Fusion
        candidates are always exempt; GTF candidates are exempt by default
        (``polya5p_exempt_gtf``), so only novel candidates are gated. No-ops if
        krill returns nothing.
        """
        from fin.scoring.polya import compute_polya

        read_seqs, read_ends = self._fetch_read_seqs_and_ends()
        polya_map = compute_polya(
            read_seqs,
            self.config.signal_path,
            pore=self.config.krill_pore,
            use_gpu=self.config.use_gpu,
        )
        if not polya_map:
            logger.warning(
                "polyA+5' filter: krill returned no polyA estimates; skipping"
            )
            return aggregated

        drop_ids = polya5p_drops(
            aggregated,
            polya_map,
            read_ends,
            window=self.config.polya5p_window_bp,
            min_polya_len=self.config.min_polya_length,
            min_reads=self.config.min_polya5p_reads,
            exempt_gtf=self.config.polya5p_exempt_gtf,
        )
        if not drop_ids:
            return aggregated

        n_gtf = sum(
            1 for cid in drop_ids if aggregated[cid].source == "gtf"
        )
        n_novel = len(drop_ids) - n_gtf
        aggregated = {
            cid: qr for cid, qr in aggregated.items() if cid not in drop_ids
        }
        logger.info(
            "Dropped %d gtf + %d novel transcripts failing polyA+5' (>= %d reads)",
            n_gtf,
            n_novel,
            self.config.min_polya5p_reads,
        )
        return aggregated

    def _adaptive_sigma(self, dist_read_to_tx: np.ndarray) -> float:
        """Data-adaptive EM sigma: median per-read range of the read×tx distance,
        clipped to [em_sigma_min, em_sigma_max]. Falls back to em_sigma.

        With per-event-normalized krill distances the absolute scale is dataset-
        dependent (~0.5-5 nats); a fixed sigma=1.0 often collapses R to one-hot.
        """
        n_reads, n_tx = dist_read_to_tx.shape
        if n_tx >= 2 and n_reads > 0:
            d_max = dist_read_to_tx.max(axis=1)
            d_min = dist_read_to_tx.min(axis=1)
            adaptive = float(np.median(d_max - d_min))
            return float(np.clip(
                adaptive if adaptive > 0 else self.config.em_sigma,
                getattr(self.config, "em_sigma_min", 0.05),
                getattr(self.config, "em_sigma_max", 50.0),
            ))
        return self.config.em_sigma

    def _eff_lengths(self, cand_list) -> Optional[np.ndarray]:
        """Per-candidate spliced effective length in cand_list column order, for
        Salmon-style abundance-feedback length normalization. Returns None unless
        abundance_feedback AND abundance_length_norm are both on. Lengths are
        floored at 1.0 so a degenerate (zero-length) candidate cannot divide-by-0.
        """
        if not (self.config.abundance_feedback and self.config.abundance_length_norm):
            return None
        lengths = np.array(
            [
                max(1.0, float(sum(e - s for s, e in _exons_from_candidate(c))))
                for c in cand_list
            ],
            dtype=float,
        )
        return lengths

    def _quant_m1_em(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> Optional[List[QuantResult]]:
        """quant_mode='m1_em': EM seeded by the M1 mappy AS-gap distance (β=0,
        no signal coherence). Pure-alignment soft assignment over krill-free
        mappy distances."""
        from fin.ablation.mappy_argmax import mappy_multimap_responsibilities

        cand_list = list(candidate_set.candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        reads_iter = [(rid, read_sequences.get(rid, "")) for rid in read_ids]
        reads_iter = [(rid, seq) for rid, seq in reads_iter if seq]
        R_mm, kept_read_ids = mappy_multimap_responsibilities(reads_iter, cand_list)
        if R_mm.size == 0:
            return []

        read_seqs = {rid: seq for rid, seq in reads_iter}
        n_reads_em = len(kept_read_ids)
        n_cands_em = len(cand_list)
        max_iter_em = (
            self.config.em_max_iter_override
            if self.config.em_max_iter_override is not None
            else self.config.em_max_iter
        )

        m1 = compute_mappy_distance(read_seqs, cand_list, kept_read_ids)
        m2_dummy = np.zeros((n_reads_em, n_cands_em), dtype=np.float32)
        m3_dummy = np.zeros((n_reads_em, n_reads_em), dtype=np.float32)
        dist_read_to_tx, _, _ = build_em_matrices(
            "m1", m1, m2_dummy, m3_dummy, em_beta=0.0,
        )
        dist_read_to_read = np.zeros((n_reads_em, n_reads_em), dtype=np.float32)
        R, hard_assignments, _ = em_with_coherence(
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_read=dist_read_to_read,
            sigma=self.config.em_sigma,
            beta=0.0,
            max_iter=max_iter_em,
            tol=self.config.em_tol,
            verbose=False,
            use_gpu=self.config.use_gpu,
            abundance_feedback=self.config.abundance_feedback,
            abundance_length_norm=self.config.abundance_length_norm,
            eff_lengths=self._eff_lengths(cand_list),
        )

        if self.config.krill_tiebreak:
            R = krill_tiebreak(
                R=R, read_ids=kept_read_ids, read_seqs=read_seqs,
                candidates=cand_list, signal_path=self.config.signal_path,
                pore=self.config.krill_pore,
                ambig_threshold=self.config.tiebreak_ambig_threshold,
                use_gpu=self.config.use_gpu,
            )
            hard_assignments = R.argmax(axis=1)

        quant_results = quantify_transcripts(
            R, hard_assignments, cand_list, kept_read_ids
        )
        for j, qr in enumerate(quant_results):
            qr.max_R = float(R[:, j].max()) if R.shape[0] > 0 else 0.0
        logger.info(
            "m1_em interval %s: %d reads -> %d candidates",
            interval.region_string, len(read_ids), len(cand_list),
        )
        return quant_results

    def _tie_nll(self, kept_read_ids, read_seqs, cand_list, aligners, raw):
        """Per-read best-AS tie junction-NLL (the m2_resolve_tie niche).

        For each read whose raw mappy AS is simultaneously best across >=2
        candidates, build the class junction-discrimination window over that tie
        set and score each tied candidate with the per-event junction mean-NLL
        (non-HMM krill eventalign over the window). Only tied, window-scorable
        cells are touched — never the dense read×candidate matrix.

        Returns (nlls_by_read {i: {j: nll}}, ties_by_read {i: [j,...]},
        n_ties, n_refined, cover_by_read {i: bool}). ``cover_by_read`` is populated
        only when the diff-region coverage gate (config.m2_diff_cover_gate) is ON;
        otherwise it is empty.
        """
        import krill

        from fin.scoring.krill_aligner import (
            krill_thread_count,
            make_krill_aligner,
        )
        from fin.scoring.m2_junction_nll import (
            _mean_nll_in_gset,
            class_junction_window_set,
            event_genomic_positions,
            read_cand_mean_nll,
            read_straddles,
            wobble_diff_spans,
        )
        from fin.scoring.mappy_score import score_hit

        gate = bool(getattr(self.config, "m2_diff_cover_gate", False))
        num_threads = krill_thread_count()
        n_c = len(cand_list)
        nlls_by_read: Dict[int, Dict[int, float]] = {}
        ties_by_read: Dict[int, List[int]] = {}
        cover_by_read: Dict[int, bool] = {}
        krill_aligner, eff_gpu = make_krill_aligner(
            krill, self.config.krill_pore, self.config.use_gpu,
            hmm_confidence=False, num_thread=num_threads,
        )
        if krill_aligner is None:
            # No krill signal backend. Gate OFF: preserve the legacy behavior
            # (empty ties -> _quant_m2_em leaves rows MISSING, byte-identical).
            # Gate ON: still expose the AS tie sets (derivable from ``raw`` alone,
            # no signal needed) so the gate's unique-best / drop logic is
            # well-defined; without them every row would be all-MISSING and the
            # EM would normalize each read to a spurious uniform vote.
            if gate:
                for i in range(len(kept_read_ids)):
                    row = raw[i]
                    finite = np.isfinite(row)
                    if not finite.any():
                        continue
                    best_as = float(row[finite].max())
                    ties_by_read[i] = [
                        j for j in range(n_c) if finite[j] and row[j] >= best_as - 1e-9
                    ]
            return nlls_by_read, ties_by_read, 0, 0, cover_by_read

        k = self.config.m2_tiebreak_junction_k
        pore = self.config.krill_pore
        sig_path = self.config.signal_path

        # --- Pass 1: tie detection + per-(read, candidate) mappy best-hit slices.
        # Each tied cell is reduced to the SAME (sliced sequence, start offset)
        # the singular read_cand_mean_nll would have eventaligned, so the batched
        # eventalign below is byte-identical to the per-call path -- only the
        # Python/GIL/per-call overhead is removed.
        gset_by_read: Dict[int, set] = {}
        spans_by_read: Dict[int, List[Tuple[int, int]]] = {}  # diff-cover gate spans
        reads_variants: Dict[str, Dict[str, str]] = {}  # rid -> {cand_id: seq}
        starts: Dict[str, Dict[str, int]] = {}          # rid -> {cand_id: r_st}
        cidj_by_read: Dict[str, Dict[str, int]] = {}     # rid -> {cand_id: col j}
        n_ties = 0
        for i, rid in enumerate(kept_read_ids):
            row = raw[i]
            finite = np.isfinite(row)
            if not finite.any():
                continue
            best_as = float(row[finite].max())
            tie = [j for j in range(n_c) if finite[j] and row[j] >= best_as - 1e-9]
            ties_by_read[i] = tie
            if len(tie) < 2:
                continue
            n_ties += 1
            gset = class_junction_window_set(
                [cand_list[j] for j in tie], flank=_M2_EM_FLANK, k=k
            )
            if not gset:
                continue
            if gate:
                spans_by_read[i] = wobble_diff_spans(
                    [cand_list[j] for j in tie], flank=_M2_EM_FLANK
                )
            seq = read_seqs[rid]
            per_seqs: Dict[str, str] = {}
            per_starts: Dict[str, int] = {}
            per_cidj: Dict[str, int] = {}
            for j in tie:
                aln = aligners[j]
                if aln is None:
                    continue
                best_hit = None
                best_hit_as = None
                for h in aln.map(seq):
                    v = score_hit(h)
                    if v is None:
                        continue
                    if best_hit_as is None or v > best_hit_as:
                        best_hit_as, best_hit = v, h
                if best_hit is None:
                    continue
                cand = cand_list[j]
                per_seqs[cand.candidate_id] = cand.sequence[best_hit.r_st:best_hit.r_en]
                per_starts[cand.candidate_id] = best_hit.r_st
                per_cidj[cand.candidate_id] = j
            if not per_seqs:
                continue
            gset_by_read[i] = gset
            reads_variants[rid] = per_seqs
            starts[rid] = per_starts
            cidj_by_read[rid] = per_cidj

        if not reads_variants:
            return nlls_by_read, ties_by_read, n_ties, 0, cover_by_read

        # --- Pass 2: ONE batched eventalign over every tied pair (per-pair start
        # via the {label: offset} form); per-read singular fallback on absence or
        # batch failure. Both branches score with the same _mean_nll_in_gset.
        rid_to_i = {rid: i for i, rid in enumerate(kept_read_ids)}
        n_refined = 0
        use_batch = hasattr(krill, "align_reads_variants")
        batch_out = None
        if use_batch:
            try:
                batch_out = krill.align_reads_variants(
                    sig_path, reads_variants, pore=pore,
                    use_gpu=eff_gpu, num_thread=num_threads,
                    aligner=krill_aligner, start=starts,
                )
            except Exception as exc:  # noqa: BLE001 - krill raises broad errors
                logger.warning("M2 batch eventalign failed (%s); per-read fallback", exc)
                use_batch = False

        if use_batch and batch_out is not None:
            for rid, res_list in batch_out.items():
                i = rid_to_i.get(rid)
                if i is None:
                    continue
                gset = gset_by_read.get(i)
                if gset is None:
                    continue
                by_label = {x.get("variant_label"): x for x in res_list}
                nlls: Dict[int, float] = {}
                spans = spans_by_read.get(i) if gate else None
                span_cov = [False] * len(spans) if spans else None
                for cid, j in cidj_by_read.get(rid, {}).items():
                    res = by_label.get(cid)
                    if res is None or res.get("status", -1) != 0:
                        continue
                    nll, n_ev = _mean_nll_in_gset(res, cand_list[j], gset)
                    if n_ev > 0 and np.isfinite(nll):
                        nlls[j] = float(nll)
                    if span_cov is not None:
                        egp = event_genomic_positions(res, cand_list[j])
                        if egp:
                            for s_idx, (lo, hi) in enumerate(spans):
                                if not span_cov[s_idx] and read_straddles(egp, lo, hi):
                                    span_cov[s_idx] = True
                if nlls:
                    n_refined += 1
                    nlls_by_read[i] = nlls
                    if gate:
                        # covered iff every wobbling diff span is straddled by the
                        # read on >=1 candidate (empty spans -> vacuously covered).
                        cover_by_read[i] = all(span_cov) if span_cov else True
            return nlls_by_read, ties_by_read, n_ties, n_refined, cover_by_read

        for rid, per_cidj in cidj_by_read.items():
            i = rid_to_i.get(rid)
            if i is None:
                continue
            gset = gset_by_read.get(i)
            seq = read_seqs[rid]
            nlls = {}
            for j in per_cidj.values():
                nll, n_ev = read_cand_mean_nll(
                    rid, seq, cand_list[j], [], krill_aligner, aligners[j],
                    sig_path, pore, gset=gset, use_gpu=eff_gpu,
                    num_thread=num_threads,
                )
                if n_ev > 0 and np.isfinite(nll):
                    nlls[j] = float(nll)
            if nlls:
                n_refined += 1
                nlls_by_read[i] = nlls
        # Per-read fallback path does not expose eventalign positions, so the
        # coverage map is left empty here; the d_tx builder defaults missing reads
        # to covered=True (recall-safe: never drops on the rare fallback path).
        return nlls_by_read, ties_by_read, n_ties, n_refined, cover_by_read

    def _observed_junctions(self, interval: GenomicInterval):
        """Strand-keyed {strand: Counter{(donor, acceptor): n_reads}} of intron
        junctions directly observed in the interval's primary-read CIGARs.

        Delegates to ``evidence.compute_observed_junctions`` (the Evidence layer) with a
        single-entry per-interval memo: this map is consumed by up to three gates within one
        interval (cluster-recheck GTF guard, Lever-2 per-junction support gate,
        junction-dominance / guided support gate), which used to each re-open the BAM;
        computing it once here is byte-identical. Returns None on empty/unreadable BAM so every
        support gate self-disables (fail-open). See fin/pipeline/evidence.py for the caveat.
        """
        from fin.pipeline.evidence import compute_observed_junctions

        # Memo key = region only: the map fetches chrom:start-end and buckets by READ
        # strand, so it is independent of the interval's own strand (two same-region
        # intervals share it). region_string is present on every interval incl. test stubs.
        key = interval.region_string
        memo = getattr(self, "_obs_junc_memo", None)
        if memo is not None and memo[0] == key:
            return memo[1]
        value = compute_observed_junctions(self.config.bam_path, interval)
        self._obs_junc_memo = (key, value)
        return value

    def _quant_m2_em(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> Optional[List[QuantResult]]:
        """quant_mode='m2_em' (production default): pure tie-break junction-NLL EM.

        M1/AS is used ONLY as a hard selector — it picks each read's best-AS tie
        set and masks mappability — and the per-event junction NLL is the SOLE
        graded distance over that tie set (m2_resolve_tie semantics, NOT the dense
        read×candidate matrix, which collapses). Optionally adds the M3 read×read
        junction-window DTW coherence when ``config.m3_coherence`` is True
        (β=em_beta); OFF by default because the pairwise DTW is expensive. All
        signal scoring is in-memory krill — no f5c CLI. (``m4_source`` does not
        affect this path.)
        """
        import mappy

        from fin.scoring.m2_junction_nll import MISSING
        from fin.scoring.mappy_preset import get_m1_preset
        from fin.scoring.mappy_score import score_hit

        cand_list = list(candidate_set.candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        reads_iter = [(rid, read_sequences.get(rid, "")) for rid in read_ids]
        reads_iter = [(rid, seq) for rid, seq in reads_iter if seq]
        if not reads_iter:
            return []

        read_seqs = {rid: seq for rid, seq in reads_iter}
        n_c = len(cand_list)
        max_iter_em = (
            self.config.em_max_iter_override
            if self.config.em_max_iter_override is not None
            else self.config.em_max_iter
        )

        # --- Per-candidate mappy aligners + raw AS matrix (tie detection). This
        #     single AS pass also yields the kept-read set, so the previously
        #     separate mappy_multimap_responsibilities pass (which recomputed the
        #     identical AS only to derive kept reads, never feeding the EM) is
        #     removed: ~40% of the per-interval cost on dense loci. ---
        preset = get_m1_preset()
        aligners = [
            mappy.Aligner(seq=c.sequence, preset=preset) if c.sequence else None
            for c in cand_list
        ]
        all_ids = [rid for rid, _ in reads_iter]
        raw_all = np.full((len(all_ids), n_c), -np.inf, dtype=np.float64)
        for i, rid in enumerate(all_ids):
            seq = read_seqs[rid]
            for j, aln in enumerate(aligners):
                if aln is None:
                    continue
                best = None
                for h in aln.map(seq):
                    v = score_hit(h)
                    if v is None:
                        continue
                    if best is None or v > best:
                        best = v
                if best is not None:
                    raw_all[i, j] = best

        # Keep reads with >=1 candidate at AS>0 (matches the removed
        # mappy_multimap_responsibilities keep rule: row.sum()>0 over best>=0 AS).
        # Mirror its MAPPY_R1_MIN_AS env knob so the kept set stays identical when
        # that R1 tuning threshold is set.
        import os as _os
        _min_as = float(_os.environ.get("MAPPY_R1_MIN_AS", "0") or "0")
        keep_rows = ((raw_all > 0.0) & (raw_all >= _min_as)).any(axis=1)
        if not keep_rows.any():
            return []
        kept_read_ids = [rid for rid, k in zip(all_ids, keep_rows) if k]
        raw = raw_all[keep_rows]
        n_r = len(kept_read_ids)

        # --- Junction NLL on the per-read AS-tie cells only. ---
        nlls_by_read, ties_by_read, n_ties, n_refined, cover_by_read = self._tie_nll(
            kept_read_ids, read_seqs, cand_list, aligners, raw
        )

        # --- d_tx skeleton. Gate OFF: the pure soft NLL-graded skeleton (M1/AS
        #     picks the tie set; junction-NLL is the only graded distance; cells
        #     outside the tie stay MISSING). Gate ON: the diff-region coverage
        #     decision matrix with proportional redistribution of ambiguous reads
        #     (no read is ever dropped -> recall-safe). ---
        gate = bool(getattr(self.config, "m2_diff_cover_gate", False))
        d_tx = np.full((n_r, n_c), MISSING, dtype=np.float64)
        if not gate:
            for i in range(n_r):
                tie = ties_by_read.get(i)
                if not tie:
                    continue
                if len(tie) < 2:
                    d_tx[i, tie[0]] = 0.0
                    continue
                nlls = nlls_by_read.get(i)
                if not nlls:
                    # nothing scorable -> flat tie (EM 1/K split among ties)
                    for j in tie:
                        d_tx[i, j] = 0.0
                    continue
                lo = min(nlls.values())
                hi = max(nlls.values())
                miss = (hi - lo) + _M2_EM_PAD
                for j in tie:
                    d_tx[i, j] = (nlls[j] - lo) if j in nlls else miss
            small = (d_tx > 0) & (d_tx < 1e5)
            diffs = d_tx[small]
            sigma2 = max(float(np.median(diffs)) if diffs.size else 1.0, 1e-3)
            d_tx[small] = d_tx[small] / sigma2
        else:
            # Diff-region coverage gate. margin = runner-up NLL - best NLL.
            #
            # Pass A: classify each tied read and accumulate a per-candidate prior
            # (covered_vote) from ONLY the covered + distinguishing reads -- the
            # gold-standard evidence of the locus isoform ratio. Ambiguous reads
            # (no junction signal, or covered-but-indistinguishable that did not
            # straddle every diff span) are deferred to Pass B.
            margin_thr = float(getattr(self.config, "m2_diff_cover_margin", 0.0))
            covered_vote = np.zeros(n_c, dtype=np.float64)
            fuzzy: List[int] = []
            for i in range(n_r):
                tie = ties_by_read.get(i)
                if not tie:
                    continue
                if len(tie) < 2:
                    d_tx[i, tie[0]] = 0.0
                    continue
                nlls = nlls_by_read.get(i)
                if not nlls:
                    fuzzy.append(i)  # no junction signal -> defer (was: drop)
                    continue
                ordered = sorted(nlls.items(), key=lambda kv: kv[1])
                best_j = ordered[0][0]
                margin = (ordered[1][1] - ordered[0][1]) if len(ordered) > 1 else float("inf")
                # The redistribution PRIOR (covered_vote) must be fed only by reads
                # with EXPLICIT coverage: the per-read _tie_nll fallback leaves
                # cover_by_read empty (coverage uncomputable), so a missing entry
                # must NOT seed covered_vote (else the fallback path biases fuzzy
                # reads off unverified coverage).
                covered_for_vote = cover_by_read.get(i) is True
                if margin >= margin_thr:
                    # M2 distinguishes -> hard assign to the lowest-NLL candidate.
                    d_tx[i, best_j] = 0.0  # every other tie cell stays MISSING
                    if covered_for_vote:
                        # Only covered+distinguishing reads define the locus ratio.
                        covered_vote[best_j] += 1.0
                else:
                    # Indistinguishable (margin < thr), whether or not covered:
                    # defer to Pass B and redistribute by the covered-read ratio
                    # (flat 1/K only when the tie has no covered prior). The margin
                    # is thus the SOLE decider of hard-assign vs ratio-follow.
                    fuzzy.append(i)

            # Pass B: redistribute each ambiguous read across its tie in proportion
            # to the covered_vote ratio of those candidates (one-shot prior, not an
            # EM feedback loop). d_tx = -log(p) so softmax(-d_tx) reproduces p at
            # sigma=1; a zero-prior candidate stays MISSING. If NONE of the tie
            # candidates earned covered votes, fall back to the flat 1/K split
            # (== legacy behavior -> signal-dead loci are never starved).
            for i in fuzzy:
                tie = ties_by_read[i]
                masses = np.array([covered_vote[j] for j in tie], dtype=np.float64)
                s = masses.sum()
                if s > 0:
                    p = masses / s
                    for idx, j in enumerate(tie):
                        if p[idx] > 0:
                            d_tx[i, j] = -float(np.log(p[idx]))
                        # p==0 -> leave MISSING (no covered support for this cand)
                else:
                    for j in tie:
                        d_tx[i, j] = 0.0  # flat 1/K (recall-safe fallback)
            # ON-path d_tx values are already calibrated (0 / -log(p) / MISSING) for
            # sigma=1; no sigma2 normalization (that is the OFF skeleton's step).
        d_tx = d_tx.astype(np.float32)

        # --- M3 read×read junction-window DTW coherence (opt-in; DTW is costly).
        #     Each read anchored to its d_tx-argmin (best) candidate; σ3 (median
        #     non-zero DTW) calibrates it to O(1) so β mixes cleanly. ---
        if self.config.m3_coherence:
            from fin.scoring.m3_junction_coherence import build_m3_coherence

            winner_col = np.asarray(d_tx).argmin(axis=1).astype(np.int64)
            no_data = np.asarray(d_tx).min(axis=1) >= 1e5  # all-MISSING rows
            winner_col[no_data] = -1
            dist_read_to_read = build_m3_coherence(
                kept_read_ids, read_seqs, cand_list, winner_col,
                self.config.signal_path, pore=self.config.krill_pore,
                junction_k=self.config.m2_tiebreak_junction_k, flank=_M2_EM_FLANK,
                use_gpu=self.config.use_gpu,
            )
            nz = dist_read_to_read[dist_read_to_read > 0]
            sigma3 = max(float(np.median(nz)) if nz.size else 1.0, 1e-3)
            dist_read_to_read = (dist_read_to_read / sigma3).astype(np.float32)
            beta_use = self.config.em_beta
        else:
            dist_read_to_read = np.zeros((n_r, n_r), dtype=np.float32)
            beta_use = 0.0

        R, hard_assignments, _ = em_with_coherence(
            dist_read_to_tx=d_tx,
            dist_read_to_read=dist_read_to_read,
            sigma=1.0,
            beta=beta_use,
            max_iter=max_iter_em,
            tol=self.config.em_tol,
            verbose=False,
            use_gpu=self.config.use_gpu,
            abundance_feedback=self.config.abundance_feedback,
            abundance_length_norm=self.config.abundance_length_norm,
            eff_lengths=self._eff_lengths(cand_list),
        )

        if self.config.krill_tiebreak:
            R = krill_tiebreak(
                R=R, read_ids=kept_read_ids, read_seqs=read_seqs,
                candidates=cand_list, signal_path=self.config.signal_path,
                pore=self.config.krill_pore,
                ambig_threshold=self.config.tiebreak_ambig_threshold,
                use_gpu=self.config.use_gpu,
            )
            hard_assignments = R.argmax(axis=1)

        # --- Cluster-internal wobble recheck (post-EM, data-layer shadow kill).
        #     Cluster ALL multi-exon candidates by structure (same intron count +
        #     every junction within +-bp); within a cluster the highest-EM-abundance
        #     candidate anchors (often the true isoform, frequently a GTF
        #     passthrough), and a NOVEL sibling whose abundance/anchor < fraction is
        #     a wobble shadow -> drop. GTF/fusion never dropped (they only anchor).
        #     Pure abundance evidence: GTF participates in the abundance race but is
        #     NOT used as a correctness oracle, so this stays robust when the
        #     annotation is wrong (and self-disables with no GTF: a real novel
        #     anchors its own cluster). OFF -> no drops (byte-identical). ---
        drop_cols: set = set()
        # Strand-keyed {strand: Counter{(donor,acceptor): n_reads}} of directly-
        # observed read junctions; built lazily (reopens BAM) and shared by the
        # cluster-recheck GTF guard and the Lever-2 junction-support gate.
        observed_jct = None
        if getattr(self.config, "m2_cluster_recheck", False):
            from fin.scoring.m2_junction_nll import structural_wobble_clusters

            frac = self.config.m2_cluster_recheck_fraction
            if frac <= 0.0:
                frac = self.config.min_isoform_fraction
            if frac > 0.0:
                em_ab = np.asarray(R).sum(axis=0)
                bp = self.config.m2_cluster_recheck_bp
                cassette_bp = getattr(
                    self.config, "m2_cluster_recheck_cassette_max_exon_bp", 0
                )
                clusters = structural_wobble_clusters(
                    cand_list, bp, cassette_max_exon_bp=cassette_bp
                )
                from fin.scoring.m2_junction_nll import (
                    gtf_guard_needed,
                    wobble_shadow_drops,
                )

                gtf_drop = bool(getattr(
                    self.config, "m2_cluster_recheck_novel_displaces_gtf", True
                ))
                # Build the directly-observed read junctions for the GTF read-support
                # guard ONLY when a clustered low-abundance GTF sibling could actually
                # be dropped — this reopens the BAM, so skip it for no-GTF intervals
                # (and avoids touching the BAM for mocked/placeholder-path callers).
                if gtf_drop and gtf_guard_needed(cand_list, em_ab, clusters, frac):
                    observed_jct = self._observed_junctions(interval)
                drop_cols |= wobble_shadow_drops(
                    cand_list, em_ab, clusters, frac,
                    gtf_drop_enabled=gtf_drop,
                    observed_jct=observed_jct,
                    jct_tol=int(getattr(self.config, "m2_cluster_recheck_jct_tol", 2)),
                    gtf_min_jct_reads=int(getattr(
                        self.config, "m2_cluster_recheck_gtf_min_jct_reads", 1
                    )),
                )

        # --- Lever 2: per-junction read-support gate. Drop a NOVEL multi-exon
        #     candidate if ANY of its junctions is spliced by fewer than
        #     novel_junction_min_reads directly-observed reads (CIGAR introns,
        #     strand-keyed, matched within novel_junction_reads_tol bp) — a novel
        #     junction must be carried by >= N reads, not just 1. gtf/fusion/mono
        #     exempt. <=1 disables (byte-identical). Reuses observed_jct if the
        #     cluster-recheck already built it; otherwise builds it here. ---
        _njmr = int(getattr(self.config, "novel_junction_min_reads", 0))
        if _njmr > 1:
            from fin.scoring.m2_junction_nll import junction_support_drops
            if observed_jct is None:
                observed_jct = self._observed_junctions(interval)
            drop_cols |= junction_support_drops(
                cand_list, observed_jct, min_reads=_njmr,
                tol=int(getattr(self.config, "novel_junction_reads_tol", 2)),
            )

        # --- GUIDED junction-support gate (EXPERIMENT, default-off). Mirror of
        #     Lever 2 for GTF-passthrough candidates: drop a guided multi-exon
        #     candidate if ANY of its junctions lacks >= guided_junction_min_reads
        #     directly-observed reads (CIGAR introns within guided_junction_reads
        #     _tol bp). Extends the coordinate-EXACT read-support check — which
        #     Lever 2 applies only to novel candidates — to GTF junctions, so a
        #     jitter-corrupted annotation junction (coords > tol from the true
        #     site, thus zero exact support) is dropped instead of surviving the
        #     coordinate-inexact M2/M1 gate. Reuses observed_jct. 0 disables
        #     (byte-identical). NOT recall-safe (also drops low-coverage annotated
        #     junctions); default-off pending real-data validation. ---
        _gjmr = int(getattr(self.config, "guided_junction_min_reads", 0))
        if _gjmr > 0:
            from fin.scoring.m2_junction_nll import guided_junction_support_drops
            if observed_jct is None:
                observed_jct = self._observed_junctions(interval)
            drop_cols |= guided_junction_support_drops(
                cand_list, observed_jct, min_reads=_gjmr,
                tol=int(getattr(self.config, "guided_junction_reads_tol", 2)),
            )

        # --- M2/M1 read-support gate. A multi-exon candidate (GTF or novel;
        #     fusion/mono exempt) must earn >=1 read's support: either it is some
        #     read's M1 SOLE best-AS (tie set == just this candidate), or it is
        #     some read's M2-best (lowest junction-NLL in that read's tie). A
        #     jitter-corrupted GTF junction wins neither (loses the M2 contest to
        #     the true-junction candidate and never holds a read's sole AS), so it
        #     is dropped; a genuine isoform always earns one or the other. ---
        if getattr(self.config, "m2_support_gate", False):
            from fin.scoring.m2_junction_nll import support_gate_drops

            drop_cols |= support_gate_drops(
                cand_list, ties_by_read, nlls_by_read,
                tie_ok=bool(getattr(self.config, "m2_support_gate_tie", True)),
            )

        # --- Containment-cluster drop: drop a NOVEL candidate whose intron chain is a
        #     contiguous SUB-CHAIN of a longer candidate (a truncation / exon-skip shadow
        #     the same-intron-count wobble cluster never groups) when it is a low-support
        #     shadow by BOTH EM abundance AND supporting-read count. Runs AFTER all the
        #     structural/support gates so ``exclude=drop_cols`` covers every already-doomed
        #     parent — a shadow is never folded into a parent that is itself dropped. The
        #     read-support guard keeps most genuine low-abundance short/alt-TSS isoforms;
        #     gtf/fusion never dropped. DEFAULT-ON (--no-containment-cluster disables). ---
        if getattr(self.config, "containment_cluster", False):
            from fin.scoring.m2_junction_nll import containment_cluster_drops
            em_ab_c = np.asarray(R).sum(axis=0)
            read_counts = [len(getattr(c, "supporting_read_ids", ()) or ())
                           for c in cand_list]
            drop_cols |= containment_cluster_drops(
                cand_list, em_ab_c, read_counts,
                wobble_bp=int(getattr(self.config, "containment_cluster_wobble_bp", 6)),
                min_ab_ratio=float(getattr(
                    self.config, "containment_cluster_min_ab_ratio", 0.3)),
                min_read_ratio=float(getattr(
                    self.config, "containment_cluster_min_read_ratio", 0.3)),
                max_shadow_reads=int(getattr(
                    self.config, "containment_cluster_max_shadow_reads", 0)),
                exclude=set(drop_cols),
            )

        quant_results = quantify_transcripts(
            R, hard_assignments, cand_list, kept_read_ids
        )
        # Stamp max_R BEFORE any filtering: quant_results is column-aligned to R
        # here (enumerate index == R column). Filtering would break that alignment.
        for j, qr in enumerate(quant_results):
            qr.max_R = float(R[:, j].max()) if R.shape[0] > 0 else 0.0
        # --- Lever 1: containment / 5'-truncation collapse (post-EM mass-fold).
        #     Fold a NOVEL candidate whose intron chain is a pure 3' suffix of a
        #     longer candidate (a 5'-truncation shadow) INTO that parent: the
        #     shadow's hard reads (read-id union) + soft EM mass move to the parent,
        #     then the shadow joins drop_cols. Parents already in drop_cols
        #     (wobble/support) are EXCLUDED, so a fold never targets an
        #     ALREADY-dropped interval candidate; a folded parent CAN still be
        #     removed by a later _finalize_and_write filter (isoform-fraction /
        #     soft-mass / fulllen / polyA), which would lose the absorbed reads --
        #     an accepted limitation given this lever is default-off. Chained
        #     containment resolves to the terminal longest parent. OFF (default)
        #     -> no folds (byte-identical). NOT recall-safe by
        #     construction (a 3'-suffix also fits a real alt-TSS isoform); see
        #     containment_shadow_drops -> kept default-off pending real-data tuning.
        if getattr(self.config, "containment_collapse", False):
            from fin.scoring.m2_junction_nll import containment_shadow_drops

            em_ab_c = np.asarray(R).sum(axis=0)
            fold = containment_shadow_drops(
                cand_list, em_ab_c,
                tol_bp=int(self.config.containment_3p_tol_bp),
                min_ratio=float(self.config.containment_min_abundance_ratio),
                exclude=drop_cols,
            )
            if fold:
                qr_by_id = {q.candidate_id: q for q in quant_results}
                n_folded = 0
                for shadow_col, parent_col in fold.items():
                    sq = qr_by_id.get(cand_list[shadow_col].candidate_id)
                    pq = qr_by_id.get(cand_list[parent_col].candidate_id)
                    if sq is None or pq is None:
                        continue
                    # Hard reads are disjoint by construction (EM argmax mask), so
                    # the union count == sum, but union is robust either way. Soft
                    # mass is added so aggregate_across_intervals (single interval:
                    # unique/weight ratio == 1) reports parent + shadow abundance.
                    union = set(pq.assigned_read_ids) | set(sq.assigned_read_ids)
                    pq.assigned_read_ids = tuple(sorted(union))
                    pq.num_assigned_reads = len(union)
                    pq.abundance += sq.abundance
                    # Keep the parent's max EM responsibility consistent with the
                    # absorbed reads (max_R is a reporting/FP-analysis metric, not a
                    # filter input; confidence is left as the parent's own mean).
                    pq.max_R = max(pq.max_R, sq.max_R)
                    drop_cols.add(shadow_col)
                    n_folded += 1
                if n_folded:
                    logger.info(
                        "m2_em interval %s: containment folded %d "
                        "5'-truncation shadows", interval.region_string, n_folded,
                    )

        # --- Post-EM mono-exon resolution (mono_resolve_post_em). Re-resolve each
        #     SURVIVING single-exon candidate's reads against the SURVIVING multi-exon
        #     candidates by strict strand-aware exonic containment: a read wholly inside
        #     one exon of exactly one surviving multi folds into it (a degradation
        #     fragment); inside several -> the highest-EM-abundance one; uncovered reads
        #     stay on the mono candidate, which is dropped if it retains fewer than
        #     mono_resolve_min_reads. Runs AFTER every structural gate, so it can only
        #     target survivors (guard 1: never fold into a doomed candidate). Containment
        #     is strict + strand-aware with terminal slop (guard 2). Multi-cover assigns
        #     the whole hard read to the top-abundance candidate, not a fractional split
        #     (guard 3); the uncovered read-support floor gates genuine mono calls. ---
        if getattr(self.config, "mono_resolve_post_em", False):
            from fin.scoring.m2_junction_nll import mono_resolve_drops
            mono_drops = mono_resolve_drops(
                cand_list, quant_results, np.asarray(R).sum(axis=0),
                getattr(candidate_set, "read_spans", None) or {},
                exclude=drop_cols,
                slop_bp=int(getattr(self.config, "mono_resolve_slop_bp", 10)),
                min_reads=int(getattr(self.config, "mono_resolve_min_reads", 2)),
            )
            if mono_drops:
                drop_cols |= mono_drops
                logger.info(
                    "m2_em interval %s: mono-resolve dropped %d under-supported "
                    "single-exon candidates", interval.region_string, len(mono_drops),
                )
        if drop_cols:
            drop_ids = {cand_list[j].candidate_id for j in drop_cols}
            quant_results = [q for q in quant_results if q.candidate_id not in drop_ids]
            logger.info(
                "m2_em interval %s: post-EM gates dropped %d candidates "
                "(cluster-recheck / support / containment)",
                interval.region_string, len(drop_ids),
            )
        logger.info(
            "m2_em interval %s: %d reads -> %d candidates "
            "(ties=%d refined=%d m3=%s beta=%.2f)",
            interval.region_string, len(read_ids), len(cand_list),
            n_ties, n_refined, self.config.m3_coherence, beta_use,
        )
        return quant_results

    def _get_fusion_genome_aligner(self):
        """Lazily build and cache the genome-wide mappy aligner for fusion arms."""
        if not self._fusion_aligner_built:
            from fin.fusion import build_genome_aligner

            self._fusion_aligner_built = True
            if self.config.genome_fasta_path:
                self._fusion_genome_aligner = build_genome_aligner(
                    self.config.genome_fasta_path
                )
        return self._fusion_genome_aligner

    def _augment_with_fusion_candidates(
        self, candidate_set: CandidateSet, interval: GenomicInterval
    ) -> CandidateSet:
        """Run the read-driven fusion sub-pipeline and merge into candidate_set.

        Collects chimeric (soft-clip) reads in the interval, re-aligns their
        soft-clip arms to find partners (Stage F1), infers per-arm splice
        structure (F2), and stitches cross-breakpoint fusion candidates (F3).
        """
        from fin.candidates.canonical import parse_motifs
        from fin.fusion import detect_fusion_candidates
        from fin.io.io_bam import BamReader

        aligner = self._get_fusion_genome_aligner()
        if aligner is None:
            return candidate_set

        with BamReader(self.config.bam_path) as bam:
            read_dicts = bam.get_reads_in_region(interval.region_string)
        if not read_dicts:
            return candidate_set

        fusion_cands = detect_fusion_candidates(
            read_dicts,
            aligner,
            self._genome_fasta or {},
            gtf_reader=self._gtf_reader,
            motif_set=parse_motifs(self.config.canonical_motifs),
            max_internal_gap_bp=self.config.fusion_max_internal_gap_bp,
            max_dist=self.config.fusion_max_dist,
            search_bp=self.config.canonical_search_bp,
            max_chains_per_read=self.config.max_chains_per_read,
            min_support=self.config.fusion_min_support,
        )
        if not fusion_cands:
            return candidate_set

        return merge_fusion_candidates(candidate_set, fusion_cands)

    def _apply_canonical_gate(
        self, candidate_set: CandidateSet, chrom_seq: str
    ) -> None:
        """Stage B: drop NOVEL multi-exon candidates with non-canonical junctions.

        Mutates ``candidate_set.candidates`` in place, keeping a candidate iff:
          - it is GTF-passthrough (source="gtf") or fusion (source="fusion") — EXEMPT;
          - it is single-exon (empty intron chain) — trivially canonical;
          - every internal junction's (donor,acceptor) motif is in canonical_motifs.

        ``candidate_set.read_ids`` is left untouched so reads off a dropped
        candidate are simply re-competed in quantification (they will land on a
        surviving candidate by argmin).
        """
        if not chrom_seq:
            # No genome sequence → cannot evaluate motifs; skip the gate rather
            # than drop everything.
            return
        motif_set = parse_motifs(self.config.canonical_motifs)
        kept: List = []
        dropped = 0
        for c in candidate_set.candidates:
            if c.source in ("gtf", "fusion"):
                kept.append(c)
                continue
            introns = c.intron_chain.introns
            if not introns:
                kept.append(c)  # mono: no internal junction
                continue
            if chain_all_canonical(introns, chrom_seq, c.strand, motif_set):
                kept.append(c)
            else:
                dropped += 1
        if dropped:
            logger.info(
                "Canonical gate: dropped %d/%d novel candidates (interval %s)",
                dropped,
                len(candidate_set.candidates),
                candidate_set.interval.region_string,
            )
        candidate_set.candidates = kept

    def _apply_junction_dominance_gate(
        self, candidate_set: CandidateSet, interval: GenomicInterval
    ) -> None:
        """Pre-EM junction-first gate: drop NOVEL candidates with a weak or
        non-dominant junction (mutates candidate_set.candidates in place).

        Builds the directly-observed read junctions (reopens BAM) and drops a
        novel multi-exon candidate if any junction has < min_reads reads or is
        dominated by a stronger DIFFERENT junction within the window. read_ids is
        left untouched (reads off a dropped candidate re-compete in quant).
        """
        from fin.scoring.m2_junction_nll import dominant_junction_drops

        observed = self._observed_junctions(interval)
        drops = dominant_junction_drops(
            candidate_set.candidates, observed,
            min_reads=int(self.config.junction_dominance_min_reads),
            window=int(self.config.junction_dominance_window_bp),
            tol=int(getattr(self.config, "junction_dominance_tol_bp", 2)),
        )
        if drops:
            kept = [c for i, c in enumerate(candidate_set.candidates)
                    if i not in drops]
            logger.info(
                "Junction-dominance gate: dropped %d/%d novel candidates "
                "(interval %s)", len(drops), len(candidate_set.candidates),
                interval.region_string,
            )
            candidate_set.candidates = kept

    def cleanup(self):
        """Close file handles."""
        if self._gtf_reader:
            self._gtf_reader.close()
        if self._signal_reader:
            self._signal_reader.close()

    def _load_genome_fasta(self, path: str) -> Dict[str, str]:
        """Load genome FASTA into a dict of chrom -> sequence."""
        from fin.io.io_fasta import FASTAReader

        seqs = {}
        with FASTAReader(path) as reader:
            for record in reader.iterate_records():
                seqs[record.id] = record.sequence
        return seqs

    def _open_signal_reader(self):
        """Open the appropriate signal reader."""
        if self.config.signal_format == "pod5":
            from fin.io.io_pod5 import Pod5Reader

            reader = Pod5Reader(self.config.signal_path)
            reader.open()
            return reader
        else:
            from fin.io.io_slow5 import Slow5Reader

            reader = Slow5Reader(self.config.signal_path)
            reader.open()
            return reader
