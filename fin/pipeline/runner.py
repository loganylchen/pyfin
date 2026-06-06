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


class PipelineRunner:
    """Orchestrates the full pyfin pipeline."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self._gtf_reader = None
        self._genome_fasta = None
        self._signal_reader = None

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

        # Process each interval
        all_quant_results: List[List[QuantResult]] = []
        for i, interval in enumerate(intervals):
            logger.info("Processing interval %d/%d: %s", i + 1, len(intervals), interval.region_string)
            quant = self.process_interval(interval)
            if quant:
                all_quant_results.append(quant)

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

        # A4: post-EM abundance filter. GTF-sourced and fusion transcripts are
        # exempt so users can still see zero-abundance annotated entries
        # (useful for debugging coverage gaps) and fusion calls are not
        # silently dropped by NOVEL-targeted filters. Only NOVEL candidates
        # are dropped here.
        if _score_filter_on and self.config.min_abundance > 0.0:
            before = len(aggregated)
            aggregated = {
                cid: qr
                for cid, qr in aggregated.items()
                if qr.source in ("gtf", "fusion")
                or qr.abundance >= self.config.min_abundance
            }
            dropped = before - len(aggregated)
            if dropped:
                logger.info(
                    "Dropped %d novel transcripts with abundance < %.3f",
                    dropped,
                    self.config.min_abundance,
                )

        # Phase A Tier-2: max_R filter. EM responsibility is a strong FP
        # discriminator (Cohen's d=1.34 vs combined_score's 0.70) -- but that
        # statistic was measured on the EM quant modes. NOTE: under the
        # production default quant_mode='argmax', max_R is the max single-read
        # *mappy* assignment weight (set in _quant_argmax_keep), NOT an EM
        # responsibility, so the d=1.34 calibration does not transfer; only
        # m1_em / m2_em populate max_R with a true EM responsibility.
        # GTF-sourced and fusion transcripts are exempt to preserve annotation
        # visibility and avoid silently dropping fusion calls.
        if _score_filter_on and self.config.min_max_r > 0.0:
            before = len(aggregated)
            aggregated = {
                cid: qr
                for cid, qr in aggregated.items()
                if qr.source in ("gtf", "fusion")
                or qr.max_R >= self.config.min_max_r
            }
            dropped = before - len(aggregated)
            if dropped:
                logger.info(
                    "Dropped %d novel transcripts with max_R < %.3f",
                    dropped,
                    self.config.min_max_r,
                )

        # Step 3: novel combined_score filter (FP-by-scoring). combined_score is
        # the geometric mean of coherence and discrimination (Cohen's d=0.70).
        # The F1-optimal threshold from profile_fp.md is 0.288 for with-GTF and
        # 0.428 for no-GTF runs. GTF-sourced and fusion transcripts are exempt.
        #
        # WARNING: combined_score is ONLY populated by the composite scorer,
        # which the assembly pipeline (this runner) never invokes -- it is wired
        # solely into the `quantify` subcommand. In every assembly quant_mode
        # (argmax / m1_em / m2_em) combined_score stays at its 0.0 default, so a
        # threshold > 0 here would drop EVERY novel candidate. Guard against that
        # foot-gun: skip the filter (with a loud warning) unless some novel
        # candidate actually carries a non-zero combined_score.
        if _score_filter_on and self.config.min_novel_combined_score > 0.0:
            if not any(
                qr.source not in ("gtf", "fusion") and qr.combined_score > 0.0
                for qr in aggregated.values()
            ):
                logger.warning(
                    "min_novel_combined_score=%.3f requested but no novel "
                    "candidate has a non-zero combined_score (the composite "
                    "scorer is not run in assembly mode); skipping this filter "
                    "to avoid dropping all novel transcripts.",
                    self.config.min_novel_combined_score,
                )
            else:
                before = len(aggregated)
                aggregated = {
                    cid: qr
                    for cid, qr in aggregated.items()
                    if qr.source in ("gtf", "fusion")
                    or qr.combined_score >= self.config.min_novel_combined_score
                }
                dropped = before - len(aggregated)
                if dropped:
                    logger.info(
                        "Dropped %d novel transcripts with combined_score < %.3f",
                        dropped,
                        self.config.min_novel_combined_score,
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

        # Full-length end-coherence filter (FLAIR/TALON-style full-length read
        # support). Drops NOVEL multi-exon (>=2 intron) transcripts whose
        # fraction of full-length assigned reads (read genomic 5' AND 3' both
        # within fulllen_window_bp of the candidate's ends) is below
        # min_fulllen_fraction. The fulllen_frac METRIC itself is signal-free
        # (BAM primary-alignment spans only), but it is computed per-interval
        # over the production argmax assignment population (which may use the
        # M2-tiebreak signal); candidates with fulllen_frac < 0 (unreachable
        # or never scored, e.g. the legacy EM path) are EXEMPT and never
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
        # 5' end. UNLIKE the other filters this also gates GTF-sourced
        # candidates (fusion stays exempt). Requires signal; no-ops when absent.
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
            canonical_search_bp=self.config.canonical_search_bp,
            max_chains_per_read=self.config.max_chains_per_read,
            canonical_motifs=self.config.canonical_motifs,
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
        quant_mode = getattr(self.config, "quant_mode", "argmax")
        if quant_mode == "argmax":
            results = self._quant_argmax_keep(candidate_set, read_ids, interval)
        elif quant_mode == "m1_em":
            results = self._quant_m1_em(candidate_set, read_ids, interval)
        elif quant_mode == "m2_em":
            results = self._quant_m2_em(candidate_set, read_ids, interval)
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
        if m2_on:
            try:
                import krill

                from fin.scoring.krill_aligner import make_krill_aligner

                m2_krill, m2_gpu = make_krill_aligner(
                    krill, self.config.krill_pore, self.config.use_gpu,
                    hmm_confidence=False,
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
                    return_scored=True, use_gpu=m2_gpu,
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
            )
            qr.max_R = max_weight[j]
            quant_results.append(qr)

        logger.info(
            "R1 M1-keep(split) interval %s: %d reads -> %d candidates "
            "(M2 overrides=%d)",
            interval.region_string, len(read_ids), len(cand_list), n_m2_override,
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
        """Drop candidates failing the polyA + 5'-proximity gate (gates GTF too).

        Runs a krill whole-read polyA pass over all reads, then removes any
        novel/gtf candidate with fewer than ``min_polya5p_reads`` reads that both
        have a confident polyA tail and map 5'-flush to the candidate. Fusion
        candidates are exempt. No-ops if krill returns nothing.
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

    def _quant_m2_em(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> Optional[List[QuantResult]]:
        """quant_mode='m2_em': EM seeded by the M2 krill junction distance plus
        the M3 read×read krill DTW coherence (β=em_beta when m4_source != 'none').
        All signal scoring is in-memory krill — no f5c CLI."""
        from fin.ablation.mappy_argmax import mappy_multimap_responsibilities
        from fin.scoring.krill_tiebreak import _build_m2_krill
        from fin.scoring.m3_junction_coherence import build_m3_coherence

        cand_list = list(candidate_set.candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        reads_iter = [(rid, read_sequences.get(rid, "")) for rid in read_ids]
        reads_iter = [(rid, seq) for rid, seq in reads_iter if seq]
        R_mm, kept_read_ids = mappy_multimap_responsibilities(reads_iter, cand_list)
        if R_mm.size == 0:
            return []

        read_seqs = {rid: seq for rid, seq in reads_iter}
        n_reads_em = len(kept_read_ids)
        max_iter_em = (
            self.config.em_max_iter_override
            if self.config.em_max_iter_override is not None
            else self.config.em_max_iter
        )

        # M2: per-(read, candidate) krill junction distance (lower = better).
        dist_read_to_tx = _build_m2_krill(
            kept_read_ids, read_seqs, cand_list,
            self.config.signal_path, self.config.krill_pore, as_distance=True,
            use_gpu=self.config.use_gpu,
        )

        # M3: read×read junction-window DTW coherence, each read anchored to its
        # M2-best candidate. m4_source='none' disables coherence (β=0).
        m4_src = getattr(self.config, "m4_source", "diff_region")
        if m4_src == "none":
            dist_read_to_read = np.zeros((n_reads_em, n_reads_em), dtype=np.float32)
            beta_use = 0.0
        else:
            if m4_src == "whole_read":
                logger.warning(
                    "m2_em: m4_source='whole_read' is not supported on krill; "
                    "using junction-window coherence (diff_region)."
                )
            winner_col = np.asarray(dist_read_to_tx).argmin(axis=1).astype(np.int64)
            # Reads with no krill signal have an all-default row; mark uncoupled.
            no_data = dist_read_to_tx.min(axis=1) >= 0.999
            winner_col[no_data] = -1
            dist_read_to_read = build_m3_coherence(
                kept_read_ids, read_seqs, cand_list, winner_col,
                self.config.signal_path, pore=self.config.krill_pore,
                junction_k=self.config.m2_tiebreak_junction_k,
                use_gpu=self.config.use_gpu,
            )
            beta_use = self.config.em_beta

        sigma_use = self._adaptive_sigma(dist_read_to_tx)
        R, hard_assignments, _ = em_with_coherence(
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_read=dist_read_to_read,
            sigma=sigma_use,
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

        quant_results = quantify_transcripts(
            R, hard_assignments, cand_list, kept_read_ids
        )
        for j, qr in enumerate(quant_results):
            qr.max_R = float(R[:, j].max()) if R.shape[0] > 0 else 0.0
        logger.info(
            "m2_em interval %s: %d reads -> %d candidates (m4=%s, beta=%.2f)",
            interval.region_string, len(read_ids), len(cand_list),
            m4_src, beta_use,
        )
        return quant_results

    def _augment_with_fusion_candidates(
        self, candidate_set: CandidateSet, interval: GenomicInterval
    ) -> CandidateSet:
        """Detect fusion breakpoints for the interval and merge into candidate_set."""
        from fin.fusion import (
            build_fusion_candidates,
            cluster_breakpoints,
            parse_sa_tags,
        )

        region = interval.region_string
        raw_bps = parse_sa_tags(
            self.config.bam_path, region=region, min_mapq=10
        )
        clusters = cluster_breakpoints(
            raw_bps,
            max_dist=self.config.fusion_max_dist,
            min_support=self.config.fusion_min_support,
        )
        if not clusters or build_fusion_candidates is None:
            return candidate_set

        fusion_cands = build_fusion_candidates(
            clusters,
            self._genome_fasta or {},
            flank_bp=self.config.fusion_flank_bp,
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
