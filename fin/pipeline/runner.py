"""Pipeline orchestrator: interval -> candidates -> scoring -> EM -> quantification."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pysam

from fin.analysis.abundance_refit import (
    ResponsibilityLedger,
    annotate_selection_metadata,
    build_responsibility_ledger,
    refit_survivor_abundance,
    write_abundance_refit_diagnostics,
)
from fin.analysis.assignments import em_with_coherence
from fin.analysis.quantification import (
    QuantResult,
    _exons_from_candidate,
    aggregate_across_intervals,
    compute_fulllen_frac,
    polya5p_drops,
    quantify_transcripts,
)
from fin.candidates.dataclasses import CandidateSet
from fin.candidates.discovery import discover_candidates, merge_fusion_candidates
from fin.io.interval_manager import GenomicInterval, generate_isolated_intervals
from fin.pipeline.assignment import Assigner
from fin.pipeline.config import PipelineConfig
from fin.pipeline.finalize import finalize_outputs, write_unfiltered_diagnostic
from fin.pipeline.junction_snap import apply_junction_snap
from fin.pipeline.selection import (
    canonical_gate_select,
    junction_dominance_select,
    select_global,
    select_m2_interval,
)
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


IntervalOutput = Tuple[List[QuantResult], Optional[ResponsibilityLedger]]


def _order_interval_outputs(
    outputs: List[IntervalOutput], refit_effective: bool
) -> List[IntervalOutput]:
    """Canonicalize aggregation order only for the refit-enabled path."""
    if not refit_effective:
        return outputs
    if any(ledger is None for _, ledger in outputs):
        raise RuntimeError(
            "post-selection refit requires a responsibility ledger for every "
            "quantified interval"
        )
    # Genomic lexical order is sufficient; interval_key includes strand and is
    # unique for the isolated interval construction.
    return sorted(outputs, key=lambda pair: pair[1].interval_key)


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
        # m2_em assignment layer (mappy AS -> tie_nll -> d_tx -> EM -> quantify).
        self._assigner = Assigner(config)

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
        interval_outputs: List[IntervalOutput] = []
        if self.config.threads <= 1:
            for i, interval in enumerate(intervals):
                logger.info("Processing interval %d/%d: %s", i + 1, len(intervals), interval.region_string)
                quant = self.process_interval(interval)
                if quant is not None:
                    interval_outputs.append(quant)
        else:
            from fin.pipeline.parallel import run_parallel

            log_level = logging.getLevelName(logging.getLogger("fin").getEffectiveLevel())
            interval_outputs = run_parallel(
                self.config,
                intervals,
                self.config.threads,
                self.config.gpu_workers,
                log_level,
            )

        refit_effective = bool(getattr(
            self.config, "post_selection_refit_effective", False
        ))
        interval_outputs = _order_interval_outputs(
            interval_outputs, refit_effective
        )
        all_quant_results: List[List[QuantResult]] = []
        responsibility_ledgers: List[ResponsibilityLedger] = []
        for results, ledger in interval_outputs:
            if results:
                all_quant_results.append(results)
            if ledger is not None:
                responsibility_ledgers.append(ledger)

        # Aggregate across intervals. Responsibility ledgers remain separate so
        # the historical aggregation/selection evidence is not changed by refit.
        aggregated = aggregate_across_intervals(all_quant_results)
        return self._finalize_and_write(
            aggregated,
            output_gtf=self.config.output_gtf,
            output_tsv=self.config.output_tsv,
            responsibility_ledgers=responsibility_ledgers,
        )

    def _finalize_and_write(
        self,
        aggregated: Dict[str, QuantResult],
        output_gtf: Optional[str],
        output_tsv: Optional[str],
        responsibility_ledgers: Optional[List[ResponsibilityLedger]] = None,
    ) -> Dict[str, QuantResult]:
        """Finalize: diagnostic dump -> GLOBAL selection -> gene-id resolution + writers.

        Thin orchestrator wiring three independent layers: the pre-filter diagnostic and the
        output writers live in fin.pipeline.finalize; the GLOBAL selection cascade lives in
        fin.pipeline.selection. Kept as a method so tests can call this seam directly."""
        write_unfiltered_diagnostic(
            self.config, aggregated, output_tsv, self._gtf_reader
        )
        aggregated, _global_outcomes = select_global(
            self.config, aggregated, self._apply_polya5p_filter,
        )
        want_evidence = bool(getattr(self.config, "candidate_evidence", False))
        want_ranking = (
            getattr(self.config, "ranking_mode", "off") == "filter"
        )
        if want_evidence or want_ranking:
            evidence_rows = self._compute_candidate_evidence_rows(
                aggregated, strict=want_ranking,
            )
            if want_evidence and evidence_rows is not None:
                from fin.analysis.candidate_evidence import (
                    write_candidate_evidence,
                )
                write_candidate_evidence(
                    Path(self.config.work_dir) / "candidate_evidence.tsv",
                    evidence_rows,
                )
            if want_ranking:
                if evidence_rows is None:
                    # An explicitly requested filter must never silently
                    # degrade into a no-op baseline run.
                    raise RuntimeError(
                        "ranking_mode=filter requires complete candidate "
                        "evidence, but the evidence computation failed"
                    )
                from fin.analysis.candidate_ranking import ranking_filter

                aggregated, _scores, dropped = ranking_filter(
                    aggregated,
                    evidence_rows,
                    threshold=getattr(
                        self.config, "ranking_threshold", None
                    ),
                )
                if dropped:
                    logger.info(
                        "Ranking filter removed %d candidates before "
                        "finalization", len(dropped),
                    )
        refit_effective = bool(getattr(
            self.config, "post_selection_refit_effective", False
        ))
        endpoint_routes: dict = {}
        endpoint_primary: dict = {}
        if refit_effective and getattr(self.config, "endpoint_refine", False):
            from fin.analysis.endpoint_refine import (
                apply_endpoint_splits,
                plan_endpoint_splits,
            )

            bam_ev = getattr(self, "_ranking_bam_evidence", None)
            if bam_ev is None:
                from fin.analysis.candidate_evidence import (
                    collect_ranking_bam_evidence,
                )
                bam_ev = collect_ranking_bam_evidence(self.config.bam_path)
                self._ranking_bam_evidence = bam_ev
            if not bam_ev.complete:
                raise RuntimeError(
                    "endpoint_refine requires a complete BAM end-evidence "
                    f"scan ({bam_ev.error})"
                )
            polya_read_ids = self._endpoint_polya_read_ids()
            tss_mode = getattr(self.config, "tss_evidence_mode", "off")
            background = 0.0
            if tss_mode != "off":
                background = self._estimate_degradation_background(
                    aggregated, bam_ev.read_ends,
                    bin_bp=int(getattr(self.config, "endpoint_window_bp", 25)),
                )
            plan = plan_endpoint_splits(
                aggregated,
                bam_ev.read_ends,
                tss_evidence_mode=tss_mode,
                background_hazard=background,
                window_bp=int(getattr(self.config, "endpoint_window_bp", 25)),
                min_end_reads=int(getattr(self.config, "endpoint_min_reads", 3)),
                min_pair_frac=float(
                    getattr(self.config, "endpoint_min_pair_frac", 0.15)
                ),
                max_splits=int(getattr(self.config, "endpoint_max_splits", 2)),
                polya_read_ids=polya_read_ids,
            )
            import json as _json
            diag_path = Path(self.config.work_dir) / "endpoint_refine.json"
            tmp = diag_path.with_suffix(".json.tmp")
            tmp.write_text(_json.dumps({
                "schema_version": 1,
                "polya_available": polya_read_ids is not None,
                "n_polya_reads": (
                    len(polya_read_ids) if polya_read_ids is not None else None
                ),
                "tss_evidence_mode": tss_mode,
                "degradation_background_hazard": background,
                "models_split": len(plan.replacements),
                "endpoint_states": plan.n_splits,
                "details": plan.details,
                "tss_verdicts": plan.tss_verdicts,
            }, indent=2, sort_keys=True) + "\n")
            tmp.replace(diag_path)
            if plan.replacements:
                aggregated = apply_endpoint_splits(aggregated, plan)
                endpoint_routes = plan.read_routes
                endpoint_primary = plan.primary
                logger.info(
                    "EndpointRefine split %d models into %d endpoint states "
                    "(polyA %s)",
                    len(plan.replacements), plan.n_splits,
                    "available" if polya_read_ids is not None else "unavailable",
                )
        if refit_effective:
            aggregated, snap_redirects = apply_junction_snap(
                self.config,
                aggregated,
                getattr(self, "_genome_fasta", None),
                return_redirects=True,
            )
            ledgers = responsibility_ledgers or []
            if aggregated and not ledgers:
                raise RuntimeError(
                    "post-selection refit is enabled but no responsibility "
                    "ledgers reached global finalization"
                )
            aggregated, diagnostics = refit_survivor_abundance(
                aggregated, ledgers, snap_redirects=snap_redirects,
                split_routes=endpoint_routes, split_primary=endpoint_primary,
            )
            write_abundance_refit_diagnostics(
                Path(self.config.work_dir) / "abundance_refit.json",
                diagnostics,
            )
        else:
            aggregated = apply_junction_snap(
                self.config, aggregated, getattr(self, "_genome_fasta", None)
            )
            if (
                getattr(self.config, "post_selection_refit", False)
                and getattr(
                    self.config, "post_selection_refit_disable_reason", None
                )
            ):
                write_abundance_refit_diagnostics(
                    Path(self.config.work_dir) / "abundance_refit.json",
                    {
                        "schema_version": 1,
                        "mode": "post_selection_survivor_renormalization",
                        "effective": False,
                        "disable_reason": self.config.post_selection_refit_disable_reason,
                        "structural_identity": "not_applicable",
                    },
                )
        return finalize_outputs(
            self.config, aggregated, output_gtf, output_tsv, self._gtf_reader,
        )

    def _estimate_degradation_background(
        self, aggregated: Dict[str, QuantResult],
        read_ends: Mapping[str, tuple], *, bin_bp: int = 25,
    ) -> float:
        """Pooled 5'-termination hazard from models with no contained sibling.

        The null must describe degradation, so it is estimated only from
        models that no other surviving model is nested inside: a locus that
        already holds a shorter sibling would contaminate the background with
        that sibling's genuine TSS.
        """
        from fin.analysis.tss_evidence import (
            build_hazard_profile,
            genomic_to_offset,
            pooled_background_hazard,
            spliced_length,
        )

        def chain(ex):
            ex = sorted(ex)
            return tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))

        by_key: Dict[tuple, list] = {}
        for cid, qr in aggregated.items():
            by_key.setdefault((qr.chrom, qr.strand), []).append(cid)
        has_nested = set()
        for _key, ids in by_key.items():
            for a in ids:
                ca = chain(aggregated[a].exons)
                for b in ids:
                    if a == b:
                        continue
                    qb = aggregated[b]
                    qa = aggregated[a]
                    if qb.start <= qa.start and qa.end <= qb.end and \
                            (qa.end - qa.start) < (qb.end - qb.start):
                        has_nested.add(b)
        profiles = []
        for cid, qr in aggregated.items():
            if cid in has_nested or len(qr.exons) < 2:
                continue
            offs = []
            plus = qr.strand == "+"
            for rid in getattr(qr, "assigned_read_ids", ()) or ():
                sp = read_ends.get(rid)
                if sp is None:
                    continue
                r5 = sp[0] if plus else sp[1]
                o = genomic_to_offset(r5, qr.exons, qr.strand)
                if o is not None:
                    offs.append(o)
            if len(offs) >= 20:
                profiles.append(build_hazard_profile(
                    offs, spliced_length(qr.exons), bin_bp=bin_bp))
        bg = pooled_background_hazard(profiles)
        logger.info(
            "TSS evidence: degradation background %.6f per %dbp bin from %d "
            "nesting-free models", bg, bin_bp, len(profiles),
        )
        return bg

    def _endpoint_polya_read_ids(self) -> Optional[set]:
        """Poly(A)-confident read IDs for EndpointRefine TES support.

        Reuses the same krill whole-read polyA pass and confidence predicate
        as the polyA+5' filter. Returns None (documented degraded mode: TES
        modes rely on end sharpness alone) when no signal file is configured
        or krill yields nothing; EndpointRefine logs which mode ran.
        """
        if not getattr(self.config, "signal_path", None):
            return None
        try:
            from fin.analysis.quantification import polya_read_passes
            from fin.scoring.polya import compute_polya

            read_seqs, _read_ends = self._fetch_read_seqs_and_ends()
            polya_map = compute_polya(
                read_seqs,
                self.config.signal_path,
                pore=self.config.krill_pore,
                use_gpu=self.config.use_gpu,
            )
            if not polya_map:
                return None
            min_len = float(getattr(self.config, "min_polya_length", 15.0))
            return {
                rid for rid in polya_map
                if polya_read_passes(rid, polya_map, min_len)
            }
        except Exception:
            logger.exception(
                "endpoint polyA pass failed; TES support degrades to "
                "end-sharpness only"
            )
            return None

    def _compute_candidate_evidence_rows(
        self, aggregated: Dict[str, QuantResult], *, strict: bool = False
    ) -> Optional[list]:
        """Observable-feature rows for the post-selection survivor set.

        Computed BEFORE junction snapping and the abundance refit, so rows
        show exactly the feature view the selection-stage scorer sees. One
        whole-BAM pass supplies junction support and read ends together.

        ``strict=False`` (pure audit): failures and incomplete BAM scans are
        logged and tolerated (partial evidence may still be written for
        inspection). ``strict=True`` (the ranking filter): an incomplete BAM
        scan or any failure returns None so the caller refuses to filter -
        undercounted junction support must never depress ranking scores.
        """
        try:
            from fin.analysis.candidate_evidence import (
                collect_ranking_bam_evidence,
                compute_candidate_evidence,
            )
            from fin.candidates.canonical import parse_motifs

            bam_ev = collect_ranking_bam_evidence(self.config.bam_path)
            self._ranking_bam_evidence = bam_ev
            if strict and not bam_ev.complete:
                logger.error(
                    "ranking evidence BAM scan incomplete (%s); refusing to "
                    "rank on partial junction support", bam_ev.error,
                )
                return None
            return compute_candidate_evidence(
                aggregated,
                junction_support=bam_ev.junction_support,
                read_ends=bam_ev.read_ends,
                genome=self._genome_fasta,
                canonical_motifs=parse_motifs(
                    getattr(self.config, "canonical_motifs", ()) or ()
                ),
                end_window_bp=int(getattr(self.config, "fulllen_window_bp", 25)),
            )
        except Exception:
            logger.exception("candidate-evidence computation failed")
            return None

    def process_interval(
        self, interval: GenomicInterval
    ) -> Optional[Tuple[List[QuantResult], Optional[ResponsibilityLedger]]]:
        """Discover candidates and quantify a single interval via ``quant_mode``.

        Returns:
            ``(results, ledger)`` for a quantified interval. The ledger is
            populated only when post-selection refit is effectively enabled.
            Returns ``None`` when the interval has no candidates or no reads.
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
                # Defer the generation-time mono fold ONLY when post-EM resolution will
                # actually run (mono_resolve_post_em is m2_em-only). For argmax/m1_em/cluster
                # there is no post-EM mono step, so keep the generation fold to suppress FPs.
                and not (getattr(self.config, "mono_resolve_post_em", False)
                         and getattr(self.config, "quant_mode", "m2_em") == "m2_em")),
            chain_cluster_fold_span_guard=getattr(
                self.config, "chain_cluster_fold_span_guard", False),
            clustering=getattr(self.config, "clustering", "read_chains"),
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
        ledger: Optional[ResponsibilityLedger] = None
        if quant_mode == "argmax":
            results = self._quant_argmax_keep(candidate_set, read_ids, interval)
        elif quant_mode == "m1_em":
            results = self._quant_m1_em(candidate_set, read_ids, interval)
        elif quant_mode == "m2_em":
            results, ledger = self._quant_m2_em(
                candidate_set, read_ids, interval
            )
        elif quant_mode == "cluster":
            results = self._quant_cluster(candidate_set, read_ids, interval)
        else:
            raise ValueError(f"unknown quant_mode: {quant_mode!r}")
        if not results and ledger is None:
            return None

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
        return results, ledger

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
                family_id=cand.family_id,
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
                family_id=cand.family_id,
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
                    # Honor the configured tight-window metric so the
                    # experimental sqrt_count variant behaves consistently in
                    # cluster mode; anything else keeps the validated summed
                    # contrast this resolver was built around.
                    _metric = (
                        self.config.m2_metric
                        if getattr(self.config, "m2_metric", "summed_llr")
                        in ("summed_llr", "sqrt_count_mean_llr")
                        else "summed_llr"
                    )
                    out = m2_resolve_tie(
                        rid, read_sequences[rid], tied_cands,
                        self.config.signal_path, pore=self.config.krill_pore,
                        krill_aligner=krill_aligner, mappy_aligners=tied_aligners,
                        use_gpu=eff_gpu, num_thread=krill_threads,
                        metric=_metric,
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

                        short_chain = IntronChain(introns=schain)
                        seq = _build_spliced_sequence(
                            genome_seq, start, end, short_chain, strand,
                        )
                        if not seq:
                            continue

                        new_id = _generate_novel_id(
                            cand.chrom, strand, start, end, short_chain
                        )
                        new_cand = TranscriptCandidate(
                            candidate_id=new_id,
                            intron_chain=short_chain,
                            three_prime_pos=(end if strand == "+" else start),
                            sequence=seq,
                            source="novel",
                            supporting_read_ids=set(sreads),
                            chrom=cand.chrom,
                            strand=strand,
                            start=start,
                            end=end,
                            family_id=cand.family_id,
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
                family_id=cand.family_id,
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
        dist_read_to_tx, _, _ = build_em_matrices("m1", m1, m2_dummy)
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
        """Per-read best-AS tie junction-NLL. Delegates to the assignment layer
        (fin.pipeline.assignment.tie_nll); kept as a thin method so tests can patch
        this mock seam. Returns (nlls_by_read, ties_by_read, n_ties, n_refined,
        cover_by_read)."""
        from fin.pipeline.assignment import tie_nll
        return tie_nll(
            self.config, kept_read_ids, read_seqs, cand_list, aligners, raw)

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
    ) -> Tuple[List[QuantResult], Optional[ResponsibilityLedger]]:
        """quant_mode='m2_em' (production default): pure tie-break junction-NLL EM.

        M1/AS is used ONLY as a hard selector — it picks each read's best-AS tie
        set and masks mappability — and the per-event junction NLL is the SOLE
        graded distance over that tie set (m2_resolve_tie semantics, NOT the dense
        read×candidate matrix, which collapses). Production uses no read×read
        coherence term. All signal scoring is in-memory krill; there is no f5c
        CLI path.
        """
        qo = self._assigner.assign(
            candidate_set, read_ids,
            tie_nll_fn=self._tie_nll, eff_lengths_fn=self._eff_lengths,
        )
        if qo is None:
            return [], None
        ledger = None
        if getattr(self.config, "post_selection_refit_effective", False):
            ledger = build_responsibility_ledger(
                qo.R,
                qo.kept_read_ids,
                [candidate.candidate_id for candidate in qo.cand_list],
                input_read_ids=read_ids,
                interval_key=(
                    f"{interval.chrom}:{interval.start}-{interval.end}:"
                    f"{interval.strand}"
                ),
            )
        quant_results, outcomes = select_m2_interval(
            self.config, qo, candidate_set, interval, self._observed_junctions,
        )
        if ledger is not None:
            annotate_selection_metadata(
                ledger,
                candidates=qo.cand_list,
                read_ids=qo.kept_read_ids,
                hard_assignments=qo.hard_assignments,
                surviving_results=quant_results,
                outcomes=outcomes,
                mono_resolve_applied=bool(getattr(
                    self.config, "mono_resolve_post_em", False
                )),
            )
        logger.info(
            "m2_em interval %s: %d reads -> %d candidates "
            "(ties=%d refined=%d)",
            interval.region_string, len(read_ids), len(qo.cand_list),
            qo.n_ties, qo.n_refined,
        )
        return quant_results, ledger

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
        """Stage B canonical-junction gate. Delegates to the selection layer
        (fin.pipeline.selection.canonical_gate_select); kept as a method so tests can
        call this seam directly. Mutates ``candidate_set.candidates`` in place."""
        canonical_gate_select(self.config, candidate_set, chrom_seq)

    def _apply_junction_dominance_gate(
        self, candidate_set: CandidateSet, interval: GenomicInterval
    ) -> None:
        """Pre-EM junction-dominance gate. Delegates to the selection layer
        (fin.pipeline.selection.junction_dominance_select); kept as a method so tests can
        call this seam directly. Mutates ``candidate_set.candidates`` in place."""
        junction_dominance_select(
            self.config, candidate_set, interval, self._observed_junctions
        )

    def cleanup(self):
        """Close file handles."""
        if self._gtf_reader:
            self._gtf_reader.close()
        if self._signal_reader:
            self._signal_reader.close()
        genome = getattr(self, "_genome_fasta", None)
        if genome is not None and hasattr(genome, "close"):
            genome.close()

    def _load_genome_fasta(self, path: str):
        """Open the genome as a mapping (lazy indexed access by default).

        Memory attribution showed the eager whole-genome dict dominated
        worker RSS (~3.1 GB per spawn worker); the lazy mapping keeps the
        same ``chrom -> sequence`` contract while holding only a bounded
        chromosome cache per process. ``lazy_genome=False`` restores the
        historical eager dict.
        """
        from fin.io.lazy_genome import open_genome

        return open_genome(
            path,
            lazy=bool(getattr(self.config, "lazy_genome", True)),
            cache_chroms=int(getattr(self.config, "genome_cache_chroms", 2)),
        )

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
