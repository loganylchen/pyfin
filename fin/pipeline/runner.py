"""Pipeline orchestrator: interval -> candidates -> scoring -> EM -> quantification."""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pysam

from fin.analysis.assignments import em_with_coherence
from fin.analysis.quantification import (
    QuantResult,
    _exons_from_candidate,
    aggregate_across_intervals,
    quantify_transcripts,
)
from fin.candidates.canonical import chain_all_canonical, parse_motifs
from fin.candidates.dataclasses import CandidateSet
from fin.candidates.discovery import discover_candidates, merge_fusion_candidates
from fin.io.interval_manager import GenomicInterval, generate_isolated_intervals
from fin.pipeline.config import PipelineConfig
from fin.scoring.composite import (
    derive_prior_weights,
    populate_quant_scores,
    score_candidates_composite,
    subsample_reads_for_dtw,
)
from fin.scoring.diff_region_dtw import compute_diff_region_m4
from fin.scoring.em_inputs import build_em_matrices
from fin.scoring.eventalign_parser import (
    build_distance_matrix,
    parse_eventalign_tsv,
)
from fin.scoring.external_tools import ExternalToolPaths, ExternalToolRunner
from fin.scoring.mappy_distance import compute_mappy_distance
from fin.scoring.signal_dtw import compute_read_to_read_dtw, extract_signal_segments
from fin.scoring.signal_tiebreak import signal_tiebreak
from fin.scoring.fake_fastq_tiebreak import fake_fastq_tiebreak
from fin.scoring.f5c_rna_tiebreak import f5c_rna_tiebreak

logger = logging.getLogger(__name__)


def _project_responsibilities_full(
    R_sub: np.ndarray,
    sub_read_ids: List[str],
    full_read_ids: List[str],
    dist_read_to_tx_full: np.ndarray,
    sigma: float,
    prior_weights: Optional[np.ndarray],
) -> tuple:
    """Project EM responsibilities from subsampled to full read set.

    Subsampled reads keep their EM responsibilities (which include coherence).
    Non-subsampled reads get a coherence-free softmax over the full d_tx using
    the same sigma + (optional) prior. This restores correct quantification
    when DTW subsampling is in effect.
    """
    n_full = len(full_read_ids)
    n_tx = R_sub.shape[1]
    sub_index = {rid: i for i, rid in enumerate(sub_read_ids)}

    R_full = np.zeros((n_full, n_tx), dtype=R_sub.dtype)

    # Coherence-free softmax for non-subsampled reads (numerically stable).
    d_tx = dist_read_to_tx_full
    d_min = d_tx.min(axis=1, keepdims=True)
    R_softmax = np.exp(-(d_tx - d_min) / max(sigma, 1e-6))
    if prior_weights is not None:
        R_softmax = R_softmax * np.asarray(prior_weights).reshape(1, -1)
    row_sums = np.maximum(R_softmax.sum(axis=1, keepdims=True), 1e-10)
    R_softmax = R_softmax / row_sums

    for i, rid in enumerate(full_read_ids):
        si = sub_index.get(rid)
        if si is not None:
            R_full[i, :] = R_sub[si, :]
        else:
            R_full[i, :] = R_softmax[i, :]

    hard_full = np.argmax(R_full, axis=1)
    return R_full, hard_full


@dataclass
class _IntervalBuildState:
    """Intermediate artifacts that are independent of em_matrix_subset.

    Built once per interval; reused across all EM-subset variants so the
    expensive Phase 1-3 work (candidate discovery, eventalign, DTW, mappy)
    is not repeated 7x for a 7-subset ablation sweep.
    """

    interval: GenomicInterval
    work_dir: Path
    candidate_set: CandidateSet
    read_ids: List[str]
    candidate_ids: List[str]
    dtw_read_ids: List[str]
    dist_read_to_tx: np.ndarray         # m2 (eventalign), rows=dtw_read_ids
    dist_read_to_tx_full: np.ndarray    # m2 over full read set (for projection)
    dist_read_to_read: np.ndarray       # m3 (DTW or diff-region)
    m1: Optional[np.ndarray]            # mappy distance, built lazily
    composite_scores: list
    combined_scores_arr: np.ndarray
    sigma_use: float
    em_max_iter_use: int
    scores_by_pair: Optional[dict] = None  # (read, cand) → ReadCandidateScore


class PipelineRunner:
    """Orchestrates the full pyfin pipeline."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self._gtf_reader = None
        self._genome_fasta = None
        self._signal_reader = None
        self._tool_runner = None

    def setup(self):
        """Validate tools, open file handles, load references."""
        # Validate external tools
        tool_paths = ExternalToolPaths(f5c=self.config.f5c_path)
        missing = tool_paths.validate()
        if missing:
            raise RuntimeError(f"Missing external tools: {', '.join(missing)}")

        # Create tool runner and build f5c index ONCE with all reads
        self._tool_runner = ExternalToolRunner(
            fastq_path=self.config.fastq_path,
            signal_path=self.config.signal_path,
            signal_format=self.config.signal_format,
            work_dir=self.config.work_dir,
            tools=tool_paths,
        )
        self._tool_runner.build_f5c_index()

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

    def run_subsets(
        self,
        subsets: List[str],
        output_paths: Dict[str, Tuple[Optional[str], Optional[str]]],
    ) -> Dict[str, Dict[str, QuantResult]]:
        """Run all intervals once, branch EM/quant per subset (matrix reuse).

        For each interval the expensive Phase 1-4 work (candidates, eventalign,
        DTW, mappy, composite scoring) is computed once via
        ``_build_interval_state``. Then for each subset in ``subsets`` we run
        only Phases 5-6 (``_em_score_and_quantify``). After all intervals are
        processed, per-subset results are aggregated, filtered, and written.

        Args:
            subsets: list of em_matrix_subset values, e.g. ['m2+m3', 'all'].
            output_paths: subset -> (gtf_path or None, tsv_path or None).

        Returns:
            subset -> aggregated quant results dict.
        """
        result = generate_isolated_intervals(
            self.config.bam_path,
            gtf_path=self.config.gtf_path,
            max_gap=self.config.max_gap,
            max_reads=self.config.max_reads,
        )
        intervals = result["intervals"]
        logger.info(
            "run_subsets: %d intervals × %d subsets = %d EM passes",
            len(intervals), len(subsets), len(intervals) * len(subsets),
        )

        per_subset_lists: Dict[str, List[List[QuantResult]]] = {
            s: [] for s in subsets
        }

        for i, interval in enumerate(intervals):
            logger.info(
                "Building state for interval %d/%d: %s",
                i + 1, len(intervals), interval.region_string,
            )
            state = self._build_interval_state(interval)
            if state is None:
                continue
            if isinstance(state, list):
                # R1 mappy-argmax bypass: same result for every subset.
                for s in subsets:
                    per_subset_lists[s].append(state)
                continue
            for s in subsets:
                logger.info(
                    "  subset=%s EM/quant for %s",
                    s, interval.region_string,
                )
                qr = self._em_score_and_quantify(
                    state, s, r_suffix=s.replace("+", "_")
                )
                if qr:
                    per_subset_lists[s].append(qr)

        finalized: Dict[str, Dict[str, QuantResult]] = {}
        for s in subsets:
            agg = aggregate_across_intervals(per_subset_lists[s])
            out_gtf, out_tsv = output_paths.get(s, (None, None))
            finalized[s] = self._finalize_and_write(
                agg, output_gtf=out_gtf, output_tsv=out_tsv
            )
        return finalized

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

        # Phase A Tier-2: max_R filter (FP-by-EM). EM responsibility is a
        # strong FP discriminator (Cohen's d=1.34 vs combined_score's 0.70).
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
        if _score_filter_on and self.config.min_novel_combined_score > 0.0:
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
        """Process a single interval; returns QuantResults for the configured subset."""
        state = self._build_interval_state(interval)
        if state is None:
            return None
        if isinstance(state, list):
            # R1 mappy-argmax path bypasses EM entirely.
            return state
        subset = getattr(self.config, "em_matrix_subset", "m2+m3")
        return self._em_score_and_quantify(state, subset)

    def _build_interval_state(
        self, interval: GenomicInterval
    ) -> "Optional[_IntervalBuildState | List[QuantResult]]":
        """Run Phases 1-4 (subset-independent) and return reusable state.

        Returns:
            - _IntervalBuildState for the normal EM path.
            - List[QuantResult] when R1 mappy-argmax bypass is engaged
              (enable_signal=False) — caller should treat it as a finished result.
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

        # R1 ablation: enable_signal=False → skip signal/EM, use mappy argmax.
        if not getattr(self.config, "enable_signal", True):
            return self._process_interval_mappy_argmax(
                candidate_set, read_ids, interval
            )

        # --- Phase 2: Per-candidate scoring (mappy + f5c eventalign) ---
        tsv_paths = self._tool_runner.score_candidates(candidate_set, work_dir)

        # Merge per-candidate eventalign results
        candidate_lengths = {
            c.candidate_id: len(c.sequence) for c in candidate_set.candidates
        }
        all_scores = []
        for tsv_path in tsv_paths:
            scores = parse_eventalign_tsv(str(tsv_path), candidate_lengths)
            all_scores.extend(scores)

        # Build scores lookup for signal tiebreak (keyed by (read, cand))
        scores_by_pair = {
            (s.read_name, s.candidate_id): s for s in all_scores
        }

        # DTW subsampling (m1): uniformly subsample reads when above cap.
        dtw_read_ids = subsample_reads_for_dtw(
            read_ids, self.config.max_reads_per_interval_for_dtw
        )
        if len(dtw_read_ids) != len(read_ids):
            logger.info(
                "Interval %s: subsampling DTW reads %d -> %d",
                interval.region_string,
                len(read_ids),
                len(dtw_read_ids),
            )

        # Read-to-tx distances (aligned with dtw_read_ids for m3 consistency).
        dist_read_to_tx = build_distance_matrix(
            all_scores, dtw_read_ids, candidate_ids,
            metric=self.config.m2_metric,
        )
        # Also build the full-read d_tx so quantification can include
        # non-subsampled reads (P0-1: prevent abundance underestimation when
        # max_reads_per_interval_for_dtw caps the EM input).
        dist_read_to_tx_full = build_distance_matrix(
            all_scores, read_ids, candidate_ids,
            metric=self.config.m2_metric,
        )

        # --- Phase 3: Signal DTW (read-to-read) ---
        # m4_source controls which read-to-read matrix is used (AC1).
        m4_src = getattr(self.config, "m4_source", "whole_read")
        if m4_src == "none":
            # AC5: zero matrix (no coherence); not None so EM shape checks pass.
            dist_read_to_read = np.zeros(
                (len(dtw_read_ids), len(dtw_read_ids)), dtype=np.float32
            )
        elif m4_src == "diff_region":
            # Collect per-event records needed by diff-region DTW (AC7-pre).
            scores_with_events = []
            cand_lengths_map = {
                c.candidate_id: len(c.sequence) for c in candidate_set.candidates
            }
            for tsv_path in tsv_paths:
                ev_scores = parse_eventalign_tsv(
                    str(tsv_path), cand_lengths_map, collect_events=True
                )
                scores_with_events.extend(ev_scores)
            scores_by_pair = {
                (s.read_name, s.candidate_id): s for s in scores_with_events
            }
            from fin.ablation.runner import _nan_to_row_mean

            m4_raw = compute_diff_region_m4(
                read_ids=dtw_read_ids,
                candidates=candidate_set.candidates,
                scores_by_pair=scores_by_pair,
                signal_reader=self._signal_reader,
                interval_start=interval.start,
                interval_end=interval.end,
                signal_format=self.config.signal_format,
                use_gpu=self.config.use_gpu,
                normalize=self.config.signal_normalize,
            )
            # AC8: sanitize NaN rows before EM.
            dist_read_to_read = _nan_to_row_mean(m4_raw)
        else:
            # Default: whole-read signal DTW.
            segments = extract_signal_segments(
                all_scores,
                self._signal_reader,
                signal_format=self.config.signal_format,
                normalize=self.config.signal_normalize,
            )
            dist_read_to_read = compute_read_to_read_dtw(
                segments, dtw_read_ids, use_gpu=self.config.use_gpu
            )

        # m3 alignment check: both distance matrices share the same read axis.
        if not (
            dist_read_to_tx.shape[0]
            == dist_read_to_read.shape[0]
            == len(dtw_read_ids)
        ):
            raise ValueError(
                f"Read axis mismatch: dist_read_to_tx rows={dist_read_to_tx.shape[0]}, "
                f"dist_read_to_read rows={dist_read_to_read.shape[0]}, "
                f"dtw_read_ids={len(dtw_read_ids)}"
            )
        if dist_read_to_read.shape[0] != dist_read_to_read.shape[1]:
            raise ValueError("dist_read_to_read must be square")

        n_reads = len(dtw_read_ids)
        n_tx = len(candidate_ids)

        # --- Phase 4: Composite scoring (uniform R seed) ---
        R_uniform = np.full((n_reads, n_tx), 1.0 / max(n_tx, 1))
        composite_scores = score_candidates_composite(
            candidates=candidate_set.candidates,
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_read=dist_read_to_read,
            R=R_uniform,
            alpha=self.config.score_alpha,
            use_gpu=self.config.use_gpu,
        )
        combined_scores_arr = np.array(
            [s.combined for s in composite_scores], dtype=float
        )

        # Data-adaptive sigma (P0-α): with per-event normalized d_tx the
        # absolute scale is dataset-dependent (~0.5-5 nats). A fixed sigma=1.0
        # often collapses R to one-hot, defeating EM. Use the median per-read
        # range as a robust scale; clip into [min_sigma, em_sigma_max].
        if n_tx >= 2 and n_reads > 0:
            d_max = dist_read_to_tx.max(axis=1)
            d_min = dist_read_to_tx.min(axis=1)
            adaptive = float(np.median(d_max - d_min))
            sigma_use = float(np.clip(
                adaptive if adaptive > 0 else self.config.em_sigma,
                getattr(self.config, "em_sigma_min", 0.05),
                getattr(self.config, "em_sigma_max", 50.0),
            ))
        else:
            sigma_use = self.config.em_sigma

        # R2 ablation: em_max_iter_override=1 forces single-step EM.
        em_max_iter_use = (
            self.config.em_max_iter_override
            if getattr(self.config, "em_max_iter_override", None) is not None
            else self.config.em_max_iter
        )

        return _IntervalBuildState(
            interval=interval,
            work_dir=work_dir,
            candidate_set=candidate_set,
            read_ids=read_ids,
            candidate_ids=list(candidate_ids),
            dtw_read_ids=dtw_read_ids,
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_tx_full=dist_read_to_tx_full,
            dist_read_to_read=dist_read_to_read,
            m1=None,
            composite_scores=composite_scores,
            combined_scores_arr=combined_scores_arr,
            sigma_use=sigma_use,
            em_max_iter_use=em_max_iter_use,
            scores_by_pair=scores_by_pair,
        )

    def _em_score_and_quantify(
        self,
        state: _IntervalBuildState,
        subset: str,
        r_suffix: Optional[str] = None,
    ) -> List[QuantResult]:
        """Phase 5+6: run EM with the chosen M-subset, then quantify.

        Reuses everything in ``state`` so candidate discovery, eventalign,
        DTW, and mappy are not recomputed across subsets.
        """
        # Prior is subset-independent (depends only on combined_scores_arr).
        prior_weights: Optional[np.ndarray] = None
        if self.config.use_prior:
            prior_weights = derive_prior_weights(
                state.combined_scores_arr,
                len(state.candidate_ids),
                self.config.prior_weight_cap,
            )

        # M-subset selector: m1 cached on state so a 7-subset sweep computes it once.
        if "m1" in subset:
            if state.m1 is None:
                read_sequences = (
                    getattr(state.candidate_set, "read_sequences", {}) or {}
                )
                state.m1 = compute_mappy_distance(
                    read_sequences,
                    state.candidate_set.candidates,
                    state.dtw_read_ids,
                )
            m1 = state.m1
        else:
            m1 = np.zeros_like(state.dist_read_to_tx)

        d_tx_em, d_rr_em, beta_use = build_em_matrices(
            subset=subset,
            m1=m1,
            m2=state.dist_read_to_tx,
            m3=state.dist_read_to_read,
            em_beta=self.config.em_beta,
            m2_norm=self.config.m2_norm,
            m1m2_fusion=self.config.m1m2_fusion,
        )

        R, hard_assignments, _log_likelihoods = em_with_coherence(
            dist_read_to_tx=d_tx_em,
            dist_read_to_read=d_rr_em,
            sigma=state.sigma_use,
            beta=beta_use,
            max_iter=state.em_max_iter_use,
            tol=self.config.em_tol,
            verbose=False,
            use_gpu=self.config.use_gpu,
            prior_weights=prior_weights,
        )

        # --- Phase 5.5: Signal tiebreak (optional) ---
        if (
            self.config.signal_tiebreak
            and self.config.enable_signal
            and state.scores_by_pair
        ):
            R = signal_tiebreak(
                R=R,
                read_ids=state.dtw_read_ids,
                candidates=state.candidate_set.candidates,
                scores_by_pair=state.scores_by_pair,
                ambig_threshold=self.config.tiebreak_ambig_threshold,
            )
            hard_assignments = R.argmax(axis=1)

        # --- Phase 5.6: Fake-FASTQ tiebreak (optional) ---
        if (
            self.config.fake_fastq_tiebreak
            and self.config.enable_signal
        ):
            # Build read_seqs from BAM for the tiebreak
            _tb_read_seqs = {}
            with pysam.AlignmentFile(self.config.bam_path, "rb") as _tb_bam:
                for _r in _tb_bam.fetch(
                    state.interval.chrom, state.interval.start, state.interval.end
                ):
                    if _r.is_secondary or _r.is_supplementary:
                        continue
                    if _r.query_sequence and _r.query_name not in _tb_read_seqs:
                        _tb_read_seqs[_r.query_name] = _r.query_sequence

            R = fake_fastq_tiebreak(
                R=R,
                read_ids=state.dtw_read_ids,
                read_seqs=_tb_read_seqs,
                candidates=state.candidate_set.candidates,
                signal_path=self.config.signal_path,
                f5c_path=self.config.f5c_path,
                ambig_threshold=self.config.tiebreak_ambig_threshold,
            )
            hard_assignments = R.argmax(axis=1)

        # --- Phase 5.7: f5c_rna in-memory tiebreak (optional) ---
        if self.config.f5c_rna_tiebreak:
            _tb_read_seqs = {}
            with pysam.AlignmentFile(self.config.bam_path, "rb") as _tb_bam:
                for _r in _tb_bam.fetch(
                    state.interval.chrom, state.interval.start, state.interval.end
                ):
                    if _r.is_secondary or _r.is_supplementary:
                        continue
                    if _r.query_sequence and _r.query_name not in _tb_read_seqs:
                        _tb_read_seqs[_r.query_name] = _r.query_sequence

            R = f5c_rna_tiebreak(
                R=R,
                read_ids=state.dtw_read_ids,
                read_seqs=_tb_read_seqs,
                candidates=state.candidate_set.candidates,
                signal_path=self.config.signal_path,
                pore=self.config.f5c_rna_pore,
                ambig_threshold=self.config.tiebreak_ambig_threshold,
            )
            hard_assignments = R.argmax(axis=1)

        # --- Phase 6: Quantification + score field population ---
        if len(state.dtw_read_ids) == len(state.read_ids):
            R_quant = R
            hard_quant = hard_assignments
            quant_read_ids = state.dtw_read_ids
        else:
            R_quant, hard_quant = _project_responsibilities_full(
                R_sub=R,
                sub_read_ids=state.dtw_read_ids,
                full_read_ids=state.read_ids,
                dist_read_to_tx_full=state.dist_read_to_tx_full,
                sigma=state.sigma_use,
                prior_weights=prior_weights,
            )
            quant_read_ids = state.read_ids

        quant_results = quantify_transcripts(
            R_quant, hard_quant, state.candidate_set.candidates, quant_read_ids
        )
        populate_quant_scores(quant_results, state.composite_scores)

        for j, qr in enumerate(quant_results):
            qr.max_R = float(R_quant[:, j].max()) if R_quant.shape[0] > 0 else 0.0

        # --- Post-hoc M2 validation (option 5) ---
        # After M1-based EM, check if M2 agrees with each read assignment.
        # Discount abundance for candidates where M2 doesn't support.
        if (
            self.config.m2_posthoc
            and self.config.enable_signal
            and hasattr(state, "dist_read_to_tx_full")
            and state.dist_read_to_tx_full is not None
        ):
            m2_full = state.dist_read_to_tx_full
            top_k = self.config.m2_posthoc_top_k
            n_reads_q = len(quant_read_ids)
            for j, qr in enumerate(quant_results):
                if qr.abundance <= 0:
                    continue
                assigned = np.where(hard_quant == j)[0]
                if len(assigned) == 0:
                    continue
                n_supported = 0
                for i in assigned:
                    if i >= m2_full.shape[0]:
                        n_supported += 1  # safety: count as supported
                        continue
                    m2_row = m2_full[i]
                    m2_j = m2_row[j]
                    if m2_j >= 1e5:
                        continue  # no M2 data for this pair
                    # Rank: how many valid candidates have lower M2 distance?
                    valid_mask = m2_row < 1e5
                    rank = int(np.sum(m2_row[valid_mask] < m2_j)) + 1
                    n_valid = int(valid_mask.sum())
                    k_use = min(top_k, max(1, n_valid // 3))
                    if rank <= k_use:
                        n_supported += 1
                support_frac = n_supported / len(assigned) if len(assigned) > 0 else 1.0
                old_ab = qr.abundance
                qr.abundance *= support_frac
                if support_frac < 1.0:
                    logger.debug(
                        "M2 posthoc: %s support=%.2f (%d/%d) abundance %.2f→%.2f",
                        qr.candidate_id, support_frac, n_supported,
                        len(assigned), old_ab, qr.abundance,
                    )

        if self.config.persist_R_matrix:
            if r_suffix:
                r_name = f"R_{r_suffix}.npy"
                meta_name = f"R_{r_suffix}_meta.json"
            else:
                r_name = "R.npy"
                meta_name = "R_meta.json"
            np.save(str(state.work_dir / r_name), R_quant.astype(np.float32))
            meta = {
                "read_ids": list(quant_read_ids),
                "candidate_ids": list(state.candidate_ids),
                "interval_region": state.interval.region_string,
                "sigma_used": state.sigma_use,
                "subset": subset,
                "was_subsampled": len(state.dtw_read_ids) != len(state.read_ids),
                "n_reads_subsampled": len(state.dtw_read_ids),
                "n_reads_full": len(state.read_ids),
            }
            with open(state.work_dir / meta_name, "w") as _f:
                json.dump(meta, _f)

        return quant_results

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
        # f5c_rna aligner; auto-skip the whole leg if signal is absent.
        m2_f5c = None
        m2_on = bool(getattr(self.config, "m2_tiebreak", False)) and bool(
            self.config.signal_path
        )
        if m2_on:
            try:
                import f5c_rna

                m2_f5c = f5c_rna.Aligner(
                    pore=self.config.f5c_rna_pore, use_gpu=False,
                    hmm_confidence=False,
                )
            except Exception as exc:  # signal stack unavailable -> keep 1/K split
                logger.warning("M2 tiebreak disabled (f5c_rna init failed): %s", exc)
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
                best_local, margin = m2_resolve_tie(
                    rid, seq, tied_cands, self.config.signal_path,
                    pore=self.config.f5c_rna_pore,
                    junction_k=self.config.m2_tiebreak_junction_k,
                    f5c_aligner=m2_f5c, mappy_aligners=tied_aligners,
                )
                if best_local is not None and margin >= self.config.m2_tiebreak_margin:
                    j = tied[best_local]
                    counts[j] += 1.0
                    assigned[j].append(rid)
                    if 1.0 > max_weight[j]:
                        max_weight[j] = 1.0
                    n_m2_override += 1
                    continue
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

    def _process_interval_mappy_argmax(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> Optional[List[QuantResult]]:
        """R1 ablation: AS-weighted multi-mapping baseline (no signal, no EM).

        Each read distributes responsibility proportional to mappy alignment
        score across all hit candidates (salmon/NanoCount-style). Returns a
        list of QuantResult with abundance = fractional column-sum.
        """
        variant = getattr(self.config, "r1_variant", "argmax_keep")

        # M1-keep split (production default): assign each read to ALL its
        # simultaneously-best-AS candidates, abundance split 1/K (read mass
        # conserved). Mirrors harness _quant_tie(mode="split") == M1-split.
        if variant == "argmax_keep":
            return self._quant_argmax_keep(candidate_set, read_ids, interval)

        # M1-first hard argmin: single best-AS pick per read, ties → lowest
        # candidate index (GTF prior).
        if variant == "argmax_first":
            return self._quant_argmax_first(candidate_set, read_ids, interval)

        from fin.ablation.mappy_argmax import (
            mappy_multimap_responsibilities,
            per_tx_counts_from_responsibilities,
        )

        candidates = candidate_set.candidates
        cand_list = list(candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        reads_iter = [
            (rid, read_sequences.get(rid, "")) for rid in read_ids
        ]
        reads_iter = [(rid, seq) for rid, seq in reads_iter if seq]

        R_mm, kept_read_ids = mappy_multimap_responsibilities(reads_iter, cand_list)

        if variant == "argmax_em" and R_mm.size > 0:
            # R1c: real m1 mappy-edit-distance matrix as EM dist_read_to_tx
            # (β=0, dist_read_to_read=zeros). Pure-alignment EM — no signal
            # involved. m1 carries strictly more information than R_mm at
            # this point: it preserves per-read AS gaps in score units rather
            # than collapsing to normalized shares.
            from fin.scoring.em_inputs import build_em_matrices
            from fin.scoring.mappy_distance import compute_mappy_distance

            read_seqs = {rid: seq for rid, seq in reads_iter}
            n_reads_em = len(kept_read_ids)
            n_cands_em = len(cand_list)
            max_iter_em = (
                self.config.em_max_iter_override
                if self.config.em_max_iter_override is not None
                else self.config.em_max_iter
            )

            scoring = getattr(self.config, "r1_scoring", "mappy")
            use_em = (max_iter_em > 1)

            if scoring == "f5c_rna":
                from fin.scoring.f5c_rna_tiebreak import _build_m2_f5c_rna

                if use_em:
                    # R2c: build M2 distance, run EM
                    dist_read_to_tx = _build_m2_f5c_rna(
                        kept_read_ids, read_seqs, cand_list,
                        self.config.signal_path, self.config.f5c_rna_pore,
                        as_distance=True,
                    )
                    dist_read_to_read = np.zeros(
                        (n_reads_em, n_reads_em), dtype=np.float32)
                    R, hard_assignments, _ = em_with_coherence(
                        dist_read_to_tx=dist_read_to_tx,
                        dist_read_to_read=dist_read_to_read,
                        sigma=self.config.em_sigma,
                        beta=0.0,
                        max_iter=max_iter_em,
                        tol=self.config.em_tol,
                        verbose=False,
                        use_gpu=self.config.use_gpu,
                    )
                else:
                    # R2a: normalized scores as R directly (no EM)
                    R = _build_m2_f5c_rna(
                        kept_read_ids, read_seqs, cand_list,
                        self.config.signal_path, self.config.f5c_rna_pore,
                        as_distance=False,
                    )
                    hard_assignments = R.argmax(axis=1)

                # f5c_rna tiebreak (if enabled, further refine)
                if self.config.f5c_rna_tiebreak:
                    R = f5c_rna_tiebreak(
                        R=R,
                        read_ids=kept_read_ids,
                        read_seqs=read_seqs,
                        candidates=cand_list,
                        signal_path=self.config.signal_path,
                        pore=self.config.f5c_rna_pore,
                        ambig_threshold=self.config.tiebreak_ambig_threshold,
                    )
                    hard_assignments = R.argmax(axis=1)

                quant_results = quantify_transcripts(
                    R, hard_assignments, cand_list, kept_read_ids
                )
                for j, qr in enumerate(quant_results):
                    qr.max_R = float(R[:, j].max()) if R.shape[0] > 0 else 0.0
                logger.info(
                    "R2 f5c_rna interval %s: %d reads -> %d candidates",
                    interval.region_string, len(read_ids), len(candidates),
                )
                return quant_results
            else:
                # M1: mappy AS-gap distance → normalized scores (sum=1)
                m1 = compute_mappy_distance(read_seqs, cand_list, kept_read_ids)
                if use_em:
                    # M1 + EM: use distance matrix through EM softmax (R1c)
                    m2_dummy = np.zeros((n_reads_em, n_cands_em), dtype=np.float32)
                    m3_dummy = np.zeros((n_reads_em, n_reads_em), dtype=np.float32)
                    dist_read_to_tx, _, _ = build_em_matrices(
                        "m1", m1, m2_dummy, m3_dummy, em_beta=0.0,
                    )
                    dist_read_to_read = np.zeros(
                        (n_reads_em, n_reads_em), dtype=np.float32)
                    R, hard_assignments, _ = em_with_coherence(
                        dist_read_to_tx=dist_read_to_tx,
                        dist_read_to_read=dist_read_to_read,
                        sigma=self.config.em_sigma,
                        beta=0.0,
                        max_iter=max_iter_em,
                        tol=self.config.em_tol,
                        verbose=False,
                        use_gpu=self.config.use_gpu,
                    )
                else:
                    # M1 normalized scores: 1/(1+dist), row-sum=1
                    m1_scores = 1.0 / (1.0 + m1.astype(np.float64))
                    row_sums = m1_scores.sum(axis=1, keepdims=True)
                    row_sums = np.where(row_sums > 0, row_sums, 1.0)
                    R = (m1_scores / row_sums).astype(np.float32)
                    hard_assignments = R.argmax(axis=1)

                # f5c_rna tiebreak in R1 path
                if self.config.f5c_rna_tiebreak:
                    R = f5c_rna_tiebreak(
                        R=R,
                        read_ids=kept_read_ids,
                        read_seqs=read_seqs,
                        candidates=cand_list,
                        signal_path=self.config.signal_path,
                        pore=self.config.f5c_rna_pore,
                        ambig_threshold=self.config.tiebreak_ambig_threshold,
                    )
                    hard_assignments = R.argmax(axis=1)

                quant_results = quantify_transcripts(
                    R, hard_assignments, cand_list, kept_read_ids
                )
                for j, qr in enumerate(quant_results):
                    qr.max_R = float(R[:, j].max()) if R.shape[0] > 0 else 0.0
                logger.info(
                    "R1 mappy interval %s: %d reads -> %d candidates",
                    interval.region_string, len(read_ids), len(candidates),
                )
                return quant_results
            R, hard_assignments, _ = em_with_coherence(
                dist_read_to_tx=dist_read_to_tx,
                dist_read_to_read=dist_read_to_read,
                sigma=self.config.em_sigma,
                beta=beta_use,
                max_iter=max_iter_em,
                tol=self.config.em_tol,
                verbose=False,
                use_gpu=self.config.use_gpu,
            )

            # f5c_rna tiebreak in R1 path (no signal pipeline needed)
            if self.config.f5c_rna_tiebreak:
                R = f5c_rna_tiebreak(
                    R=R,
                    read_ids=kept_read_ids,
                    read_seqs=read_seqs,
                    candidates=cand_list,
                    signal_path=self.config.signal_path,
                    pore=self.config.f5c_rna_pore,
                    ambig_threshold=self.config.tiebreak_ambig_threshold,
                )
                hard_assignments = R.argmax(axis=1)

            quant_results = quantify_transcripts(
                R, hard_assignments, cand_list, kept_read_ids
            )
            for j, qr in enumerate(quant_results):
                qr.max_R = float(R[:, j].max()) if R.shape[0] > 0 else 0.0
            logger.info(
                "R1 mappy+EM interval %s: %d reads -> %d candidates",
                interval.region_string,
                len(read_ids),
                len(candidates),
            )
            return quant_results

        # variant == "argmax_only": existing AS-weighted multimap baseline.
        counts = per_tx_counts_from_responsibilities(R_mm, cand_list)

        # Hard argmax per read for reporting `assigned_read_ids`.
        argmax_by_read: Dict[str, str] = {}
        if R_mm.size > 0:
            arg = np.argmax(R_mm, axis=1)
            for i, rid in enumerate(kept_read_ids):
                argmax_by_read[rid] = cand_list[int(arg[i])].candidate_id

        quant_results: List[QuantResult] = []
        for j, cand in enumerate(cand_list):
            cnt = counts.get(cand.candidate_id, 0.0)
            assigned = [
                rid for rid, cid in argmax_by_read.items()
                if cid == cand.candidate_id
            ]
            col_max = float(R_mm[:, j].max()) if R_mm.size > 0 else 0.0
            qr = QuantResult(
                candidate_id=cand.candidate_id,
                abundance=cnt,
                confidence=col_max,
                num_assigned_reads=len(assigned),
                source=cand.source,
                chrom=cand.chrom,
                strand=cand.strand,
                start=cand.start,
                end=cand.end,
                exons=_exons_from_candidate(cand),
                assigned_read_ids=tuple(assigned),
            )
            qr.max_R = col_max
            quant_results.append(qr)

        logger.info(
            "R1 mappy multimap interval %s: %d reads -> %d candidates",
            interval.region_string,
            len(read_ids),
            len(candidates),
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
