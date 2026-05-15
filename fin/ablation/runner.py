"""Ablation runner: sweeps Stage 1-5 configurations over a fixed interval set.

Each ablation row (R1..R5) maps to a PipelineConfig override that isolates one
component of the pipeline. The runner processes a pre-built list of intervals
using an already-opened signal reader and tool runner, so the expensive setup
(index build, file open) happens once per benchmark run.

Row definitions (AC1):
  R1  enable_signal=False          mappy argmax only, no EM/DTW
  R2  em_max_iter_override=1       single-step EM (warm start from prior)
  R3  m4_source="none"             EM with coherence disabled (beta irrelevant)
  R4  m4_source="diff_region"      diff-region DTW as m4 source
  R5  m4_source="diff_region",     diff-region DTW + score filters on
      enable_score_filter=True     (production-equivalent)

The default (no override) is the full production pipeline (R5 semantics without
the ablation label).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from fin.analysis.assignments import em_with_coherence
from fin.analysis.quantification import (
    QuantResult,
    aggregate_across_intervals,
    quantify_transcripts,
)
from fin.candidates.dataclasses import CandidateSet
from fin.io.interval_manager import GenomicInterval
from fin.pipeline.config import PipelineConfig
from fin.scoring.composite import (
    derive_prior_weights,
    populate_quant_scores,
    score_candidates_composite,
    subsample_reads_for_dtw,
)
from fin.scoring.diff_region_dtw import compute_diff_region_m4
from fin.scoring.eventalign_parser import (
    ReadCandidateScore,
    build_distance_matrix,
    parse_eventalign_tsv,
)
from fin.scoring.signal_dtw import compute_read_to_read_dtw, extract_signal_segments

logger = logging.getLogger(__name__)


@dataclass
class AblationRowConfig:
    """One row in the ablation sweep."""

    row_id: str          # e.g. "R1", "R2", ...
    label: str           # human-readable description
    enable_signal: bool = True
    em_max_iter_override: Optional[int] = None
    m4_source: str = "whole_read"   # "whole_read" | "diff_region" | "none"
    enable_score_filter: bool = True


# Canonical Stage 1-5 rows (AC1).
ABLATION_ROWS: List[AblationRowConfig] = [
    AblationRowConfig(
        row_id="R1",
        label="mappy_argmax_only",
        enable_signal=False,
        em_max_iter_override=None,
        m4_source="whole_read",
        enable_score_filter=False,
    ),
    AblationRowConfig(
        row_id="R2",
        label="single_step_em",
        enable_signal=True,
        em_max_iter_override=1,
        m4_source="none",
        enable_score_filter=False,
    ),
    AblationRowConfig(
        row_id="R3",
        label="em_no_coherence",
        enable_signal=True,
        em_max_iter_override=None,
        m4_source="none",
        enable_score_filter=False,
    ),
    AblationRowConfig(
        row_id="R4",
        label="diff_region_dtw",
        enable_signal=True,
        em_max_iter_override=None,
        m4_source="diff_region",
        enable_score_filter=False,
    ),
    AblationRowConfig(
        row_id="R5",
        label="diff_region_dtw_scored",
        enable_signal=True,
        em_max_iter_override=None,
        m4_source="diff_region",
        enable_score_filter=True,
    ),
]


@dataclass
class AblationIntervalResult:
    """Per-interval result for one ablation row."""

    row_id: str
    interval_region: str
    quant_results: List[QuantResult]
    n_reads: int = 0
    n_candidates: int = 0


@dataclass
class AblationResult:
    """Aggregated result for one ablation row across all intervals."""

    row_id: str
    label: str
    aggregated: Dict[str, QuantResult] = field(default_factory=dict)
    interval_results: List[AblationIntervalResult] = field(default_factory=list)


def _nan_to_row_mean(m: np.ndarray) -> np.ndarray:
    """AC8: replace NaN cells with the row mean; if entire row is NaN use 0."""
    out = m.copy()
    for i in range(out.shape[0]):
        row = out[i]
        nan_mask = np.isnan(row)
        if nan_mask.all():
            out[i] = 0.0
        elif nan_mask.any():
            row_mean = float(np.nanmean(row))
            out[i, nan_mask] = row_mean
    return out


def _make_m4(
    row_cfg: AblationRowConfig,
    dtw_read_ids: List[str],
    candidates,
    scores_by_pair: Dict[Tuple[str, str], ReadCandidateScore],
    signal_reader,
    interval: GenomicInterval,
    all_scores: List[ReadCandidateScore],
    config: PipelineConfig,
) -> np.ndarray:
    """Build the read-to-read distance matrix (m4) for the given row config."""
    n = len(dtw_read_ids)
    src = row_cfg.m4_source

    if src == "none":
        # AC5: zero matrix, not None
        return np.zeros((n, n), dtype=np.float32)

    if src == "diff_region":
        m4 = compute_diff_region_m4(
            read_ids=dtw_read_ids,
            candidates=list(candidates),
            scores_by_pair=scores_by_pair,
            signal_reader=signal_reader,
            interval_start=interval.start,
            interval_end=interval.end,
            signal_format=config.signal_format,
            use_gpu=config.use_gpu,
            normalize=config.signal_normalize,
        )
        # AC8: sanitize NaN cells before passing to EM
        return _nan_to_row_mean(m4)

    # "whole_read" — standard signal DTW path
    segments = extract_signal_segments(
        all_scores,
        signal_reader,
        signal_format=config.signal_format,
        normalize=config.signal_normalize,
    )
    m4 = compute_read_to_read_dtw(segments, dtw_read_ids, use_gpu=config.use_gpu)
    return m4


def run_ablation_row(
    row_cfg: AblationRowConfig,
    intervals: List[GenomicInterval],
    candidate_sets: Dict[str, CandidateSet],
    tsv_paths_by_interval: Dict[str, List[str]],
    signal_reader,
    config: PipelineConfig,
) -> AblationResult:
    """Run one ablation row over all intervals.

    Args:
        row_cfg: The ablation row configuration.
        intervals: List of GenomicInterval to process.
        candidate_sets: interval.region_string -> CandidateSet (pre-built).
        tsv_paths_by_interval: interval.region_string -> [tsv_path, ...].
        signal_reader: Open Pod5Reader or Slow5Reader (shared, read-only).
        config: Base PipelineConfig; ablation fields are overridden per row.

    Returns:
        AblationResult with per-interval and aggregated quant results.
    """
    result = AblationResult(row_id=row_cfg.row_id, label=row_cfg.label)

    all_quant_lists: List[List[QuantResult]] = []
    for interval in intervals:
        region = interval.region_string
        candidate_set = candidate_sets.get(region)
        if candidate_set is None or candidate_set.num_candidates == 0:
            continue
        if not candidate_set.read_ids:
            continue

        tsv_paths = tsv_paths_by_interval.get(region, [])

        interval_result = _process_interval_for_row(
            row_cfg=row_cfg,
            interval=interval,
            candidate_set=candidate_set,
            tsv_paths=tsv_paths,
            signal_reader=signal_reader,
            config=config,
        )
        if interval_result is not None:
            result.interval_results.append(interval_result)
            all_quant_lists.append(interval_result.quant_results)

    # Apply post-EM filters only if enable_score_filter is on (AC7)
    aggregated = aggregate_across_intervals(all_quant_lists)
    if row_cfg.enable_score_filter:
        aggregated = _apply_score_filters(aggregated, config)

    result.aggregated = aggregated
    return result


def _process_interval_for_row(
    row_cfg: AblationRowConfig,
    interval: GenomicInterval,
    candidate_set: CandidateSet,
    tsv_paths: List[str],
    signal_reader,
    config: PipelineConfig,
) -> Optional[AblationIntervalResult]:
    """Process one interval for one ablation row.

    Returns AblationIntervalResult or None if nothing to do.
    """
    read_ids = sorted(candidate_set.read_ids)
    candidate_ids = candidate_set.candidate_ids()
    candidates = candidate_set.candidates
    n_tx = len(candidate_ids)

    # R1: skip signal entirely — AS-weighted multi-mapping baseline (AC1).
    # Each read distributes responsibility proportional to mappy AS across all
    # hit candidates (salmon/NanoCount-style alignment-only baseline). This
    # preserves alignment uncertainty unlike hard argmax.
    if not row_cfg.enable_signal:
        from fin.ablation.mappy_argmax import (
            mappy_multimap_responsibilities,
            per_tx_counts_from_responsibilities,
        )

        reads_iter = [
            (rid, _seq_for_read(rid, candidate_set)) for rid in read_ids
        ]
        reads_iter = [(rid, seq) for rid, seq in reads_iter if seq]
        cand_list = list(candidates)
        R_mm, kept_read_ids = mappy_multimap_responsibilities(reads_iter, cand_list)
        counts = per_tx_counts_from_responsibilities(R_mm, cand_list)

        # Hard argmax (per row) for assigned_read_ids / num_assigned_reads
        # so downstream reporting still has a "primary candidate per read".
        argmax_by_read: Dict[str, str] = {}
        if R_mm.size > 0:
            arg = np.argmax(R_mm, axis=1)
            for i, rid in enumerate(kept_read_ids):
                argmax_by_read[rid] = cand_list[int(arg[i])].candidate_id

        quant_results = []
        for j, cand in enumerate(cand_list):
            cnt = counts.get(cand.candidate_id, 0.0)
            assigned = [rid for rid, cid in argmax_by_read.items() if cid == cand.candidate_id]
            # Confidence = max responsibility this candidate received from any
            # read (analogous to max_R in the EM path). 0.0 when no reads hit.
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
                assigned_read_ids=tuple(assigned),
            )
            # Populate max_R so downstream filters (R5 enable_score_filter) can
            # treat R1 output uniformly with the EM rows.
            qr.max_R = col_max
            quant_results.append(qr)

        return AblationIntervalResult(
            row_id=row_cfg.row_id,
            interval_region=interval.region_string,
            quant_results=quant_results,
            n_reads=len(read_ids),
            n_candidates=n_tx,
        )

    # Signal path (R2-R5)
    # Parse eventalign TSVs — collect_events for diff_region m4 source (AC7-pre)
    collect_events = (row_cfg.m4_source == "diff_region")
    candidate_lengths = {c.candidate_id: len(c.sequence) for c in candidates}
    all_scores: List[ReadCandidateScore] = []
    for tsv_path in tsv_paths:
        scores = parse_eventalign_tsv(
            str(tsv_path),
            candidate_lengths,
            collect_events=collect_events,
        )
        all_scores.extend(scores)

    # DTW subsampling
    dtw_read_ids = subsample_reads_for_dtw(
        read_ids, config.max_reads_per_interval_for_dtw
    )

    dist_read_to_tx = build_distance_matrix(all_scores, dtw_read_ids, candidate_ids)
    dist_read_to_tx_full = build_distance_matrix(all_scores, read_ids, candidate_ids)

    # Build scores_by_pair lookup for diff_region path
    scores_by_pair: Dict[Tuple[str, str], ReadCandidateScore] = {}
    if collect_events:
        for s in all_scores:
            scores_by_pair[(s.read_name, s.candidate_id)] = s

    # m4 matrix
    dist_read_to_read = _make_m4(
        row_cfg=row_cfg,
        dtw_read_ids=dtw_read_ids,
        candidates=candidates,
        scores_by_pair=scores_by_pair,
        signal_reader=signal_reader,
        interval=interval,
        all_scores=all_scores,
        config=config,
    )

    n_reads = len(dtw_read_ids)
    n_cands = len(candidate_ids)

    # Composite scoring
    R_uniform = np.full((n_reads, n_cands), 1.0 / max(n_cands, 1))
    composite_scores = score_candidates_composite(
        candidates=list(candidates),
        dist_read_to_tx=dist_read_to_tx,
        dist_read_to_read=dist_read_to_read,
        R=R_uniform,
        alpha=config.score_alpha,
        use_gpu=config.use_gpu,
    )
    combined_scores_arr = np.array(
        [s.combined for s in composite_scores], dtype=float
    )

    # Prior
    prior_weights: Optional[np.ndarray] = None
    if config.use_prior:
        prior_weights = derive_prior_weights(
            combined_scores_arr, n_cands, config.prior_weight_cap
        )

    # Adaptive sigma
    if n_cands >= 2 and n_reads > 0:
        d_max = dist_read_to_tx.max(axis=1)
        d_min = dist_read_to_tx.min(axis=1)
        adaptive = float(np.median(d_max - d_min))
        sigma_use = float(np.clip(
            adaptive if adaptive > 0 else config.em_sigma,
            getattr(config, "em_sigma_min", 0.05),
            getattr(config, "em_sigma_max", 50.0),
        ))
    else:
        sigma_use = config.em_sigma

    # EM max_iter override (R2)
    max_iter = (
        row_cfg.em_max_iter_override
        if row_cfg.em_max_iter_override is not None
        else config.em_max_iter
    )

    R, hard_assignments, _ = em_with_coherence(
        dist_read_to_tx=dist_read_to_tx,
        dist_read_to_read=dist_read_to_read,
        sigma=sigma_use,
        beta=config.em_beta,
        max_iter=max_iter,
        tol=config.em_tol,
        verbose=False,
        use_gpu=config.use_gpu,
        prior_weights=prior_weights,
    )

    # Project back to full read set if subsampled
    if len(dtw_read_ids) == len(read_ids):
        R_quant = R
        hard_quant = hard_assignments
        quant_read_ids = dtw_read_ids
    else:
        from fin.pipeline.runner import _project_responsibilities_full

        R_quant, hard_quant = _project_responsibilities_full(
            R_sub=R,
            sub_read_ids=dtw_read_ids,
            full_read_ids=read_ids,
            dist_read_to_tx_full=dist_read_to_tx_full,
            sigma=sigma_use,
            prior_weights=prior_weights,
        )
        quant_read_ids = read_ids

    quant_results = quantify_transcripts(
        R_quant, hard_quant, list(candidates), quant_read_ids
    )
    populate_quant_scores(quant_results, composite_scores)

    for j, qr in enumerate(quant_results):
        qr.max_R = float(R_quant[:, j].max()) if R_quant.shape[0] > 0 else 0.0

    return AblationIntervalResult(
        row_id=row_cfg.row_id,
        interval_region=interval.region_string,
        quant_results=quant_results,
        n_reads=len(read_ids),
        n_candidates=n_cands,
    )


def _apply_score_filters(
    aggregated: Dict[str, QuantResult],
    config: PipelineConfig,
) -> Dict[str, QuantResult]:
    """Apply min_abundance, min_max_r, min_novel_combined_score filters."""
    if config.min_abundance > 0.0:
        aggregated = {
            cid: qr
            for cid, qr in aggregated.items()
            if qr.source in ("gtf", "fusion") or qr.abundance >= config.min_abundance
        }
    if config.min_max_r > 0.0:
        aggregated = {
            cid: qr
            for cid, qr in aggregated.items()
            if qr.source in ("gtf", "fusion") or qr.max_R >= config.min_max_r
        }
    if config.min_novel_combined_score > 0.0:
        aggregated = {
            cid: qr
            for cid, qr in aggregated.items()
            if qr.source in ("gtf", "fusion")
            or qr.combined_score >= config.min_novel_combined_score
        }
    return aggregated


def _seq_for_read(read_id: str, candidate_set: CandidateSet) -> str:
    """Extract sequence for a read_id from the candidate set read sequences.

    CandidateSet stores read sequences in read_sequences if available.
    Falls back to empty string (read will be dropped by mappy_argmax).
    """
    seq = getattr(candidate_set, "read_sequences", {})
    if isinstance(seq, dict):
        return seq.get(read_id, "")
    return ""
