"""Ablation runner: sweeps quant-mode configurations over a fixed interval set.

Each ablation row (R1..R5) maps to a PipelineConfig override that isolates one
component of the pipeline. The benchmark driver (``benchmarks/ablation_5rows.py``)
runs the full ``PipelineRunner`` once per row, applying the row's overrides via
``dataclasses.replace``; all signal scoring is in-memory krill (no f5c CLI).

Row definitions (AC1):
  R1a quant_mode="argmax"           mappy argmax + M2 tiebreak, no EM
  R1b quant_mode="argmax"           + score filters on
  R1c quant_mode="m1_em"            mappy-distance EM (β=0), score filters on
  R2  quant_mode="m2_em",           single-step krill EM + DTW coherence
      em_max_iter_override=1,        (isolates iteration count vs R4)
      m3_coherence=True
  R3  quant_mode="m2_em",           krill EM, coherence disabled (β=0)
      m3_coherence=False
  R4  quant_mode="m2_em",           krill EM + junction-window DTW coherence
      m3_coherence=True
  R5  quant_mode="m2_em",           R4 + score filters on (production-equivalent)
      m3_coherence=True,
      enable_score_filter=True

Note: the m2_em assembly path gates M3 read×read DTW coherence on
``PipelineConfig.m3_coherence`` (NOT ``m4_source``). ``m4_source`` is retained
because the ablation rows below still use it to select the coherence source.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from fin.analysis.quantification import QuantResult
from fin.pipeline.config import PipelineConfig


@dataclass
class AblationRowConfig:
    """One row in the ablation sweep."""

    row_id: str          # e.g. "R1a", "R2", ...
    label: str           # human-readable description
    quant_mode: str = "m2_em"   # "argmax" | "m1_em" | "m2_em"
    em_max_iter_override: Optional[int] = None
    m4_source: str = "diff_region"   # "diff_region" | "none" (ablation coherence src)
    m3_coherence: bool = False  # gate read×read DTW coherence in the m2_em path
    enable_score_filter: bool = True


# Canonical ablation rows (AC1).
ABLATION_ROWS: List[AblationRowConfig] = [
    AblationRowConfig(
        row_id="R1a",
        label="mappy_argmax_only",
        quant_mode="argmax",
        em_max_iter_override=None,
        m4_source="none",
        enable_score_filter=False,
    ),
    AblationRowConfig(
        row_id="R1b",
        label="mappy_argmax_filtered",
        quant_mode="argmax",
        em_max_iter_override=None,
        m4_source="none",
        enable_score_filter=True,
    ),
    AblationRowConfig(
        row_id="R1c",
        label="mappy_em_filtered",
        quant_mode="m1_em",
        em_max_iter_override=None,
        m4_source="none",
        enable_score_filter=True,
    ),
    AblationRowConfig(
        row_id="R2",
        label="single_step_em",
        quant_mode="m2_em",
        em_max_iter_override=1,
        m4_source="diff_region",
        m3_coherence=True,
        enable_score_filter=False,
    ),
    AblationRowConfig(
        row_id="R3",
        label="em_no_coherence",
        quant_mode="m2_em",
        em_max_iter_override=None,
        m4_source="none",
        m3_coherence=False,
        enable_score_filter=False,
    ),
    AblationRowConfig(
        row_id="R4",
        label="diff_region_dtw",
        quant_mode="m2_em",
        em_max_iter_override=None,
        m4_source="diff_region",
        m3_coherence=True,
        enable_score_filter=False,
    ),
    AblationRowConfig(
        row_id="R5",
        label="diff_region_dtw_scored",
        quant_mode="m2_em",
        em_max_iter_override=None,
        m4_source="diff_region",
        m3_coherence=True,
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
