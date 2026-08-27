"""Ablation runner for active quantification configurations.

Each row maps to a ``PipelineConfig`` override applied by
``benchmarks/ablation_5rows.py``. The former R4 read-by-read DTW coherence row
was retired with M3; its prototype is preserved under
``experiments/m3_coherence``.

Active rows:
  R1a  mappy argmax + M2 tiebreak, filters off
  R1b  mappy argmax + M2 tiebreak, filters on
  R1c  mappy-distance EM, filters on
  R2   single-step M2-seeded EM, filters off
  R3   converged M2-seeded EM, filters off
  R5   converged M2-seeded EM, filters on
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

    row_id: str
    label: str
    quant_mode: str = "m2_em"
    em_max_iter_override: Optional[int] = None
    enable_score_filter: bool = True


ABLATION_ROWS: List[AblationRowConfig] = [
    AblationRowConfig(
        row_id="R1a",
        label="mappy_argmax_only",
        quant_mode="argmax",
        enable_score_filter=False,
    ),
    AblationRowConfig(
        row_id="R1b",
        label="mappy_argmax_filtered",
        quant_mode="argmax",
        enable_score_filter=True,
    ),
    AblationRowConfig(
        row_id="R1c",
        label="mappy_em_filtered",
        quant_mode="m1_em",
        enable_score_filter=True,
    ),
    AblationRowConfig(
        row_id="R2",
        label="single_step_em",
        quant_mode="m2_em",
        em_max_iter_override=1,
        enable_score_filter=False,
    ),
    AblationRowConfig(
        row_id="R3",
        label="m2_em_unfiltered",
        quant_mode="m2_em",
        enable_score_filter=False,
    ),
    AblationRowConfig(
        row_id="R5",
        label="m2_em_filtered",
        quant_mode="m2_em",
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
    """AC8: replace NaN cells with the row mean; all-NaN rows become zero."""
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
    """Apply the min_abundance novel/fusion-exempt floor filter."""
    if config.min_abundance > 0.0:
        aggregated = {
            cid: qr
            for cid, qr in aggregated.items()
            if qr.source in ("gtf", "fusion") or qr.abundance >= config.min_abundance
        }
    return aggregated
