"""Build M1/M2 EM input matrices for alignment and signal ablations.

Production M1/M2 modes do not use read-by-read coherence, so this helper always
returns a zero read-distance matrix and ``beta=0``.
"""
from __future__ import annotations

import logging
from typing import Tuple

import numpy as np

logger = logging.getLogger(__name__)


VALID_SUBSETS = ("m1", "m2", "m1+m2")


def _row_zscore(M: np.ndarray) -> np.ndarray:
    """Per-row z-score, treating sentinel 1e6 and non-finite as missing."""
    arr = np.asarray(M, dtype=np.float32)
    bad = (~np.isfinite(arr)) | (arr >= 1e5)
    clean = np.where(bad, np.nan, arr)
    mu = np.nanmean(clean, axis=1, keepdims=True)
    sd = np.nanstd(clean, axis=1, keepdims=True)
    sd = np.where(sd > 0, sd, 1.0)
    z = (clean - mu) / sd
    return np.where(np.isnan(z), 0.0, z).astype(np.float32)


def _zscore_mean(M1: np.ndarray, M2: np.ndarray) -> np.ndarray:
    """Per-row z-score each matrix then take element-wise mean.

    Output is shifted to be non-negative (subtract row-min) so it can be
    consumed by ``em_with_coherence`` which requires d >= 0.
    """
    z1 = _row_zscore(M1)
    z2 = _row_zscore(M2)
    mean = 0.5 * (z1 + z2)
    # Shift each row to start at 0 (preserves rank, satisfies EM contract).
    row_min = mean.min(axis=1, keepdims=True)
    return (mean - row_min).astype(np.float32)


def _row_center(M: np.ndarray) -> np.ndarray:
    """Per-row min subtraction, treating sentinel 1e6 and non-finite as missing.

    After centering the best candidate per row has distance 0, matching the
    semantic structure of M1 (mappy_distance.py:85-99).
    """
    arr = np.asarray(M, dtype=np.float32)
    bad = (~np.isfinite(arr)) | (arr >= 1e5)
    clean = np.where(bad, np.nan, arr)
    row_min = np.nanmin(clean, axis=1, keepdims=True)
    centered = clean - row_min
    return np.where(np.isnan(centered), 1e6, centered).astype(np.float32)


def _rank_fusion(M1: np.ndarray, M2: np.ndarray) -> np.ndarray:
    """Reciprocal rank fusion of M1 and M2, per row.

    For each read, rank candidates by M1 and M2 distances (lower distance =
    rank 1 = best).  Combined score = 1/rank_m1 + 1/rank_m2.  Convert back
    to a non-negative distance by subtracting from the row max so the best
    candidate gets distance 0.  Completely scale-invariant.
    """
    n_reads, n_cands = M1.shape
    result = np.zeros((n_reads, n_cands), dtype=np.float32)
    for i in range(n_reads):
        # argsort-of-argsort gives 0-based rank; +1 for 1-based
        rank_m1 = np.argsort(np.argsort(M1[i])).astype(np.float64) + 1.0
        rank_m2 = np.argsort(np.argsort(M2[i])).astype(np.float64) + 1.0
        rrf = 1.0 / rank_m1 + 1.0 / rank_m2
        # Best (highest rrf) → distance 0
        result[i] = (rrf.max() - rrf).astype(np.float32)
    return result


def build_em_matrices(
    subset: str,
    m1: np.ndarray,
    m2: np.ndarray,
    m2_norm: str = "none",
    m1m2_fusion: str = "zscore_mean",
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Assemble ``(d_tx, zero_d_rr, 0.0)`` for ``em_with_coherence``.

    ``subset`` may be ``m1``, ``m2``, or ``m1+m2``. The remaining arguments
    control M2 normalization and M1/M2 fusion for historical ablations.
    """
    if subset not in VALID_SUBSETS:
        raise ValueError(
            f"em_matrix_subset must be one of {VALID_SUBSETS}, got {subset!r}"
        )

    n_reads = m2.shape[0]
    zeros_rr = np.zeros((n_reads, n_reads), dtype=np.float32)

    # Apply M2 normalization when M2 is used standalone.
    if m2_norm == "center":
        m2_normed = _row_center(m2)
    elif m2_norm == "zscore":
        m2_normed = _row_zscore(m2)
        # Shift non-negative (EM contract)
        rmin = m2_normed.min(axis=1, keepdims=True)
        m2_normed = (m2_normed - rmin).astype(np.float32)
    else:
        m2_normed = m2.astype(np.float32)

    if subset == "m1":
        d_tx = m1.astype(np.float32)
    elif subset == "m2":
        d_tx = m2_normed
    elif m1m2_fusion == "rank":
        d_tx = _rank_fusion(m1, m2)
    else:
        d_tx = _zscore_mean(m1, m2)
    d_rr, beta_use = zeros_rr, 0.0

    # Defensive: EM requires non-negative distances. Floor any small
    # numerical drift below zero (can happen after z-scoring).
    if d_tx.min() < 0:
        d_tx = d_tx - d_tx.min()
    if d_rr.min() < 0:
        d_rr = d_rr - d_rr.min()
    return d_tx.astype(np.float32), d_rr.astype(np.float32), beta_use
