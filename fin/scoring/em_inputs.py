"""Build EM input matrices from M1/M2/M3 according to em_matrix_subset.

This module assembles ``dist_read_to_tx`` and ``dist_read_to_read`` for
``em_with_coherence`` given the seven supported ablation subsets and the
already-computed M1/M2/M3 matrices.

Public function:
    build_em_matrices(subset, m1, m2, m3, em_beta) -> (d_tx, d_rr, beta_use)
"""
from __future__ import annotations

import logging
from typing import Tuple

import numpy as np

logger = logging.getLogger(__name__)


VALID_SUBSETS = ("m1", "m2", "m3", "m1+m2", "m1+m3", "m2+m3", "all")


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


def build_em_matrices(
    subset: str,
    m1: np.ndarray,
    m2: np.ndarray,
    m3: np.ndarray,
    em_beta: float,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Assemble ``(d_tx, d_rr, beta_use)`` for ``em_with_coherence``.

    Args:
        subset: One of :data:`VALID_SUBSETS`. Selects which of M1/M2/M3
            enter EM. Subsets containing "m3" pass ``em_beta`` through;
            subsets without "m3" force ``beta=0`` (no coherence term).
        m1: (n_reads, n_tx) mappy distance matrix. Pass zeros if unused.
        m2: (n_reads, n_tx) eventalign distance matrix. Pass zeros if unused.
        m3: (n_reads, n_reads) read-read distance matrix. Pass zeros if
            unused.
        em_beta: Coherence weight (used only when subset includes "m3").

    Returns:
        d_tx: (n_reads, n_tx) non-negative read×tx distance.
        d_rr: (n_reads, n_reads) non-negative read×read distance.
        beta_use: Effective beta (0 when M3 not selected).

    Raises:
        ValueError: If ``subset`` is not in :data:`VALID_SUBSETS`.
    """
    if subset not in VALID_SUBSETS:
        raise ValueError(
            f"em_matrix_subset must be one of {VALID_SUBSETS}, got {subset!r}"
        )

    n_reads, n_tx = m2.shape
    zeros_tx = np.zeros((n_reads, n_tx), dtype=np.float32)
    zeros_rr = np.zeros((n_reads, n_reads), dtype=np.float32)

    if subset == "m1":
        d_tx, d_rr, beta_use = m1.astype(np.float32), zeros_rr, 0.0
    elif subset == "m2":
        d_tx, d_rr, beta_use = m2.astype(np.float32), zeros_rr, 0.0
    elif subset == "m3":
        # No read×tx information; rely entirely on coherence.
        d_tx, d_rr, beta_use = zeros_tx, m3.astype(np.float32), float(em_beta)
    elif subset == "m1+m2":
        d_tx, d_rr, beta_use = _zscore_mean(m1, m2), zeros_rr, 0.0
    elif subset == "m1+m3":
        d_tx, d_rr, beta_use = m1.astype(np.float32), m3.astype(np.float32), float(em_beta)
    elif subset == "m2+m3":
        d_tx, d_rr, beta_use = m2.astype(np.float32), m3.astype(np.float32), float(em_beta)
    else:  # "all"
        d_tx, d_rr, beta_use = _zscore_mean(m1, m2), m3.astype(np.float32), float(em_beta)

    # Defensive: EM requires non-negative distances. Floor any small
    # numerical drift below zero (can happen after z-scoring).
    if d_tx.min() < 0:
        d_tx = d_tx - d_tx.min()
    if d_rr.min() < 0:
        d_rr = d_rr - d_rr.min()
    return d_tx.astype(np.float32), d_rr.astype(np.float32), beta_use
