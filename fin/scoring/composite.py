"""Composite scoring functions for transcript candidates.

Computes coherence, discrimination, and combined scores from EM soft
assignment matrices and distance matrices.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import List, Optional

from fin.candidates.dataclasses import TranscriptCandidate

# Optional CuPy support — same pattern as fin/analysis/assignments.py
try:
    import cupy as cp  # type: ignore

    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False


def _get_array_module(use_gpu: bool = True):
    """Return (xp, to_device, to_host) tuple for numpy or cupy."""
    if use_gpu and CUPY_AVAILABLE:
        return cp, cp.asarray, cp.asnumpy
    return np, lambda x: x, lambda x: x


# ---------------------------------------------------------------------------
# Dataclass
# ---------------------------------------------------------------------------


@dataclass
class CompositeScore:
    """Composite scoring result for a single transcript candidate."""

    candidate_id: str
    coherence: float
    discrimination: float
    combined: float


# ---------------------------------------------------------------------------
# Individual scoring functions
# ---------------------------------------------------------------------------


def compute_coherence_score(
    dist_read_to_read: np.ndarray,
    R: np.ndarray,
    candidate_idx: int,
) -> float:
    """Measure how similar the reads supporting a candidate are to each other.

    Args:
        dist_read_to_read: (n_reads, n_reads) pairwise read distance matrix.
        R: (n_reads, n_tx) soft assignment matrix from EM.
        candidate_idx: Column index for the candidate of interest.

    Returns:
        Coherence score in [0, 1].  Higher means tighter cluster.
    """
    w_i = R[:, candidate_idx].astype(np.float64)  # shape (n_reads,)
    w_sum = float(w_i.sum())

    if w_sum < 1e-10:
        return 0.0

    # Weighted mean pairwise distance: w^T D w / (w_sum^2)
    d_mean = float(w_i @ dist_read_to_read @ w_i) / (w_sum * w_sum + 1e-10)

    # Reference scale: median of the full distance matrix
    d_ref = float(max(np.median(dist_read_to_read), 1e-6))

    coherence = float(np.exp(-d_mean / d_ref))
    coherence = float(np.clip(coherence, 0.0, 1.0))
    return coherence


def compute_discrimination_score(
    dist_read_to_tx: np.ndarray,
    candidate_idx: int,
    tau: Optional[float] = None,
) -> float:
    """Measure how clearly a candidate beats its alternatives for each read.

    Args:
        dist_read_to_tx: (n_reads, n_tx) distance matrix.
        candidate_idx: Column index for the candidate of interest.
        tau: Temperature for sigmoid smoothing.  Defaults to std(delta).

    Returns:
        Discrimination score in [0, 1].  0.5 = tied; >0.5 = this candidate wins.
    """
    n_reads, n_tx = dist_read_to_tx.shape

    # Guard: single candidate — no alternative exists
    if n_tx == 1:
        return 0.5

    d_self = dist_read_to_tx[:, candidate_idx]  # (n_reads,)

    # next-best distance for each read: min over all other candidates
    other_cols = [j for j in range(n_tx) if j != candidate_idx]
    d_others = dist_read_to_tx[:, other_cols]  # (n_reads, n_tx-1)
    d_next_best = d_others.min(axis=1)  # (n_reads,)

    # delta_i > 0 means this candidate is closer (better) than the next best
    delta = (d_next_best - d_self).astype(np.float64)

    if tau is None:
        tau = float(max(float(np.std(delta)), 1e-6))

    def _sigmoid(x: np.ndarray) -> np.ndarray:
        return 1.0 / (1.0 + np.exp(-x))

    discrimination = float(np.mean(_sigmoid(delta / tau)))
    discrimination = float(np.clip(discrimination, 0.0, 1.0))
    return discrimination


def compute_combined_score(
    coherence: float,
    discrimination: float,
    alpha: float = 0.5,
) -> float:
    """Geometric mean of coherence and discrimination.

    Args:
        coherence: Coherence score in [0, 1].
        discrimination: Discrimination score in [0, 1].
        alpha: Weight for coherence (1-alpha for discrimination).

    Returns:
        Combined score in [0, 1].
    """
    coh = float(np.clip(coherence, 0.0, 1.0))
    dis = float(np.clip(discrimination, 0.0, 1.0))
    combined = float((coh ** alpha) * (dis ** (1.0 - alpha)))
    return float(np.clip(combined, 0.0, 1.0))


# ---------------------------------------------------------------------------
# Batch scoring
# ---------------------------------------------------------------------------


def score_candidates_composite(
    candidates: List[TranscriptCandidate],
    dist_read_to_tx: np.ndarray,
    dist_read_to_read: np.ndarray,
    R: np.ndarray,
    alpha: float = 0.5,
    tau: Optional[float] = None,
    use_gpu: bool = False,
) -> List[CompositeScore]:
    """Score all candidates with coherence, discrimination, and combined metrics.

    Args:
        candidates: Ordered list of TranscriptCandidate objects.
        dist_read_to_tx: (n_reads, n_tx) distance matrix.
        dist_read_to_read: (n_reads, n_reads) pairwise read distances.
        R: (n_reads, n_tx) soft assignment matrix from EM.
        alpha: Coherence weight in combined score.
        tau: Optional temperature for discrimination sigmoid.
        use_gpu: Hint to use GPU when CuPy is available (falls back to numpy).

    Returns:
        List of CompositeScore, one per candidate, aligned with input order.
    """
    # use_gpu is accepted as a forward-compatibility hint; numpy impl used
    # unless cupy ops are explicitly desired in future extension.
    _ = use_gpu  # acknowledged; numpy fallback always safe

    scores: List[CompositeScore] = []
    for j, cand in enumerate(candidates):
        coh = compute_coherence_score(dist_read_to_read, R, j)
        dis = compute_discrimination_score(dist_read_to_tx, j, tau=tau)
        com = compute_combined_score(coh, dis, alpha=alpha)
        scores.append(
            CompositeScore(
                candidate_id=cand.candidate_id,
                coherence=coh,
                discrimination=dis,
                combined=com,
            )
        )
    return scores


# ---------------------------------------------------------------------------
# Shared helpers for pipeline runners
# ---------------------------------------------------------------------------


def derive_prior_weights(
    combined_scores: np.ndarray,
    n_tx: int,
    cap: float,
) -> Optional[np.ndarray]:
    """Derive EM prior weights from composite combined scores.

    Normalizes combined scores to a probability vector, clips the per-candidate
    boost relative to a uniform prior to [1/cap, cap], then re-normalizes so
    the result sums to 1.

    Args:
        combined_scores: Array of combined scores, one per candidate.
        n_tx: Number of candidates (length of combined_scores).
        cap: Maximum multiplicative boost/penalty relative to uniform prior.

    Returns:
        Prior weight vector of shape (n_tx,), or None if n_tx == 0.
    """
    if n_tx <= 0:
        return None
    uniform_weight = 1.0 / n_tx
    cap = max(float(cap), 1.0)
    raw = np.asarray(combined_scores, dtype=float) + 1e-12
    normalized = raw / raw.sum()
    ratio = normalized / uniform_weight
    ratio = np.clip(ratio, 1.0 / cap, cap)
    shaped = ratio * uniform_weight
    return shaped / shaped.sum()


def populate_quant_scores(
    quant_results,
    composite_scores: List[CompositeScore],
) -> None:
    """Populate coherence/discrimination/combined fields on QuantResult objects.

    Args:
        quant_results: Iterable of QuantResult (mutated in place).
        composite_scores: Composite scores from score_candidates_composite.
    """
    score_by_cid = {s.candidate_id: s for s in composite_scores}
    for qr in quant_results:
        s = score_by_cid.get(qr.candidate_id)
        if s is not None:
            qr.coherence_score = s.coherence
            qr.discrimination_score = s.discrimination
            qr.combined_score = s.combined


def subsample_reads_for_dtw(
    read_ids: List[str],
    max_dtw: Optional[int],
) -> List[str]:
    """Uniformly subsample read_ids when count exceeds max_dtw.

    Args:
        read_ids: Sorted list of read IDs.
        max_dtw: Maximum number of reads to keep for DTW; None/0 disables.

    Returns:
        Subsampled list (or the original list when no subsampling is needed).
    """
    if not max_dtw or len(read_ids) <= max_dtw:
        return read_ids
    idx = np.linspace(0, len(read_ids) - 1, max_dtw).astype(int)
    return [read_ids[i] for i in idx]
