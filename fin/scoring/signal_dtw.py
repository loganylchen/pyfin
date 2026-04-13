"""Signal segment extraction and pairwise DTW computation.

Uses GPU-accelerated DTW when available, with scipy/numpy CPU fallback.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

from fin.scoring.eventalign_parser import ReadCandidateScore

logger = logging.getLogger(__name__)


def extract_signal_segments(
    scores: List[ReadCandidateScore],
    signal_reader,
    signal_format: str = "pod5",
) -> Dict[str, np.ndarray]:
    """Extract raw signal segments for each read using eventalign index bounds.

    For each read, uses the signal_start_idx and signal_end_idx from its
    best-scoring candidate to extract the relevant signal segment.

    Args:
        scores: Parsed ReadCandidateScore list.
        signal_reader: Opened Pod5Reader or Slow5Reader instance.
        signal_format: "pod5" or "slow5" to select the right API.

    Returns:
        Dict mapping read_name -> signal segment (np.ndarray of float).
    """
    # Find best-scoring candidate per read (highest log-likelihood)
    best_per_read: Dict[str, ReadCandidateScore] = {}
    for s in scores:
        existing = best_per_read.get(s.read_name)
        if existing is None or s.total_log_likelihood > existing.total_log_likelihood:
            best_per_read[s.read_name] = s

    segments: Dict[str, np.ndarray] = {}
    for read_name, best in best_per_read.items():
        start = best.signal_start_idx
        end = best.signal_end_idx

        if signal_format == "pod5":
            result = signal_reader.get_calibrated_signal(read_name)
        else:
            result = signal_reader.get_picoamp_signal(read_name)

        if result is None:
            logger.warning("No signal found for read %s", read_name)
            continue

        signal, _meta = result
        signal = np.asarray(signal, dtype=np.float32)

        # Clamp indices
        start = max(0, start)
        end = min(len(signal), end)

        if end <= start:
            logger.warning(
                "Invalid signal range [%d, %d) for read %s", start, end, read_name
            )
            continue

        segments[read_name] = signal[start:end]

    logger.info("Extracted signal segments for %d / %d reads", len(segments), len(best_per_read))
    return segments


def compute_read_to_read_dtw(
    segments: Dict[str, np.ndarray],
    read_ids: List[str],
    use_gpu: bool = True,
) -> np.ndarray:
    """Compute pairwise DTW distances between read signal segments.

    Uses GPU-accelerated DTW when available, falls back to scipy/numpy.

    Args:
        segments: Dict mapping read_name -> signal segment.
        read_ids: Ordered list of read IDs for matrix indexing.
        use_gpu: Whether to attempt GPU acceleration.

    Returns:
        np.ndarray of shape (n_reads, n_reads) with pairwise DTW distances.
        Reads without segments get distance 0 to all others.
    """
    n = len(read_ids)
    dist = np.zeros((n, n), dtype=np.float64)

    # Collect segments in order
    available = []
    available_indices = []
    for i, rid in enumerate(read_ids):
        if rid in segments:
            available.append(segments[rid])
            available_indices.append(i)

    if len(available) < 2:
        return dist

    # Try GPU batch DTW first (equal-length only)
    if use_gpu:
        try:
            from fin._dtw import dtw_pairwise, dtw_distance, is_available

            if is_available():
                lengths = [len(s) for s in available]
                if len(set(lengths)) == 1:
                    # Equal length: use efficient batch pairwise
                    batch = np.stack(available)
                    pairwise = dtw_pairwise(batch, use_open_end=True)
                    for ii, gi in enumerate(available_indices):
                        for jj, gj in enumerate(available_indices):
                            dist[gi, gj] = pairwise[ii, jj]
                    return dist
                else:
                    # Variable length: use GPU dtw_distance per pair
                    return _pairwise_gpu_variable(
                        available, available_indices, n
                    )
        except (ImportError, RuntimeError) as e:
            logger.info("GPU DTW unavailable, falling back to CPU: %s", e)

    # CPU fallback using scipy
    return _pairwise_cpu(available, available_indices, n)


def _pairwise_gpu_variable(
    segments: List[np.ndarray],
    indices: List[int],
    n: int,
) -> np.ndarray:
    """Compute pairwise DTW for variable-length segments using GPU."""
    from fin._dtw import dtw_distance

    dist = np.zeros((n, n), dtype=np.float64)
    for i in range(len(segments)):
        for j in range(i + 1, len(segments)):
            d = dtw_distance(segments[i], segments[j], use_open_end=True)
            dist[indices[i], indices[j]] = d
            dist[indices[j], indices[i]] = d
    return dist


def _pairwise_cpu(
    segments: List[np.ndarray],
    indices: List[int],
    n: int,
) -> np.ndarray:
    """Compute pairwise DTW using scipy or optimized numpy (CPU fallback)."""
    dist = np.zeros((n, n), dtype=np.float64)
    for i in range(len(segments)):
        for j in range(i + 1, len(segments)):
            d = _cpu_dtw(segments[i], segments[j])
            dist[indices[i], indices[j]] = d
            dist[indices[j], indices[i]] = d
    return dist


def _cpu_dtw(seq1: np.ndarray, seq2: np.ndarray) -> float:
    """DTW distance using scipy if available, else vectorized numpy.

    Open-end variant: takes minimum of last row.
    """
    try:
        from scipy.spatial.distance import cdist

        n, m = len(seq1), len(seq2)
        # Compute pairwise cost matrix
        cost_matrix = (seq1[:, None] - seq2[None, :]) ** 2

        # Accumulate DTW cost matrix row by row (vectorized per row)
        dtw_cost = np.full((n + 1, m + 1), np.inf)
        dtw_cost[0, 0] = 0.0

        for i in range(1, n + 1):
            prev_min = np.minimum(dtw_cost[i - 1, :-1], dtw_cost[i - 1, 1:])
            prev_min = np.minimum(prev_min, dtw_cost[i, :-1])
            dtw_cost[i, 1:] = cost_matrix[i - 1] + prev_min

        # Open-end: minimum of last row
        return float(np.min(dtw_cost[n, 1:]))
    except ImportError:
        # Pure numpy fallback (row-vectorized, not per-cell Python loop)
        return _numpy_dtw(seq1, seq2)


def _numpy_dtw(seq1: np.ndarray, seq2: np.ndarray) -> float:
    """Row-vectorized DTW in pure numpy (no scipy dependency)."""
    n, m = len(seq1), len(seq2)
    # Previous and current rows only (memory efficient)
    prev_row = np.full(m + 1, np.inf)
    prev_row[0] = 0.0

    for i in range(1, n + 1):
        curr_row = np.full(m + 1, np.inf)
        costs = (seq1[i - 1] - seq2) ** 2  # vectorized cost for row i

        for j in range(1, m + 1):
            curr_row[j] = costs[j - 1] + min(prev_row[j], prev_row[j - 1], curr_row[j - 1])

        prev_row = curr_row

    # Open-end: minimum of last row
    return float(np.min(prev_row[1:]))
