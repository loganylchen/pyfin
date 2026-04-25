"""Signal segment extraction and pairwise DTW computation.

Uses GPU-accelerated DTW when available, with scipy/numpy CPU fallback.
"""

from __future__ import annotations

import logging
from typing import Dict, List

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

    if use_gpu:
        try:
            from fin._dtw import is_available

            if is_available():
                return _pairwise_gpu_chunked(available, available_indices, n)
        except (ImportError, RuntimeError) as e:
            logger.info("GPU DTW unavailable: %s", e)

    return _pairwise_cpu(available, available_indices, n)


def _pairwise_gpu_chunked(
    segments: List[np.ndarray],
    indices: List[int],
    n: int,
) -> np.ndarray:
    """Compute pairwise DTW with memory-aware GPU chunking."""
    from fin._dtw import dtw_pairwise_varlen, estimate_gpu_memory, get_free_gpu_memory

    max_len = max(len(s) for s in segments)
    num_seqs = len(segments)

    estimated = estimate_gpu_memory(num_seqs, max_len)
    free_mem = get_free_gpu_memory()
    budget = int(free_mem * 0.8) if free_mem > 0 else 0

    if budget > 0 and estimated <= budget:
        pairwise = dtw_pairwise_varlen(segments, use_open_end=True)
        return _map_to_full(pairwise, indices, n)

    # Find max chunk size that fits GPU memory
    chunk_size = _find_chunk_size(max_len, budget) if budget > 0 else 32
    logger.info(
        "DTW chunking: %d seqs, max_len=%d, chunk_size=%d (est=%s, free=%s)",
        num_seqs, max_len, chunk_size,
        _fmt_bytes(estimated), _fmt_bytes(free_mem),
    )

    dist = np.zeros((n, n), dtype=np.float64)

    for i in range(0, num_seqs, chunk_size):
        chunk_a = list(range(i, min(i + chunk_size, num_seqs)))

        # Intra-chunk pairwise
        segs_a = [segments[k] for k in chunk_a]
        idx_a = [indices[k] for k in chunk_a]
        if len(segs_a) >= 2:
            pw = dtw_pairwise_varlen(segs_a, use_open_end=True)
            _fill_dist(dist, pw, idx_a)

        # Cross-chunk: pair chunk_a with each later chunk
        for j in range(i + chunk_size, num_seqs, chunk_size):
            chunk_b = list(range(j, min(j + chunk_size, num_seqs)))
            combined_segs = segs_a + [segments[k] for k in chunk_b]
            pw = dtw_pairwise_varlen(combined_segs, use_open_end=True)
            _fill_cross(dist, pw, idx_a, [indices[k] for k in chunk_b], len(segs_a))

    return dist


def _find_chunk_size(max_len: int, budget: int) -> int:
    """Binary search for largest num_seqs that fits in GPU memory budget."""
    from fin._dtw import estimate_gpu_memory

    lo, hi = 2, 10000
    while lo < hi:
        mid = (lo + hi + 1) // 2
        if estimate_gpu_memory(mid, max_len) <= budget:
            lo = mid
        else:
            hi = mid - 1
    return max(lo, 2)


def _map_to_full(pairwise: np.ndarray, indices: List[int], n: int) -> np.ndarray:
    """Map small pairwise matrix back to full n x n matrix."""
    dist = np.zeros((n, n), dtype=np.float64)
    for ii, gi in enumerate(indices):
        for jj, gj in enumerate(indices):
            dist[gi, gj] = pairwise[ii, jj]
    return dist


def _fill_dist(dist: np.ndarray, pw: np.ndarray, idx: List[int]) -> None:
    """Fill distance matrix from pairwise result."""
    for ii, gi in enumerate(idx):
        for jj, gj in enumerate(idx):
            dist[gi, gj] = pw[ii, jj]


def _fill_cross(
    dist: np.ndarray, pw: np.ndarray, idx_a: List[int], idx_b: List[int], len_a: int
) -> None:
    """Fill only cross-chunk distances (skip intra-chunk already filled)."""
    for ii, gi in enumerate(idx_a):
        for jj, gj in enumerate(idx_b):
            d = pw[ii, len_a + jj]
            dist[gi, gj] = d
            dist[gj, gi] = d


def _fmt_bytes(n: int) -> str:
    """Format byte count for logging."""
    if n >= 1 << 30:
        return f"{n / (1 << 30):.1f} GB"
    if n >= 1 << 20:
        return f"{n / (1 << 20):.1f} MB"
    return f"{n / (1 << 10):.1f} KB"


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
        import scipy.spatial.distance  # noqa: F401  # availability probe

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
