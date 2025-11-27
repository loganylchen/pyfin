"""
GPU-accelerated DTW (Dynamic Time Warping) using OpenDBA.

This module provides high-performance DTW distance computation using CUDA,
based on the OpenDBA implementation by Paul G. (https://github.com/nodrogluap/OpenDBA).

Features:
- Single-pair DTW distance computation
- Pairwise distance matrix for multiple sequences
- Open-start and open-end alignment options
- Automatic GPU/CPU fallback
"""

import numpy as np
from typing import List, Union, Optional
import warnings

# Try to import CUDA extension
try:
    from fin._opendba import opendba_cuda
    OPENDBA_AVAILABLE = True
except ImportError:
    warnings.warn(
        "OpenDBA CUDA extension not available. Install with: pip install -e .\n"
        "Will use CPU fallback for DTW calculations.",
        ImportWarning
    )
    OPENDBA_AVAILABLE = False


def dtw_distance(
    seq1: np.ndarray,
    seq2: np.ndarray,
    open_start: bool = False,
    open_end: bool = False,
    use_cuda: bool = True
) -> float:
    """
    Compute DTW distance between two sequences using GPU acceleration.

    Args:
        seq1: First sequence (1D numpy array)
        seq2: Second sequence (1D numpy array)
        open_start: Allow open start (skip beginning of sequences). Default: False
        open_end: Allow open end (skip ending of sequences). Default: False
        use_cuda: Use GPU acceleration if available. Default: True

    Returns:
        float: DTW distance between sequences

    Examples:
        >>> import numpy as np
        >>> seq1 = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)
        >>> seq2 = np.array([1.1, 2.1, 3.1, 4.1], dtype=np.float32)
        >>> dist = dtw_distance(seq1, seq2)
        >>> print(f"DTW distance: {dist:.2f}")

        With open start/end for time series with different lengths:
        >>> seq1 = np.array([0.0, 0.0, 1.0, 2.0, 3.0], dtype=np.float32)
        >>> seq2 = np.array([1.0, 2.0, 3.0, 0.0, 0.0], dtype=np.float32)
        >>> dist = dtw_distance(seq1, seq2, open_start=True, open_end=True)
    """
    # Convert to float32 if needed
    if seq1.dtype != np.float32:
        seq1 = seq1.astype(np.float32)
    if seq2.dtype != np.float32:
        seq2 = seq2.astype(np.float32)

    # Ensure contiguous arrays
    seq1 = np.ascontiguousarray(seq1)
    seq2 = np.ascontiguousarray(seq2)

    if OPENDBA_AVAILABLE:
        return opendba_cuda.dtw_distance(
            seq1, seq2,
            open_start=open_start,
            open_end=open_end,
            use_cuda=use_cuda
        )
    else:
        # CPU fallback - simple DTW implementation
        return _dtw_distance_cpu(seq1, seq2, open_start, open_end)


def dtw_pairwise_matrix(
    sequences: List[np.ndarray],
    open_start: bool = False,
    open_end: bool = False,
    use_cuda: bool = True
) -> np.ndarray:
    """
    Compute pairwise DTW distance matrix for list of sequences.

    This is highly optimized for batch processing using GPU parallelization.

    Args:
        sequences: List of sequences (each as 1D numpy array)
        open_start: Allow open start for alignments. Default: False
        open_end: Allow open end for alignments. Default: False
        use_cuda: Use GPU acceleration. Default: True

    Returns:
        numpy.ndarray: Square distance matrix (n_sequences x n_sequences)
                      with DTW distances. Diagonal is 0.

    Examples:
        >>> # Generate some test sequences
        >>> np.random.seed(42)
        >>> sequences = [np.random.randn(100).astype(np.float32) for _ in range(10)]
        >>> dist_matrix = dtw_pairwise_matrix(sequences)
        >>> print(f"Distance matrix shape: {dist_matrix.shape}")
        >>> print(f"First few distances:\n{dist_matrix[:3, :3]}")

        >>> # With open alignment options
        >>> sequences = [np.sin(np.linspace(0, 4*np.pi, 100)).astype(np.float32),
        ...              np.sin(np.linspace(0.5, 4.5*np.pi, 110)).astype(np.float32)]
        >>> dist_matrix = dtw_pairwise_matrix(sequences, open_start=True, open_end=True)
    """
    # Ensure all sequences are float32
    sequences = [np.asarray(seq, dtype=np.float32) for seq in sequences]

    if OPENDBA_AVAILABLE:
        return opendba_cuda.dtw_pairwise_matrix(
            sequences,
            open_start=open_start,
            open_end=open_end,
            use_cuda=use_cuda
        )
    else:
        # CPU fallback
        return _dtw_pairwise_matrix_cpu(sequences, open_start, open_end)


def dtw_nearest_neighbors(
    query: np.ndarray,
    references: List[np.ndarray],
    k: int = 1,
    open_start: bool = False,
    open_end: bool = False,
    use_cuda: bool = True
) -> List[tuple[int, float]]:
    """
    Find k nearest neighbors to query sequence using DTW distance.

    Args:
        query: Query sequence
        references: List of reference sequences
        k: Number of nearest neighbors to return. Default: 1
        open_start: Allow open start. Default: False
        open_end: Allow open end. Default: False
        use_cuda: Use GPU acceleration. Default: True

    Returns:
        List of (index, distance) tuples for k nearest neighbors,
        sorted by distance (closest first)
    """
    # Compute distances to all references
    distances = []

    for idx, ref in enumerate(references):
        dist = dtw_distance(
            query, ref,
            open_start=open_start,
            open_end=open_end,
            use_cuda=use_cuda
        )
        distances.append((idx, dist))

    # Sort by distance and return top k
    distances.sort(key=lambda x: x[1])

    return distances[:k]


def _dtw_distance_cpu(seq1: np.ndarray, seq2: np.ndarray,
                      open_start: bool, open_end: bool) -> float:
    """CPU fallback implementation of DTW distance."""
    len1, len2 = len(seq1), len(seq2)

    # Initialize DP matrix
    dp = np.full((len1 + 1, len2 + 1), np.inf, dtype=np.float32)
    dp[0, 0] = 0.0

    # Open start allows skipping beginning
    if open_start:
        dp[1:, 0] = 0.0
        dp[0, 1:] = 0.0

    # Fill DP matrix
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            cost = abs(seq1[i-1] - seq2[j-1])

            min_prev = min(dp[i-1, j-1], dp[i-1, j], dp[i, j-1])
            dp[i, j] = cost + min_prev

    # Normal DTW: return bottom-right
    result = dp[len1, len2]

    # Open end: find minimum in last row or column
    if open_end:
        min_end = min(np.min(dp[:, len2]), np.min(dp[len1, :]))
        result = min_end

    return float(result)


def _dtw_pairwise_matrix_cpu(sequences: List[np.ndarray],
                             open_start: bool, open_end: bool) -> np.ndarray:
    """CPU fallback for pairwise DTW distance matrix."""
    n = len(sequences)
    dist_matrix = np.zeros((n, n), dtype=np.float32)

    for i in range(n):
        for j in range(i + 1, n):
            dist = _dtw_distance_cpu(
                sequences[i], sequences[j],
                open_start, open_end
            )
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    return dist_matrix


class DTWBatchProcessor:
    """
    Batch processor for computing DTW distances efficiently.

    Maintains GPU context and settings for repeated DTW computations.
    """

    def __init__(self, use_cuda: bool = True, open_start: bool = False, open_end: bool = False):
        """
        Initialize batch processor.

        Args:
            use_cuda: Use GPU acceleration. Default: True
            open_start: Default open start setting. Default: False
            open_end: Default open end setting. Default: False
        """
        self.use_cuda = use_cuda and OPENDBA_AVAILABLE
        self.default_open_start = open_start
        self.default_open_end = open_end

        if not self.use_cuda:
            warnings.warn("CUDA not available, using CPU implementation")

    def distance(self, seq1: np.ndarray, seq2: np.ndarray,
                 open_start: Optional[bool] = None, open_end: Optional[bool] = None) -> float:
        """
        Compute DTW distance with processor settings.

        Args:
            seq1: First sequence
            seq2: Second sequence
            open_start: Override open_start setting
            open_end: Override open_end setting

        Returns:
            float: DTW distance
        """
        os = self.default_open_start if open_start is None else open_start
        oe = self.default_open_end if open_end is None else open_end

        return dtw_distance(seq1, seq2, open_start=os, open_end=oe, use_cuda=self.use_cuda)

    def pairwise_matrix(self, sequences: List[np.ndarray],
                        open_start: Optional[bool] = None,
                        open_end: Optional[bool] = None) -> np.ndarray:
        """
        Compute pairwise distance matrix with processor settings.

        Args:
            sequences: List of sequences
            open_start: Override open_start setting
            open_end: Override open_end setting

        Returns:
            numpy.ndarray: Distance matrix
        """
        os = self.default_open_start if open_start is None else open_start
        oe = self.default_open_end if open_end is None else open_end

        return dtw_pairwise_matrix(sequences, open_start=os, open_end=oe, use_cuda=self.use_cuda)

    def cluster(self, sequences: List[np.ndarray], threshold: float = 0.1) -> List[List[int]]:
        """
        Simple clustering of sequences based on DTW distances.

        Args:
            sequences: List of sequences
            threshold: Distance threshold for clustering

        Returns:
            List of clusters (each cluster is list of indices)
        """
        dist_matrix = self.pairwise_matrix(sequences)

        # Simple hierarchical clustering with threshold
        n = len(sequences)
        clusters = [[i] for i in range(n)]

        # Find pairs to merge
        while True:
            min_dist = float('inf')
            merge_i, merge_j = -1, -1

            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    # Find minimum distance between clusters
                    min_cluster_dist = float('inf')
                    for idx_i in clusters[i]:
                        for idx_j in clusters[j]:
                            if dist_matrix[idx_i, idx_j] < min_cluster_dist:
                                min_cluster_dist = dist_matrix[idx_i, idx_j]

                    if min_cluster_dist < min_dist:
                        min_dist = min_cluster_dist
                        merge_i, merge_j = i, j

            if min_dist <= threshold and merge_i != -1:
                # Merge clusters
                clusters[merge_i].extend(clusters[merge_j])
                clusters.pop(merge_j)
            else:
                break

        return clusters


# Convenience functions

def dtw_similarity(seq1: np.ndarray, seq2: np.ndarray,
                     open_start: bool = False, open_end: bool = False,
                     use_cuda: bool = True) -> float:
    """
    Compute DTW similarity (1 / (1 + distance)) between two sequences.

    This is useful for applications where you need a similarity score
    in range (0, 1] instead of a distance.

    Args:
        seq1: First sequence
        seq2: Second sequence
        open_start: Allow open start
        open_end: Allow open end
        use_cuda: Use GPU acceleration

    Returns:
        float: Similarity score between 0 and 1

    Examples:
        >>> seq1 = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)
        >>> seq2 = np.array([1.1, 2.1, 3.1, 4.1], dtype=np.float32)
        >>> sim = dtw_similarity(seq1, seq2)
        >>> print(f"Similarity: {sim:.3f}")
    """
    dist = dtw_distance(seq1, seq2, open_start, open_end, use_cuda)
    return 1.0 / (1.0 + dist)
