"""
CUDA-accelerated Dynamic Time Warping (DTW) module

This module provides GPU-accelerated DTW distance calculation.
"""

import subprocess
from typing import Union

import numpy as np

# Try to import the CUDA extension
try:
    from ._cuda_dtw import dtw_distance as _dtw_distance_cuda
    from ._cuda_dtw import dtw_pairwise as _dtw_pairwise_cuda
    from ._cuda_dtw import dtw_pairwise_varlen as _dtw_pairwise_varlen_cuda
    from ._cuda_dtw import get_free_gpu_memory as _get_free_gpu_memory_cuda
    from ._cuda_dtw import cleanup as _cuda_cleanup

    CUDA_AVAILABLE = True
except ImportError as e:
    CUDA_AVAILABLE = False
    _import_error = str(e)

# Log backend availability upon module import
if CUDA_AVAILABLE:
    print("🚀 DTW: GPU (CUDA) acceleration ENABLED")
else:
    print("💻 DTW: CPU implementation (CUDA not available)")


def dtw_distance(
    seq1: Union[np.ndarray, list],
    seq2: Union[np.ndarray, list],
    use_open_start: bool = False,
    use_open_end: bool = False,
) -> float:
    """
    Compute DTW distance between two sequences using CUDA acceleration.

    Parameters
    ----------
    seq1 : array-like
        First sequence (will be converted to float32 numpy array)
    seq2 : array-like
        Second sequence (will be converted to float32 numpy array)
    use_open_start : bool, optional
        Enable open start boundary condition (default: False)
    use_open_end : bool, optional
        Enable open end boundary condition (default: False)

    Returns
    -------
    float
        DTW distance between seq1 and seq2

    Raises
    ------
    RuntimeError
        If CUDA extension is not available
    ValueError
        If input sequences are invalid

    Examples
    --------
    >>> import numpy as np
    >>> from fin._dtw import dtw_distance
    >>> seq1 = np.random.randn(100).astype(np.float32)
    >>> seq2 = np.random.randn(100).astype(np.float32)
    >>> distance = dtw_distance(seq1, seq2)
    >>> print(f"DTW distance: {distance}")
    """
    if not CUDA_AVAILABLE:
        raise RuntimeError(
            f"CUDA DTW extension is not available.\n"
            f"Import error: {_import_error}\n\n"
            f"The extension was not built during installation. This can happen if:\n"
            f"  1. CUDA Toolkit is not installed (check with: nvcc --version)\n"
            f"  2. The build process skipped the CUDA extension\n"
            f"  3. The package was installed in a different environment\n\n"
            f"To build with CUDA support:\n"
            f"  1. Install CUDA Toolkit from NVIDIA\n"
            f"  2. Ensure nvcc is in PATH\n"
            f"  3. Reinstall: pip uninstall py-fin && pip install -e .\n\n"
            f"Check availability with: fin._dtw.is_available()"
        )

    # Convert inputs to numpy arrays if needed
    if not isinstance(seq1, np.ndarray):
        seq1 = np.array(seq1, dtype=np.float32)
    else:
        seq1 = np.asarray(seq1, dtype=np.float32)

    if not isinstance(seq2, np.ndarray):
        seq2 = np.array(seq2, dtype=np.float32)
    else:
        seq2 = np.asarray(seq2, dtype=np.float32)

    # Ensure arrays are contiguous
    if not seq1.flags["C_CONTIGUOUS"]:
        seq1 = np.ascontiguousarray(seq1)
    if not seq2.flags["C_CONTIGUOUS"]:
        seq2 = np.ascontiguousarray(seq2)

    # Validate shapes
    if seq1.ndim != 1:
        raise ValueError(f"seq1 must be 1-dimensional, got shape {seq1.shape}")
    if seq2.ndim != 1:
        raise ValueError(f"seq2 must be 1-dimensional, got shape {seq2.shape}")

    if len(seq1) == 0 or len(seq2) == 0:
        raise ValueError("Sequences cannot be empty")

    # Call the CUDA function
    return _dtw_distance_cuda(
        seq1, seq2, use_open_start=int(use_open_start), use_open_end=int(use_open_end)
    )


def dtw_pairwise(
    sequences: Union[np.ndarray, list], use_open_start: bool = False, use_open_end: bool = False
) -> np.ndarray:
    """
    Compute pairwise DTW distances for a batch of sequences using CUDA.

    This is significantly more efficient than calling dtw_distance() in a loop,
    as it:
    - Transfers all sequences to GPU in one batch
    - Computes multiple DTW pairs in parallel
    - Amortizes memory allocation/deallocation overhead

    Parameters
    ----------
    sequences : array-like
        2D array of sequences with shape (num_sequences, seq_length)
        All sequences must have the same length
        Will be converted to float32 if needed
    use_open_start : bool, optional
        Enable open start boundary condition (default: False)
    use_open_end : bool, optional
        Enable open end boundary condition (default: False)

    Returns
    -------
    np.ndarray
        Distance matrix of shape (num_sequences, num_sequences)
        Matrix is symmetric with zeros on the diagonal

    Raises
    ------
    RuntimeError
        If CUDA extension is not available
    ValueError
        If input sequences are invalid

    Examples
    --------
    >>> import numpy as np
    >>> from fin._dtw import dtw_pairwise
    >>> # Generate 10 sequences of length 100
    >>> sequences = np.random.randn(10, 100).astype(np.float32)
    >>> distance_matrix = dtw_pairwise(sequences)
    >>> print(f"Distance matrix shape: {distance_matrix.shape}")
    >>> # distance_matrix[i, j] is the DTW distance between sequences[i] and sequences[j]
    """
    if not CUDA_AVAILABLE:
        raise RuntimeError(
            f"CUDA DTW extension is not available.\n"
            f"Import error: {_import_error}\n\n"
            f"Check availability with: fin._dtw.is_available()"
        )

    # Convert to numpy array if needed
    if not isinstance(sequences, np.ndarray):
        sequences = np.array(sequences, dtype=np.float32)
    else:
        sequences = np.asarray(sequences, dtype=np.float32)

    # Validate input
    if sequences.ndim != 2:
        raise ValueError(f"sequences must be 2D array, got shape {sequences.shape}")

    if sequences.shape[0] < 2:
        raise ValueError(f"Need at least 2 sequences, got {sequences.shape[0]}")

    if sequences.shape[1] == 0:
        raise ValueError("Sequence length cannot be 0")

    # GPU memory guard: estimate required memory and check availability
    num_seq, seq_len = sequences.shape
    # Memory estimate: cost matrix = seq_len * (num_seq-1) * 4 bytes + sequences + distance matrix
    estimated_bytes = (seq_len * (num_seq - 1) * 4) + (num_seq * seq_len * 4) + (num_seq * num_seq * 4)
    try:
        import subprocess as _sp
        _mem_result = _sp.run(
            ["nvidia-smi", "--query-gpu=memory.free", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=5,
        )
        if _mem_result.returncode == 0:
            free_mb = int(_mem_result.stdout.strip().split("\n")[0])
            free_bytes = free_mb * 1024 * 1024
            if estimated_bytes > free_bytes * 0.9:  # 90% safety margin
                raise MemoryError(
                    f"Insufficient GPU memory for DTW pairwise batch: "
                    f"estimated {estimated_bytes / 1e9:.1f} GB needed, "
                    f"{free_bytes / 1e9:.1f} GB available. "
                    f"Reduce batch size or use dtw_distance() in a loop."
                )
    except (FileNotFoundError, subprocess.TimeoutExpired, ValueError):
        pass  # Can't check memory, proceed and let CUDA handle it

    # Call the CUDA function
    return _dtw_pairwise_cuda(
        sequences, use_open_start=int(use_open_start), use_open_end=int(use_open_end)
    )


def cleanup():
    """
    Reset CUDA device and free all GPU resources.

    This function should be called when you're done using the CUDA DTW
    functionality to properly clean up GPU resources.

    Examples
    --------
    >>> from fin._dtw import cleanup
    >>> # After computing many DTW distances...
    >>> cleanup()
    """
    if not CUDA_AVAILABLE:
        return  # Nothing to cleanup if CUDA is not available

    _cuda_cleanup()


def is_available() -> bool:
    """
    Check if CUDA DTW extension is available.

    Returns
    -------
    bool
        True if CUDA extension is available, False otherwise
    """
    return CUDA_AVAILABLE


def dtw_pairwise_varlen(
    segments,
    use_open_start: bool = False,
    use_open_end: bool = False,
) -> np.ndarray:
    """Pairwise DTW for variable-length sequences (single GPU batch call).

    Accepts a list of variable-length 1D arrays. Handles padding internally.

    Parameters
    ----------
    segments : list of array-like
        Variable-length 1D sequences.
    use_open_start : bool, optional
        Enable open start boundary (default: False).
    use_open_end : bool, optional
        Enable open end boundary (default: False).

    Returns
    -------
    np.ndarray
        Distance matrix (num_sequences, num_sequences).
    """
    if not CUDA_AVAILABLE:
        raise RuntimeError(
            f"CUDA DTW extension is not available.\n"
            f"Import error: {_import_error}\n"
            f"Check availability with: fin._dtw.is_available()"
        )

    lengths = np.array([len(s) for s in segments], dtype=np.int64)
    max_len = int(lengths.max())

    padded = np.zeros((len(segments), max_len), dtype=np.float32)
    for i, s in enumerate(segments):
        s = np.asarray(s, dtype=np.float32)
        padded[i, : len(s)] = s

    return _dtw_pairwise_varlen_cuda(
        padded, lengths,
        use_open_start=int(use_open_start),
        use_open_end=int(use_open_end),
    )


def estimate_gpu_memory(num_sequences: int, max_length: int) -> int:
    """Estimate GPU bytes needed for pairwise varlen DTW.

    Parameters
    ----------
    num_sequences : int
        Number of sequences.
    max_length : int
        Maximum sequence length (padding width).

    Returns
    -------
    int
        Estimated GPU memory in bytes (with 20% headroom).
    """
    input_bytes = num_sequences * max_length * 4
    lengths_bytes = num_sequences * 8
    num_pairs = num_sequences * (num_sequences - 1) // 2
    pairs_bytes = num_pairs * 4
    max_parallel = num_sequences - 1
    cost_bytes = max_length * max_parallel * 4 * 2
    return int((input_bytes + lengths_bytes + pairs_bytes + cost_bytes) * 1.2)


def get_free_gpu_memory() -> int:
    """Query free GPU memory in bytes.

    Returns
    -------
    int
        Free GPU memory in bytes, or 0 if CUDA unavailable.
    """
    if not CUDA_AVAILABLE:
        return 0
    return _get_free_gpu_memory_cuda()


__all__ = [
    "dtw_distance",
    "dtw_pairwise",
    "dtw_pairwise_varlen",
    "estimate_gpu_memory",
    "get_free_gpu_memory",
    "cleanup",
    "is_available",
    "CUDA_AVAILABLE",
]
