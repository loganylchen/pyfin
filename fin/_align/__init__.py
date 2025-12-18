"""
F5C original source code - Event-to-sequence alignment for nanopore data

This module provides Python wrappers for the original f5c eventalign implementation.
Both CPU and GPU (CUDA) versions are supported if available.
"""

import numpy as np

# Try to import CPU version
try:
    from ._align import eventalign as _eventalign_cpu
    from ._align import profile_hmm_eventalign as _profile_hmm_eventalign_cpu

    _ALIGN_CPU_AVAILABLE = True
except ImportError as e:
    _ALIGN_CPU_AVAILABLE = False
    _align_cpu_import_error = str(e)

# Try to import CUDA version
try:
    from ._align_cuda import eventalign as _eventalign_cuda
    from ._align_cuda import profile_hmm_eventalign as _profile_hmm_eventalign_cuda

    _ALIGN_CUDA_AVAILABLE = True
except ImportError as e:
    _ALIGN_CUDA_AVAILABLE = False
    _align_cuda_import_error = str(e)

# Log backend availability
if _ALIGN_CUDA_AVAILABLE:
    print("🚀 F5C Eventalign: GPU (CUDA) acceleration ENABLED")
    _default_backend = "cuda"
elif _ALIGN_CPU_AVAILABLE:
    print("💻 F5C Eventalign: Using CPU implementation")
    _default_backend = "cpu"
else:
    print("⚠️  F5C Eventalign: No backend available - please rebuild the package")
    _default_backend = None


def eventalign(
    raw_signal: np.ndarray,
    sequence: str,
    kmer_size: int = 5,
    model: dict = None,
    use_gpu: bool = None,
) -> dict:
    """
    Align detected RNA events to a reference sequence using f5c algorithm.

    This uses the original f5c source code for event-to-sequence alignment.
    Events are automatically reversed internally to match RNA 5'->3' direction.

    Args:
        raw_signal: 1D numpy array of raw nanopore signal (float32).
        sequence: Reference RNA sequence in 5'->3' direction.
        kmer_size: K-mer size (5 or 9, default: 5).
        model: Optional k-mer model dict (uses built-in RNA model if None).
        use_gpu: Force GPU (True) or CPU (False). Auto-detect if None.

    Returns:
        dict with keys:
            - base_to_event_map: List of dicts per k-mer with 'start', 'stop', 'kmer'
            - scaling: Dict with 'scale' and 'shift'
            - n_events: Total number of events
            - n_aligned_pairs: Number of alignments

    Raises:
        RuntimeError: If no backend available or alignment fails.
        TypeError: If inputs have incorrect types.
        ValueError: If sequence is too short.

    Example:
        >>> signal = np.random.randn(10000).astype(np.float32)
        >>> seq = "ACGUACGUACGU"
        >>> result = eventalign(signal, seq, kmer_size=5)
    """
    # Auto-detect backend
    if use_gpu is None:
        use_gpu = _default_backend == "cuda"

    # Select backend
    if use_gpu:
        if not _ALIGN_CUDA_AVAILABLE:
            raise RuntimeError(
                f"GPU backend requested but not available. " f"Error: {_align_cuda_import_error}"
            )
        backend_func = _eventalign_cuda
    else:
        if not _ALIGN_CPU_AVAILABLE:
            raise RuntimeError(f"CPU backend not available. " f"Error: {_align_cpu_import_error}")
        backend_func = _eventalign_cpu

    # Validate inputs
    if not isinstance(raw_signal, np.ndarray):
        raise TypeError(f"raw_signal must be numpy array (got {type(raw_signal)})")
    if raw_signal.ndim != 1:
        raise TypeError(f"raw_signal must be 1D (got {raw_signal.ndim}D)")
    if raw_signal.dtype != np.float32:
        raw_signal = raw_signal.astype(np.float32)
        print("Warning: raw_signal converted to float32")

    if not isinstance(sequence, str):
        raise TypeError(f"sequence must be string (got {type(sequence)})")
    if len(sequence) < kmer_size:
        raise ValueError(f"sequence length ({len(sequence)}) < kmer_size ({kmer_size})")

    # Call backend (RNA-only)
    return backend_func(raw_signal, sequence, model, kmer_size)


def profile_hmm_eventalign(
    raw_signal: np.ndarray,
    sequence: str,
    kmer_size: int = 5,
    events_per_base: float = 3.0,
    use_gpu: bool = None,
) -> dict:
    """
    Full f5c Profile HMM eventalign with detailed alignment output.

    This is the complete f5c eventalign using Viterbi HMM, matching
    the original f5c command-line tool output format.

    Args:
        raw_signal: 1D numpy array of raw signal (float32).
        sequence: Reference RNA sequence in 5'->3' direction.
        kmer_size: K-mer size (5 or 9, default: 5).
        events_per_base: Expected events per base (default: 3.0).
        use_gpu: Force GPU (True) or CPU (False). Auto-detect if None.

    Returns:
        dict with keys:
            - alignment: List of dicts with full event_alignment_t data:
                * ref_position: Reference position (0-based)
                * ref_kmer: Reference k-mer string
                * event_idx: Event index (-1 for kmer skips)
                * hmm_state: 'M' (match), 'K' (kmer_skip), 'B' (bad_event)
                * event_mean, event_stdv, event_duration: Event stats
                * model_mean, model_stdv: Model parameters
                * scaled_model_mean, scaled_model_stdv: Scaled parameters
            - scaling: Dict with 'scale' and 'shift'
            - n_events: Total events detected
            - n_aligned: Number of alignment records
            - events_per_base: Ratio used

    Raises:
        RuntimeError: If backend not available.

    Example:
        >>> signal = np.random.randn(10000).astype(np.float32)
        >>> seq = "ACGUACGUACGU"
        >>> result = profile_hmm_eventalign(signal, seq)
        >>> print(f"Found {result['n_aligned']} alignments")
    """
    # Auto-detect backend
    if use_gpu is None:
        use_gpu = _default_backend == "cuda"

    # Select backend
    if use_gpu:
        if not _ALIGN_CUDA_AVAILABLE:
            raise RuntimeError(
                f"GPU backend requested but not available. " f"Error: {_align_cuda_import_error}"
            )
        backend_func = _profile_hmm_eventalign_cuda
    else:
        if not _ALIGN_CPU_AVAILABLE:
            raise RuntimeError(f"CPU backend not available. " f"Error: {_align_cpu_import_error}")
        backend_func = _profile_hmm_eventalign_cpu

    # Validate inputs
    if not isinstance(raw_signal, np.ndarray):
        raise TypeError(f"raw_signal must be numpy array")
    if raw_signal.ndim != 1:
        raise TypeError(f"raw_signal must be 1D")
    if raw_signal.dtype != np.float32:
        raw_signal = raw_signal.astype(np.float32)

    if not isinstance(sequence, str):
        raise TypeError(f"sequence must be string")
    if len(sequence) < kmer_size:
        raise ValueError(f"sequence too short")

    # Call backend (RNA-only)
    return backend_func(raw_signal, sequence, kmer_size, events_per_base)


__all__ = ["eventalign", "profile_hmm_eventalign"]
