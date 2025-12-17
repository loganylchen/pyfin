import numpy as np

# Try to import the compiled C extension
try:
    from ._event import get_events

    _EVENT_AVAILABLE = True
except ImportError as e:
    _EVENT_AVAILABLE = False
    _import_error = str(e)

# Try to import eventalign extension (CPU)
try:
    from ._eventalign import eventalign as _eventalign_c
    from ._eventalign import profile_hmm_eventalign as _profile_hmm_eventalign_c

    _EVENTALIGN_AVAILABLE = True
except ImportError as e:
    _EVENTALIGN_AVAILABLE = False
    _eventalign_import_error = str(e)

# Try to import CUDA eventalign extension (GPU)
try:
    from ._eventalign_cuda import eventalign as _eventalign_cuda
    from ._eventalign_cuda import profile_hmm_eventalign as _profile_hmm_eventalign_cuda

    _EVENTALIGN_CUDA_AVAILABLE = True
except ImportError as e:
    _EVENTALIGN_CUDA_AVAILABLE = False
    _eventalign_cuda_import_error = str(e)

# Log backend availability upon module import
if _EVENTALIGN_AVAILABLE:
    if _EVENTALIGN_CUDA_AVAILABLE:
        print("🚀 Eventalign: GPU (CUDA) acceleration ENABLED")
    else:
        print("💻 Eventalign: Using CPU implementation (CUDA not available)")
elif _EVENT_AVAILABLE:
    print("ℹ️  Event detection available (eventalign not built)")


def detect_events(raw_signal: np.ndarray) -> list[dict]:
    """
    High-level wrapper for nanopore RNA event detection (C-backed).

    Returns events in RAW SIGNAL order (no reversal):
    - events[0] = first event in raw signal (corresponds to 3' end of RNA)
    - events[n-1] = last event in raw signal (corresponds to 5' end of RNA)

    This matches the order you would see when iterating through the raw signal.
    The event_idx values in eventalign results use this same indexing.

    Args:
        raw_signal: 1D numpy array of raw nanopore signal values (must be float32).

    Returns:
        List of event dicts (in raw signal order), each with:
            - mean: Event mean current (float)
            - stdv: Event standard deviation (float)
            - start: Start index of the event in raw signal (int)
            - length: Length of the event (float)

    Raises:
        TypeError: If raw_signal is not a 1D float32 numpy array.
        RuntimeError: If event detection fails (no events found) or extension not available.
    """
    if not _EVENT_AVAILABLE:
        raise RuntimeError(
            f"f5c event detection extension is not available. "
            f"Import error: {_import_error}\n"
            f"Please ensure the package was built correctly with: pip install -e ."
        )

    # Validate input (add user-friendly checks)
    if not isinstance(raw_signal, np.ndarray):
        raise TypeError(f"raw_signal must be a numpy array (got {type(raw_signal)})")
    if raw_signal.ndim != 1:
        raise TypeError(f"raw_signal must be 1D (got {raw_signal.ndim}D)")
    if raw_signal.dtype != np.float32:
        # Convert to float32 if needed (common user mistake)
        raw_signal = raw_signal.astype(np.float32)
        print("Warning: raw_signal converted to float32 (required for C backend)")

    # Call C function (RNA mode with reversed events)
    return get_events(raw_signal)


def eventalign(
    raw_signal: np.ndarray,
    sequence: str,
    kmer_size: int = 5,
    model: dict = None,
) -> dict:
    """
    Align detected RNA events to a reference sequence.

    RNA-only: Events are automatically reversed internally to match the
    5'->3' sequence direction. You do NOT need to reverse your sequence.

    This function performs event-to-sequence alignment, creating a mapping
    between k-mers in the sequence and events in the raw signal.

    Args:
        raw_signal: 1D numpy array of raw nanopore signal (float32).
        sequence: Reference RNA sequence string in standard 5'->3' direction.
                  Do NOT reverse the sequence - event reversal is handled internally.
        kmer_size: Size of k-mers for alignment (5 or 9, default: 5).
        model: Optional k-mer model dictionary (default: uses built-in RNA model).

    Returns:
        dict with keys:
            - base_to_event_map: List of dicts, one per k-mer, with:
                * 'start': First event index for this k-mer
                * 'stop': Last event index for this k-mer
                * 'kmer': The k-mer sequence string
            - scaling: Dict with 'scale' and 'shift' calibration parameters
            - n_events: Total number of detected events
            - n_aligned_pairs: Number of event-to-kmer alignments

    Raises:
        RuntimeError: If eventalign extension is not available or alignment fails.
        TypeError: If inputs have incorrect types.
        ValueError: If sequence is too short for kmer_size.

    Example:
        >>> signal = np.random.randn(10000).astype(np.float32)
        >>> sequence = "ACGUACGUACGU"  # Use original 5'->3' sequence
        >>> result = eventalign(signal, sequence, kmer_size=5)
        >>> print(f"Detected {result['n_events']} events")
    """
    if not _EVENTALIGN_AVAILABLE:
        raise RuntimeError(
            f"eventalign extension is not available. "
            f"Import error: {_eventalign_import_error}\n"
            f"Please ensure the package was built correctly with: pip install -e ."
        )

    # Validate inputs
    if not isinstance(raw_signal, np.ndarray):
        raise TypeError(f"raw_signal must be a numpy array (got {type(raw_signal)})")
    if raw_signal.ndim != 1:
        raise TypeError(f"raw_signal must be 1D (got {raw_signal.ndim}D)")
    if raw_signal.dtype != np.float32:
        raw_signal = raw_signal.astype(np.float32)
        print("Warning: raw_signal converted to float32")

    if not isinstance(sequence, str):
        raise TypeError(f"sequence must be a string (got {type(sequence)})")

    if len(sequence) < kmer_size:
        raise ValueError(f"sequence length ({len(sequence)}) must be >= kmer_size ({kmer_size})")

    # Call C extension (RNA-only)
    return _eventalign_c(raw_signal, sequence, model, kmer_size)


def profile_hmm_eventalign(
    raw_signal: np.ndarray,
    sequence: str,
    kmer_size: int = 5,
    events_per_base: float = 3.0,
) -> dict:
    """
    Full f5c Profile HMM eventalign with detailed alignment output.

    RNA-only: Events are automatically reversed internally to match the
    5'->3' sequence direction. You do NOT need to reverse your sequence.

    This is the true f5c eventalign implementation using Viterbi HMM.
    Returns detailed event_alignment_t structures with HMM states.

    Args:
        raw_signal: 1D numpy array of raw nanopore signal (float32).
        sequence: Reference RNA sequence string in standard 5'->3' direction.
                  Do NOT reverse the sequence - event reversal is handled internally.
        kmer_size: Size of k-mers (5 or 9, default: 5).
        events_per_base: Expected events per base (default: 3.0).

    Returns:
        dict with keys:
            - alignment: List of dicts with full event_alignment_t data:
                * ref_position: Reference position (0-based)
                * ref_kmer: Reference k-mer string
                * event_idx: Event index (-1 for kmer skips)
                * hmm_state: 'M' (match), 'K' (kmer_skip), 'B' (bad_event)
                * event_mean, event_stdv, event_duration: Observed event stats
                * model_mean, model_stdv: Expected model stats
                * scaled_model_mean, scaled_model_stdv: Scaled model stats
            - scaling: dict with 'scale' and 'shift'
            - n_events: Total number of events detected
            - n_aligned: Number of alignment records
            - events_per_base: Events per base ratio used

    Raises:
        RuntimeError: If eventalign extension is not available.
    """
    if not _EVENTALIGN_AVAILABLE:
        raise RuntimeError(
            f"eventalign extension is not available. "
            f"Import error: {_eventalign_import_error}\n"
            f"Please ensure the package was built correctly with: pip install -e ."
        )

    # Validate inputs
    if not isinstance(raw_signal, np.ndarray):
        raise TypeError(f"raw_signal must be a numpy array (got {type(raw_signal)})")
    if raw_signal.ndim != 1:
        raise TypeError(f"raw_signal must be 1D (got {raw_signal.ndim}D)")
    if raw_signal.dtype != np.float32:
        raw_signal = raw_signal.astype(np.float32)

    if not isinstance(sequence, str):
        raise TypeError(f"sequence must be a string (got {type(sequence)})")

    if len(sequence) < kmer_size:
        raise ValueError(f"sequence length ({len(sequence)}) must be >= kmer_size ({kmer_size})")

    # Call C extension (RNA-only)
    return _profile_hmm_eventalign_c(raw_signal, sequence, None, kmer_size, events_per_base)


def is_available() -> bool:
    """
    Check if the f5c event detection extension is available.

    Returns:
        bool: True if extension is available, False otherwise
    """
    return _EVENT_AVAILABLE


def eventalign_is_available() -> bool:
    """
    Check if the eventalign extension (CPU) is available.

    Returns:
        bool: True if extension is available, False otherwise
    """
    return _EVENTALIGN_AVAILABLE


def eventalign_cuda_is_available() -> bool:
    """
    Check if the CUDA (GPU) eventalign extension is available.

    Returns:
        bool: True if CUDA extension is available, False otherwise
    """
    return _EVENTALIGN_CUDA_AVAILABLE


def profile_hmm_eventalign_cuda(
    raw_signal: np.ndarray,
    sequence: str,
    kmer_size: int = 5,
    events_per_base: float = 3.0,
) -> dict:
    """
    GPU-accelerated Profile HMM eventalign using CUDA.

    RNA-only: Events are automatically reversed internally to match the
    5'->3' sequence direction. You do NOT need to reverse your sequence.

    This is the CUDA-accelerated version of profile_hmm_eventalign.
    Provides significant speedup for large signals on NVIDIA GPUs.

    Args:
        raw_signal: 1D numpy array of raw nanopore signal (float32).
        sequence: Reference RNA sequence string in standard 5'->3' direction.
                  Do NOT reverse the sequence - event reversal is handled internally.
        kmer_size: Size of k-mers (5 or 9, default: 5).
        events_per_base: Expected events per base (default: 3.0).

    Returns:
        Same as profile_hmm_eventalign

    Raises:
        RuntimeError: If CUDA eventalign extension is not available.
    """
    if not _EVENTALIGN_CUDA_AVAILABLE:
        raise RuntimeError(
            f"CUDA eventalign extension is not available. "
            f"Import error: {_eventalign_cuda_import_error}\n"
            f"Build with CUDA: CUDA_HOME=/usr/local/cuda pip install -e ."
        )

    # Validate inputs
    if not isinstance(raw_signal, np.ndarray):
        raise TypeError(f"raw_signal must be a numpy array (got {type(raw_signal)})")
    if raw_signal.ndim != 1:
        raise TypeError(f"raw_signal must be 1D (got {raw_signal.ndim}D)")
    if raw_signal.dtype != np.float32:
        raw_signal = raw_signal.astype(np.float32)

    if not isinstance(sequence, str):
        raise TypeError(f"sequence must be a string (got {type(sequence)})")

    if len(sequence) < kmer_size:
        raise ValueError(f"sequence length ({len(sequence)}) must be >= kmer_size ({kmer_size})")

    # Call CUDA extension (RNA-only)
    return _profile_hmm_eventalign_cuda(raw_signal, sequence, None, kmer_size, events_per_base)


def eventalign_cuda(
    raw_signal: np.ndarray,
    sequence: str,
    kmer_size: int = 5,
    model: dict = None,
) -> dict:
    """
    GPU-accelerated event-to-sequence alignment using CUDA.

    RNA-only: Events are automatically reversed internally to match the
    5'->3' sequence direction. You do NOT need to reverse your sequence.

    This is the CUDA-accelerated version of eventalign.
    Provides significant speedup for large signals on NVIDIA GPUs.

    Args:
        raw_signal: 1D numpy array of raw nanopore signal (float32).
        sequence: Reference RNA sequence string in standard 5'->3' direction.
                  Do NOT reverse the sequence - event reversal is handled internally.
        kmer_size: Size of k-mers (5 or 9, default: 5).
        model: Optional k-mer model dictionary.

    Returns:
        Same as eventalign

    Raises:
        RuntimeError: If CUDA eventalign extension is not available.
    """
    if not _EVENTALIGN_CUDA_AVAILABLE:
        raise RuntimeError(
            f"CUDA eventalign extension is not available. "
            f"Import error: {_eventalign_cuda_import_error}\n"
            f"Build with CUDA: CUDA_HOME=/usr/local/cuda pip install -e ."
        )

    # Validate inputs
    if not isinstance(raw_signal, np.ndarray):
        raise TypeError(f"raw_signal must be a numpy array (got {type(raw_signal)})")
    if raw_signal.ndim != 1:
        raise TypeError(f"raw_signal must be 1D (got {raw_signal.ndim}D)")
    if raw_signal.dtype != np.float32:
        raw_signal = raw_signal.astype(np.float32)

    if not isinstance(sequence, str):
        raise TypeError(f"sequence must be a string (got {type(sequence)})")

    if len(sequence) < kmer_size:
        raise ValueError(f"sequence length ({len(sequence)}) must be >= kmer_size ({kmer_size})")

    # Call CUDA extension (RNA-only)
    return _eventalign_cuda(raw_signal, sequence, model, kmer_size)


# Export public functions
__all__ = [
    "detect_events",
    "is_available",
    "eventalign",
    "profile_hmm_eventalign",
    "eventalign_is_available",
]
