import numpy as np

# Try to import the compiled C extension
try:
    from ._event import get_events

    _EVENT_AVAILABLE = True
except ImportError as e:
    _EVENT_AVAILABLE = False
    _import_error = str(e)

# Try to import eventalign extension
try:
    from ._eventalign import eventalign as _eventalign_c
    from ._eventalign import profile_hmm_eventalign as _profile_hmm_eventalign_c

    _EVENTALIGN_AVAILABLE = True
except ImportError as e:
    _EVENTALIGN_AVAILABLE = False
    _eventalign_import_error = str(e)


def detect_events(raw_signal: np.ndarray, is_rna: bool = True) -> list[dict]:
    """
    High-level wrapper for nanopore event detection (C-backed).

    Args:
        raw_signal: 1D numpy array of raw nanopore signal values (must be float32).
        is_rna: If True, use RNA-specific detection parameters (default: DNA).

    Returns:
        List of event dicts, each with:
            - mean: Event mean current (float)
            - stdv: Event standard deviation (float)
            - start: Start index of the event (int)
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

    # Call C function (convert is_rna to int: 1/0)
    return get_events(raw_signal, int(is_rna))


def eventalign(
    raw_signal: np.ndarray,
    sequence: str,
    is_rna: bool = False,
    kmer_size: int = 5,
    model: dict = None,
) -> dict:
    """
    Align detected events to a reference sequence.

    This function performs event-to-sequence alignment, creating a mapping
    between k-mers in the sequence and events in the raw signal.

    Args:
        raw_signal: 1D numpy array of raw nanopore signal (float32).
        sequence: Reference DNA/RNA sequence string (e.g., "ACGTACGT").
        is_rna: If True, use RNA-specific detection parameters (default: False).
        kmer_size: Size of k-mers for alignment (default: 5).
        model: Optional k-mer model dictionary (default: uses built-in model).

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
        >>> sequence = "ACGTACGTACGT"
        >>> result = eventalign(signal, sequence, is_rna=False, kmer_size=5)
        >>> print(f"Detected {result['n_events']} events")
        >>> print(f"Scaling: {result['scaling']}")
        >>> for i, mapping in enumerate(result['base_to_event_map']):
        >>>     print(f"K-mer {mapping['kmer']}: events {mapping['start']}-{mapping['stop']}")
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

    # Call C extension
    return _eventalign_c(raw_signal, sequence, model, int(is_rna), kmer_size)


def profile_hmm_eventalign(
    raw_signal: np.ndarray,
    sequence: str,
    is_rna: bool = False,
    kmer_size: int = 5,
    events_per_base: float = 3.0,
) -> dict:
    """
    Full f5c Profile HMM eventalign with detailed alignment output.

    This is the true f5c eventalign implementation using Viterbi HMM.
    Returns detailed event_alignment_t structures with HMM states.

    Args:
        raw_signal: 1D numpy array of raw nanopore signal (float32).
        sequence: Reference DNA/RNA sequence string.
        is_rna: If True, use RNA-specific parameters (default: False).
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

    # Call C extension
    return _profile_hmm_eventalign_c(
        raw_signal, sequence, None, int(is_rna), kmer_size, events_per_base
    )


def is_available() -> bool:
    """
    Check if the f5c event detection extension is available.

    Returns:
        bool: True if extension is available, False otherwise
    """
    return _EVENT_AVAILABLE


def eventalign_is_available() -> bool:
    """
    Check if the eventalign extension is available.

    Returns:
        bool: True if extension is available, False otherwise
    """
    return _EVENTALIGN_AVAILABLE


# Export public functions
__all__ = [
    "detect_events",
    "is_available",
    "eventalign",
    "profile_hmm_eventalign",
    "eventalign_is_available",
]
