import numpy as np
from . import get_events  # Import the compiled C function

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
        RuntimeError: If event detection fails (no events found).
    """
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

# Export only the high-level function (hide low-level _f5c)
__all__ = ["detect_events"]