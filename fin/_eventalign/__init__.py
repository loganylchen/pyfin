"""
F5C Event Detection and Model API Module

This module provides Python bindings for the f5c event detection and model loading functions.

The module supports:
- Event detection from raw nanopore signal
- Loading pore models (RNA002 and RNA004)

Example usage:
    >>> from fin._eventalign import getevents, set_model, MODEL_RNA002
    >>> import numpy as np
    >>>
    >>> # Detect events from signal
    >>> signal = np.random.randn(10000).astype(np.float32)
    >>> events = getevents(signal)
    >>> print(f"Detected {events['n_events']} events")
    >>>
    >>> # Load pore model
    >>> model = set_model(MODEL_RNA002)
    >>> print(f"K-mer size: {model['kmer_size']}")
"""

import sys

# Try to import the C extension
try:
    from ._eventalign import getevents as _getevents_c
    from ._eventalign import set_model as _set_model_c
    from ._eventalign import MODEL_RNA002, MODEL_RNA004, MAX_KMER_SIZE, MAX_NUM_KMER
    _EVENTALIGN_AVAILABLE = True
except ImportError as e:
    _EVENTALIGN_AVAILABLE = False
    _eventalign_import_error = str(e)
    # Define constants for when module is not available
    MODEL_RNA002 = 1
    MODEL_RNA004 = 2
    MAX_KMER_SIZE = 9
    MAX_NUM_KMER = 262144

__all__ = [
    "getevents",
    "set_model",
    "MODEL_RNA002",
    "MODEL_RNA004",
    "MAX_KMER_SIZE",
    "MAX_NUM_KMER",
]


def _check_available():
    """Check if the eventalign extension is available."""
    if not _EVENTALIGN_AVAILABLE:
        raise ImportError(
            f"Eventalign C extension not available.\n"
            f"Error: {_eventalign_import_error}\n"
            "Please rebuild the package with: pip install -e ."
        )


def getevents(signal):
    """
    Detect events from raw nanopore signal.

    This function performs event detection on raw nanopore current signal,
    segmenting the continuous signal into discrete events characterized by
    their start position, length, mean, and standard deviation.

    Args:
        signal: numpy array of float32 raw signal values

    Returns:
        dict with keys:
            - n_events: number of events detected
            - starts: numpy uint64 array of event start positions (in samples)
            - lengths: numpy float32 array of event lengths (in samples)
            - means: numpy float32 array of event mean current values
            - stdvs: numpy float32 array of event standard deviations

    Example:
        >>> import numpy as np
        >>> from fin._eventalign import getevents
        >>>
        >>> # Generate synthetic signal
        >>> signal = np.random.randn(10000).astype(np.float32)
        >>>
        >>> # Detect events
        >>> events = getevents(signal)
        >>>
        >>> print(f"Detected {events['n_events']} events")
        >>> print(f"First event starts at sample {events['starts'][0]}")
        >>> print(f"First event mean: {events['means'][0]:.2f}")
    """
    _check_available()
    return _getevents_c(signal)


def set_model(model_id: int = MODEL_RNA002):
    """
    Load a pore model for nanopore signal alignment.

    This function loads the pore model parameters (mean levels and standard deviations)
    for all k-mers of a given model type.

    Args:
        model_id: Model type - 1 for RNA002 (k=5), 2 for RNA004 (k=9)

    Returns:
        dict with keys:
            - kmer_size: k-mer size (5 for RNA002, 9 for RNA004)
            - num_kmer: number of k-mers (4^kmer_size)
            - level_means: numpy float32 array of k-mer mean levels
            - level_stdvs: numpy float32 array of k-mer standard deviations

    Raises:
        ValueError: If model_id is not 1 or 2

    Example:
        >>> from fin._eventalign import set_model, MODEL_RNA002, MODEL_RNA004
        >>>
        >>> # Load RNA002 model (k=5, 1024 kmers)
        >>> model_002 = set_model(MODEL_RNA002)
        >>> print(f"RNA002: k={model_002['kmer_size']}, {model_002['num_kmer']} kmers")
        >>>
        >>> # Load RNA004 model (k=9, 262144 kmers)
        >>> model_004 = set_model(MODEL_RNA004)
        >>> print(f"RNA004: k={model_004['kmer_size']}, {model_004['num_kmer']} kmers")
    """
    _check_available()
    return _set_model_c(model_id)


# Print availability status on import
def _print_status():
    """Print the status of the eventalign extension."""
    if _EVENTALIGN_AVAILABLE:
        print("[fin._eventalign] Event detection and model API loaded successfully")
    else:
        print(f"[fin._eventalign] Warning: C extension not available - {_eventalign_import_error}")


_print_status()
