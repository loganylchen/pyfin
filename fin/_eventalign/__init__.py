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
import os
import traceback

# Try to import the C extension
_eventalign_import_error = None
_eventalign_diagnostics = None
_EVENTALIGN_AVAILABLE = False

try:
    from ._eventalign import getevents as _getevents_c
    from ._eventalign import set_model as _set_model_c
    from ._eventalign import init_db_from_python as _init_db_from_python_c
    from ._eventalign import free_db as _free_db_c
    from ._eventalign import MODEL_RNA002, MODEL_RNA004, MAX_KMER_SIZE, MAX_NUM_KMER
    _EVENTALIGN_AVAILABLE = True
except ImportError as e:
    _EVENTALIGN_AVAILABLE = False
    _eventalign_import_error = str(e)

    # Collect diagnostic information
    import sysconfig
    diagnostics = []

    # Python version info
    diagnostics.append(f"Python version: {sys.version}")
    diagnostics.append(f"Platform: {sys.platform}")

    # Check if extension file exists
    ext_path = os.path.join(os.path.dirname(__file__), "_eventalign.so")
    if os.path.exists(ext_path):
        diagnostics.append(f"Extension file exists: {ext_path}")
        diagnostics.append(f"Extension file size: {os.path.getsize(ext_path)} bytes")
    else:
        diagnostics.append(f"Extension file NOT found: {ext_path}")

    # Check for common missing dependencies
    try:
        import numpy
        diagnostics.append(f"NumPy available: {numpy.__version__}")
    except ImportError:
        diagnostics.append("NumPy: NOT INSTALLED")

    # Compiler info
    diagnostics.append(f"ABI flags: {sysconfig.get_config_var('ABI')}")
    diagnostics.append(f"SOABI: {sysconfig.get_config_var('SOABI')}")

    # Full traceback
    diagnostics.append("\n--- Full traceback ---")
    diagnostics.append(traceback.format_exc())

    _eventalign_diagnostics = "\n  ".join(diagnostics)

    # Define constants for when module is not available
    MODEL_RNA002 = 1
    MODEL_RNA004 = 2
    MAX_KMER_SIZE = 9
    MAX_NUM_KMER = 262144

__all__ = [
    "getevents",
    "set_model",
    "init_db_from_python",
    "free_db",
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
            f"\nDiagnostics:\n  {_eventalign_diagnostics}\n"
            "\nTo fix this issue, try rebuilding the package:\n"
            "  pip install -e .\n"
            "Or clean build and reinstall:\n"
            "  pip uninstall -y py-fin && rm -rf build/fin/_eventalign && pip install -e ."
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


def init_db_from_python(
    read_ids,
    read_seqs,
    ref_seqs,
    ref_names,
    ref_lens,
    signals,
    signal_drifts=None,
    signal_scales=None,
    signal_shifts=None,
):
    """
    Initialize a db_t structure from Python-provided data.

    This function creates a db_t structure (used internally by f5c for event alignment)
    populated with data provided directly from Python, rather than reading from BAM/SLOW5 files.

    Args:
        read_ids: list of read identifier strings
        read_seqs: list of read sequence strings
        ref_seqs: list of reference sequence strings
        ref_names: list of reference name strings
        ref_lens: list of reference sequence lengths (int)
        signals: list of 1D float32 numpy arrays (raw signal data)
        signal_drifts: list of float drift values (optional, default=None)
        signal_scales: list of float scale values (optional, default=None)
        signal_shifts: list of float shift values (optional, default=None)

    Returns:
        int: pointer to db_t structure (as Python int). Use this pointer with other
             functions that operate on db_t, and call free_db() when done.

    Example:
        >>> from fin._eventalign import init_db_from_python, free_db
        >>> import numpy as np
        >>>
        >>> read_ids = ['read1', 'read2']
        >>> read_seqs = ['ACGTACGT', 'TGCATGCA']
        >>> ref_seqs = ['ACGTACGT', 'TGCATGCA']
        >>> ref_names = ['chr1', 'chr1']
        >>> ref_lens = [8, 8]
        >>> signals = [
        ...     np.random.randn(10000).astype(np.float32),
        ...     np.random.randn(12000).astype(np.float32),
        ... ]
        >>>
        >>> db_ptr = init_db_from_python(
        ...     read_ids, read_seqs, ref_seqs, ref_names, ref_lens, signals
        ... )
        >>> # ... use db_ptr with other functions ...
        >>> free_db(db_ptr)
    """
    _check_available()
    return _init_db_from_python_c(
        read_ids,
        read_seqs,
        ref_seqs,
        ref_names,
        ref_lens,
        signals,
        signal_drifts,
        signal_scales,
        signal_shifts,
    )


def free_db(db_ptr):
    """
    Free a db_t structure created by init_db_from_python.

    Args:
        db_ptr: pointer to db_t structure (as returned by init_db_from_python)

    Example:
        >>> from fin._eventalign import init_db_from_python, free_db
        >>> db_ptr = init_db_from_python(...)
        >>> free_db(db_ptr)
    """
    _check_available()
    _free_db_c(db_ptr)


# Print availability status on import
def _print_status():
    """Print the status of the eventalign extension."""
    if _EVENTALIGN_AVAILABLE:
        print("[fin._eventalign] Event detection and model API loaded successfully")
    else:
        import logging
        logger = logging.getLogger(__name__)
        logger.warning(
            f"f5c extension not available - event detection disabled\n"
            f"Error: {_eventalign_import_error}\n"
            f"Diagnostics:\n  {_eventalign_diagnostics}"
        )
        print(
            f"[fin._eventalign] WARNING: C extension not available - event detection disabled\n"
            f"  Error: {_eventalign_import_error}\n"
            f"  Diagnostics:\n    {_eventalign_diagnostics}\n"
            f"  To fix: pip install -e ."
        )


_print_status()
