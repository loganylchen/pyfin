"""
F5C Event Alignment Module

This module provides Python bindings for the f5c event alignment algorithm,
which aligns nanopore signal events to reference k-mers.

The module supports:
- Single read to single reference alignment
- Batch alignment of multiple reads to multiple references
- RNA002 (k=5) and RNA004 (k=9) models
- Both CPU and GPU (CUDA) acceleration (if available)

Example usage:
    >>> from fin._eventalign import EventAligner
    >>> import numpy as np
    >>>
    >>> # Initialize with RNA002 model (k=5)
    >>> aligner = EventAligner(model=1)
    >>>
    >>> # Align a single read
    >>> signal = np.random.randn(10000).astype(np.float32)
    >>> result = aligner.align_read_single_ref(
    ...     signal=signal,
    ...     read_name="read001",
    ...     ref_sequence="ACGUACGUACGU",
    ...     ref_name="chr1",
    ...     sample_rate=4000.0
    ... )
    >>> print(f"Aligned {result['n_alignments']} pairs")
"""

import sys

# Try to import CPU version
try:
    from ._eventalign_cpu import EventAligner as _EventAlignerCPU
    from ._eventalign_cpu import MODEL_RNA002, MODEL_RNA004
    _EVENTALIGN_CPU_AVAILABLE = True
except ImportError as e:
    _EVENTALIGN_CPU_AVAILABLE = False
    _eventalign_cpu_import_error = str(e)
    _EventAlignerCPU = None
    MODEL_RNA002 = 1
    MODEL_RNA004 = 2

# Try to import GPU version
try:
    from ._eventalign_gpu import EventAligner as _EventAlignerGPU
    _EVENTALIGN_GPU_AVAILABLE = True
except ImportError as e:
    _EVENTALIGN_GPU_AVAILABLE = False
    _eventalign_gpu_import_error = str(e)
    _EventAlignerGPU = None

__all__ = ["EventAligner", "MODEL_RNA002", "MODEL_RNA004"]


def _check_available():
    """Check if any eventalign extension is available."""
    if not _EVENTALIGN_CPU_AVAILABLE and not _EVENTALIGN_GPU_AVAILABLE:
        raise ImportError(
            f"No EventAlign C extension available.\n"
            f"CPU error: {_eventalign_cpu_import_error}\n"
            f"GPU error: {_eventalign_gpu_import_error}\n"
            "Please rebuild the package with: pip install -e ."
        )


# Create a unified EventAligner class that uses GPU if available, falls back to CPU
class EventAligner:
    """
    Unified EventAligner class that automatically uses GPU acceleration if available.

    This class wraps both CPU and GPU implementations, automatically selecting
    the best available backend.

    Args:
        model: Model type (1=RNA002 with k=5, 2=RNA004 with k=9)
        use_gpu: Force GPU (True) or CPU (False). Auto-detect if None.

    Attributes:
        model_id: The model type being used
        kmer_size: The k-mer size (5 for RNA002, 9 for RNA004)
        using_gpu: Whether GPU acceleration is being used

    Example:
        >>> aligner = EventAligner(model=1)  # Auto-detect GPU
        >>> result = aligner.align_read_single_ref(signal, "read1", "ACGUACGU", "chr1")
    """

    def __init__(self, model: int = 1, use_gpu: bool = None):
        _check_available()

        if model != 1 and model != 2:
            raise ValueError("model must be 1 (RNA002) or 2 (RNA004)")

        self.model_id = model
        self.kmer_size = 5 if model == 1 else 9

        # Determine which backend to use
        if use_gpu is None:
            # Auto-detect: prefer GPU if available
            self._using_gpu = _EVENTALIGN_GPU_AVAILABLE
        else:
            self._using_gpu = use_gpu

        if self._using_gpu and not _EVENTALIGN_GPU_AVAILABLE:
            raise ImportError(
                f"GPU backend requested but not available. Error: {_eventalign_gpu_import_error}"
            )

        # Create the appropriate backend
        if self._using_gpu:
            self._backend = _EventAlignerGPU(model)
        else:
            if not _EVENTALIGN_CPU_AVAILABLE:
                raise ImportError(
                    f"CPU backend not available. Error: {_eventalign_cpu_import_error}"
                )
            self._backend = _EventAlignerCPU(model)

    @property
    def using_gpu(self) -> bool:
        """Whether GPU acceleration is being used."""
        return self._using_gpu

    def align_read_single_ref(self, signal, read_name: str, ref_sequence: str,
                             ref_name: str, sample_rate: float = 4000.0) -> dict:
        """
        Align a single read to a single reference sequence.

        Args:
            signal: numpy array of float32 raw signal
            read_name: name identifier for the read
            ref_sequence: reference sequence string
            ref_name: name identifier for the reference
            sample_rate: signal sample rate in Hz (default 4000)

        Returns:
            dict with alignment results including:
                - read_name: name of the read
                - ref_name: name of the reference
                - success: whether alignment succeeded
                - n_events: number of events detected
                - n_alignments: number of aligned pairs
                - events_per_base: events per base ratio
                - ref_positions: numpy array of reference positions
                - read_positions: numpy array of event positions
        """
        return self._backend.align_read_single_ref(
            signal, read_name, ref_sequence, ref_name, sample_rate
        )

    def align_batch(self, signals, read_names, ref_sequences, ref_names,
                    ref_lengths, sample_rate: float = 4000.0) -> dict:
        """
        Align multiple reads to multiple reference sequences.

        Args:
            signals: list of numpy float32 arrays
            read_names: list of read name strings
            ref_sequences: list of reference sequence strings
            ref_names: list of reference name strings
            ref_lengths: list of reference sequence lengths
            sample_rate: signal sample rate in Hz (default 4000)

        Returns:
            dict mapping (read_name, ref_name) tuples to alignment result dicts
        """
        return self._backend.align_batch(
            signals, read_names, ref_sequences, ref_names, ref_lengths, sample_rate
        )


# Module-level convenience functions
def create_aligner(model: int = 1, use_gpu: bool = None):
    """
    Create an EventAligner instance.

    Args:
        model: Model type (1=RNA002 with k=5, 2=RNA004 with k=9)
        use_gpu: Force GPU (True) or CPU (False). Auto-detect if None.

    Returns:
        EventAligner instance

    Raises:
        ImportError: If no C extension is available
        ValueError: If model is not 1 or 2
    """
    return EventAligner(model, use_gpu)


def align_single(signal, read_name: str, ref_sequence: str, ref_name: str,
                 model: int = 1, sample_rate: float = 4000.0, use_gpu: bool = None) -> dict:
    """
    Convenience function to align a single read to a single reference.

    Args:
        signal: numpy array of float32 raw signal
        read_name: name identifier for the read
        ref_sequence: reference sequence string
        ref_name: name identifier for the reference
        model: model type (1=RNA002, 2=RNA004)
        sample_rate: signal sample rate in Hz
        use_gpu: Force GPU (True) or CPU (False). Auto-detect if None.

    Returns:
        dict with alignment results

    Example:
        >>> import numpy as np
        >>> from fin._eventalign import align_single
        >>> signal = np.random.randn(10000).astype(np.float32)
        >>> result = align_single(signal, "read1", "ACGUACGU", "chr1")
        >>> print(f"Success: {result['success']}, Alignments: {result['n_alignments']}")
    """
    aligner = create_aligner(model, use_gpu)
    return aligner.align_read_single_ref(signal, read_name, ref_sequence,
                                        ref_name, sample_rate)


def align_batch(signals, read_names, ref_sequences, ref_names, ref_lengths,
                model: int = 1, sample_rate: float = 4000.0, use_gpu: bool = None) -> dict:
    """
    Convenience function to align multiple reads to multiple references.

    Args:
        signals: list of numpy float32 arrays
        read_names: list of read name strings
        ref_sequences: list of reference sequence strings
        ref_names: list of reference name strings
        ref_lengths: list of reference sequence lengths
        model: model type (1=RNA002, 2=RNA004)
        sample_rate: signal sample rate in Hz
        use_gpu: Force GPU (True) or CPU (False). Auto-detect if None.

    Returns:
        dict mapping (read_name, ref_name) tuples to alignment result dicts

    Example:
        >>> import numpy as np
        >>> from fin._eventalign import align_batch
        >>> signals = [np.random.randn(10000).astype(np.float32) for _ in range(2)]
        >>> results = align_batch(
        ...     signals,
        ...     ["read1", "read2"],
        ...     ["ACGUACGU", "UGCAUGCA"],
        ...     ["chr1", "chr2"],
        ...     [8, 8]
        ... )
        >>> for key, result in results.items():
        ...     print(f"{key}: {result['n_alignments']} alignments")
    """
    aligner = create_aligner(model, use_gpu)
    return aligner.align_batch(signals, read_names, ref_sequences,
                             ref_names, ref_lengths, sample_rate)


# Print availability status on import
def _print_status():
    """Print the status of available backends."""
    if _EVENTALIGN_GPU_AVAILABLE:
        print("🚀 F5C EventAlign: GPU (CUDA) acceleration ENABLED")
    elif _EVENTALIGN_CPU_AVAILABLE:
        print("💻 F5C EventAlign: Using CPU implementation")
    else:
        print("⚠️  F5C EventAlign: No backend available - please rebuild the package")


_print_status()

# Add convenience functions to __all__
__all__ += ["create_aligner", "align_single", "align_batch"]
