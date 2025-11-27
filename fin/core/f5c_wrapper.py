"""
Python wrapper for f5c C extension module.

This module provides a Pythonic interface to the f5c event alignment functionality.
The actual implementation calls C functions from the f5c core library.
"""

import numpy as np
from typing import List, Tuple, Optional
import warnings

# Try to import the compiled f5c module
try:
    from fin._f5c import f5c_python
    F5C_AVAILABLE = True
except ImportError:
    warnings.warn(
        "f5c C extension not available. Install with: pip install -e .\n"
        "Falling back to pure Python implementation.",
        ImportWarning
    )
    F5C_AVAILABLE = False


@dataclass
class Event:
    """Nanopore signal event."""
    mean: float
    stdv: float
    start: int
    length: int


@dataclass
class AlignedEvent:
    """Event aligned to reference."""
    ref_position: int
    ref_kmer: str
    event_idx: int
    event_mean: float
    event_stdv: float
    event_length: int
    strand: str
    alignment_score: float


class F5CWrapper:
    """
    Python wrapper for f5c event detection and alignment.

    This class provides a clean interface to the f5c C library functions,
    handling the conversion between Python objects and C data structures.
    """

    def __init__(self, use_c_extension: bool = True):
        """
        Initialize f5c wrapper.

        Args:
            use_c_extension: Whether to use the C extension (if available)
        """
        self.use_c_extension = use_c_extension and F5C_AVAILABLE

        if not self.use_c_extension:
            warnings.warn(
                "Using pure Python fallback implementation. "
                "For better performance, install the C extension.",
                RuntimeWarning
            )
            # Initialize fallback implementation
            self._init_fallback()

    def _init_fallback(self):
        """Initialize pure Python fallback."""
        from .event_detection import EventDetector
        from .eventalign import AdaptiveBandedAligner

        self.event_detector = EventDetector()
        self.aligner = AdaptiveBandedAligner()

    def detect_events(
        self,
        signal: np.ndarray,
        outlier_threshold: float = 3.0,
        min_event_length: int = 3,
        sampling_rate: float = 4000.0
    ) -> List[Event]:
        """
        Detect events from raw nanopore signal using f5c algorithm.

        Args:
            signal: Raw ionic current signal (1D array)
            outlier_threshold: Threshold for outlier removal
            min_event_length: Minimum samples per event
            sampling_rate: Sampling rate in Hz (default: 4000 for MinION)

        Returns:
            List of detected events
        """
        if self.use_c_extension:
            return self._detect_events_c(
                signal, outlier_threshold, min_event_length
            )
        else:
            return self._detect_events_python(
                signal, outlier_threshold, min_event_length, sampling_rate
            )

    def _detect_events_c(
        self,
        signal: np.ndarray,
        outlier_threshold: float,
        min_event_length: int
    ) -> List[Event]:
        """
        Call f5c C function to detect events.

        Note: This is a placeholder that calls the C extension.
        In the actual implementation, this would call f5c's detect_events().
        """
        # Ensure contiguous float32 array
        signal_f32 = np.asarray(signal, dtype=np.float32, order='C')

        # Call C function (simplified implementation in current C file)
        result = f5c_python.detect_events(
            signal_f32, outlier_threshold, min_event_length
        )

        # Parse result
        means, stdvs, starts, lengths = result

        events = []
        for i in range(len(means)):
            events.append(Event(
                mean=means[i],
                stdv=stdvs[i],
                start=starts[i],
                length=lengths[i]
            ))

        return events

    def _detect_events_python(
        self,
        signal: np.ndarray,
        outlier_threshold: float,
        min_event_length: int,
        sampling_rate: float
    ) -> List[Event]:
        """
        Pure Python fallback for event detection.
        """
        events_dummy = self.event_detector.detect_events(signal)

        # Convert to Event dataclass
        return [
            Event(mean=e.mean, stdv=e.stdv, start=e.start, length=e.length)
            for e in events_dummy
        ]

    def align_to_sequence(
        self,
        events: List[Event],
        sequence: str,
        is_rna: bool = False,
        band_width: int = 100
    ) -> List[AlignedEvent]:
        """
        Align events to reference sequence using f5c banded DP algorithm.

        Args:
            events: List of detected events
            sequence: Reference sequence
            is_rna: Whether this is RNA (triggers reverse complement)
            band_width: Band width for alignment

        Returns:
            List of aligned events with reference positions
        """
        if self.use_c_extension:
            return self._align_to_sequence_c(
                events, sequence, is_rna, band_width
            )
        else:
            return self._align_to_sequence_python(
                events, sequence, is_rna, band_width
            )

    def _align_to_sequence_c(
        self,
        events: List[Event],
        sequence: str,
        is_rna: bool,
        band_width: int
    ) -> List[AlignedEvent]:
        """Call f5c C function for alignment."""
        # Convert events to numpy arrays
        n_events = len(events)

        means = np.array([e.mean for e in events], dtype=np.float32)
        stdvs = np.array([e.stdv for e in events], dtype=np.float32)
        starts = np.array([e.start for e in events], dtype=np.int32)
        lengths = np.array([e.length for e in events], dtype=np.int32)

        # Call C function
        result = f5c_python.align_events_to_sequence(
            means, stdvs, starts, lengths,
            sequence, is_rna, band_width
        )

        # Parse result
        ref_positions, ref_kmers, event_indices, scores = result

        aligned_events = []
        for i in range(len(ref_positions)):
            event_idx = event_indices[i]
            event = events[event_idx]

            aligned_events.append(AlignedEvent(
                ref_position=ref_positions[i],
                ref_kmer=ref_kmers[i],
                event_idx=event_idx,
                event_mean=event.mean,
                event_stdv=event.stdv,
                event_length=event.length,
                strand="-" if is_rna else "+",
                alignment_score=scores[i]
            ))

        return aligned_events

    def _align_to_sequence_python(
        self,
        events: List[Event],
        sequence: str,
        is_rna: bool,
        band_width: int
    ) -> List[AlignedEvent]:
        """Pure Python fallback for alignment."""
        aligned_dummy = self.aligner.align(
            events, sequence, is_rna=is_rna
        )

        # Convert to AlignedEvent dataclass
        return [
            AlignedEvent(
                ref_position=ae.ref_position,
                ref_kmer=ae.ref_kmer,
                event_idx=ae.event_idx,
                event_mean=ae.event_mean,
                event_stdv=ae.event_stdv,
                event_length=ae.event_length,
                strand=ae.strand,
                alignment_score=ae.alignment_score
            )
            for ae in aligned_dummy
        ]

    def call_eventalign(
        self,
        raw_signal: np.ndarray,
        reference_sequence: str,
        is_rna: bool = False
    ) -> tuple[List[Event], List[AlignedEvent]]:
        """
        Complete eventalign workflow: detect events and align to sequence.

        Args:
            raw_signal: Raw nanopore signal
            reference_sequence: Reference sequence
            is_rna: Whether this is RNA

        Returns:
            Tuple of (events, aligned_events)
        """
        # Detect events
        events = self.detect_events(raw_signal)

        if not events:
            return [], []

        # Align to sequence
        aligned_events = self.align_to_sequence(
            events, reference_sequence, is_rna=is_rna
        )

        return events, aligned_events


# Convenience functions

def f5c_detect_events(
    signal: np.ndarray,
    **kwargs
) -> List[Event]:
    """
    Convenience function: detect events using f5c algorithm.

    Args:
        signal: Raw ionic current signal
        **kwargs: Additional arguments for F5CWrapper.detect_events

    Returns:
        List of detected events
    """
    wrapper = F5CWrapper()
    return wrapper.detect_events(signal, **kwargs)


def f5c_eventalign(
    signal: np.ndarray,
    sequence: str,
    is_rna: bool = False
) -> tuple[List[Event], List[AlignedEvent]]:
    """
    Convenience function: perform eventalign using f5c algorithm.

    Args:
        signal: Raw nanopore signal
        sequence: Reference sequence
        is_rna: Whether this is RNA

    Returns:
        Tuple of (events, aligned_events)
    """
    wrapper = F5CWrapper()
    return wrapper.call_eventalign(signal, sequence, is_rna=is_rna)
