"""
Signal processing module for nanopore DRS raw signals.

This module handles:
- Loading and preprocessing of nanopore signal data
- Signal normalization and quality control
- Signal alignment and segmentation
- Event detection and alignment to sequences
"""

import numpy as np
from typing import List, Optional, Dict, Any
from .event_detection import EventDetector
from .eventalign import AdaptiveBandedAligner, Event


class SignalProcessor:
    """
    Process nanopore raw signals from fast5/blow5/pod5 files.

    Integrates event detection and alignment capabilities.
    """

    def __init__(
        self,
        band_width: int = 100,
        kmer_size: int = 5,
        window_size: int = 5,
        min_event_length: int = 3
    ):
        """
        Initialize signal processor.

        Args:
            band_width: Band width for event alignment
            kmer_size: K-mer size for pore model (5 for DNA/R9.4, 6 for RNA004)
            window_size: Smoothing window size for event detection
            min_event_length: Minimum samples per event
        """
        self.event_detector = EventDetector(
            window_size=window_size,
            min_event_length=min_event_length
        )
        self.aligner = AdaptiveBandedAligner(
            band_width=band_width,
            kmer_size=kmer_size
        )
        self.normalization_params = None

    def load_signal(self, signal_file: str) -> Dict[str, Any]:
        """
        Load nanopore signal from file.

        Args:
            signal_file: Path to signal file (fast5/blow5/pod5)

        Returns:
            dict: Raw signal data with metadata
        """
        raise NotImplementedError("Signal loading not yet implemented")

    def normalize_signal(self, signal: np.ndarray) -> np.ndarray:
        """
        Normalize nanopore signal to correct for technical variation.

        Args:
            signal: Raw signal array

        Returns:
            array: Normalized signal
        """
        # Scale signal by median and MAD
        median = np.median(signal)
        mad = np.median(np.abs(signal - median))
        # MAD to std: std = MAD * 1.4826
        scale = mad * 1.4826 if mad > 0 else 1.0

        normalized = (signal - median) / scale
        return normalized

    def denoise_signal(self, signal: np.ndarray, method: str = 'wavelet') -> np.ndarray:
        """
        Denoise signal using various methods.

        Args:
            signal: Signal to denoise
            method: 'wavelet' or 'median'

        Returns:
            Denoised signal
        """
        if method == 'median':
            # Apply median filter
            from scipy.signal import medfilt
            return medfilt(signal, kernel_size=self.event_detector.window_size)
        elif method == 'wavelet':
            # Simple thresholding approximation
            # In practice, use PyWavelets
            return self._wavelet_denoise(signal)
        else:
            return signal

    def _wavelet_denoise(self, signal: np.ndarray, threshold: float = 0.1) -> np.ndarray:
        """
        Simple wavelet denoising approximation.

        Args:
            signal: Signal to denoise
            threshold: Threshold for denoising

        Returns:
            Denoised signal
        """
        # Use FFT-based approximation for denoising
        # This is a placeholder - use PyWavelets for real wavelet denoising
        freq = np.fft.fft(signal)
        freq_abs = np.abs(freq)

        # Simple thresholding
        mask = freq_abs > (threshold * np.max(freq_abs))
        freq_filtered = freq * mask

        denoised = np.real(np.fft.ifft(freq_filtered))
        return denoised

    def detect_events(self, signal: np.ndarray) -> List[Event]:
        """
        Detect events from raw signal.

        Args:
            signal: Raw ionic current signal

        Returns:
            List of detected events
        """
        return self.event_detector.detect_events(signal)

    def align_to_sequence(
        self,
        events: List[Event],
        sequence: str,
        is_rna: bool = False
    ) -> List:
        """
        Align events to reference sequence.

        Args:
            events: Detected events
            sequence: Reference sequence
            is_rna: Whether this is RNA (requires reverse complement)

        Returns:
            List of aligned events
        """
        return self.aligner.align(events, sequence, is_rna=is_rna)

    def call_events_and_align(
        self,
        raw_signal: np.ndarray,
        reference_sequence: str,
        is_rna: bool = False
    ) -> tuple[List[Event], List]:
        """
        Convenience method: detect events and align to sequence.

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
        aligned_events = self.align_to_sequence(events, reference_sequence, is_rna=is_rna)

        return events, aligned_events

    def summarize_alignment(self, aligned_events: List) -> Dict[str, Any]:
        """
        Generate summary statistics from aligned events.

        Args:
            aligned_events: List of aligned events

        Returns:
            Dictionary of summary statistics
        """
        if not aligned_events:
            return {
                'n_aligned': 0,
                'mean_score': 0.0,
                'coverage': 0.0,
                'mean_dwell_time': 0.0
            }

        n_aligned = len(aligned_events)
        scores = [e.alignment_score for e in aligned_events]
        dwell_times = [e.event_length for e in aligned_events]

        summary = {
            'n_aligned': n_aligned,
            'mean_score': np.mean(scores),
            'std_score': np.std(scores),
            'min_score': np.min(scores),
            'max_score': np.max(scores),
            'mean_dwell_time': np.mean(dwell_times),
            'std_dwell_time': np.std(dwell_times),
            'mean_dwell_ms': np.mean(dwell_times) / 4.0,  # Assuming 4000 Hz sampling
        }

        # Calculate coverage of reference
        if aligned_events:
            positions = [e.ref_position for e in aligned_events]
            summary['coverage'] = len(set(positions))

        return summary

