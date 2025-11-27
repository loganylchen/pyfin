"""
Event detection module for segmenting nanopore raw signals.

Extracts discrete events from continuous ionic current signals using
segmentation algorithms inspired by Scrappie and f5c.
"""

import numpy as np
from typing import List, Optional, Tuple
from scipy import signal as scipy_signal

from .eventalign import Event


class EventDetector:
    """
    Detect events from nanopore raw ionic current signals.

    Implements segmentation algorithms to identify discrete events representing
    nucleotides or groups of nucleotides passing through the pore.
    """

    def __init__(
        self,
        window_size: int = 5,
        min_event_length: int = 3,
        outlier_threshold: float = 3.0
    ):
        """
        Initialize event detector.

        Args:
            window_size: Size of smoothing window
            min_event_length: Minimum samples per event
            outlier_threshold: Std dev threshold for outlier removal
        """
        self.window_size = window_size
        self.min_event_length = min_event_length
        self.outlier_threshold = outlier_threshold

    def detect_events(self, raw_signal: np.ndarray) -> List[Event]:
        """
        Detect events from raw nanopore signal.

        Args:
            raw_signal: Raw ionic current signal (1D array)

        Returns:
            List of detected events
        """
        # Step 1: Preprocess signal
        smoothed = self._preprocess_signal(raw_signal)

        # Step 2: Segment signal into discrete events
        event_boundaries = self._segment_events(smoothed)

        # Step 3: Extract event statistics
        events = self._extract_events(smoothed, event_boundaries)

        # Step 4: Filter events
        filtered_events = self._filter_events(events)

        return filtered_events

    def _preprocess_signal(self, signal: np.ndarray) -> np.ndarray:
        """
        Preprocess raw signal: smooth and remove outliers.

        Args:
            signal: Raw signal

        Returns:
            Preprocessed signal
        """
        # Remove extreme outliers (likely artifacts)
        median = np.median(signal)
        mad = np.median(np.abs(signal - median))
        std_est = mad * 1.4826  # Convert MAD to std dev

        # Clip extreme values
        if std_est > 0:
            lower = median - self.outlier_threshold * std_est
            upper = median + self.outlier_threshold * std_est
            signal_clipped = np.clip(signal, lower, upper)
        else:
            signal_clipped = signal

        # Apply low-pass filter (moving average)
        if len(signal_clipped) > self.window_size:
            window = np.ones(self.window_size) / self.window_size
            smoothed = scipy_signal.convolve(signal_clipped, window, mode='valid')
            # Pad to match original length
            pad_len = len(signal_clipped) - len(smoothed)
            smoothed = np.pad(smoothed, (pad_len // 2, pad_len - pad_len // 2), mode='edge')
        else:
            smoothed = signal_clipped

        return smoothed

    def _segment_events(self, signal: np.ndarray) -> List[Tuple[int, int]]:
        """
        Segment signal into discrete events using change point detection.

        Uses a simplified version of the algorithm from Scrappie:
        - Detect sudden changes in signal level
        - Merge short segments
        - Ensure minimum event length

        Args:
            signal: Preprocessed signal

        Returns:
            List of (start, end) boundaries for each event
        """
        if len(signal) < self.min_event_length:
            return [(0, len(signal))]

        # Calculate signal differences (first derivative)
        diff = np.diff(signal)

        # Find change points where derivative exceeds threshold
        # Threshold based on signal gradient magnitude
        threshold = np.std(diff) * 0.5

        change_points = [0]  # Start at beginning

        for i in range(len(diff)):
            if abs(diff[i]) > threshold:
                # Potential change point
                # Verify it persists (not transient noise)
                if i + 3 < len(signal):
                    # Check if signal stabilizes at new level
                    next_segment = signal[i+1:i+4]
                    prev_segment = signal[max(0, i-3):i]

                    if (np.std(next_segment) < np.std(prev_segment) * 0.7 and
                        abs(np.mean(next_segment) - np.mean(prev_segment)) > threshold):
                        change_points.append(i + 1)

        change_points.append(len(signal))  # End at signal length

        # Create segments
        segments = []
        for i in range(len(change_points) - 1):
            start = change_points[i]
            end = change_points[i + 1]

            # Ensure minimum length
            if end - start >= self.min_event_length:
                segments.append((start, end))

        # Merge small segments
        merged_segments = self._merge_small_segments(segments, signal)

        return merged_segments

    def _merge_small_segments(
        self,
        segments: List[Tuple[int, int]],
        signal: np.ndarray
    ) -> List[Tuple[int, int]]:
        """
        Merge small segments that are too short.

        Args:
            segments: Initial segments
            signal: Signal array

        Returns:
            Merged segments
        """
        if len(segments) <= 1:
            return segments

        merged = [segments[0]]

        for i in range(1, len(segments)):
            prev_segment = merged[-1]
            curr_segment = segments[i]

            prev_length = prev_segment[1] - prev_segment[0]
            curr_length = curr_segment[1] - curr_segment[0]

            # If either segment is too short, merge them
            if (prev_length < self.min_event_length or
                curr_length < self.min_event_length):

                # Check if levels are similar (should be merged)
                prev_mean = np.mean(signal[prev_segment[0]:prev_segment[1]])
                curr_mean = np.mean(signal[curr_segment[0]:curr_segment[1]])

                # Merge if levels are close (within 2 std dev)
                combined_std = np.std(signal[prev_segment[0]:curr_segment[1]])
                if abs(prev_mean - curr_mean) < 2 * combined_std:
                    # Merge segments
                    merged[-1] = (prev_segment[0], curr_segment[1])
                    continue

            merged.append(curr_segment)

        return merged

    def _extract_events(
        self,
        signal: np.ndarray,
        boundaries: List[Tuple[int, int]]
    ) -> List[Event]:
        """
        Extract event statistics from signal segments.

        Args:
            signal: Signal array
            boundaries: Event boundaries [(start, end), ...]

        Returns:
            List of Event objects
        """
        events = []

        for i, (start, end) in enumerate(boundaries):
            event_signal = signal[start:end]

            if len(event_signal) == 0:
                continue

            # Calculate statistics
            mean = np.mean(event_signal)
            stdv = np.std(event_signal) if len(event_signal) > 1 else 0.0
            length = len(event_signal)

            # Create event
            event = Event(
                mean=mean,
                stdv=stdv,
                start=start,
                length=length,
                samples=event_signal.copy()
            )
            events.append(event)

        return events

    def _filter_events(self, events: List[Event]) -> List[Event]:
        """
        Filter out invalid or low-quality events.

        Args:
            events: Detected events

        Returns:
            Filtered events
        """
        filtered = []

        for event in events:
            # Check valid signal range
            if event.mean < 0 or event.mean > 200:
                continue

            # Check min duration
            if event.length < self.min_event_length:
                continue

            # Check reasonable standard deviation
            if event.stdv > event.mean * 0.5:
                continue

            filtered.append(event)

        return filtered

    def detect_events_hmm(
        self,
        raw_signal: np.ndarray,
        sampling_rate: float = 4000.0
    ) -> List[Event]:
        """
        Alternative event detection using Hidden Markov Model approach.

        This is a simplified HMM-based detector that models transitions
        between states (nucleotides).

        Args:
            raw_signal: Raw ionic current signal
            sampling_rate: Sampling rate in Hz (default: 4000 for MinION)

        Returns:
            List of detected events
        """
        # Preprocess
        smoothed = self._preprocess_signal(raw_signal)

        # Simple HMM with 3 states: same, up, down
        # This is a simplified version - a full HMM would use pore models
        states = []
        current_state = 0  # 0=same, 1=up, 2=down

        min_event_samples = int(sampling_rate * 0.002)  # 2ms minimum

        for i in range(1, len(smoothed)):
            diff = smoothed[i] - smoothed[i-1]

            if abs(diff) < 1.0:  # No significant change
                next_state = 0
            elif diff > 0:
                next_state = 1  # Level increase
            else:
                next_state = 2  # Level decrease

            if next_state != current_state:
                # Potential state change
                # Verify it's not noise by checking duration
                if len(states) > 0 and (i - states[-1]) < min_event_samples:
                    # Too short, stay in current state
                    continue

                # State transition
                current_state = next_state
                states.append(i)

        states.append(len(smoothed))  # End

        # Create boundaries
        boundaries = [(states[i], states[i+1]) for i in range(len(states)-1)]

        # Extract events
        events = self._extract_events(smoothed, boundaries)
        filtered = self._filter_events(events)

        return filtered
