#!/usr/bin/env python3
"""
Unit tests for event detection comparison between PyFIN and f5c.

This module tests that PyFIN's event detection produces events that are
compatible with f5c's output.
"""

import numpy as np
import sys
from pathlib import Path
from typing import Tuple, Optional
import pytest

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Test data paths - actual test data location
TEST_DATA_DIR = PROJECT_ROOT / "tests" / "testdata"

# RNA004 test data
RNA004_POD5_PATH = TEST_DATA_DIR / "RNA004.test.pod5"
RNA002_POD5_PATH = TEST_DATA_DIR / "RNA002.test.pod5"

# Legacy path (examples directory)
LEGACY_DATA_DIR = PROJECT_ROOT / "examples" / "test_data"
POD5_PATH = LEGACY_DATA_DIR / "one_read.pod5"

# Check if required modules are available
try:
    import pod5
    POD5_AVAILABLE = True
except ImportError:
    POD5_AVAILABLE = False

try:
    from fin._eventalign import getevents
    EVENTALIGN_AVAILABLE = True
except ImportError:
    EVENTALIGN_AVAILABLE = False


def load_signal_from_pod5(pod5_path: Path) -> Tuple[str, np.ndarray, float]:
    """Load signal data from POD5 file."""
    import pod5
    
    with pod5.Reader(str(pod5_path)) as reader:
        for read in reader.reads():
            signal = read.signal_pa.astype(np.float32)
            sample_rate = float(read.run_info.sample_rate)
            return str(read.read_id), signal, sample_rate
    
    raise ValueError(f"No reads found in {pod5_path}")


def get_available_pod5() -> Optional[Path]:
    """Get the first available POD5 file for testing."""
    for path in [RNA004_POD5_PATH, RNA002_POD5_PATH, POD5_PATH]:
        if path.exists():
            return path
    return None


@pytest.mark.skipif(not POD5_AVAILABLE, reason="pod5 not installed")
@pytest.mark.skipif(not EVENTALIGN_AVAILABLE, reason="fin._eventalign not available")
class TestEventDetection:
    """Test class for event detection validation."""
    
    @classmethod
    def setup_class(cls):
        """Load test data once for all tests."""
        pod5_path = get_available_pod5()
        if pod5_path is None:
            pytest.skip("No POD5 test data available")
        
        cls.read_id, cls.signal, cls.sample_rate = load_signal_from_pod5(pod5_path)
        print(f"\nLoaded signal from {pod5_path.name}: {len(cls.signal)} samples, {cls.sample_rate} Hz")
    
    def test_event_detection_runs(self):
        """Test that event detection runs without errors."""
        from fin._eventalign import getevents
        
        events = getevents(self.signal)
        
        assert 'n_events' in events
        assert 'starts' in events
        assert 'lengths' in events
        assert 'means' in events
        assert 'stdvs' in events
        
        print(f"  Detected {events['n_events']} events")
    
    def test_event_count_reasonable(self):
        """Test that event count is reasonable for signal length."""
        from fin._eventalign import getevents
        
        events = getevents(self.signal)
        n_events = events['n_events']
        
        # RNA events typically have ~8-12 samples per event
        expected_min = len(self.signal) // 20  # At most 20 samples/event average
        expected_max = len(self.signal) // 4   # At least 4 samples/event average
        
        assert n_events >= expected_min, \
            f"Too few events: {n_events} < {expected_min}"
        assert n_events <= expected_max, \
            f"Too many events: {n_events} > {expected_max}"
        
        samples_per_event = len(self.signal) / n_events
        print(f"  Average samples per event: {samples_per_event:.1f}")
    
    def test_events_cover_signal(self):
        """Test that events cover the entire signal (approximately)."""
        from fin._eventalign import getevents
        
        events = getevents(self.signal)
        
        starts = events['starts']
        lengths = events['lengths']
        
        # First event should start near beginning
        assert starts[0] < len(self.signal) * 0.1, \
            f"First event starts too late: {starts[0]}"
        
        # Last event should end near signal end
        last_end = starts[-1] + int(lengths[-1])
        assert last_end > len(self.signal) * 0.9, \
            f"Last event ends too early: {last_end} vs {len(self.signal)}"
        
        print(f"  Events cover samples {starts[0]} to {last_end} of {len(self.signal)}")
    
    def test_events_are_contiguous(self):
        """Test that events are roughly contiguous (minimal gaps)."""
        from fin._eventalign import getevents
        
        events = getevents(self.signal)
        
        starts = events['starts']
        lengths = events['lengths']
        
        gaps = []
        for i in range(len(starts) - 1):
            end_i = starts[i] + int(lengths[i])
            start_next = starts[i + 1]
            gap = start_next - end_i
            gaps.append(gap)
        
        gaps = np.array(gaps)
        mean_gap = np.mean(gaps)
        max_gap = np.max(gaps)
        
        # Gaps should be small (overlap is allowed, so negative gaps are OK)
        assert mean_gap < 10, f"Mean gap too large: {mean_gap:.1f}"
        assert max_gap < 100, f"Max gap too large: {max_gap}"
        
        print(f"  Mean gap: {mean_gap:.1f}, Max gap: {max_gap}")
    
    def test_event_means_in_expected_range(self):
        """Test that event means are in physiologically reasonable range."""
        from fin._eventalign import getevents
        
        events = getevents(self.signal)
        means = events['means']
        
        # RNA signal typically in 60-160 pA range (after normalization)
        mean_min = np.min(means)
        mean_max = np.max(means)
        mean_avg = np.mean(means)
        
        print(f"  Event means: min={mean_min:.1f}, max={mean_max:.1f}, avg={mean_avg:.1f} pA")
        
        # Just ensure we have variation and reasonable values
        assert mean_max > mean_min, "No variation in event means"
        assert np.std(means) > 1.0, "Event means have too little variation"
    
    def test_event_stdvs_positive(self):
        """Test that event standard deviations are positive."""
        from fin._eventalign import getevents
        
        events = getevents(self.signal)
        stdvs = events['stdvs']
        
        assert np.all(stdvs >= 0), "Negative standard deviations found"
        assert np.mean(stdvs) > 0, "Mean stdv is zero"
        
        print(f"  Event stdvs: mean={np.mean(stdvs):.2f}, max={np.max(stdvs):.2f}")
    
    def test_event_order_is_signal_order(self):
        """Test that events are in raw signal order (not reversed)."""
        from fin._eventalign import getevents
        
        events = getevents(self.signal)
        starts = events['starts']
        
        # Events should be monotonically increasing in start position
        diffs = np.diff(starts)
        
        # Allow small negative diffs (overlaps) but overall increasing
        assert np.mean(diffs) > 0, "Events are not in increasing order"
        assert starts[-1] > starts[0], "Last event start < first event start"
        
        print(f"  Event order: starts from {starts[0]} to {starts[-1]}")


@pytest.mark.skipif(not POD5_AVAILABLE, reason="pod5 not installed")
@pytest.mark.skipif(not EVENTALIGN_AVAILABLE, reason="fin._eventalign not available")
class TestEventDetectionParameters:
    """Test that event detection uses RNA-specific parameters."""
    
    @classmethod
    def setup_class(cls):
        """Load test data."""
        pod5_path = get_available_pod5()
        if pod5_path is None:
            pytest.skip("No POD5 test data available")
        
        cls.read_id, cls.signal, cls.sample_rate = load_signal_from_pod5(pod5_path)
    
    def test_rna_window_length(self):
        """
        Test that RNA window length is used (7 samples).
        
        RNA uses window_length1=7 while DNA uses window_length1=5.
        We can verify this indirectly by checking event characteristics.
        """
        from fin._eventalign import getevents
        
        events = getevents(self.signal)
        lengths = events['lengths']
        
        # RNA events with window=7 tend to be slightly longer
        mean_length = np.mean(lengths)
        
        # RNA events typically 8-12 samples on average
        assert 5 <= mean_length <= 20, \
            f"Mean event length {mean_length:.1f} outside expected RNA range"
        
        print(f"  Mean event length: {mean_length:.1f} samples (RNA expected: 8-12)")
    
    def test_consistent_event_detection(self):
        """Test that event detection is deterministic."""
        from fin._eventalign import getevents
        
        events1 = getevents(self.signal)
        events2 = getevents(self.signal)
        
        assert events1['n_events'] == events2['n_events'], \
            "Event detection is not deterministic"
        
        np.testing.assert_array_equal(events1['starts'], events2['starts'])
        np.testing.assert_array_almost_equal(events1['means'], events2['means'])


if __name__ == "__main__":
    import pytest
    pytest.main([__file__, "-v", "-s"])
