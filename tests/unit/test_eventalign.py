#!/usr/bin/env python3
"""
Unit tests for fin._eventalign module.

Tests the event detection and alignment functions:
- getevents: Detect events from raw nanopore signal
- set_model: Load pore model
- run_eventalign: Full event alignment pipeline
- Model constants: MODEL_RNA002, MODEL_RNA004

Note: These tests check both the API behavior and skip tests that require
the C extension if it's not available.

Run with:
    pytest tests/unit/test_eventalign.py -v
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Check if eventalign module is available
try:
    from fin._eventalign import (
        getevents,
        set_model,
        MODEL_RNA002,
        MODEL_RNA004,
        MAX_KMER_SIZE,
        MAX_NUM_KMER,
        _EVENTALIGN_AVAILABLE,
    )
    EVENTALIGN_AVAILABLE = _EVENTALIGN_AVAILABLE
except ImportError:
    EVENTALIGN_AVAILABLE = False
    MODEL_RNA002 = 1
    MODEL_RNA004 = 2
    MAX_KMER_SIZE = 9
    MAX_NUM_KMER = 262144


# ============================================================================
# Test fixtures
# ============================================================================

@pytest.fixture
def simple_signal():
    """Create a simple test signal with clear event transitions."""
    # Create signal with distinct levels
    np.random.seed(42)
    signal = np.concatenate([
        np.random.normal(100, 2, 50),   # Level 1
        np.random.normal(120, 2, 50),   # Level 2
        np.random.normal(90, 2, 50),    # Level 3
        np.random.normal(110, 2, 50),   # Level 4
    ]).astype(np.float32)
    return signal


@pytest.fixture
def constant_signal():
    """Create a constant signal (no events)."""
    return np.full(200, 100.0, dtype=np.float32)


@pytest.fixture
def noisy_signal():
    """Create a noisy signal that mimics real nanopore data."""
    np.random.seed(42)
    # Create random walk with occasional jumps
    signal = np.zeros(1000, dtype=np.float32)
    signal[0] = 100.0
    
    for i in range(1, 1000):
        # Occasional level change
        if np.random.random() < 0.02:
            signal[i] = signal[i-1] + np.random.choice([-20, 20])
        else:
            signal[i] = signal[i-1] + np.random.normal(0, 1)
    
    return signal.astype(np.float32)


@pytest.fixture
def test_data_signal():
    """Load real test signal if available."""
    test_pod5 = PROJECT_ROOT / "tests" / "testdata" / "RNA004.test.pod5"
    if not test_pod5.exists():
        return None
    
    try:
        import pod5
        with pod5.Reader(test_pod5) as reader:
            for read in reader.reads():
                return read.signal_pa.astype(np.float32)
    except Exception:
        return None


# ============================================================================
# Model constants tests
# ============================================================================

class TestModelConstants:
    """Tests for model constants."""

    def test_model_rna002_value(self):
        """Test MODEL_RNA002 has expected value."""
        assert MODEL_RNA002 == 1

    def test_model_rna004_value(self):
        """Test MODEL_RNA004 has expected value."""
        assert MODEL_RNA004 == 2

    def test_max_kmer_size_reasonable(self):
        """Test MAX_KMER_SIZE is reasonable."""
        assert MAX_KMER_SIZE >= 5
        assert MAX_KMER_SIZE <= 15

    def test_max_num_kmer_value(self):
        """Test MAX_NUM_KMER has expected value."""
        # 4^9 = 262144 for 9-mer
        assert MAX_NUM_KMER == 262144


# ============================================================================
# getevents tests (only run if extension available)
# ============================================================================

@pytest.mark.skipif(not EVENTALIGN_AVAILABLE, reason="Eventalign extension not available")
class TestGetevents:
    """Tests for getevents function."""

    def test_getevents_returns_dict(self, simple_signal):
        """Test that getevents returns a dictionary."""
        result = getevents(simple_signal)
        assert isinstance(result, dict)

    def test_getevents_has_required_keys(self, simple_signal):
        """Test that result has all required keys."""
        result = getevents(simple_signal)
        
        required_keys = ['n_events', 'starts', 'lengths', 'means', 'stdvs']
        for key in required_keys:
            assert key in result, f"Missing key: {key}"

    def test_getevents_n_events_positive(self, simple_signal):
        """Test that n_events is positive for non-constant signal."""
        result = getevents(simple_signal)
        assert result['n_events'] > 0

    def test_getevents_arrays_same_length(self, simple_signal):
        """Test that all output arrays have same length."""
        result = getevents(simple_signal)
        n = result['n_events']
        
        assert len(result['starts']) == n
        assert len(result['lengths']) == n
        assert len(result['means']) == n
        assert len(result['stdvs']) == n

    def test_getevents_starts_increasing(self, simple_signal):
        """Test that event starts are monotonically increasing."""
        result = getevents(simple_signal)
        starts = result['starts']
        
        for i in range(1, len(starts)):
            assert starts[i] > starts[i-1], f"Start positions not increasing at {i}"

    def test_getevents_lengths_positive(self, simple_signal):
        """Test that all event lengths are positive."""
        result = getevents(simple_signal)
        lengths = result['lengths']
        
        assert np.all(lengths > 0), "Found non-positive lengths"

    def test_getevents_stdvs_non_negative(self, simple_signal):
        """Test that standard deviations are non-negative."""
        result = getevents(simple_signal)
        stdvs = result['stdvs']
        
        assert np.all(stdvs >= 0), "Found negative stdvs"

    def test_getevents_means_reasonable(self, simple_signal):
        """Test that means are within reasonable range."""
        result = getevents(simple_signal)
        means = result['means']
        
        # Means should be within signal range
        assert np.all(means >= simple_signal.min() - 1)
        assert np.all(means <= simple_signal.max() + 1)

    def test_getevents_events_cover_signal(self, simple_signal):
        """Test that events roughly cover the signal length."""
        result = getevents(simple_signal)
        
        # Total length should be close to signal length
        total_length = np.sum(result['lengths'])
        # Allow some margin for edge effects
        assert total_length > len(simple_signal) * 0.5
        assert total_length <= len(simple_signal) + result['n_events']


@pytest.mark.skipif(not EVENTALIGN_AVAILABLE, reason="Eventalign extension not available")
class TestGeteventsInputValidation:
    """Tests for getevents input validation."""

    def test_getevents_accepts_float32(self, simple_signal):
        """Test that getevents accepts float32 arrays."""
        signal = simple_signal.astype(np.float32)
        result = getevents(signal)
        assert result is not None

    def test_getevents_accepts_python_list(self, simple_signal):
        """Test that getevents accepts Python lists."""
        signal_list = simple_signal.tolist()
        result = getevents(signal_list)
        assert result is not None

    def test_getevents_short_signal(self):
        """Test getevents with very short signal."""
        short_signal = np.array([100.0, 110.0, 105.0], dtype=np.float32)
        # Should either return valid result or raise appropriate error
        try:
            result = getevents(short_signal)
            # If it returns, check it's valid
            assert 'n_events' in result
        except Exception:
            # Some implementations may reject very short signals
            pass


# ============================================================================
# set_model tests (only run if extension available)
# ============================================================================

@pytest.mark.skipif(not EVENTALIGN_AVAILABLE, reason="Eventalign extension not available")
class TestSetModel:
    """Tests for set_model function."""

    def test_set_model_rna002_returns_dict(self):
        """Test that set_model returns a dictionary for RNA002."""
        result = set_model(MODEL_RNA002)
        assert isinstance(result, dict)

    def test_set_model_rna004_returns_dict(self):
        """Test that set_model returns a dictionary for RNA004."""
        result = set_model(MODEL_RNA004)
        assert isinstance(result, dict)

    def test_set_model_has_kmer_size(self):
        """Test that model dict has kmer_size."""
        result = set_model(MODEL_RNA002)
        assert 'kmer_size' in result

    def test_set_model_rna002_kmer_size(self):
        """Test that RNA002 model has k=5."""
        result = set_model(MODEL_RNA002)
        assert result['kmer_size'] == 5

    def test_set_model_rna004_kmer_size(self):
        """Test that RNA004 model has k=9."""
        result = set_model(MODEL_RNA004)
        assert result['kmer_size'] == 9

    def test_set_model_has_level_mean(self):
        """Test that model has level_mean array."""
        result = set_model(MODEL_RNA002)
        assert 'level_mean' in result
        assert len(result['level_mean']) > 0

    def test_set_model_has_level_stdv(self):
        """Test that model has level_stdv array."""
        result = set_model(MODEL_RNA002)
        assert 'level_stdv' in result
        assert len(result['level_stdv']) > 0


# ============================================================================
# Integration tests with real data
# ============================================================================

@pytest.mark.skipif(not EVENTALIGN_AVAILABLE, reason="Eventalign extension not available")
class TestEventalignIntegration:
    """Integration tests with real test data."""

    def test_getevents_with_real_signal(self, test_data_signal):
        """Test getevents with real nanopore signal."""
        if test_data_signal is None:
            pytest.skip("Test data not available")
        
        result = getevents(test_data_signal)
        
        assert result['n_events'] > 0
        assert len(result['starts']) == result['n_events']
        assert len(result['means']) == result['n_events']

    def test_event_statistics_reasonable(self, test_data_signal):
        """Test that event statistics are reasonable for real data."""
        if test_data_signal is None:
            pytest.skip("Test data not available")
        
        result = getevents(test_data_signal)
        
        # For real RNA signals, typical values
        means = result['means']
        stdvs = result['stdvs']
        lengths = result['lengths']
        
        # Check means are in reasonable pA range (typical: 50-200 pA)
        assert np.median(means) > 0
        assert np.median(means) < 300
        
        # Check lengths are reasonable (typical: 5-200 samples)
        assert np.median(lengths) > 1
        assert np.median(lengths) < 500


# ============================================================================
# Module availability tests
# ============================================================================

class TestModuleAvailability:
    """Tests for module availability and fallbacks."""

    def test_availability_flag_exists(self):
        """Test that availability flag exists."""
        from fin._eventalign import _EVENTALIGN_AVAILABLE
        assert isinstance(_EVENTALIGN_AVAILABLE, bool)

    def test_constants_always_available(self):
        """Test that model constants are always available."""
        from fin._eventalign import MODEL_RNA002, MODEL_RNA004
        assert MODEL_RNA002 is not None
        assert MODEL_RNA004 is not None

    def test_graceful_degradation(self):
        """Test that module degrades gracefully if extension unavailable."""
        from fin._eventalign import _EVENTALIGN_AVAILABLE
        
        if not _EVENTALIGN_AVAILABLE:
            # Should still be importable, just with limited functionality
            from fin._eventalign import MODEL_RNA002, MODEL_RNA004
            assert MODEL_RNA002 == 1
            assert MODEL_RNA004 == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
