#!/usr/bin/env python3
"""
Unit tests for fin._dtw module.

Tests the DTW (Dynamic Time Warping) functions:
- dtw_distance: Compute DTW distance between two sequences
- dtw_pairwise: Compute pairwise DTW distances for batch of sequences
- is_available: Check if CUDA DTW is available
- cleanup: Clean up CUDA resources

Note: These tests check both the API behavior and skip CUDA-specific tests
if CUDA is not available.

Run with:
    pytest tests/unit/test_dtw.py -v
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Check if DTW module is available
try:
    from fin._dtw import (
        dtw_distance,
        dtw_pairwise,
        cleanup,
        is_available,
        CUDA_AVAILABLE,
    )
    DTW_AVAILABLE = True
except ImportError:
    DTW_AVAILABLE = False
    CUDA_AVAILABLE = False


# Skip all tests if DTW module not available
pytestmark = pytest.mark.skipif(
    not DTW_AVAILABLE,
    reason="fin._dtw module not available (CUDA not installed)"
)


# ============================================================================
# Test fixtures
# ============================================================================

@pytest.fixture
def simple_sequences():
    """Create simple test sequences."""
    seq1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
    seq2 = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
    return seq1, seq2


@pytest.fixture
def different_sequences():
    """Create two different sequences."""
    seq1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
    seq2 = np.array([2.0, 3.0, 4.0, 5.0, 6.0], dtype=np.float32)
    return seq1, seq2


@pytest.fixture
def longer_sequences():
    """Create longer test sequences."""
    np.random.seed(42)
    seq1 = np.random.randn(100).astype(np.float32)
    seq2 = np.random.randn(100).astype(np.float32)
    return seq1, seq2


@pytest.fixture
def batch_sequences():
    """Create batch of sequences for pairwise computation."""
    np.random.seed(42)
    # 5 sequences of length 50
    return np.random.randn(5, 50).astype(np.float32)


# ============================================================================
# is_available tests
# ============================================================================

class TestIsAvailable:
    """Tests for is_available function."""

    def test_is_available_returns_bool(self):
        """Test that is_available returns a boolean."""
        result = is_available()
        assert isinstance(result, bool)

    def test_is_available_consistent(self):
        """Test that is_available returns consistent results."""
        result1 = is_available()
        result2 = is_available()
        assert result1 == result2


# ============================================================================
# dtw_distance tests (only run if CUDA available)
# ============================================================================

@pytest.mark.skipif(not CUDA_AVAILABLE, reason="CUDA not available")
class TestDtwDistance:
    """Tests for dtw_distance function."""

    def test_identical_sequences_zero_distance(self, simple_sequences):
        """Test that identical sequences have zero DTW distance."""
        seq1, seq2 = simple_sequences
        distance = dtw_distance(seq1, seq2)
        
        assert distance >= 0.0
        # Identical sequences should have very small distance
        assert distance < 1e-5

    def test_distance_non_negative(self, different_sequences):
        """Test that DTW distance is non-negative."""
        seq1, seq2 = different_sequences
        distance = dtw_distance(seq1, seq2)
        
        assert distance >= 0.0

    def test_distance_symmetric(self, different_sequences):
        """Test that DTW distance is symmetric."""
        seq1, seq2 = different_sequences
        d1 = dtw_distance(seq1, seq2)
        d2 = dtw_distance(seq2, seq1)
        
        assert abs(d1 - d2) < 1e-5

    def test_distance_returns_float(self, simple_sequences):
        """Test that dtw_distance returns a float."""
        seq1, seq2 = simple_sequences
        distance = dtw_distance(seq1, seq2)
        
        assert isinstance(distance, (float, np.floating))

    def test_distance_with_list_input(self):
        """Test dtw_distance accepts Python lists."""
        seq1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        seq2 = [1.0, 2.0, 3.0, 4.0, 5.0]
        
        distance = dtw_distance(seq1, seq2)
        assert distance >= 0.0

    def test_distance_converts_dtype(self):
        """Test that dtw_distance converts non-float32 arrays."""
        seq1 = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        seq2 = np.array([1.0, 2.0, 3.0], dtype=np.float64)
        
        distance = dtw_distance(seq1, seq2)
        assert distance >= 0.0

    def test_distance_with_open_start(self, different_sequences):
        """Test dtw_distance with open_start boundary."""
        seq1, seq2 = different_sequences
        
        d_closed = dtw_distance(seq1, seq2, use_open_start=False)
        d_open = dtw_distance(seq1, seq2, use_open_start=True)
        
        # Open start should allow skipping beginning, potentially lower distance
        assert d_open <= d_closed + 1e-5

    def test_distance_with_open_end(self, different_sequences):
        """Test dtw_distance with open_end boundary."""
        seq1, seq2 = different_sequences
        
        d_closed = dtw_distance(seq1, seq2, use_open_end=False)
        d_open = dtw_distance(seq1, seq2, use_open_end=True)
        
        # Open end should allow skipping end, potentially lower distance
        assert d_open <= d_closed + 1e-5

    def test_distance_longer_sequences(self, longer_sequences):
        """Test dtw_distance with longer sequences."""
        seq1, seq2 = longer_sequences
        distance = dtw_distance(seq1, seq2)
        
        assert distance >= 0.0
        assert np.isfinite(distance)


@pytest.mark.skipif(not CUDA_AVAILABLE, reason="CUDA not available")
class TestDtwDistanceValidation:
    """Tests for dtw_distance input validation."""

    def test_empty_sequence_raises(self):
        """Test that empty sequences raise error."""
        seq1 = np.array([], dtype=np.float32)
        seq2 = np.array([1.0, 2.0], dtype=np.float32)
        
        with pytest.raises(ValueError):
            dtw_distance(seq1, seq2)

    def test_2d_array_raises(self):
        """Test that 2D arrays raise error."""
        seq1 = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
        seq2 = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)
        
        with pytest.raises(ValueError):
            dtw_distance(seq1, seq2)


# ============================================================================
# dtw_pairwise tests (only run if CUDA available)
# ============================================================================

@pytest.mark.skipif(not CUDA_AVAILABLE, reason="CUDA not available")
class TestDtwPairwise:
    """Tests for dtw_pairwise function."""

    def test_pairwise_returns_square_matrix(self, batch_sequences):
        """Test that dtw_pairwise returns square distance matrix."""
        distances = dtw_pairwise(batch_sequences)
        
        n = batch_sequences.shape[0]
        assert distances.shape == (n, n)

    def test_pairwise_diagonal_near_zero(self, batch_sequences):
        """Test that diagonal elements are near zero."""
        distances = dtw_pairwise(batch_sequences)
        
        diagonal = np.diag(distances)
        assert np.allclose(diagonal, 0.0, atol=1e-5)

    def test_pairwise_symmetric(self, batch_sequences):
        """Test that distance matrix is symmetric."""
        distances = dtw_pairwise(batch_sequences)
        
        assert np.allclose(distances, distances.T, atol=1e-5)

    def test_pairwise_non_negative(self, batch_sequences):
        """Test that all distances are non-negative."""
        distances = dtw_pairwise(batch_sequences)
        
        assert np.all(distances >= -1e-5)

    def test_pairwise_with_open_boundaries(self, batch_sequences):
        """Test dtw_pairwise with open boundaries."""
        d_closed = dtw_pairwise(batch_sequences, use_open_start=False, use_open_end=False)
        d_open = dtw_pairwise(batch_sequences, use_open_start=True, use_open_end=True)
        
        # Both should be valid distance matrices
        assert d_closed.shape == d_open.shape
        assert np.all(d_closed >= -1e-5)
        assert np.all(d_open >= -1e-5)


# ============================================================================
# cleanup tests
# ============================================================================

@pytest.mark.skipif(not CUDA_AVAILABLE, reason="CUDA not available")
class TestCleanup:
    """Tests for cleanup function."""

    def test_cleanup_succeeds(self):
        """Test that cleanup runs without error."""
        # This should not raise an exception
        cleanup()

    def test_cleanup_idempotent(self):
        """Test that cleanup can be called multiple times."""
        cleanup()
        cleanup()
        cleanup()
        # Should not raise


# ============================================================================
# Module availability tests
# ============================================================================

class TestModuleAvailability:
    """Tests for module availability."""

    def test_cuda_available_constant(self):
        """Test that CUDA_AVAILABLE constant is boolean."""
        assert isinstance(CUDA_AVAILABLE, bool)

    def test_functions_exist(self):
        """Test that expected functions exist in module."""
        from fin import _dtw
        
        # Check function availability
        assert hasattr(_dtw, 'dtw_distance') or not CUDA_AVAILABLE
        assert hasattr(_dtw, 'dtw_pairwise') or not CUDA_AVAILABLE
        assert hasattr(_dtw, 'cleanup') or not CUDA_AVAILABLE
        assert hasattr(_dtw, 'is_available')


# ============================================================================
# Integration tests
# ============================================================================

@pytest.mark.skipif(not CUDA_AVAILABLE, reason="CUDA not available")
class TestDtwIntegration:
    """Integration tests for DTW module."""

    def test_pairwise_matches_individual(self, batch_sequences):
        """Test that pairwise distances match individual dtw_distance calls."""
        # Get pairwise distances
        pairwise = dtw_pairwise(batch_sequences)
        
        # Compare with individual calls
        for i in range(batch_sequences.shape[0]):
            for j in range(batch_sequences.shape[0]):
                individual = dtw_distance(batch_sequences[i], batch_sequences[j])
                # Allow small numerical differences
                assert abs(pairwise[i, j] - individual) < 0.1

    def test_typical_nanopore_signal_lengths(self):
        """Test with typical nanopore signal lengths."""
        np.random.seed(42)
        # Typical event-level signals might have 100-10000 samples
        seq1 = np.random.randn(1000).astype(np.float32)
        seq2 = np.random.randn(1000).astype(np.float32)
        
        distance = dtw_distance(seq1, seq2)
        
        assert np.isfinite(distance)
        assert distance >= 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
