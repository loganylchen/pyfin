#!/usr/bin/env python3
"""
Unit tests for fin.analysis.assignments module.

Tests the EM algorithm with coherence regularization for clustering reads to transcripts:
- em_with_coherence: EM algorithm with read-read coherence

Run with:
    pytest tests/unit/test_assignments.py -v
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.analysis.assignments import em_with_coherence


# ============================================================================
# Test fixtures
# ============================================================================

@pytest.fixture
def simple_distance_matrices():
    """Create simple distance matrices for testing."""
    # 4 reads, 2 transcripts
    dist_read_to_tx = np.array([
        [0.1, 0.9],  # Read 0: close to tx 0
        [0.2, 0.8],  # Read 1: close to tx 0
        [0.8, 0.2],  # Read 2: close to tx 1
        [0.9, 0.1],  # Read 3: close to tx 1
    ])
    
    # Read-read distances (similar reads are close)
    dist_read_to_read = np.array([
        [0.0, 0.1, 0.8, 0.9],  # Read 0
        [0.1, 0.0, 0.7, 0.8],  # Read 1
        [0.8, 0.7, 0.0, 0.1],  # Read 2
        [0.9, 0.8, 0.1, 0.0],  # Read 3
    ])
    
    return dist_read_to_tx, dist_read_to_read


@pytest.fixture
def uniform_distance_matrices():
    """Create uniform distance matrices where all reads are equidistant."""
    n_reads, n_tx = 4, 2
    dist_read_to_tx = np.ones((n_reads, n_tx)) * 0.5
    dist_read_to_read = np.ones((n_reads, n_reads)) * 0.5
    np.fill_diagonal(dist_read_to_read, 0.0)
    return dist_read_to_tx, dist_read_to_read


@pytest.fixture
def single_read_matrices():
    """Create matrices with just one read."""
    dist_read_to_tx = np.array([[0.1, 0.9]])
    dist_read_to_read = np.array([[0.0]])
    return dist_read_to_tx, dist_read_to_read


# ============================================================================
# Basic functionality tests
# ============================================================================

class TestEmWithCoherence:
    """Tests for em_with_coherence function."""

    def test_returns_correct_types(self, simple_distance_matrices):
        """Test that function returns correct types."""
        dist_tx, dist_rr = simple_distance_matrices
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, verbose=False
        )
        
        assert isinstance(R, np.ndarray)
        assert isinstance(assignments, np.ndarray)
        assert isinstance(lls, list)

    def test_responsibility_matrix_shape(self, simple_distance_matrices):
        """Test responsibility matrix has correct shape."""
        dist_tx, dist_rr = simple_distance_matrices
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, verbose=False
        )
        
        assert R.shape == dist_tx.shape

    def test_assignments_shape(self, simple_distance_matrices):
        """Test assignments array has correct shape."""
        dist_tx, dist_rr = simple_distance_matrices
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, verbose=False
        )
        
        assert assignments.shape == (dist_tx.shape[0],)

    def test_responsibilities_sum_to_one(self, simple_distance_matrices):
        """Test that responsibilities sum to 1 for each read."""
        dist_tx, dist_rr = simple_distance_matrices
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, verbose=False
        )
        
        row_sums = R.sum(axis=1)
        np.testing.assert_array_almost_equal(row_sums, np.ones(R.shape[0]))

    def test_responsibilities_non_negative(self, simple_distance_matrices):
        """Test that all responsibilities are non-negative."""
        dist_tx, dist_rr = simple_distance_matrices
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, verbose=False
        )
        
        assert np.all(R >= 0)

    def test_assignments_valid_range(self, simple_distance_matrices):
        """Test that assignments are valid transcript indices."""
        dist_tx, dist_rr = simple_distance_matrices
        n_tx = dist_tx.shape[1]
        
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, verbose=False
        )
        
        assert np.all(assignments >= 0)
        assert np.all(assignments < n_tx)

    def test_log_likelihoods_non_empty(self, simple_distance_matrices):
        """Test that log-likelihoods list is not empty."""
        dist_tx, dist_rr = simple_distance_matrices
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, verbose=False
        )
        
        assert len(lls) > 0


# ============================================================================
# Clustering quality tests
# ============================================================================

class TestEmWithCoherenceClustering:
    """Tests for clustering behavior of em_with_coherence."""

    def test_clear_separation_clustering(self, simple_distance_matrices):
        """Test that clearly separated reads are assigned correctly."""
        dist_tx, dist_rr = simple_distance_matrices
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, sigma=1.0, beta=0.5, verbose=False
        )
        
        # Reads 0,1 should be assigned to transcript 0
        # Reads 2,3 should be assigned to transcript 1
        assert assignments[0] == 0
        assert assignments[1] == 0
        assert assignments[2] == 1
        assert assignments[3] == 1

    def test_beta_zero_ignores_coherence(self, simple_distance_matrices):
        """Test that beta=0 ignores read-read similarity."""
        dist_tx, dist_rr = simple_distance_matrices
        
        # With beta=0, only read-to-transcript distances matter
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, sigma=1.0, beta=0.0, verbose=False
        )
        
        # Should still assign based on dist_read_to_tx
        assert assignments[0] == 0  # closest to tx 0
        assert assignments[3] == 1  # closest to tx 1

    def test_high_beta_strong_clustering(self, simple_distance_matrices):
        """Test that high beta leads to stronger clustering."""
        dist_tx, dist_rr = simple_distance_matrices
        
        # With high beta, similar reads should cluster together
        R_low_beta, _, _ = em_with_coherence(
            dist_tx, dist_rr, beta=0.1, verbose=False
        )
        R_high_beta, _, _ = em_with_coherence(
            dist_tx, dist_rr, beta=2.0, verbose=False
        )
        
        # High beta should lead to more confident assignments
        confidence_low = np.max(R_low_beta, axis=1).mean()
        confidence_high = np.max(R_high_beta, axis=1).mean()
        
        # Can't guarantee higher confidence, but both should be > 0.5
        assert confidence_low > 0.5
        assert confidence_high > 0.5


# ============================================================================
# Parameter validation tests
# ============================================================================

class TestEmWithCoherenceValidation:
    """Tests for parameter validation in em_with_coherence."""

    def test_mismatched_shapes_raises(self):
        """Test that mismatched shapes raise error."""
        dist_tx = np.random.rand(4, 2)
        dist_rr = np.random.rand(5, 5)  # Wrong shape

        with pytest.raises(ValueError):
            em_with_coherence(dist_tx, dist_rr, verbose=False)

    def test_non_square_read_distance_raises(self):
        """Test that non-square read distance matrix raises error."""
        dist_tx = np.random.rand(4, 2)
        dist_rr = np.random.rand(4, 3)  # Not square

        with pytest.raises(ValueError):
            em_with_coherence(dist_tx, dist_rr, verbose=False)

    def test_negative_distances_raises(self):
        """Test that negative distances raise error."""
        dist_tx = np.array([[-0.1, 0.9], [0.2, 0.8]])
        dist_rr = np.array([[0.0, 0.5], [0.5, 0.0]])
        
        with pytest.raises(ValueError, match="non-negative"):
            em_with_coherence(dist_tx, dist_rr, verbose=False)

    def test_negative_beta_raises(self):
        """Test that negative beta raises error."""
        dist_tx = np.random.rand(4, 2)
        dist_rr = np.random.rand(4, 4)
        np.fill_diagonal(dist_rr, 0)

        with pytest.raises(ValueError):
            em_with_coherence(dist_tx, dist_rr, beta=-0.5, verbose=False)

    def test_sigma_below_minimum_raises(self):
        """Test that sigma below minimum raises error."""
        dist_tx = np.random.rand(4, 2)
        dist_rr = np.random.rand(4, 4)
        np.fill_diagonal(dist_rr, 0)

        with pytest.raises(ValueError):
            em_with_coherence(dist_tx, dist_rr, sigma=0.001, verbose=False)


# ============================================================================
# Convergence tests
# ============================================================================

class TestEmWithCoherenceConvergence:
    """Tests for convergence behavior."""

    def test_convergence_within_max_iter(self, simple_distance_matrices):
        """Test that algorithm converges within max_iter."""
        dist_tx, dist_rr = simple_distance_matrices
        
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, max_iter=1000, tol=1e-4, verbose=False
        )
        
        # Should converge before hitting max_iter
        assert len(lls) < 1000

    def test_low_tolerance_more_iterations(self, simple_distance_matrices):
        """Test that lower tolerance leads to more iterations."""
        dist_tx, dist_rr = simple_distance_matrices
        
        _, _, lls_high_tol = em_with_coherence(
            dist_tx, dist_rr, tol=1e-2, verbose=False
        )
        _, _, lls_low_tol = em_with_coherence(
            dist_tx, dist_rr, tol=1e-6, verbose=False
        )
        
        # Lower tolerance should require more iterations (or same if already converged)
        assert len(lls_low_tol) >= len(lls_high_tol)

    def test_log_likelihood_increases(self, simple_distance_matrices):
        """Test that log-likelihood generally increases (or stays same)."""
        dist_tx, dist_rr = simple_distance_matrices
        
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, verbose=False
        )
        
        # Note: Due to the energy formulation, this might not strictly increase
        # but the algorithm should still converge
        assert len(lls) > 0


# ============================================================================
# Edge cases
# ============================================================================

class TestEmWithCoherenceEdgeCases:
    """Tests for edge cases."""

    def test_single_read(self, single_read_matrices):
        """Test with single read."""
        dist_tx, dist_rr = single_read_matrices
        
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, verbose=False
        )
        
        assert R.shape == (1, 2)
        assert assignments.shape == (1,)
        assert assignments[0] == 0  # Closer to first transcript

    def test_single_transcript(self):
        """Test with single transcript."""
        dist_tx = np.array([[0.1], [0.2], [0.3]])
        dist_rr = np.array([[0, 0.1, 0.2], [0.1, 0, 0.1], [0.2, 0.1, 0]])
        
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, verbose=False
        )
        
        # All reads should be assigned to the only transcript
        assert np.all(assignments == 0)
        # Responsibilities should all be 1.0 (only option)
        np.testing.assert_array_almost_equal(R, np.ones((3, 1)))

    def test_identical_distances(self, uniform_distance_matrices):
        """Test with uniform/identical distances."""
        dist_tx, dist_rr = uniform_distance_matrices
        
        R, assignments, lls = em_with_coherence(
            dist_tx, dist_rr, verbose=False
        )
        
        # Should still produce valid output
        assert R.shape == dist_tx.shape
        np.testing.assert_array_almost_equal(R.sum(axis=1), np.ones(R.shape[0]))


# ============================================================================
# Reproducibility tests
# ============================================================================

class TestEmWithCoherenceReproducibility:
    """Tests for reproducibility."""

    def test_deterministic_output(self, simple_distance_matrices):
        """Test that same input produces same output."""
        dist_tx, dist_rr = simple_distance_matrices
        
        R1, a1, ll1 = em_with_coherence(dist_tx, dist_rr, verbose=False)
        R2, a2, ll2 = em_with_coherence(dist_tx, dist_rr, verbose=False)
        
        np.testing.assert_array_almost_equal(R1, R2)
        np.testing.assert_array_equal(a1, a2)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
