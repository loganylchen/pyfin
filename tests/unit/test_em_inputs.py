"""Tests for fin.scoring.em_inputs.build_em_matrices."""

from __future__ import annotations

import numpy as np
import pytest

from fin.scoring.em_inputs import VALID_SUBSETS, build_em_matrices


def _mats(n_reads=3, n_tx=4):
    rng = np.random.default_rng(42)
    m1 = rng.uniform(0, 5, size=(n_reads, n_tx)).astype(np.float32)
    m2 = rng.uniform(0, 5, size=(n_reads, n_tx)).astype(np.float32)
    return m1, m2


class TestSubsetBehavior:
    def test_m1_only_uses_m1_and_zeros_rr(self):
        m1, m2 = _mats()
        d_tx, d_rr, beta = build_em_matrices("m1", m1, m2)
        np.testing.assert_array_equal(d_tx, m1)
        np.testing.assert_array_equal(
            d_rr, np.zeros((m1.shape[0], m1.shape[0]), dtype=np.float32)
        )
        assert beta == 0.0

    def test_m2_only_uses_m2_and_zeros_rr(self):
        m1, m2 = _mats()
        d_tx, d_rr, beta = build_em_matrices("m2", m1, m2)
        np.testing.assert_array_equal(d_tx, m2)
        np.testing.assert_array_equal(
            d_rr, np.zeros((m1.shape[0], m1.shape[0]), dtype=np.float32)
        )
        assert beta == 0.0

    def test_m1_plus_m2_zscore_mean_nonnegative(self):
        m1, m2 = _mats()
        d_tx, d_rr, beta = build_em_matrices("m1+m2", m1, m2)
        assert d_tx.min() >= 0
        np.testing.assert_array_equal(
            d_rr, np.zeros((m1.shape[0], m1.shape[0]), dtype=np.float32)
        )
        assert beta == 0.0
        assert d_tx.min(axis=1).max() < 1e-5


class TestValidation:
    def test_unknown_subset_raises(self):
        m1, m2 = _mats()
        with pytest.raises(ValueError, match="em_matrix_subset"):
            build_em_matrices("bogus", m1, m2)

    def test_all_subsets_are_coherence_free(self):
        m1, m2 = _mats()
        for subset in VALID_SUBSETS:
            d_tx, d_rr, beta = build_em_matrices(subset, m1, m2)
            assert d_tx.min() >= 0
            assert d_rr.min() >= 0
            assert beta == 0.0

    def test_output_dtype_float32(self):
        m1, m2 = _mats()
        d_tx, d_rr, _ = build_em_matrices("m1+m2", m1, m2)
        assert d_tx.dtype == np.float32
        assert d_rr.dtype == np.float32


class TestZscoreFusion:
    def test_zscore_mean_handles_sentinel(self):
        m1 = np.array([[0.0, 1.0, 1e6]], dtype=np.float32)
        m2 = np.array([[0.0, 2.0, 4.0]], dtype=np.float32)
        d_tx, _, _ = build_em_matrices("m1+m2", m1, m2)
        assert np.isfinite(d_tx).all()
        assert d_tx.min() >= 0
