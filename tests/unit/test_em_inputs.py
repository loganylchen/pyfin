"""Tests for fin.scoring.em_inputs.build_em_matrices."""

from __future__ import annotations

import numpy as np
import pytest

from fin.scoring.em_inputs import VALID_SUBSETS, build_em_matrices


def _mats(n_reads=3, n_tx=4):
    rng = np.random.default_rng(42)
    m1 = rng.uniform(0, 5, size=(n_reads, n_tx)).astype(np.float32)
    m2 = rng.uniform(0, 5, size=(n_reads, n_tx)).astype(np.float32)
    m3 = rng.uniform(0, 5, size=(n_reads, n_reads)).astype(np.float32)
    return m1, m2, m3


class TestSubsetBehavior:
    def test_m1_only_uses_m1_and_zeros_rr(self):
        m1, m2, m3 = _mats()
        d_tx, d_rr, beta = build_em_matrices("m1", m1, m2, m3, em_beta=0.5)
        np.testing.assert_array_equal(d_tx, m1)
        np.testing.assert_array_equal(d_rr, np.zeros_like(m3))
        assert beta == 0.0

    def test_m2_only_uses_m2_and_zeros_rr(self):
        m1, m2, m3 = _mats()
        d_tx, d_rr, beta = build_em_matrices("m2", m1, m2, m3, em_beta=0.5)
        np.testing.assert_array_equal(d_tx, m2)
        np.testing.assert_array_equal(d_rr, np.zeros_like(m3))
        assert beta == 0.0

    def test_m3_only_zeros_tx_keeps_beta(self):
        m1, m2, m3 = _mats()
        d_tx, d_rr, beta = build_em_matrices("m3", m1, m2, m3, em_beta=0.5)
        assert d_tx.shape == m2.shape
        np.testing.assert_array_equal(d_tx, np.zeros_like(m2))
        np.testing.assert_array_equal(d_rr, m3)
        assert beta == 0.5

    def test_m1_plus_m2_zscore_mean_nonneg(self):
        m1, m2, m3 = _mats()
        d_tx, d_rr, beta = build_em_matrices("m1+m2", m1, m2, m3, em_beta=0.5)
        assert d_tx.min() >= 0
        np.testing.assert_array_equal(d_rr, np.zeros_like(m3))
        assert beta == 0.0
        # row-min ≈ 0 after shift
        assert d_tx.min(axis=1).max() < 1e-5

    def test_m1_plus_m3(self):
        m1, m2, m3 = _mats()
        d_tx, d_rr, beta = build_em_matrices("m1+m3", m1, m2, m3, em_beta=0.7)
        np.testing.assert_array_equal(d_tx, m1)
        np.testing.assert_array_equal(d_rr, m3)
        assert beta == 0.7

    def test_m2_plus_m3_default(self):
        m1, m2, m3 = _mats()
        d_tx, d_rr, beta = build_em_matrices("m2+m3", m1, m2, m3, em_beta=0.7)
        np.testing.assert_array_equal(d_tx, m2)
        np.testing.assert_array_equal(d_rr, m3)
        assert beta == 0.7

    def test_all_uses_zscore_mean_and_m3(self):
        m1, m2, m3 = _mats()
        d_tx, d_rr, beta = build_em_matrices("all", m1, m2, m3, em_beta=0.3)
        assert d_tx.min() >= 0
        np.testing.assert_array_equal(d_rr, m3)
        assert beta == 0.3


class TestValidation:
    def test_unknown_subset_raises(self):
        m1, m2, m3 = _mats()
        with pytest.raises(ValueError, match="em_matrix_subset"):
            build_em_matrices("bogus", m1, m2, m3, em_beta=0.5)

    def test_all_subsets_covered(self):
        m1, m2, m3 = _mats()
        for s in VALID_SUBSETS:
            d_tx, d_rr, beta = build_em_matrices(s, m1, m2, m3, em_beta=0.5)
            assert d_tx.min() >= 0
            assert d_rr.min() >= 0
            assert beta in (0.0, 0.5)

    def test_output_dtype_float32(self):
        m1, m2, m3 = _mats()
        d_tx, d_rr, _ = build_em_matrices("all", m1, m2, m3, em_beta=0.5)
        assert d_tx.dtype == np.float32
        assert d_rr.dtype == np.float32


class TestZscoreFusion:
    def test_zscore_mean_handles_sentinel(self):
        # 1e6 sentinel should be treated as missing.
        m1 = np.array([[0.0, 1.0, 1e6]], dtype=np.float32)
        m2 = np.array([[0.0, 2.0, 4.0]], dtype=np.float32)
        m3 = np.zeros((1, 1), dtype=np.float32)
        d_tx, _, _ = build_em_matrices("m1+m2", m1, m2, m3, em_beta=0.0)
        assert np.isfinite(d_tx).all()
        assert d_tx.min() >= 0
