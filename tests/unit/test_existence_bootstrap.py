"""Calibration properties of the existence null (spec D3 audit gate), pinned at
a Monte-Carlo scale small enough for a unit test. The full study lives in
experiments/prod_validation/existence_null_calibration.py."""
from __future__ import annotations

import warnings

import numpy as np
import pytest

from fin.analysis.existence_bootstrap import (
    NULL_METHOD, balanced_folds, bootstrap_component, by_adjust, crossfit_batched,
    fit_mixture, holm_adjust, predictive_ll, required_replicates,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)


def _h0(rng, n, k, sparsity, noise):
    W = rng.normal(0, noise, (n, k))
    W[rng.random((n, k)) < sparsity] = -np.inf
    empty = ~np.isfinite(W).any(axis=1)
    W[empty, 0] = 0.0
    return W


class TestMixture:
    def test_fit_recovers_weights(self):
        rng = np.random.default_rng(0)
        n = 4000
        z = rng.random(n) < 0.7                      # 70% from col 0
        W = np.zeros((n, 2))
        W[z, 1] = -3.0; W[~z, 0] = -3.0              # each read strongly prefers its origin
        theta, ok = fit_mixture(W)
        assert ok and abs(theta[0] - 0.7) < 0.03

    def test_absent_column_gets_zero_weight(self):
        W = np.array([[0.0, -np.inf], [0.0, -np.inf]])
        theta, _ = fit_mixture(W)
        assert theta[1] == 0.0 and theta[0] == pytest.approx(1.0)

    def test_predictive_ll_marginalises(self):
        W = np.array([[0.0, -1.0]])
        theta = np.array([0.5, 0.5])
        assert predictive_ll(W, theta)[0] == pytest.approx(np.log(0.5 + 0.5 * np.exp(-1)))


class TestFolds:
    def test_balanced_and_nonempty(self):
        rng = np.random.default_rng(1)
        for n in (2, 3, 7, 40):
            f, k = balanced_folds(n, 5, rng)
            assert k == min(5, n)
            counts = np.bincount(f, minlength=k)
            assert counts.min() >= 1 and counts.max() - counts.min() <= 1
        assert balanced_folds(1, 5, rng)[1] == 0


class TestCrossfit:
    def test_indistinguishable_columns_score_zero(self):
        W = np.zeros((30, 2))
        folds, k = balanced_folds(30, 5, np.random.default_rng(0))
        T = crossfit_batched(W[None], np.ones((1, 30), bool), folds[None], k, [0, 1])
        assert abs(T[0][0]) < 1e-9 and abs(T[1][0]) < 1e-9

    def test_systematic_advantage_scores_positive(self):
        W = np.zeros((30, 2)); W[:, 1] = -1.0
        folds, k = balanced_folds(30, 5, np.random.default_rng(0))
        T = crossfit_batched(W[None], np.ones((1, 30), bool), folds[None], k, [0, 1])
        assert T[0][0] > 0 and T[1][0] <= 1e-9

    def test_ineligible_reads_do_not_score(self):
        W = np.zeros((30, 2)); W[:, 1] = -1.0
        folds, k = balanced_folds(30, 5, np.random.default_rng(0))
        T = crossfit_batched(W[None], np.zeros((1, 30), bool), folds[None], k, [0])
        assert T[0][0] == 0.0


class TestNull:
    def test_h0_is_calibrated_or_conservative(self):
        rng = np.random.default_rng(7)
        ps, pm = [], []
        for s in range(40):
            W = _h0(rng, 30, 3, 0.3, 0.5)
            r = bootstrap_component(W, np.ones(30, bool), [0], max_folds=5,
                                    replicates=99, seed=s, component_tag="h0")
            if r.degenerate:
                continue
            ps.append(r.p_plus[0]); pm.append(r.p_minus[0])
        ps, pm = np.array(ps), np.array(pm)
        assert len(ps) >= 30
        # 40 sets: at nominal 0.10 the MC 2-sigma band is ~+0.10
        assert (ps <= 0.10).mean() <= 0.20 and (pm <= 0.10).mean() <= 0.20

    def test_h1_has_power(self):
        rng = np.random.default_rng(3)
        hits = 0
        for s in range(10):
            W = np.zeros((20, 2)); W[:, 1] = -1.0 + rng.normal(0, 0.3, 20)
            r = bootstrap_component(W, np.ones(20, bool), [0], max_folds=5,
                                    replicates=199, seed=s, component_tag="h1")
            hits += (not r.degenerate) and r.p_plus[0] <= 0.025
        assert hits >= 8

    def test_duplicate_of_alternative_is_degenerate_not_significant(self):
        rng = np.random.default_rng(2)
        W = np.zeros((20, 2)); W[:, 0] = rng.normal(0, 0.2, 20); W[:, 1] = W[:, 0]
        r = bootstrap_component(W, np.ones(20, bool), [0, 1], max_folds=5,
                                replicates=99, seed=0, component_tag="dup")
        assert r.degenerate and r.p_plus == {}

    def test_single_molecule_is_degenerate(self):
        W = np.zeros((1, 2)); W[:, 1] = -5.0
        r = bootstrap_component(W, np.ones(1, bool), [0], max_folds=5,
                                replicates=49, seed=0, component_tag="t")
        assert r.degenerate and r.p_plus == {}

    def test_very_small_components_cannot_reach_significance(self):
        # n=2/3 are resampleable (LOO), but the permutation null has too few
        # configurations for a small p: the verdict must stay unresolved.
        for n in (2, 3):
            W = np.zeros((n, 2)); W[:, 1] = -5.0
            r = bootstrap_component(W, np.ones(n, bool), [0], max_folds=5,
                                    replicates=199, seed=0, component_tag="t")
            assert r.degenerate or r.p_plus[0] > 0.025

    def test_deterministic(self):
        W = np.zeros((15, 2)); W[:, 1] = -0.7 + np.random.default_rng(5).normal(0, 0.3, 15)
        a = bootstrap_component(W, np.ones(15, bool), [0], max_folds=5, replicates=63, seed=9, component_tag="x")
        b = bootstrap_component(W, np.ones(15, bool), [0], max_folds=5, replicates=63, seed=9, component_tag="x")
        assert a.p_plus == b.p_plus and a.t_obs == b.t_obs
        assert NULL_METHOD == "label_permutation_molecule_bootstrap"


class TestCorrections:
    def test_required_replicates_grows_with_family(self):
        assert required_replicates(1, 1) == 4095
        assert required_replicates(40, 6) == 8191
        assert required_replicates(400, 6) > 8191

    def test_by_and_holm(self):
        q = by_adjust([0.001, 0.02, 0.2, 0.6])
        assert all(b >= a - 1e-12 for a, b in zip(q, q[1:]))
        assert q[0] > 0.001                                # penalised
        assert holm_adjust([0.01, 0.02, 0.5]) == pytest.approx([0.03, 0.04, 0.5])
        assert by_adjust([]) == [] and holm_adjust([]) == []
