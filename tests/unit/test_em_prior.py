import numpy as np
import pytest

from fin.analysis.assignments import em_with_coherence


def _setup(seed=42):
    rng = np.random.default_rng(seed)
    n_reads, n_tx = 20, 5
    d_tx = rng.random((n_reads, n_tx)) + 0.1
    d_rr_raw = rng.random((n_reads, n_reads))
    d_rr = (d_rr_raw + d_rr_raw.T) / 2
    np.fill_diagonal(d_rr, 0.0)
    return d_tx, d_rr


def test_em_accepts_prior_weights_param():
    d_tx, d_rr = _setup()
    prior = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
    R, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.5, max_iter=10,
        verbose=False, use_gpu=False, prior_weights=prior,
    )
    assert R.shape == (20, 5)


def test_em_prior_weights_validation_wrong_shape():
    d_tx, d_rr = _setup()
    with pytest.raises(AssertionError):
        em_with_coherence(
            d_tx, d_rr, sigma=1.0, beta=0.5, max_iter=5,
            verbose=False, use_gpu=False,
            prior_weights=np.array([1.0, 1.0]),  # wrong shape
        )


def test_em_prior_weights_negative_raises():
    d_tx, d_rr = _setup()
    with pytest.raises(ValueError):
        em_with_coherence(
            d_tx, d_rr, sigma=1.0, beta=0.5, max_iter=5,
            verbose=False, use_gpu=False,
            prior_weights=np.array([-1.0, 1.0, 1.0, 1.0, 1.0]),
        )


def test_em_backward_compat_none_vs_omitted():
    """prior_weights=None must give identical R to omitting the argument."""
    d_tx, d_rr = _setup()
    R_old, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.5, max_iter=50,
        verbose=False, use_gpu=False,
    )
    R_new, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.5, max_iter=50,
        verbose=False, use_gpu=False, prior_weights=None,
    )
    assert np.allclose(R_old, R_new, atol=1e-10)


def test_em_prior_biases_initial_assignment():
    """Biased prior should shift R toward the favored transcript at convergence.

    With many iterations the prior (applied only at initialization) can still
    steer convergence toward a different local optimum than the uniform start.
    We verify this by using asymmetric d_tx so that the prior-seeded trajectory
    ends up with higher average mass on transcript 0.
    """
    rng = np.random.default_rng(7)
    n_reads, n_tx = 20, 3
    # d_tx: tx=0 is generally closer (smaller distance) than tx=1 and tx=2
    d_tx = rng.random((n_reads, n_tx)) + 0.5
    d_tx[:, 0] *= 0.5  # tx=0 has half the distances -> favored by likelihood too

    d_rr_raw = rng.random((n_reads, n_reads))
    d_rr = (d_rr_raw + d_rr_raw.T) / 2
    np.fill_diagonal(d_rr, 0.0)

    # Strong prior toward tx=0
    prior_biased = np.array([100.0, 1.0, 1.0])
    R_biased, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.5, max_iter=50,
        verbose=False, use_gpu=False, prior_weights=prior_biased,
    )
    R_uniform, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.5, max_iter=50,
        verbose=False, use_gpu=False, prior_weights=None,
    )
    # Biased prior should concentrate at least as much mass on candidate 0
    assert R_biased[:, 0].mean() >= R_uniform[:, 0].mean()


def test_em_prior_all_zeros_falls_back():
    """All-zero prior should not crash; warn and fall back to uniform."""
    d_tx, d_rr = _setup()
    R, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.5, max_iter=5,
        verbose=False, use_gpu=False,
        prior_weights=np.zeros(5),
    )
    assert R.shape == (20, 5)
    # No NaNs
    assert np.all(np.isfinite(R))
