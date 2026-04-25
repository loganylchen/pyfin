"""GPU/CPU equivalence tests for EM quantification (US-019).

These tests verify that CPU and GPU code paths produce numerically equivalent
results from em_with_coherence. All tests are skipped when no GPU / CuPy is
available so CI environments without CUDA continue to pass.
"""

from __future__ import annotations

import numpy as np
import pytest

try:
    import cupy as cp

    GPU_AVAILABLE = bool(cp.cuda.is_available())
except Exception:
    GPU_AVAILABLE = False


@pytest.mark.skipif(not GPU_AVAILABLE, reason="No GPU / CuPy available")
def test_gpu_cpu_equivalence_em() -> None:
    """AC: CPU and GPU em_with_coherence converge to the same R within atol=1e-6."""
    from fin.analysis.assignments import em_with_coherence

    rng = np.random.default_rng(42)
    n_reads, n_tx = 20, 5
    dist_tx = rng.random((n_reads, n_tx)) * 5
    dist_rr = rng.random((n_reads, n_reads))
    dist_rr = (dist_rr + dist_rr.T) / 2
    np.fill_diagonal(dist_rr, 0)

    R_cpu, _, _ = em_with_coherence(
        dist_tx,
        dist_rr,
        sigma=1.0,
        beta=0.5,
        max_iter=50,
        tol=1e-6,
        verbose=False,
        use_gpu=False,
    )
    R_gpu, _, _ = em_with_coherence(
        dist_tx,
        dist_rr,
        sigma=1.0,
        beta=0.5,
        max_iter=50,
        tol=1e-6,
        verbose=False,
        use_gpu=True,
    )

    np.testing.assert_allclose(R_cpu, R_gpu, atol=1e-6, rtol=1e-4)


@pytest.mark.skipif(not GPU_AVAILABLE, reason="No GPU / CuPy available")
def test_gpu_cpu_equivalence_with_prior() -> None:
    """Same prior_weights must produce equivalent R on CPU vs GPU."""
    from fin.analysis.assignments import em_with_coherence

    rng = np.random.default_rng(7)
    n_reads, n_tx = 15, 3
    dist_tx = rng.random((n_reads, n_tx))
    dist_rr = rng.random((n_reads, n_reads))
    dist_rr = (dist_rr + dist_rr.T) / 2
    np.fill_diagonal(dist_rr, 0)
    prior = np.array([3.0, 1.0, 2.0])

    R_cpu, _, _ = em_with_coherence(
        dist_tx,
        dist_rr,
        sigma=1.0,
        beta=0.5,
        max_iter=30,
        tol=1e-6,
        verbose=False,
        use_gpu=False,
        prior_weights=prior,
    )
    R_gpu, _, _ = em_with_coherence(
        dist_tx,
        dist_rr,
        sigma=1.0,
        beta=0.5,
        max_iter=30,
        tol=1e-6,
        verbose=False,
        use_gpu=True,
        prior_weights=prior,
    )

    np.testing.assert_allclose(R_cpu, R_gpu, atol=1e-6, rtol=1e-4)
