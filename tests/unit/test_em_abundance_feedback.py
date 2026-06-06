"""Unit tests for RSEM/Salmon-style abundance feedback in ``em_with_coherence``.

Covers ``fin.analysis.assignments.em_with_coherence`` with the new
``abundance_feedback`` / ``eff_lengths`` / ``abundance_length_norm`` /
``abundance_eps`` parameters:
  - feedback OFF is identical to not passing the new kwargs (production safety).
  - feedback ON pulls a read shared between two candidates toward the more
    abundant one (the 0.5/0.5 -> 0.99/0.01 behaviour the user asked for).
  - non-identifiable abundances (all reads equidistant) stay 0.5/0.5 even ON.
  - effective-length normalization down-weights long transcripts.
  - theta underflow is eps-guarded (no NaN/inf).
  - the EM still converges with feedback on.
  - CPU/GPU parity (marked ``gpu``).
And a config-default test (feedback OFF by default).
"""
from __future__ import annotations

import numpy as np
import pytest

from fin.analysis.assignments import CUPY_AVAILABLE, em_with_coherence


def _shared_read_scenario():
    """9 reads that strongly prefer candidate A (col 0), plus 1 read that is
    exactly equidistant to A and B (col 1). Read-to-read is zero (beta=0)."""
    n_reads, n_tx = 10, 2
    d_tx = np.zeros((n_reads, n_tx), dtype=np.float64)
    d_tx[:9] = [0.0, 10.0]   # reads 0..8: clearly candidate A
    d_tx[9] = [1.0, 1.0]     # read 9: shared / equidistant
    d_rr = np.zeros((n_reads, n_reads), dtype=np.float64)
    return d_tx, d_rr


# --------------------------------------------------------------------------
# OFF path == defaults (production safety)
# --------------------------------------------------------------------------
def test_off_is_byte_identical():
    d_tx, d_rr = _shared_read_scenario()
    R_implicit, h_implicit, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False
    )
    R_explicit, h_explicit, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False,
        abundance_feedback=False, abundance_length_norm=False, eff_lengths=None,
    )
    assert np.array_equal(R_implicit, R_explicit)
    assert np.array_equal(h_implicit, h_explicit)


def test_off_shared_read_stays_split():
    d_tx, d_rr = _shared_read_scenario()
    R, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False,
        abundance_feedback=False,
    )
    # The shared read (row 9) must remain ~0.5/0.5 without feedback.
    assert R[9, 0] == pytest.approx(0.5, abs=1e-3)
    assert R[9, 1] == pytest.approx(0.5, abs=1e-3)


# --------------------------------------------------------------------------
# ON path: abundance pull
# --------------------------------------------------------------------------
def test_feedback_pulls_shared_read():
    d_tx, d_rr = _shared_read_scenario()
    R, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False,
        abundance_feedback=True,
    )
    # With ~9.5 vs ~0.5 mass, the shared read is pulled strongly toward A.
    assert R[9, 0] > 0.9
    assert R[9, 0] + R[9, 1] == pytest.approx(1.0, abs=1e-9)


def test_non_identifiability_stays_split():
    # Every read equidistant to both candidates -> abundances non-identifiable.
    n = 8
    d_tx = np.full((n, 2), 3.0, dtype=np.float64)
    d_rr = np.zeros((n, n), dtype=np.float64)
    R, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False,
        abundance_feedback=True,
    )
    # theta stays uniform -> split stays 0.5/0.5 (documented RSEM behaviour).
    assert np.allclose(R, 0.5, atol=1e-6)


# --------------------------------------------------------------------------
# effective-length normalization
# --------------------------------------------------------------------------
def test_eff_length_norm_path():
    # 5 reads prefer A, 5 prefer B, 1 shared -> raw counts ~equal (theta~0.5/0.5).
    n_reads = 11
    d_tx = np.zeros((n_reads, 2), dtype=np.float64)
    d_tx[:5] = [0.0, 10.0]    # prefer A
    d_tx[5:10] = [10.0, 0.0]  # prefer B
    d_tx[10] = [1.0, 1.0]     # shared
    d_rr = np.zeros((n_reads, n_reads), dtype=np.float64)

    R_plain, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False,
        abundance_feedback=True, abundance_length_norm=False,
    )
    R_norm, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False,
        abundance_feedback=True, abundance_length_norm=True,
        eff_lengths=np.array([1.0, 10.0]),  # B is 10x longer -> down-weighted
    )
    # Raw counts ~equal -> shared read stays ~0.5 without length norm.
    assert R_plain[10, 0] == pytest.approx(0.5, abs=0.05)
    # Length norm penalizes the long candidate B -> shared read tilts to A.
    assert R_norm[10, 0] > R_plain[10, 0]
    assert R_norm[10, 0] > 0.5


def test_eff_length_norm_requires_eff_lengths():
    d_tx, d_rr = _shared_read_scenario()
    with pytest.raises(ValueError, match="requires eff_lengths"):
        em_with_coherence(
            d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False,
            abundance_feedback=True, abundance_length_norm=True, eff_lengths=None,
        )


def test_eff_length_wrong_shape_rejected():
    d_tx, d_rr = _shared_read_scenario()
    with pytest.raises(ValueError, match="eff_lengths must be"):
        em_with_coherence(
            d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False,
            abundance_feedback=True, abundance_length_norm=True,
            eff_lengths=np.array([1.0, 2.0, 3.0]),  # n_tx is 2
        )


def test_eff_length_nonpositive_rejected():
    d_tx, d_rr = _shared_read_scenario()
    with pytest.raises(ValueError, match="finite and strictly positive"):
        em_with_coherence(
            d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False,
            abundance_feedback=True, abundance_length_norm=True,
            eff_lengths=np.array([1.0, 0.0]),
        )


# --------------------------------------------------------------------------
# numerical robustness + convergence
# --------------------------------------------------------------------------
def test_theta_underflow_no_nan():
    # All reads overwhelmingly prefer A -> theta_B underflows toward 0.
    n = 12
    d_tx = np.tile([0.0, 1000.0], (n, 1)).astype(np.float64)
    d_rr = np.zeros((n, n), dtype=np.float64)
    R, hard, lls = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False,
        abundance_feedback=True,
    )
    assert np.all(np.isfinite(R))
    assert np.all(np.isfinite(lls))
    assert np.all(hard == 0)


def test_convergence_with_feedback():
    d_tx, d_rr = _shared_read_scenario()
    R, _, lls = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.0, max_iter=1000, tol=1e-4,
        verbose=False, use_gpu=False, abundance_feedback=True,
    )
    assert len(lls) < 1000          # converged before hitting the cap
    assert np.all(np.isfinite(lls))
    assert np.allclose(R.sum(axis=1), 1.0)  # rows are valid distributions


# --------------------------------------------------------------------------
# GPU/CPU parity
# --------------------------------------------------------------------------
@pytest.mark.gpu
def test_gpu_cpu_parity():
    if not CUPY_AVAILABLE:
        pytest.skip("CuPy not available")
    d_tx, d_rr = _shared_read_scenario()
    R_cpu, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=False,
        abundance_feedback=True,
    )
    R_gpu, _, _ = em_with_coherence(
        d_tx, d_rr, sigma=1.0, beta=0.0, verbose=False, use_gpu=True,
        abundance_feedback=True,
    )
    assert np.allclose(R_cpu, R_gpu, atol=1e-5)


# --------------------------------------------------------------------------
# config default
# --------------------------------------------------------------------------
def test_config_defaults_off():
    from fin.pipeline.config import PipelineConfig

    cfg = PipelineConfig(bam_path="x.bam")
    assert cfg.abundance_feedback is False
    assert cfg.abundance_length_norm is False
