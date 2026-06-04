"""Tests for shared helpers in fin/scoring/composite.py."""

from __future__ import annotations

import numpy as np
import pytest

from fin.scoring.composite import (
    derive_prior_weights,
    populate_quant_scores,
    subsample_reads_for_dtw,
    CompositeScore,
)


# ---------------------------------------------------------------------------
# derive_prior_weights
# ---------------------------------------------------------------------------


def test_derive_prior_weights_returns_none_for_zero_tx():
    assert derive_prior_weights(np.array([]), n_tx=0, cap=10.0) is None


def test_derive_prior_weights_n_tx_one_collapses_to_singleton():
    pw = derive_prior_weights(np.array([0.7]), n_tx=1, cap=10.0)
    assert pw is not None
    assert pw.shape == (1,)
    assert pw[0] == pytest.approx(1.0)


def test_derive_prior_weights_all_zero_scores_yields_uniform():
    """All-zero scores must NOT produce NaNs; expect (near) uniform output."""
    pw = derive_prior_weights(np.zeros(4), n_tx=4, cap=10.0)
    assert pw is not None
    assert pw.shape == (4,)
    assert np.all(np.isfinite(pw))
    assert pw.sum() == pytest.approx(1.0)
    np.testing.assert_allclose(pw, np.full(4, 0.25))


def test_derive_prior_weights_respects_cap():
    """A dominant score must be capped relative to uniform.

    After re-normalization the lower bound can drift below uniform/cap, but the
    upper bound is always honored: no candidate may exceed uniform*cap.
    """
    cap = 5.0
    n_tx = 4
    uniform = 1.0 / n_tx
    # One score is heavily dominant.
    scores = np.array([1.0, 0.001, 0.001, 0.001])
    pw = derive_prior_weights(scores, n_tx=n_tx, cap=cap)
    assert pw.sum() == pytest.approx(1.0)
    assert np.all(pw > 0)
    assert np.all(pw <= uniform * cap + 1e-9)
    # Without clipping, the dominant weight would be ~0.997; clipping keeps it
    # at most uniform * cap = 1.25 -> after renorm strictly below 1.
    assert pw[0] < 0.95


def test_derive_prior_weights_clamps_cap_below_one():
    """cap < 1.0 must be clamped to 1.0 (= no shaping)."""
    pw = derive_prior_weights(np.array([0.5, 0.5]), n_tx=2, cap=0.1)
    np.testing.assert_allclose(pw, np.full(2, 0.5))


# ---------------------------------------------------------------------------
# populate_quant_scores
# ---------------------------------------------------------------------------


def test_populate_quant_scores_writes_fields():
    from fin.analysis.quantification import QuantResult

    qr_a = QuantResult(
        candidate_id="A", abundance=1.0, confidence=0.0,
        num_assigned_reads=0, source="gtf",
    )
    qr_b = QuantResult(
        candidate_id="B", abundance=2.0, confidence=0.0,
        num_assigned_reads=0, source="gtf",
    )
    scores = [
        CompositeScore(candidate_id="A", coherence=0.1, discrimination=0.2, combined=0.3),
        CompositeScore(candidate_id="B", coherence=0.4, discrimination=0.5, combined=0.6),
    ]
    populate_quant_scores([qr_a, qr_b], scores)
    assert qr_a.coherence_score == 0.1
    assert qr_a.discrimination_score == 0.2
    assert qr_a.combined_score == 0.3
    assert qr_b.coherence_score == 0.4
    assert qr_b.discrimination_score == 0.5
    assert qr_b.combined_score == 0.6


def test_populate_quant_scores_skips_unmatched_ids():
    from fin.analysis.quantification import QuantResult

    qr = QuantResult(
        candidate_id="A", abundance=1.0, confidence=0.0,
        num_assigned_reads=0, source="gtf",
    )
    populate_quant_scores(
        [qr],
        [CompositeScore(candidate_id="Z", coherence=1.0, discrimination=1.0, combined=1.0)],
    )
    # Defaults preserved when no composite score matches.
    assert qr.coherence_score == 0.0
    assert qr.discrimination_score == 0.0
    assert qr.combined_score == 0.0


# ---------------------------------------------------------------------------
# subsample_reads_for_dtw
# ---------------------------------------------------------------------------


def test_subsample_reads_for_dtw_below_cap_returns_input():
    reads = ["r0", "r1", "r2"]
    assert subsample_reads_for_dtw(reads, max_dtw=10) is reads


def test_subsample_reads_for_dtw_zero_or_none_disables():
    reads = [f"r{i}" for i in range(50)]
    assert subsample_reads_for_dtw(reads, max_dtw=None) is reads
    assert subsample_reads_for_dtw(reads, max_dtw=0) is reads


def test_subsample_reads_for_dtw_uniform_indices():
    reads = [f"r{i:03d}" for i in range(100)]
    sub = subsample_reads_for_dtw(reads, max_dtw=10)
    assert len(sub) == 10
    # Endpoints are preserved by np.linspace.
    assert sub[0] == "r000"
    assert sub[-1] == "r099"
    # Strictly increasing original ordering.
    assert sub == sorted(sub)
