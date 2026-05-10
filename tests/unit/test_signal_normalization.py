"""Tests for per-read robust z-score normalization in signal_dtw.

Verifies that with `normalize=True`, two signals that differ only by an
affine offset+scale produce a much smaller DTW distance than when raw.
"""

from __future__ import annotations

import numpy as np

from fin.scoring.eventalign_parser import ReadCandidateScore
from fin.scoring.signal_dtw import _cpu_dtw, extract_signal_segments


class _FakeSignalReader:
    """Minimal stub matching the slow5/pod5 reader interface used by signal_dtw."""

    def __init__(self, signals):
        self._signals = signals

    def get_picoamp_signal(self, read_name):
        sig = self._signals.get(read_name)
        if sig is None:
            return None
        return sig, {}

    def get_calibrated_signal(self, read_name):
        return self.get_picoamp_signal(read_name)


def _shape(n=200, seed=0):
    rng = np.random.default_rng(seed)
    base = np.cumsum(rng.normal(0, 1, size=n)).astype(np.float32)
    # Smooth it a touch so DTW has structure to follow.
    return base


def test_normalization_collapses_offset_scale_drift():
    """Two reads with same signal *shape* but different offset/scale should
    have near-zero DTW after normalization, much larger DTW without.
    """
    base = _shape(n=200, seed=7)
    sig_a = base.copy()
    sig_b = (base * 5.0 + 100.0).astype(np.float32)  # affine transform

    scores = [
        ReadCandidateScore(
            read_name="r_a",
            candidate_id="c1",
            total_log_likelihood=-1.0,
            signal_start_idx=0,
            signal_end_idx=len(sig_a),
        ),
        ReadCandidateScore(
            read_name="r_b",
            candidate_id="c1",
            total_log_likelihood=-1.0,
            signal_start_idx=0,
            signal_end_idx=len(sig_b),
        ),
    ]
    reader = _FakeSignalReader({"r_a": sig_a, "r_b": sig_b})

    norm_segs = extract_signal_segments(scores, reader, signal_format="slow5", normalize=True)
    raw_segs = extract_signal_segments(scores, reader, signal_format="slow5", normalize=False)

    d_norm = _cpu_dtw(norm_segs["r_a"], norm_segs["r_b"])
    d_raw = _cpu_dtw(raw_segs["r_a"], raw_segs["r_b"])

    # Normalized DTW should be dramatically smaller (orders of magnitude).
    assert d_norm < d_raw, f"normalized DTW {d_norm} not smaller than raw {d_raw}"
    assert d_norm < d_raw * 0.05, (
        f"normalization should remove most affine drift; "
        f"got d_norm={d_norm:.2f}, d_raw={d_raw:.2f}"
    )


def test_normalization_disabled_passes_through_raw_signal():
    """With normalize=False, segments should match raw input slice exactly."""
    sig = np.array([10.0, 20.0, 30.0, 40.0], dtype=np.float32)
    score = ReadCandidateScore(
        read_name="r",
        candidate_id="c",
        total_log_likelihood=0.0,
        signal_start_idx=0,
        signal_end_idx=4,
    )
    reader = _FakeSignalReader({"r": sig})

    segs = extract_signal_segments([score], reader, signal_format="slow5", normalize=False)
    np.testing.assert_array_equal(segs["r"], sig)


def test_normalization_handles_constant_signal():
    """A constant segment has MAD=0; normalization must not divide by zero."""
    sig = np.full(50, 7.5, dtype=np.float32)
    score = ReadCandidateScore(
        read_name="r",
        candidate_id="c",
        total_log_likelihood=0.0,
        signal_start_idx=0,
        signal_end_idx=50,
    )
    reader = _FakeSignalReader({"r": sig})

    segs = extract_signal_segments([score], reader, signal_format="slow5", normalize=True)
    out = segs["r"]
    assert np.all(np.isfinite(out))
    # All same value -> median equals value -> z-scores are all zero.
    assert np.allclose(out, 0.0)
