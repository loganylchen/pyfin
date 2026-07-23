"""Unit tests for 5'-TSS short-isoform recovery (peak-vs-ramp)."""
from __future__ import annotations

from fin.candidates.isoform_recovery import (
    recover_5p_peaks,
    evaluate_tss_peak,
)


def _sharp_peak(pos, n, flank_n):
    """n reads exactly at pos (in +-window), flank_n scattered in the +-3w flanks."""
    ends = [(pos, 1.0) for _ in range(n)]
    for i in range(flank_n):
        ends.append((pos + 40 + i, 1.0))   # 40 > window(20), <= 60 (3w) -> flank
    return ends


def test_sharp_peak_accepted():
    ends = _sharp_peak(1000, 30, 4)
    r = evaluate_tss_peak(ends, 1000)
    assert r is not None
    assert r.support == 30
    assert r.oe > 4 and r.concentration > 0.6
    assert r.excess > 25             # peak mass well above flank background


def test_flat_ramp_rejected():
    # reads spread evenly every 5 bp across +-60 -> low concentration, low O/E.
    ends = [(940 + 5 * k, 1.0) for k in range(25)]   # 940..1060
    assert evaluate_tss_peak(ends, 1000) is None


def test_below_support_floor_rejected():
    # a clean, concentrated peak but only 5 reads (< min_reads 8) -> rejected.
    ends = _sharp_peak(1000, 5, 0)
    assert evaluate_tss_peak(ends, 1000) is None


def test_min_frac_floor_scales_with_family_mass():
    # 10 reads in-window would pass the absolute floor (>=8), but with a large family
    # mass the 3% fractional floor (0.03 * 400 = 12) rejects it.
    ends = _sharp_peak(1000, 10, 2)
    small = evaluate_tss_peak(ends, 1000, total_mass=12)     # 3% -> 0.36, floor=8
    assert small is not None
    big = evaluate_tss_peak(ends, 1000, total_mass=400)      # 3% -> 12 > 10
    assert big is None


def test_concentration_gate_rejects_broad_bump():
    # equal mass inside and outside +-window -> concentration 0.5 < 0.6 -> rejected
    # even though absolute counts are high.
    ends = [(1000, 1.0) for _ in range(20)] + [(1040 + i, 1.0) for i in range(20)]
    r = evaluate_tss_peak(ends, 1000)
    assert r is None


def test_recover_dedups_near_proposals():
    # two proposals 10 bp apart both land on the same peak -> merged to one.
    ends = _sharp_peak(1000, 30, 4)
    peaks = recover_5p_peaks(ends, [1000, 1010])
    assert len(peaks) == 1
    assert abs(peaks[0].pos - 1000) <= 10


def test_recover_keeps_distinct_peaks():
    ends = _sharp_peak(1000, 30, 2) + _sharp_peak(5000, 25, 2)
    peaks = recover_5p_peaks(ends, [1000, 5000])
    assert {p.pos for p in peaks} == {1000, 5000}


def test_recover_empty_when_no_proposal_passes():
    ends = [(940 + 5 * k, 1.0) for k in range(25)]
    assert recover_5p_peaks(ends, [1000]) == []
