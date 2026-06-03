"""Unit tests for the full-length end-coherence filter.

Covers ``fin.analysis.quantification``:
  - ``compute_fulllen_frac`` — strand-aware full-length read fraction over a
    candidate's assigned reads (unreachable -> -1.0 sentinel).
  - ``fulllen_fraction_drops`` — the production drop set (FLAIR/TALON-style
    full-length read support; orthogonal to ``isoform_fraction_drops``).
  - ``aggregate_across_intervals`` — threads ``fulllen_frac`` through (max of
    reachable values, recall-safe).
"""
from __future__ import annotations

from fin.analysis.quantification import (
    QuantResult,
    aggregate_across_intervals,
    compute_fulllen_frac,
    fulllen_fraction_drops,
)


def _qr(
    cid,
    *,
    source="novel",
    chrom="chr1",
    strand="+",
    start=0,
    end=1000,
    exons=None,
    assigned=(),
    fulllen_frac=-1.0,
    abundance=1.0,
):
    if exons is None:
        exons = ((start, 300), (500, 700), (800, end))  # 3 exons => 2 introns
    return QuantResult(
        candidate_id=cid,
        abundance=abundance,
        confidence=1.0,
        num_assigned_reads=len(assigned),
        source=source,
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
        exons=exons,
        assigned_read_ids=tuple(assigned),
        fulllen_frac=fulllen_frac,
    )


def _drops(results, thr):
    return fulllen_fraction_drops({qr.candidate_id: qr for qr in results}, thr)


# --------------------------------------------------------------------------
# compute_fulllen_frac
# --------------------------------------------------------------------------
def test_compute_all_full_length_plus_strand():
    qr = _qr("c", start=100, end=900, assigned=("r1", "r2", "r3", "r4"))
    ends = {
        "r1": (100, 900),
        "r2": (110, 890),  # within 25 of both ends
        "r3": (90, 910),
        "r4": (105, 895),
    }
    assert compute_fulllen_frac(qr, ends, window=25, min_reads=4) == 1.0


def test_compute_partial_full_length():
    qr = _qr("c", start=100, end=900, assigned=("r1", "r2", "r3", "r4"))
    ends = {
        "r1": (100, 900),   # full
        "r2": (100, 900),   # full
        "r3": (500, 900),   # 5' off by 400 -> not full
        "r4": (100, 500),   # 3' off by 400 -> not full
    }
    assert compute_fulllen_frac(qr, ends, window=25, min_reads=4) == 0.5


def test_compute_minus_strand_uses_swapped_ends():
    # On '-' strand the genomic 3' end is start, genomic 5' end is end.
    qr = _qr("c", strand="-", start=100, end=900, assigned=("r1", "r2"))
    ends = {
        "r1": (100, 900),  # 5'=900, 3'=100 vs cand 5'=900, 3'=100 -> full
        "r2": (102, 898),  # within 25 -> full
    }
    assert compute_fulllen_frac(qr, ends, window=25, min_reads=2) == 1.0


def test_compute_unreachable_returns_sentinel():
    qr = _qr("c", start=100, end=900, assigned=("r1", "r2", "r3", "r4"))
    ends = {"r1": (100, 900)}  # only 1 read has a span; min_reads=4
    assert compute_fulllen_frac(qr, ends, window=25, min_reads=4) == -1.0


def test_compute_reads_without_span_are_not_counted():
    qr = _qr("c", start=100, end=900, assigned=("r1", "r2", "r3"))
    ends = {"r1": (100, 900), "r2": (100, 900)}  # r3 has no span
    # n_used = 2 >= min_reads 2 -> 2/2 = 1.0
    assert compute_fulllen_frac(qr, ends, window=25, min_reads=2) == 1.0


def test_compute_window_boundary_inclusive():
    qr = _qr("c", start=100, end=900, assigned=("r1", "r2"))
    ends = {
        "r1": (125, 875),  # exactly 25 off both ends -> inclusive, full
        "r2": (126, 900),  # 26 off 5' -> not full
    }
    assert compute_fulllen_frac(qr, ends, window=25, min_reads=2) == 0.5


# --------------------------------------------------------------------------
# fulllen_fraction_drops
# --------------------------------------------------------------------------
def test_drop_below_threshold():
    low = _qr("low", fulllen_frac=0.05)
    high = _qr("high", fulllen_frac=0.80)
    assert _drops([low, high], 0.1) == {"low"}


def test_at_threshold_is_kept():
    qr = _qr("c", fulllen_frac=0.1)  # 0.1 < 0.1 is False -> kept
    assert _drops([qr], 0.1) == set()


def test_unreachable_sentinel_never_dropped():
    qr = _qr("c", fulllen_frac=-1.0)
    assert _drops([qr], 0.1) == set()


def test_zero_fulllen_frac_is_dropped():
    qr = _qr("c", fulllen_frac=0.0)
    assert _drops([qr], 0.1) == {"c"}


def test_gtf_and_fusion_exempt():
    g = _qr("g", source="gtf", fulllen_frac=0.0)
    f = _qr("f", source="fusion", fulllen_frac=0.0)
    assert _drops([g, f], 0.1) == set()


def test_mono_and_single_intron_exempt():
    mono = _qr("mono", exons=((0, 1000),), fulllen_frac=0.0)
    two_exon = _qr("two", exons=((0, 400), (600, 1000)), fulllen_frac=0.0)
    # both have < 3 exons (< 2 introns) -> exempt
    assert _drops([mono, two_exon], 0.1) == set()


def test_zero_threshold_disables_filter():
    qr = _qr("c", fulllen_frac=0.0)
    assert _drops([qr], 0.0) == set()


def test_higher_threshold_drops_more():
    qr = _qr("c", fulllen_frac=0.3)
    assert _drops([qr], 0.2) == set()
    assert _drops([qr], 0.4) == {"c"}


# --------------------------------------------------------------------------
# aggregate_across_intervals threading
# --------------------------------------------------------------------------
def test_aggregate_preserves_max_fulllen_frac():
    a = _qr("c", assigned=("r1",), fulllen_frac=0.2)
    b = _qr("c", assigned=("r2",), fulllen_frac=0.6)
    agg = aggregate_across_intervals([[a], [b]])
    assert agg["c"].fulllen_frac == 0.6


def test_aggregate_unreachable_when_all_intervals_unreachable():
    a = _qr("c", assigned=("r1",), fulllen_frac=-1.0)
    b = _qr("c", assigned=("r2",), fulllen_frac=-1.0)
    agg = aggregate_across_intervals([[a], [b]])
    assert agg["c"].fulllen_frac == -1.0


def test_aggregate_reachable_beats_unreachable():
    a = _qr("c", assigned=("r1",), fulllen_frac=-1.0)
    b = _qr("c", assigned=("r2",), fulllen_frac=0.05)
    agg = aggregate_across_intervals([[a], [b]])
    assert agg["c"].fulllen_frac == 0.05
