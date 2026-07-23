"""Unit tests for fin/scoring/m2_junction_nll.py pure-logic helpers.

The signal path (``read_cand_mean_nll`` GTF<NOV behaviour) is validated
end-to-end against real krill + blow5 data in the Docker ablation run; these
tests cover the coordinate / diff-window logic that has no external deps.
"""
from __future__ import annotations

import numpy as np

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.scoring.m2_junction_nll import (
    _build_exons,
    _internal_bounds,
    _mean_nll_in_gset,
    _mean_nll_in_window,
    _tx2genome,
    diff_junction_windows,
    internal_diff_regions,
    tx2genome_array,
)


def _cand(cid, introns, start, end, strand="+", seq="A"):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=end if strand == "+" else start,
        sequence=seq,
        source="gtf",
        supporting_read_ids=set(),
        chrom="chrT",
        strand=strand,
        start=start,
        end=end,
    )


def test_build_exons_plus_strand():
    c = _cand("a", [(200, 300), (400, 500)], 100, 600)
    assert _build_exons(c) == [(100, 200), (300, 400), (500, 600)]


def test_internal_bounds_drops_terminal_boundaries():
    c = _cand("a", [(200, 300), (400, 500)], 100, 600)
    # flat = [200,300,400,500]; internal = [300,400]
    assert _internal_bounds(c) == [300, 400]


def test_tx2genome_plus_strand():
    c = _cand("a", [(200, 300)], 100, 400)  # exons (100,200),(300,400)
    assert _tx2genome(c, 0) == 100
    assert _tx2genome(c, 50) == 150
    assert _tx2genome(c, 100) == 300  # second exon, offset 0


def test_tx2genome_array_matches_scalar_plus():
    # Vectorized projection must equal the scalar _tx2genome for every position,
    # with -1 where the scalar returns None (out of range).
    c = _cand("a", [(200, 300), (400, 500)], 100, 600)  # tx len = 300
    positions = list(range(-2, 305))
    vec = tx2genome_array(c, np.array(positions, dtype=np.int64))
    for p, gv in zip(positions, vec.tolist()):
        scalar = _tx2genome(c, p)
        assert gv == (scalar if scalar is not None else -1)


def test_tx2genome_array_matches_scalar_minus():
    c = _cand("m", [(200, 300)], 100, 400, strand="-")  # tx len = 200
    positions = list(range(-1, 205))
    vec = tx2genome_array(c, np.array(positions, dtype=np.int64))
    for p, gv in zip(positions, vec.tolist()):
        scalar = _tx2genome(c, p)
        assert gv == (scalar if scalar is not None else -1)


def test_tx2genome_array_empty():
    c = _cand("a", [(200, 300)], 100, 400)
    assert tx2genome_array(c, np.array([], dtype=np.int64)).tolist() == []


def test_tx2genome_minus_strand():
    c = _cand("a", [(200, 300)], 100, 400, strand="-")
    # transcript 5'->3' order reversed: (300,400) then (100,200)
    assert _tx2genome(c, 0) == 399
    assert _tx2genome(c, 100) == 199


def test_internal_diff_region_for_internal_wobble():
    # Same class: internal donor of the 2nd intron shifted 400 -> 410.
    a = _cand("a", [(200, 300), (400, 500)], 100, 600)
    b = _cand("b", [(200, 300), (410, 500)], 100, 600)
    regions = internal_diff_regions([a, b], flank=2)
    # sym diff sliver is (400,410) inside internal span [300,410]; padded +/-2.
    assert regions == [(398, 412)]


def test_terminal_differences_excluded():
    # Only 5' start and 3' end differ; intron chain identical -> no internal diff.
    a = _cand("a", [(200, 300), (400, 500)], 100, 600)
    d = _cand("d", [(200, 300), (400, 500)], 120, 650)
    assert internal_diff_regions([a, d], flank=2) == []


def test_internal_diff_region_singleton_empty():
    a = _cand("a", [(200, 300), (400, 500)], 100, 600)
    assert internal_diff_regions([a], flank=2) == []


# --- diff_junction_windows: tight +-flank windows around DIFFERING boundaries ----
def test_diff_junction_windows_acceptor_wobble():
    # acceptor shifts 300 -> 304; donor 200 shared. As whole-intron tuples both the
    # shared donor (200) and the two acceptors (300,304) are collected; +-6 merges
    # the two nearby acceptor windows.
    a = _cand("a", [(200, 300)], 100, 400)
    b = _cand("b", [(200, 304)], 100, 400)
    assert diff_junction_windows([a, b], flank=6) == [(194, 206), (294, 310)]


def test_diff_junction_windows_unanimous_empty():
    a = _cand("a", [(200, 300)], 100, 400)
    b = _cand("b", [(200, 300)], 100, 500)  # identical chain, only 3' end differs
    assert diff_junction_windows([a, b], flank=6) == []


def test_diff_junction_windows_singleton_empty():
    a = _cand("a", [(200, 300)], 100, 400)
    assert diff_junction_windows([a], flank=6) == []


# --- reduce="sum" is exactly n * mean (summed-LLR building block) ----------------
def _res(positions, event_means, model_mean=100.0, model_stdv=1.0, kmer="AAAAA"):
    n = len(positions)
    return {
        "position": list(positions),
        "reference_kmer": [kmer] * n,
        "model_kmer": [kmer] * n,
        "event_level_mean": list(event_means),
        "model_mean": [model_mean] * n,
        "model_stdv": [model_stdv] * n,
        "status": 0,
    }


def test_reduce_sum_is_n_times_mean_gset():
    # cand exons (100,200),(300,400); tx 10->g110, 20->g120, 150->g350.
    c = _cand("a", [(200, 300)], 100, 400, seq="A" * 200)
    res = _res([10, 20, 150], [101.0, 102.0, 103.0])  # z=1,2,3 -> NLL 0.5,2.0,4.5
    gset = {110, 120}  # excludes g350
    mean, n_m = _mean_nll_in_gset(res, c, gset, reduce="mean")
    total, n_s = _mean_nll_in_gset(res, c, gset, reduce="sum")
    assert n_m == n_s == 2
    assert total == mean * 2
    assert abs(total - (0.5 + 2.0)) < 1e-9  # only the two in-window events


def test_reduce_sum_is_n_times_mean_window():
    c = _cand("a", [(200, 300)], 100, 400, seq="A" * 200)
    res = _res([10, 20, 150], [101.0, 102.0, 103.0])
    windows = [(105, 125)]  # genomic; captures g110,g120, not g350
    mean, n_m = _mean_nll_in_window(res, c, windows, reduce="mean")
    total, n_s = _mean_nll_in_window(res, c, windows, reduce="sum")
    assert n_m == n_s == 2
    assert total == mean * 2


def test_reduce_default_is_mean():
    c = _cand("a", [(200, 300)], 100, 400, seq="A" * 200)
    res = _res([10, 20], [101.0, 103.0])  # z=1,3 -> NLL 0.5,4.5 -> mean 2.5
    mean_default, _ = _mean_nll_in_gset(res, c, {110, 120})
    assert abs(mean_default - 2.5) < 1e-9
