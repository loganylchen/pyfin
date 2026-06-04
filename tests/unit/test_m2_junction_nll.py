"""Unit tests for fin/scoring/m2_junction_nll.py pure-logic helpers.

The signal path (``read_cand_mean_nll`` GTF<NOV behaviour) is validated
end-to-end against real krill + blow5 data in the Docker ablation run; these
tests cover the coordinate / diff-window logic that has no external deps.
"""
from __future__ import annotations

import numpy as np

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.scoring.m2_junction_nll import (
    MISSING,
    _build_exons,
    _internal_bounds,
    _tx2genome,
    build_m2_distance,
    internal_diff_regions,
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


def test_build_m2_distance_no_window_degrades_to_m1():
    # Two structurally distinct candidates (different intron counts) -> each is
    # its own class -> no discrimination window -> M2 == M1 on accepted cells,
    # MISSING on rejected cells.
    a = _cand("a", [(200, 300)], 100, 400)
    b = _cand("b", [(150, 250), (300, 350)], 100, 400)
    read_ids = ["r0", "r1"]
    read_seqs = {"r0": "ACGT", "r1": "ACGT"}
    mask = np.array([[True, False], [False, True]])
    dist_m1 = np.array([[0.0, MISSING], [MISSING, 3.0]])
    out = build_m2_distance(
        read_ids, read_seqs, [a, b], mask, dist_m1,
        sig_path="/nonexistent.blow5", pore="rna002", flank=2, chain_wobble=20,
    )
    # No windowed class -> f5c never invoked; falls back to M1 exactly.
    assert out[0, 0] == 0.0
    assert out[0, 1] == MISSING
    assert out[1, 0] == MISSING
    assert out[1, 1] == 3.0
