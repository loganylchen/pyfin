"""Prototype tests for the retired production M3 coherence implementation.

The signal path (winner-anchored eventalign -> read x read DTW) is validated
end-to-end against real krill + blow5 data in the Docker ablation run; these
tests cover the window-signal projection and the no-window short-circuit, which
have no external deps.
"""
from __future__ import annotations

import numpy as np

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from experiments.m3_coherence.m3_junction_coherence import (
    _window_signal,
    build_m3_coherence,
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


def test_window_signal_orders_by_genome_and_normalizes():
    # exons (100,200),(300,400); tx pos 0->g100, 50->g150, 100->g300.
    c = _cand("a", [(200, 300)], 100, 400)
    res = {
        "position": [100, 0, 50],          # deliberately out of genomic order
        "event_level_mean": [10.0, 2.0, 6.0],
    }
    gset = {100, 150, 300}
    seg = _window_signal(res, c, gset)
    assert seg is not None
    # genomic order: g100(elm2), g150(elm6), g300(elm10) -> median 6, MAD 4.
    np.testing.assert_allclose(seg, np.array([-1.0, 0.0, 1.0], dtype=np.float32))


def test_window_signal_none_when_no_event_in_window():
    c = _cand("a", [(200, 300)], 100, 400)
    res = {"position": [0, 50], "event_level_mean": [2.0, 6.0]}
    assert _window_signal(res, c, gset=set()) is None


def test_window_signal_skips_nonfinite():
    c = _cand("a", [(200, 300)], 100, 400)
    res = {"position": [0, 50], "event_level_mean": [np.nan, 6.0]}
    seg = _window_signal(res, c, gset={100, 150})
    # only g150 survives -> single value, median==value -> 0.
    assert seg is not None and seg.shape == (1,)
    assert seg[0] == 0.0


def test_build_m3_no_window_returns_zeros():
    # Two structurally distinct candidates -> each its own class -> no window.
    a = _cand("a", [(200, 300)], 100, 400)
    b = _cand("b", [(150, 250), (300, 350)], 100, 400)
    read_ids = ["r0", "r1"]
    read_seqs = {"r0": "ACGT", "r1": "ACGT"}
    winner_col = np.array([0, 1], dtype=np.int64)
    m3 = build_m3_coherence(
        read_ids, read_seqs, [a, b], winner_col,
        sig_path="/nonexistent.blow5", pore="rna002",
        flank=2, chain_wobble=20, junction_k=10,
    )
    # No windowed class -> f5c never invoked; all-zero coherence.
    assert m3.shape == (2, 2)
    assert np.array_equal(m3, np.zeros((2, 2)))


def test_build_m3_empty():
    m3 = build_m3_coherence(
        [], {}, [], np.array([], dtype=np.int64),
        sig_path="/nonexistent.blow5",
    )
    assert m3.shape == (0, 0)
