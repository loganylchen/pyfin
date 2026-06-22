"""Unit tests for the M2/M1 read-support gate.

Covers ``fin.scoring.m2_junction_nll.support_gate_drops`` — a multi-exon
candidate is kept iff it is some read's M1 sole best-AS OR some read's M2-best
(lowest junction-NLL); otherwise dropped. Catches corrupted-annotation junctions
that win neither, without误删 genuine isoforms.
"""
from __future__ import annotations

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.scoring.m2_junction_nll import support_gate_drops


def _cand(cid, source="novel", introns=((100, 200),)):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=300,
        sequence="ACGT" * 50,
        source=source,
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=0,
        end=300,
    )


def test_unsupported_candidate_dropped():
    # c0 is no read's sole-AS and no read's M2-best -> dropped.
    cands = [_cand("c0"), _cand("c1")]
    ties = {0: [0, 1]}                 # both tied at M1 (neither is sole)
    nlls = {0: {0: 5.0, 1: 1.0}}       # c1 wins M2; c0 loses
    assert support_gate_drops(cands, ties, nlls) == {0}


def test_m1_sole_best_keeps():
    # c0 is read 0's SOLE M1 best-AS -> kept even with no M2 score.
    cands = [_cand("c0"), _cand("c1")]
    ties = {0: [0], 1: [1]}            # each read uniquely picks its candidate
    nlls = {}
    assert support_gate_drops(cands, ties, nlls) == set()


def test_m2_best_keeps():
    # c0 wins M2 for read 0 -> kept; c1 wins for read 1 -> kept.
    cands = [_cand("c0"), _cand("c1")]
    ties = {0: [0, 1], 1: [0, 1]}      # tied at M1
    nlls = {0: {0: 1.0, 1: 5.0}, 1: {0: 5.0, 1: 1.0}}
    assert support_gate_drops(cands, ties, nlls) == set()


def test_tie_accept_keeps_both():
    # read 0 ties c0/c1 at the lowest NLL; tie_ok=True keeps both.
    cands = [_cand("c0"), _cand("c1")]
    ties = {0: [0, 1]}
    nlls = {0: {0: 1.0, 1: 1.0}}
    assert support_gate_drops(cands, ties, nlls, tie_ok=True) == set()


def test_tie_strict_drops_both():
    # same tie, tie_ok=False: no strict unique winner -> neither supported -> both drop.
    cands = [_cand("c0"), _cand("c1")]
    ties = {0: [0, 1]}
    nlls = {0: {0: 1.0, 1: 1.0}}
    assert support_gate_drops(cands, ties, nlls, tie_ok=False) == {0, 1}


def test_fusion_exempt():
    cands = [_cand("f", source="fusion")]
    assert support_gate_drops(cands, {}, {}) == set()


def test_mono_exon_exempt():
    cands = [_cand("m", introns=())]   # single-exon: no internal junction
    assert support_gate_drops(cands, {}, {}) == set()


def test_gtf_not_exempt():
    # A GTF candidate winning no read's support IS dropped (the corrupted-GTF case).
    cands = [_cand("gtf_bad", source="gtf"), _cand("c1")]
    ties = {0: [0, 1]}
    nlls = {0: {0: 9.0, 1: 1.0}}       # gtf_bad loses M2, never sole-AS
    assert support_gate_drops(cands, ties, nlls) == {0}


def test_empty_inputs():
    assert support_gate_drops([], {}, {}) == set()


def test_no_evidence_drops_nothing():
    # No ties AND no nlls (e.g. no signal backend + diff-cover gate off): the gate
    # has nothing to judge on, so it must NOT drop any candidate.
    cands = [_cand("c0"), _cand("c1", source="gtf")]
    assert support_gate_drops(cands, {}, {}) == set()
