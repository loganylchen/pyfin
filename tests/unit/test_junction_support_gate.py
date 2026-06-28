"""Unit tests for Lever 2 — per-junction read-support gate.

``junction_support_drops`` drops a NOVEL multi-exon candidate if ANY of its
junctions is spliced by fewer than ``min_reads`` directly-observed reads
(strand-keyed ``observed_jct`` matched within ``tol`` bp). gtf/fusion/mono
exempt; min_reads<=1 or empty observed_jct disables.
"""
from __future__ import annotations

from collections import Counter

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.scoring.m2_junction_nll import junction_support_drops


def _cand(introns, *, source="novel", strand="+", chrom="chr1", cid="c"):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=(introns[-1][1] + 100) if introns else 200,
        sequence="ACGT" * 50,
        source=source,
        supporting_read_ids=set(),
        chrom=chrom,
        strand=strand,
        start=1,
        end=(introns[-1][1] + 100) if introns else 200,
    )


def _obs(d, pairs):  # {strand: {(donor,acceptor): n}}
    return {d: Counter(pairs)}


class TestDisabled:
    def test_min_reads_le_1_noop(self):
        c = _cand([(100, 200)])
        assert junction_support_drops([c], _obs("+", {(100, 200): 0}),
                                      min_reads=1, tol=0) == set()

    def test_empty_observed_noop(self):
        c = _cand([(100, 200)])
        assert junction_support_drops([c], None, min_reads=2, tol=0) == set()
        assert junction_support_drops([c], {}, min_reads=2, tol=0) == set()


class TestGate:
    def test_under_supported_junction_dropped(self):
        c = _cand([(100, 200), (400, 500)])
        # first junction has 5 reads, second only 1 -> drop (needs >=2)
        obs = _obs("+", {(100, 200): 5, (400, 500): 1})
        assert junction_support_drops([c], obs, min_reads=2, tol=0) == {0}

    def test_all_junctions_supported_kept(self):
        c = _cand([(100, 200), (400, 500)])
        obs = _obs("+", {(100, 200): 3, (400, 500): 2})
        assert junction_support_drops([c], obs, min_reads=2, tol=0) == set()

    def test_tolerance_counts_nearby(self):
        c = _cand([(100, 200)])
        obs = _obs("+", {(102, 201): 2})  # 2bp off the candidate junction
        # tol=2 -> the 2 nearby reads count -> kept (no drop)
        assert junction_support_drops([c], obs, min_reads=2, tol=2) == set()
        # tol=0 -> strict, nothing matches at exactly (100,200) -> 0 reads -> drop
        assert junction_support_drops([c], obs, min_reads=2, tol=0) == {0}


class TestStrand:
    def test_antisense_reads_do_not_support(self):
        c = _cand([(100, 200)], strand="+")
        # support only on the '-' strand -> '+' candidate sees 0 -> drop
        obs = _obs("-", {(100, 200): 5})
        assert junction_support_drops([c], obs, min_reads=2, tol=0) == {0}


class TestExemptions:
    def test_gtf_exempt(self):
        c = _cand([(100, 200)], source="gtf", cid="ENST1")
        assert junction_support_drops([c], _obs("+", {(100, 200): 0}),
                                      min_reads=2, tol=0) == set()

    def test_fusion_exempt(self):
        c = _cand([(100, 200)], source="fusion")
        assert junction_support_drops([c], _obs("+", {(100, 200): 0}),
                                      min_reads=2, tol=0) == set()

    def test_mono_exempt(self):
        c = _cand([], source="novel")  # single-exon: no junctions
        assert junction_support_drops([c], {}, min_reads=2, tol=0) == set()
