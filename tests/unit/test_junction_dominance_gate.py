"""Unit tests for the pre-EM junction-dominance gate (junction-first idea).

``dominant_junction_drops`` drops a NOVEL multi-exon candidate if any junction is
(a) supported by < min_reads observed reads, OR (b) not locally dominant — a
DIFFERENT observed junction within `window` bp (beyond `tol`) has strictly more
reads. gtf/fusion/mono exempt; min_reads<=0 or empty observed disables.
"""
from __future__ import annotations

from collections import Counter

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.scoring.m2_junction_nll import dominant_junction_drops


def _cand(introns, *, source="novel", strand="+", cid="c"):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=(introns[-1][1] + 100) if introns else 200,
        sequence="ACGT" * 50,
        source=source,
        supporting_read_ids=set(),
        chrom="chr1",
        strand=strand,
        start=1,
        end=(introns[-1][1] + 100) if introns else 200,
    )


def _obs(strand, pairs):
    return {strand: Counter(pairs)}


def _drop(cands, obs, **kw):
    kw.setdefault("tol", 2)
    return dominant_junction_drops(cands, obs, **kw)


class TestDisabled:
    def test_min_reads_zero(self):
        c = _cand([(100, 200)])
        assert _drop([c], _obs("+", {(100, 200): 0}), min_reads=0, window=20) == set()

    def test_empty_observed(self):
        c = _cand([(100, 200)])
        assert _drop([c], {}, min_reads=2, window=20) == set()
        assert _drop([c], None, min_reads=2, window=20) == set()


class TestDominance:
    def test_dominant_junction_kept(self):
        # the candidate's junction is the strongest in its neighborhood
        c = _cand([(100, 200)])
        obs = _obs("+", {(100, 200): 5, (108, 205): 1})  # neighbor weaker
        assert _drop([c], obs, min_reads=2, window=20) == set()

    def test_dominated_shadow_dropped(self):
        # the candidate IS the weak shadow; a different junction 8bp away is stronger
        c = _cand([(108, 205)])
        obs = _obs("+", {(100, 200): 9, (108, 205): 2})  # true junction dominates
        assert _drop([c], obs, min_reads=2, window=20) == {0}

    def test_tie_kept(self):
        # equal-support neighbor does NOT demote (needs STRICTLY more)
        c = _cand([(108, 205)])
        obs = _obs("+", {(100, 200): 3, (108, 205): 3})
        assert _drop([c], obs, min_reads=2, window=20) == set()

    def test_stronger_junction_beyond_window_does_not_demote(self):
        c = _cand([(108, 205)])
        # the stronger junction is 100bp away -> outside window -> no demotion
        obs = _obs("+", {(300, 405): 9, (108, 205): 3})
        assert _drop([c], obs, min_reads=2, window=20) == set()


class TestCount:
    def test_weak_junction_dropped(self):
        c = _cand([(100, 200)])
        obs = _obs("+", {(100, 200): 1})  # below min_reads=2, no competitor
        assert _drop([c], obs, min_reads=2, window=20) == {0}


class TestExemptions:
    def test_gtf_exempt(self):
        c = _cand([(108, 205)], source="gtf", cid="ENST1")
        obs = _obs("+", {(100, 200): 9, (108, 205): 1})
        assert _drop([c], obs, min_reads=2, window=20) == set()

    def test_fusion_exempt(self):
        c = _cand([(108, 205)], source="fusion")
        obs = _obs("+", {(100, 200): 9, (108, 205): 1})
        assert _drop([c], obs, min_reads=2, window=20) == set()

    def test_mono_exempt(self):
        c = _cand([], source="novel")
        assert _drop([c], {}, min_reads=2, window=20) == set()
