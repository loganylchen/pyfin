"""Unit tests for the GUIDED junction-support gate (experiment, default-off).

``guided_junction_support_drops`` is the mirror of Lever 2 for GTF-passthrough
candidates: it drops a ``source == "gtf"`` multi-exon candidate if ANY of its
junctions is spliced by fewer than ``min_reads`` directly-observed reads
(strand-keyed ``observed_jct`` matched within ``tol`` bp). novel/fusion/mono are
exempt here; ``min_reads <= 0`` or empty ``observed_jct`` disables.

The decisive case is jitter: a GTF junction whose coordinates sit > ``tol`` bp
from the true (read-observed) site has ZERO exact support and is dropped — the
coordinate-exact check Lever 2 never applies to guided candidates.
"""
from __future__ import annotations

from collections import Counter

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.scoring.m2_junction_nll import guided_junction_support_drops


def _cand(introns, *, source="gtf", strand="+", chrom="chr1", cid="ENST1"):
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
    def test_min_reads_le_0_noop(self):
        c = _cand([(100, 200)])
        assert guided_junction_support_drops([c], _obs("+", {(100, 200): 0}),
                                             min_reads=0, tol=2) == set()

    def test_empty_or_none_observed_failopen(self):
        # None (fetch failure / not built) and {} (empty, indistinguishable from a
        # swallowed fetch failure) both fail-open — drop nothing. Recall-safe.
        c = _cand([(100, 200)])
        assert guided_junction_support_drops([c], None, min_reads=2, tol=2) == set()
        assert guided_junction_support_drops([c], {}, min_reads=2, tol=2) == set()


class TestGate:
    def test_under_supported_gtf_junction_dropped(self):
        c = _cand([(100, 200), (400, 500)])
        # first junction 5 reads, second only 1 -> drop (needs >=2)
        obs = _obs("+", {(100, 200): 5, (400, 500): 1})
        assert guided_junction_support_drops([c], obs, min_reads=2, tol=0) == {0}

    def test_all_junctions_supported_kept(self):
        c = _cand([(100, 200), (400, 500)])
        obs = _obs("+", {(100, 200): 3, (400, 500): 2})
        assert guided_junction_support_drops([c], obs, min_reads=2, tol=0) == set()

    def test_jitter_junction_beyond_tol_dropped(self):
        # THE case: GTF junction at (100,200) but reads splice (110,210) — a +10bp
        # jitter. With tol=2 the shifted GTF junction has 0 exact support -> dropped.
        c = _cand([(100, 200)])
        obs = _obs("+", {(110, 210): 8})
        assert guided_junction_support_drops([c], obs, min_reads=2, tol=2) == {0}
        # if reads DID splice the annotated site, it would be kept
        obs_true = _obs("+", {(100, 200): 8})
        assert guided_junction_support_drops([c], obs_true, min_reads=2, tol=2) == set()

    def test_tolerance_counts_nearby(self):
        c = _cand([(100, 200)])
        obs = _obs("+", {(102, 201): 2})  # 2bp off
        assert guided_junction_support_drops([c], obs, min_reads=2, tol=2) == set()
        assert guided_junction_support_drops([c], obs, min_reads=2, tol=0) == {0}


class TestStrand:
    def test_antisense_reads_do_not_support(self):
        c = _cand([(100, 200)], strand="+")
        obs = _obs("-", {(100, 200): 5})
        assert guided_junction_support_drops([c], obs, min_reads=2, tol=0) == {0}


class TestExemptions:
    def test_novel_exempt(self):
        # novel candidates go through Lever 2, not this gate
        c = _cand([(100, 200)], source="novel", cid="c")
        assert guided_junction_support_drops([c], _obs("+", {(100, 200): 0}),
                                             min_reads=2, tol=0) == set()

    def test_fusion_exempt(self):
        c = _cand([(100, 200)], source="fusion")
        assert guided_junction_support_drops([c], _obs("+", {(100, 200): 0}),
                                             min_reads=2, tol=0) == set()

    def test_mono_exempt(self):
        c = _cand([], source="gtf")  # single-exon: no junctions
        assert guided_junction_support_drops([c], {}, min_reads=2, tol=0) == set()


class TestMixed:
    def test_only_gtf_gated_novel_untouched(self):
        gtf_bad = _cand([(100, 200)], source="gtf", cid="ENST_bad")
        novel_bad = _cand([(100, 200)], source="novel", cid="novel_bad")
        obs = _obs("+", {(300, 400): 9})  # neither is supported at its own site
        # only the GTF candidate (index 0) is dropped; novel (index 1) is exempt here
        assert guided_junction_support_drops([gtf_bad, novel_bad], obs,
                                             min_reads=2, tol=0) == {0}
