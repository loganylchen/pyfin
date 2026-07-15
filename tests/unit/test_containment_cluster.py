"""Unit tests for containment_cluster_drops — the recall-safER, read-guarded,
wobble-tolerant sub-chain shadow drop.

A NOVEL candidate whose intron chain is a contiguous SUB-CHAIN (within wobble_bp) of
a longer candidate is dropped ONLY when it is a low-support shadow by BOTH EM
abundance (<= parent*min_ab_ratio) AND supporting-read count (<= parent*min_read_ratio).
The read-support guard is what keeps a genuine low-abundance short isoform (reads
comparable to the parent) from being folded away.
"""
from __future__ import annotations

import uuid

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.scoring.m2_junction_nll import containment_cluster_drops

J1, J2, J3 = (100, 200), (400, 500), (700, 800)


def _cand(introns, *, start, end, source="novel", strand="+", chrom="chr1"):
    return TranscriptCandidate(
        candidate_id="novel_" + uuid.uuid4().hex[:8],
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=end, sequence="ACGT" * 50, source=source,
        supporting_read_ids=set(), chrom=chrom, strand=strand, start=start, end=end)


def _drop(cands, em_ab, reads, **kw):
    kw.setdefault("wobble_bp", 6)
    kw.setdefault("min_ab_ratio", 0.3)
    kw.setdefault("min_read_ratio", 0.3)
    return containment_cluster_drops(cands, em_ab, reads, **kw)


class TestContainmentClusterDrop:
    def test_low_support_subchain_dropped(self):
        parent = _cand([J1, J2, J3], start=50, end=850)
        shadow = _cand([J2, J3], start=350, end=850)   # 3'-suffix sub-chain
        assert _drop([parent, shadow], [10.0, 2.0], [13, 1]) == {1}

    def test_real_short_isoform_kept_by_read_guard(self):
        # SAME low abundance, but comparable read support -> a genuine short isoform.
        parent = _cand([J1, J2, J3], start=50, end=850)
        shadow = _cand([J2, J3], start=350, end=850)
        assert _drop([parent, shadow], [10.0, 2.0], [13, 10]) == set()

    def test_high_abundance_subchain_kept(self):
        parent = _cand([J1, J2, J3], start=50, end=850)
        shadow = _cand([J2, J3], start=350, end=850)
        assert _drop([parent, shadow], [10.0, 5.0], [13, 1]) == set()

    def test_wobble_tolerant_match(self):
        # shadow junctions within 6bp of the parent's suffix junctions.
        parent = _cand([J1, J2, J3], start=50, end=850)
        shadow = _cand([(402, 498), (700, 800)], start=350, end=850)
        assert _drop([parent, shadow], [10.0, 2.0], [13, 1]) == {1}

    def test_beyond_wobble_not_matched(self):
        parent = _cand([J1, J2, J3], start=50, end=850)
        shadow = _cand([(430, 470), (700, 800)], start=350, end=850)  # 30bp off
        assert _drop([parent, shadow], [10.0, 2.0], [13, 1]) == set()

    def test_internal_subchain_dropped(self):
        parent = _cand([J1, J2, J3], start=50, end=850)
        shadow = _cand([J2], start=350, end=550)   # middle single-intron sub-chain
        assert _drop([parent, shadow], [10.0, 2.0], [13, 1]) == {1}

    def test_gtf_shadow_never_dropped(self):
        parent = _cand([J1, J2, J3], start=50, end=850)
        shadow = _cand([J2, J3], start=350, end=850, source="gtf")
        assert _drop([parent, shadow], [10.0, 2.0], [13, 1]) == set()

    def test_equal_length_not_subchain(self):
        # same intron count -> handled by the wobble cluster, not this function.
        a = _cand([J1, J2], start=50, end=550)
        b = _cand([J1, J2], start=50, end=550)
        assert _drop([a, b], [10.0, 2.0], [13, 1]) == set()

    def test_chained_containment_both_dropped(self):
        big = _cand([J1, J2, J3], start=50, end=850)
        mid = _cand([J2, J3], start=350, end=850)
        small = _cand([J3], start=650, end=850)
        # both mid and small are low-support sub-chains -> both dropped, big kept.
        assert _drop([big, mid, small], [10.0, 2.0, 1.0], [20, 2, 1]) == {1, 2}


class TestGuardsAndExemptions:
    def test_gtf_parent_absorbs_novel_shadow(self):
        # a gtf parent is a valid parent; the novel sub-chain shadow still drops.
        parent = _cand([J1, J2, J3], start=50, end=850, source="gtf")
        shadow = _cand([J2, J3], start=350, end=850)   # novel
        assert _drop([parent, shadow], [10.0, 2.0], [13, 1]) == {1}

    def test_fusion_never_parent_nor_shadow(self):
        # fusion excluded from the bucket entirely: a novel sub-chain of a FUSION is
        # NOT dropped (fusion is not a valid parent), and a fusion sub-chain is kept.
        fusion = _cand([J1, J2, J3], start=50, end=850, source="fusion")
        novel_sub = _cand([J2, J3], start=350, end=850)
        assert _drop([fusion, novel_sub], [10.0, 2.0], [13, 1]) == set()
        parent = _cand([J1, J2, J3], start=50, end=850)
        fusion_sub = _cand([J2, J3], start=350, end=850, source="fusion")
        assert _drop([parent, fusion_sub], [10.0, 2.0], [13, 1]) == set()

    def test_exclude_skips_shadow_and_parent(self):
        parent = _cand([J1, J2, J3], start=50, end=850)
        shadow = _cand([J2, J3], start=350, end=850)
        # excluding the parent (col 0) -> no valid parent -> shadow kept.
        assert _drop([parent, shadow], [10.0, 2.0], [13, 1], exclude={0}) == set()
        # excluding the shadow (col 1) -> shadow not considered -> nothing dropped.
        assert _drop([parent, shadow], [10.0, 2.0], [13, 1], exclude={1}) == set()

    def test_ratio_equality_is_inclusive(self):
        # em_ab[a] == parent*ratio and reads[a] == parent*ratio -> dropped (<=).
        parent = _cand([J1, J2, J3], start=50, end=850)
        shadow = _cand([J2, J3], start=350, end=850)
        assert _drop([parent, shadow], [10.0, 3.0], [10, 3]) == {1}   # 3<=10*.3, 3<=10*.3

    def test_deterministic_under_shuffled_input(self):
        # same candidate set in two orders -> the dropped candidate_ids match.
        big = _cand([J1, J2, J3], start=50, end=850)
        mid = _cand([J2, J3], start=350, end=850)
        cands = [big, mid]
        d1 = _drop(cands, [10.0, 2.0], [13, 1])
        cands_rev = [mid, big]
        d2 = _drop(cands_rev, [2.0, 10.0], [1, 13])
        # d1 drops index 1 (mid); d2 drops index 0 (mid) -> same candidate.
        assert {cands[i].candidate_id for i in d1} == \
               {cands_rev[i].candidate_id for i in d2} == {mid.candidate_id}

    def test_absolute_cap_protects_high_read_shadow(self):
        # parent 100 reads, shadow 20 reads & 20% abundance: passes BOTH ratios
        # (20 <= 100*0.3=30, 4 <= 20*0.3) but exceeds the absolute cap -> kept.
        parent = _cand([J1, J2, J3], start=50, end=850)
        shadow = _cand([J2, J3], start=350, end=850)
        assert _drop([parent, shadow], [20.0, 4.0], [100, 20],
                     max_shadow_reads=10) == set()
        # same, cap disabled (0) -> ratio-only -> dropped.
        assert _drop([parent, shadow], [20.0, 4.0], [100, 20],
                     max_shadow_reads=0) == {1}
        # a small shadow under the cap is still dropped.
        assert _drop([parent, shadow], [20.0, 2.0], [100, 3],
                     max_shadow_reads=10) == {1}

    def test_zero_reads_parent_not_divzero(self):
        parent = _cand([J1, J2, J3], start=50, end=850)
        shadow = _cand([J2, J3], start=350, end=850)
        # parent 0 reads -> read guard 0<=0*.3=0 only if shadow also 0; no crash.
        assert _drop([parent, shadow], [10.0, 0.0], [0, 0]) == {1}
        assert _drop([parent, shadow], [10.0, 2.0], [0, 1]) == set()  # 1 > 0*.3
