"""Unit tests for the cassette-skip extension of ``structural_wobble_clusters``.

The cassette extension unions a K-intron candidate with a (K-1)-intron candidate
when the chains align as: one small extra exon (< cassette_max_exon_bp) in the
K-intron chain, replacing one spanning intron in the (K-1)-intron chain, with
every other junction within ``bp``. This targets minimap2's small-exon-skip
failure mode where the same biological isoform splits across two candidates.
"""
from __future__ import annotations

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.scoring.m2_junction_nll import (
    structural_wobble_clusters,
    wobble_shadow_drops,
)


def _cand(introns, *, source="novel", strand="+", chrom="chr1"):
    seq = "ACGT" * 200
    ex = []
    if introns:
        ex.append((1, introns[0][0]))
        for i, (s, e) in enumerate(introns):
            nxt = introns[i + 1][0] if i + 1 < len(introns) else e + 200
            ex.append((e, nxt))
    return TranscriptCandidate(
        candidate_id="c",
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=ex[-1][1] if ex else 200,
        sequence=seq,
        source=source,
        supporting_read_ids=set(),
        chrom=chrom,
        strand=strand,
        start=1,
        end=ex[-1][1] if ex else 200,
    )


def _classes(clusters):
    """Return classes of size >= 2 as a list of sorted-tuple cluster member lists."""
    return sorted(
        tuple(sorted(cols)) for cols in clusters.values() if len(cols) >= 2
    )


class TestCassetteDisabled:
    def test_cassette_disabled_by_default(self):
        # K=2 (with 64bp cassette exon at 4003-4067) vs K=1 (spanning intron).
        a = _cand([(2005, 4003), (4067, 6057)])
        b = _cand([(2005, 6057)])
        clusters = structural_wobble_clusters([a, b], bp=20)
        assert _classes(clusters) == []  # cassette_max_exon_bp default 0 -> no join

    def test_cassette_disabled_explicitly(self):
        a = _cand([(2005, 4003), (4067, 6057)])
        b = _cand([(2005, 6057)])
        clusters = structural_wobble_clusters(
            [a, b], bp=20, cassette_max_exon_bp=0
        )
        assert _classes(clusters) == []


class TestCassetteEnabled:
    def test_cassette_exact_at_position_0(self):
        # SIRV-style: A has 2 introns flanking a 64bp cassette exon at 4003-4067.
        # B has 1 spanning intron 2005-6057. Should cluster.
        a = _cand([(2005, 4003), (4067, 6057)])
        b = _cand([(2005, 6057)])
        clusters = structural_wobble_clusters(
            [a, b], bp=20, cassette_max_exon_bp=70
        )
        assert _classes(clusters) == [(0, 1)]

    def test_cassette_with_outer_wobble_within_bp(self):
        # outer ends differ by 5bp / 3bp (within bp=20). Should still cluster.
        a = _cand([(2005, 4003), (4067, 6057)])
        b = _cand([(2010, 6060)])
        clusters = structural_wobble_clusters(
            [a, b], bp=20, cassette_max_exon_bp=70
        )
        assert _classes(clusters) == [(0, 1)]

    def test_cassette_too_long_rejected(self):
        # cassette exon is 200bp (4003-4203), exceeds threshold 70.
        a = _cand([(2005, 4003), (4203, 6057)])
        b = _cand([(2005, 6057)])
        clusters = structural_wobble_clusters(
            [a, b], bp=20, cassette_max_exon_bp=70
        )
        assert _classes(clusters) == []

    def test_cassette_outer_wobble_too_big(self):
        # outer ends differ by 25bp (exceeds bp=20).
        a = _cand([(2005, 4003), (4067, 6057)])
        b = _cand([(2030, 6057)])
        clusters = structural_wobble_clusters(
            [a, b], bp=20, cassette_max_exon_bp=70
        )
        assert _classes(clusters) == []

    def test_cassette_middle_of_chain(self):
        # 4-intron A vs 3-intron B; cassette at middle position (i=1).
        # A's middle cassette exon spans 400..460 (60 bp), which < 70 → join.
        # A: [(100,200), (260,400), (460,600), (650,750)]
        # B: [(100,200),            (260,600), (650,750)]
        a = _cand([(100, 200), (260, 400), (460, 600), (650, 750)])
        b = _cand([(100, 200), (260, 600), (650, 750)])
        clusters = structural_wobble_clusters(
            [a, b], bp=20, cassette_max_exon_bp=70
        )
        assert _classes(clusters) == [(0, 1)]

    def test_cassette_at_boundary_exclusive(self):
        # cassette_size = 70 is rejected (the threshold is strict <).
        a = _cand([(100, 200), (260, 400), (470, 600), (650, 750)])
        b = _cand([(100, 200), (260, 600), (650, 750)])
        clusters = structural_wobble_clusters(
            [a, b], bp=20, cassette_max_exon_bp=70
        )
        assert _classes(clusters) == []

    def test_non_cassette_difference_not_joined(self):
        # A and B differ by intron count 1 but the chains don't align as cassette
        # (the differing region is structural, not a small extra exon).
        a = _cand([(100, 200), (300, 500)])      # exons end at 100, then 200..300, 500..end
        b = _cand([(100, 200)])                  # only 1 intron
        clusters = structural_wobble_clusters(
            [a, b], bp=20, cassette_max_exon_bp=70
        )
        # A has intron at the END that B lacks (not a cassette swap inside the chain)
        # -> should NOT cluster.
        assert _classes(clusters) == []

    def test_cassette_different_chrom_strand_not_joined(self):
        a = _cand([(2005, 4003), (4067, 6057)], chrom="chr1", strand="+")
        b = _cand([(2005, 6057)], chrom="chr2", strand="+")
        c = _cand([(2005, 6057)], chrom="chr1", strand="-")
        clusters = structural_wobble_clusters(
            [a, b, c], bp=20, cassette_max_exon_bp=70
        )
        assert _classes(clusters) == []

    def test_cassette_three_way_cluster_via_transitive_union(self):
        # Anchor B (1 intron, the long-span). A1 and A2 are two cassette variants
        # of B (slightly different wobble), both <70bp cassette. All three should
        # end up in one cluster via union-find transitivity.
        b = _cand([(2005, 6057)])
        a1 = _cand([(2005, 4003), (4067, 6057)])   # cassette 64bp
        a2 = _cand([(2010, 4007), (4070, 6060)])   # similar cassette 63bp, slight wobble
        clusters = structural_wobble_clusters(
            [b, a1, a2], bp=20, cassette_max_exon_bp=70
        )
        assert _classes(clusters) == [(0, 1, 2)]


class TestSameCountWobbleUnaffected:
    """The classic same-intron-count wobble path must not regress when
    cassette is enabled."""

    def test_pure_wobble_still_clusters(self):
        a = _cand([(200, 300)])
        b = _cand([(208, 305)])
        clusters = structural_wobble_clusters(
            [a, b], bp=20, cassette_max_exon_bp=70
        )
        assert _classes(clusters) == [(0, 1)]

    def test_pure_wobble_beyond_bp_no_cluster(self):
        a = _cand([(200, 300)])
        b = _cand([(225, 300)])  # 25 bp donor wobble, exceeds bp=20
        clusters = structural_wobble_clusters(
            [a, b], bp=20, cassette_max_exon_bp=70
        )
        assert _classes(clusters) == []


def _src_cand(source, introns=((200, 300), (400, 500))):
    return TranscriptCandidate(
        candidate_id=f"{source}_c",
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=600,
        sequence="ACGT" * 200,
        source=source,
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=100,
        end=600,
    )


class TestWobbleShadowDrops:
    """Shadow-drop logic: novel by abundance; GTF by direct read-support guard."""

    def _one_cluster(self, n):
        return {0: list(range(n))}

    # anchor introns ((200,300),(400,500)); GTF sibling differs at junction 2 donor
    # (414 vs 400 = 14bp) -> distinguishing junction is (414,500).
    def _anchor(self, source="novel"):
        return _src_cand(source, introns=((200, 300), (400, 500)))

    def _gtf_sibling(self):
        return _src_cand("gtf", introns=((200, 300), (414, 500)))

    def test_gtf_dropped_when_distinguishing_jct_has_no_reads(self):
        cands = [self._anchor("novel"), self._gtf_sibling()]
        em_ab = [30.0, 2.0]  # gtf/anchor = 0.067 < 0.15
        observed = {"+": {(200, 300): 50, (400, 500): 40}}  # NO read at (414,500)
        drops = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                    gtf_drop_enabled=True, observed_jct=observed,
                                    jct_tol=0, gtf_min_jct_reads=1)
        assert drops == {1}  # phantom jittered GTF dropped

    def test_gtf_kept_when_distinguishing_jct_has_reads(self):
        # SIRV606-style: the GTF's own junction (414,500) is directly observed.
        cands = [self._anchor("novel"), self._gtf_sibling()]
        em_ab = [30.0, 2.0]
        observed = {"+": {(200, 300): 50, (400, 500): 40, (414, 500): 9}}
        drops = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                    gtf_drop_enabled=True, observed_jct=observed,
                                    jct_tol=0, gtf_min_jct_reads=1)
        assert drops == set()  # real isoform with its own reads survives

    def test_gtf_kept_when_anchor_is_gtf_but_now_judged_by_own_support(self):
        # anchor source no longer matters: a GTF sibling with read support is kept...
        cands = [self._anchor("gtf"), self._gtf_sibling()]
        em_ab = [30.0, 2.0]
        observed = {"+": {(414, 500): 5}}
        drops = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                    gtf_drop_enabled=True, observed_jct=observed,
                                    jct_tol=0, gtf_min_jct_reads=1)
        assert drops == set()

    def test_gtf_dropped_even_under_gtf_anchor_if_no_support(self):
        # ...and a zero-support GTF sibling IS dropped even when the anchor is GTF.
        cands = [self._anchor("gtf"), self._gtf_sibling()]
        em_ab = [30.0, 2.0]
        observed = {"+": {(200, 300): 9}}  # nothing at (414,500)
        drops = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                    gtf_drop_enabled=True, observed_jct=observed,
                                    jct_tol=0, gtf_min_jct_reads=1)
        assert drops == {1}

    def test_antisense_read_does_not_support_gtf(self):
        # [P2 fix] a read with the SAME junction coords but on the '-' strand must
        # NOT count as support for a '+' strand GTF sibling -> still dropped.
        cands = [self._anchor("novel"), self._gtf_sibling()]  # both '+'
        em_ab = [30.0, 2.0]
        observed = {"-": {(414, 500): 99}}  # antisense reads at the GTF junction
        drops = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                    gtf_drop_enabled=True, observed_jct=observed,
                                    jct_tol=0, gtf_min_jct_reads=1)
        assert drops == {1}  # opposite-strand support ignored

    def test_gtf_same_chain_diff_ends_kept(self):
        # [P2 fix] GTF shares EVERY junction with the anchor (alt-TSS/TES, no
        # distinguishing junction) -> no failed-support evidence -> keep.
        cands = [self._anchor("novel"),
                 _src_cand("gtf", introns=((200, 300), (400, 500)))]  # identical chain
        em_ab = [30.0, 2.0]
        observed = {"+": {}}  # even with zero observed reads, must NOT drop
        drops = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                    gtf_drop_enabled=True, observed_jct=observed,
                                    jct_tol=0, gtf_min_jct_reads=1)
        assert drops == set()

    def test_loose_tol_lets_gtf_borrow_neighbour_support(self):
        # at jct_tol=2 a read at (412,500) "supports" the GTF's (414,500) -> kept;
        # at jct_tol=0 it does not -> dropped. (Why strict tol matters.)
        cands = [self._anchor("novel"), self._gtf_sibling()]
        em_ab = [30.0, 2.0]
        observed = {"+": {(412, 500): 9}}  # 2bp from the GTF junction
        kept = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                   gtf_drop_enabled=True, observed_jct=observed,
                                   jct_tol=2, gtf_min_jct_reads=1)
        dropped = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                      gtf_drop_enabled=True, observed_jct=observed,
                                      jct_tol=0, gtf_min_jct_reads=1)
        assert kept == set() and dropped == {1}

    def test_gtf_kept_when_drop_disabled(self):
        cands = [self._anchor("novel"), self._gtf_sibling()]
        em_ab = [30.0, 2.0]
        observed = {(200, 300): 9}
        drops = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                    gtf_drop_enabled=False, observed_jct=observed)
        assert drops == set()

    def test_gtf_kept_when_no_observed_jct(self):
        # recall-safe fallback: no read info -> never drop a GTF
        cands = [self._anchor("novel"), self._gtf_sibling()]
        em_ab = [30.0, 2.0]
        drops = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                    gtf_drop_enabled=True, observed_jct=None)
        assert drops == set()

    def test_high_support_gtf_above_frac_not_dropped(self):
        # GTF above frac*anchor is never a shadow regardless of read support
        cands = [self._anchor("novel"), self._gtf_sibling()]
        em_ab = [30.0, 10.0]  # 0.33 > 0.15
        drops = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                    gtf_drop_enabled=True, observed_jct={}, jct_tol=0)
        assert drops == set()

    def test_novel_sibling_dropped_on_abundance_alone(self):
        # novel low-ab sibling dropped by abundance, no read guard needed
        cands = [self._anchor("novel"), _src_cand("novel", introns=((200, 300), (414, 500)))]
        em_ab = [30.0, 2.0]
        drops = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                    gtf_drop_enabled=True, observed_jct=None)
        assert drops == {1}

    def test_fusion_never_dropped(self):
        cands = [self._anchor("novel"), _src_cand("fusion", introns=((200, 300), (414, 500)))]
        em_ab = [30.0, 1.0]
        drops = wobble_shadow_drops(cands, em_ab, self._one_cluster(2), frac=0.15,
                                    gtf_drop_enabled=True, observed_jct={})
        assert drops == set()

    def test_singleton_cluster_no_drop(self):
        cands = [self._anchor("novel")]
        em_ab = [5.0]
        drops = wobble_shadow_drops(cands, em_ab, {0: [0]}, frac=0.15,
                                    gtf_drop_enabled=True, observed_jct={})
        assert drops == set()


class TestGtfGuardNeeded:
    """Pre-check that gates the (BAM-reading) read-support guard."""

    def _one_cluster(self, n):
        return {0: list(range(n))}

    def _anchor(self, source="novel"):
        return _src_cand(source, introns=((200, 300), (400, 500)))

    def _gtf_sibling(self):
        return _src_cand("gtf", introns=((200, 300), (414, 500)))

    def test_needed_when_low_ab_gtf_clustered(self):
        from fin.scoring.m2_junction_nll import gtf_guard_needed
        cands = [self._anchor("novel"), self._gtf_sibling()]
        assert gtf_guard_needed(cands, [30.0, 2.0], self._one_cluster(2), 0.15) is True

    def test_not_needed_when_no_gtf_sibling(self):
        from fin.scoring.m2_junction_nll import gtf_guard_needed
        cands = [self._anchor("novel"), _src_cand("novel", introns=((200, 300), (414, 500)))]
        assert gtf_guard_needed(cands, [30.0, 2.0], self._one_cluster(2), 0.15) is False

    def test_not_needed_when_gtf_above_frac(self):
        from fin.scoring.m2_junction_nll import gtf_guard_needed
        cands = [self._anchor("novel"), self._gtf_sibling()]
        assert gtf_guard_needed(cands, [30.0, 10.0], self._one_cluster(2), 0.15) is False

    def test_not_needed_for_singleton(self):
        from fin.scoring.m2_junction_nll import gtf_guard_needed
        cands = [self._gtf_sibling()]
        assert gtf_guard_needed(cands, [5.0], {0: [0]}, 0.15) is False
