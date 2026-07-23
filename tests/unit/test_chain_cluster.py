"""Unit tests for generation-side chain clustering (chain_cluster).

Group reads by EXACT chain (3' ignored); fold EXACT sub-chains into their container;
union-find the survivors into clusters by wobble / cassette / containment. Members
(distinct candidates) are kept; no snap, no shadow generation.
"""
from __future__ import annotations

from fin.candidates.dataclasses import IntronChain
from fin.candidates.chain_cluster import (
    cluster_read_chains, _exact_subchain, _wobble, _contains, _cassette, _related,
)

J1, J2, J3 = (100, 200), (300, 400), (500, 600)


def _rc(pairs):
    return [({"query_name": rid}, IntronChain(introns=introns)) for rid, introns in pairs]


def _member_chains(cluster):
    return {m.chain.introns for m in cluster.members}


class TestPredicates:
    def test_exact_subchain(self):
        assert _exact_subchain((J2, J3), (J1, J2, J3))
        assert not _exact_subchain(((302, 398), J3), (J1, J2, J3))  # wobbled -> not exact
        assert not _exact_subchain((J1, J2, J3), (J1, J2, J3))       # equal length

    def test_wobble(self):
        assert _wobble((J1, J2), ((102, 198), (301, 402)), bp=6)
        assert not _wobble((J1, J2), (J1, J2, J3), bp=6)

    def test_contains_tolerant(self):
        assert _contains(((302, 398), J3), (J1, J2, J3), bp=6)   # wobbled suffix
        assert not _contains(((330, 370), J3), (J1, J2, J3), bp=6)

    def test_cassette(self):
        long_ = (J1, (300, 350), (370, 600))    # small exon 350..370
        short = (J1, (300, 600))
        assert _cassette(long_, short, bp=6, max_exon_bp=70)
        assert not _cassette(long_, short, bp=6, max_exon_bp=10)


class TestClustering:
    def test_low_support_exact_subchain_collapses(self):
        # exact suffix with FEW reads (1) vs a well-supported container (4 reads):
        # folded away (a truncation shadow).
        cl = cluster_read_chains(_rc(
            [("a1", (J1, J2, J3)), ("a2", (J1, J2, J3)), ("a3", (J1, J2, J3)),
             ("a4", (J1, J2, J3)), ("b", (J2, J3))]))
        assert len(cl) == 1
        assert _member_chains(cl[0]) == {(J1, J2, J3)}   # only the container survives
        assert cl[0].read_ids == {"a1", "a2", "a3", "a4", "b"}

    def test_exact_subchain_always_folds_regardless_of_support(self):
        # UNCONDITIONAL fold: an exact contiguous sub-chain is folded into its
        # container even when well-supported. Whether it is a genuine short isoform
        # is NOT decided here -- that is a separate downstream recovery step.
        cl = cluster_read_chains(_rc(
            [("a1", (J1, J2, J3)), ("a2", (J1, J2, J3)),
             ("b1", (J2, J3)), ("b2", (J2, J3))]))
        assert len(cl) == 1
        assert _member_chains(cl[0]) == {(J1, J2, J3)}   # sub-chain folded away
        assert cl[0].read_ids == {"a1", "a2", "b1", "b2"}  # its reads pooled in

    def test_folded_subchain_retained_as_shadow(self):
        # The folded sub-chain is kept as a SHADOW (provenance) on its container:
        # its own chain + its own reads, for the downstream recovery step.
        cl = cluster_read_chains(_rc(
            [("a1", (J1, J2, J3)), ("a2", (J1, J2, J3)),
             ("b1", (J2, J3)), ("b2", (J2, J3))]))
        (m,) = cl[0].members
        assert m.chain.introns == (J1, J2, J3)
        assert len(m.folded) == 1
        assert m.folded[0].chain.introns == (J2, J3)
        assert m.folded[0].read_ids == {"b1", "b2"}      # shadow's own reads
        assert m.read_ids == {"a1", "a2", "b1", "b2"}    # container still pools all

    def test_wobble_variants_same_cluster_kept(self):
        # 2bp-wobbled sibling: NOT an exact subchain -> both KEPT, one cluster.
        cl = cluster_read_chains(
            _rc([("a", (J1, J2, J3)), ("b", (J1, (302, 398), J3))]))
        assert len(cl) == 1
        assert _member_chains(cl[0]) == {(J1, J2, J3), (J1, (302, 398), J3)}
        assert cl[0].read_ids == {"a", "b"}

    def test_cassette_same_cluster_kept(self):
        long_ = (J1, (300, 350), (370, 600))
        short = (J1, (300, 600))
        cl = cluster_read_chains(_rc([("a", long_), ("b", short)]))
        assert len(cl) == 1
        assert _member_chains(cl[0]) == {long_, short}   # both kept

    def test_tolerant_containment_same_cluster_kept(self):
        # a wobbled (non-exact) sub-chain -> clustered but KEPT (not collapsed).
        cl = cluster_read_chains(
            _rc([("a", (J1, J2, J3)), ("b", ((302, 398), J3))]))
        assert len(cl) == 1
        assert _member_chains(cl[0]) == {(J1, J2, J3), ((302, 398), J3)}

    def test_distinct_chains_separate(self):
        far = ((1000, 1100), (1300, 1400))
        cl = cluster_read_chains(_rc([("a", (J1, J2, J3)), ("b", far)]))
        assert len(cl) == 2

    def test_three_prime_ignored_same_chain_pools(self):
        cl = cluster_read_chains(
            _rc([("a", (J1, J2)), ("b", (J1, J2)), ("c", (J1, J2))]))
        assert len(cl) == 1
        assert len(cl[0].members) == 1
        assert cl[0].read_ids == {"a", "b", "c"}

    def test_mono_exon_singleton(self):
        cl = cluster_read_chains(_rc([("a", ()), ("b", ()), ("m", (J1, J2))]))
        by = {tuple(c.members[0].chain.introns): c for c in cl}
        assert by[()].read_ids == {"a", "b"}
        assert by[(J1, J2)].read_ids == {"m"}
        assert len(cl) == 2

    def test_representative_is_longest(self):
        # tolerant-containment cluster: 3-intron chain vs a wobbled 2-intron sub-chain
        # (kept, not collapsed) -> representative is the longer (3-intron) member.
        cl = cluster_read_chains(
            _rc([("a", (J1, J2, J3)), ("b", ((302, 398), J3))]))
        assert cl[0].representative.chain.introns == (J1, J2, J3)


def _rcs(entries):
    """entries: (rid, introns, ref_start, ref_end) -> read_chains with spans."""
    return [({"query_name": rid, "reference_start": s, "reference_end": e},
             IntronChain(introns=introns)) for rid, introns, s, e in entries]


class TestMonoexonFold:
    """fold_monoexon_contained: single-exon reads inside a multi-exon candidate's exon
    are folded into it; intronic / uncontained / junction-crossing mono reads stay."""

    # multi candidate (J1,J2): exons [50,100) [200,300) [400,500); introns [100,200) [300,400)
    def _rc(self, mono):
        return _rcs([("m1", (J1, J2), 50, 500), ("m2", (J1, J2), 50, 500)] + mono)

    def test_monoexon_in_exon_folds(self):
        cl = cluster_read_chains(self._rc([("s1", (), 210, 290)]),
                                 fold_monoexon_contained=True)
        assert len(cl) == 1                                   # no standalone mono
        assert _member_chains(cl[0]) == {(J1, J2)}
        assert "s1" in cl[0].read_ids                         # folded into the multi

    def test_monoexon_in_intron_stays(self):
        cl = cluster_read_chains(self._rc([("s1", (), 120, 180)]),  # inside intron J1
                                 fold_monoexon_contained=True)
        assert any(m.chain.introns == () for c in cl for m in c.members)

    def test_monoexon_outside_span_stays(self):
        cl = cluster_read_chains(self._rc([("s1", (), 600, 700)]),
                                 fold_monoexon_contained=True)
        assert any(m.chain.introns == () for c in cl for m in c.members)

    def test_monoexon_crossing_junction_stays(self):
        cl = cluster_read_chains(self._rc([("s1", (), 250, 350)]),  # spans intron J2
                                 fold_monoexon_contained=True)
        assert any(m.chain.introns == () for c in cl for m in c.members)

    def test_default_off_emits_standalone_mono(self):
        cl = cluster_read_chains(self._rc([("s1", (), 210, 290)]))  # flag off
        assert any(m.chain.introns == () for c in cl for m in c.members)


class TestFoldSpanGuard:
    """fold_span_guard: a read whose span runs exonically across a container's EXTRA
    intron (retained-intron / alt isoform) is kept as its own candidate, not folded."""

    # container (J1,J2,J3): extra intron vs (J1,J2) is J3=(500,600).
    def _rc(self):
        return _rcs([
            ("big1", (J1, J2, J3), 50, 650),
            ("big2", (J1, J2, J3), 50, 650),
            ("trunc", (J1, J2), 50, 450),     # ends before J3 -> truncation (fold)
            ("retain", (J1, J2), 50, 650),    # spans across J3 -> retained (keep)
        ])

    def test_guard_keeps_retained_intron_folds_truncation(self):
        cl = cluster_read_chains(self._rc(), fold_span_guard=True)
        members = {m.chain.introns: m for c in cl for m in c.members}
        assert (J1, J2, J3) in members
        assert (J1, J2) in members                          # retained-intron kept
        assert "retain" in members[(J1, J2)].read_ids
        assert "trunc" in members[(J1, J2, J3)].read_ids     # truncation folded in
        assert "retain" not in members[(J1, J2, J3)].read_ids

    def test_default_off_folds_everything(self):
        cl = cluster_read_chains(self._rc())                # guard off (default)
        members = {m.chain.introns for c in cl for m in c.members}
        assert (J1, J2) not in members                       # unconditional fold
        assert (J1, J2, J3) in members
