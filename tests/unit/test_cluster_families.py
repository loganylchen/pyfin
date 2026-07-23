"""Unit tests for the clustering-only primitive ``cluster_families``.

This is the clean, single-responsibility clustering step (clustering redesign): group reads
into multi-exon structural families by wobble/cassette/containment, pool their reads, and
bucket every chain-less read into ``mono_reads`` UNTOUCHED. It must NOT fold exact sub-chains,
NOT process mono, and NOT rewrite coordinates.
"""
from fin.candidates.chain_cluster import ChainFamily, ClusterResult, cluster_families
from fin.candidates.dataclasses import IntronChain


def _rc(qname, introns):
    """(read_dict, IntronChain) with a minimal read_dict (only query_name is read)."""
    return ({"query_name": qname}, IntronChain(introns=tuple(introns)))


def _families_by_variant_count(res):
    return sorted(len(f.variants) for f in res.families)


class TestGroupingAndPool:
    def test_identical_chains_pool_into_one_variant(self):
        res = cluster_families([
            _rc("r1", [(1000, 2000), (3000, 4000)]),
            _rc("r2", [(1000, 2000), (3000, 4000)]),
            _rc("r3", [(1000, 2000), (3000, 4000)]),
        ])
        assert len(res.families) == 1
        fam = res.families[0]
        assert fam.variants == [((1000, 2000), (3000, 4000))]     # one distinct structure
        assert fam.read_pool == {"r1", "r2", "r3"}               # reads pooled at family
        assert res.mono_reads == set()

    def test_unrelated_chains_are_separate_families(self):
        res = cluster_families([
            _rc("a", [(1000, 2000)]),
            _rc("b", [(90000, 91000)]),   # far away -> not related
        ])
        assert len(res.families) == 2
        assert all(len(f.variants) == 1 for f in res.families)


class TestWobbleKeptNotMerged:
    def test_wobble_siblings_one_family_two_variants(self):
        # 2bp-wobbled junctions: same family (wobble), but BOTH variants kept (no merge,
        # no coordinate snapping) and reads pooled.
        res = cluster_families([
            _rc("r1", [(1000, 2000), (3000, 4000)]),
            _rc("r2", [(1002, 1998), (3000, 4000)]),
        ], wobble_bp=6)
        assert len(res.families) == 1
        fam = res.families[0]
        assert set(fam.variants) == {
            ((1000, 2000), (3000, 4000)),
            ((1002, 1998), (3000, 4000)),
        }
        assert fam.read_pool == {"r1", "r2"}
        # coordinates untouched: both original chains still present verbatim
        assert ((1000, 2000), (3000, 4000)) in fam.variant_reads
        assert fam.variant_reads[((1002, 1998), (3000, 4000))] == {"r2"}

    def test_beyond_wobble_tolerance_splits(self):
        # 12bp apart with tol=6 and NO bridge -> two separate families.
        res = cluster_families([
            _rc("r1", [(1000, 2000)]),
            _rc("r2", [(1012, 1988)]),
        ], wobble_bp=6)
        assert len(res.families) == 2


class TestNoFold:
    def test_exact_subchain_is_not_folded_but_clusters_via_containment(self):
        # (3000,4000) is an EXACT sub-chain of the 2-intron chain. cluster_read_chains
        # would FOLD it away; cluster_families must NOT fold -> it stays a distinct variant,
        # joined into the SAME family by tolerant containment.
        res = cluster_families([
            _rc("long", [(1000, 2000), (3000, 4000)]),
            _rc("short", [(3000, 4000)]),
        ])
        assert len(res.families) == 1
        fam = res.families[0]
        assert set(fam.variants) == {
            ((1000, 2000), (3000, 4000)),
            ((3000, 4000),),
        }
        assert fam.read_pool == {"long", "short"}


class TestCassette:
    def test_cassette_sibling_same_family(self):
        # Inclusion isoform (K=3 introns) splits the base's second intron (3000,4000) into
        # (3000,3500) + (3550,4000) with a 50bp (< 70) cassette exon [3500,3550]; outer
        # donor/acceptor match the base intron -> cassette -> same family.
        base = [(1000, 2000), (3000, 4000)]
        cassette = [(1000, 2000), (3000, 3500), (3550, 4000)]
        res = cluster_families([
            _rc("b", base),
            _rc("c", cassette),
        ], wobble_bp=6, cassette_max_exon_bp=70)
        assert len(res.families) == 1
        assert _families_by_variant_count(res) == [2]


class TestMonoBucket:
    def test_monoexon_reads_go_to_bucket_untouched(self):
        res = cluster_families([
            _rc("m1", []),
            _rc("m2", []),
            _rc("multi", [(1000, 2000), (3000, 4000)]),
        ])
        assert res.mono_reads == {"m1", "m2"}     # bucketed, not folded, not split
        assert len(res.families) == 1              # only the multi-exon read forms a family
        assert res.families[0].variants == [((1000, 2000), (3000, 4000))]

    def test_all_mono_no_families(self):
        res = cluster_families([_rc("m1", []), _rc("m2", [])])
        assert res.families == []
        assert res.mono_reads == {"m1", "m2"}


class TestSingleLinkageAndDeterminism:
    def test_single_linkage_bridges_three(self):
        # A~B (2bp) and B~C (2bp) but A vs C is 4bp on the first donor -> still all one
        # family via single-linkage transitivity (accepted by design).
        res = cluster_families([
            _rc("a", [(1000, 2000)]),
            _rc("b", [(1002, 2000)]),
            _rc("c", [(1004, 2000)]),
        ], wobble_bp=2)
        assert len(res.families) == 1
        assert len(res.families[0].variants) == 3
        assert res.families[0].read_pool == {"a", "b", "c"}

    def test_deterministic_across_input_order(self):
        reads_a = [
            _rc("r1", [(1000, 2000), (3000, 4000)]),
            _rc("r2", [(1002, 1998), (3000, 4000)]),
            _rc("r3", [(5000, 6000)]),
        ]
        reads_b = list(reversed(reads_a))
        ra = cluster_families(reads_a)
        rb = cluster_families(reads_b)
        # same families (variant lists coordinate-sorted) and same mono bucket
        assert [f.variants for f in ra.families] == [f.variants for f in rb.families]
        assert ra.mono_reads == rb.mono_reads


class TestSkipsUnnamedReads:
    def test_missing_query_name_skipped(self):
        res = cluster_families([
            ({"query_name": None}, IntronChain(introns=((1000, 2000),))),
            _rc("ok", [(1000, 2000)]),
        ])
        assert len(res.families) == 1
        assert res.families[0].read_pool == {"ok"}


class TestEmptyAndRepresentative:
    def test_empty_input(self):
        res = cluster_families([])
        assert res.families == []
        assert res.mono_reads == set()

    def test_representative_is_longest_no_snap(self):
        # containment family: short (3000,4000) + long 2-intron chain -> representative is
        # the longest variant, with ORIGINAL coordinates (no consensus/snap).
        res = cluster_families([
            _rc("s", [(3000, 4000)]),
            _rc("l", [(1000, 2000), (3000, 4000)]),
        ])
        assert len(res.families) == 1
        assert res.families[0].representative == ((1000, 2000), (3000, 4000))

    def test_representative_tiebreak_by_support(self):
        # two wobble variants of equal length -> representative is the higher-support one.
        res = cluster_families([
            _rc("a1", [(1000, 2000)]),
            _rc("a2", [(1000, 2000)]),   # (1000,2000) has 2 reads
            _rc("b1", [(1003, 2000)]),   # (1003,2000) has 1 read
        ], wobble_bp=6)
        assert len(res.families) == 1
        assert res.families[0].representative == ((1000, 2000),)


class TestGtfAttach:
    def _one(self, res):
        assert len(res.families) == 1
        return res.families[0]

    def test_gtf_wobble_attaches_zero_read(self):
        # GTF wobble-matches the read family -> attached as a zero-read gtf_member; the read
        # pool / variants are UNCHANGED (GTF adds no read mass, no structure).
        gtf = [("gA", ((1002, 1998), (3000, 4000)))]
        res = cluster_families([
            _rc("r1", [(1000, 2000), (3000, 4000)]),
            _rc("r2", [(1000, 2000), (3000, 4000)]),
        ], gtf_variants=gtf)
        fam = self._one(res)
        assert fam.variants == [((1000, 2000), (3000, 4000))]     # read structure unchanged
        assert fam.read_pool == {"r1", "r2"}                      # GTF added no reads
        assert fam.gtf_members == [("gA", ((1002, 1998), (3000, 4000)))]

    def test_gtf_unmatched_becomes_gtf_only_family(self):
        gtf = [("gFar", ((90000, 91000),))]
        res = cluster_families([
            _rc("r1", [(1000, 2000)]),
        ], gtf_variants=gtf)
        assert len(res.families) == 2
        read_fam = [f for f in res.families if f.variants]
        gtf_fam = [f for f in res.families if not f.variants]
        assert len(read_fam) == 1 and len(gtf_fam) == 1
        assert gtf_fam[0].read_pool == set()                     # GTF-only: empty pool
        assert gtf_fam[0].gtf_members == [("gFar", ((90000, 91000),))]
        assert gtf_fam[0].representative is None                 # no read variant

    def test_gtf_does_not_bridge_two_read_families(self):
        # GTF (1000,2000)+(8000,9000) is containment-related to BOTH single-intron read
        # families, which are NOT related to each other. It must attach to exactly ONE and
        # NEVER merge the two families.
        gtf = [("gBridge", ((1000, 2000), (8000, 9000)))]
        res = cluster_families([
            _rc("rA", [(1000, 2000)]),
            _rc("rB", [(8000, 9000)]),
        ], gtf_variants=gtf)
        assert len(res.families) == 2                            # NOT merged
        total_gtf = sum(len(f.gtf_members) for f in res.families)
        assert total_gtf == 1                                    # attached once, not duplicated
        # read families keep their own single reads
        pools = sorted((tuple(sorted(f.read_pool)) for f in res.families))
        assert pools == [("rA",), ("rB",)]

    def test_gtf_none_is_backward_compatible(self):
        res = cluster_families([_rc("r1", [(1000, 2000), (3000, 4000)])], gtf_variants=None)
        assert res.families[0].gtf_members == []

    def test_mono_gtf_is_skipped(self):
        # empty-chain (single-exon) GTF is deferred to the mono finalizer, not attached here.
        res = cluster_families(
            [_rc("r1", [(1000, 2000)])],
            gtf_variants=[("gMono", ())],
        )
        assert len(res.families) == 1
        assert res.families[0].gtf_members == []

    def test_exact_match_gtf_attaches_not_merges_reads(self):
        # a GTF whose chain EXACTLY equals the read chain attaches as a 0-read member (the new
        # design: GTF competes on evidence) -- it does NOT silently absorb the reads at generation.
        gtf = [("gExact", ((1000, 2000), (3000, 4000)))]
        res = cluster_families([
            _rc("r1", [(1000, 2000), (3000, 4000)]),
        ], gtf_variants=gtf)
        fam = self._one(res)
        assert fam.read_pool == {"r1"}
        assert fam.variant_reads == {((1000, 2000), (3000, 4000)): {"r1"}}
        assert fam.gtf_members == [("gExact", ((1000, 2000), (3000, 4000)))]

    def test_deterministic_with_gtf(self):
        reads = [
            _rc("r1", [(1000, 2000), (3000, 4000)]),
            _rc("r2", [(5000, 6000)]),
        ]
        gtf = [("gA", ((1002, 1998), (3000, 4000))), ("gFar", ((90000, 91000),))]
        ra = cluster_families(list(reads), gtf_variants=list(gtf))
        rb = cluster_families(list(reversed(reads)), gtf_variants=list(reversed(gtf)))
        norm = lambda res: [(f.variants, sorted(f.gtf_members)) for f in res.families]
        assert norm(ra) == norm(rb)

    def test_two_gtf_attach_to_same_family(self):
        # two annotation isoforms both wobble the one read family -> both attach as members.
        gtf = [
            ("gA", ((1002, 1998), (3000, 4000))),
            ("gB", ((999, 2001), (3000, 4000))),
        ]
        res = cluster_families([
            _rc("r1", [(1000, 2000), (3000, 4000)]),
        ], gtf_variants=gtf)
        fam = self._one(res)
        assert sorted(fam.gtf_members) == sorted(gtf)
        assert fam.read_pool == {"r1"}

    def test_gtf_only_family_does_not_attract_later_gtf(self):
        # g1 matches no read -> GTF-only family. g2 is a wobble sibling of g1 but ALSO matches
        # no read family; it must NOT attach to g1's GTF-only family -> two GTF-only families.
        gtf = [
            ("g1", ((50000, 51000),)),
            ("g2", ((50002, 50998),)),   # wobble of g1, unrelated to the read
        ]
        res = cluster_families([
            _rc("r1", [(1000, 2000)]),
        ], gtf_variants=gtf)
        assert len(res.families) == 3                      # 1 read family + 2 GTF-only
        gtf_only = [f for f in res.families if not f.variants]
        assert len(gtf_only) == 2
        assert all(len(f.gtf_members) == 1 for f in gtf_only)
