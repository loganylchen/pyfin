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
