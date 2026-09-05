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


class TestGtfVariantKeys:
    """Annotation enters cluster_families as ORDINARY zero-read variant keys (no flag,
    the default): it is listed in ``variants`` and competes like any chain, but it can
    never BRIDGE two read-derived components, and an annotation-only chain related to
    several components joins none of them. ``gtf_members`` records which keys are
    annotated."""

    def _one(self, res):
        assert len(res.families) == 1
        return res.families[0]

    def test_gtf_wobble_joins_as_zero_read_variant(self):
        gtf = [("gA", ((1002, 1998), (3000, 4000)))]
        res = cluster_families([
            _rc("r1", [(1000, 2000), (3000, 4000)]),
            _rc("r2", [(1000, 2000), (3000, 4000)]),
        ], gtf_variants=gtf)
        fam = self._one(res)
        assert fam.variants == [((1000, 2000), (3000, 4000)), ((1002, 1998), (3000, 4000))]
        assert fam.read_pool == {"r1", "r2"}                      # GTF added no reads
        assert fam.variant_reads[((1002, 1998), (3000, 4000))] == set()
        assert fam.gtf_members == [("gA", ((1002, 1998), (3000, 4000)))]

    def test_gtf_unmatched_becomes_gtf_only_family(self):
        gtf = [("gFar", ((90000, 91000),))]
        res = cluster_families([
            _rc("r1", [(1000, 2000)]),
        ], gtf_variants=gtf)
        assert len(res.families) == 2
        gtf_fam = [f for f in res.families if f.read_pool == set()]
        assert len(gtf_fam) == 1
        assert gtf_fam[0].variants == [((90000, 91000),)]      # a real (zero-read) key
        assert gtf_fam[0].gtf_members == [("gFar", ((90000, 91000),))]

    def test_supported_annotation_may_bridge(self):
        # every junction of the annotated chain is observed in reads -> the data
        # agrees the connection exists, so it may merge the two components.
        gtf = [("gBridge", ((1000, 2000), (8000, 9000)))]
        res = cluster_families([
            _rc("rA", [(1000, 2000)]),
            _rc("rB", [(8000, 9000)]),
        ], gtf_variants=gtf,
            observed_junctions={(1000, 2000), (8000, 9000)})
        assert len(res.families) == 1                       # bridged
        fam = res.families[0]
        assert fam.read_pool == {"rA", "rB"}
        assert fam.gtf_members == [("gBridge", ((1000, 2000), (8000, 9000)))]

    def test_uncontradicted_annotation_may_bridge(self):
        # the annotation's extra junction is unobserved but nothing near it was
        # observed either: the data is silent, not disagreeing -> it may bridge.
        gtf = [("gQuiet", ((1000, 2000), (5000, 6000), (8000, 9000)))]
        res = cluster_families([
            _rc("rA", [(1000, 2000)]),
            _rc("rB", [(8000, 9000)]),
        ], gtf_variants=gtf, observed_junctions={(1000, 2000), (8000, 9000)})
        assert len(res.families) == 1
        assert res.families[0].read_pool == {"rA", "rB"}

    def test_nearby_different_junction_contradicts(self):
        gtf = [("gJit", ((1004, 2000),))]
        res = cluster_families([
            _rc("rA", [(1000, 2000)]),
        ], gtf_variants=gtf, observed_junctions={(1000, 2000)})
        fam = next(f for f in res.families if f.gtf_members)
        assert ((1004, 2000),) not in fam.uncontradicted_gtf

    def test_silent_region_is_not_contradiction(self):
        gtf = [("gFar", ((50000, 51000),))]
        res = cluster_families([
            _rc("rA", [(1000, 2000)]),
        ], gtf_variants=gtf, observed_junctions={(1000, 2000)})
        fam = next(f for f in res.families if f.gtf_members)
        assert ((50000, 51000),) in fam.uncontradicted_gtf

    def test_contradicted_junction_may_not_bridge(self):
        # a read puts a DIFFERENT junction at the annotation's donor (a
        # jitter/merge signature) -> the annotation may not reshape the topology.
        gtf = [("gJitter", ((1000, 2000), (5000, 6000), (8000, 9000)))]
        res = cluster_families([
            _rc("rA", [(1000, 2000)]),
            _rc("rB", [(8000, 9000)]),
            _rc("rC", [(5000, 6400)]),          # same donor, different acceptor
        ], gtf_variants=gtf,
            observed_junctions={(1000, 2000), (8000, 9000), (5000, 6400)})
        fam = next(f for f in res.families if f.gtf_members)
        assert ((1000, 2000), (5000, 6000), (8000, 9000)) not in fam.uncontradicted_gtf

    def test_fabricated_junction_may_not_bridge(self):
        # the annotation carries a junction no read ever produced (a merge/jitter
        # corruption): it must not weld the two real families together.
        gtf = [("gMerged", ((1000, 2000), (5000, 6000), (8000, 9000)))]
        res = cluster_families([
            _rc("rA", [(1000, 2000)]),
            _rc("rB", [(8000, 9000)]),
        ], gtf_variants=gtf,
            observed_junctions={(1000, 2000), (8000, 9000), (5003, 6002)})
        assert len(res.families) >= 2
        assert sorted(sorted(f.read_pool) for f in res.families) == [[], ["rA"], ["rB"]]

    def test_bridging_is_off_without_a_registry(self):
        gtf = [("gBridge", ((1000, 2000), (8000, 9000)))]
        res = cluster_families([
            _rc("rA", [(1000, 2000)]),
            _rc("rB", [(8000, 9000)]),
        ], gtf_variants=gtf)                    # no observed_junctions
        assert len(res.families) == 3

    def test_bridging_is_order_invariant(self):
        gtf = [("gBridge", ((1000, 2000), (8000, 9000)))]
        obs = {(1000, 2000), (8000, 9000)}
        reads = [_rc("rA", [(1000, 2000)]), _rc("rB", [(8000, 9000)])]
        key = lambda r: sorted((tuple(sorted(f.read_pool)),
                                tuple(sorted(g for g, _ in f.gtf_members)))
                               for f in r.families)
        a = cluster_families(list(reads), gtf_variants=gtf, observed_junctions=obs)
        b = cluster_families(list(reversed(reads)), gtf_variants=gtf,
                             observed_junctions=obs)
        assert key(a) == key(b)

    def test_gtf_never_bridges_two_read_families(self):
        # (1000,2000)+(8000,9000) is containment-related to BOTH single-intron read
        # chains, which are unrelated to each other. A zero-read annotation must not
        # merge them: it attaches to exactly one component and the data's partition
        # is unchanged.
        gtf = [("gBridge", ((1000, 2000), (8000, 9000)))]
        res = cluster_families([
            _rc("rA", [(1000, 2000)]),
            _rc("rB", [(8000, 9000)]),
        ], gtf_variants=gtf)
        assert len(res.families) == 3                      # rA, rB, and the annotation
        assert sorted(sorted(f.read_pool) for f in res.families) == [[], ["rA"], ["rB"]]
        assert sum(len(f.gtf_members) for f in res.families) == 1   # counted once
        host = next(f for f in res.families if f.gtf_members)
        assert ((1000, 2000), (8000, 9000)) in host.variants        # a real variant key

    def test_ambiguous_annotation_attaches_to_no_component(self):
        # related to BOTH read components -> nothing shows which it belongs to, so it
        # stays its own zero-read component instead of being placed by a tie-break.
        gtf = [("gB", ((1000, 2000), (8000, 9000)))]
        reads = [_rc("rA", [(1000, 2000)]), _rc("rB", [(8000, 9000)])]
        res = cluster_families(list(reads), gtf_variants=gtf)
        assert len(res.families) == 3
        host = next(f for f in res.families if f.gtf_members)
        assert host.read_pool == set()
        assert host.variants == [((1000, 2000), (8000, 9000))]
        # ...and the answer does not depend on input order
        rev = cluster_families(list(reversed(reads)), gtf_variants=gtf)
        key = lambda r: sorted((tuple(sorted(f.read_pool)),
                                tuple(sorted(g for g, _ in f.gtf_members)))
                               for f in r.families)
        assert key(res) == key(rev)

    def test_unambiguous_annotation_joins_its_single_component(self):
        gtf = [("gA", ((1002, 1998), (3000, 4000)))]
        res = cluster_families([
            _rc("r1", [(1000, 2000), (3000, 4000)]),
            _rc("r2", [(5000, 6000), (7000, 8000)]),
        ], gtf_variants=gtf)
        host = next(f for f in res.families if f.gtf_members)
        assert host.read_pool == {"r1"}
        assert ((1002, 1998), (3000, 4000)) in host.variants

    def test_two_gtf_only_chains_do_not_form_a_bridge_chain(self):
        # g1 relates to rA, g2 relates to rB, and g1~g2: allowing annotation-annotation
        # edges would merge rA and rB one level down. They must stay separate.
        gtf = [("g1", ((1000, 2000), (5000, 6000))),
               ("g2", ((5000, 6000), (8000, 9000)))]
        res = cluster_families([
            _rc("rA", [(1000, 2000)]),
            _rc("rB", [(8000, 9000)]),
        ], gtf_variants=gtf)
        assert len(res.families) == 2
        assert sorted(sorted(f.read_pool) for f in res.families) == [["rA"], ["rB"]]

    def test_gtf_none_is_backward_compatible(self):
        res = cluster_families([_rc("r1", [(1000, 2000), (3000, 4000)])], gtf_variants=None)
        assert res.families[0].gtf_members == []

    def test_mono_gtf_is_skipped(self):
        # empty-chain (single-exon) GTF is deferred to the mono finalizer, not a key here.
        res = cluster_families(
            [_rc("r1", [(1000, 2000)])],
            gtf_variants=[("gMono", ())],
        )
        assert len(res.families) == 1
        assert res.families[0].gtf_members == []
        assert res.families[0].variants == [((1000, 2000),)]

    def test_exact_match_gtf_is_the_same_key(self):
        # a GTF whose chain EXACTLY equals a read chain is the SAME variant key: the
        # key carries the reads and is annotated -- no duplicate, no read absorption.
        gtf = [("gExact", ((1000, 2000), (3000, 4000)))]
        res = cluster_families([
            _rc("r1", [(1000, 2000), (3000, 4000)]),
        ], gtf_variants=gtf)
        fam = self._one(res)
        assert fam.variants == [((1000, 2000), (3000, 4000))]
        assert fam.variant_reads == {((1000, 2000), (3000, 4000)): {"r1"}}
        assert fam.gtf_members == [("gExact", ((1000, 2000), (3000, 4000)))]

    def test_two_gtf_models_share_one_key(self):
        # two annotated transcripts with the same chain (alternative TSS/TES) -> one key,
        # two gtf_members entries.
        gtf = [("gA", ((1000, 2000),)), ("gB", ((1000, 2000),))]
        res = cluster_families([_rc("r1", [(1000, 2000)])], gtf_variants=gtf)
        fam = self._one(res)
        assert fam.variants == [((1000, 2000),)]
        assert fam.gtf_members == [("gA", ((1000, 2000),)), ("gB", ((1000, 2000),))]

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

    def test_two_gtf_wobble_siblings_join_same_family(self):
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
        assert len(fam.variants) == 3

    def test_unrelated_gtf_only_chains_stay_separate(self):
        # neither annotation relates to a read chain: each is its own zero-read
        # component (they are not joined to each other either).
        gtf = [("g1", ((50000, 51000),)), ("g2", ((50002, 50998),))]
        res = cluster_families([
            _rc("r1", [(1000, 2000)]),
        ], gtf_variants=gtf)
        gtf_only = [f for f in res.families if f.read_pool == set()]
        assert len(gtf_only) == 2
        assert all(len(f.gtf_members) == 1 for f in gtf_only)


class TestContradictionRegistryCoversEveryReadChain:
    """The contradiction evidence must come from every read chain that enters
    clustering. Using only the MAPQ>=1 SUPPORT registry would let a low-MAPQ read
    shape the read-derived topology while being invisible to the gate, so its
    junction would count as 'the reads are silent' when they in fact disagree.
    """

    def test_a_read_junction_absent_from_the_registry_still_contradicts(self):
        # the registry omits (5000,6000) -> the annotation would look uncontradicted
        gtf = [("gMerge", ((1000, 2000), (5000, 6000), (8000, 9000)))]
        reads = [_rc("rA", [(1000, 2000)]), _rc("rB", [(8000, 9000)]),
                 _rc("rC", [(5000, 6400)])]
        # registry built from EVERY read chain (what discovery now passes)
        obs = {j for _, ch in [(None, type("C", (), {"introns": c})()) for c in
                               ([(1000, 2000)], [(8000, 9000)], [(5000, 6400)])]
               for j in ch.introns}
        res = cluster_families(reads, gtf_variants=gtf, observed_junctions=obs)
        fam = next(f for f in res.families if f.gtf_members)
        assert ((1000, 2000), (5000, 6000), (8000, 9000)) not in fam.uncontradicted_gtf
