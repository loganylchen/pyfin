"""Unit tests for the collapse step (single-family exact-subchain FOLD).

collapse(family) folds every EXACT contiguous sub-chain into its longest exact container
(absorb-only, fold-all), keeping folded shadows as dormant provenance, and leaves wobble /
cassette / non-exact siblings as separate members. It reproduces the old cluster_read_chains
per-member fold/shadow contract for one family, so the downstream emission, shadow_map and
TSS-recovery contracts are unchanged. EM SCOPE may differ (an exact sub-chain can bridge two
maximals into one family; see TestBridgeAndRobustness) -- an intended behavior change, validated
by metrics.
"""
from fin.candidates.chain_cluster import ChainFamily, collapse

# introns, genomically ordered with positive exons between
I1 = (1000, 2000)
I2 = (3000, 4000)
I3 = (5000, 6000)
I4 = (7000, 8000)


def _fam(variant_reads):
    """variant_reads: {chain_tuple: set(read_ids)} -> a ChainFamily (read_pool = union)."""
    pool = set()
    for r in variant_reads.values():
        pool |= r
    return ChainFamily(variants=sorted(variant_reads), read_pool=pool,
                       variant_reads=dict(variant_reads))


def _chain(member):
    return tuple(member.chain.introns)


class TestExactSubchainFold:
    def test_suffix_subchain_folds_into_longest(self):
        # 5' degradation ladder: (I2,I3) and (I3,) are exact suffix sub-chains of (I1,I2,I3).
        fam = _fam({(I1, I2, I3): {"a"}, (I2, I3): {"b"}, (I3,): {"c"}})
        cl = collapse(fam)
        assert [_chain(m) for m in cl.members] == [(I1, I2, I3)]   # only the maximal survives
        m = cl.members[0]
        assert m.read_ids == {"a", "b", "c"}                       # shadows' reads pooled in
        folded = {tuple(f.chain.introns): f.read_ids for f in m.folded}
        assert folded == {(I2, I3): {"b"}, (I3,): {"c"}}           # shadows kept as provenance

    def test_internal_subchain_folds(self):
        # (I2,) is an exact internal/contiguous sub-chain of (I1,I2,I3).
        fam = _fam({(I1, I2, I3): {"a"}, (I2,): {"b"}})
        cl = collapse(fam)
        assert [_chain(m) for m in cl.members] == [(I1, I2, I3)]
        assert cl.members[0].read_ids == {"a", "b"}

    def test_folds_into_longest_not_intermediate(self):
        # (I3,) is a sub-chain of BOTH (I2,I3) and (I1,I2,I3); must fold into the LONGEST.
        fam = _fam({(I1, I2, I3): {"a"}, (I2, I3): {"b"}, (I3,): {"c"}})
        cl = collapse(fam)
        m = cl.members[0]
        # every non-maximal chain is a shadow of the single survivor (longest container)
        assert {tuple(f.chain.introns) for f in m.folded} == {(I2, I3), (I3,)}


class TestNonExactKept:
    def test_wobble_sibling_not_folded(self):
        # a 5bp-wobbled full chain is NOT an exact sub-chain -> kept as a separate member.
        wob = ((1005, 2000), (3000, 4000), (5000, 6000))
        fam = _fam({(I1, I2, I3): {"a"}, wob: {"b"}})
        cl = collapse(fam)
        chains = sorted(_chain(m) for m in cl.members)
        assert chains == sorted([(I1, I2, I3), wob])               # both survive, no fold

    def test_cassette_sibling_not_folded(self):
        # a cassette-skip sibling (differs by one small exon) is not an exact sub-chain.
        cas = ((1000, 2000), (5000, 6000))                          # skips the I2..I3 exon
        fam = _fam({(I1, I2, I3): {"a"}, cas: {"b"}})
        cl = collapse(fam)
        assert len(cl.members) == 2

    def test_disjoint_variants_all_kept(self):
        fam = _fam({(I1, I2): {"a"}, (I3, I4): {"b"}})
        cl = collapse(fam)
        assert len(cl.members) == 2


class TestInvariants:
    def test_absorb_only_no_new_chains(self):
        fam = _fam({(I1, I2, I3): {"a"}, (I2, I3): {"b"}})
        cl = collapse(fam)
        for m in cl.members:
            assert _chain(m) in fam.variants                        # every survivor is an input chain

    def test_reads_conserved_and_pooled(self):
        fam = _fam({(I1, I2, I3): {"a", "x"}, (I2, I3): {"b"}, (I3,): {"c"}})
        cl = collapse(fam)
        seen = set()
        for m in cl.members:
            seen |= m.read_ids
        assert seen == {"a", "x", "b", "c"} == fam.read_pool       # no read lost or duplicated
        assert cl.read_ids == fam.read_pool

    def test_deterministic_across_input_order(self):
        d1 = {(I1, I2, I3): {"a"}, (I2, I3): {"b"}, (I3,): {"c"}}
        d2 = {(I3,): {"c"}, (I1, I2, I3): {"a"}, (I2, I3): {"b"}}

        def sig(cl):
            return sorted((tuple(m.chain.introns), tuple(sorted(m.read_ids)),
                           tuple(sorted(tuple(f.chain.introns) for f in m.folded)))
                          for m in cl.members)
        assert sig(collapse(_fam(d1))) == sig(collapse(_fam(d2)))

    def test_does_not_mutate_family(self):
        fam = _fam({(I1, I2, I3): {"a"}, (I2, I3): {"b"}})
        before_variants = list(fam.variants)
        before_reads = {c: set(r) for c, r in fam.variant_reads.items()}
        collapse(fam)
        assert fam.variants == before_variants
        assert fam.variant_reads == before_reads

    def test_empty_gtf_only_family(self):
        fam = ChainFamily(variants=[], read_pool=set(), variant_reads={},
                          gtf_members=[("ENST1", (I1, I2))])
        cl = collapse(fam)
        assert cl.members == [] and cl.read_ids == set()


class TestBridgeAndRobustness:
    def test_exact_subchain_bridge_keeps_maximals_coscoped(self):
        # (I2,) is an exact sub-chain of BOTH (I1,I2) and (I2,I3) -> cluster_families bridges the
        # two maximal chains into one family. collapse folds the bridge sub-chain but does NOT
        # re-cluster: both maximals stay as members of the SAME cluster (intended behavior change
        # vs the old fold-first-then-cluster path; Codex gpt-5.6-sol).
        fam = _fam({(I1, I2): {"a"}, (I2, I3): {"b"}, (I2,): {"s"}})
        cl = collapse(fam)
        chains = sorted(_chain(m) for m in cl.members)
        assert chains == [(I1, I2), (I2, I3)]                 # both maximals survive, one cluster
        # the bridge sub-chain folded into exactly one of them (deterministic longest/coord-first)
        holders = [m for m in cl.members if any(f.chain.introns == (I2,) for f in m.folded)]
        assert len(holders) == 1 and "s" in holders[0].read_ids

    def test_deterministic_with_unsorted_variants(self):
        # build ChainFamily directly with variants in a NON-sorted order (bypass _fam's sort).
        vr = {(I1, I2, I3): {"a"}, (I2, I3): {"b"}, (I3,): {"c"}}
        pool = {"a", "b", "c"}
        f_sorted = ChainFamily(variants=sorted(vr), read_pool=pool, variant_reads=dict(vr))
        f_rev = ChainFamily(variants=list(reversed(sorted(vr))), read_pool=pool,
                            variant_reads=dict(vr))

        def sig(cl):
            return sorted((tuple(m.chain.introns), tuple(sorted(m.read_ids)))
                          for m in cl.members)
        assert sig(collapse(f_sorted)) == sig(collapse(f_rev))

    def test_missing_variant_reads_entry_no_crash(self):
        # a maximal variant absent from variant_reads (e.g. a future explore variant whose reads
        # the adapter has not registered) contributes no reads and does not crash.
        fam = ChainFamily(variants=[(I1, I2), (I3, I4)], read_pool={"a"},
                          variant_reads={(I1, I2): {"a"}})   # (I3,I4) has no entry
        cl = collapse(fam)
        by = {_chain(m): m.read_ids for m in cl.members}
        assert by == {(I1, I2): {"a"}, (I3, I4): set()}

    def test_no_read_in_two_members_when_disjoint(self):
        # under the partition precondition, every read lands in exactly one member.
        fam = _fam({(I1, I2, I3): {"a", "b"}, (I4,): {"c"}})   # disjoint reads, disjoint chains
        cl = collapse(fam)
        counts = {}
        for m in cl.members:
            for r in m.read_ids:
                counts[r] = counts.get(r, 0) + 1
        assert all(v == 1 for v in counts.values())


class TestSpanGuard:
    # container (I1,I2,I3); sub-chain (I3,) is a 5' suffix. A read supporting (I3,) that reads
    # exonically across the container's extra intron I2 is a retained-intron / alt isoform and
    # must be KEPT; a 5'-degradation read (starts downstream of I2) folds.
    C = (I1, I2, I3)
    SUB = (I3,)                       # extra introns vs container = I1, I2

    def _fam_sub(self):
        return ChainFamily(
            variants=[self.C, self.SUB],
            read_pool={"full", "deg", "ret"},
            variant_reads={self.C: {"full"}, self.SUB: {"deg", "ret"}})

    def test_span_guard_keeps_retention_read_folds_degradation(self):
        spans = {"full": (900, 6100), "deg": (4500, 6500), "ret": (2500, 6500)}
        cl = collapse(self._fam_sub(), span_guard=True, read_spans=spans)
        by = {tuple(m.chain.introns): m for m in cl.members}
        assert set(by) == {self.C, self.SUB}          # sub-chain survives (retention kept)
        assert by[self.C].read_ids == {"full", "deg"}  # degradation folded into container
        assert by[self.SUB].read_ids == {"ret"}        # retention read kept on the sub-chain
        # only the degradation read is recorded as a folded shadow
        folded = {tuple(f.chain.introns): f.read_ids for f in by[self.C].folded}
        assert folded == {self.SUB: {"deg"}}

    def test_span_guard_off_folds_everything(self):
        spans = {"full": (900, 6100), "deg": (4500, 6500), "ret": (2500, 6500)}
        cl = collapse(self._fam_sub(), span_guard=False, read_spans=spans)
        assert [tuple(m.chain.introns) for m in cl.members] == [self.C]
        assert cl.members[0].read_ids == {"full", "deg", "ret"}

    def test_span_guard_without_spans_is_noop(self):
        # span_guard=True but no read_spans -> cannot judge -> folds unconditionally.
        cl = collapse(self._fam_sub(), span_guard=True, read_spans=None)
        assert [tuple(m.chain.introns) for m in cl.members] == [self.C]
        assert cl.members[0].read_ids == {"full", "deg", "ret"}
