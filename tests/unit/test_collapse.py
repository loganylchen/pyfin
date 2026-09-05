"""Unit tests for the collapse step (single-family exact-subchain FOLD).

collapse(family) folds every EXACT contiguous sub-chain into its longest exact container
(absorb-only, fold-all), keeping folded shadows as dormant provenance, and leaves wobble /
cassette / non-exact siblings as separate members. It reproduces the old cluster_read_chains
per-member fold/shadow contract for one family, so the downstream emission, shadow_map and
TSS-recovery contracts are unchanged. EM SCOPE may differ (an exact sub-chain can bridge two
maximals into one family; see TestBridgeAndRobustness) -- an intended behavior change, validated
by metrics.
"""
from fin.candidates.chain_cluster import (
    ChainFamily, _cigar_exonically_covers, collapse, read_alignment_from_cigar)

# introns, genomically ordered with positive exons between
I1 = (1000, 2000)
I2 = (3000, 4000)
I3 = (5000, 6000)
I4 = (7000, 8000)


def _aln(ref_start, cigar):
    return read_alignment_from_cigar(cigar, ref_start)


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

    def test_gtf_only_family_keeps_zero_read_key(self):
        # an annotated chain with no reads is an ordinary (zero-read) variant: it is
        # kept as a member so discovery can emit the annotation model from it.
        fam = ChainFamily(variants=[(I1, I2)], read_pool=set(),
                          variant_reads={(I1, I2): set()},
                          gtf_members=[("ENST1", (I1, I2))])
        cl = collapse(fam)
        assert [tuple(m.chain.introns) for m in cl.members] == [(I1, I2)]
        assert cl.members[0].read_ids == set() and cl.read_ids == set()

    def test_annotated_subchain_is_never_a_foldee(self):
        # an annotated exact sub-chain (zero reads OR with reads) stays a member: the
        # fold is a read-level degradation call a true short annotated isoform can
        # never refute. Measured -3.0 F1 when folded (SIRV4 grid 2026-09-03).
        fam = ChainFamily(variants=[(I1, I2, I3), (I2, I3)], read_pool={"a"},
                          variant_reads={(I1, I2, I3): {"a"}, (I2, I3): set()},
                          gtf_members=[("ENST_short", (I2, I3))])
        cl = collapse(fam, span_guard=True,
                      read_alignments={"a": _aln(900, [(0, 8100)])})
        assert [tuple(m.chain.introns) for m in cl.members] == [(I1, I2, I3), (I2, I3)]
        assert cl.members[0].folded == []
        fam2 = ChainFamily(variants=[(I1, I2, I3), (I2, I3)], read_pool={"a", "d"},
                           variant_reads={(I1, I2, I3): {"a"}, (I2, I3): {"d"}},
                           gtf_members=[("ENST_short", (I2, I3))])
        cl2 = collapse(fam2, span_guard=True,
                       read_alignments={"a": _aln(900, [(0, 8100)]),
                                        "d": _aln(4500, [(0, 500), (3, 1000), (0, 500)])})
        by = {tuple(m.chain.introns): m for m in cl2.members}
        assert by[(I2, I3)].read_ids == {"d"} and by[(I1, I2, I3)].read_ids == {"a"}

    def test_uncontradicted_annotation_may_contain(self):
        # the reads do not contradict the annotated chain -> it may absorb the
        # 5'-truncated read sub-chain (those reads are fragments of it).
        fam = ChainFamily(variants=[(I1, I2, I3), (I2, I3)], read_pool={"b"},
                          variant_reads={(I1, I2, I3): set(), (I2, I3): {"b"}},
                          gtf_members=[("ENST_long", (I1, I2, I3))],
                          uncontradicted_gtf={(I1, I2, I3)})
        cl = collapse(fam)
        assert [tuple(m.chain.introns) for m in cl.members] == [(I1, I2, I3)]
        assert cl.members[0].read_ids == {"b"}

    def test_annotation_only_key_is_not_a_fold_container(self):
        # a zero-read annotated container must NOT absorb a read-supported sub-chain:
        # a hypothesis with no reads cannot delete one that has them before any
        # evidence is weighed. Both survive and compete in assignment.
        fam = ChainFamily(variants=[(I1, I2, I3), (I2, I3)], read_pool={"b"},
                          variant_reads={(I1, I2, I3): set(), (I2, I3): {"b"}},
                          gtf_members=[("ENST_long", (I1, I2, I3))])
        cl = collapse(fam)
        assert [tuple(m.chain.introns) for m in cl.members] == [(I1, I2, I3), (I2, I3)]
        assert {tuple(m.chain.introns): m.read_ids for m in cl.members}[(I2, I3)] == {"b"}

    def test_annotated_key_with_its_own_reads_may_contain(self):
        # the annotated chain is also a READ chain here -> an ordinary key, so the
        # exact sub-chain folds into it as usual.
        fam = ChainFamily(variants=[(I1, I2, I3), (I2, I3)], read_pool={"a", "b"},
                          variant_reads={(I1, I2, I3): {"a"}, (I2, I3): {"b"}},
                          gtf_members=[("ENST_long", (I1, I2, I3))])
        cl = collapse(fam)
        assert [tuple(m.chain.introns) for m in cl.members] == [(I1, I2, I3)]
        assert cl.members[0].read_ids == {"a", "b"}


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

    def _alns(self):
        # full: exact container read. deg: 5'-degraded, starts inside exon 3 (after I2).
        # ret: reads exonically through I2 = (3000,4000): continuous M from 2500 to 6500
        # apart from the I3 splice.
        return {
            "full": _aln(900, [(0, 100), (3, 1000), (0, 1000), (3, 1000), (0, 1000),
                               (3, 1000), (0, 100)]),
            "deg": _aln(4500, [(0, 500), (3, 1000), (0, 500)]),
            "ret": _aln(2500, [(0, 2500), (3, 1000), (0, 500)]),
        }

    def test_span_guard_keeps_retention_read_folds_degradation(self):
        cl = collapse(self._fam_sub(), span_guard=True, read_alignments=self._alns())
        by = {tuple(m.chain.introns): m for m in cl.members}
        assert set(by) == {self.C, self.SUB}          # sub-chain survives (retention kept)
        assert by[self.C].read_ids == {"full", "deg"}  # degradation folded into container
        assert by[self.SUB].read_ids == {"ret"}        # retention read kept on the sub-chain
        # only the degradation read is recorded as a folded shadow
        folded = {tuple(f.chain.introns): f.read_ids for f in by[self.C].folded}
        assert folded == {self.SUB: {"deg"}}

    def test_span_guard_off_folds_everything(self):
        cl = collapse(self._fam_sub(), span_guard=False, read_alignments=self._alns())
        assert [tuple(m.chain.introns) for m in cl.members] == [self.C]
        assert cl.members[0].read_ids == {"full", "deg", "ret"}

    def test_span_guard_without_alignments_is_noop(self):
        # span_guard=True but no read_alignments -> cannot judge -> folds unconditionally.
        cl = collapse(self._fam_sub(), span_guard=True, read_alignments=None)
        assert [tuple(m.chain.introns) for m in cl.members] == [self.C]
        assert cl.members[0].read_ids == {"full", "deg", "ret"}

    def test_span_guard_reference_span_alone_is_not_read_through(self):
        # `ret` brackets I2 by reference span but jumps it with a 900 bp deletion -> fold.
        alns = self._alns()
        alns["ret"] = _aln(2500, [(0, 550), (2, 900), (0, 1050), (3, 1000), (0, 500)])
        cl = collapse(self._fam_sub(), span_guard=True, read_alignments=alns)
        assert [tuple(m.chain.introns) for m in cl.members] == [self.C]


class TestCigarExonicallyCovers:
    """Unit rules of the CIGAR read-through test on intron [3000,4000), W=[2990,4010)."""
    D, A = 3000, 4000

    def _ok(self, aln):
        return _cigar_exonically_covers(aln, self.D, self.A)

    def test_continuous_match_keeps(self):
        assert self._ok(_aln(2500, [(0, 2000)]))

    def test_small_deletions_keep(self):
        for dlen in (1, 2, 3):
            assert self._ok(_aln(2500, [(0, 800), (2, dlen), (0, 1200)])), dlen

    def test_deletion_of_four_folds(self):
        assert not self._ok(_aln(2500, [(0, 800), (2, 4), (0, 1200)]))

    def test_many_small_deletions_over_total_fold(self):
        # 3 x 3 bp = 9 bp > 6 bp total
        assert not self._ok(_aln(2500, [(0, 700), (2, 3), (0, 200), (2, 3), (0, 200),
                                        (2, 3), (0, 700)]))
        # 2 x 3 bp = 6 bp is still allowed
        assert self._ok(_aln(2500, [(0, 700), (2, 3), (0, 300), (2, 3), (0, 800)]))

    def test_skip_across_intron_folds(self):
        assert not self._ok(_aln(2500, [(0, 500), (3, 1000), (0, 500)]))
        # an N that only clips the window edge also disqualifies
        assert not self._ok(_aln(2500, [(0, 495), (3, 10), (0, 1500)]))

    def test_large_deletion_with_bracketing_span_folds(self):
        assert not self._ok(_aln(2500, [(0, 400), (2, 1200), (0, 400)]))

    def test_soft_clip_leaving_upstream_flank_unaligned_folds(self):
        # reference_start 2995 -> only 5 bp of the 10 bp upstream flank are matched
        assert not self._ok(_aln(2995, [(4, 300), (0, 1500)]))
        # the footprint must cover ALL of W: starting 2 bp inside the flank still folds
        assert not self._ok(_aln(2992, [(0, 1500)]))
        # ...but a 2 bp deletion INSIDE the flank (footprint intact, 8 matched) keeps
        assert self._ok(_aln(2990, [(0, 4), (2, 2), (0, 1500)]))

    def test_one_sided_flank_folds(self):
        assert not self._ok(_aln(2500, [(0, 1505)]))      # ends at 4005: downstream flank short
        assert not self._ok(_aln(3005, [(0, 1500)]))      # starts inside the intron

    def test_multiple_extra_introns_any_read_through_keeps(self):
        fam = ChainFamily(
            variants=[(I1, I2, I3), (I3,)], read_pool={"c", "r"},
            variant_reads={(I1, I2, I3): {"c"}, (I3,): {"r"}})
        # r splices I1 like the container but reads through I2 -> genuinely distinct
        alns = {"r": _aln(900, [(0, 100), (3, 1000), (0, 3000), (3, 1000), (0, 100)])}
        cl = collapse(fam, span_guard=True, read_alignments=alns)
        assert {tuple(m.chain.introns) for m in cl.members} == {(I1, I2, I3), (I3,)}
        assert next(m for m in cl.members if tuple(m.chain.introns) == (I3,)).read_ids == {"r"}
