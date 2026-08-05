"""Unit tests for cluster-level diff-region extraction (structural, candidates-only).

Synthetic candidates are built over a fixed genome so each candidate's spliced sequence is exact and
the extracted alleles can be checked against known sub-sequences. The invariants under test:
  - a wobbling junction yields one small diff region with two differing, context-bearing alleles;
  - a skipped (cassette) exon still gives the SKIPPER a non-empty allele (shared context), never
    "nothing covering the region";
  - a pure transcript-end difference (no opposing splice) is NOT a diff region;
  - a candidate that does not reach a locus gets a None allele and marks the region terminal;
  - minus-strand candidates extract consistent (reverse-complemented) alleles.
"""
import random

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.candidates.diff_regions import cluster_diff_regions

_RC = str.maketrans("ACGT", "TGCA")


def _rc(s: str) -> str:
    return s.translate(_RC)[::-1]


_RNG = random.Random(7)                                  # ONE instance: bases actually vary
GENOME = "".join(_RNG.choice("ACGT") for _ in range(300))
assert len(set(GENOME)) == 4                              # guard: a degenerate genome hides flip bugs


def _cand(cid: str, exons, strand: str = "+") -> TranscriptCandidate:
    spliced = "".join(GENOME[s:e] for s, e in exons)
    seq = _rc(spliced) if strand == "-" else spliced
    introns = tuple((exons[k][1], exons[k + 1][0]) for k in range(len(exons) - 1))
    return TranscriptCandidate(
        candidate_id=cid, intron_chain=IntronChain(introns),
        three_prime_pos=(exons[0][0] if strand == "-" else exons[-1][1]),
        sequence=seq, source="novel", supporting_read_ids=set(),
        chrom="chr1", strand=strand, start=exons[0][0], end=exons[-1][1])


class TestWobble:
    def test_single_small_region_two_alleles(self):
        a = _cand("A", [(10, 40), (60, 90), (110, 140)])
        b = _cand("B", [(10, 40), (64, 90), (110, 140)])   # exon2 acceptor 60 -> 64 (4bp wobble)
        regs = cluster_diff_regions([a, b], context=3)
        assert len(regs) == 1
        r = regs[0]
        assert (r.g_lo, r.g_hi) == (60, 64) and r.kind == "wobble"
        # shared 3bp context on both sides; A carries the 4 wobble bp, B splices across them
        assert r.alleles["A"] == GENOME[37:40] + GENOME[60:67]      # ctx + 4bp diff + ctx
        assert r.alleles["B"] == GENOME[37:40] + GENOME[64:67]      # ctx + ctx (no diff bp)
        assert [sorted(g) for g in r.groups] == [["A"], ["B"]]


class TestCassetteSkipHasContext:
    def test_skipper_allele_is_nonempty_shared_context(self):
        # A,B include middle exon [60,90]; C skips it (intron [40,110]).
        a = _cand("A", [(10, 40), (60, 90), (110, 140)])
        b = _cand("B", [(10, 40), (60, 90), (110, 140)])
        c = _cand("C", [(10, 40), (110, 140)])
        regs = cluster_diff_regions([a, b, c], context=3)
        assert len(regs) == 1
        r = regs[0]
        assert (r.g_lo, r.g_hi) == (60, 90) and r.kind == "cassette"
        # THE invariant: the skipper is NOT empty -- it carries the SHARED junction context.
        # C shares exon1 (3bp tail) and exon3 (3bp head) with the includers -> allele = tail + head.
        assert r.alleles["C"] == GENOME[37:40] + GENOME[110:113]
        # includers carry the whole middle exon between the SAME shared flanks and differ from C
        assert r.alleles["A"] == r.alleles["B"] and r.alleles["A"] != r.alleles["C"]
        assert r.alleles["A"] == GENOME[37:40] + GENOME[60:90] + GENOME[110:113]
        assert [sorted(g) for g in r.groups] == [["A", "B"], ["C"]]


class TestEndpointExcluded:
    def test_ragged_end_is_not_a_diff_region(self):
        # identical splicing; B just extends 20bp further at the 3' end -> no opposing splice.
        a = _cand("A", [(10, 40), (60, 90)])
        b = _cand("B", [(10, 40), (60, 110)])
        assert cluster_diff_regions([a, b], context=3) == []


class TestStructuralOutlierExcluded:
    """A candidate that does NOT share both flanking anchors (a truncation / large multi-exon skip) is
    not comparable at the locus and is dropped from the fine region -- it is distinguished by gross
    whole-block presence elsewhere, not by these narrow context-anchored regions."""

    def test_absent_candidate_dropped(self):
        # A/B differ at the last-exon acceptor (110 vs 120); T ends at 90 so it never reaches that
        # locus -> no shared right flank -> excluded. The region is the clean A-vs-B difference.
        a = _cand("A", [(10, 40), (60, 90), (110, 140)])
        b = _cand("B", [(10, 40), (60, 90), (120, 140)])
        t = _cand("T", [(10, 40), (60, 90)])
        regs = cluster_diff_regions([a, b, t], context=3)
        assert len(regs) == 1
        r = regs[0]
        assert (r.g_lo, r.g_hi) == (110, 120)
        assert set(r.alleles) == {"A", "B"}                    # T not comparable -> excluded
        assert r.alleles["A"] != r.alleles["B"]

    def test_truncation_inside_region_dropped(self):
        # A/B differ at exon2's donor (90 vs 75); P ends at 80, strictly inside [75,90], so its right
        # flank is empty -> P does not share flanks -> excluded (not folded into the fine region).
        a = _cand("A", [(10, 40), (60, 90), (110, 140)])
        b = _cand("B", [(10, 40), (60, 75), (110, 140)])
        p = _cand("P", [(10, 40), (60, 80)])
        regs = cluster_diff_regions([a, b, p], context=3)
        assert len(regs) == 1
        assert set(regs[0].alleles) == {"A", "B"}              # P excluded (structural outlier)


class TestMinusStrand:
    def test_cassette_skip_minus_strand(self):
        a = _cand("A", [(10, 40), (60, 90), (110, 140)], strand="-")
        c = _cand("C", [(10, 40), (110, 140)], strand="-")
        regs = cluster_diff_regions([a, c], context=3)
        assert len(regs) == 1
        r = regs[0]
        assert (r.g_lo, r.g_hi) == (60, 90) and r.kind == "cassette"
        # minus-strand alleles = reverse-complement of the plus-strand ones (same shared context)
        assert r.alleles["C"] == _rc(GENOME[37:40] + GENOME[110:113])
        assert r.alleles["A"] == _rc(GENOME[37:40] + GENOME[60:90] + GENOME[110:113])


class TestStrandParity:
    """The engine must give the same STRUCTURE on both strands: a locus on the minus strand yields the
    reverse-complement of the plus-strand allele, with identical grouping."""

    def test_plus_and_minus_are_reverse_complements(self):
        exons = [[(10, 40), (60, 90), (110, 140)], [(10, 40), (64, 90), (110, 140)],
                 [(10, 40), (110, 140)]]
        plus = [_cand(n, e, "+") for n, e in zip("ABC", exons)]
        minus = [_cand(n, e, "-") for n, e in zip("ABC", exons)]
        rp = cluster_diff_regions(plus, context=3)
        rm = cluster_diff_regions(minus, context=3)
        assert len(rp) == len(rm) and len(rp) >= 1
        for a, b in zip(rp, rm):
            assert (a.g_lo, a.g_hi) == (b.g_lo, b.g_hi)          # same genomic loci
            assert a.groups == b.groups                          # same partition of candidates
            for cid in a.alleles:                                # minus allele == rc(plus allele)
                assert b.alleles[cid] == _rc(a.alleles[cid])


class TestDeterministicOrder:
    def test_region_order_independent_of_input_order(self):
        # two flank groups (A/B share exon1 start 10, C/D share start 20) both emit a wobble at the
        # SAME span 60-64 -> the output order must not depend on candidate input order.
        cands = [_cand("A", [(10, 40), (60, 90), (110, 140)]),
                 _cand("B", [(10, 40), (64, 90), (110, 140)]),
                 _cand("C", [(20, 50), (60, 90), (110, 140)]),
                 _cand("D", [(20, 50), (64, 90), (110, 140)])]
        key = lambda regs: [(r.g_lo, r.g_hi, tuple(sorted(r.alleles))) for r in regs]
        forward = key(cluster_diff_regions(cands, context=3))
        assert forward == key(cluster_diff_regions(list(reversed(cands)), context=3))


class TestDegenerate:
    def test_fewer_than_two_multiexon(self):
        a = _cand("A", [(10, 40), (60, 90)])
        assert cluster_diff_regions([a], context=3) == []
        assert cluster_diff_regions([], context=3) == []

    def test_identical_candidates_no_regions(self):
        a = _cand("A", [(10, 40), (60, 90), (110, 140)])
        b = _cand("B", [(10, 40), (60, 90), (110, 140)])
        assert cluster_diff_regions([a, b], context=3) == []
