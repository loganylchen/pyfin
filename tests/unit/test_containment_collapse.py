"""Unit tests for Lever 1 — containment / 5'-truncation collapse.

``containment_shadow_drops`` maps a NOVEL candidate whose intron chain is a pure,
exact 3' SUFFIX of a longer candidate (a 5'-truncation shadow) to the longer
parent it should fold into. Strand sets which genomic end is the 3' terminus:
"+" -> high coord (suffix = last introns), "-" -> low coord (suffix = first
introns). The match is strict, so exon-skipping / alt-internal-splicing isoforms
never fold. Only novels are foldable; a gtf candidate may be a parent but is never
folded away. A shadow folds only when its EM abundance <= parent * min_ratio.

The aggregate-layer test confirms a single-interval fold (parent.abundance +=
shadow + read-id union) survives ``aggregate_across_intervals`` unchanged.
"""
from __future__ import annotations

import uuid

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.analysis.quantification import QuantResult, aggregate_across_intervals
from fin.scoring.m2_junction_nll import containment_shadow_drops


def _cand(introns, *, start, end, source="novel", strand="+", chrom="chr1",
          cid=None):
    return TranscriptCandidate(
        candidate_id=cid or ("novel_" + uuid.uuid4().hex[:8]),
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=end,
        sequence="ACGT" * 50,
        source=source,
        supporting_read_ids=set(),
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
    )


class TestPureTruncationFolds:
    def test_plus_strand_5p_truncation(self):
        # parent: 2 introns spanning 100..700; shadow: drops the 5' intron, keeps
        # the 3' (high-coord) intron + same 3' end; 5' start interior.
        parent = _cand([(100, 200), (400, 500)], start=50, end=700)
        shadow = _cand([(400, 500)], start=300, end=700)
        em_ab = [10.0, 3.0]  # parent, shadow
        fold = containment_shadow_drops(
            [parent, shadow], em_ab, tol_bp=20, min_ratio=1.0
        )
        assert fold == {1: 0}

    def test_minus_strand_5p_truncation(self):
        # "-" strand: 3' = low coord. shadow shares the FIRST (low-coord) intron
        # + same 3' (start) end; shadow's 5' (high-coord, end) interior to parent.
        parent = _cand([(100, 200), (400, 500)], start=50, end=700, strand="-")
        shadow = _cand([(100, 200)], start=50, end=350, strand="-")
        em_ab = [10.0, 3.0]
        fold = containment_shadow_drops(
            [parent, shadow], em_ab, tol_bp=20, min_ratio=1.0
        )
        assert fold == {1: 0}


class TestNonFolds:
    def test_exon_skip_not_folded(self):
        # shadow's single intron (100..500) spans the parent's two introns -> NOT
        # an exact suffix of the parent chain -> never folded.
        parent = _cand([(100, 200), (400, 500)], start=50, end=700)
        shadow = _cand([(100, 500)], start=50, end=700)
        em_ab = [10.0, 3.0]
        assert containment_shadow_drops(
            [parent, shadow], em_ab, tol_bp=20, min_ratio=1.0
        ) == {}

    def test_alt_3p_end_not_folded(self):
        # exact suffix intron, but 3' termini differ by > tol -> not a truncation
        # of THIS parent (different 3' end = alt-TES, keep).
        parent = _cand([(100, 200), (400, 500)], start=50, end=700)
        shadow = _cand([(400, 500)], start=300, end=900)  # 3' end 200bp off
        em_ab = [10.0, 3.0]
        assert containment_shadow_drops(
            [parent, shadow], em_ab, tol_bp=20, min_ratio=1.0
        ) == {}

    def test_unstranded_never_folded(self):
        parent = _cand([(100, 200), (400, 500)], start=50, end=700, strand=".")
        shadow = _cand([(400, 500)], start=300, end=700, strand=".")
        em_ab = [10.0, 3.0]
        assert containment_shadow_drops(
            [parent, shadow], em_ab, tol_bp=20, min_ratio=1.0
        ) == {}

    def test_higher_abundance_fragment_not_folded(self):
        # shadow more abundant than parent -> not the minor member -> keep.
        parent = _cand([(100, 200), (400, 500)], start=50, end=700)
        shadow = _cand([(400, 500)], start=300, end=700)
        em_ab = [3.0, 10.0]  # shadow >> parent
        assert containment_shadow_drops(
            [parent, shadow], em_ab, tol_bp=20, min_ratio=1.0
        ) == {}


class TestSourceRules:
    def test_gtf_parent_absorbs_novel_child(self):
        parent = _cand([(100, 200), (400, 500)], start=50, end=700,
                       source="gtf", cid="ENST1")
        shadow = _cand([(400, 500)], start=300, end=700, source="novel")
        em_ab = [10.0, 3.0]
        assert containment_shadow_drops(
            [parent, shadow], em_ab, tol_bp=20, min_ratio=1.0
        ) == {1: 0}

    def test_gtf_child_never_folded(self):
        # the truncation is a gtf candidate -> never a foldable shadow.
        parent = _cand([(100, 200), (400, 500)], start=50, end=700, source="novel")
        shadow = _cand([(400, 500)], start=300, end=700, source="gtf", cid="ENST2")
        em_ab = [10.0, 3.0]
        assert containment_shadow_drops(
            [parent, shadow], em_ab, tol_bp=20, min_ratio=1.0
        ) == {}


class TestChainingAndExclude:
    def test_chained_resolves_to_terminal_parent(self):
        # c (3 introns) > b (2, suffix of c) > a (1, suffix of b). a and b both
        # fold to the TERMINAL longest parent c (col 0).
        c = _cand([(100, 200), (400, 500), (700, 800)], start=50, end=1000)
        b = _cand([(400, 500), (700, 800)], start=300, end=1000)
        a = _cand([(700, 800)], start=600, end=1000)
        em_ab = [10.0, 5.0, 2.0]  # c, b, a
        fold = containment_shadow_drops([c, b, a], em_ab, tol_bp=20, min_ratio=1.0)
        assert fold == {1: 0, 2: 0}

    def test_excluded_parent_not_targeted(self):
        # the only possible parent (col 0) is excluded -> no fold.
        parent = _cand([(100, 200), (400, 500)], start=50, end=700)
        shadow = _cand([(400, 500)], start=300, end=700)
        em_ab = [10.0, 3.0]
        assert containment_shadow_drops(
            [parent, shadow], em_ab, tol_bp=20, min_ratio=1.0, exclude={0}
        ) == {}

    def test_empty_when_nothing_contained(self):
        a = _cand([(100, 200), (400, 500)], start=50, end=700)
        b = _cand([(1000, 1100), (1300, 1400)], start=950, end=1600)
        assert containment_shadow_drops(
            [a, b], [10.0, 8.0], tol_bp=20, min_ratio=1.0
        ) == {}


class TestDeterminism:
    def test_parent_choice_stable_across_random_ids(self):
        # Two equal-abundance parents (same structure) + a shadow. The structural
        # tie-break must pick the SAME parent regardless of the random uuid ids.
        results = set()
        for _ in range(8):
            p1 = _cand([(100, 200), (400, 500)], start=50, end=700)
            p2 = _cand([(100, 200), (400, 500)], start=50, end=700)
            shadow = _cand([(400, 500)], start=300, end=700)
            em_ab = [5.0, 5.0, 2.0]  # p1, p2 tied
            fold = containment_shadow_drops(
                [p1, p2, shadow], em_ab, tol_bp=20, min_ratio=1.0
            )
            results.add(fold.get(2))
        # identical structural key for p1/p2 -> tie resolves to the same column
        # every run (no dependence on the random candidate_id).
        assert results == {0}


class TestFoldAbundanceSurvivesAggregate:
    def test_single_interval_fold_preserves_abundance(self):
        # Simulate the runner fold on a single interval, then aggregate: the
        # parent's reported abundance must equal parent + shadow soft mass, and
        # num_assigned_reads the read-id union.
        parent = QuantResult(
            candidate_id="P", abundance=10.0, confidence=1.0,
            num_assigned_reads=4, source="novel", chrom="chr1", strand="+",
            start=50, end=700, exons=((50, 100), (200, 400), (500, 700)),
            assigned_read_ids=("r1", "r2", "r3", "r4"),
        )
        shadow_reads = ("r5", "r6")
        # runner fold: union read ids, count, add soft mass.
        union = tuple(sorted(set(parent.assigned_read_ids) | set(shadow_reads)))
        parent.assigned_read_ids = union
        parent.num_assigned_reads = len(union)
        parent.abundance += 3.0  # shadow soft mass
        agg = aggregate_across_intervals([[parent]])
        out = agg["P"]
        assert out.num_assigned_reads == 6
        assert abs(out.abundance - 13.0) < 1e-6
