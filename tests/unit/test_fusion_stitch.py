"""Unit tests for Stage F3 fusion stitching (fin/fusion/stitch.py).

Covers:
  - stitched sequence = spliced(arm A) ++ spliced(arm B) with per-arm strand;
  - distinct splice combinations -> distinct candidates;
  - read x read support = arm read-set intersection, gated by min_support;
  - annotation-containing combos are kept regardless of support;
  - dedup by (chainA, chainB) unions supporting reads;
  - breakpoint / fusion_junction fields populated.
"""
from __future__ import annotations

from fin.candidates.dataclasses import IntronChain
from fin.fusion.arm_assembly import ArmVariant, FusionPairCluster
from fin.fusion.chimeric import ArmAlignment, ChimericRead
from fin.fusion.stitch import build_fusion_candidates_v2, stitch_cluster


# --------------------------------------------------------------------------
# Builders
# --------------------------------------------------------------------------
def _variant(chrom, start, end, introns=(), strand="+", source="read", reads=()):
    return ArmVariant(
        chrom=chrom, strand=strand, start=start, end=end,
        intron_chain=IntronChain(tuple(introns)), source=source,
        supporting_read_ids=set(reads),
    )


def _read(qname):
    arm = ArmAlignment("chr1", 0, 1, "+", 0, 1, ())
    return ChimericRead(qname, arm, arm, ("chr1", 1, "+"), ("chr2", 5000, "+"), 0)


def _cluster(a_vars, b_vars, read_names=("r1", "r2"),
             bp_a=("chr1", 1500, "+"), bp_b=("chr2", 5000, "+")):
    return FusionPairCluster(
        signature=(bp_a[0], bp_a[2], bp_b[0], bp_b[2]),
        breakpoint_a=bp_a,
        breakpoint_b=bp_b,
        reads=[_read(n) for n in read_names],
        arm_a_variants=list(a_vars),
        arm_b_variants=list(b_vars),
    )


# A genome: chrom -> sequence. Use distinct letters so we can verify concat.
_GENOME = {
    "chr1": "A" * 2000,
    "chr2": "C" * 7000,
}


# --------------------------------------------------------------------------
# Sequence stitching
# --------------------------------------------------------------------------
def test_stitch_concatenates_arm_sequences():
    a = _variant("chr1", 1000, 1100, reads=("r1", "r2"))   # 100 A's
    b = _variant("chr2", 5000, 5150, reads=("r1", "r2"))   # 150 C's
    cands = stitch_cluster(_cluster([a], [b]), _GENOME, min_support=2)
    assert len(cands) == 1
    seq = cands[0].sequence
    assert seq == "A" * 100 + "C" * 150


def test_stitch_spliced_arm_drops_intron():
    # arm A [1000,1300) with intron (1100,1250) -> exon lengths 100 + 50 = 150
    a = _variant("chr1", 1000, 1300, introns=((1100, 1250),), reads=("r1", "r2"))
    b = _variant("chr2", 5000, 5100, reads=("r1", "r2"))
    cands = stitch_cluster(_cluster([a], [b]), _GENOME, min_support=2)
    assert len(cands[0].sequence) == 150 + 100


# --------------------------------------------------------------------------
# Variant combinatorics
# --------------------------------------------------------------------------
def test_distinct_splice_combos_distinct_candidates():
    a1 = _variant("chr1", 1000, 1100, reads=("r1", "r2"))
    a2 = _variant("chr1", 1000, 1200, introns=((1050, 1150),), reads=("r1", "r2"))
    b = _variant("chr2", 5000, 5100, reads=("r1", "r2"))
    cands = stitch_cluster(_cluster([a1, a2], [b]), _GENOME, min_support=2)
    assert len(cands) == 2
    ids = {c.candidate_id for c in cands}
    assert len(ids) == 2  # distinct, deterministic ids


# --------------------------------------------------------------------------
# Support model
# --------------------------------------------------------------------------
def test_read_combo_support_is_intersection():
    a = _variant("chr1", 1000, 1100, reads=("r1", "r2", "r3"))
    b = _variant("chr2", 5000, 5100, reads=("r2", "r3", "r4"))
    cands = stitch_cluster(_cluster([a], [b]), _GENOME, min_support=2)
    assert cands[0].supporting_read_ids == {"r2", "r3"}


def test_read_combo_dropped_below_min_support():
    a = _variant("chr1", 1000, 1100, reads=("r1", "r2"))
    b = _variant("chr2", 5000, 5100, reads=("r2", "r9"))   # intersection {r2} = 1
    cands = stitch_cluster(_cluster([a], [b]), _GENOME, min_support=2)
    assert cands == []


def test_annotation_combo_kept_regardless_of_support():
    # arm A annotation (no reads), arm B read-derived
    a = _variant("chr1", 1000, 1100, source="annotation", reads=())
    b = _variant("chr2", 5000, 5100, source="read", reads=("r1",))
    cands = stitch_cluster(_cluster([a], [b]), _GENOME, min_support=5)
    assert len(cands) == 1
    assert cands[0].supporting_read_ids == {"r1"}   # the read-derived arm's reads


def test_both_annotation_uses_cluster_reads():
    a = _variant("chr1", 1000, 1100, source="annotation", reads=())
    b = _variant("chr2", 5000, 5100, source="annotation", reads=())
    cands = stitch_cluster(_cluster([a], [b], read_names=("x", "y", "z")),
                           _GENOME, min_support=2)
    assert len(cands) == 1
    assert cands[0].supporting_read_ids == {"x", "y", "z"}


# --------------------------------------------------------------------------
# Dedup
# --------------------------------------------------------------------------
def test_dedup_unions_reads_for_same_chains():
    # Two arm-A variants with the SAME chain but different read sets collapse.
    a1 = _variant("chr1", 1000, 1100, introns=(), reads=("r1",))
    a2 = _variant("chr1", 1000, 1100, introns=(), reads=("r2",))
    b = _variant("chr2", 5000, 5100, reads=("r1", "r2"))
    cands = stitch_cluster(_cluster([a1, a2], [b]), _GENOME, min_support=1)
    assert len(cands) == 1
    assert cands[0].supporting_read_ids == {"r1", "r2"}


# --------------------------------------------------------------------------
# Fields + entry point
# --------------------------------------------------------------------------
def test_breakpoint_fields_populated():
    a = _variant("chr1", 1000, 1100, reads=("r1", "r2"))
    b = _variant("chr2", 5000, 5100, reads=("r1", "r2"))
    cl = _cluster([a], [b], bp_a=("chr1", 1500, "+"), bp_b=("chr2", 5000, "-"))
    c = stitch_cluster(cl, _GENOME, min_support=2)[0]
    assert c.source == "fusion"
    assert c.chrom == "chr1::chr2"
    assert c.fusion_junction == (1500, 5000)
    assert c.breakpoint_left == ("chr1", 1500, "+")
    assert c.breakpoint_right == ("chr2", 5000, "-")
    assert c.start == 1500 and c.end == 5000


def test_build_entry_point_over_multiple_clusters():
    a = _variant("chr1", 1000, 1100, reads=("r1", "r2"))
    b = _variant("chr2", 5000, 5100, reads=("r1", "r2"))
    c1 = _cluster([a], [b], bp_a=("chr1", 1500, "+"), bp_b=("chr2", 5000, "+"))
    c2 = _cluster([a], [b], bp_a=("chr1", 9000, "+"), bp_b=("chr2", 8000, "+"))
    out = build_fusion_candidates_v2([c1, c2], _GENOME, min_support=2)
    assert len(out) == 2


def test_empty_arm_pool_yields_nothing():
    a = _variant("chr1", 1000, 1100, reads=("r1", "r2"))
    assert stitch_cluster(_cluster([a], []), _GENOME, min_support=2) == []


def test_missing_arm_chrom_drops_candidate():
    # arm B is on a contig absent from the genome -> seq_b == "" -> no half-fusion
    a = _variant("chr1", 1000, 1100, reads=("r1", "r2"))
    b = _variant("chrUNKNOWN", 5000, 5100, reads=("r1", "r2"))
    cl = _cluster([a], [b], bp_b=("chrUNKNOWN", 5000, "+"))
    assert stitch_cluster(cl, _GENOME, min_support=2) == []
