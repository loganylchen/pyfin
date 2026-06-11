"""Unit tests for Stage F2 fusion per-arm splice inference (arm_assembly.py).

Covers:
  - cluster_chimeric_reads: signature bucketing, breakpoint-proximity merge,
    consensus breakpoint, determinism;
  - read-derived chain variants from CIGAR introns + min_support filter;
  - annotation chains as parallel additive candidates (strand-filtered, exempt
    from min_support);
  - infer_arm_variants / assemble_fusion_arms wiring.
"""
from __future__ import annotations

from fin.candidates.dataclasses import IntronChain
from fin.fusion.arm_assembly import (
    assemble_fusion_arms,
    cluster_chimeric_reads,
    infer_arm_variants,
)
from fin.fusion.chimeric import ArmAlignment, ChimericRead


# --------------------------------------------------------------------------
# Builders
# --------------------------------------------------------------------------
def _arm(chrom, ref_start, ref_end, strand="+", cigar=((0, 100),)):
    return ArmAlignment(
        chrom=chrom, ref_start=ref_start, ref_end=ref_end, strand=strand,
        q_start=0, q_end=100, cigartuples=tuple(cigar),
    )


def _cread(qname, *, a_chrom="chr1", a_pos=1100, a_strand="+",
           b_chrom="chr2", b_pos=5000, b_strand="+",
           arm_a=None, arm_b=None):
    return ChimericRead(
        query_name=qname,
        arm_a=arm_a or _arm(a_chrom, a_pos - 100, a_pos, a_strand),
        arm_b=arm_b or _arm(b_chrom, b_pos, b_pos + 100, b_strand),
        breakpoint_a=(a_chrom, a_pos, a_strand),
        breakpoint_b=(b_chrom, b_pos, b_strand),
        internal_gap=0,
    )


# Fake GTF transcript + reader
class _FakeTx:
    def __init__(self, tid, strand, exons):
        self.transcript_id = tid
        self.strand = strand
        self.exons = list(exons)
        self.start = min(e[0] for e in exons)
        self.end = max(e[1] for e in exons)

    def sort_features(self):
        self.exons.sort()


class _FakeGtf:
    def __init__(self, txs):
        self._txs = txs

    def get_transcripts_in_region(self, chrom, start, end):
        return [t for t in self._txs
                if any(s < end and e > start for s, e in t.exons)]


# --------------------------------------------------------------------------
# cluster_chimeric_reads
# --------------------------------------------------------------------------
def test_cluster_merges_nearby_same_signature():
    reads = [
        _cread("r1", a_pos=1100, b_pos=5000),
        _cread("r2", a_pos=1110, b_pos=5005),   # within max_dist of r1
    ]
    clusters = cluster_chimeric_reads(reads, max_dist=500)
    assert len(clusters) == 1
    c = clusters[0]
    assert len(c.reads) == 2
    assert c.breakpoint_a[1] == round((1100 + 1110) / 2)
    assert c.breakpoint_b[1] == round((5000 + 5005) / 2)


def test_cluster_splits_distant_breakpoints():
    reads = [
        _cread("r1", a_pos=1100, b_pos=5000),
        _cread("r2", a_pos=9000, b_pos=5000),   # far on arm A
    ]
    clusters = cluster_chimeric_reads(reads, max_dist=500)
    assert len(clusters) == 2


def test_cluster_splits_different_signature():
    reads = [
        _cread("r1", b_chrom="chr2", b_strand="+"),
        _cread("r2", b_chrom="chr3", b_strand="+"),  # different partner chrom
    ]
    clusters = cluster_chimeric_reads(reads, max_dist=500)
    assert len(clusters) == 2


def test_cluster_empty():
    assert cluster_chimeric_reads([], max_dist=500) == []


# --------------------------------------------------------------------------
# read-derived variants (wobble off via search_bp=0)
# --------------------------------------------------------------------------
def _spliced_arm(chrom, ref_start, intron, strand="+"):
    # one intron of given (s,e); CIGAR M (to intron start) N (intron) M (rest)
    s, e = intron
    m1 = s - ref_start
    n = e - s
    cig = ((0, m1), (3, n), (0, 100))
    return _arm(chrom, ref_start, e + 100, strand, cig)


def test_read_variant_distinct_chains_and_min_support():
    # 2 reads share intron (1200,1400); 1 read has a different intron (1200,1500)
    a1 = _spliced_arm("chr1", 1000, (1200, 1400))
    a2 = _spliced_arm("chr1", 1000, (1200, 1400))
    a3 = _spliced_arm("chr1", 1000, (1200, 1500))
    reads = [
        _cread("r1", a_pos=1500, arm_a=a1),
        _cread("r2", a_pos=1500, arm_a=a2),
        _cread("r3", a_pos=1600, arm_a=a3),
    ]
    cluster = cluster_chimeric_reads(reads, max_dist=500)[0]
    variants = infer_arm_variants(
        cluster, "a", genome_by_chrom={}, gtf_reader=None,
        motif_set=None, search_bp=0, min_support=2,
    )
    chains = {v.intron_chain for v in variants}
    # only the 2-read chain survives min_support=2
    assert IntronChain(((1200, 1400),)) in chains
    assert IntronChain(((1200, 1500),)) not in chains
    assert all(v.source == "read" for v in variants)


def test_read_variant_span_is_min_max():
    a1 = _spliced_arm("chr1", 1000, (1200, 1400))
    a2 = _spliced_arm("chr1", 900, (1200, 1400))   # starts earlier
    reads = [_cread("r1", a_pos=1500, arm_a=a1), _cread("r2", a_pos=1500, arm_a=a2)]
    cluster = cluster_chimeric_reads(reads, max_dist=500)[0]
    v = infer_arm_variants(cluster, "a", {}, None, None, 0, min_support=2)[0]
    assert v.start == 900            # min ref_start
    assert v.end == 1500             # max ref_end (1400+100)


# --------------------------------------------------------------------------
# annotation variants (parallel, additive)
# --------------------------------------------------------------------------
def test_annotation_variants_added_and_strand_filtered():
    a1 = _spliced_arm("chr1", 1000, (1200, 1400))
    reads = [_cread("r1", a_pos=1500, arm_a=a1), _cread("r2", a_pos=1500, arm_a=a1)]
    cluster = cluster_chimeric_reads(reads, max_dist=500)[0]
    gtf = _FakeGtf([
        _FakeTx("plus_tx", "+", [(1000, 1180), (1450, 1600)]),   # same strand → kept
        _FakeTx("minus_tx", "-", [(1000, 1180), (1450, 1600)]),  # other strand → dropped
    ])
    variants = infer_arm_variants(
        cluster, "a", genome_by_chrom={}, gtf_reader=gtf,
        motif_set=None, search_bp=0, min_support=2,
    )
    annot = [v for v in variants if v.source == "annotation"]
    assert len(annot) == 1
    assert annot[0].intron_chain == IntronChain(((1180, 1450),))
    # read-derived variant still present alongside annotation (parallel, not snapped)
    assert any(v.source == "read" and v.intron_chain == IntronChain(((1200, 1400),))
               for v in variants)


def test_annotation_variant_exempt_from_min_support():
    # single read (below min_support) → no read variant, but annotation still offered
    a1 = _spliced_arm("chr1", 1000, (1200, 1400))
    reads = [_cread("r1", a_pos=1500, arm_a=a1)]
    cluster = cluster_chimeric_reads(reads, max_dist=500)[0]
    gtf = _FakeGtf([_FakeTx("tx", "+", [(1000, 1180), (1450, 1600)])])
    variants = infer_arm_variants(cluster, "a", {}, gtf, None, 0, min_support=2)
    assert [v.source for v in variants] == ["annotation"]


# --------------------------------------------------------------------------
# assemble_fusion_arms end-to-end
# --------------------------------------------------------------------------
def test_assemble_populates_both_arms():
    a = _spliced_arm("chr1", 1000, (1200, 1400))
    b = _spliced_arm("chr2", 5000, (5200, 5400))
    reads = [
        _cread("r1", a_pos=1500, b_pos=5500, arm_a=a, arm_b=b),
        _cread("r2", a_pos=1500, b_pos=5500, arm_a=a, arm_b=b),
    ]
    clusters = assemble_fusion_arms(reads, genome_by_chrom={}, min_support=2, search_bp=0)
    assert len(clusters) == 1
    c = clusters[0]
    assert c.arm_a_variants and c.arm_b_variants
    assert c.arm_a_variants[0].intron_chain == IntronChain(((1200, 1400),))
    assert c.arm_b_variants[0].intron_chain == IntronChain(((5200, 5400),))
