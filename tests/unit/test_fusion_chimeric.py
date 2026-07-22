"""Unit tests for Stage F1 fusion chimeric-read collection (fin/fusion/chimeric.py).

Covers:
  - soft-clip segment extraction (front/back) from read dicts;
  - internal-gap computation between the two arms;
  - the adapter-bridged chimera guard (drop when internal gap >= threshold);
  - arm re-alignment wiring via an injected fake aligner;
  - graceful no-op when the genome aligner is absent.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple

from fin.fusion.chimeric import (
    _clip_segments,
    _internal_gap,
    collect_chimeric_reads,
)


# --------------------------------------------------------------------------
# Fake mappy hit + aligner
# --------------------------------------------------------------------------
@dataclass
class _FakeHit:
    ctg: str
    r_st: int
    r_en: int
    strand: int
    q_st: int
    q_en: int
    mlen: int
    cigar: Tuple[Tuple[int, int], ...] = ()


class _FakeAligner:
    """Returns a fixed list of hits for any segment (ignores the sequence)."""

    def __init__(self, hits: List[_FakeHit]):
        self._hits = hits

    def map(self, seq):  # noqa: A003 - mimic mappy API
        return list(self._hits)


def _read(
    *,
    qname="r1",
    chrom="chr1",
    ref_start=1000,
    ref_end=1100,
    is_reverse=False,
    front_clip=0,
    back_clip=0,
    aln_len=100,
):
    """Build a primary-alignment read dict with optional soft-clips.

    query_sequence length = front_clip + aln_len + back_clip; the aligned span is
    [front_clip, front_clip + aln_len).
    """
    total = front_clip + aln_len + back_clip
    cigar = []
    if front_clip:
        cigar.append((4, front_clip))
    cigar.append((0, aln_len))
    if back_clip:
        cigar.append((4, back_clip))
    return {
        "query_name": qname,
        "reference_name": chrom,
        "reference_start": ref_start,
        "reference_end": ref_end,
        "cigartuples": cigar,
        "query_sequence": "A" * total,
        "query_alignment_start": front_clip,
        "query_alignment_end": front_clip + aln_len,
        "is_reverse": is_reverse,
        "is_supplementary": False,
    }


# --------------------------------------------------------------------------
# _clip_segments
# --------------------------------------------------------------------------
def test_clip_segments_back():
    r = _read(front_clip=0, aln_len=100, back_clip=80)
    segs = _clip_segments(r)
    assert segs == [("back", 100, 180)]


def test_clip_segments_front():
    r = _read(front_clip=70, aln_len=100, back_clip=0)
    segs = _clip_segments(r)
    assert segs == [("front", 0, 70)]


def test_clip_segments_none_when_no_softclip():
    r = _read(front_clip=0, aln_len=100, back_clip=0)
    assert _clip_segments(r) == []


# --------------------------------------------------------------------------
# _internal_gap
# --------------------------------------------------------------------------
def test_internal_gap_back_clean_junction():
    # primary ends at 100; partner segment starts at read 100, hit q_st=0 -> gap 0
    assert _internal_gap("back", 0, 100, seg_q_start=100, hit_q_st=0, hit_q_en=80) == 0


def test_internal_gap_back_with_adapter():
    # partner alignment starts 40 bp into the clip -> 40 bp internal unmapped gap
    assert _internal_gap("back", 0, 100, seg_q_start=100, hit_q_st=40, hit_q_en=80) == 40


def test_internal_gap_front_clean_junction():
    # front clip [0,70); partner ends at q_en=70 -> abuts primary start 70 -> gap 0
    assert _internal_gap("front", 70, 170, seg_q_start=0, hit_q_st=0, hit_q_en=70) == 0


def test_internal_gap_front_with_adapter():
    # partner ends at 30, primary starts at 70 -> 40 bp gap
    assert _internal_gap("front", 70, 170, seg_q_start=0, hit_q_st=0, hit_q_en=30) == 40


# --------------------------------------------------------------------------
# collect_chimeric_reads
# --------------------------------------------------------------------------
def test_collect_emits_clean_back_clip_fusion():
    r = _read(qname="f1", chrom="chr1", ref_start=1000, ref_end=1100, back_clip=260)
    aligner = _FakeAligner([_FakeHit("chr2", 5000, 5080, 1, 0, 80, mlen=80)])
    out = collect_chimeric_reads([r], aligner, max_internal_gap_bp=30)
    assert len(out) == 1
    cr = out[0]
    assert cr.query_name == "f1"
    assert cr.arm_a.chrom == "chr1" and cr.arm_b.chrom == "chr2"
    assert cr.breakpoint_a == ("chr1", 1100, "+")   # primary reference_end
    assert cr.breakpoint_b == ("chr2", 5000, "+")   # partner r_st
    assert cr.internal_gap == 0


def test_collect_drops_adapter_chimera():
    r = _read(qname="art", ref_end=1100, back_clip=260)
    # partner aligns only after a 50 bp unmapped prefix in the clip -> gap 50
    aligner = _FakeAligner([_FakeHit("chr2", 5000, 5070, 1, 50, 120, mlen=70)])
    out = collect_chimeric_reads([r], aligner, max_internal_gap_bp=30)
    assert out == []


def test_collect_keeps_when_gap_below_threshold():
    r = _read(qname="ok", ref_end=1100, back_clip=260)
    aligner = _FakeAligner([_FakeHit("chr2", 5000, 5090, 1, 20, 110, mlen=90)])
    out = collect_chimeric_reads([r], aligner, max_internal_gap_bp=30)
    assert len(out) == 1 and out[0].internal_gap == 20


def test_collect_gap_guard_disabled_when_zero():
    r = _read(qname="big", ref_end=1100, back_clip=260)
    aligner = _FakeAligner([_FakeHit("chr2", 5000, 5100, 1, 100, 200, mlen=100)])
    out = collect_chimeric_reads([r], aligner, max_internal_gap_bp=0)
    assert len(out) == 1 and out[0].internal_gap == 100


def test_collect_skips_non_fusion_reads():
    r = _read(qname="plain", front_clip=0, aln_len=100, back_clip=0)  # no soft-clip
    aligner = _FakeAligner([_FakeHit("chr2", 5000, 5080, 1, 0, 80, mlen=80)])
    assert collect_chimeric_reads([r], aligner) == []


def test_collect_skips_when_segment_unmapped():
    r = _read(qname="nohit", back_clip=260)
    aligner = _FakeAligner([])  # partner segment maps nowhere
    assert collect_chimeric_reads([r], aligner) == []


def test_collect_none_aligner_returns_empty():
    r = _read(back_clip=80)
    assert collect_chimeric_reads([r], None) == []


def test_collect_minus_strand_partner_back_clip():
    # back clip + '-' partner: junction edge is the partner's HIGH genomic end.
    r = _read(qname="rev", ref_end=1100, back_clip=260)
    aligner = _FakeAligner([_FakeHit("chr3", 200, 280, -1, 0, 80, mlen=80)])
    out = collect_chimeric_reads([r], aligner, max_internal_gap_bp=30)
    assert len(out) == 1
    assert out[0].arm_b.strand == "-"
    assert out[0].breakpoint_b == ("chr3", 280, "-")


def test_collect_minus_strand_primary_back_clip():
    # pysam normalizes the primary cigar to genome-forward regardless of
    # is_reverse, so a back clip's junction is ref_end even on a '-' read.
    r = _read(qname="rp", ref_start=1000, ref_end=1100, is_reverse=True, back_clip=260)
    aligner = _FakeAligner([_FakeHit("chr2", 5000, 5080, 1, 0, 80, mlen=80)])
    out = collect_chimeric_reads([r], aligner, max_internal_gap_bp=30)
    assert len(out) == 1
    assert out[0].breakpoint_a == ("chr1", 1100, "-")   # ref_end (no strand flip)
    assert out[0].breakpoint_b == ("chr2", 5000, "+")   # '+' partner -> r_st


def test_collect_front_clip_plus_strands():
    # front clip on '+' primary with '+' partner: arm A junction = ref_start,
    # arm B junction = partner r_en.
    r = _read(qname="fc", ref_start=1000, ref_end=1100, front_clip=260, aln_len=100)
    # partner maps the full 260bp front clip -> r_en 5070 (junction edge), gap 0
    aligner = _FakeAligner([_FakeHit("chr2", 4810, 5070, 1, 0, 260, mlen=260)])
    out = collect_chimeric_reads([r], aligner, max_internal_gap_bp=30)
    assert len(out) == 1
    assert out[0].breakpoint_a == ("chr1", 1000, "+")   # front -> ref_start on '+'
    assert out[0].breakpoint_b == ("chr2", 5070, "+")   # front -> r_en on '+'
