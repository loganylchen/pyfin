"""Unit tests for fin.fusion.candidates.build_fusion_candidates (US-005)."""

from __future__ import annotations


from fin.fusion.candidates import Breakpoint, _revcomp, build_fusion_candidates
from fin.candidates.dataclasses import IntronChain


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _make_bp(
    chromA="chr1", posA=5000, strandA="+",
    chromB="chr2", posB=3000, strandB="+",
    support_count=1, read_ids=None,
):
    return Breakpoint(
        chromA=chromA, posA=posA, strandA=strandA,
        chromB=chromB, posB=posB, strandB=strandB,
        support_count=support_count,
        supporting_read_ids=set(read_ids or []),
    )


# ---------------------------------------------------------------------------
# _revcomp helper
# ---------------------------------------------------------------------------

def test_revcomp_simple():
    assert _revcomp("ACGT") == "ACGT"


def test_revcomp_asymmetric():
    assert _revcomp("AAAA") == "TTTT"
    assert _revcomp("CCCC") == "GGGG"


def test_revcomp_n():
    assert _revcomp("ANGT") == "ACNT"


# ---------------------------------------------------------------------------
# AC-04 basic test
# ---------------------------------------------------------------------------

def test_build_fusion_candidates_basic():
    genome = {"chr1": "A" * 10000, "chr2": "T" * 10000}
    bp = _make_bp(chromA="chr1", posA=5000, strandA="+",
                  chromB="chr2", posB=3000, strandB="+",
                  support_count=3)

    results = build_fusion_candidates([bp], genome_fasta=genome, flank_bp=500)

    assert len(results) == 1
    tc = results[0]

    assert tc.source == "fusion"
    assert tc.candidate_id.startswith("fusion_")
    assert tc.chrom == "chr1::chr2"
    assert tc.intron_chain.introns == ()
    assert tc.fusion_junction == (5000, 3000)
    assert tc.sequence == "A" * 500 + "T" * 500
    assert len(tc.sequence) == 1000
    assert tc.breakpoint_left == ("chr1", 5000, "+")
    assert tc.breakpoint_right == ("chr2", 3000, "+")


# ---------------------------------------------------------------------------
# Reverse-complement on strandA == '-'
# ---------------------------------------------------------------------------

def test_build_fusion_candidates_revcomp():
    # chr1: repeating ACGT pattern; chr2: all T
    genome = {"chr1": "ACGT" * 2500, "chr2": "T" * 10000}
    bp = _make_bp(chromA="chr1", posA=500, strandA="-",
                  chromB="chr2", posB=0, strandB="+")

    results = build_fusion_candidates([bp], genome_fasta=genome, flank_bp=500)
    assert len(results) == 1
    tc = results[0]

    # Extract the raw upstream slice and compute expected revcomp
    raw_left = "ACGT" * 2500  # genome_fasta["chr1"]
    raw_slice = raw_left[0:500]          # max(0, 500-500) : 500
    expected_left = _revcomp(raw_slice)

    assert tc.sequence.startswith(expected_left)
    assert tc.breakpoint_left == ("chr1", 500, "-")


# ---------------------------------------------------------------------------
# Multiple clusters → multiple candidates with distinct IDs
# ---------------------------------------------------------------------------

def test_build_fusion_candidates_multiple():
    genome = {"chr1": "A" * 10000, "chr2": "T" * 10000, "chr3": "G" * 10000}
    bp1 = _make_bp(chromA="chr1", posA=1000, strandA="+",
                   chromB="chr2", posB=2000, strandB="+")
    bp2 = _make_bp(chromA="chr1", posA=3000, strandA="+",
                   chromB="chr3", posB=4000, strandB="+")

    results = build_fusion_candidates([bp1, bp2], genome_fasta=genome, flank_bp=200)

    assert len(results) == 2
    ids = {tc.candidate_id for tc in results}
    assert len(ids) == 2
    assert "fusion_chr1_1000_chr2_2000" in ids
    assert "fusion_chr1_3000_chr3_4000" in ids


# ---------------------------------------------------------------------------
# Flank clipping: posA < flank_bp → left_seq shorter than flank_bp
# ---------------------------------------------------------------------------

def test_build_fusion_candidates_flank_clipping():
    genome = {"chr1": "A" * 10000, "chr2": "T" * 10000}
    # posA=100 with flank_bp=500 → only 100 bases available upstream
    bp = _make_bp(chromA="chr1", posA=100, strandA="+",
                  chromB="chr2", posB=5000, strandB="+")

    results = build_fusion_candidates([bp], genome_fasta=genome, flank_bp=500)
    assert len(results) == 1
    tc = results[0]

    # Left side clipped to 100 bases, right side full 500 bases
    assert tc.sequence == "A" * 100 + "T" * 500
    assert len(tc.sequence) == 600


# ---------------------------------------------------------------------------
# Empty intron chain: must be () not a tuple with zero-width entries
# ---------------------------------------------------------------------------

def test_build_fusion_candidate_empty_intron_chain():
    genome = {"chr1": "A" * 10000, "chr2": "T" * 10000}
    bp = _make_bp()

    results = build_fusion_candidates([bp], genome_fasta=genome, flank_bp=100)
    assert len(results) == 1
    tc = results[0]

    assert isinstance(tc.intron_chain, IntronChain)
    assert tc.intron_chain.introns == ()
    assert len(tc.intron_chain.introns) == 0


# ---------------------------------------------------------------------------
# supporting_read_ids propagated correctly
# ---------------------------------------------------------------------------

def test_build_fusion_candidates_supporting_reads():
    genome = {"chr1": "A" * 10000, "chr2": "T" * 10000}
    bp = _make_bp(read_ids=["read_a", "read_b", "read_c"])

    results = build_fusion_candidates([bp], genome_fasta=genome, flank_bp=100)
    tc = results[0]

    assert tc.supporting_read_ids == {"read_a", "read_b", "read_c"}


# ---------------------------------------------------------------------------
# Empty cluster list → empty result
# ---------------------------------------------------------------------------

def test_build_fusion_candidates_empty_input():
    results = build_fusion_candidates([], genome_fasta={}, flank_bp=500)
    assert results == []
