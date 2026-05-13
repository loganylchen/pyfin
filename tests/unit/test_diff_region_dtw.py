"""Tests for fin/scoring/diff_region_dtw.py.

Covers:
  - extract_diff_regions: basic diff detection, single candidate, merged regions
  - genomic_region_to_cdna: plus-strand, minus-strand, intronic overlap
  - cdna_region_to_signal_range: event lookup
  - compute_diff_region_m4: zero matrix on no diff regions (AC5), NaN handling
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval
from fin.scoring.diff_region_dtw import (
    cdna_region_to_signal_range,
    compute_diff_region_m4,
    extract_diff_regions,
    genomic_region_to_cdna,
)
from fin.scoring.eventalign_parser import ReadCandidateScore


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_candidate(
    cid: str,
    chrom: str,
    start: int,
    end: int,
    strand: str = "+",
    introns: Tuple[Tuple[int, int], ...] = (),
    seq: str = "",
) -> TranscriptCandidate:
    interval = GenomicInterval(chrom=chrom, start=start, end=end, strand=strand)
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=introns),
        three_prime_pos=end,
        sequence=seq or ("A" * (end - start)),
        source="novel",
        supporting_read_ids=frozenset(),
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
    )


def _make_score(
    read_name: str,
    cid: str,
    events: List[Tuple[int, int, int]],
) -> ReadCandidateScore:
    s = ReadCandidateScore(read_name=read_name, candidate_id=cid)
    s.events = events
    return s


# ---------------------------------------------------------------------------
# extract_diff_regions
# ---------------------------------------------------------------------------

class TestExtractDiffRegions:
    def test_fewer_than_two_candidates_returns_empty(self):
        c = _make_candidate("c1", "chr1", 100, 200)
        assert extract_diff_regions([c]) == []
        assert extract_diff_regions([]) == []

    def test_identical_single_exon_no_diff(self):
        """Two single-exon candidates with same span → no diff bases."""
        c1 = _make_candidate("c1", "chr1", 100, 300)
        c2 = _make_candidate("c2", "chr1", 100, 300)
        assert extract_diff_regions([c1, c2]) == []

    def test_one_has_intron_other_does_not(self):
        """c1: exon 100-300 (no intron); c2: exon 100-150, intron 150-250, exon 250-300.
        Bases 150-250 are exon in c1 and intron in c2 → diff region [150, 250)."""
        c1 = _make_candidate("c1", "chr1", 100, 300)
        c2 = _make_candidate("c2", "chr1", 100, 300, introns=((150, 250),))
        regions = extract_diff_regions([c1, c2])
        assert len(regions) == 1
        assert regions[0] == (150, 250)

    def test_multiple_diff_regions(self):
        """Two diff introns → two separate diff regions."""
        c1 = _make_candidate("c1", "chr1", 0, 600)
        # c2 has two introns: [100,200) and [400,500)
        c2 = _make_candidate("c2", "chr1", 0, 600, introns=((100, 200), (400, 500)))
        regions = extract_diff_regions([c1, c2])
        assert len(regions) == 2
        assert regions[0] == (100, 200)
        assert regions[1] == (400, 500)

    def test_no_diff_if_all_candidates_have_same_intron(self):
        """Both candidates have the same intron → no exon/intron split."""
        c1 = _make_candidate("c1", "chr1", 0, 500, introns=((100, 200),))
        c2 = _make_candidate("c2", "chr1", 0, 500, introns=((100, 200),))
        assert extract_diff_regions([c1, c2]) == []


# ---------------------------------------------------------------------------
# genomic_region_to_cdna
# ---------------------------------------------------------------------------

class TestGenomicRegionToCdna:
    def test_single_exon_full_overlap(self):
        """Single exon 100-300, region 100-300 → cDNA 0-200."""
        c = _make_candidate("c1", "chr1", 100, 300)
        result = genomic_region_to_cdna(c, (100, 300))
        assert result == (0, 200)

    def test_single_exon_partial_overlap(self):
        """Single exon 100-300, region 150-250 → cDNA 50-150."""
        c = _make_candidate("c1", "chr1", 100, 300)
        result = genomic_region_to_cdna(c, (150, 250))
        assert result == (50, 150)

    def test_no_overlap_returns_none(self):
        c = _make_candidate("c1", "chr1", 100, 300)
        assert genomic_region_to_cdna(c, (400, 500)) is None

    def test_intronic_region_returns_none(self):
        """Region entirely within intron → no exonic bases → None."""
        c = _make_candidate("c1", "chr1", 0, 500, introns=((100, 400),))
        assert genomic_region_to_cdna(c, (150, 350)) is None

    def test_two_exon_plus_strand(self):
        """Exons: 0-100, 200-300; region 0-100 → cDNA 0-100."""
        c = _make_candidate("c1", "chr1", 0, 300, introns=((100, 200),))
        result = genomic_region_to_cdna(c, (0, 100))
        assert result == (0, 100)

    def test_two_exon_plus_strand_second_exon(self):
        """Exons: 0-100, 200-300; region 200-300 → cDNA 100-200."""
        c = _make_candidate("c1", "chr1", 0, 300, introns=((100, 200),))
        result = genomic_region_to_cdna(c, (200, 300))
        assert result == (100, 200)

    def test_minus_strand_single_exon(self):
        """Minus strand single exon 0-100; region 0-100.
        cDNA is reverse-complemented: spliced_len=100, c_start_rc=100-100=0, c_end_rc=100-0=100."""
        c = _make_candidate("c1", "chr1", 0, 100, strand="-")
        result = genomic_region_to_cdna(c, (0, 100))
        assert result == (0, 100)

    def test_minus_strand_partial(self):
        """Minus strand single exon 0-100; region 60-100.
        Genomic c_start=60, c_end=100. RC: c_start_rc=100-100=0, c_end_rc=100-60=40."""
        c = _make_candidate("c1", "chr1", 0, 100, strand="-")
        result = genomic_region_to_cdna(c, (60, 100))
        assert result == (0, 40)


# ---------------------------------------------------------------------------
# cdna_region_to_signal_range
# ---------------------------------------------------------------------------

class TestCdnaRegionToSignalRange:
    def test_no_events_returns_none(self):
        s = ReadCandidateScore(read_name="r1", candidate_id="c1")
        assert cdna_region_to_signal_range(s, (0, 100)) is None

    def test_events_in_range(self):
        s = _make_score("r1", "c1", [(10, 100, 200), (20, 200, 300), (50, 300, 400)])
        result = cdna_region_to_signal_range(s, (0, 60))
        assert result == (100, 400)

    def test_events_outside_range_returns_none(self):
        s = _make_score("r1", "c1", [(100, 0, 50)])
        assert cdna_region_to_signal_range(s, (0, 50)) is None

    def test_boundary_exclusive_end(self):
        """pos == c_end is excluded (half-open interval)."""
        s = _make_score("r1", "c1", [(50, 0, 100)])
        assert cdna_region_to_signal_range(s, (0, 50)) is None

    def test_single_event_in_range(self):
        s = _make_score("r1", "c1", [(5, 10, 20)])
        result = cdna_region_to_signal_range(s, (0, 10))
        assert result == (10, 20)

    def test_zero_width_signal_range_returns_none(self):
        """sig_lo == sig_hi → None."""
        s = _make_score("r1", "c1", [(5, 50, 50)])
        assert cdna_region_to_signal_range(s, (0, 10)) is None


# ---------------------------------------------------------------------------
# compute_diff_region_m4
# ---------------------------------------------------------------------------

class TestComputeDiffRegionM4:
    def test_empty_read_ids_returns_empty_matrix(self):
        m4 = compute_diff_region_m4(
            read_ids=[],
            candidates=[],
            scores_by_pair={},
            signal_reader=None,
            interval_start=0,
            interval_end=1000,
        )
        assert m4.shape == (0, 0)

    def test_no_diff_regions_returns_zero_matrix(self):
        """AC5: when no diff regions, must return zeros (not NaN)."""
        # Two identical single-exon candidates → no diff regions.
        c1 = _make_candidate("c1", "chr1", 0, 500)
        c2 = _make_candidate("c2", "chr1", 0, 500)
        read_ids = ["r1", "r2"]
        m4 = compute_diff_region_m4(
            read_ids=read_ids,
            candidates=[c1, c2],
            scores_by_pair={},
            signal_reader=None,
            interval_start=0,
            interval_end=500,
        )
        assert m4.shape == (2, 2)
        assert np.all(m4 == 0.0), f"Expected all zeros, got:\n{m4}"

    def test_diagonal_is_zero(self):
        """Self-distance must always be 0."""
        c1 = _make_candidate("c1", "chr1", 0, 500)
        c2 = _make_candidate("c2", "chr1", 0, 500, introns=((100, 400),))
        # No signal reader → all segments will be empty → off-diagonal NaN.
        m4 = compute_diff_region_m4(
            read_ids=["r1", "r2"],
            candidates=[c1, c2],
            scores_by_pair={},
            signal_reader=None,
            interval_start=0,
            interval_end=500,
        )
        assert m4[0, 0] == 0.0
        assert m4[1, 1] == 0.0

    def test_single_read_returns_1x1_zero(self):
        c1 = _make_candidate("c1", "chr1", 0, 200)
        c2 = _make_candidate("c2", "chr1", 0, 200, introns=((50, 150),))
        m4 = compute_diff_region_m4(
            read_ids=["r1"],
            candidates=[c1, c2],
            scores_by_pair={},
            signal_reader=None,
            interval_start=0,
            interval_end=200,
        )
        assert m4.shape == (1, 1)
        assert m4[0, 0] == 0.0
