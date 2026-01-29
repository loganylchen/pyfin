#!/usr/bin/env python3
"""
Unit tests for fin.io.interval_manager module.

Tests the genomic interval management functions:
- GenomicInterval: Data class for genomic intervals
- intervals_to_bed: Export intervals to BED format
- is_fusion_read: Check for fusion read candidates
- extract_strand_from_read: Extract strand from alignment
- Other helper functions

Run with:
    pytest tests/unit/test_interval_manager.py -v
"""

import pytest
import tempfile
import os
import sys
from pathlib import Path
from io import StringIO

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.io.interval_manager import (
    GenomicInterval,
    intervals_to_bed,
    is_fusion_read,
    extract_strand_from_read,
    extract_three_prime_pos,
)


# ============================================================================
# GenomicInterval tests
# ============================================================================

class TestGenomicInterval:
    """Tests for GenomicInterval dataclass."""

    def test_interval_initialization_basic(self):
        """Test basic initialization of GenomicInterval."""
        interval = GenomicInterval(chrom="chr1", start=100, end=200)
        assert interval.chrom == "chr1"
        assert interval.start == 100
        assert interval.end == 200

    def test_interval_initialization_with_strand(self):
        """Test initialization with strand."""
        interval = GenomicInterval(chrom="chr1", start=100, end=200, strand="+")
        assert interval.strand == "+"

    def test_interval_initialization_with_read_count(self):
        """Test initialization with read count."""
        interval = GenomicInterval(chrom="chr1", start=100, end=200, read_count=50)
        assert interval.read_count == 50

    def test_interval_initialization_with_id(self):
        """Test initialization with interval_id."""
        interval = GenomicInterval(
            chrom="chr1", start=100, end=200, interval_id="interval_1"
        )
        assert interval.interval_id == "interval_1"

    def test_interval_default_values(self):
        """Test default values are set correctly."""
        interval = GenomicInterval(chrom="chr1", start=100, end=200)
        assert interval.strand is None
        assert interval.read_count == 0
        assert interval.interval_id is None
        assert interval.attrs == []
        assert interval.three_prime_pos == []

    def test_interval_tuple_property(self):
        """Test interval_tuple property returns (chrom, start, end)."""
        interval = GenomicInterval(chrom="chr1", start=100, end=200)
        assert interval.interval_tuple == ("chr1", 100, 200)

    def test_interval_region_string_property(self):
        """Test region_string property."""
        interval = GenomicInterval(chrom="chr1", start=100, end=200)
        assert interval.region_string == "chr1:100-200"

    def test_interval_overlaps_true(self):
        """Test overlaps method returns True for overlapping intervals."""
        interval1 = GenomicInterval(chrom="chr1", start=100, end=200)
        interval2 = GenomicInterval(chrom="chr1", start=150, end=250)
        assert interval1.overlaps(interval2) is True
        assert interval2.overlaps(interval1) is True

    def test_interval_overlaps_false_no_overlap(self):
        """Test overlaps method returns False for non-overlapping intervals."""
        interval1 = GenomicInterval(chrom="chr1", start=100, end=200)
        interval2 = GenomicInterval(chrom="chr1", start=300, end=400)
        assert interval1.overlaps(interval2) is False

    def test_interval_overlaps_false_different_chrom(self):
        """Test overlaps returns False for different chromosomes."""
        interval1 = GenomicInterval(chrom="chr1", start=100, end=200)
        interval2 = GenomicInterval(chrom="chr2", start=100, end=200)
        assert interval1.overlaps(interval2) is False

    def test_interval_overlaps_adjacent(self):
        """Test adjacent intervals do not overlap."""
        interval1 = GenomicInterval(chrom="chr1", start=100, end=200)
        interval2 = GenomicInterval(chrom="chr1", start=200, end=300)
        assert interval1.overlaps(interval2) is False

    def test_interval_overlaps_contained(self):
        """Test interval contained within another overlaps."""
        interval1 = GenomicInterval(chrom="chr1", start=100, end=400)
        interval2 = GenomicInterval(chrom="chr1", start=200, end=300)
        assert interval1.overlaps(interval2) is True
        assert interval2.overlaps(interval1) is True

    def test_interval_merge(self):
        """Test merging two overlapping intervals."""
        interval1 = GenomicInterval(
            chrom="chr1", start=100, end=200, strand="+", read_count=10
        )
        interval2 = GenomicInterval(
            chrom="chr1", start=150, end=250, strand="+", read_count=15
        )
        
        merged = interval1.merge(interval2)
        
        assert merged.chrom == "chr1"
        assert merged.start == 100
        assert merged.end == 250
        assert merged.strand == "+"
        assert merged.read_count == 25

    def test_interval_merge_different_strands_raises(self):
        """Test merging intervals with different strands raises error."""
        interval1 = GenomicInterval(chrom="chr1", start=100, end=200, strand="+")
        interval2 = GenomicInterval(chrom="chr1", start=150, end=250, strand="-")
        
        with pytest.raises(ValueError, match="different strands"):
            interval1.merge(interval2)

    def test_interval_merge_preserves_attrs(self):
        """Test merge combines attrs lists."""
        interval1 = GenomicInterval(
            chrom="chr1", start=100, end=200, strand="+", attrs=["attr1"]
        )
        interval2 = GenomicInterval(
            chrom="chr1", start=150, end=250, strand="+", attrs=["attr2"]
        )
        
        merged = interval1.merge(interval2)
        assert "attr1" in merged.attrs
        assert "attr2" in merged.attrs


# ============================================================================
# intervals_to_bed tests
# ============================================================================

class TestIntervalsToBed:
    """Tests for intervals_to_bed function."""

    def test_intervals_to_bed_empty(self):
        """Test intervals_to_bed with empty list."""
        result = intervals_to_bed([])
        assert result == ""

    def test_intervals_to_bed_single_interval(self):
        """Test intervals_to_bed with single interval."""
        intervals = [GenomicInterval(chrom="chr1", start=100, end=200)]
        result = intervals_to_bed(intervals)
        
        lines = result.strip().split("\n")
        assert len(lines) == 1
        fields = lines[0].split("\t")
        assert fields[0] == "chr1"
        assert fields[1] == "100"
        assert fields[2] == "200"

    def test_intervals_to_bed_with_strand(self):
        """Test intervals_to_bed includes strand."""
        intervals = [GenomicInterval(chrom="chr1", start=100, end=200, strand="+")]
        result = intervals_to_bed(intervals)
        
        fields = result.strip().split("\t")
        assert fields[5] == "+"

    def test_intervals_to_bed_with_read_count(self):
        """Test intervals_to_bed includes read_count as score."""
        intervals = [
            GenomicInterval(chrom="chr1", start=100, end=200, read_count=50)
        ]
        result = intervals_to_bed(intervals)
        
        fields = result.strip().split("\t")
        assert fields[4] == "50"

    def test_intervals_to_bed_with_interval_id(self):
        """Test intervals_to_bed includes interval_id as name."""
        intervals = [
            GenomicInterval(
                chrom="chr1", start=100, end=200, interval_id="my_interval"
            )
        ]
        result = intervals_to_bed(intervals)
        
        fields = result.strip().split("\t")
        assert fields[3] == "my_interval"

    def test_intervals_to_bed_multiple_intervals(self):
        """Test intervals_to_bed with multiple intervals."""
        intervals = [
            GenomicInterval(chrom="chr1", start=100, end=200),
            GenomicInterval(chrom="chr1", start=300, end=400),
            GenomicInterval(chrom="chr2", start=500, end=600),
        ]
        result = intervals_to_bed(intervals)
        
        lines = result.strip().split("\n")
        assert len(lines) == 3

    def test_intervals_to_bed_to_file(self):
        """Test intervals_to_bed writes to file."""
        intervals = [GenomicInterval(chrom="chr1", start=100, end=200)]
        
        with tempfile.NamedTemporaryFile(
            mode='w', suffix='.bed', delete=False
        ) as f:
            fpath = f.name
        
        try:
            intervals_to_bed(intervals, output_path=fpath)
            
            with open(fpath, 'r') as f:
                content = f.read()
            
            assert "chr1" in content
            assert "100" in content
            assert "200" in content
        finally:
            os.unlink(fpath)

    def test_intervals_to_bed_to_file_object(self):
        """Test intervals_to_bed writes to file-like object."""
        intervals = [GenomicInterval(chrom="chr1", start=100, end=200)]
        
        output = StringIO()
        intervals_to_bed(intervals, output_path=output)
        
        content = output.getvalue()
        assert "chr1" in content

    def test_intervals_to_bed_no_strand(self):
        """Test intervals_to_bed with None strand uses '.'."""
        intervals = [GenomicInterval(chrom="chr1", start=100, end=200, strand=None)]
        result = intervals_to_bed(intervals)
        
        fields = result.strip().split("\t")
        assert fields[5] == "."


# ============================================================================
# is_fusion_read tests
# ============================================================================

class TestIsFusionRead:
    """Tests for is_fusion_read function."""

    def test_is_fusion_read_supplementary(self):
        """Test supplementary reads are identified as fusion."""
        read_dict = {"is_supplementary": True, "cigartuples": [(0, 100)]}
        assert is_fusion_read(read_dict) is True

    def test_is_fusion_read_not_supplementary(self):
        """Test non-supplementary reads without long soft-clips."""
        read_dict = {"is_supplementary": False, "cigartuples": [(0, 100)]}
        assert is_fusion_read(read_dict) is False

    def test_is_fusion_read_long_soft_clip_start(self):
        """Test reads with long soft-clip at start are fusion."""
        # CIGAR operation 4 is soft-clip
        read_dict = {"is_supplementary": False, "cigartuples": [(4, 60), (0, 100)]}
        assert is_fusion_read(read_dict) is True

    def test_is_fusion_read_long_soft_clip_end(self):
        """Test reads with long soft-clip at end are fusion."""
        read_dict = {"is_supplementary": False, "cigartuples": [(0, 100), (4, 60)]}
        assert is_fusion_read(read_dict) is True

    def test_is_fusion_read_short_soft_clip(self):
        """Test reads with short soft-clips are not fusion."""
        read_dict = {"is_supplementary": False, "cigartuples": [(4, 30), (0, 100)]}
        assert is_fusion_read(read_dict) is False

    def test_is_fusion_read_no_cigar(self):
        """Test reads without cigartuples."""
        read_dict = {"is_supplementary": False, "cigartuples": None}
        assert is_fusion_read(read_dict) is False

    def test_is_fusion_read_empty_cigar(self):
        """Test reads with empty cigartuples."""
        read_dict = {"is_supplementary": False, "cigartuples": []}
        assert is_fusion_read(read_dict) is False


# ============================================================================
# extract_strand_from_read tests
# ============================================================================

class TestExtractStrandFromRead:
    """Tests for extract_strand_from_read function."""

    def test_extract_strand_forward(self):
        """Test forward strand extraction."""
        read_dict = {"is_forward": True}
        assert extract_strand_from_read(read_dict) == "+"

    def test_extract_strand_reverse(self):
        """Test reverse strand extraction."""
        read_dict = {"is_forward": False}
        assert extract_strand_from_read(read_dict) == "-"

    def test_extract_strand_missing_is_forward(self):
        """Test when is_forward is missing (defaults to falsy -> '-')."""
        read_dict = {}
        assert extract_strand_from_read(read_dict) == "-"


# ============================================================================
# extract_three_prime_pos tests
# ============================================================================

class TestExtractThreePrimePos:
    """Tests for extract_three_prime_pos function."""

    def test_extract_three_prime_forward(self):
        """Test 3' position for forward strand."""
        read_dict = {"is_forward": True, "reference_start": 100, "reference_end": 200}
        assert extract_three_prime_pos(read_dict) == 200

    def test_extract_three_prime_reverse(self):
        """Test 3' position for reverse strand."""
        read_dict = {"is_forward": False, "reference_start": 100, "reference_end": 200}
        assert extract_three_prime_pos(read_dict) == 100


# ============================================================================
# Integration tests with real data (if available)
# ============================================================================

class TestIntervalManagerIntegration:
    """Integration tests for interval_manager module."""

    @pytest.fixture
    def sample_intervals(self):
        """Create a set of sample intervals for testing."""
        return [
            GenomicInterval(
                chrom="chr1",
                start=1000,
                end=2000,
                strand="+",
                read_count=100,
                interval_id="interval_1",
            ),
            GenomicInterval(
                chrom="chr1",
                start=3000,
                end=4000,
                strand="-",
                read_count=50,
                interval_id="interval_2",
            ),
            GenomicInterval(
                chrom="chr2",
                start=1000,
                end=2000,
                strand="+",
                read_count=75,
                interval_id="interval_3",
            ),
        ]

    def test_multiple_intervals_to_bed(self, sample_intervals):
        """Test converting multiple intervals to BED format."""
        result = intervals_to_bed(sample_intervals)
        lines = result.strip().split("\n")
        
        assert len(lines) == 3
        
        # Check first interval
        fields1 = lines[0].split("\t")
        assert fields1[0] == "chr1"
        assert fields1[1] == "1000"
        assert fields1[2] == "2000"
        assert fields1[3] == "interval_1"
        assert fields1[4] == "100"
        assert fields1[5] == "+"

    def test_find_overlapping_intervals(self, sample_intervals):
        """Test finding overlapping intervals."""
        test_interval = GenomicInterval(chrom="chr1", start=1500, end=2500)
        
        overlapping = [
            i for i in sample_intervals if i.overlaps(test_interval)
        ]
        
        assert len(overlapping) == 1
        assert overlapping[0].interval_id == "interval_1"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
