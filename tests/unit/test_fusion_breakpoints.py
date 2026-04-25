"""Unit tests for fin.fusion.breakpoints — SA tag parsing and clustering."""

from __future__ import annotations

import os

import pysam

from fin.fusion.breakpoints import Breakpoint, cluster_breakpoints, parse_sa_tags


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_bam(tmp_path: str, reads: list[dict]) -> str:
    """Write a minimal BAM file containing the given synthetic reads.

    Each *read* dict may have:
        name (str), chrom (str), pos (int, 0-based), mapq (int),
        reverse (bool), sa_tag (str or None), seq (str).
    """
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [
                {"SN": "chr1", "LN": 100_000},
                {"SN": "chr2", "LN": 100_000},
                {"SN": "chr3", "LN": 100_000},
            ],
        }
    )

    bam_path = os.path.join(tmp_path, "test.bam")
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for r in reads:
            seg = pysam.AlignedSegment(header)
            seg.query_name = r["name"]
            seg.query_sequence = r.get("seq", "A" * 50)
            seg.flag = 0
            if r.get("reverse", False):
                seg.flag |= 0x10
            ref_id = header.get_tid(r["chrom"])
            seg.reference_id = ref_id
            seg.reference_start = r["pos"]
            seg.mapping_quality = r["mapq"]
            seg.cigar = [(0, len(seg.query_sequence))]  # simple match
            seg.query_qualities = pysam.qualitystring_to_array("I" * len(seg.query_sequence))
            if r.get("sa_tag"):
                seg.set_tag("SA", r["sa_tag"])
            bam.write(seg)

    pysam.sort("-o", bam_path + ".sorted.bam", bam_path)
    pysam.index(bam_path + ".sorted.bam")
    return bam_path + ".sorted.bam"


# ---------------------------------------------------------------------------
# Tests: parse_sa_tags
# ---------------------------------------------------------------------------

class TestParseSaTags:
    def test_parse_sa_tags_basic(self, tmp_path):
        """AC-01: 3 reads (2 with SA tags, 1 without) → exactly 2 Breakpoints."""
        reads = [
            {
                "name": "read1",
                "chrom": "chr1",
                "pos": 1000,
                "mapq": 60,
                "reverse": False,
                # SA: rname,pos,strand,CIGAR,mapQ,NM  (pos is 1-based in SA tag)
                "sa_tag": "chr2,5001,+,50M,60,0;",
            },
            {
                "name": "read2",
                "chrom": "chr1",
                "pos": 2000,
                "mapq": 40,
                "reverse": True,
                "sa_tag": "chr3,8001,-,50M,40,2;",
            },
            {
                "name": "read3",
                "chrom": "chr1",
                "pos": 3000,
                "mapq": 60,
                "reverse": False,
                # No SA tag → should produce no Breakpoint
            },
        ]
        bam_path = _write_bam(str(tmp_path), reads)
        result = parse_sa_tags(bam_path, min_mapq=10)

        assert len(result) == 2, f"Expected 2 breakpoints, got {len(result)}"

        # Verify read1 breakpoint
        bp1 = next(b for b in result if b.supporting_read_ids == {"read1"})
        assert bp1.chromA == "chr1"
        assert bp1.posA == 1000
        assert bp1.strandA == "+"
        assert bp1.chromB == "chr2"
        assert bp1.posB == 5000  # SA pos 5001 → 0-based 5000
        assert bp1.strandB == "+"

        # Verify read2 breakpoint
        bp2 = next(b for b in result if b.supporting_read_ids == {"read2"})
        assert bp2.chromA == "chr1"
        assert bp2.posA == 2000
        assert bp2.strandA == "-"
        assert bp2.chromB == "chr3"
        assert bp2.posB == 8000  # SA pos 8001 → 0-based 8000
        assert bp2.strandB == "-"

    def test_parse_sa_tags_mapq_filter(self, tmp_path):
        """Primary MAPQ=5 should be skipped when min_mapq=10."""
        reads = [
            {
                "name": "lowq",
                "chrom": "chr1",
                "pos": 500,
                "mapq": 5,
                "reverse": False,
                "sa_tag": "chr2,1001,+,50M,60,0;",
            },
        ]
        bam_path = _write_bam(str(tmp_path), reads)
        result = parse_sa_tags(bam_path, min_mapq=10)
        assert result == [], "Low-MAPQ primary read should be filtered out"

    def test_parse_sa_tags_sa_mapq_filter(self, tmp_path):
        """SA entry with mapQ below threshold should be skipped even if primary passes."""
        reads = [
            {
                "name": "read_ok_primary",
                "chrom": "chr1",
                "pos": 500,
                "mapq": 60,
                "reverse": False,
                "sa_tag": "chr2,1001,+,50M,5,0;",  # SA mapQ = 5 < 10
            },
        ]
        bam_path = _write_bam(str(tmp_path), reads)
        result = parse_sa_tags(bam_path, min_mapq=10)
        assert result == [], "Low-MAPQ SA entry should be filtered out"

    def test_parse_sa_tags_no_sa_reads(self, tmp_path):
        """BAM with no SA tags at all returns empty list."""
        reads = [
            {"name": "r1", "chrom": "chr1", "pos": 100, "mapq": 60, "reverse": False},
        ]
        bam_path = _write_bam(str(tmp_path), reads)
        assert parse_sa_tags(bam_path) == []


# ---------------------------------------------------------------------------
# Tests: cluster_breakpoints
# ---------------------------------------------------------------------------

class TestClusterBreakpoints:
    def _bp(
        self,
        chromA: str,
        posA: int,
        strandA: str,
        chromB: str,
        posB: int,
        strandB: str,
        read_id: str = "r",
    ) -> Breakpoint:
        return Breakpoint(
            chromA=chromA,
            posA=posA,
            strandA=strandA,
            chromB=chromB,
            posB=posB,
            strandB=strandB,
            support_count=1,
            supporting_read_ids={read_id},
        )

    def test_cluster_breakpoints_basic(self):
        """AC-02: 5 breakpoints → 2 clusters with support_count 3 and 2."""
        # Group A: 3 close breakpoints
        bps = [
            self._bp("chr1", 1000, "+", "chr2", 5000, "+", "rA1"),
            self._bp("chr1", 1050, "+", "chr2", 5020, "+", "rA2"),
            self._bp("chr1", 1100, "+", "chr2", 4980, "+", "rA3"),
            # Group B: 2 close breakpoints on different chroms
            self._bp("chr3", 2000, "-", "chr4", 8000, "-", "rB1"),
            self._bp("chr3", 2050, "-", "chr4", 8010, "-", "rB2"),
        ]
        result = cluster_breakpoints(bps, max_dist=500, min_support=2)

        assert len(result) == 2, f"Expected 2 clusters, got {len(result)}"
        counts = sorted(r.support_count for r in result)
        assert counts == [2, 3], f"Expected support counts [2, 3], got {counts}"

    def test_cluster_min_support_filter(self):
        """AC-03: singleton cluster dropped when min_support=2."""
        bps = [
            # Singleton
            self._bp("chr1", 100, "+", "chr2", 200, "+", "lone"),
            # Pair that should survive
            self._bp("chr1", 1000, "+", "chr2", 5000, "+", "p1"),
            self._bp("chr1", 1010, "+", "chr2", 5010, "+", "p2"),
        ]
        result = cluster_breakpoints(bps, max_dist=500, min_support=2)
        assert len(result) == 1
        assert result[0].support_count == 2

    def test_cluster_preserves_read_ids(self):
        """Merged breakpoint supporting_read_ids is the union of all members."""
        bps = [
            self._bp("chr1", 1000, "+", "chr2", 5000, "+", "readA"),
            self._bp("chr1", 1010, "+", "chr2", 5010, "+", "readB"),
            self._bp("chr1", 1020, "+", "chr2", 5020, "+", "readC"),
        ]
        result = cluster_breakpoints(bps, max_dist=500, min_support=2)
        assert len(result) == 1
        assert result[0].supporting_read_ids == {"readA", "readB", "readC"}
        assert result[0].support_count == 3

    def test_cluster_empty_input(self):
        assert cluster_breakpoints([], max_dist=500, min_support=2) == []

    def test_cluster_different_strands_not_merged(self):
        """Breakpoints with different strands must not be merged."""
        bps = [
            self._bp("chr1", 1000, "+", "chr2", 5000, "+", "r1"),
            self._bp("chr1", 1000, "-", "chr2", 5000, "+", "r2"),  # strandA differs
        ]
        # Both are singletons → both dropped at min_support=2
        result = cluster_breakpoints(bps, max_dist=500, min_support=2)
        assert result == []

    def test_cluster_mean_position(self):
        """Representative positions are rounded means of cluster members."""
        bps = [
            self._bp("chr1", 1000, "+", "chr2", 5000, "+", "r1"),
            self._bp("chr1", 1100, "+", "chr2", 5100, "+", "r2"),
        ]
        result = cluster_breakpoints(bps, max_dist=500, min_support=2)
        assert len(result) == 1
        assert result[0].posA == 1050
        assert result[0].posB == 5050
