"""Tests for multi-sample quantification components."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pysam

from fin.analysis.quantification import QuantResult, compute_tpm
from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval
from fin.pipeline.quantify_runner import IntervalAssignment, intervals_from_gtf


# ---------------------------------------------------------------------------
# compute_tpm
# ---------------------------------------------------------------------------

class TestComputeTpm:
    def test_basic_tpm_sums_to_1e6(self):
        results = {
            "tx1": QuantResult(candidate_id="tx1", abundance=100, confidence=0.9, num_assigned_reads=80, source="gtf"),
            "tx2": QuantResult(candidate_id="tx2", abundance=200, confidence=0.8, num_assigned_reads=160, source="gtf"),
        }
        lengths = {"tx1": 1000, "tx2": 2000}

        tpm = compute_tpm(results, lengths)

        assert abs(sum(tpm.values()) - 1e6) < 0.01
        assert len(tpm) == 2
        # tx1: RPK = 100/1 = 100, tx2: RPK = 200/2 = 100 → equal TPM
        assert abs(tpm["tx1"] - tpm["tx2"]) < 0.01

    def test_zero_abundance(self):
        results = {
            "tx1": QuantResult(candidate_id="tx1", abundance=0, confidence=0, num_assigned_reads=0, source="gtf"),
            "tx2": QuantResult(candidate_id="tx2", abundance=50, confidence=0.5, num_assigned_reads=40, source="gtf"),
        }
        lengths = {"tx1": 1000, "tx2": 500}

        tpm = compute_tpm(results, lengths)

        assert tpm["tx1"] == 0.0
        assert abs(tpm["tx2"] - 1e6) < 0.01

    def test_all_zero_abundance(self):
        results = {
            "tx1": QuantResult(candidate_id="tx1", abundance=0, confidence=0, num_assigned_reads=0, source="gtf"),
        }
        lengths = {"tx1": 1000}

        tpm = compute_tpm(results, lengths)

        assert tpm["tx1"] == 0.0

    def test_missing_length_treated_as_zero(self):
        results = {
            "tx1": QuantResult(candidate_id="tx1", abundance=100, confidence=0.9, num_assigned_reads=80, source="gtf"),
            "tx2": QuantResult(candidate_id="tx2", abundance=50, confidence=0.8, num_assigned_reads=40, source="gtf"),
        }
        lengths = {"tx1": 1000}  # tx2 missing

        tpm = compute_tpm(results, lengths)

        assert tpm["tx2"] == 0.0
        assert abs(tpm["tx1"] - 1e6) < 0.01


# ---------------------------------------------------------------------------
# intervals_from_gtf
# ---------------------------------------------------------------------------

def _make_mock_gene(gene_id, chrom, start, end, strand="+"):
    gene = MagicMock()
    gene.gene_id = gene_id
    gene.chrom = chrom
    gene.start = start
    gene.end = end
    gene.strand = strand
    return gene


class TestIntervalsFromGtf:
    def test_non_overlapping_genes(self):
        genes = [
            _make_mock_gene("g1", "chr1", 100, 500, "+"),
            _make_mock_gene("g2", "chr1", 1000, 2000, "+"),
        ]
        reader = MagicMock()
        reader.iterate_genes.return_value = iter(genes)

        intervals = intervals_from_gtf(reader)

        assert len(intervals) == 2
        assert intervals[0].start == 100
        assert intervals[0].end == 500
        assert intervals[1].start == 1000
        assert intervals[1].end == 2000

    def test_overlapping_genes_merge(self):
        genes = [
            _make_mock_gene("g1", "chr1", 100, 500, "+"),
            _make_mock_gene("g2", "chr1", 400, 800, "+"),
            _make_mock_gene("g3", "chr1", 2000, 3000, "+"),
        ]
        reader = MagicMock()
        reader.iterate_genes.return_value = iter(genes)

        intervals = intervals_from_gtf(reader)

        assert len(intervals) == 2
        assert intervals[0].start == 100
        assert intervals[0].end == 800
        assert intervals[1].start == 2000
        assert intervals[1].end == 3000

    def test_multi_chromosome(self):
        genes = [
            _make_mock_gene("g1", "chr1", 100, 500, "+"),
            _make_mock_gene("g2", "chr2", 200, 600, "-"),
        ]
        reader = MagicMock()
        reader.iterate_genes.return_value = iter(genes)

        intervals = intervals_from_gtf(reader)

        assert len(intervals) == 2
        chroms = {iv.chrom for iv in intervals}
        assert chroms == {"chr1", "chr2"}

    def test_empty_gtf(self):
        reader = MagicMock()
        reader.iterate_genes.return_value = iter([])

        intervals = intervals_from_gtf(reader)

        assert intervals == []

    def test_skip_invalid_gene(self):
        """Genes with start >= end should be skipped."""
        genes = [
            _make_mock_gene("g1", "chr1", 500, 500, "+"),  # zero-length
            _make_mock_gene("g2", "chr1", 100, 200, "+"),
        ]
        reader = MagicMock()
        reader.iterate_genes.return_value = iter(genes)

        intervals = intervals_from_gtf(reader)

        assert len(intervals) == 1
        assert intervals[0].start == 100


# ---------------------------------------------------------------------------
# discover_gtf_only (unit test with mocks)
# ---------------------------------------------------------------------------

class TestDiscoverGtfOnly:
    def test_returns_only_gtf_candidates(self):
        """Verify discover_gtf_only returns only GTF candidates, no novel."""
        from unittest.mock import patch

        from fin.candidates.discovery import discover_gtf_only

        interval = GenomicInterval(chrom="chr1", start=100, end=1000, strand="+")

        # Mock GTF transcript
        mock_tx = MagicMock()
        mock_tx.transcript_id = "ENST001"
        mock_tx.strand = "+"
        mock_tx.start = 100
        mock_tx.end = 500
        mock_tx.exons = [(100, 200), (300, 500)]
        mock_tx.get_spliced_sequence.return_value = "ATCG" * 50
        mock_tx.sort_features.return_value = None

        mock_gtf_reader = MagicMock()
        mock_gtf_reader.get_transcripts_in_region.return_value = [mock_tx]

        # Mock BamReader to return some reads
        mock_bam_instance = MagicMock()
        mock_alignment = MagicMock()
        mock_rd = {"is_mapped": True, "query_name": "read1", "is_forward": True}
        mock_bam_instance.fetch.return_value = [mock_alignment]
        mock_bam_instance.alignment_to_dict.return_value = mock_rd
        mock_bam_instance.__enter__ = MagicMock(return_value=mock_bam_instance)
        mock_bam_instance.__exit__ = MagicMock(return_value=False)

        with patch("fin.io.io_bam.BamReader", return_value=mock_bam_instance):
            result = discover_gtf_only(
                interval=interval,
                bam_path="fake.bam",
                gtf_reader=mock_gtf_reader,
                genome_fasta="A" * 2000,
            )

        assert isinstance(result, CandidateSet)
        assert result.num_candidates == 1
        assert result.candidates[0].source == "gtf"
        assert result.candidates[0].candidate_id == "ENST001"
        assert "read1" in result.read_ids

    def test_rejects_opposite_strand_reads(self):
        """A strand-pure interval must drop reads mapped to the other strand.

        bam.fetch returns reads on BOTH strands overlapping the region; each
        read belongs to its own mapped strand's interval. A '-' read fetched
        inside a '+' interval must not be swallowed (it belongs to the '-'
        interval), otherwise sense/antisense overlaps corrupt strand assignment.
        """
        from unittest.mock import patch

        from fin.candidates.discovery import discover_gtf_only

        interval = GenomicInterval(chrom="chr1", start=100, end=1000, strand="+")

        mock_gtf_reader = MagicMock()
        mock_gtf_reader.get_transcripts_in_region.return_value = []

        # One '+' read (kept) and one '-' read (must be rejected).
        fwd = MagicMock()
        rev = MagicMock()
        mock_bam_instance = MagicMock()
        mock_bam_instance.fetch.return_value = [fwd, rev]
        mock_bam_instance.alignment_to_dict.side_effect = [
            {"is_mapped": True, "query_name": "fwd_read", "is_forward": True},
            {"is_mapped": True, "query_name": "rev_read", "is_forward": False},
        ]
        mock_bam_instance.__enter__ = MagicMock(return_value=mock_bam_instance)
        mock_bam_instance.__exit__ = MagicMock(return_value=False)

        with patch("fin.io.io_bam.BamReader", return_value=mock_bam_instance):
            result = discover_gtf_only(
                interval=interval,
                bam_path="fake.bam",
                gtf_reader=mock_gtf_reader,
                genome_fasta="A" * 2000,
            )

        assert "fwd_read" in result.read_ids
        assert "rev_read" not in result.read_ids


# ---------------------------------------------------------------------------
# _build_assignment_bam
# ---------------------------------------------------------------------------

def _make_candidate(candidate_id, sequence="ATCGATCG" * 50):
    """Create a TranscriptCandidate with minimal fields."""
    return TranscriptCandidate(
        candidate_id=candidate_id,
        intron_chain=IntronChain(introns=()),
        three_prime_pos=len(sequence),
        sequence=sequence,
        source="gtf",
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=0,
        end=len(sequence),
    )


def _write_candidate_bam(bam_path: Path, candidate_id: str, seq_len: int, reads):
    """Write a single-reference BAM with given reads.

    reads: list of (query_name, query_sequence, ref_start)
    """
    bam_path.parent.mkdir(parents=True, exist_ok=True)
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": candidate_id, "LN": seq_len}],
        }
    )
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:
        for qname, qseq, ref_start in sorted(reads, key=lambda r: r[2]):
            a = pysam.AlignedSegment(header)
            a.query_name = qname
            a.query_sequence = qseq
            a.flag = 0
            a.reference_id = 0
            a.reference_start = ref_start
            a.mapping_quality = 60
            a.cigar = [(0, len(qseq))]  # simple M cigar
            a.query_qualities = pysam.qualitystring_to_array("I" * len(qseq))
            outf.write(a)
    pysam.index(str(bam_path))


class TestBuildAssignmentBam:
    def test_merges_assigned_reads(self, tmp_path):
        """Reads assigned to different candidates appear in output with correct refs."""
        from fin.pipeline.quantify_runner import QuantifyRunner, SampleInput

        seq = "ATCGATCG" * 50  # 400bp
        cand_a = _make_candidate("txA", seq)
        cand_b = _make_candidate("txB", seq)

        # Set up per-candidate BAMs in work_dir
        # In practice, mappy aligns ALL reads to EVERY candidate
        work_dir = tmp_path / "work"
        all_reads = [
            ("read1", "ATCGATCG", 0),
            ("read2", "ATCGATCG", 10),
            ("read3", "ATCGATCG", 20),
            ("read4", "ATCGATCG", 5),
        ]
        _write_candidate_bam(
            work_dir / "txA" / "aligned.bam", "txA", len(seq), all_reads
        )
        _write_candidate_bam(
            work_dir / "txB" / "aligned.bam", "txB", len(seq), all_reads
        )

        # read1 → txA (idx 0), read2 → txB (idx 1), read3 → txB, read4 → txA
        assignment = IntervalAssignment(
            read_ids=["read1", "read2", "read3", "read4"],
            hard_assignments=np.array([0, 1, 1, 0]),
            candidates=[cand_a, cand_b],
            work_dir=work_dir,
        )

        sample = SampleInput(name="s1", bam_path="", fastq_path="", signal_path="")
        output_bam = tmp_path / "s1.bam"

        # Create a minimal runner just to call the method
        runner = QuantifyRunner.__new__(QuantifyRunner)
        runner._build_assignment_bam(sample, [assignment], output_bam)

        # Verify output BAM
        assert output_bam.exists()
        assert (tmp_path / "s1.bam.bai").exists()

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            refs = {bam.get_reference_name(i) for i in range(bam.nreferences)}
            assert refs == {"txA", "txB"}

            reads_by_ref = {}
            for aln in bam:
                ref_name = bam.get_reference_name(aln.reference_id)
                reads_by_ref.setdefault(ref_name, set()).add(aln.query_name)

        # read1, read4 → txA; read2, read3 → txB
        assert reads_by_ref["txA"] == {"read1", "read4"}
        assert reads_by_ref["txB"] == {"read2", "read3"}

    def test_skips_unassigned_reads(self, tmp_path):
        """Reads with negative assignment index are omitted."""
        from fin.pipeline.quantify_runner import QuantifyRunner, SampleInput

        seq = "ATCGATCG" * 50
        cand = _make_candidate("tx1", seq)

        work_dir = tmp_path / "work"
        _write_candidate_bam(
            work_dir / "tx1" / "aligned.bam",
            "tx1",
            len(seq),
            [("read1", "ATCGATCG", 0), ("read2", "ATCGATCG", 10)],
        )

        # read1 assigned to tx1, read2 unassigned (-1)
        assignment = IntervalAssignment(
            read_ids=["read1", "read2"],
            hard_assignments=np.array([0, -1]),
            candidates=[cand],
            work_dir=work_dir,
        )

        sample = SampleInput(name="s1", bam_path="", fastq_path="", signal_path="")
        output_bam = tmp_path / "s1.bam"

        runner = QuantifyRunner.__new__(QuantifyRunner)
        runner._build_assignment_bam(sample, [assignment], output_bam)

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            read_names = {aln.query_name for aln in bam}

        assert read_names == {"read1"}

    def test_multiple_intervals(self, tmp_path):
        """Reads from multiple intervals are merged into a single BAM."""
        from fin.pipeline.quantify_runner import QuantifyRunner, SampleInput

        seq = "ATCGATCG" * 50
        cand_a = _make_candidate("txA", seq)
        cand_b = _make_candidate("txB", seq)

        work_dir_1 = tmp_path / "interval1"
        work_dir_2 = tmp_path / "interval2"

        _write_candidate_bam(
            work_dir_1 / "txA" / "aligned.bam",
            "txA",
            len(seq),
            [("r1", "ATCGATCG", 0)],
        )
        _write_candidate_bam(
            work_dir_2 / "txB" / "aligned.bam",
            "txB",
            len(seq),
            [("r2", "ATCGATCG", 0)],
        )

        ia1 = IntervalAssignment(
            read_ids=["r1"],
            hard_assignments=np.array([0]),
            candidates=[cand_a],
            work_dir=work_dir_1,
        )
        ia2 = IntervalAssignment(
            read_ids=["r2"],
            hard_assignments=np.array([0]),
            candidates=[cand_b],
            work_dir=work_dir_2,
        )

        sample = SampleInput(name="s1", bam_path="", fastq_path="", signal_path="")
        output_bam = tmp_path / "s1.bam"

        runner = QuantifyRunner.__new__(QuantifyRunner)
        runner._build_assignment_bam(sample, [ia1, ia2], output_bam)

        with pysam.AlignmentFile(str(output_bam), "rb") as bam:
            reads_by_ref = {}
            for aln in bam:
                ref_name = bam.get_reference_name(aln.reference_id)
                reads_by_ref.setdefault(ref_name, set()).add(aln.query_name)

        assert reads_by_ref == {"txA": {"r1"}, "txB": {"r2"}}
