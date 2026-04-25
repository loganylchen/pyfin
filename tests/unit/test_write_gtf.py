"""Tests for GTF output writer."""

import os

import pytest

from fin.analysis.quantification import QuantResult
from fin.io.io_gtf import write_gtf


def _make_result(candidate_id, gene_id, chrom, strand, start, end, exons, **kwargs):
    defaults = dict(abundance=1.0, confidence=0.9, num_assigned_reads=10, source="gtf")
    defaults.update(kwargs)
    return QuantResult(
        candidate_id=candidate_id,
        gene_id=gene_id,
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
        exons=exons,
        **defaults,
    )


@pytest.fixture
def sample_results():
    return {
        "tx1": _make_result(
            "tx1", "gene1", "chr1", "+", 100, 500,
            exons=((100, 200), (300, 500)),
            abundance=5.5, confidence=0.85, num_assigned_reads=42,
        ),
        "tx2": _make_result(
            "tx2", "gene1", "chr1", "+", 100, 600,
            exons=((100, 200), (300, 400), (500, 600)),
            abundance=3.2, confidence=0.72, num_assigned_reads=15,
        ),
        "tx3": _make_result(
            "tx3", "gene2", "chr2", "-", 1000, 2000,
            exons=((1000, 2000),),
            abundance=10.0, confidence=0.99, num_assigned_reads=100,
            source="novel",
        ),
    }


class TestWriteGtf:
    def test_creates_file(self, sample_results, tmp_path):
        out = str(tmp_path / "out.gtf")
        write_gtf(sample_results, out)
        assert os.path.exists(out)

    def test_tab_separated_columns(self, sample_results, tmp_path):
        out = str(tmp_path / "out.gtf")
        write_gtf(sample_results, out)
        with open(out) as f:
            for line in f:
                fields = line.strip().split("\t")
                assert len(fields) == 9, f"Expected 9 tab-separated fields, got {len(fields)}: {line}"

    def test_one_based_coordinates(self, sample_results, tmp_path):
        out = str(tmp_path / "out.gtf")
        write_gtf(sample_results, out)
        with open(out) as f:
            for line in f:
                fields = line.strip().split("\t")
                start = int(fields[3])
                assert start >= 1, f"GTF start should be 1-based, got {start}"

    def test_gene_transcript_exon_hierarchy(self, sample_results, tmp_path):
        out = str(tmp_path / "out.gtf")
        write_gtf(sample_results, out)
        with open(out) as f:
            features = [line.strip().split("\t")[2] for line in f]

        # gene1 has 2 transcripts: gene, tx, 2 exons, tx, 3 exons = 1+1+2+1+3 = 8
        # gene2 has 1 transcript: gene, tx, 1 exon = 3
        assert features.count("gene") == 2
        assert features.count("transcript") == 3
        assert features.count("exon") == 6  # 2 + 3 + 1

    def test_gene_line_spans_all_transcripts(self, sample_results, tmp_path):
        out = str(tmp_path / "out.gtf")
        write_gtf(sample_results, out)
        with open(out) as f:
            for line in f:
                fields = line.strip().split("\t")
                if fields[2] == "gene" and 'gene1' in fields[8]:
                    # gene1 spans tx1 (100-500) and tx2 (100-600) -> 1-based: 101..600
                    assert int(fields[3]) == 101
                    assert int(fields[4]) == 600

    def test_transcript_attributes(self, sample_results, tmp_path):
        out = str(tmp_path / "out.gtf")
        write_gtf(sample_results, out)
        with open(out) as f:
            for line in f:
                fields = line.strip().split("\t")
                if fields[2] == "transcript" and 'tx1' in fields[8]:
                    attrs = fields[8]
                    assert 'abundance "5.5000"' in attrs
                    assert 'confidence "0.8500"' in attrs
                    assert 'num_reads "42"' in attrs
                    assert 'transcript_source "gtf"' in attrs
                    return
        pytest.fail("tx1 transcript line not found")

    def test_sorted_by_chrom_start(self, sample_results, tmp_path):
        out = str(tmp_path / "out.gtf")
        write_gtf(sample_results, out)
        with open(out) as f:
            chroms = []
            for line in f:
                fields = line.strip().split("\t")
                if fields[2] == "gene":
                    chroms.append(fields[0])
        assert chroms == ["chr1", "chr2"]

    def test_exon_numbers(self, sample_results, tmp_path):
        out = str(tmp_path / "out.gtf")
        write_gtf(sample_results, out)
        with open(out) as f:
            tx2_exon_nums = []
            for line in f:
                fields = line.strip().split("\t")
                if fields[2] == "exon" and 'tx2' in fields[8]:
                    # Extract exon_number
                    for part in fields[8].split(";"):
                        part = part.strip()
                        if part.startswith("exon_number"):
                            num = part.split('"')[1]
                            tx2_exon_nums.append(int(num))
            assert tx2_exon_nums == [1, 2, 3]

    def test_single_exon_transcript(self, tmp_path):
        results = {
            "mono": _make_result(
                "mono", "gMono", "chrX", "+", 0, 100,
                exons=((0, 100),),
            ),
        }
        out = str(tmp_path / "out.gtf")
        write_gtf(results, out)
        with open(out) as f:
            lines = f.readlines()
        assert len(lines) == 3  # gene + transcript + 1 exon

    def test_creates_parent_directories(self, tmp_path):
        out = str(tmp_path / "sub" / "dir" / "out.gtf")
        results = {
            "t1": _make_result("t1", "g1", "chr1", "+", 0, 100, exons=((0, 100),)),
        }
        write_gtf(results, out)
        assert os.path.exists(out)

    def test_empty_results(self, tmp_path):
        out = str(tmp_path / "empty.gtf")
        write_gtf({}, out)
        with open(out) as f:
            assert f.read() == ""


class TestWriteGtfScoresAndFusion:
    """AC-18: score attributes and fusion chrom splitting."""

    def _write_and_read_transcript_line(self, qr: "QuantResult", tmp_path) -> str:
        """Helper: write single result, return the transcript line's attribute field."""
        out = str(tmp_path / "out.gtf")
        write_gtf({qr.candidate_id: qr}, out)
        with open(out) as f:
            for line in f:
                fields = line.strip().split("\t")
                if fields[2] == "transcript":
                    return fields[8]
        pytest.fail("No transcript line found in output")

    def test_gtf_emits_score_attributes(self, tmp_path):
        qr = _make_result(
            "tx_s", "g_s", "chr1", "+", 0, 100, exons=((0, 100),),
            coherence_score=0.1234, discrimination_score=0.5678, combined_score=0.3456,
        )
        attrs = self._write_and_read_transcript_line(qr, tmp_path)
        assert 'coherence_score "0.1234"' in attrs
        assert 'discrimination_score "0.5678"' in attrs
        assert 'combined_score "0.3456"' in attrs

    def test_gtf_score_attributes_order(self, tmp_path):
        """Scores must appear AFTER transcript_source (strict append-only, m6)."""
        qr = _make_result(
            "tx_ord", "g_ord", "chr1", "+", 0, 100, exons=((0, 100),),
            coherence_score=0.1, discrimination_score=0.2, combined_score=0.3,
        )
        attrs = self._write_and_read_transcript_line(qr, tmp_path)
        src_pos = attrs.index("transcript_source")
        coh_pos = attrs.index("coherence_score")
        disc_pos = attrs.index("discrimination_score")
        comb_pos = attrs.index("combined_score")
        assert src_pos < coh_pos < disc_pos < comb_pos, (
            f"Score attrs must follow transcript_source. Positions: "
            f"transcript_source={src_pos}, coherence={coh_pos}, "
            f"discrimination={disc_pos}, combined={comb_pos}"
        )

    def test_gtf_fusion_seqname_split(self, tmp_path):
        qr = _make_result(
            "tx_f", "g_f", "chr1::chr2", "+", 0, 100, exons=((0, 100),),
            source="fusion",
        )
        out = str(tmp_path / "out.gtf")
        write_gtf({qr.candidate_id: qr}, out)
        with open(out) as f:
            for line in f:
                fields = line.strip().split("\t")
                if fields[2] == "transcript":
                    assert fields[0] == "chr1", f"seqname should be chr1, got {fields[0]}"
                    assert 'fusion_partner_chrom "chr2"' in fields[8], (
                        f"fusion_partner_chrom not found in: {fields[8]}"
                    )
                    return
        pytest.fail("No transcript line found")

    def test_gtf_non_fusion_unchanged(self, tmp_path):
        qr = _make_result(
            "tx_nf", "g_nf", "chr5", "+", 0, 100, exons=((0, 100),),
            source="gtf",
        )
        out = str(tmp_path / "out.gtf")
        write_gtf({qr.candidate_id: qr}, out)
        with open(out) as f:
            for line in f:
                fields = line.strip().split("\t")
                if fields[2] == "transcript":
                    assert fields[0] == "chr5", f"chrom should be chr5, got {fields[0]}"
                    assert "fusion_partner_chrom" not in fields[8], (
                        "fusion_partner_chrom should not appear for non-fusion records"
                    )
                    return
        pytest.fail("No transcript line found")

    def test_gtf_default_scores_zero_format(self, tmp_path):
        qr = _make_result(
            "tx_z", "g_z", "chr1", "+", 0, 100, exons=((0, 100),),
        )
        attrs = self._write_and_read_transcript_line(qr, tmp_path)
        assert 'coherence_score "0.0000"' in attrs
        assert 'discrimination_score "0.0000"' in attrs
        assert 'combined_score "0.0000"' in attrs
