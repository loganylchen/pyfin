import pysam
from fin.io.interval_manager import generate_intervals_from_reads, generate_isolated_intervals, is_fusion_read


def _write_synthetic_bam(tmp_path):
    """Build a minimal BAM with known fusion and non-fusion reads.

    Returns path to sorted+indexed BAM.
    """
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 10000}, {"SN": "chr2", "LN": 10000}],
    }
    unsorted = tmp_path / "test.unsorted.bam"
    bam_path = tmp_path / "test.bam"
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as f:
        for i, (qn, ref_start, cigar, is_suppl, sa_tag, md_tag) in enumerate([
            ("read_norm_1", 1000, [(0, 100)], False, None, "100"),
            ("read_norm_2", 1100, [(0, 100)], False, None, "100"),
            ("read_norm_3", 1200, [(0, 100)], False, None, "100"),
            # Fusion-by-supplementary: supplementary flag set
            ("read_fus_1", 1500, [(0, 80)], True, "chr2,5000,+,80M,60,0;", "80"),
            # Fusion-by-soft-clip: 260bp soft clip at start (>= 250 threshold)
            ("read_fus_2", 2000, [(4, 260), (0, 100)], False, None, "100"),
        ]):
            a = pysam.AlignedSegment(pysam.AlignmentHeader.from_dict(header))
            a.query_name = qn
            a.query_sequence = "A" * sum(L for op, L in cigar if op in (0, 4, 7, 8))
            a.flag = 2048 if is_suppl else 0
            a.reference_id = 0
            a.reference_start = ref_start
            a.mapping_quality = 60
            a.cigar = cigar
            a.set_tag("MD", md_tag, "Z")
            if sa_tag:
                a.set_tag("SA", sa_tag, "Z")
            f.write(a)

    pysam.sort("-o", str(bam_path), str(unsorted))
    pysam.index(str(bam_path))
    return str(bam_path)


def test_generate_intervals_returns_three_tuple(tmp_path):
    bam = _write_synthetic_bam(tmp_path)
    result = generate_intervals_from_reads(bam)
    assert len(result) == 3, "must return (read_alignments, fusion_read_ids, fusion_reads)"
    read_alignments, fusion_read_ids, fusion_reads = result
    assert isinstance(read_alignments, list)
    assert isinstance(fusion_read_ids, set)
    assert isinstance(fusion_reads, list)


def test_fusion_reads_retained(tmp_path):
    bam = _write_synthetic_bam(tmp_path)
    _, fusion_read_ids, fusion_reads = generate_intervals_from_reads(bam)
    # 2 fusion reads should be captured (supplementary + long soft-clip)
    assert len(fusion_reads) == 2
    # fusion_read_ids matches fusion_reads
    assert fusion_read_ids == {r["query_name"] for r in fusion_reads}
    # Each entry in fusion_reads satisfies is_fusion_read
    for rd in fusion_reads:
        assert is_fusion_read(rd)


def test_non_fusion_count_unchanged(tmp_path):
    bam = _write_synthetic_bam(tmp_path)
    read_alignments, _, _ = generate_intervals_from_reads(bam)
    # 3 normal reads, all non-fusion
    assert len(read_alignments) == 3
    for rd in read_alignments:
        assert not is_fusion_read(rd)


def test_isolated_intervals_includes_fusion_reads_key(tmp_path):
    bam = _write_synthetic_bam(tmp_path)
    result = generate_isolated_intervals(bam)
    assert "fusion_reads" in result
    assert isinstance(result["fusion_reads"], list)
    assert len(result["fusion_reads"]) == 2
    # existing keys still present
    assert "intervals" in result
    assert "fusion_read_ids" in result
    assert "num_reads_processed" in result
