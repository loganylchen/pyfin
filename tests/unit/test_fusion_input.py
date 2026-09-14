"""Fusion-only BAM evidence must be independent of assembly intervals."""

from pathlib import Path
from unittest.mock import patch

import pysam
import pytest

from fin.fusion.annotation import GeneIndex, GeneModel
from fin.fusion.input import FusionInputError, scan_fusion_bam
from fin.fusion.paths import find_read_paths
from fin.utils.sequences import reverse_complement


def _write_bam(path, entries):
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": 10000}
                                           for chrom in ("chr1", "chr2", "chr3")]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as writer:
        for entry in entries:
            record = pysam.AlignedSegment(writer.header)
            record.query_name = entry.get("name", "fusion")
            record.flag = entry.get("flag", 0)
            record.reference_id = writer.header.get_tid(entry.get("chrom", "chr1"))
            record.reference_start = entry.get("start", 100)
            record.mapping_quality = 60
            record.cigarstring = entry.get("cigar", "100M200S")
            query_length = sum(length for code, length in record.cigartuples if code in (0, 1, 4, 7, 8))
            record.query_sequence = entry.get("sequence", "A" * query_length)
            record.set_tag("NM", 0)
            if "sa" in entry:
                record.set_tag("SA", entry["sa"])
            writer.write(record)


def test_fusion_only_locus_is_scanned_without_index_annotation_or_intervals(tmp_path):
    bam = tmp_path / "fusion.bam"
    _write_bam(bam, [{"sa": "chr2,101,+,100S100M100S,60,0;chr3,101,+,200S100M,60,0;"}])
    before = bam.read_bytes()
    evidence = scan_fusion_bam(bam)
    assert evidence.complete is True
    assert len(evidence.molecules) == 1
    molecule, = evidence.molecules
    result = find_read_paths(molecule.segments)
    assert result.status == "unique"
    assert len(result.paths[0].segments) == 3
    assert molecule.sequence == "A" * 300
    assert bam.read_bytes() == before
    assert not Path(str(bam) + ".bai").exists()


def test_supplementary_records_are_found_even_without_sa_tags(tmp_path):
    bam = tmp_path / "fusion.bam"
    _write_bam(bam, [
        {"cigar": "100M200S"},
        {"chrom": "chr2", "cigar": "100H100M100H", "flag": 2048},
        {"chrom": "chr3", "cigar": "200H100M", "flag": 2048},
    ])
    evidence = scan_fusion_bam(bam)
    molecule, = evidence.molecules
    assert len(molecule.segments) == 3
    assert find_read_paths(molecule.segments).status == "unique"
    assert evidence.record_counts["nominated_molecules"] == 1


def test_hard_clipped_primary_without_sa_is_nominated_for_sequence_recovery(tmp_path):
    bam = tmp_path / "fusion.bam"
    # the nomination threshold is 250 bp of clip on real direct RNA (adapter
    # and poly(A) carry-through sit at ~80 bp and would nominate ~9% of a
    # library at 100); a 300-bp hard clip is above it
    _write_bam(bam, [{"cigar": "100M300H"}])
    molecule, = scan_fusion_bam(bam).molecules
    assert molecule.status == "needs_sequence"
    assert molecule.segments[0].read_length == 400


def test_secondary_normal_alignment_is_retained_as_competitor(tmp_path):
    bam = tmp_path / "fusion.bam"
    _write_bam(bam, [
        {"cigar": "100M100S", "sa": "chr2,101,+,100S100M,60,0;"},
        {"chrom": "chr3", "cigar": "200M", "flag": 256},
    ])
    molecule, = scan_fusion_bam(bam).molecules
    assert any(segment.secondary for segment in molecule.segments)
    assert find_read_paths(molecule.segments).status == "ambiguous"


def test_same_chromosome_gene_transition_is_nominated_without_soft_clipping(tmp_path):
    bam = tmp_path / "readthrough.bam"
    _write_bam(bam, [{"cigar": "100M800N100M"}])
    genes = GeneIndex([
        GeneModel("geneA", "chr1", "+", ((100, 200),)),
        GeneModel("geneB", "chr1", "+", ((1000, 1100),)),
    ])
    assert len(scan_fusion_bam(bam, genes).molecules) == 1
    assert scan_fusion_bam(bam).molecules == ()


def test_reverse_primary_full_sequence_is_normalized(tmp_path):
    bam = tmp_path / "reverse.bam"
    sequence = "A" * 100 + "C" * 100
    _write_bam(bam, [{"cigar": "100S100M", "flag": 16,
                      "sequence": reverse_complement(sequence),
                      "sa": "chr2,101,+,100S100M,60,0;"}])
    molecule, = scan_fusion_bam(bam).molecules
    assert molecule.sequence == sequence


def test_hard_clipped_only_molecule_requires_external_sequence(tmp_path):
    bam = tmp_path / "hard.bam"
    _write_bam(bam, [{"flag": 2048, "cigar": "100H100M100H"}])
    molecule, = scan_fusion_bam(bam).molecules
    assert molecule.sequence is None
    assert molecule.status == "needs_sequence"


def test_malformed_sa_invalidates_molecule_instead_of_dropping_alternative(tmp_path):
    bam = tmp_path / "bad.bam"
    _write_bam(bam, [{"sa": "chr2,101,+,100S100M100S,60,0;malformed;"}])
    molecule, = scan_fusion_bam(bam).molecules
    assert molecule.status == "unresolved"
    assert molecule.segments == ()


def test_duplicate_sa_and_bam_records_do_not_trip_distinct_segment_limit(tmp_path):
    bam = tmp_path / "copies.bam"
    _write_bam(bam, [
        {"sa": "chr2,101,+,100S100M100S,60,0;"},
        *[{"chrom": "chr2", "flag": 2048, "cigar": "100H100M100H"}] * 10,
    ])
    molecule, = scan_fusion_bam(bam, max_segments_per_read=2).molecules
    assert len(molecule.segments) == 2
    assert molecule.status == "ready"


def test_candidate_limit_does_not_publish_a_prefix(tmp_path):
    bam = tmp_path / "many.bam"
    # clips above the 250-bp nomination threshold, so both reads nominate
    _write_bam(bam, [{"name": "one", "cigar": "100M300S"}, {"name": "two", "cigar": "100M300S"}])
    with pytest.raises(FusionInputError, match="limit"):
        scan_fusion_bam(bam, max_candidate_reads=1)


def test_file_scan_failure_is_not_reported_as_no_fusions(tmp_path):
    bam = tmp_path / "invalid.bam"
    bam.write_bytes(b"not a BAM")
    with pytest.raises(FusionInputError, match="incomplete"):
        scan_fusion_bam(bam)


def test_changing_input_invalidates_scan(tmp_path):
    bam = tmp_path / "changing.bam"
    _write_bam(bam, [{}])
    with patch("fin.fusion.input._snapshot", side_effect=[("before",), ("after",)]):
        with pytest.raises(FusionInputError, match="changed"):
            scan_fusion_bam(bam)


def test_conflicting_molecule_sequences_are_unresolved(tmp_path):
    bam = tmp_path / "conflict.bam"
    _write_bam(bam, [
        {"cigar": "100M100S", "sequence": "A" * 200},
        {"cigar": "100S100M", "chrom": "chr2", "flag": 2048, "sequence": "C" * 200},
    ])
    molecule, = scan_fusion_bam(bam).molecules
    assert molecule.status == "unresolved"
    assert molecule.reason == "conflicting_forward_sequences"


def test_record_order_does_not_change_evidence(tmp_path):
    first = tmp_path / "first.bam"
    second = tmp_path / "second.bam"
    entries = [
        {"cigar": "100M200S", "sa": "chr2,101,+,100S100M100S,60,0;"},
        {"chrom": "chr2", "flag": 2048, "cigar": "100H100M100H"},
        {"chrom": "chr3", "flag": 2048, "cigar": "200H100M"},
    ]
    _write_bam(first, entries)
    _write_bam(second, list(reversed(entries)))
    assert scan_fusion_bam(first) == scan_fusion_bam(second)


def test_entering_an_overlapping_gene_is_not_a_gene_transition(tmp_path):
    """`labels()` attaches a gene on ONE base of exonic overlap, so on a real
    annotation a read that merely enters an overlapping gene's territory
    ({A} -> {A,B}) splits into arms. Nominating those made 6.1% of a real
    library candidates -- more than the cap -- for no fusion evidence."""
    bam = tmp_path / "overlap.bam"
    _write_bam(bam, [{"cigar": "100M800N100M"}])
    genes = GeneIndex([
        GeneModel("geneA", "chr1", "+", ((100, 200), (1000, 1100))),
        GeneModel("geneB", "chr1", "+", ((1000, 1100),)),
    ])
    assert scan_fusion_bam(bam, genes).molecules == ()


def test_a_real_gene_transition_is_still_nominated(tmp_path):
    bam = tmp_path / "readthrough.bam"
    _write_bam(bam, [{"cigar": "100M800N100M"}])
    genes = GeneIndex([
        GeneModel("geneA", "chr1", "+", ((100, 200),)),
        GeneModel("geneB", "chr1", "+", ((1000, 1100),)),
    ])
    molecule, = scan_fusion_bam(bam, genes).molecules
    assert scan_fusion_bam(bam, genes).record_counts["nominated_by_gene_transition"] == 1


def test_leaving_annotated_territory_is_not_a_gene_transition(tmp_path):
    # {A} -> {} is novel structure; a genuine distant fusion would arrive
    # with a supplementary alignment or a large clip, not on one collinear
    # segment
    bam = tmp_path / "novel.bam"
    _write_bam(bam, [{"cigar": "100M800N100M"}])
    genes = GeneIndex([GeneModel("geneA", "chr1", "+", ((100, 200),))])
    assert scan_fusion_bam(bam, genes).molecules == ()


def test_nomination_routes_are_counted_separately(tmp_path):
    bam = tmp_path / "routes.bam"
    _write_bam(bam, [{"name": "clipped", "cigar": "100M300S"},
                     {"name": "tagged", "cigar": "100M300S", "sa": "chr2,101,+,100S100M200S,60,0;"}])
    counts = scan_fusion_bam(bam).record_counts
    assert counts["nominated_by_supplementary_or_sa"] == 1
    assert counts["nominated_by_soft_clip"] == 1
    assert counts["nominated_molecules"] == 2
