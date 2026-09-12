"""Full-query coordinates, molecule identity and multi-arm path contracts."""

from dataclasses import replace
from itertools import permutations

import pysam
import pytest

from fin.fusion.paths import ReadPath, find_read_paths
from fin.fusion.records import (
    deduplicate_segments,
    parse_cigar,
    parse_sa_tag,
    segment_from_cigar,
    segments_from_record,
)


def _segment(chrom, cigar, *, strand="+", read_id="molecule", quality=60, start=100):
    return segment_from_cigar(read_id, chrom, start, strand, cigar, quality)


def _record(chrom, cigar, *, flag=0, sa=None, read_id="molecule"):
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6"},
        "SQ": [{"SN": name, "LN": 10000} for name in ("chr1", "chr2", "chr3", "chr4")],
    })
    record = pysam.AlignedSegment(header)
    record.query_name = read_id
    record.reference_id = header.get_tid(chrom)
    record.reference_start = 100
    record.flag = flag
    record.mapping_quality = 60
    record.cigarstring = cigar
    record.query_sequence = "A" * sum(length for code, length in parse_cigar(cigar)
                                       if code in (0, 1, 4, 7, 8))
    record.set_tag("NM", 0)
    if sa is not None:
        record.set_tag("SA", sa)
    return record


@pytest.mark.parametrize("cigar", ["", "*", "0M", "10Mgarbage", "10M5S10M", "5S5H10M", "10I"])
def test_malformed_or_unaligned_cigar_is_rejected(cigar):
    with pytest.raises(ValueError):
        parse_cigar(cigar)


def test_forward_cigar_counts_hard_clips_on_full_molecule_axis():
    segment = _segment("chr1", "20H10S40M100N30M5S15H")
    assert segment.read_length == 120
    assert (segment.query_start, segment.query_end) == (30, 100)
    assert segment.exons == ((100, 140), (240, 270))
    assert segment.entry_position == 100
    assert segment.exit_position == 270


def test_reverse_cigar_reverses_blocks_and_full_query_coordinates():
    segment = _segment("chr1", "20H10S40M100N30M5S15H", strand="-")
    assert segment.read_length == 120
    assert (segment.query_start, segment.query_end) == (20, 90)
    assert segment.exons == ((240, 270), (100, 140))
    assert [(block.query_start, block.query_end) for block in segment.blocks] == [(20, 50), (50, 90)]
    assert segment.entry_position == 270
    assert segment.exit_position == 100


def test_insertions_and_deletions_do_not_create_fusion_segments():
    segment = _segment("chr1", "10S30M3I20M4D10M20S")
    assert segment.read_length == 93
    assert len(segment.blocks) == 1
    assert segment.exons == ((100, 164),)
    assert segment.aligned_bases == 60
    assert (segment.query_start, segment.query_end) == (10, 73)


def test_sa_coordinates_are_one_based_and_include_all_entries():
    segments = parse_sa_tag("molecule", "chr2,101,+,100S100M100S,60,2;chr3,201,-,100M200H,40,3;")
    assert len(segments) == 2
    assert segments[0].exons == ((100, 200),)
    assert (segments[1].query_start, segments[1].query_end) == (200, 300)
    assert segments[1].exons == ((200, 300),)


@pytest.mark.parametrize("tag", [
    "chr1,0,+,10M,60,0;",
    "chr1,1,?,10M,60,0;",
    "chr1,1,+,10M,256,0;",
    "chr1,1,+,10M,60,-1;",
    "chr1,1,+,10M,60;",
    "chr1,1,+,10M,60,0;;",
    "chr1,1,+,10M,60,0;bad_entry;",
])
def test_sa_parsing_does_not_silently_ignore_invalid_alternatives(tag):
    with pytest.raises(ValueError):
        parse_sa_tag("molecule", tag)


def test_sa_bam_duplicate_does_not_inflate_segments():
    primary = _record("chr1", "100M200S", sa="chr2,101,+,100S100M100S,60,0;")
    supplementary = _record("chr2", "100H100M100H", flag=2048)
    segments = deduplicate_segments(segments_from_record(primary) + segments_from_record(supplementary))
    assert len(segments) == 2
    assert next(segment for segment in segments if segment.chrom == "chr2").sources == (
        "bam_supplementary", "sa",
    )


def test_sa_length_mismatch_invalidates_record():
    record = _record("chr1", "100M200S", sa="chr2,101,+,100S100M,60,0;")
    with pytest.raises(ValueError, match="length"):
        segments_from_record(record)


def test_middle_primary_can_resolve_three_ordered_genomic_arms():
    primary = _record("chr2", "100S100M100S", sa=(
        "chr1,101,+,100M200S,60,0;chr3,101,+,200S100M,60,0;"
    ))
    result = find_read_paths(segments_from_record(primary))
    assert result.status == "unique"
    path = result.paths[0]
    assert [segment.chrom for segment in path.segments] == ["chr1", "chr2", "chr3"]
    assert len(path.junctions) == 2
    assert path.query_covered_bases == 300
    assert path.coverage_fraction == 1.0


def test_three_arm_path_is_invariant_to_primary_choice():
    descriptions = (
        ("chr1", "100M200S"),
        ("chr2", "100S100M100S"),
        ("chr3", "200S100M"),
    )
    signatures = []
    for primary_index in range(3):
        segments = []
        for index, (chrom, cigar) in enumerate(descriptions):
            segments.extend(segments_from_record(_record(chrom, cigar, flag=0 if index == primary_index else 2048)))
        result = find_read_paths(segments)
        signatures.append(result.paths[0].geometry_key)
    assert signatures[0] == signatures[1] == signatures[2]


def test_four_arm_path_does_not_flatten_to_one_pair():
    segments = [
        _segment("chr1", "100M300S"),
        _segment("chr2", "100S100M200S"),
        _segment("chr3", "200S100M100S"),
        _segment("chr4", "300S100M"),
    ]
    result = find_read_paths(segments)
    assert result.status == "unique"
    assert len(result.paths[0].segments) == 4
    assert len(result.paths[0].junctions) == 3


def test_reverse_middle_arm_keeps_forward_query_order():
    segments = [
        _segment("chr1", "100M200S"),
        _segment("chr2", "100H100M100H", strand="-"),
        _segment("chr3", "200S100M"),
    ]
    result = find_read_paths(segments)
    assert result.status == "unique"
    assert result.paths[0].junctions == (
        (("chr1", 200, "+"), ("chr2", 200, "-")),
        (("chr2", 100, "-"), ("chr3", 100, "+")),
    )


def test_different_molecules_cannot_be_joined_into_a_longer_path():
    segments = [
        _segment("chr1", "100M100S", read_id="ab"),
        _segment("chr2", "100S100M", read_id="ab"),
        _segment("chr2", "100M100S", read_id="bc"),
        _segment("chr3", "100S100M", read_id="bc"),
    ]
    with pytest.raises(ValueError, match="different molecules"):
        find_read_paths(segments)


def test_competing_partner_paths_are_not_broken_by_record_order():
    segments = [
        _segment("chr1", "100M100S"),
        _segment("chr2", "100S100M"),
        _segment("chr3", "100S100M"),
    ]
    signatures = set()
    for order in permutations(segments):
        result = find_read_paths(order)
        assert result.status == "ambiguous"
        signatures.add(tuple(path.geometry_key for path in result.paths))
    assert len(signatures) == 1


def test_normal_whole_read_alignment_remains_a_competing_explanation():
    result = find_read_paths([
        _segment("chr1", "100M100S"),
        _segment("chr2", "100S100M"),
        _segment("chr3", "200M"),
    ])
    assert result.status == "ambiguous"
    assert sorted(len(path.segments) for path in result.paths) == [1, 2]


def test_unknown_mapq_is_not_treated_as_extreme_confidence():
    result = find_read_paths([_segment("chr1", "100M100S"), _segment("chr2", "100S100M", quality=255)])
    assert result.status == "unresolved"
    assert result.reason == "insufficient_query_coverage"


def test_large_internal_gap_does_not_form_a_complete_path():
    result = find_read_paths([_segment("chr1", "100M140S"), _segment("chr2", "140S100M")])
    assert result.status == "unresolved"
    assert result.paths == ()


def test_excessive_overlap_does_not_form_a_complete_path():
    result = find_read_paths([_segment("chr1", "120M80S"), _segment("chr2", "80S120M")])
    assert result.status == "unresolved"


def test_microhomology_overlap_does_not_count_query_bases_twice():
    result = find_read_paths([_segment("chr1", "105M95S"), _segment("chr2", "100S100M")])
    assert result.status == "unique"
    assert result.paths[0].query_covered_bases == 200


def test_search_limits_cannot_create_a_false_unique_path():
    result = find_read_paths([
        _segment("chr1", "100M100S"),
        _segment("chr2", "100S100M"),
        _segment("chr3", "100S100M"),
    ], max_paths=1)
    assert result.status == "unresolved"
    assert result.reason == "path_limit"
    assert result.paths == ()


def test_segment_limit_is_applied_after_deduplication():
    segment = _segment("chr1", "100M")
    result = find_read_paths([segment] * 100, max_segments=1)
    assert result.status == "unique"
    assert result.segment_count == 1


def test_quality_conflicts_between_duplicate_records_are_conservative():
    segment = _segment("chr1", "100M")
    unique = deduplicate_segments([segment, replace(segment, mapping_quality=0, sources=("sa",))])
    assert unique[0].mapping_quality == 0


def test_full_read_length_conflict_is_visible_abstention():
    result = find_read_paths([_segment("chr1", "100M100S"), _segment("chr2", "100S120M")])
    assert result.status == "unresolved"
    assert result.reason == "inconsistent_read_length"


def test_read_path_cannot_reverse_query_order():
    with pytest.raises(ValueError, match="advance"):
        ReadPath("molecule", (_segment("chr2", "100S100M"), _segment("chr1", "100M100S")))
