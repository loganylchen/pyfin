"""Real minimap2 searches recover missing arms and preserve normal competitors."""

import random

import pytest

from fin.fusion.realignment import FusionGenomeAligner
from fin.fusion.records import segment_from_cigar
from fin.utils.sequences import reverse_complement


@pytest.fixture
def reference(tmp_path):
    generator = random.Random(591)
    genome = {chrom: "".join(generator.choices("ACGT", k=2000)) for chrom in "ABCD"}
    path = tmp_path / "genome.fa"
    path.write_text("".join(f">{chrom}\n{sequence}\n" for chrom, sequence in genome.items()))
    return path, genome


@pytest.mark.parametrize("count", [3, 4])
def test_real_genome_mapping_recovers_all_arms_without_sa_records(reference, count):
    path, genome = reference
    chromosomes = list(genome)[:count]
    query = "".join(genome[chrom][500:800] for chrom in chromosomes)
    result = FusionGenomeAligner(path).map_molecule("read", query)
    assert result.status == "complete"
    for index, chrom in enumerate(chromosomes):
        assert any(segment.chrom == chrom and segment.strand == "+"
                   and abs(segment.query_start - index * 300) <= 5
                   and abs(segment.query_end - (index + 1) * 300) <= 5
                   for segment in result.segments)
    assert (0, len(query)) in result.searched_windows


def test_reverse_middle_arm_uses_forward_molecule_coordinates(reference):
    path, genome = reference
    query = genome["A"][500:800] + reverse_complement(genome["B"][500:800]) + genome["C"][500:800]
    result = FusionGenomeAligner(path).map_molecule("read", query)
    assert result.status == "complete"
    assert any(segment.chrom == "B" and segment.strand == "-"
               and abs(segment.query_start - 300) <= 5 and abs(segment.query_end - 600) <= 5
               for segment in result.segments)


def test_full_length_normal_competitor_is_not_omitted(reference):
    path, genome = reference
    query = genome["A"][500:800] + genome["B"][500:800]
    with path.open("a") as handle:
        handle.write(f">normal\n{query}\n")
    existing = (segment_from_cigar("read", "A", 500, "+", "300M300S", 60),
                segment_from_cigar("read", "B", 500, "+", "300S300M", 60))
    result = FusionGenomeAligner(path).map_molecule("read", query, existing)
    assert result.status == "complete"
    assert any(segment.chrom == "normal" and segment.query_start == 0
               and segment.query_end == 600 for segment in result.segments)


def test_paralog_alternatives_survive_even_with_zero_mapq(reference):
    path, genome = reference
    with path.open("a") as handle:
        handle.write(f">A_copy\n{genome['A']}\n")
    result = FusionGenomeAligner(path).map_molecule("read", genome["A"][500:800])
    assert result.status == "complete"
    assert {segment.chrom for segment in result.segments} >= {"A", "A_copy"}
    assert any(segment.mapping_quality == 0 for segment in result.segments)


def test_hit_limit_withholds_all_partial_evidence(reference):
    path, genome = reference
    with path.open("a") as handle:
        handle.write("".join(f">copy{index}\n{genome['A']}\n" for index in range(8)))
    result = FusionGenomeAligner(path, max_hits=3).map_molecule("read", genome["A"][500:800])
    assert result.status == "unresolved"
    assert result.reason == "hit_limit"
    assert result.segments == ()


def test_window_limit_is_not_a_successful_partial_search(reference):
    path, genome = reference
    query = genome["A"][500:800] + genome["B"][500:800]
    result = FusionGenomeAligner(path, max_windows=1).map_molecule("read", query)
    assert result.status == "unresolved"
    assert result.reason == "window_limit"
    assert result.segments == ()


def test_sequence_missing_abstains(reference):
    path, _ = reference
    result = FusionGenomeAligner(path).map_molecule("read", "")
    assert result.reason == "missing_sequence"


def test_existing_records_must_belong_to_same_molecule(reference):
    path, genome = reference
    existing = (segment_from_cigar("other", "A", 500, "+", "300M", 60),)
    with pytest.raises(ValueError, match="disagree"):
        FusionGenomeAligner(path).map_molecule("read", genome["A"][500:800], existing)


def test_existing_clipped_tail_is_searched_even_if_whole_read_has_no_hits(reference):
    path, genome = reference
    query = genome["A"][500:800] + genome["B"][500:800]
    existing = (segment_from_cigar("read", "A", 500, "+", "300M300S", 60),)
    aligner = FusionGenomeAligner(path)
    backend = aligner.aligner

    class WholeReadAbstention:
        def map(self, sequence):
            return iter(()) if len(sequence) == 600 else backend.map(sequence)

    aligner.aligner = WholeReadAbstention()
    result = aligner.map_molecule("read", query, existing)
    assert result.status == "complete"
    assert (300, 600) in result.searched_windows
    assert any(segment.chrom == "B" and segment.query_start == 300 for segment in result.segments)
