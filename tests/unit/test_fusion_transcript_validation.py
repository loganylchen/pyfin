"""Full stitched validation cannot substitute separate local edge matches."""

from dataclasses import replace
import random

from fin.fusion.annotation import AnnotatedReadPath, FusionArm
from fin.fusion.records import AlignmentBlock
from fin.fusion.validation import stitched_template, validate_stitched_path


def _fixture():
    generator = random.Random(193)
    genome = {chrom: "".join(generator.choices("ACGT", k=1000)) for chrom in "ABC"}
    arms = tuple(FusionArm(chrom, "+", (AlignmentBlock(100, 400, index * 300,
                                                     (index + 1) * 300, 300),), (chrom,))
                 for index, chrom in enumerate(genome))
    path = AnnotatedReadPath("read", arms, "interchromosomal_candidate", "resolved")
    query, boundaries = stitched_template(path, genome)
    return genome, path, query, boundaries


def test_whole_three_gene_transcript_spans_both_edges():
    genome, path, query, boundaries = _fixture()
    assert boundaries == (300, 600)
    result = validate_stitched_path(path, query, genome)
    assert result.status == "validated"
    assert result.query_coverage == result.identity == 1
    assert result.score == 1800
    # anchors now count aligned template positions within the 24-bp window on
    # each side of a boundary (a nanopore read with a 2-bp deletion inside a
    # 12-bp flank is no longer rejected); an error-free read fills the window
    assert result.junction_anchors == ((24, 24), (24, 24))


def test_partial_read_cannot_validate_complete_chain():
    genome, path, query, _ = _fixture()
    result = validate_stitched_path(path, query[:600], genome)
    assert result.status != "validated"


def test_internal_adapter_insertion_is_not_hidden_by_full_length_alignment():
    genome, path, query, _ = _fixture()
    result = validate_stitched_path(path, query[:300] + "A" * 40 + query[300:], genome)
    assert result.status != "validated"


def test_wrong_middle_partner_does_not_pass_outer_arm_matches():
    genome, path, query, _ = _fixture()
    middle = replace(path.arms[1], blocks=(AlignmentBlock(600, 900, 300, 600, 300),))
    wrong = replace(path, arms=(path.arms[0], middle, path.arms[2]))
    assert validate_stitched_path(wrong, query, genome).status != "validated"


def test_invalid_reference_sequence_abstains():
    genome, path, query, _ = _fixture()
    genome["B"] = "N" * 1000
    assert validate_stitched_path(path, query, genome).reason == "ambiguous_sequence"
