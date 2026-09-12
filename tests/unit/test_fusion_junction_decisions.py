"""Coordinate correction requires anchored contrasts, not just splice motifs."""

from dataclasses import replace
import random

import pytest

from fin.fusion.annotation import FusionArm
from fin.fusion.junctions import arm_sequence, evaluate_junction, move_junction
from fin.fusion.records import AlignmentBlock


def _fixture(left_strand="+", right_strand="+", internal=False, seed=732):
    generator = random.Random(seed)
    genome = {chrom: "".join(generator.choices("ACGT", k=600))
              for chrom in ("left", "right")}
    left = FusionArm("left", left_strand, (AlignmentBlock(100, 180, 0, 80, 80),), ("A",))
    right = FusionArm("left" if internal else "right", right_strand,
                      (AlignmentBlock(300, 380, 80, 160, 80),), ("A" if internal else "B",))
    if internal and left_strand == "-":
        left = replace(left, blocks=(AlignmentBlock(300, 380, 0, 80, 80),))
        right = replace(right, blocks=(AlignmentBlock(100, 180, 80, 160, 80),))
    query = arm_sequence(left, genome) + arm_sequence(right, genome)
    return genome, left, right, query


@pytest.mark.parametrize("left_strand,right_strand", [("+", "+"), ("+", "-"), ("-", "+"), ("-", "-")])
def test_read_supported_shifted_fusion_edge_is_corrected(left_strand, right_strand):
    genome, left, right, query = _fixture(left_strand, right_strand)
    truth = left.exit_position, right.entry_position
    shifted = (truth[0] + (2 if left_strand == "+" else -2), truth[1])
    observed_left, observed_right = move_junction(left, right, shifted)
    decision = evaluate_junction(query, observed_left, observed_right, genome,
                                 supported_positions=(truth,))
    assert decision.status == "correction_proposed"
    assert decision.selected_positions == truth
    assert decision.signal_status == "not_scored"
    assert decision.query_start == 45
    assert decision.query_end == 115


@pytest.mark.parametrize("strand", ["+", "-"])
def test_internal_intron_shift_is_corrected(strand):
    genome, left, right, query = _fixture(strand, strand, internal=True)
    truth = left.exit_position, right.entry_position
    observed = move_junction(left, right, (truth[0] + (2 if strand == "+" else -2), truth[1]))
    decision = evaluate_junction(query, *observed, genome, supported_positions=(truth,),
                                 internal_intron=True)
    assert decision.status == "correction_proposed"
    assert decision.selected_positions == truth


def test_annotation_nomination_cannot_override_exact_read_sequence():
    genome, left, right, query = _fixture()
    truth = left.exit_position, right.entry_position
    decision = evaluate_junction(query, left, right, genome,
                                 annotation_positions=((truth[0] + 2, truth[1]),))
    assert decision.status == "retained"
    assert decision.selected_positions == truth


def test_identical_microhomology_templates_preserve_coordinate_ambiguity():
    genome, left, right, _ = _fixture()
    genome = {chrom: "A" * 600 for chrom in genome}
    decision = evaluate_junction("A" * 160, left, right, genome,
                                 supported_positions=((181, 301),))
    assert decision.status == "unresolved"
    assert decision.reason == "indistinguishable_coordinates"
    assert decision.selected_positions is None


def test_canonical_microhomology_tie_is_not_forced_to_supported_coordinate():
    genome, left, right, query = _fixture("-", "-", seed=731)
    truth = left.exit_position, right.entry_position
    observed = move_junction(left, right, (truth[0] - 2, truth[1]))
    decision = evaluate_junction(query, *observed, genome, supported_positions=(truth,))
    assert decision.status == "unresolved"
    assert decision.reason == "indistinguishable_coordinates"
    winners = [item.positions for item in decision.contrasts if item.alignment.score == 140]
    assert winners == [(99, 379), truth]


def test_candidate_limit_cannot_turn_truncated_search_into_unique_call():
    genome, left, right, query = _fixture()
    decision = evaluate_junction(query, left, right, genome,
                                 supported_positions=((182, 300),), max_candidates=1)
    assert decision.status == "unresolved"
    assert decision.reason == "candidate_limit"
    assert not decision.contrasts


def test_missing_reference_abstains_without_editing_observation():
    genome, left, right, query = _fixture()
    del genome["right"]
    decision = evaluate_junction(query, left, right, genome)
    assert decision.status == "unresolved"
    assert decision.reason == "invalid_observed_template"
    assert decision.observed_positions == (180, 300)
