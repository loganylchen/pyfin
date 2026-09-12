"""Whole-path correction covers internal introns and every multi-gene edge."""

from dataclasses import replace
import random
from types import SimpleNamespace

from fin.fusion.annotation import AnnotatedReadPath, FusionArm
from fin.fusion.correction import annotation_junction_registry, correct_read_path, observed_junction_registry
from fin.fusion.junctions import arm_sequence, move_junction
from fin.fusion.records import AlignmentBlock


def _fixture():
    generator = random.Random(846)
    genome = {chrom: "".join(generator.choices("ACGT", k=800))
              for chrom in ("A", "B", "C", "D")}
    arms = []
    for index, chrom in enumerate(genome):
        query_start = index * 100
        blocks = ((AlignmentBlock(100, 150, query_start, query_start + 50, 50),
                   AlignmentBlock(300, 350, query_start + 50, query_start + 100, 50))
                  if index == 1 else
                  (AlignmentBlock(100, 200, query_start, query_start + 100, 100),))
        arms.append(FusionArm(chrom, "+", blocks, (chrom,)))
    path = AnnotatedReadPath("read", tuple(arms), "interchromosomal_candidate", "resolved")
    sequence = "".join(arm_sequence(arm, genome) for arm in path.arms)
    return genome, path, sequence


def test_four_gene_chain_validates_all_edges_and_internal_intron():
    genome, path, sequence = _fixture()
    result = correct_read_path(path, sequence, genome, tolerance=0)
    assert result.status == "retained"
    assert result.proposed == path
    assert [item.kind for item in result.junctions] == ["fusion", "intron", "fusion", "fusion"]
    assert len(result.validation) == 4
    assert all(item.decision.status == "retained" for item in result.validation)


def test_joint_internal_and_fusion_correction_preserves_original_evidence():
    genome, truth, sequence = _fixture()
    arms = list(truth.arms)
    arms[0], arms[1] = move_junction(arms[0], arms[1], (202, 100))
    middle = arms[1]
    arms[1] = replace(middle, blocks=(replace(middle.blocks[0], reference_end=152),
                                    middle.blocks[1]))
    observed = replace(truth, arms=tuple(arms))
    result = correct_read_path(observed, sequence, genome,
                               supported_registry=observed_junction_registry([truth]))
    assert result.status == "correction_proposed"
    assert result.proposed == truth
    assert result.observed == observed
    assert result.observed.arms[0].exit_position == 202
    assert sum(item.decision.status == "correction_proposed" for item in result.junctions) == 2
    assert all(item.decision.status == "retained" for item in result.validation)


def test_failed_edge_withholds_whole_transaction_not_just_that_edge():
    genome, truth, sequence = _fixture()
    arms = list(truth.arms)
    arms[0], arms[1] = move_junction(arms[0], arms[1], (202, 100))
    observed = replace(truth, arms=tuple(arms))
    del genome["D"]
    result = correct_read_path(observed, sequence, genome,
                               supported_registry=observed_junction_registry([truth]))
    assert result.status == "unresolved"
    assert result.reason == "unresolved_junction"
    assert result.proposed is None
    assert result.observed == observed
    assert any(item.decision.status == "correction_proposed" for item in result.junctions)


def test_observed_registry_is_order_and_duplicate_invariant():
    _, path, _ = _fixture()
    alternate = replace(path, read_id="other")
    assert observed_junction_registry([path, alternate, path]) == observed_junction_registry([alternate, path])


def test_missing_sequence_is_audited_not_replaced_by_reference():
    genome, path, _ = _fixture()
    result = correct_read_path(path, "", genome)
    assert result.status == "unresolved"
    assert result.proposed is None
    assert all(item.decision.reason == "invalid_query_window" for item in result.junctions)


def test_annotation_nominations_are_connected_to_read_scored_correction():
    genome, truth, sequence = _fixture()
    arms = list(truth.arms)
    arms[0], arms[1] = move_junction(arms[0], arms[1], (202, 100))
    observed = replace(truth, arms=tuple(arms))
    transcripts = [SimpleNamespace(chrom=arm.chrom, strand=arm.strand, gene_id=arm.gene_ids[0],
                                   exons=arm.exons) for arm in truth.arms]
    registry = annotation_junction_registry([observed], transcripts)
    result = correct_read_path(observed, sequence, genome, annotation_registry=registry)
    assert result.status == "correction_proposed"
    assert result.proposed == truth
    assert any("annotation" in contrast.origin for junction in result.junctions
               for contrast in junction.decision.contrasts)


def test_reverse_annotation_intron_uses_transcript_order():
    _, original, _ = _fixture()
    arm = FusionArm("B", "-", (AlignmentBlock(300, 350, 0, 50, 50),
                               AlignmentBlock(100, 150, 50, 100, 50)), ("B",))
    path = replace(original, arms=(arm,))
    transcript = SimpleNamespace(chrom="B", strand="-", gene_id="B", exons=((100, 150), (300, 350)))
    registry = annotation_junction_registry([path], [transcript])
    assert registry[("intron", "B", "-", ("B",), "B", "-", ("B",))] == ((300, 150),)
