"""Full sequence recovery must not erase alignment or duplicate-record errors."""

from dataclasses import replace

import pytest

from fin.fusion.input import FusionInput, FusionInputError, MoleculeEvidence, supply_fusion_sequences
from fin.fusion.records import segment_from_cigar


def _evidence(sequence=None):
    segment = segment_from_cigar("read", "A", 100, "+", "100M200H", 60)
    molecule = MoleculeEvidence("read", (segment,), sequence,
                                "ready" if sequence else "needs_sequence", "initial")
    return FusionInput((molecule,), {"nominated_molecules": 1}, True)


def test_hard_clipped_full_sequence_is_recovered(tmp_path):
    path = tmp_path / "reads.fastq"
    path.write_text("@read description\n" + "ACG" * 100 + "\n+\n" + "I" * 300 + "\n")
    result = supply_fusion_sequences(_evidence(), path)
    assert result.molecules[0].status == "ready"
    assert result.molecules[0].sequence == "ACG" * 100
    assert result.record_counts["sequence_recovered_molecules"] == 1
    assert result.record_counts["molecules_missing_sequence"] == 0


@pytest.mark.parametrize("sequences,reason", [
    (["A" * 299], "sequence_length_mismatch"),
    (["A" * 300, "C" * 300], "conflicting_sequence_records"),
])
def test_invalid_full_sequences_remain_unresolved(tmp_path, sequences, reason):
    path = tmp_path / "reads.fa"
    path.write_text("".join(f">read\n{sequence}\n" for sequence in sequences))
    result = supply_fusion_sequences(_evidence(), path)
    assert result.molecules[0].status == "unresolved"
    assert result.molecules[0].reason == reason
    assert result.molecules[0].sequence is None


def test_bam_fastq_conflict_cannot_silently_overwrite_bam(tmp_path):
    path = tmp_path / "reads.fa"
    path.write_text(">read\n" + "C" * 300 + "\n")
    result = supply_fusion_sequences(_evidence("A" * 300), path)
    assert result.molecules[0].reason == "bam_sequence_conflict"


def test_matching_duplicates_supply_only_one_molecule(tmp_path):
    path = tmp_path / "reads.fa"
    path.write_text((">read\n" + "A" * 300 + "\n") * 2)
    result = supply_fusion_sequences(_evidence(), path)
    assert result.record_counts["sequence_validated_molecules"] == 1
    assert len(result.molecules) == 1


def test_missing_id_is_not_filled_with_another_molecule(tmp_path):
    path = tmp_path / "reads.fa"
    path.write_text(">other\n" + "A" * 300 + "\n")
    result = supply_fusion_sequences(_evidence(), path)
    assert result.molecules[0].status == "needs_sequence"


def test_sequence_cannot_rescue_invalid_alignment_evidence(tmp_path):
    path = tmp_path / "reads.fa"
    path.write_text(">read\n" + "A" * 300 + "\n")
    original = _evidence()
    invalid = replace(original, molecules=(replace(original.molecules[0], status="unresolved",
                                                   reason="malformed_sa"),))
    result = supply_fusion_sequences(invalid, path)
    assert result.molecules[0].reason == "malformed_sa"


def test_incomplete_sequence_scan_fails_the_analysis(tmp_path):
    path = tmp_path / "reads.fastq"
    path.write_text("@read\nAAAA\n+\nI\n")
    with pytest.raises(FusionInputError, match="incomplete"):
        supply_fusion_sequences(_evidence(), path)
