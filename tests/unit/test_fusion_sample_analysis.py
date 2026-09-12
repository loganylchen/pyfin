"""Real reference searches and per-molecule adjudication before support counts."""

from dataclasses import replace
import random
import json

import pytest

from fin.fusion.analysis import analyze_fusion_input, FusionSettings
from fin.fusion.annotation import GeneIndex, GeneModel
from fin.fusion.input import FusionInput, MoleculeEvidence
from fin.fusion.realignment import FusionGenomeAligner
from fin.fusion.records import segment_from_cigar
from fin.fusion.output import write_fusion_outputs


@pytest.fixture
def fixture(tmp_path):
    generator = random.Random(591)
    genome = {chrom: "".join(generator.choices("ACGT", k=2000)) for chrom in "ABCD"}
    reference = tmp_path / "genome.fa"
    reference.write_text("".join(f">{chrom}\n{sequence}\n" for chrom, sequence in genome.items()))
    genes = GeneIndex([GeneModel(chrom, chrom, "+", ((400, 900),)) for chrom in genome])
    return reference, genome, genes


def _input(genome, paths):
    molecules = []
    for read_id, chromosomes in enumerate(paths):
        name = f"read{read_id}"
        sequence = "".join(genome[chrom][500:800] for chrom in chromosomes)
        cigar = "300M" + (f"{len(sequence) - 300}S" if len(sequence) > 300 else "")
        segment = segment_from_cigar(name, chromosomes[0], 500, "+", cigar, 60)
        molecules.append(MoleculeEvidence(name, (segment,), sequence, "ready", "test"))
    return FusionInput(tuple(molecules), {}, True)


def test_three_gene_chain_has_complete_read_support(fixture):
    reference, genome, genes = fixture
    result = analyze_fusion_input(_input(genome, ["ABC", "ABC"]), genome,
                                   FusionGenomeAligner(reference), genes)
    assert len(result.chains) == 1
    chain, = result.chains
    assert chain.representative.gene_path == (("A",), ("B",), ("C",))
    assert chain.full_support_read_ids == ("read0", "read1")
    assert len(chain.edge_support_read_ids) == 2
    assert chain.status.startswith("supported_")
    assert all(read.status == "supported_path" for read in result.reads)


def test_independent_pairs_never_manufacture_a_three_gene_fusion(fixture):
    reference, genome, genes = fixture
    result = analyze_fusion_input(_input(genome, ["AB", "AB", "BC", "BC"]), genome,
                                   FusionGenomeAligner(reference), genes)
    assert result.chains
    assert all(len(chain.representative.arms) == 2 for chain in result.chains)


def test_full_length_normal_explanation_prevents_fusion_call(fixture):
    reference, genome, genes = fixture
    genome["normal"] = genome["A"][500:800] + genome["B"][500:800]
    with reference.open("a") as handle:
        handle.write(f">normal\n{genome['normal']}\n")
    result = analyze_fusion_input(_input(genome, ["AB", "AB"]), genome,
                                   FusionGenomeAligner(reference), genes)
    assert result.chains == ()
    assert all(read.reason == "normal_alignment_competitor" for read in result.reads)


def test_paralog_alternative_cannot_be_hidden_by_original_bam_mapq(fixture):
    reference, genome, genes = fixture
    genome["A_copy"] = genome["A"]
    with reference.open("a") as handle:
        handle.write(f">A_copy\n{genome['A_copy']}\n")
    result = analyze_fusion_input(_input(genome, ["AB", "AB"]), genome,
                                   FusionGenomeAligner(reference), genes)
    assert result.chains == ()
    assert all(read.status == "unresolved" for read in result.reads)


def test_single_read_cannot_reach_two_read_support(fixture):
    reference, genome, genes = fixture
    result = analyze_fusion_input(_input(genome, ["AB"]), genome,
                                   FusionGenomeAligner(reference), genes)
    assert all(chain.status == "insufficient_full_chain_support" for chain in result.chains)
    assert len(result.chains) == 1


def test_invalid_settings_and_incomplete_input_fail_before_mapping(fixture):
    reference, genome, genes = fixture
    aligner = FusionGenomeAligner(reference)
    evidence = _input(genome, ["AB"])
    with pytest.raises(ValueError, match="incomplete"):
        analyze_fusion_input(replace(evidence, complete=False), genome, aligner, genes)
    with pytest.raises(ValueError, match="settings"):
        analyze_fusion_input(evidence, genome, aligner, genes, FusionSettings(min_score_margin=float("nan")))


def test_failed_local_corrections_remain_in_published_audit(fixture, tmp_path):
    reference, genome, genes = fixture
    evidence = _input(genome, ["AB"])
    original = evidence.molecules[0]
    altered = replace(original, sequence=original.sequence[:290] + "N" * 20 + original.sequence[310:])
    result = analyze_fusion_input(replace(evidence, molecules=(altered,)), genome,
                                   FusionGenomeAligner(reference), genes)
    read, = result.reads
    assert read.status == "unresolved"
    assert read.corrections
    assert any(correction.proposed is None for correction in read.corrections)
    output = tmp_path / "audit"
    write_fusion_outputs(result, genome, output)
    corrections = [json.loads(line) for line in (output / "junction_corrections.jsonl").read_text().splitlines()]
    assert corrections
    assert any(correction["proposed"] is None for correction in corrections)
    read_record = json.loads((output / "fusion_reads.jsonl").read_text())
    assert set(read_record["correction_ids"]) == {correction["correction_id"] for correction in corrections}
