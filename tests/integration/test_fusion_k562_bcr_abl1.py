"""The real positive control: BCR-ABL1 in K562 direct RNA.

`tests/data/fusion_k562_bcr_abl1/` freezes the 95 real records overlapping the
BCR and ABL1 loci of SG-NEx K562 rep4 run1 -- unmodified sequences, SA tags
(including the summarised `175M723D` CIGAR), hard clips, indels, MAPQ -- on
two 300-kb reference slices with every BCR/ABL1 GENCODE v44 transcript kept
whole. Before the real-data fixes the analysis found NOTHING here at the
defaults, and at a relaxed coverage found one wrong, unlabelled chain. The
expectations below are what the reads actually show, fixed after inspection
of each molecule: K562 carries the e14a2 (b3a2) junction, BCR exon 14 joined
to ABL1 exon 2, on transcripts ENST00000305877/ENST00000359540 (BCR) and
ENST00000318560/ENST00000372348 (ABL1).

This is a two-locus reference: it proves the regression, not genome-wide
uniqueness of the mapping (the full-reference run is recorded separately).
"""
import csv
import json
from pathlib import Path

import pytest

from fin.fusion.analysis import FusionSettings
from fin.fusion.pipeline import run_fusion_analysis

FIXTURE = Path(__file__).resolve().parents[1] / "data" / "fusion_k562_bcr_abl1"

# the SA-tagged BCR-ABL1 molecules of this library slice
BCR_ABL1_MOLECULES = {
    "35fc0c7f-7c48-481a-8120-425a4c679269", "4736bf19-b8bd-4e49-8145-1e1201614e4b",
    "65bf7bff-aa90-4fe7-ad10-8760057ff722", "a2d444fe-834b-4595-9253-3a306c35aac2",
    "bff45671-0f39-44c9-ae2a-623499c946f5", "938666a7-48c4-4ddc-8088-9921d441f8ce",
    "97ea8286-b546-4318-b314-86e01ff971fb",
}


@pytest.fixture(scope="module")
def published(tmp_path_factory):
    out = tmp_path_factory.mktemp("k562") / "result"
    run_fusion_analysis(FIXTURE / "reads.bam", FIXTURE / "reference.fa", out,
                        gtf_path=FIXTURE / "annotation.gtf", sequence_path=FIXTURE / "reads.fastq",
                        settings=FusionSettings(), threads=2)
    return out / "fusion"


def _rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_the_canonical_e14a2_fused_transcript_is_the_one_model(published):
    models = _rows(published / "fusion_models.tsv")
    assert len(models) == 1, [m["junction_names"] for m in models]
    model = models[0]
    assert json.loads(model["gene_names"]) == [["BCR"], ["ABL1"]]
    assert json.loads(model["junction_names"]) == ["BCR:e14>ABL1:e2"]
    assert json.loads(model["arm_match_classes"]) == ["arm_fsm", "arm_fsm"]
    tx = json.loads(model["compatible_transcript_ids"])
    assert {t.split(".")[0] for t in tx[0]} == {"ENST00000305877", "ENST00000359540"}
    assert {t.split(".")[0] for t in tx[1]} == {"ENST00000318560", "ENST00000372348"}
    # four of the seven SA-tagged molecules are exact members (each holds the
    # group's alternative placements itself); a fifth settled decisively on a
    # 2-bp-shifted internal donor and is REPORTED as nearby, not counted --
    # coordinate proximity is not evidence
    unique = set(json.loads(model["unique_read_ids"]))
    assert len(unique) >= 4
    assert unique <= BCR_ABL1_MOLECULES
    assert model["support_compatible"] == "0"
    nearby = set(json.loads(model["nearby_read_ids"]))
    assert nearby and nearby <= BCR_ABL1_MOLECULES and not (nearby & unique)
    assert model["envelope_status"] == "single_molecule_spanned"
    # the 3' end reaches ABL1's annotated end (direct RNA is 3'-anchored); the
    # 5' end is whatever the least-degraded molecule reached
    assert model["three_prime_status"] == "annotated_boundary_observed"
    assert int(model["observed_length"]) > 7000


def test_the_fused_transcript_sequence_is_the_stitched_reference_exons(published):
    model = _rows(published / "fusion_models.tsv")[0]
    header, sequence = None, ""
    for line in (published / "fusion_models.fa").read_text().splitlines():
        if line.startswith(">"):
            header = line
        else:
            sequence += line
    assert header.startswith(">" + model["sequence_id"] + " ")
    assert "junction=BCR:e14>ABL1:e2" in header and "kind=evidence_constrained" in header
    assert len(sequence) == int(model["observed_length"])
    # the exon table reconstructs the FASTA from the reference exactly
    import pysam

    fasta = pysam.FastaFile(str(FIXTURE / "reference.fa"))
    exons = [e for e in _rows(published / "fusion_model_exons.tsv") if e["sequence_id"] == model["sequence_id"]]
    exons.sort(key=lambda e: int(e["tx_start0"]))
    rebuilt = "".join(fasta.fetch(e["chrom"], int(e["start0"]), int(e["end0"])).upper() for e in exons)
    assert rebuilt == sequence
    # BCR exon 14 is the last chr22 exon, ABL1 exon 2 the first chr9 exon
    chr22 = [e for e in exons if e["chrom"] == "chr22"]
    chr9 = [e for e in exons if e["chrom"] == "chr9"]
    assert int(chr22[-1]["end0"]) == 23290413 - 23100000
    assert int(chr9[0]["start0"]) == 130854063 - 130600000


def test_the_summarised_sa_geometry_was_superseded_not_selected(published):
    reads = [json.loads(line) for line in (published / "fusion_reads.jsonl").read_text().splitlines()]
    by_id = {r["read_id"]: r for r in reads}
    supported = {rid for rid, r in by_id.items() if r["status"] == "supported_path"}
    assert len(supported & BCR_ABL1_MOLECULES) >= 5
    # every supported BCR-ABL1 molecule selected a SPLICED chr22 arm, never
    # the single intron-spanning block the SA summary describes
    for rid in supported & BCR_ABL1_MOLECULES:
        for path in by_id[rid]["selected"]:
            chr22 = [a for a in path["arms"] if a["chrom"] == "chr22"]
            assert chr22 and len(chr22[0]["blocks"]) >= 2
            assert all(b["reference_end"] - b["reference_start"] < 2000 for b in chr22[0]["blocks"])
        assert any(e["event"] == "sa_geometry_superseded" for e in by_id[rid]["sa_supersessions"])


def test_chain_level_artifacts_still_report_the_over_split_molecules(published):
    # the chain layer keys on exact topology, so degraded molecules of one
    # fusion are separate chains there; the transcript layer is what joins
    # them. Both views are published; neither is silently rewritten.
    chains = _rows(published / "fusions.tsv")
    bcr_abl1 = [c for c in chains if "ENSG00000186716" in c["gene_path"] and "ENSG00000097007" in c["gene_path"]]
    assert len(bcr_abl1) >= 4
    assert sum(int(c["full_read_support"]) for c in bcr_abl1) >= 5
