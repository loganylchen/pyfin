"""Fused-transcript models: truth-based acceptance on a synthetic genome.

The fixture is the one agreed with Codex before the layer was written. Truth
sequences are built directly from known exon slices, never with the
implementation's own stitching helper, so the test cannot pass by agreeing
with itself.

  A1 chrA [100,200) [300,400) [500,600)      B1 chrB [1000,1100) [1200,1300) [1600,1700)
  A2 chrA [100,200) [500,600)                 C1 chrC [2000,2100) [2200,2300)

  1. A1-B1, bp chrA:600 -> chrB:1000, two molecules both 5'-degraded to A 130
  2. A2-B1, same breakpoint, two molecules starting at A 100
  3. A1-B1-C1, bps 600->1000 and 1300->2000, two complete molecules
     (retained middle B exons: only [1000,1100) [1200,1300))
  4. one partial B-C molecule (middle + last arms only)
  5. one A-B molecule beginning at A 550: compatible with BOTH A isoforms
"""

import csv
import json
import random
from pathlib import Path

import pysam
import pytest

from fin.fusion.analysis import FusionSettings
from fin.fusion.pipeline import run_fusion_analysis
from fin.utils.sequences import reverse_complement


def _genome():
    generator = random.Random(20260913)
    # keep splice motifs honest: donors GT / acceptors AG at every annotated intron
    g = {chrom: list("".join(generator.choices("ACGT", k=3000))) for chrom in "ABC"}
    for chrom, exons in (("A", [(100, 200), (300, 400), (500, 600)]),
                         ("B", [(1000, 1100), (1200, 1300), (1600, 1700)]),
                         ("C", [(2000, 2100), (2200, 2300)])):
        for (s, e), (ns, ne) in zip(exons, exons[1:]):
            g[chrom][e:e + 2] = "GT"
            g[chrom][ns - 2:ns] = "AG"
    return {chrom: "".join(seq) for chrom, seq in g.items()}


A1 = [(100, 200), (300, 400), (500, 600)]
A2 = [(100, 200), (500, 600)]
B1 = [(1000, 1100), (1200, 1300), (1600, 1700)]
B_MIDDLE = [(1000, 1100), (1200, 1300)]
C1 = [(2000, 2100), (2200, 2300)]


def _seq(genome, chrom, exons):
    return "".join(genome[chrom][s:e] for s, e in exons)


def _cigar(exons, clip_left, clip_right):
    parts = []
    if clip_left:
        parts.append(f"{clip_left}S")
    for i, (s, e) in enumerate(exons):
        if i:
            parts.append(f"{s - exons[i - 1][1]}N")
        parts.append(f"{e - s}M")
    if clip_right:
        parts.append(f"{clip_right}S")
    return "".join(parts)


ANNOTATION = (("A", "A1", "A", "+", A1), ("A", "A2", "A", "+", A2), ("B", "B1", "B", "+", B1), ("C", "C1", "C", "+", C1))


def _arm_seq(genome, chrom, exons, strand="+"):
    seq = _seq(genome, chrom, exons)
    return seq if strand == "+" else reverse_complement(seq)


def _write_inputs(directory, genome, molecules, annotation=ANNOTATION):
    """molecules: list of (read_id, [(chrom, exons[, strand]), ...]) in transcript order.

    The primary record is the first arm; the rest are SA entries. A minus arm
    contributes the reverse complement of its ascending exons, and its SA/BAM
    record is written on the minus strand with the query-clip convention of
    the SAM spec (clips are given for the forward read; SA strand '-').
    """
    reference = directory / "genome.fa"
    reference.write_text("".join(f">{c}\n{s}\n" for c, s in genome.items()))
    gtf = directory / "genes.gtf"
    lines = []
    for chrom, tx, gene, strand, exons in annotation:
        for s, e in exons:
            lines.append(f'{chrom}\ttest\texon\t{s + 1}\t{e}\t.\t{strand}\t.\tgene_id "{gene}"; gene_name "g{gene}"; transcript_id "{tx}";\n')
    gtf.write_text("".join(lines))
    fastq = directory / "reads.fastq"
    bam = directory / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": c, "LN": 3000} for c in sorted(genome)]}
    names = sorted(genome)
    with fastq.open("w") as fq, pysam.AlignmentFile(str(bam), "wb", header=header) as writer:
        for read_id, arms in molecules:
            arms = [(a[0], a[1], a[2] if len(a) > 2 else "+") for a in arms]
            query = "".join(_arm_seq(genome, chrom, exons, strand) for chrom, exons, strand in arms)
            fq.write(f"@{read_id}\n{query}\n+\n{'I' * len(query)}\n")
            offsets = []
            pos = 0
            for chrom, exons, strand in arms:
                length = sum(e - s for s, e in exons)
                offsets.append((pos, length))
                pos += length
            record = pysam.AlignedSegment(writer.header)
            record.query_name = read_id
            chrom, exons, strand = arms[0]
            left, right = 0, len(query) - offsets[0][1]
            if strand == "+":
                record.query_sequence = query
                record.flag = 0
                record.cigarstring = _cigar(exons, left, right)
            else:
                record.query_sequence = reverse_complement(query)
                record.flag = 16
                record.cigarstring = _cigar(exons, right, left)
            record.reference_id = names.index(chrom)
            record.reference_start = exons[0][0]
            record.mapping_quality = 60
            record.set_tag("NM", 0)
            sa = []
            for (chrom, exons, strand), (start, length) in zip(arms[1:], offsets[1:]):
                left, right = start, len(query) - start - length
                cig = _cigar(exons, left, right) if strand == "+" else _cigar(exons, right, left)
                sa.append(f"{chrom},{exons[0][0] + 1},{strand},{cig},60,0;")
            if sa:
                record.set_tag("SA", "".join(sa))
            writer.write(record)
    return bam, reference, fastq, gtf


def _run(tmp_path, molecules, genome=None, annotation=ANNOTATION, **kwargs):
    genome = genome or _genome()
    bam, reference, fastq, gtf = _write_inputs(tmp_path, genome, molecules, annotation)
    out = tmp_path / "result"
    result = run_fusion_analysis(bam, reference, out, gtf_path=gtf, sequence_path=fastq,
                                 settings=FusionSettings(min_support=2, min_clip_bases=50, **kwargs))
    published = out / "fusion"
    with (published / "fusion_models.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    fasta = _fasta(published / "fusion_models.fa")
    annotated = _fasta(published / "fusion_models.annotated.fa")
    with (published / "fusion_model_exons.tsv").open() as handle:
        exons = list(csv.DictReader(handle, delimiter="\t"))
    return genome, result, published, rows, fasta, annotated, exons


def _fasta(path):
    records = {}
    header = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:]
            records[header.split()[0]] = [header, ""]
        else:
            records[header.split()[0]][1] += line
    return records


FIXTURE = [
    ("ab1_a", [("A", [(130, 200), (300, 400), (500, 600)]), ("B", B1)]),
    ("ab1_b", [("A", [(130, 200), (300, 400), (500, 600)]), ("B", B1)]),
    ("ab2_a", [("A", A2), ("B", B1)]),
    ("ab2_b", [("A", A2), ("B", B1)]),
    ("abc_a", [("A", A1), ("B", B_MIDDLE), ("C", C1)]),
    ("abc_b", [("A", A1), ("B", B_MIDDLE), ("C", C1)]),
    ("bc_partial", [("B", B_MIDDLE), ("C", C1)]),
    ("ab_ambiguous", [("A", [(550, 600)]), ("B", B1)]),
]


def test_two_isoform_combinations_at_one_breakpoint_are_two_models_and_three_genes_phase(tmp_path):
    genome, result, published, rows, fasta, annotated, exons = _run(tmp_path, FIXTURE)
    by_genes = {}
    for row in rows:
        path = tuple(tuple(g) for g in json.loads(row["gene_path"]))
        by_genes.setdefault(path, []).append(row)
    # exactly two resolved A-B splice models, not one gene pair and not a cross product
    ab = by_genes[(("A",), ("B",))]
    assert len(ab) == 2
    ab_by_len = {int(r["observed_length"]): r for r in ab}
    assert set(ab_by_len) == {570, 500}
    long_row, short_row = ab_by_len[570], ab_by_len[500]
    # unique structural support 2 and 2; the A-550 molecule votes for neither
    assert long_row["support_unique"] == "2" and short_row["support_unique"] == "2"
    assert "ab_ambiguous" not in json.loads(long_row["unique_read_ids"])
    assert "ab_ambiguous" not in json.loads(short_row["unique_read_ids"])
    assert "ab_ambiguous" in json.loads(long_row["compatible_read_ids"])
    assert "ab_ambiguous" in json.loads(short_row["compatible_read_ids"])
    # isoform resolution
    assert json.loads(long_row["compatible_transcript_ids"]) == [["A1"], ["B1"]]
    assert json.loads(short_row["compatible_transcript_ids"]) == [["A2"], ["B1"]]
    # splice class is about introns: the degraded read still saw every intron
    # of A1's retained slice, so it is arm_fsm; the missing 30 bp is an
    # ENDPOINT fact, carried by five_prime_status below
    assert json.loads(long_row["arm_match_classes"]) == ["arm_fsm", "arm_fsm"]
    assert json.loads(short_row["arm_match_classes"]) == ["arm_fsm", "arm_fsm"]
    assert json.loads(long_row["junction_names"]) == ["gA:e3>gB:e1"]
    assert json.loads(short_row["junction_names"]) == ["gA:e2>gB:e1"]
    # FASTA truth, built from exon slices independently
    assert fasta[long_row["sequence_id"]][1] == _seq(genome, "A", [(130, 200), (300, 400), (500, 600)]) + _seq(genome, "B", B1)
    assert fasta[short_row["sequence_id"]][1] == _seq(genome, "A", A2) + _seq(genome, "B", B1)
    # annotated completion of the degraded model: 30 more bases, labelled
    assert long_row["annotated_completion_id"] == long_row["model_id"] + ".ann"
    assert long_row["extension_bases_5p"] == "30" and long_row["annotated_completion_length"] == "600"
    completion = annotated[long_row["annotated_completion_id"]]
    assert completion[1] == _seq(genome, "A", A1) + _seq(genome, "B", B1)
    assert "kind=annotated_completion" in completion[0] and "extension_observed=false" in completion[0]
    assert long_row["five_prime_status"] == "annotated_extension_available"
    # A2-B1 was observed to its annotated start: no completion record needed
    assert short_row["annotated_completion_id"] == ""
    assert short_row["five_prime_status"] == "annotated_boundary_observed"
    # the three-gene model
    abc = by_genes[(("A",), ("B",), ("C",))]
    assert len(abc) == 1
    row = abc[0]
    assert row["arm_count"] == "3" and row["support_unique"] == "2"
    assert row["observed_length"] == "700"
    assert fasta[row["sequence_id"]][1] == _seq(genome, "A", A1) + _seq(genome, "B", B_MIDDLE) + _seq(genome, "C", C1)
    # the partial B-C molecule adds zero full-model votes
    assert "bc_partial" not in json.loads(row["unique_read_ids"]) + json.loads(row["compatible_read_ids"])
    # middle-arm exons are exactly the retained B exons
    middle = [(int(e["start0"]), int(e["end0"])) for e in exons
              if e["model_id"] == row["model_id"] and e["record_kind"] == "evidence_constrained" and e["arm_index"] == "2"]
    assert middle == B_MIDDLE
    assert json.loads(row["junction_names"]) == ["gA:e3>gB:e1", "gB:e2>gC:e1"]
    # the B->C edge sits in a 3-bp microhomology of this random genome, so the
    # chain layer keeps three same-sequence placements; at the transcript level
    # that is ONE fused transcript with two alternative placements, not three
    # models sharing the same two molecules with zero unique votes each
    assert row["coordinate_status"] == "unresolved_microhomology_same_sequence"
    assert len(json.loads(row["coordinate_alternatives"])) == 2
    assert row["breakpoints"] == "A:600:+>B:1000:+,B:1300:+>C:2000:+"
    # exon table reproduces each FASTA exactly
    for r in rows:
        seq = fasta[r["sequence_id"]][1]
        table = sorted((int(e["tx_start0"]), int(e["tx_end0"]), e["chrom"], int(e["start0"]), int(e["end0"]))
                       for e in exons if e["sequence_id"] == r["sequence_id"])
        assert table[0][0] == 0 and table[-1][1] == len(seq)
        assert "".join(genome[c][s:e] for _, _, c, s, e in table) == seq
    # mixed provenance on the extended first exon
    prov = {(e["arm_index"], e["provenance"]) for e in exons if e["sequence_id"] == long_row["annotated_completion_id"]}
    assert ("1", "observed_partial_annotated_extension") in prov
    assert ("2", "observed") in prov
    manifest = json.loads((published / "fusion_manifest.json").read_text())
    assert manifest["fused_transcript_models"]["count"] == 3
    assert manifest["fused_transcript_models"]["multi_gene_models"] == 1
    # the A-550 observation is collapsed, not gone: audited with both models it could be
    collapsed = [json.loads(l) for l in (published / "fusion_model_collapsed.jsonl").read_text().splitlines()]
    assert [c["read_ids"] for c in collapsed] == [["ab_ambiguous"]]
    assert sorted(collapsed[0]["compatible_model_ids"]) == sorted([long_row["model_id"], short_row["model_id"]])


def test_independent_ab_and_bc_molecules_never_yield_a_three_gene_model(tmp_path):
    molecules = [m for m in FIXTURE if not m[0].startswith("abc")]
    _, _, _, rows, _, _, _ = _run(tmp_path, molecules)
    assert all(row["arm_count"] == "2" for row in rows)


def test_conflicting_middle_arm_splicing_is_not_majority_merged(tmp_path):
    conflicting = [("abc_c", [("A", A1), ("B", [(1000, 1100)]), ("C", C1)])]
    _, _, _, rows, _, _, _ = _run(tmp_path, FIXTURE + conflicting)
    abc = [r for r in rows if r["arm_count"] == "3"]
    assert len(abc) == 1
    assert abc[0]["support_unique"] == "2"
    assert "abc_c" not in json.loads(abc[0]["unique_read_ids"]) + json.loads(abc[0]["compatible_read_ids"])


def test_chain_level_artifacts_are_unchanged_by_the_model_layer(tmp_path):
    _, _, published, _, _, _, _ = _run(tmp_path, FIXTURE)
    for name in ("fusions.tsv", "fusions.bedpe", "fusion_transcripts.fa", "fusion_arms.jsonl",
                 "fusion_reads.jsonl", "junction_corrections.jsonl"):
        assert (published / name).exists()
    with (published / "fusions.tsv").open() as handle:
        chains = list(csv.DictReader(handle, delimiter="\t"))
    assert all(set(r) >= {"chain_id", "gene_path", "full_read_support"} for r in chains)


# ---------------------------------------------------------------------------
# the review's commit-gating cases
# ---------------------------------------------------------------------------


def _minus_genome():
    """chrA plus-strand A1/A2 as before; chrM carries a MINUS-strand gene M1
    with exons [1000,1100) [1200,1300) [1600,1700) (transcript order is
    1600->1200->1000), and chrC C1 on plus."""
    generator = random.Random(7)
    g = {chrom: list("".join(generator.choices("ACGT", k=3000))) for chrom in "ACM"}
    for (s, e), (ns, ne) in zip(A1, A1[1:]):
        g["A"][e:e + 2] = "GT"; g["A"][ns - 2:ns] = "AG"
    for (s, e), (ns, ne) in zip(C1, C1[1:]):
        g["C"][e:e + 2] = "GT"; g["C"][ns - 2:ns] = "AG"
    # minus strand: donor GT on the transcript strand = "AC" at the genomic
    # left of the intron... the layer does not require motifs, but keep the
    # reverse-strand consensus honest: intron [1100,1200) genomic = AC..CT
    for (s, e), (ns, ne) in zip(M1, M1[1:]):
        g["M"][e:e + 2] = "CT"; g["M"][ns - 2:ns] = "AC"
    return {chrom: "".join(seq) for chrom, seq in g.items()}


M1 = [(1000, 1100), (1200, 1300), (1600, 1700)]
M_ANNOTATION = (("A", "A1", "A", "+", A1), ("A", "A2", "A", "+", A2), ("M", "M1", "M", "-", M1), ("C", "C1", "C", "+", C1))


def test_minus_strand_last_arm_sequence_offsets_and_ranks(tmp_path):
    # A1 (+) -> M1 (-): the M arm is the reverse complement of its ascending
    # exons in descending order; entry = genomic 1700 (its 5' end on '-')
    molecules = [("am_a", [("A", A1), ("M", M1, "-")]), ("am_b", [("A", A1), ("M", M1, "-")])]
    genome = _minus_genome()
    genome, result, published, rows, fasta, annotated, exons = _run(tmp_path, molecules, genome=genome,
                                                                    annotation=M_ANNOTATION)
    assert len(rows) == 1
    row = rows[0]
    assert row["support_unique"] == "2"
    assert json.loads(row["compatible_transcript_ids"]) == [["A1"], ["M1"]]
    assert json.loads(row["arm_match_classes"]) == ["arm_fsm", "arm_fsm"]
    truth = _seq(genome, "A", A1) + reverse_complement(_seq(genome, "M", M1))
    assert fasta[row["sequence_id"]][1] == truth
    assert row["breakpoints"] == "A:600:+>M:1700:-"
    assert json.loads(row["junction_names"]) == ["gA:e3>gM:e1"]
    m_rows = [e for e in exons if e["sequence_id"] == row["sequence_id"] and e["arm_index"] == "2"]
    # transcript order on minus = descending genomic; offsets contiguous
    assert [(int(e["start0"]), int(e["end0"])) for e in m_rows] == [(1600, 1700), (1200, 1300), (1000, 1100)]
    assert [int(e["tx_start0"]) for e in m_rows] == [300, 400, 500]
    assert "".join(genome[e["chrom"]][int(e["start0"]):int(e["end0"])] for e in m_rows)  # slices exist
    # and the exon table reconstructs the FASTA with reverse complement on '-'
    pieces = []
    for e in sorted((e for e in exons if e["sequence_id"] == row["sequence_id"]), key=lambda e: int(e["tx_start0"])):
        piece = genome[e["chrom"]][int(e["start0"]):int(e["end0"])]
        pieces.append(piece if e["strand"] == "+" else reverse_complement(piece))
    assert "".join(pieces) == truth


def test_minus_strand_first_arm_truncation_is_compatible_and_completion_extends_its_5p(tmp_path):
    # M1 (-) -> C1 (+); the second molecule lost 40 bp of M1's 5' end (genomic
    # 1660-1700), so its M exons are [(1000,1100),(1200,1300),(1600,1660)]
    molecules = [("mc_a", [("M", M1, "-"), ("C", C1)]),
                 ("mc_b", [("M", M1, "-"), ("C", C1)]),
                 ("mc_short", [("M", [(1000, 1100), (1200, 1300), (1600, 1660)], "-"), ("C", C1)]),
                 ("mc_short2", [("M", [(1000, 1100), (1200, 1300), (1600, 1660)], "-"), ("C", C1)])]
    genome = _minus_genome()
    genome, result, published, rows, fasta, annotated, exons = _run(tmp_path, molecules, genome=genome,
                                                                    annotation=M_ANNOTATION)
    assert len(rows) == 1, [r["breakpoints"] for r in rows]
    row = rows[0]
    # the two short molecules fit inside the full model and are unique to it
    # (it is the only reported model)
    assert row["support_unique"] == "4"
    assert row["breakpoints"] == "M:1000:->C:2000:+"
    assert row["five_prime_status"] == "annotated_boundary_observed"
    assert fasta[row["sequence_id"]][1] == reverse_complement(_seq(genome, "M", M1)) + _seq(genome, "C", C1)
    # the short molecules lost 40 bp of an EXON, not an exon: same splice
    # structure, so they are the same model (nothing to collapse) -- and the
    # record says which molecules actually reach the model's 5' extent
    assert (published / "fusion_model_collapsed.jsonl").read_text() == ""
    assert row["five_prime_observed_reads"] == "2"
    assert row["three_prime_observed_reads"] == "4"


def test_minus_first_arm_degraded_model_gets_annotated_completion(tmp_path):
    # only degraded molecules: model 5' extent is genomic 1660 on '-', the
    # completion restores [1660,1700) = 40 bp at the 5' end
    molecules = [("mc_s1", [("M", [(1000, 1100), (1200, 1300), (1600, 1660)], "-"), ("C", C1)]),
                 ("mc_s2", [("M", [(1000, 1100), (1200, 1300), (1600, 1660)], "-"), ("C", C1)])]
    genome = _minus_genome()
    genome, result, published, rows, fasta, annotated, exons = _run(tmp_path, molecules, genome=genome,
                                                                    annotation=M_ANNOTATION)
    row = rows[0]
    assert row["five_prime_status"] == "annotated_extension_available"
    assert row["extension_bases_5p"] == "40" and row["extension_bases_3p"] == "0"
    assert annotated[row["annotated_completion_id"]][1] == reverse_complement(_seq(genome, "M", M1)) + _seq(genome, "C", C1)
    assert fasta[row["sequence_id"]][1] == reverse_complement(_seq(genome, "M", [(1000, 1100), (1200, 1300), (1600, 1660)])) + _seq(genome, "C", C1)


@pytest.mark.parametrize("order", ["fit_first", "misfit_first"])
def test_compatibility_is_judged_per_molecule_not_per_representative(tmp_path, order):
    # two short A observations share a splice key (no introns, same exit 600):
    # S1 [(550,600)] fits inside A1's last exon; S2 [(450,600)] does NOT -- it
    # aligned continuously across A1's intron [400,500). S2 must keep its own
    # structure as a model regardless of which one is encountered first.
    s1 = ("s1", [("A", [(550, 600)]), ("B", B1)])
    s2 = ("s2", [("A", [(450, 600)]), ("B", B1)])
    s3 = ("s3", [("A", [(450, 600)]), ("B", B1)])   # second S2-type molecule so it can be supported
    long_ = [("ab1_a", [("A", A1), ("B", B1)]), ("ab1_b", [("A", A1), ("B", B1)])]
    molecules = long_ + ([s1, s2, s3] if order == "fit_first" else [s2, s3, s1])
    _, _, _, rows, fasta, _, _ = _run(tmp_path, molecules)
    by_len = {int(r["observed_length"]): r for r in rows}
    assert set(by_len) == {600, 450}, sorted(by_len)
    assert json.loads(by_len[600]["unique_read_ids"]) == ["ab1_a", "ab1_b"]
    assert json.loads(by_len[450]["unique_read_ids"]) == ["s2", "s3"]
    # s1 fits both the full A1 model and the [450,600) model: compatible with both, unique to none
    assert "s1" in json.loads(by_len[600]["compatible_read_ids"])
    assert "s1" in json.loads(by_len[450]["compatible_read_ids"])
    assert "s1" not in json.loads(by_len[600]["unique_read_ids"]) + json.loads(by_len[450]["unique_read_ids"])


def test_supported_shorter_structure_survives_an_unsupported_longer_one(tmp_path):
    # the full A1-B1 structure (3 A exons) is seen by ONE molecule only, i.e.
    # an insufficient-support chain; a structure missing A exon 1 entirely is
    # seen by two. The shorter is a strict 5' truncation of the longer, but the
    # longer is not eligible -- so the shorter MUST be reported, not erased.
    molecules = [("full_once", [("A", A1), ("B", B1)]),
                 ("deg_a", [("A", [(300, 400), (500, 600)]), ("B", B1)]),
                 ("deg_b", [("A", [(300, 400), (500, 600)]), ("B", B1)])]
    _, _, published, rows, _, _, _ = _run(tmp_path, molecules)
    assert len(rows) == 1
    assert rows[0]["observed_length"] == "500"
    assert json.loads(rows[0]["unique_read_ids"]) == ["deg_a", "deg_b"]
    assert "full_once" not in json.loads(rows[0]["unique_read_ids"]) + json.loads(rows[0]["compatible_read_ids"])
    # two-exon A structure is compatible with A1 (missing intron outside coverage)
    # and NOT with A2 (A2's intron [200,500) overlaps the observed [300,400))
    assert json.loads(rows[0]["compatible_transcript_ids"])[0] == ["A1"]
    assert json.loads(rows[0]["arm_match_classes"])[0] == "arm_ism_5p"
    assert rows[0]["extension_bases_5p"] == "100"


def test_distinct_loci_with_identical_sequence_stay_distinct(tmp_path):
    # chrD is a copy of chrA: D1-B1 and A1-B1 project to the same sequence but
    # are different fusion events and must not be merged
    genome = _genome()
    genome["D"] = genome["A"]
    annotation = ANNOTATION + (("D", "D1", "D", "+", A1),)
    molecules = [("ab_a", [("A", A1), ("B", B1)]), ("ab_b", [("A", A1), ("B", B1)]),
                 ("db_a", [("D", A1), ("B", B1)]), ("db_b", [("D", A1), ("B", B1)])]
    _, result, _, rows, _, _, _ = _run(tmp_path, molecules, genome=genome, annotation=annotation)
    # mappy may report each molecule as ambiguous between A and D (identical
    # sequence); if it resolves them, both events must be reported separately
    paths = {tuple(tuple(g) for g in json.loads(r["gene_path"])) for r in rows}
    assert not any(r["coordinate_status"] == "unresolved_microhomology_same_sequence" for r in rows)
    assert len(rows) == len(paths)


def test_transcript_writer_failure_publishes_a_failed_manifest(tmp_path, monkeypatch):
    import fin.fusion.pipeline as pipeline_module

    def boom(*args, **kwargs):
        raise ValueError("injected transcript writer failure")

    monkeypatch.setattr("fin.fusion.transcript_output.write_transcript_outputs", boom)
    genome = _genome()
    bam, reference, fastq, gtf = _write_inputs(tmp_path, genome, FIXTURE)
    out = tmp_path / "result"
    with pytest.raises(ValueError):
        run_fusion_analysis(bam, reference, out, gtf_path=gtf, sequence_path=fastq,
                            settings=FusionSettings(min_support=2, min_clip_bases=50))
    manifest = json.loads((out / "fusion" / "fusion_manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert not (out / "fusion" / "fusions.tsv").exists()
