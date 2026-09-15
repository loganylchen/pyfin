"""Independent file-to-output fusion analysis using real BAM and minimap2."""

import csv
import hashlib
import json
import random
import subprocess
import sys
from pathlib import Path

import pysam
import pytest

from fin.fusion.analysis import FusionSettings
from fin.fusion.output import FusionPublication
from fin.fusion.pipeline import run_fusion_analysis
from fin.utils.sequences import reverse_complement


def _inputs(directory, *, count=3, reverse_records=False, unique_edges=False):
    generator = random.Random(591)
    chromosomes = "ABCD"[:count]
    genome = {chrom: "".join(generator.choices("ACGT", k=2000)) for chrom in chromosomes}
    if unique_edges:
        for left, right in zip(chromosomes, chromosomes[1:]):
            left_sequence, right_sequence = list(genome[left]), list(genome[right])
            right_sequence[499] = next(base for base in "ACGT" if base != left_sequence[799])
            left_sequence[800] = next(base for base in "ACGT" if base != right_sequence[500])
            genome[left], genome[right] = "".join(left_sequence), "".join(right_sequence)
    reference = directory / "genome.fa"
    reference.write_text("".join(f">{chrom}\n{sequence}\n" for chrom, sequence in genome.items()))
    query = "".join(sequence[500:800] for sequence in genome.values())
    fastq = directory / "reads.fastq"
    fastq.write_text("".join(f"@read{index}\n{query}\n+\n{'I' * len(query)}\n" for index in range(2)))
    bam = directory / "reads.bam"
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": chrom, "LN": 2000} for chrom in chromosomes]}
    with pysam.AlignmentFile(str(bam), "wb", header=header) as writer:
        for index in ([1, 0] if reverse_records else [0, 1]):
            record = pysam.AlignedSegment(writer.header)
            record.query_name = f"read{index}"
            record.query_sequence = query
            record.flag = 0
            record.reference_id = 0
            record.reference_start = 500
            record.mapping_quality = 60
            record.cigarstring = f"300M{len(query) - 300}S"
            record.set_tag("NM", 0)
            writer.write(record)
    gtf = directory / "genes.gtf"
    gtf.write_text("".join(f'{chrom}\ttest\texon\t401\t900\t.\t+\t.\tgene_id "{chrom}"; transcript_id "{chrom}.tx";\n'
                           for chrom in chromosomes))
    return bam, reference, fastq, gtf


def test_four_gene_file_pipeline_emits_all_edges_and_provenance(tmp_path):
    bam, reference, fastq, gtf = _inputs(tmp_path, count=4)
    output = tmp_path / "result"
    output.mkdir()
    sentinels = {name: f"unchanged {name}\n" for name in ("assembly.gtf", "scores.tsv", "abundance_refit.json")}
    for name, content in sentinels.items():
        (output / name).write_text(content)
    result = run_fusion_analysis(bam, reference, output, gtf_path=gtf, sequence_path=fastq)
    assert len(result.chains) == 1
    assert len(result.chains[0].representative.arms) == 4
    published = output / "fusion"
    manifest = json.loads((published / "fusion_manifest.json").read_text())
    assert manifest["status"] == "complete"
    assert manifest["signal_status"] == "not_scored"
    for name, digest in manifest["files"].items():
        assert hashlib.sha256((published / name).read_bytes()).hexdigest() == digest
    with (published / "fusions.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert json.loads(rows[0]["gene_path"]) == [["A"], ["B"], ["C"], ["D"]]
    assert rows[0]["full_read_support"] == "2"
    with (published / "fusions.bedpe").open() as handle:
        edges = list(csv.DictReader(handle, delimiter="\t"))
    assert [(edge["#chrom1"], edge["chrom2"]) for edge in edges] == [("A", "B"), ("B", "C"), ("C", "D")]
    assert [edge["edge_index"] for edge in edges] == ["1", "2", "3"]
    assert (published / "fusion_transcripts.fa").read_text().startswith(">fusion_chain_")
    assert (published / "junction_corrections.jsonl").stat().st_size > 0
    for name, content in sentinels.items():
        assert (output / name).read_text() == content


def test_pipeline_reads_final_assembly_evidence_and_audits_signal_abstention(tmp_path, monkeypatch):
    from fin.fusion.signal import KrillFusionBackend

    bam, reference, fastq, gtf = _inputs(tmp_path)
    assembly = tmp_path / "assembly.gtf"
    assembly.write_text('A\ttest\texon\t501\t800\t.\t+\t.\tgene_id "A"; transcript_id "A.final";\n'
                        'A\ttest\texon\t1001\t1100\t.\t+\t.\tgene_id "A"; transcript_id "A.final";\n')
    support = tmp_path / "junction_support.tsv"
    support.write_text("chrom\tstrand\tdonor\tacceptor\traw_reads\tindel_clean_reads\nA\t+\t800\t1000\t6\t4\n")
    signal = tmp_path / "unavailable.slow5"
    signal.write_text("test backend explicitly returns missing read records\n")
    monkeypatch.setattr(KrillFusionBackend, "__call__", lambda *args: ())
    result = run_fusion_analysis(bam, reference, tmp_path / "result", gtf_path=gtf, sequence_path=fastq,
                                 signal_path=signal, assembly_gtf_path=assembly, assembly_support_path=support)
    assert result.signal_status == "abstained"
    assert len(result.chains) == 1
    audits = result.reads[0].assembly_junctions
    assert any(audit["junctions"][0]["donor_clean_reads"] == 4 for audit in audits)
    published = tmp_path / "result" / "fusion"
    manifest = json.loads((published / "fusion_manifest.json").read_text())
    assert manifest["assembly_junction_evidence"]["junction_count"] == 1
    assert {"assembly_gtf", "assembly_junction_support", "signal"} <= manifest["inputs"].keys()
    assert (published / "junction_signal.jsonl").stat().st_size > 0
    with (published / "fusions.tsv").open() as handle:
        assert next(csv.DictReader(handle, delimiter="\t"))["signal_status"] == "sequence_only_signal_abstained"


def test_failed_rerun_does_not_expose_previous_fusion_calls(tmp_path):
    output = tmp_path / "output"
    with FusionPublication(output) as publication:
        (publication.work_dir / "fusions.tsv").write_text("old calls\n")
        publication.complete({})
    assert (output / "fusion" / "fusions.tsv").exists()
    with pytest.raises(RuntimeError, match="synthetic failure"):
        with FusionPublication(output) as publication:
            assert not (output / "fusion" / "fusions.tsv").exists()
            (publication.work_dir / "fusions.tsv").write_text("partial calls\n")
            raise RuntimeError("synthetic failure")
    manifest = json.loads((output / "fusion" / "fusion_manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert not (output / "fusion" / "fusions.tsv").exists()


def test_concurrent_publication_cannot_replace_an_active_run(tmp_path):
    with FusionPublication(tmp_path):
        with pytest.raises(BlockingIOError):
            with FusionPublication(tmp_path):
                pytest.fail("second writer acquired output lock")


def test_unmanaged_output_directory_is_not_deleted(tmp_path):
    existing = tmp_path / "fusion"
    existing.mkdir()
    (existing / "important.txt").write_text("keep")
    with pytest.raises(ValueError, match="managed"):
        with FusionPublication(tmp_path):
            pytest.fail("replaced unmanaged directory")
    assert (existing / "important.txt").read_text() == "keep"


def test_invalid_settings_publish_failure_not_old_success(tmp_path):
    bam, reference, fastq, gtf = _inputs(tmp_path)
    with pytest.raises(ValueError, match="settings"):
        run_fusion_analysis(bam, reference, tmp_path / "result", gtf_path=gtf, sequence_path=fastq,
                            settings=FusionSettings(min_support=0))
    manifest = json.loads((tmp_path / "result" / "fusion" / "fusion_manifest.json").read_text())
    assert manifest["status"] == "failed"


def _add_ordinary_assembly_inputs(directory):
    with pysam.FastaFile(str(directory / "genome.fa")) as reference:
        query = reference.fetch("A", 900, 1050) + reference.fetch("A", 1350, 1500)
    bam = directory / "reads.bam"
    unsorted = directory / "unsorted.bam"
    with pysam.AlignmentFile(str(bam), "rb") as reader:
        with pysam.AlignmentFile(str(unsorted), "wb", header=reader.header) as writer:
            for record in reader:
                writer.write(record)
            for index in range(4):
                record = pysam.AlignedSegment(writer.header)
                record.query_name = f"normal{index}"
                record.query_sequence = query
                record.flag = 0
                record.reference_id = 0
                record.reference_start = 900
                record.mapping_quality = 60
                record.cigarstring = "150M300N150M"
                record.set_tag("NM", 0)
                writer.write(record)
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    with (directory / "reads.fastq").open("a") as handle:
        for index in range(4):
            handle.write(f"@normal{index}\n{query}\n+\n{'I' * 300}\n")
    with (directory / "genes.gtf").open("a") as handle:
        for start, end in ((901, 1050), (1351, 1500)):
            handle.write(f'A\ttest\texon\t{start}\t{end}\t.\t+\t.\tgene_id "normal"; transcript_id "normal.tx";\n')
    signal = directory / "signal.slow5"
    signal.write_text(
        "#slow5_version\t1.0.0\n#num_read_groups\t1\n"
        "#char*\tuint32_t\tdouble\tdouble\tdouble\tdouble\tuint64_t\tint16_t*\n"
        "#read_id\tread_group\tdigitisation\toffset\trange\tsampling_rate\tlen_raw_signal\traw_signal\n"
        "read0\t0\t8192\t0\t1467\t4000\t8\t500,501,502,503,504,505,506,507\n"
    )
    return signal


def test_real_cli_assembly_bytes_match_prechange_baseline_with_fusion_off_and_on(tmp_path):
    """Golden hashes measured from HEAD 945fdae with M2 and endpoint channels off."""
    pytest.importorskip("pyslow5", reason="real CLI fixture requires the SLOW5 reader")
    bam, reference, fastq, gtf = _inputs(tmp_path)
    signal = _add_ordinary_assembly_inputs(tmp_path)
    source_root = Path(__file__).resolve().parents[2]
    command = [sys.executable, "-m", "fin.cli", "--bam", str(bam), "--genome", str(reference),
               "--gtf", str(gtf), "--fastq", str(fastq), "--signal", str(signal),
               "--signal-format", "slow5", "--no-gpu", "--m2-metric", "off",
               "--assignment-channels", "none", "--folded-tss-recovery-mode", "off"]
    for mode in ("off", "on"):
        output = tmp_path / mode
        execution = subprocess.run(command + ["--output-dir", str(output)] + (["--fusion"] if mode == "on" else []),
                                   cwd=source_root, capture_output=True, text=True, timeout=90)
        assert execution.returncode == 0, execution.stdout + execution.stderr
    baseline = {
        "assembly.gtf": "b3731bf89bd91d817b9813461cd3a544c2378383c5634cde2c295c843222b3f3",
        "scores.tsv": "bcffdbd5c5399898fe73333cc2993acfa908224b24370506d6084ceec980f872",
        "abundance_refit.json": "1cfc99f3b6538d5acd44a5de55bd087fef29191f446fe1b5c69b753d4f4b6a0e",
    }
    for name, expected in baseline.items():
        off = (tmp_path / "off" / name).read_bytes()
        on = (tmp_path / "on" / name).read_bytes()
        assert off == on
        assert hashlib.sha256(off).hexdigest() == expected
    with (tmp_path / "on" / "fusion" / "fusions.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["arm_count"] == "3"
    assert rows[0]["full_read_support"] == "2"
    assert not (tmp_path / "off" / "fusion").exists()


def test_shifted_fusion_junction_is_corrected_in_file_to_output_path(tmp_path):
    bam, reference, fastq, gtf = _inputs(tmp_path, unique_edges=True)
    shifted = tmp_path / "shifted.bam"
    with pysam.AlignmentFile(str(bam), "rb") as reader:
        with pysam.AlignmentFile(str(shifted), "wb", header=reader.header) as writer:
            for record in reader:
                record.cigarstring = "302M598S"
                writer.write(record)
    result = run_fusion_analysis(shifted, reference, tmp_path / "result", gtf_path=gtf, sequence_path=fastq)
    assert len(result.chains) == 1
    assert result.chains[0].full_support_read_ids == ("read0", "read1")
    for read in result.reads:
        assert read.status == "supported_path"
        assert {path.arms[0].exit_position for path in read.selected} == {800}
        assert any(item.observed.arms[0].exit_position == 802
                   and item.proposed.arms[0].exit_position == 800 for item in read.evaluated)


def test_record_order_does_not_change_published_biological_artifacts(tmp_path):
    outputs = []
    for reverse_records in (False, True):
        directory = tmp_path / str(reverse_records)
        directory.mkdir()
        bam, reference, fastq, gtf = _inputs(directory, reverse_records=reverse_records, unique_edges=True)
        run_fusion_analysis(bam, reference, directory / "result", gtf_path=gtf, sequence_path=fastq)
        outputs.append(directory / "result" / "fusion")
    for name in ("fusions.tsv", "fusions.bedpe", "fusion_transcripts.fa", "fusion_arms.jsonl",
                 "fusion_reads.jsonl", "junction_corrections.jsonl"):
        assert (outputs[0] / name).read_bytes() == (outputs[1] / name).read_bytes()


def test_reverse_middle_primary_and_hard_clip_recovery_preserve_complete_chain(tmp_path):
    bam, reference, fastq, gtf = _inputs(tmp_path)
    with pysam.FastaFile(str(reference)) as genome:
        query = genome.fetch("A", 500, 800) + reverse_complement(genome.fetch("B", 500, 800)) + genome.fetch("C", 500, 800)
    fastq.write_text("".join(f"@read{index}\n{query}\n+\n{'I' * len(query)}\n" for index in range(2)))
    gtf.write_text(gtf.read_text().replace("B\ttest\texon\t401\t900\t.\t+", "B\ttest\texon\t401\t900\t.\t-"))
    with pysam.AlignmentFile(str(bam), "rb") as reader:
        header = reader.header
    tables = []
    for mode in ("first_primary", "middle_primary", "hard_only"):
        candidate_bam = tmp_path / f"{mode}.bam"
        with pysam.AlignmentFile(str(candidate_bam), "wb", header=header) as writer:
            for index in range(2):
                for arm_index, chrom in enumerate("ABC"):
                    primary = arm_index == (0 if mode == "first_primary" else 1)
                    hard = not primary or mode == "hard_only"
                    front = arm_index * 300
                    back = 900 - front - 300
                    if chrom == "B":
                        front, back = back, front
                    clip = "H" if hard else "S"
                    record = pysam.AlignedSegment(writer.header)
                    record.query_name = f"read{index}"
                    record.flag = (0 if primary else 2048) | (16 if chrom == "B" else 0)
                    record.reference_id = writer.header.get_tid(chrom)
                    record.reference_start = 500
                    record.mapping_quality = 60
                    record.cigarstring = (f"{front}{clip}" if front else "") + "300M" + (f"{back}{clip}" if back else "")
                    sequence = query[arm_index * 300:(arm_index + 1) * 300] if hard else query
                    record.query_sequence = reverse_complement(sequence) if chrom == "B" else sequence
                    record.set_tag("NM", 0)
                    writer.write(record)
        output = tmp_path / mode
        result = run_fusion_analysis(candidate_bam, reference, output, gtf_path=gtf, sequence_path=fastq)
        assert len(result.chains) == 1
        assert result.chains[0].representative.gene_path == (("A",), ("B",), ("C",))
        assert tuple(arm.strand for arm in result.chains[0].representative.arms) == ("+", "-", "+")
        assert result.chains[0].full_support_read_ids == ("read0", "read1")
        tables.append((output / "fusion" / "fusions.tsv").read_bytes())
    assert tables[0] == tables[1] == tables[2]


def test_internal_splice_and_fusion_edges_are_jointly_corrected_end_to_end(tmp_path):
    bam, reference, fastq, gtf = _inputs(tmp_path, unique_edges=True)
    with pysam.FastaFile(str(reference)) as reader:
        genome = {chrom: reader.fetch(chrom) for chrom in "ABC"}
    middle = list(genome["B"])
    middle[800:802] = "GT"
    middle[1098:1100] = "AG"
    middle[1400] = next(base for base in "ACGT" if base != genome["C"][500])
    genome["B"] = "".join(middle)
    reference.write_text("".join(f">{chrom}\n{sequence}\n" for chrom, sequence in genome.items()))
    query = genome["A"][500:800] + genome["B"][500:800] + genome["B"][1100:1400] + genome["C"][500:800]
    fastq.write_text("".join(f"@read{index}\n{query}\n+\n{'I' * len(query)}\n" for index in range(2)))
    gtf.write_text(
        'A\ttest\texon\t501\t800\t.\t+\t.\tgene_id "A"; transcript_id "A.tx";\n'
        'B\ttest\texon\t501\t800\t.\t+\t.\tgene_id "B"; transcript_id "B.tx";\n'
        'B\ttest\texon\t1101\t1400\t.\t+\t.\tgene_id "B"; transcript_id "B.tx";\n'
        'C\ttest\texon\t501\t800\t.\t+\t.\tgene_id "C"; transcript_id "C.tx";\n'
    )
    with pysam.AlignmentFile(str(bam), "rb") as reader:
        header = reader.header
    with pysam.AlignmentFile(str(bam), "wb", header=header) as writer:
        for index in range(2):
            record = pysam.AlignedSegment(writer.header)
            record.query_name = f"read{index}"
            record.query_sequence = query
            record.flag = 0
            record.reference_id = 0
            record.reference_start = 500
            record.mapping_quality = 60
            record.cigarstring = "300M900S"
            record.set_tag("SA", "B,501,+,300S302M298N298M300S,60,4;C,501,+,900S300M,60,0;")
            record.set_tag("NM", 0)
            writer.write(record)
            # The shifted intron must be an OBSERVED geometry for the corrector
            # to act on. An SA entry is only a locator (minimap2 writes SA
            # CIGARs as one M plus one I/D), so it is superseded once the
            # realignment describes the same placement in full; a
            # supplementary RECORD carrying the shifted CIGAR is detailed
            # evidence and stays observed.
            supplementary = pysam.AlignedSegment(writer.header)
            supplementary.query_name = f"read{index}"
            supplementary.query_sequence = query
            supplementary.flag = 2048
            supplementary.reference_id = 1
            supplementary.reference_start = 500
            supplementary.mapping_quality = 60
            supplementary.cigarstring = "300S302M298N298M300S"
            supplementary.set_tag("SA", "A,501,+,300M900S,60,0;C,501,+,900S300M,60,0;")
            supplementary.set_tag("NM", 4)
            writer.write(supplementary)
    result = run_fusion_analysis(bam, reference, tmp_path / "result", gtf_path=gtf, sequence_path=fastq)
    assert len(result.chains) == 1
    assert result.chains[0].representative.gene_path == (("A",), ("B",), ("C",))
    assert result.chains[0].full_support_read_ids == ("read0", "read1")
    for read in result.reads:
        assert read.status == "supported_path"
        assert any(path.arms[1].exons == ((500, 800), (1100, 1400)) for path in read.selected)
        assert any(item.correction is not None and any(
            junction.kind == "intron" and junction.decision.observed_positions == (802, 1100)
            and any(contrast.positions == (800, 1100) and "annotation" in contrast.origin
                    for contrast in junction.decision.contrasts)
            for junction in item.correction.junctions) for item in read.evaluated)


def test_annotation_free_file_pipeline_reports_loci_without_inventing_gene_names(tmp_path):
    bam, reference, fastq, _ = _inputs(tmp_path, unique_edges=True)
    result = run_fusion_analysis(bam, reference, tmp_path / "result", sequence_path=fastq)
    assert result.annotation_status == "unavailable_loci_only"
    assert len(result.chains) == 1
    assert result.chains[0].representative.gene_path == ((), (), ())
    assert tuple(arm.chrom for arm in result.chains[0].representative.arms) == ("A", "B", "C")
    manifest = json.loads((tmp_path / "result" / "fusion" / "fusion_manifest.json").read_text())
    assert manifest["annotation_status"] == "unavailable_loci_only"


def test_a_remount_is_not_a_changed_input_but_a_rewrite_is_named(tmp_path):
    """The publication check killed two whole-library runs after every artifact
    was computed. The first said only "fusion input changed before publication"
    and was undiagnosable; the second, once the message named the input and the
    differing fields, reported `st_dev 144->138` on an /autofs fastq whose
    inode, size, mtime and ctime were unchanged -- the automounter had expired
    the share and remounted it mid-run.

    So: a device number must not fail the run, a real rewrite still must, and
    the message must name the field that actually differs. That last clause is
    why `SNAPSHOT_FIELDS` is shared rather than spelled twice -- labelling the
    tuple positionally from a stale list would have reported an inode change as
    a device change, resurrecting the exact message this fix removes.
    """
    import os
    import pathlib

    from fin.fusion.input import SNAPSHOT_FIELDS

    bam, reference, fastq, gtf = _inputs(tmp_path)
    assert "st_dev" not in SNAPSHOT_FIELDS

    real = pathlib.Path.stat
    remounted_paths = {str(fastq)}

    def remounted(self, *args, **kwargs):
        original = real(self, *args, **kwargs)
        if str(self) not in remounted_paths:
            return original
        fields = list(original)
        fields[2] += 6  # st_dev
        return os.stat_result(fields, {"st_atime_ns": original.st_atime_ns,
                                       "st_mtime_ns": original.st_mtime_ns,
                                       "st_ctime_ns": original.st_ctime_ns})

    original_stat = pathlib.Path.stat
    pathlib.Path.stat = remounted
    try:
        result = run_fusion_analysis(bam, reference, tmp_path / "remounted",
                                     gtf_path=gtf, sequence_path=fastq)
    finally:
        pathlib.Path.stat = original_stat
    assert result.chains

    # A genuine rewrite of the same path is still rejected, and named. Rewrite
    # it from inside the run, after the inputs were snapshotted -- which is
    # exactly when the real automount drift was observed.
    import fin.fusion.pipeline as pipeline_module

    original_writer = pipeline_module.write_fusion_outputs

    def rewrite_then_write(analysis, genome, work_dir):
        fastq.write_text(fastq.read_text() + "@extra\nACGT\n+\nIIII\n")
        return original_writer(analysis, genome, work_dir)

    pipeline_module.write_fusion_outputs = rewrite_then_write
    try:
        with pytest.raises(Exception) as failure:
            run_fusion_analysis(bam, reference, tmp_path / "rewritten",
                                gtf_path=gtf, sequence_path=fastq)
    finally:
        pipeline_module.write_fusion_outputs = original_writer
    message = str(failure.value)
    assert "fusion input changed before publication" in message
    assert "sequence" in message and "st_size" in message
    assert "st_dev" not in message
