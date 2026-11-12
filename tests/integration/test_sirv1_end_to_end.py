"""The whole command, on real reads with their real signal.

``tests/data/sirv1_mini`` is one SIRV locus of an SG-NEx direct-RNA run: 91
reads, their raw signal, the contig and its eight annotated transcripts. Under
the shipped profile it has exact alignment-score ties that only the signal can
break, so these tests reach the code a sequence-only fixture cannot.

They need the krill signal backend. CI installs it and fails if it is missing;
elsewhere a missing backend skips them.
"""
from __future__ import annotations

import importlib.util
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
DATA = REPO / "tests" / "data" / "sirv1_mini"

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("krill") is None, reason="the krill signal backend is not installed"
)

_spec = importlib.util.spec_from_file_location(
    "v1_struct_compare", REPO / "experiments" / "prod_validation" / "v1_struct_compare.py"
)
struct_compare = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(struct_compare)


def _pyfin(out: Path, *extra, data: Path = DATA, gtf: bool = True):
    args = [
        sys.executable, "-m", "fin.cli",
        "--bam", str(data / "reads.bam"),
        "--genome", str(data / "reference.fa"),
        "--signal", str(data / "reads.blow5"),
        "--output-dir", str(out),
        "--no-gpu",
        *extra,
    ]
    if gtf:
        args += ["--gtf", str(data / "annotation.gtf")]
    return subprocess.run(args, capture_output=True, text=True, cwd=str(REPO))


def _transcripts(gtf: Path):
    rows = []
    for line in gtf.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) > 8 and fields[2] == "transcript":
            rows.append(dict(re.findall(r'(\S+) "([^"]*)"', fields[8])))
    return rows


@pytest.fixture(scope="module")
def default_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("default")
    result = _pyfin(out)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
    return out, result.stdout + result.stderr


def test_the_default_command_assembles_the_locus(default_run):
    out, _log = default_run
    models = _transcripts(out / "assembly.gtf")
    annotated = {m["transcript_id"].removesuffix(".ext") for m in models if m["transcript_source"] == "gtf"}
    catalogue = {m["transcript_id"] for m in _transcripts(DATA / "annotation.gtf")}

    assert (out / "scores.tsv").exists() and (out / "run_manifest.json").exists()
    assert annotated <= catalogue
    assert len(annotated) >= 5, annotated
    assert all(float(m["abundance"]) > 0 for m in models)


def test_the_signal_was_actually_scored(default_run):
    _out, log = default_run
    counters = [(int(t), int(r)) for t, r in re.findall(r"\(ties=(\d+) refined=(\d+)\)", log)]
    assert len(counters) == 2  # two intervals
    assert sum(t for t, _ in counters) >= 10
    assert sum(r for _, r in counters) >= 3  # contrasts that the signal scored


def test_the_format_and_pore_are_read_from_the_inputs(default_run):
    out, log = default_run
    config = json.loads((out / "run_manifest.json").read_text())["config"]
    assert config["signal_format"] == "slow5"
    assert config["krill_pore"] == "rna002"
    assert "pore=rna002 (sequencing kit 'sqk-rna002')" in log


@pytest.mark.parametrize("extra", [("--m2-batch-pairs", "1"), ("--threads", "2")])
def test_resource_settings_do_not_move_the_assembly(default_run, tmp_path, extra):
    out, _log = default_run
    result = _pyfin(tmp_path / "out", *extra)
    assert result.returncode == 0, result.stderr[-3000:]
    assert struct_compare.records(str(out / "assembly.gtf")) == struct_compare.records(
        str(tmp_path / "out" / "assembly.gtf")
    )


def test_it_runs_without_an_annotation(tmp_path):
    result = _pyfin(tmp_path / "out", gtf=False)
    assert result.returncode == 0, result.stderr[-3000:]
    models = _transcripts(tmp_path / "out" / "assembly.gtf")
    assert models and all(m["transcript_source"] != "gtf" for m in models)


def _copy(tmp_path: Path, *names: str) -> Path:
    data = tmp_path / "data"
    data.mkdir()
    for name in names:
        shutil.copy(DATA / name, data / name)
    return data


ALL = ("reads.bam", "reads.bam.bai", "reads.blow5", "reads.blow5.idx", "reference.fa",
       "reference.fa.fai", "annotation.gtf")


def test_a_bam_without_its_index_is_refused_with_the_fix(tmp_path):
    data = _copy(tmp_path, *(n for n in ALL if n != "reads.bam.bai"))
    result = _pyfin(tmp_path / "out", data=data)
    assert result.returncode == 2
    assert "has no index" in result.stderr and "samtools index" in result.stderr
    assert not (tmp_path / "out" / "assembly.gtf").exists()


def test_an_annotation_with_other_contig_names_is_refused(tmp_path):
    data = _copy(tmp_path, *ALL)
    gtf = data / "annotation.gtf"
    gtf.write_text(gtf.read_text().replace("SIRV1\t", "chrSIRV1\t"))
    result = _pyfin(tmp_path / "out", data=data)
    assert result.returncode == 2
    assert "share no contig name" in result.stderr and "chrSIRV1" in result.stderr


def test_a_genome_of_another_build_is_refused(tmp_path):
    data = _copy(tmp_path, *(n for n in ALL if n != "reference.fa.fai"))
    fasta = data / "reference.fa"
    lines = fasta.read_text().splitlines()
    fasta.write_text("\n".join(lines[:-1]) + "\n")  # the same contig, one line shorter
    result = _pyfin(tmp_path / "out", data=data)
    assert result.returncode == 2 and "different assembly builds" in result.stderr


def test_signal_of_other_reads_is_refused(tmp_path):
    pysam = pytest.importorskip("pysam")
    data = _copy(tmp_path, *(n for n in ALL if not n.startswith("reads.bam")))
    with pysam.AlignmentFile(str(DATA / "reads.bam")) as src, \
            pysam.AlignmentFile(str(data / "reads.bam"), "wb", header=src.header) as dst:
        for record in src:
            record.query_name = "renamed-" + record.query_name
            dst.write(record)
    pysam.index(str(data / "reads.bam"))
    result = _pyfin(tmp_path / "out", data=data)
    assert result.returncode == 2
    assert "no read of the BAM is present in the signal file" in result.stderr
    assert "different runs" in result.stderr and "renamed-" in result.stderr
    assert not (tmp_path / "out").exists()
