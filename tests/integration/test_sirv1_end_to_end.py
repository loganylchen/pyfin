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
    out, log = default_run
    counters = [(int(t), int(r)) for t, r in re.findall(r"\(ties=(\d+) refined=(\d+)\)", log)]
    assert len(counters) == 2  # two intervals
    assert sum(t for t, _ in counters) >= 10
    assert sum(r for _, r in counters) >= 3  # contrasts that the signal scored
    # and the same is on record without reading a log
    signal = json.loads((out / "run_summary.json").read_text())["signal"]
    assert signal == {"scoring": "summed_llr", "pore": "rna002",
                      "tied_reads": sum(t for t, _ in counters),
                      "scored_contrasts": sum(r for _, r in counters)}


def test_the_summary_describes_the_run(default_run):
    out, _log = default_run
    summary = json.loads((out / "run_summary.json").read_text())
    models = _transcripts(out / "assembly.gtf")
    assert summary["models"] == {
        "total": len(models),
        "annotated": sum(m["transcript_source"] == "gtf" for m in models),
        "novel": sum(m["transcript_source"] != "gtf" for m in models),
    }
    assert summary["intervals"] == {"total": 2, "assigned": 2}
    assert 0 < summary["reads"]["assigned"] <= summary["reads"]["in_assigned_intervals"] == 91
    assert summary["quant_mode"] == "m2_em"
    assert summary["warnings"] == [] and summary["peak_rss_gib"] > 0 and summary["elapsed_seconds"] > 0
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert manifest["libraries"] == summary["libraries"] and manifest["libraries"]["krill"]
    assert manifest["pyfin_version"] == summary["pyfin_version"]


def test_the_output_directory_holds_only_results(default_run):
    out, _log = default_run
    assert sorted(p.name for p in out.iterdir() if p.is_dir()) == []
    assert {"assembly.gtf", "scores.tsv", "run_manifest.json", "run_summary.json"} <= {
        p.name for p in out.iterdir()}


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
    # workers write their interval counters as shards; the totals must not depend on who ran what
    summaries = [json.loads((d / "run_summary.json").read_text()) for d in (out, tmp_path / "out")]
    for key in ("intervals", "reads", "signal", "models"):
        assert summaries[0][key] == summaries[1][key]


def _scores(path: Path):
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def _spliced_lengths(gtf: Path):
    lengths = {}
    for line in gtf.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) > 8 and fields[2] == "exon":
            tid = re.search(r'transcript_id "([^"]*)"', fields[8]).group(1)
            lengths[tid] = lengths.get(tid, 0) + int(fields[4]) - int(fields[3]) + 1
    return lengths


def _all_but_tpm(rows):
    """Every column except tpm, with the spelling of minted novel ids masked."""
    return sorted(
        tuple((k, struct_compare._mask(v)) for k, v in row.items() if k != "tpm") for row in rows
    )


def test_tpm_is_a_share_of_abundance_and_length_normalising_changes_only_that_column(default_run, tmp_path):
    out, _log = default_run
    rows = _scores(out / "scores.tsv")
    total = sum(float(r["abundance"]) for r in rows)
    for row in rows:
        assert float(row["tpm"]) == pytest.approx(1e6 * float(row["abundance"]) / total, rel=1e-4)

    other = tmp_path / "out"
    result = _pyfin(other, "--tpm-length-normalize")
    assert result.returncode == 0, result.stderr[-3000:]
    assert struct_compare.records(str(out / "assembly.gtf")) == struct_compare.records(
        str(other / "assembly.gtf"))
    by_length = _scores(other / "scores.tsv")
    assert _all_but_tpm(rows) == _all_but_tpm(by_length)        # nothing but tpm moved

    lengths = _spliced_lengths(other / "assembly.gtf")
    density = {r["candidate_id"]: float(r["abundance"]) / lengths[r["candidate_id"]] for r in by_length}
    for row in by_length:
        expected = 1e6 * density[row["candidate_id"]] / sum(density.values())
        assert float(row["tpm"]) == pytest.approx(expected, rel=1e-3)
    assert sorted(float(r["tpm"]) for r in rows) != sorted(float(r["tpm"]) for r in by_length)


def test_an_engine_that_does_not_score_ties_reports_its_signal_counters_as_unknown(tmp_path):
    result = _pyfin(tmp_path / "out", "--quant-mode", "argmax")
    assert result.returncode == 0, result.stderr[-3000:]
    summary = json.loads((tmp_path / "out" / "run_summary.json").read_text())
    assert summary["quant_mode"] == "argmax"
    assert summary["intervals"] == {"total": 2, "assigned": 2}
    assert summary["reads"] == {"in_assigned_intervals": 91, "assigned": None}
    assert summary["signal"]["tied_reads"] is None and summary["signal"]["scored_contrasts"] is None
    assert summary["models"]["total"] == len(_transcripts(tmp_path / "out" / "assembly.gtf")) > 0


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
