"""Tests for fin CLI (assembly-only; fusion and quantify subcommands removed)."""

from __future__ import annotations

import os
import subprocess
import sys


def run_cli(*args):
    """Invoke the fin CLI via python -m fin.cli."""
    return subprocess.run(
        [sys.executable, "-m", "fin.cli", *args],
        capture_output=True,
        text=True,
    )


# ---------------------------------------------------------------------------
# Default command (assembly) — fin --help exposes assembly + fusion flags
# ---------------------------------------------------------------------------


def test_default_help_shows_assembly_flags():
    r = run_cli("--help")
    assert r.returncode == 0, f"Expected exit 0, got {r.returncode}\n{r.stderr}"
    for flag in [
        "--bam",
        "--gtf",
        "--genome",
        "--fastq",
        "--signal",
        "--output-dir",
        "--no-gpu",
        "--fusion",
    ]:
        assert flag in r.stdout, f"Missing flag in fin --help: {flag}"


def test_default_help_shows_fusion_flags():
    r = run_cli("--help")
    assert r.returncode == 0
    for flag in ["--min-support", "--max-dist", "--flank-bp"]:
        assert flag in r.stdout, f"Missing fusion flag in fin --help: {flag}"


def test_default_help_shows_containment_flags():
    r = run_cli("--help")
    assert r.returncode == 0
    for flag in [
        "--containment-collapse",
        "--containment-3p-tol-bp",
        "--containment-min-abundance-ratio",
    ]:
        assert flag in r.stdout, f"Missing containment flag in fin --help: {flag}"


def test_default_help_shows_mono_exon_flags():
    r = run_cli("--help")
    assert r.returncode == 0
    for flag in [
        "--drop-mono-exon-novel",
        "--min-mono-exon-reads",
        "--min-mono-exon-length",
    ]:
        assert flag in r.stdout, f"Missing mono-exon flag in fin --help: {flag}"


def test_default_help_shows_junction_support_flags():
    r = run_cli("--help")
    assert r.returncode == 0
    for flag in ["--novel-junction-min-reads", "--novel-junction-reads-tol"]:
        assert flag in r.stdout, f"Missing junction-support flag in fin --help: {flag}"


def test_default_help_shows_description():
    r = run_cli("--help")
    assert r.returncode == 0
    assert "assembly" in r.stdout.lower() or "assemble" in r.stdout.lower()


def test_assemble_subcommand_removed():
    r = run_cli("assemble", "--help")
    assert r.returncode != 0


def test_fusion_subcommand_removed():
    r = run_cli("fusion", "--help")
    assert r.returncode != 0


# ---------------------------------------------------------------------------
# quantify subcommand removed (pyfin is assembly-only)
# ---------------------------------------------------------------------------


def test_quantify_subcommand_removed():
    r = run_cli("quantify", "--help")
    assert r.returncode != 0


# ---------------------------------------------------------------------------
# General CLI behaviour
# ---------------------------------------------------------------------------


def test_no_args_shows_help_or_error():
    r = run_cli()
    # No args → missing required flags for default assembly → non-zero exit
    assert r.returncode != 0 or "usage" in (r.stdout + r.stderr).lower(), (
        f"Expected non-zero exit or usage text; got returncode={r.returncode}"
    )


def test_top_level_help_omits_quantify():
    r = run_cli("--help")
    assert r.returncode == 0
    assert "quantify" not in r.stdout, "Removed subcommand 'quantify' still in top-level --help"


def test_unknown_subcommand_errors():
    r = run_cli("nonexistent-command")
    assert r.returncode != 0


# ---------------------------------------------------------------------------
# Interval-level parallelism flags (--threads / --gpu-workers)
# ---------------------------------------------------------------------------

_TD = os.path.join(os.path.dirname(__file__), "..", "..", "testdata")
_REQUIRED = [
    "--bam", os.path.join(_TD, "mapped.bam"),
    "--genome", os.path.join(_TD, "SIRV.genome.fa"),
    "--fastq", os.path.join(_TD, "mapped.fq.gz"),
    "--signal", os.path.join(_TD, "mapped.blow5"),
    "--output-dir", "/tmp/_fin_cli_validate",
]


def test_help_shows_parallel_flags():
    r = run_cli("--help")
    assert r.returncode == 0
    for flag in ["--threads", "--gpu-workers"]:
        assert flag in r.stdout, f"Missing parallel flag in fin --help: {flag}"


def test_threads_zero_rejected():
    # Validation fires before the pipeline runs -> fast, non-zero exit.
    r = run_cli(*_REQUIRED, "--threads", "0")
    assert r.returncode != 0
    assert "threads" in (r.stdout + r.stderr).lower()


def test_gpu_workers_negative_rejected():
    r = run_cli(*_REQUIRED, "--threads", "2", "--gpu-workers=-1")
    assert r.returncode != 0
    assert "gpu-workers" in (r.stdout + r.stderr).lower()


def test_gpu_workers_exceeds_threads_rejected():
    r = run_cli(*_REQUIRED, "--gpu", "--threads", "2", "--gpu-workers", "3")
    assert r.returncode != 0
    assert "gpu-workers" in (r.stdout + r.stderr).lower()
