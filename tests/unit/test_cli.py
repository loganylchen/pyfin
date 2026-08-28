"""Tests for fin CLI (assembly-only; fusion and quantify subcommands removed)."""

from __future__ import annotations

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
        "--profile",
        "--no-gpu",
        "--fusion",
        "--m2-metric",
        "--m2-summed-llr-margin",
        "--strict-novel-abundance-floor",
        "--isoform-fraction-locus",
        "--post-selection-refit",
    ]:
        assert flag in r.stdout, f"Missing flag in fin --help: {flag}"


def test_isoform_fraction_locus_rejects_unknown_mode():
    r = run_cli("--isoform-fraction-locus", "gene")
    assert r.returncode != 0
    assert "Invalid value" in r.stderr


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


def test_default_help_shows_junction_snap_flags():
    r = run_cli("--help")
    assert r.returncode == 0
    for flag in [
        "--junction-snap",
        "--junction-snap-tolerance",
        "--junction-snap-min-support",
        "--junction-snap-min-ratio",
    ]:
        assert flag in r.stdout, f"Missing junction-snap flag in fin --help: {flag}"


def test_default_help_shows_junction_support_flags():
    r = run_cli("--help")
    assert r.returncode == 0
    for flag in ["--novel-junction-min-reads", "--novel-junction-reads-tol"]:
        assert flag in r.stdout, f"Missing junction-support flag in fin --help: {flag}"


def test_default_help_shows_junction_dominance_flags():
    r = run_cli("--help")
    assert r.returncode == 0
    for flag in [
        "--junction-dominance-filter",
        "--junction-dominance-min-reads",
        "--junction-dominance-window-bp",
        "--junction-dominance-tol-bp",
    ]:
        assert flag in r.stdout, f"Missing junction-dominance flag in fin --help: {flag}"


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

def _required_args(tmp_path):
    paths = []
    for flag, name in (
        ("--bam", "mapped.bam"),
        ("--genome", "SIRV.genome.fa"),
        ("--fastq", "mapped.fq.gz"),
        ("--signal", "mapped.blow5"),
    ):
        path = tmp_path / name
        path.touch()
        paths.extend((flag, str(path)))
    paths.extend(("--output-dir", str(tmp_path / "out")))
    return paths


def test_help_shows_parallel_flags():
    r = run_cli("--help")
    assert r.returncode == 0
    for flag in ["--threads", "--gpu-workers"]:
        assert flag in r.stdout, f"Missing parallel flag in fin --help: {flag}"


def test_threads_zero_rejected(tmp_path):
    # Validation fires before the pipeline runs -> fast, non-zero exit.
    r = run_cli(*_required_args(tmp_path), "--threads", "0")
    assert r.returncode != 0
    assert "threads" in (r.stdout + r.stderr).lower()


def test_gpu_workers_negative_rejected(tmp_path):
    r = run_cli(*_required_args(tmp_path), "--threads", "2", "--gpu-workers=-1")
    assert r.returncode != 0
    assert "gpu-workers" in (r.stdout + r.stderr).lower()


def test_gpu_workers_exceeds_threads_rejected(tmp_path):
    r = run_cli(
        *_required_args(tmp_path),
        "--gpu", "--threads", "2", "--gpu-workers", "3",
    )
    assert r.returncode != 0
    assert "gpu-workers" in (r.stdout + r.stderr).lower()
