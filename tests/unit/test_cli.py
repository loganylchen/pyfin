"""Tests for fin CLI (default assembly command + quantify subcommand)."""

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
        "--use-prior",
        "--no-prior",
        "--no-gpu",
        "--alpha",
        "--fusion",
    ]:
        assert flag in r.stdout, f"Missing flag in fin --help: {flag}"


def test_default_help_shows_fusion_flags():
    r = run_cli("--help")
    assert r.returncode == 0
    for flag in ["--min-support", "--max-dist", "--flank-bp"]:
        assert flag in r.stdout, f"Missing fusion flag in fin --help: {flag}"


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
# quantify subcommand backward compat
# ---------------------------------------------------------------------------


def test_quantify_help():
    r = run_cli("quantify", "--help")
    assert r.returncode == 0, f"Expected exit 0, got {r.returncode}\n{r.stderr}"
    for flag in ["--gtf", "--genome", "--sample", "--output-dir", "--use-gpu", "--no-gpu"]:
        assert flag in r.stdout, f"Missing existing flag in quantify --help: {flag}"
    for flag in ["--use-prior", "--no-prior"]:
        assert flag in r.stdout, f"Missing flag in quantify --help: {flag}"


# ---------------------------------------------------------------------------
# General CLI behaviour
# ---------------------------------------------------------------------------


def test_no_args_shows_help_or_error():
    r = run_cli()
    # No args → missing required flags for default assembly → non-zero exit
    assert r.returncode != 0 or "usage" in (r.stdout + r.stderr).lower(), (
        f"Expected non-zero exit or usage text; got returncode={r.returncode}"
    )


def test_top_level_help_lists_quantify():
    r = run_cli("--help")
    assert r.returncode == 0
    assert "quantify" in r.stdout, "Subcommand 'quantify' missing from top-level --help"


def test_unknown_subcommand_errors():
    r = run_cli("nonexistent-command")
    assert r.returncode != 0
