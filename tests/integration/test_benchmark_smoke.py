"""Smoke tests for benchmark harness. SKIP when external tools absent."""
import pytest
import shutil
import subprocess
import sys
from pathlib import Path

COMPARE_PY = Path(__file__).resolve().parent.parent.parent / "benchmarks" / "compare_results.py"


# The two harness smoke tests that lived here were retired on 2026-09-21 with
# benchmarks/run_benchmark.sh itself. That script wrote {"status": "ran"}
# markers without ever invoking a tool, and the tests covered its help output,
# exit status and [SKIP]/[RUN] messages -- real assertions about the skeleton,
# but never about whether any external tool ran. Retiring them alongside the
# script loses no tool-execution coverage. The compare_results.py tests below
# exercise a utility that does real work and are kept.


def test_compare_results_produces_tsv(tmp_path):
    """compare_results.py writes a valid TSV with header and tool rows."""
    if not COMPARE_PY.exists():
        pytest.skip("compare_results.py not found")

    input_dir = tmp_path / "in"
    input_dir.mkdir()

    # Create minimal per-tool output
    tool_dir = input_dir / "pyfin"
    tool_dir.mkdir()
    (tool_dir / "result.json").write_text('{"status": "ran"}')

    out = tmp_path / "cmp.tsv"
    r = subprocess.run(
        [sys.executable, str(COMPARE_PY), "--input-dir", str(input_dir), "--output", str(out)],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    assert out.exists()
    content = out.read_text()
    assert "tool\tmetric\tvalue" in content
    assert "pyfin" in content


def test_compare_results_handles_empty_input_dir(tmp_path):
    """compare_results.py exits 0 and writes header-only TSV for empty input dir."""
    if not COMPARE_PY.exists():
        pytest.skip("compare_results.py not found")

    input_dir = tmp_path / "empty"
    input_dir.mkdir()
    out = tmp_path / "empty.tsv"

    r = subprocess.run(
        [sys.executable, str(COMPARE_PY), "--input-dir", str(input_dir), "--output", str(out)],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    assert out.exists()
    content = out.read_text()
    assert "tool\tmetric\tvalue" in content


def test_compare_results_counts_gtf_transcripts(tmp_path):
    """num_transcripts metric is populated from *.gtf files."""
    if not COMPARE_PY.exists():
        pytest.skip("compare_results.py not found")

    input_dir = tmp_path / "in"
    tool_dir = input_dir / "pyfin"
    tool_dir.mkdir(parents=True)
    (tool_dir / "result.json").write_text('{"status": "ran"}')

    # Write a minimal GTF with 2 transcript lines
    gtf_lines = [
        "# comment line",
        'chr1\t.\ttranscript\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
        'chr1\t.\texon\t100\t200\t.\t+\t.\tgene_id "G1"; transcript_id "T1";',
        'chr1\t.\ttranscript\t300\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T2";',
    ]
    (tool_dir / "output.gtf").write_text("\n".join(gtf_lines) + "\n")

    out = tmp_path / "cmp.tsv"
    r = subprocess.run(
        [sys.executable, str(COMPARE_PY), "--input-dir", str(input_dir), "--output", str(out)],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    content = out.read_text()
    assert "num_transcripts\t2" in content
