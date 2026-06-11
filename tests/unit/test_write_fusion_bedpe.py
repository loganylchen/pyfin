"""Tests for fin.io.io_bedpe.write_fusion_bedpe (AC-20)."""

from __future__ import annotations

import os
import tempfile


from fin.analysis.quantification import QuantResult
from fin.io.io_bedpe import write_fusion_bedpe


def _make_fusion(
    candidate_id: str,
    num_assigned_reads: int = 10,
    breakpoint_left: tuple | None = ("chr1", 1000, "+"),
    breakpoint_right: tuple | None = ("chr2", 2000, "-"),
) -> QuantResult:
    return QuantResult(
        candidate_id=candidate_id,
        abundance=float(num_assigned_reads),
        confidence=0.9,
        num_assigned_reads=num_assigned_reads,
        source="fusion",
        breakpoint_left=breakpoint_left,
        breakpoint_right=breakpoint_right,
    )


def _make_gtf(candidate_id: str) -> QuantResult:
    return QuantResult(
        candidate_id=candidate_id,
        abundance=5.0,
        confidence=0.8,
        num_assigned_reads=5,
        source="gtf",
        combined_score=0.7,
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_and_read(results):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".bedpe", delete=False) as f:
        path = f.name
    try:
        write_fusion_bedpe(results, path)
        with open(path) as fh:
            lines = [ln for ln in fh.read().splitlines() if ln]
        return lines
    finally:
        os.unlink(path)


# ---------------------------------------------------------------------------
# AC-20 tests
# ---------------------------------------------------------------------------

def test_bedpe_fusion_only():
    """2 fusion + 1 gtf results → exactly 2 output lines."""
    results = {
        "f1": _make_fusion("f1"),
        "f2": _make_fusion("f2", breakpoint_left=("chr3", 500, "+"), breakpoint_right=("chr4", 600, "+")),
        "g1": _make_gtf("g1"),
    }
    lines = _write_and_read(results)
    assert len(lines) == 2


def test_bedpe_column_count():
    """Each output line has exactly 10 tab-separated fields."""
    results = {"f1": _make_fusion("f1")}
    lines = _write_and_read(results)
    assert len(lines) == 1
    fields = lines[0].split("\t")
    assert len(fields) == 10


def test_bedpe_score_is_read_support():
    """score column == num_assigned_reads (the fusion's read evidence)."""
    results = {"f1": _make_fusion("f1", num_assigned_reads=7)}
    lines = _write_and_read(results)
    fields = lines[0].split("\t")
    assert fields[7] == "7"


def test_bedpe_score_clamping():
    """num_assigned_reads > 1000 clamps to 1000; 0 stays 0."""
    results = {
        "high": _make_fusion("high", num_assigned_reads=5000),
        "low": _make_fusion(
            "low",
            num_assigned_reads=0,
            breakpoint_left=("chrZ", 9999, "+"),
            breakpoint_right=("chrZ", 10000, "+"),
        ),
    }
    lines = _write_and_read(results)
    scores = {ln.split("\t")[6]: ln.split("\t")[7] for ln in lines}
    assert scores["high"] == "1000"
    assert scores["low"] == "0"


def test_bedpe_empty_no_fusion():
    """Only non-fusion results → output file has 0 lines."""
    results = {"g1": _make_gtf("g1"), "g2": _make_gtf("g2")}
    lines = _write_and_read(results)
    assert lines == []


def test_bedpe_missing_breakpoint_skipped():
    """Fusion entry with breakpoint_left=None is skipped without crash."""
    results = {
        "bad": _make_fusion("bad", breakpoint_left=None),
        "good": _make_fusion("good"),
    }
    lines = _write_and_read(results)
    assert len(lines) == 1
    assert lines[0].split("\t")[6] == "good"


def test_bedpe_sorted_deterministic():
    """Two fusions in reverse order → output sorted by (chromA, startA)."""
    results = {
        "b": _make_fusion("b", breakpoint_left=("chrB", 5000, "+"), breakpoint_right=("chrC", 100, "+")),
        "a": _make_fusion("a", breakpoint_left=("chrA", 100, "+"), breakpoint_right=("chrC", 200, "+")),
    }
    lines = _write_and_read(results)
    assert len(lines) == 2
    first_chrom = lines[0].split("\t")[0]
    second_chrom = lines[1].split("\t")[0]
    assert first_chrom == "chrA"
    assert second_chrom == "chrB"


def test_bedpe_coordinates():
    """startA=posA, endA=posA+1; same for B."""
    qr = _make_fusion("f1", breakpoint_left=("chr1", 999, "+"), breakpoint_right=("chr2", 1999, "-"))
    lines = _write_and_read({"f1": qr})
    fields = lines[0].split("\t")
    assert fields[0] == "chr1"
    assert fields[1] == "999"
    assert fields[2] == "1000"
    assert fields[3] == "chr2"
    assert fields[4] == "1999"
    assert fields[5] == "2000"
    assert fields[8] == "+"
    assert fields[9] == "-"


def test_bedpe_no_header():
    """Output file must not contain a header line."""
    results = {"f1": _make_fusion("f1")}
    lines = _write_and_read(results)
    for line in lines:
        assert not line.lower().startswith("chrom")


def test_bedpe_accepts_iterable():
    """write_fusion_bedpe accepts an iterable (not only a dict)."""
    qrs = [_make_fusion("f1"), _make_gtf("g1")]
    lines = _write_and_read(qrs)
    assert len(lines) == 1
