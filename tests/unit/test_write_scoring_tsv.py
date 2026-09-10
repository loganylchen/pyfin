"""Tests for fin.io.io_tsv.write_scoring_tsv (US-015 / AC-19)."""

from __future__ import annotations

import pytest

import os
import tempfile
from typing import Dict


from fin.analysis.quantification import QuantResult
from fin.io.io_tsv import write_scoring_tsv

_HEADER = (
    "candidate_id\tgene_id\tchrom\tstrand\tstart\tend\tsource\t"
    "abundance\tconfidence\tcoherence_score\tdiscrimination_score\t"
    "combined_score\tnum_reads\ttpm\tmax_R"
)


def _make_qr(
    candidate_id: str,
    abundance: float,
    source: str = "gtf",
    gene_id: str = "GENE1",
    chrom: str = "chr1",
    strand: str = "+",
    start: int = 1000,
    end: int = 2000,
    num_assigned_reads: int = 10,
    coherence_score: float = 0.8,
    discrimination_score: float = 0.7,
    combined_score: float = 0.75,
    confidence: float = 0.9,
) -> QuantResult:
    return QuantResult(
        candidate_id=candidate_id,
        abundance=abundance,
        confidence=confidence,
        num_assigned_reads=num_assigned_reads,
        source=source,
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
        gene_id=gene_id,
        coherence_score=coherence_score,
        discrimination_score=discrimination_score,
        combined_score=combined_score,
    )


def _write_and_read(results: Dict[str, QuantResult], lengths: Dict[str, int]):
    """Helper: write to a temp file, return (header_line, data_lines)."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
        path = tmp.name
    try:
        write_scoring_tsv(results, lengths, path)
        with open(path) as f:
            lines = f.read().splitlines()
    finally:
        os.unlink(path)
    return lines[0], lines[1:]


def test_write_scoring_tsv_basic():
    results = {
        "tx_gtf": _make_qr("tx_gtf", abundance=100.0, source="gtf"),
        "tx_novel": _make_qr("tx_novel", abundance=50.0, source="novel"),
        "tx_fusion": _make_qr("tx_fusion", abundance=25.0, source="fusion"),
    }
    lengths = {"tx_gtf": 1000, "tx_novel": 1000, "tx_fusion": 1000}

    header, data_lines = _write_and_read(results, lengths)

    # Header must match exactly
    assert header == _HEADER

    # 3 data rows
    assert len(data_lines) == 3

    # TPM column (index 13) sums to ~1e6
    tpm_sum = sum(float(line.split("\t")[13]) for line in data_lines)
    assert abs(tpm_sum - 1e6) < 1.0

    # Float columns use 4-decimal format
    for line in data_lines:
        fields = line.split("\t")
        for col_idx in (7, 8, 9, 10, 11, 13, 14):  # abundance..tpm, max_R
            val = fields[col_idx]
            assert "." in val
            assert len(val.split(".")[1]) == 4, f"Expected 4 decimals in '{val}'"


def test_write_scoring_tsv_empty_results():
    header, data_lines = _write_and_read({}, {})
    assert header == _HEADER
    assert data_lines == []


def test_write_scoring_tsv_zero_abundance():
    results = {"tx1": _make_qr("tx1", abundance=0.0)}
    lengths = {"tx1": 1000}

    header, data_lines = _write_and_read(results, lengths)
    assert len(data_lines) == 1
    fields = data_lines[0].split("\t")
    assert fields[7] == "0.0000"   # abundance
    assert fields[13] == "0.0000"  # tpm


def test_write_scoring_tsv_missing_length():
    results = {"tx_no_len": _make_qr("tx_no_len", abundance=50.0)}
    lengths: Dict[str, int] = {}  # no entry for tx_no_len

    header, data_lines = _write_and_read(results, lengths)
    assert len(data_lines) == 1
    fields = data_lines[0].split("\t")
    assert fields[13] == "0.0000"  # TPM falls back to 0


def test_write_scoring_tsv_column_order():
    results = {"tx1": _make_qr("tx1", abundance=10.0)}
    lengths = {"tx1": 500}

    header, _ = _write_and_read(results, lengths)
    cols = header.split("\t")
    expected = [
        "candidate_id", "gene_id", "chrom", "strand", "start", "end",
        "source", "abundance", "confidence", "coherence_score",
        "discrimination_score", "combined_score", "num_reads", "tpm", "max_R",
    ]
    assert cols == expected


def test_write_scoring_tsv_sorted_by_id():
    results = {
        "zzz_tx": _make_qr("zzz_tx", abundance=10.0),
        "aaa_tx": _make_qr("aaa_tx", abundance=20.0),
        "mmm_tx": _make_qr("mmm_tx", abundance=15.0),
    }
    lengths = {"zzz_tx": 1000, "aaa_tx": 1000, "mmm_tx": 1000}

    _, data_lines = _write_and_read(results, lengths)
    ids = [line.split("\t")[0] for line in data_lines]
    assert ids == sorted(ids)


class TestTpmConvention:
    """The direct-RNA TPM convention: no length normalisation.

    A whole molecule yields at most one read regardless of length, so the count
    does not carry the length factor that short-read normalisation removes, and
    dividing by length gives an INVERSE-length-weighted quantity. Measured
    against the Sequin spike-in's known input molarity, the un-normalised value
    agrees better (see analysis.quantification.compute_tpm).
    """

    def _results(self):
        from fin.analysis.quantification import QuantResult
        # same abundance, 10x different lengths
        return {
            "short": QuantResult("short", 10.0, 1.0, 10, "gtf", exons=((0, 1000),)),
            "long": QuantResult("long", 10.0, 1.0, 10, "gtf", exons=((0, 10000),)),
        }, {"short": 1000, "long": 10000}

    def test_length_normalized_favours_the_short_transcript(self):
        from fin.analysis.quantification import compute_tpm
        r, lens = self._results()
        t = compute_tpm(r, lens)
        assert t["short"] > t["long"]
        assert abs(t["short"] / t["long"] - 10.0) < 1e-6
        assert abs(sum(t.values()) - 1e6) < 1e-3

    def test_molar_convention_treats_equal_counts_equally(self):
        from fin.analysis.quantification import compute_tpm
        r, lens = self._results()
        t = compute_tpm(r, lens, length_normalize=False)
        assert abs(t["short"] - t["long"]) < 1e-6
        assert abs(sum(t.values()) - 1e6) < 1e-3

    def test_molar_convention_ignores_the_length_table_entirely(self):
        from fin.analysis.quantification import compute_tpm
        r, _ = self._results()
        assert (compute_tpm(r, {}, length_normalize=False)
                == compute_tpm(r, {"short": 5, "long": 999999}, length_normalize=False))

    def test_zero_total_abundance_is_not_a_division_error(self):
        from fin.analysis.quantification import QuantResult, compute_tpm
        r = {"a": QuantResult("a", 0.0, 0.0, 0, "gtf", exons=((0, 100),))}
        assert compute_tpm(r, {"a": 100}, length_normalize=False) == {"a": 0.0}

    def test_writer_threads_the_flag(self, tmp_path):
        import csv
        from fin.io.io_tsv import write_scoring_tsv
        r, lens = self._results()
        for norm, expect_equal in ((True, False), (False, True)):
            p = tmp_path / f"s{norm}.tsv"
            write_scoring_tsv(r, lens, str(p), length_normalize=norm)
            got = {row["candidate_id"]: float(row["tpm"])
                   for row in csv.DictReader(open(p), delimiter="\t")}
            assert (abs(got["short"] - got["long"]) < 1e-6) is expect_equal

    def test_non_finite_or_negative_abundance_is_rejected(self):
        from fin.analysis.quantification import QuantResult, compute_tpm
        for bad in (float("nan"), float("inf"), -1.0):
            r = {"a": QuantResult("a", bad, 0.0, 0, "gtf", exons=((0, 100),)),
                 "b": QuantResult("b", 10.0, 1.0, 10, "gtf", exons=((0, 100),))}
            with pytest.raises(ValueError, match="finite"):
                compute_tpm(r, {"a": 100, "b": 100}, length_normalize=False)

    def test_empty_result_set(self):
        from fin.analysis.quantification import compute_tpm
        assert compute_tpm({}, {}, length_normalize=False) == {}

    def test_refit_diagnostics_follow_the_same_convention(self):
        """`tpm` must not mean two different things in two files of one run."""
        import inspect
        from fin.analysis import abundance_refit as A
        assert "length_normalize" in inspect.signature(A.refit_survivor_abundance).parameters
        src = inspect.getsource(A.refit_survivor_abundance)
        assert src.count("length_normalize=length_normalize") == 2, (
            "both the before and after TPM in the diagnostics must honour the flag")
