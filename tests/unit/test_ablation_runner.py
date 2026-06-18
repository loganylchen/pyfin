"""Tests for fin/ablation/runner.py and fin/ablation/io.py.

Covers:
  - ABLATION_ROWS definitions (row_id, label, field values per AC1)
  - AblationRowConfig defaults
  - _nan_to_row_mean: AC8 NaN adapter
  - _apply_score_filters: AC7 score-filter bypass
  - write_ablation_summary / write_per_row_tsv: IO correctness
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path
from typing import Dict, List

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.ablation.runner import (
    ABLATION_ROWS,
    AblationResult,
    AblationRowConfig,
    _apply_score_filters,
    _nan_to_row_mean,
)
from fin.analysis.quantification import QuantResult
from fin.pipeline.config import PipelineConfig


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _qr(cid: str, source: str = "novel", abundance: float = 1.0,
        max_R: float = 0.5, combined_score: float = 0.5) -> QuantResult:
    return QuantResult(
        candidate_id=cid,
        abundance=abundance,
        confidence=1.0,
        num_assigned_reads=int(abundance),
        source=source,
        max_R=max_R,
        combined_score=combined_score,
    )


def _ablation_result(row_id: str = "R3", label: str = "test",
                     aggregated: Dict[str, QuantResult] = None) -> AblationResult:
    return AblationResult(
        row_id=row_id,
        label=label,
        aggregated=aggregated or {},
    )


# ---------------------------------------------------------------------------
# ABLATION_ROWS definitions (AC1)
# ---------------------------------------------------------------------------

class TestAblationRowsDefinitions:
    def test_seven_rows_defined(self):
        assert len(ABLATION_ROWS) == 7

    def test_row_ids_ordered(self):
        ids = [r.row_id for r in ABLATION_ROWS]
        assert ids == ["R1a", "R1b", "R1c", "R2", "R3", "R4", "R5"]

    def test_r1a_argmax_unfiltered(self):
        r1a = next(r for r in ABLATION_ROWS if r.row_id == "R1a")
        assert r1a.quant_mode == "argmax"
        assert r1a.m4_source == "none"
        assert r1a.em_max_iter_override is None
        assert r1a.enable_score_filter is False

    def test_r1b_argmax_filtered(self):
        r1b = next(r for r in ABLATION_ROWS if r.row_id == "R1b")
        assert r1b.quant_mode == "argmax"
        assert r1b.m4_source == "none"
        assert r1b.enable_score_filter is True

    def test_r1c_mappy_em_filtered(self):
        r1c = next(r for r in ABLATION_ROWS if r.row_id == "R1c")
        assert r1c.quant_mode == "m1_em"
        assert r1c.m4_source == "none"
        assert r1c.em_max_iter_override is None
        assert r1c.enable_score_filter is True

    def test_r2_single_step_em(self):
        r2 = next(r for r in ABLATION_ROWS if r.row_id == "R2")
        assert r2.quant_mode == "m2_em"
        assert r2.em_max_iter_override == 1
        # R2 isolates "single-step EM" as the only variable vs R4 (coherence on).
        assert r2.m4_source == "diff_region"
        assert r2.m3_coherence is True
        assert r2.enable_score_filter is False

    def test_r3_no_coherence(self):
        """AC5: R3 disables coherence (m3_coherence=False -> EM beta=0)."""
        r3 = next(r for r in ABLATION_ROWS if r.row_id == "R3")
        assert r3.quant_mode == "m2_em"
        assert r3.m4_source == "none"
        assert r3.m3_coherence is False
        assert r3.em_max_iter_override is None
        assert r3.enable_score_filter is False

    def test_r4_diff_region(self):
        r4 = next(r for r in ABLATION_ROWS if r.row_id == "R4")
        assert r4.quant_mode == "m2_em"
        assert r4.m4_source == "diff_region"
        assert r4.m3_coherence is True
        assert r4.enable_score_filter is False

    def test_r5_diff_region_scored(self):
        r5 = next(r for r in ABLATION_ROWS if r.row_id == "R5")
        assert r5.quant_mode == "m2_em"
        assert r5.m4_source == "diff_region"
        assert r5.m3_coherence is True
        assert r5.enable_score_filter is True


# ---------------------------------------------------------------------------
# AblationRowConfig defaults
# ---------------------------------------------------------------------------

class TestAblationRowConfigDefaults:
    def test_defaults(self):
        cfg = AblationRowConfig(row_id="X", label="test")
        assert cfg.quant_mode == "m2_em"
        assert cfg.em_max_iter_override is None
        assert cfg.m4_source == "diff_region"
        assert cfg.m3_coherence is False
        assert cfg.enable_score_filter is True


# ---------------------------------------------------------------------------
# _nan_to_row_mean (AC8)
# ---------------------------------------------------------------------------

class TestNanToRowMean:
    def test_no_nans_unchanged(self):
        m = np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32)
        out = _nan_to_row_mean(m)
        np.testing.assert_array_equal(out, m)

    def test_nan_replaced_by_row_mean(self):
        m = np.array([[1.0, np.nan, 3.0]], dtype=np.float32)
        out = _nan_to_row_mean(m)
        # Row mean of non-NaN values: (1+3)/2 = 2.0
        assert out[0, 1] == pytest.approx(2.0)
        assert out[0, 0] == 1.0
        assert out[0, 2] == 3.0

    def test_all_nan_row_becomes_zeros(self):
        m = np.array([[np.nan, np.nan]], dtype=np.float32)
        out = _nan_to_row_mean(m)
        assert np.all(out == 0.0)

    def test_partial_nan_multiple_rows(self):
        m = np.array([
            [1.0, np.nan],
            [np.nan, np.nan],
            [2.0, 4.0],
        ], dtype=np.float32)
        out = _nan_to_row_mean(m)
        assert out[0, 1] == pytest.approx(1.0)   # row mean = 1.0
        assert out[1, 0] == 0.0                   # all-NaN → 0
        assert out[1, 1] == 0.0
        assert out[2, 0] == 2.0
        assert out[2, 1] == 4.0

    def test_input_not_mutated(self):
        m = np.array([[np.nan, 1.0]], dtype=np.float32)
        _nan_to_row_mean(m)
        assert np.isnan(m[0, 0])  # original unchanged


# ---------------------------------------------------------------------------
# _apply_score_filters (AC7)
# ---------------------------------------------------------------------------

class TestApplyScoreFilters:
    def _config(self, min_ab=0.0):
        return PipelineConfig(
            bam_path="/tmp/x.bam",
            min_abundance=min_ab,
        )

    def test_no_filters_returns_all(self):
        agg = {"c1": _qr("c1"), "c2": _qr("c2")}
        cfg = self._config()
        out = _apply_score_filters(agg, cfg)
        assert set(out) == {"c1", "c2"}

    def test_min_abundance_drops_novel(self):
        agg = {
            "n1": _qr("n1", source="novel", abundance=0.2),
            "g1": _qr("g1", source="gtf", abundance=0.2),
        }
        cfg = self._config(min_ab=0.5)
        out = _apply_score_filters(agg, cfg)
        # GTF exempt, novel dropped
        assert "g1" in out
        assert "n1" not in out

    def test_min_abundance_keeps_novel_above_threshold(self):
        agg = {"n1": _qr("n1", source="novel", abundance=1.0)}
        cfg = self._config(min_ab=0.5)
        out = _apply_score_filters(agg, cfg)
        assert "n1" in out

    def test_fusion_exempt_from_filters(self):
        agg = {"f1": _qr("f1", source="fusion", abundance=0.0)}
        cfg = self._config(min_ab=1.0)
        out = _apply_score_filters(agg, cfg)
        assert "f1" in out


# ---------------------------------------------------------------------------
# IO: write_ablation_summary and write_per_row_tsv
# ---------------------------------------------------------------------------

class TestAblationIO:
    def _make_result(self, row_id: str, label: str) -> AblationResult:
        return _ablation_result(
            row_id=row_id,
            label=label,
            aggregated={
                "tx1": _qr("tx1", source="gtf", abundance=3.5),
                "tx2": _qr("tx2", source="novel", abundance=1.2),
            },
        )

    def test_write_ablation_summary_creates_file(self):
        from fin.ablation.io import write_ablation_summary

        results = [self._make_result("R1", "mappy_argmax_only")]
        with tempfile.TemporaryDirectory() as tmp:
            path = write_ablation_summary(results, tmp)
            assert Path(path).exists()

    def test_write_ablation_summary_header(self):
        from fin.ablation.io import write_ablation_summary
        import csv

        results = [self._make_result("R3", "em_no_coherence")]
        with tempfile.TemporaryDirectory() as tmp:
            path = write_ablation_summary(results, tmp)
            with open(path, newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)
        assert len(rows) == 2  # tx1 + tx2
        assert rows[0]["row_id"] == "R3"
        assert rows[0]["label"] == "em_no_coherence"
        # Sorted by candidate_id
        assert rows[0]["candidate_id"] == "tx1"
        assert rows[1]["candidate_id"] == "tx2"

    def test_write_ablation_summary_multiple_rows(self):
        from fin.ablation.io import write_ablation_summary
        import csv

        results = [
            self._make_result("R1", "mappy_argmax_only"),
            self._make_result("R2", "single_step_em"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = write_ablation_summary(results, tmp)
            with open(path, newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)
        # 2 transcripts × 2 rows
        assert len(rows) == 4
        row_ids = {r["row_id"] for r in rows}
        assert row_ids == {"R1", "R2"}

    def test_write_per_row_tsv_creates_file(self):
        from fin.ablation.io import write_per_row_tsv

        result = self._make_result("R4", "diff_region_dtw")
        with tempfile.TemporaryDirectory() as tmp:
            path = write_per_row_tsv(result, tmp)
            assert Path(path).exists()
            assert "R4" in Path(path).name
            assert "diff_region_dtw" in Path(path).name

    def test_write_per_row_tsv_content(self):
        from fin.ablation.io import write_per_row_tsv
        import csv

        result = self._make_result("R5", "diff_region_dtw_scored")
        with tempfile.TemporaryDirectory() as tmp:
            path = write_per_row_tsv(result, tmp)
            with open(path, newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)
        assert len(rows) == 2
        assert rows[0]["candidate_id"] == "tx1"
        assert float(rows[0]["abundance"]) == pytest.approx(3.5, abs=1e-3)

    def test_write_ablation_summary_empty_results(self):
        from fin.ablation.io import write_ablation_summary

        with tempfile.TemporaryDirectory() as tmp:
            path = write_ablation_summary([], tmp)
            assert Path(path).exists()
            lines = Path(path).read_text().splitlines()
            # Only header line
            assert len(lines) == 1
