"""Tests for PipelineRunner.run() output writers (US-013)."""

from __future__ import annotations

from typing import Dict
from unittest.mock import patch

from fin.analysis.quantification import QuantResult
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_quant_result(cid: str, exons=((0, 100),)) -> QuantResult:
    """Return a minimal QuantResult sufficient for output writers."""
    return QuantResult(
        candidate_id=cid,
        gene_id=cid,
        chrom="chr1",
        strand="+",
        start=0,
        end=100,
        source="gtf",
        abundance=1.0,
        confidence=1.0,
        coherence_score=0.5,
        discrimination_score=0.5,
        combined_score=0.5,
        num_assigned_reads=10,
        exons=exons,
    )


def _make_aggregated() -> Dict[str, QuantResult]:
    return {"tx1": _make_quant_result("tx1", exons=((0, 100), (200, 300)))}


def _make_runner(
    output_gtf=None,
    output_tsv=None,
    output_bedpe=None,
    fusion_enabled=False,
) -> PipelineRunner:
    config = PipelineConfig(
        bam_path="/fake/reads.bam",
        output_gtf=output_gtf,
        output_tsv=output_tsv,
        output_bedpe=output_bedpe,
        fusion_enabled=fusion_enabled,
    )
    runner = PipelineRunner(config)
    return runner


def _patch_run_internals(runner: PipelineRunner, aggregated: Dict[str, QuantResult]):
    """Patch out everything in run() except the output section."""
    runner_mod = "fin.pipeline.runner"

    patches = [
        patch(f"{runner_mod}.generate_isolated_intervals", return_value={"intervals": []}),
        patch(f"{runner_mod}.aggregate_across_intervals", return_value=aggregated),
    ]
    return patches


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestScoringTsv:
    def test_run_calls_scoring_tsv_when_configured(self, tmp_path):
        tsv_path = str(tmp_path / "out.tsv")
        runner = _make_runner(output_tsv=tsv_path)
        aggregated = _make_aggregated()

        with patch("fin.pipeline.runner.generate_isolated_intervals",
                   return_value={"intervals": []}), \
             patch("fin.pipeline.runner.aggregate_across_intervals",
                   return_value=aggregated), \
             patch("fin.io.io_tsv.write_scoring_tsv") as mock_tsv:

            runner.run()

        mock_tsv.assert_called_once()
        call_args = mock_tsv.call_args
        assert call_args[0][0] == aggregated
        # transcript_lengths: tx1 has exons (0,100),(200,300) -> 100+100=200
        assert call_args[0][1] == {"tx1": 200}
        assert call_args[0][2] == tsv_path

    def test_run_skips_scoring_tsv_when_none(self):
        runner = _make_runner(output_tsv=None)
        aggregated = _make_aggregated()

        with patch("fin.pipeline.runner.generate_isolated_intervals",
                   return_value={"intervals": []}), \
             patch("fin.pipeline.runner.aggregate_across_intervals",
                   return_value=aggregated), \
             patch("fin.io.io_tsv.write_scoring_tsv") as mock_tsv:

            runner.run()

        mock_tsv.assert_not_called()


class TestFusionBedpe:
    def test_run_calls_fusion_bedpe_when_enabled(self, tmp_path):
        bedpe_path = str(tmp_path / "out.bedpe")
        runner = _make_runner(fusion_enabled=True, output_bedpe=bedpe_path)
        aggregated = _make_aggregated()

        with patch("fin.pipeline.runner.generate_isolated_intervals",
                   return_value={"intervals": []}), \
             patch("fin.pipeline.runner.aggregate_across_intervals",
                   return_value=aggregated), \
             patch("fin.io.io_bedpe.write_fusion_bedpe") as mock_bedpe:

            runner.run()

        mock_bedpe.assert_called_once_with(aggregated, bedpe_path)

    def test_run_skips_bedpe_when_fusion_disabled(self, tmp_path):
        bedpe_path = str(tmp_path / "out.bedpe")
        runner = _make_runner(fusion_enabled=False, output_bedpe=bedpe_path)
        aggregated = _make_aggregated()

        with patch("fin.pipeline.runner.generate_isolated_intervals",
                   return_value={"intervals": []}), \
             patch("fin.pipeline.runner.aggregate_across_intervals",
                   return_value=aggregated), \
             patch("fin.io.io_bedpe.write_fusion_bedpe") as mock_bedpe:

            runner.run()

        mock_bedpe.assert_not_called()

    def test_run_skips_bedpe_when_path_none(self):
        runner = _make_runner(fusion_enabled=True, output_bedpe=None)
        aggregated = _make_aggregated()

        with patch("fin.pipeline.runner.generate_isolated_intervals",
                   return_value={"intervals": []}), \
             patch("fin.pipeline.runner.aggregate_across_intervals",
                   return_value=aggregated), \
             patch("fin.io.io_bedpe.write_fusion_bedpe") as mock_bedpe:

            runner.run()

        mock_bedpe.assert_not_called()


class TestGtfRegressionUnchanged:
    def test_run_existing_gtf_output_unchanged(self, tmp_path):
        gtf_path = str(tmp_path / "out.gtf")
        runner = _make_runner(output_gtf=gtf_path)
        aggregated = _make_aggregated()

        with patch("fin.pipeline.runner.generate_isolated_intervals",
                   return_value={"intervals": []}), \
             patch("fin.pipeline.runner.aggregate_across_intervals",
                   return_value=aggregated), \
             patch("fin.io.io_gtf.write_gtf") as mock_gtf:

            runner.run()

        mock_gtf.assert_called_once_with(aggregated, gtf_path)

    def test_run_gtf_not_called_when_path_none(self):
        runner = _make_runner(output_gtf=None)
        aggregated = _make_aggregated()

        with patch("fin.pipeline.runner.generate_isolated_intervals",
                   return_value={"intervals": []}), \
             patch("fin.pipeline.runner.aggregate_across_intervals",
                   return_value=aggregated), \
             patch("fin.io.io_gtf.write_gtf") as mock_gtf:

            runner.run()

        mock_gtf.assert_not_called()


# ---------------------------------------------------------------------------
# T8: max_R column in scoring TSV
# ---------------------------------------------------------------------------

class TestMaxRInScoringTsv:
    """Verify max_R column appears in scoring TSV."""

    def test_max_R_column_in_header(self, tmp_path):
        from fin.io.io_tsv import write_scoring_tsv
        qr = QuantResult(
            candidate_id="tx1",
            gene_id="tx1",
            chrom="chr1",
            strand="+",
            start=0,
            end=100,
            source="gtf",
            abundance=2.5,
            confidence=0.8,
            coherence_score=0.6,
            discrimination_score=0.7,
            combined_score=0.65,
            num_assigned_reads=3,
            max_R=0.75,
        )
        tsv_path = str(tmp_path / "scores.tsv")
        write_scoring_tsv({"tx1": qr}, {"tx1": 100}, tsv_path)

        with open(tsv_path) as f:
            lines = f.readlines()

        header = lines[0].strip().split("\t")
        assert "max_R" in header, f"max_R not in header: {header}"

    def test_max_R_value_written(self, tmp_path):
        from fin.io.io_tsv import write_scoring_tsv
        qr = QuantResult(
            candidate_id="tx1",
            gene_id="tx1",
            chrom="chr1",
            strand="+",
            start=0,
            end=100,
            source="gtf",
            abundance=2.5,
            confidence=0.8,
            coherence_score=0.6,
            discrimination_score=0.7,
            combined_score=0.65,
            num_assigned_reads=3,
            max_R=0.75,
        )
        tsv_path = str(tmp_path / "scores.tsv")
        write_scoring_tsv({"tx1": qr}, {"tx1": 100}, tsv_path)

        with open(tsv_path) as f:
            lines = f.readlines()

        header = lines[0].strip().split("\t")
        row = lines[1].strip().split("\t")
        max_r_idx = header.index("max_R")
        assert abs(float(row[max_r_idx]) - 0.75) < 1e-4

    def test_max_R_defaults_to_zero_when_unset(self, tmp_path):
        from fin.io.io_tsv import write_scoring_tsv
        # QuantResult without explicit max_R (defaults to 0.0)
        qr = QuantResult(
            candidate_id="tx1",
            gene_id="tx1",
            chrom="chr1",
            strand="+",
            start=0,
            end=100,
            source="gtf",
            abundance=1.0,
            confidence=1.0,
            coherence_score=0.5,
            discrimination_score=0.5,
            combined_score=0.5,
            num_assigned_reads=1,
        )
        tsv_path = str(tmp_path / "scores.tsv")
        write_scoring_tsv({"tx1": qr}, {"tx1": 100}, tsv_path)

        with open(tsv_path) as f:
            lines = f.readlines()

        header = lines[0].strip().split("\t")
        row = lines[1].strip().split("\t")
        max_r_idx = header.index("max_R")
        assert float(row[max_r_idx]) == 0.0


# ---------------------------------------------------------------------------
# _observed_junctions: recall-safe self-disable on missing/unreadable BAM
# (Lever 2 is default-on, so this runs on every m2_em interval).
# ---------------------------------------------------------------------------

import types  # noqa: E402


def _stub_interval():
    return types.SimpleNamespace(
        chrom="chr1", start=0, end=100, region_string="chr1:1-100"
    )


def test_observed_junctions_none_on_missing_bam():
    r = PipelineRunner.__new__(PipelineRunner)
    r.config = PipelineConfig(bam_path="/nonexistent/does_not_exist.bam")
    assert r._observed_junctions(_stub_interval()) is None


def test_observed_junctions_none_on_unreadable_bam():
    # An existing-but-invalid BAM (e.g. /dev/null) must self-disable, not raise.
    r = PipelineRunner.__new__(PipelineRunner)
    r.config = PipelineConfig(bam_path="/dev/null")
    assert r._observed_junctions(_stub_interval()) is None
