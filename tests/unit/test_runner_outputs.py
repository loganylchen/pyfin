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
