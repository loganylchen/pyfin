"""Tests for Step 3 novel combined_score filter in PipelineRunner.run()."""

from __future__ import annotations

from typing import Dict
from unittest.mock import patch

from fin.analysis.quantification import QuantResult
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner


def _qr(cid: str, source: str, combined: float) -> QuantResult:
    return QuantResult(
        candidate_id=cid,
        gene_id=cid,
        chrom="chr1",
        strand="+",
        start=0,
        end=100,
        source=source,
        abundance=1.0,
        confidence=1.0,
        coherence_score=0.5,
        discrimination_score=0.5,
        combined_score=combined,
        num_assigned_reads=1,
        exons=((0, 100),),
    )


def _runner(min_novel_combined_score: float) -> PipelineRunner:
    cfg = PipelineConfig(
        bam_path="/fake/reads.bam",
        min_novel_combined_score=min_novel_combined_score,
    )
    return PipelineRunner(cfg)


def _run_with(aggregated: Dict[str, QuantResult], runner: PipelineRunner):
    with patch("fin.pipeline.runner.generate_isolated_intervals",
               return_value={"intervals": []}), \
         patch("fin.pipeline.runner.aggregate_across_intervals",
               return_value=aggregated):
        return runner.run()


class TestNovelCombinedScoreFilter:
    def test_disabled_when_threshold_zero(self):
        runner = _runner(min_novel_combined_score=0.0)
        agg = {
            "n1": _qr("n1", "novel", 0.1),
            "g1": _qr("g1", "gtf", 0.1),
        }
        result = _run_with(agg, runner)
        assert "n1" in result
        assert "g1" in result

    def test_drops_novel_below_threshold(self):
        runner = _runner(min_novel_combined_score=0.3)
        agg = {
            "novel_low": _qr("novel_low", "novel", 0.25),
            "novel_hi": _qr("novel_hi", "novel", 0.40),
        }
        result = _run_with(agg, runner)
        assert "novel_low" not in result
        assert "novel_hi" in result

    def test_gtf_exempt_from_filter(self):
        runner = _runner(min_novel_combined_score=0.5)
        agg = {
            "gtf_low": _qr("gtf_low", "gtf", 0.1),
            "gtf_zero": _qr("gtf_zero", "gtf", 0.0),
        }
        result = _run_with(agg, runner)
        assert "gtf_low" in result
        assert "gtf_zero" in result

    def test_boundary_inclusive(self):
        runner = _runner(min_novel_combined_score=0.3)
        agg = {"n_eq": _qr("n_eq", "novel", 0.3)}
        result = _run_with(agg, runner)
        assert "n_eq" in result

    def test_mixed_set_only_drops_qualifying_novel(self):
        runner = _runner(min_novel_combined_score=0.288)
        agg = {
            "g1": _qr("g1", "gtf", 0.10),
            "g2": _qr("g2", "gtf", 0.50),
            "n1": _qr("n1", "novel", 0.20),  # drop
            "n2": _qr("n2", "novel", 0.35),  # keep
            "n3": _qr("n3", "novel", 0.288),  # keep (boundary)
        }
        result = _run_with(agg, runner)
        assert set(result.keys()) == {"g1", "g2", "n2", "n3"}
