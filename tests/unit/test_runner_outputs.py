"""Tests for PipelineRunner.run() output writers (US-013)."""

from __future__ import annotations

import json
from typing import Dict
from unittest.mock import MagicMock, patch

import numpy as np

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
# T8: R-matrix persistence tests
# ---------------------------------------------------------------------------

def _make_interval_mock(region="chr1:1000-2000"):
    """Return a minimal mock GenomicInterval."""
    m = MagicMock()
    m.region_string = region
    m.chrom = "chr1"
    return m


class TestRMatrixPersistence:
    """Unit tests for R.npy + R_meta.json writeout in process_interval."""

    def _run_process_interval_with_mocks(self, tmp_path, persist_R_matrix=True):
        """
        Call process_interval() on a PipelineRunner with everything below EM
        stubbed out. Returns (work_dir, quant_results).
        """
        from fin.pipeline.runner import PipelineRunner
        from fin.pipeline.config import PipelineConfig
        from fin.candidates.dataclasses import CandidateSet

        cfg = PipelineConfig(
            bam_path="/fake/reads.bam",
            work_dir=str(tmp_path),
            persist_R_matrix=persist_R_matrix,
            enable_signal=True,  # R.npy persistence is an EM-path feature
        )
        runner = PipelineRunner(cfg)

        # Synthetic interval
        interval = _make_interval_mock("chr1_1000_2000")

        # Synthetic inputs: 4 reads x 3 transcripts
        n_reads, n_tx = 4, 3
        read_ids = [f"read{i}" for i in range(n_reads)]
        cand_ids = [f"tx{j}" for j in range(n_tx)]

        # Build fake R matrix (rows sum to 1)
        R_fake = np.array([
            [0.7, 0.2, 0.1],
            [0.1, 0.8, 0.1],
            [0.2, 0.2, 0.6],
            [0.4, 0.4, 0.2],
        ], dtype=np.float32)
        hard_fake = np.argmax(R_fake, axis=1)

        # Minimal QuantResult list (what quantify_transcripts would return)
        from fin.analysis.quantification import QuantResult
        quant_fake = [
            QuantResult(
                candidate_id=cand_ids[j],
                abundance=float(R_fake[:, j].sum()),
                confidence=0.7,
                num_assigned_reads=1,
                source="gtf",
            )
            for j in range(n_tx)
        ]

        from fin.scoring.composite import CompositeScore
        composite_fake = [
            CompositeScore(candidate_id=cid, coherence=0.5, discrimination=0.5, combined=0.5)
            for cid in cand_ids
        ]

        runner_mod = "fin.pipeline.runner"
        with patch(f"{runner_mod}.discover_candidates") as mock_disc, \
             patch(f"{runner_mod}.ExternalToolRunner") as mock_tr_cls, \
             patch(f"{runner_mod}.parse_eventalign_tsv", return_value=[]), \
             patch(f"{runner_mod}.build_distance_matrix",
                   return_value=np.ones((n_reads, n_tx), dtype=np.float64)), \
             patch(f"{runner_mod}.subsample_reads_for_dtw",
                   return_value=read_ids), \
             patch(f"{runner_mod}.extract_signal_segments", return_value=[]), \
             patch(f"{runner_mod}.compute_read_to_read_dtw",
                   return_value=np.zeros((n_reads, n_reads), dtype=np.float64)), \
             patch(f"{runner_mod}.score_candidates_composite",
                   return_value=composite_fake), \
             patch(f"{runner_mod}.derive_prior_weights", return_value=None), \
             patch(f"{runner_mod}.em_with_coherence",
                   return_value=(R_fake, hard_fake, [])), \
             patch(f"{runner_mod}.quantify_transcripts",
                   return_value=quant_fake), \
             patch(f"{runner_mod}.populate_quant_scores"):

            # Build a fake CandidateSet
            mock_cs = MagicMock()
            mock_cs.num_candidates = n_tx
            mock_cs.read_ids = set(read_ids)
            mock_cs.candidate_ids.return_value = cand_ids
            # candidates list for n_tx transcripts
            mock_cs.candidates = [MagicMock(candidate_id=cid) for cid in cand_ids]
            mock_disc.return_value = mock_cs

            # fake tool runner
            mock_tr = MagicMock()
            mock_tr.score_candidates.return_value = []
            mock_tr_cls.return_value = mock_tr
            runner._tool_runner = mock_tr

            work_dir = tmp_path / interval.region_string.replace(":", "_").replace("-", "_")
            result = runner.process_interval(interval)

        return work_dir, result, R_fake, cand_ids, read_ids

    def test_R_npy_written(self, tmp_path):
        work_dir, result, R_fake, cand_ids, read_ids = \
            self._run_process_interval_with_mocks(tmp_path)
        r_path = work_dir / "R.npy"
        assert r_path.exists(), f"R.npy not found at {r_path}"

    def test_R_npy_shape_matches_meta(self, tmp_path):
        work_dir, result, R_fake, cand_ids, read_ids = \
            self._run_process_interval_with_mocks(tmp_path)
        R_loaded = np.load(str(work_dir / "R.npy"))
        with open(work_dir / "R_meta.json") as f:
            meta = json.load(f)
        assert R_loaded.shape == (len(meta["read_ids"]), len(meta["candidate_ids"]))
        assert R_loaded.shape == R_fake.shape

    def test_R_meta_json_fields(self, tmp_path):
        work_dir, result, R_fake, cand_ids, read_ids = \
            self._run_process_interval_with_mocks(tmp_path)
        with open(work_dir / "R_meta.json") as f:
            meta = json.load(f)
        assert "read_ids" in meta
        assert "candidate_ids" in meta
        assert "interval_region" in meta
        assert "sigma_used" in meta
        assert "was_subsampled" in meta
        assert "n_reads_subsampled" in meta
        assert "n_reads_full" in meta
        assert meta["candidate_ids"] == cand_ids
        assert meta["read_ids"] == read_ids

    def test_R_npy_not_written_when_disabled(self, tmp_path):
        work_dir, result, R_fake, cand_ids, read_ids = \
            self._run_process_interval_with_mocks(tmp_path, persist_R_matrix=False)
        assert not (work_dir / "R.npy").exists()
        assert not (work_dir / "R_meta.json").exists()

    def test_max_R_populated_on_quant_results(self, tmp_path):
        work_dir, result, R_fake, cand_ids, read_ids = \
            self._run_process_interval_with_mocks(tmp_path)
        assert result is not None
        for j, qr in enumerate(result):
            expected = float(R_fake[:, j].max())
            assert abs(qr.max_R - expected) < 1e-6, \
                f"max_R mismatch for col {j}: got {qr.max_R}, expected {expected}"


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
