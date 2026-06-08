"""Tests for the baked-in --min-abundance=3 abundance floor.

Covers:
1. The `fin` assembly CLI defaults --min-abundance to 3.0 and
   --floor-gtf-abundance to False (GTF exempt by default).
2. The runner A4 abundance filter drops low-abundance NOVEL transcripts but
   leaves GTF (and fusion) untouched by default.
3. With floor_gtf_abundance=True the floor additionally drops low-abundance GTF
   (reproduces the SIRV gffcompare T>=3 with-GTF operating point); fusion stays
   exempt regardless.

The dataclass field defaults stay min_abundance=0.0 / floor_gtf_abundance=False
(programmatic/quantify/ablation callers unchanged) — asserted too.
"""

from __future__ import annotations

from fin.analysis.quantification import QuantResult
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner


def _qr(cid: str, source: str, abundance: float) -> QuantResult:
    return QuantResult(
        candidate_id=cid,
        abundance=abundance,
        confidence=0.0,
        num_assigned_reads=0,
        source=source,
    )


class TestCliDefault:
    def test_cli_min_abundance_default_is_3(self):
        """The `fin` assembly command must default --min-abundance to 3.0."""
        from fin.cli import main

        opt = next(p for p in main.params if p.name == "min_abundance")
        assert opt.default == 3.0

    def test_cli_floor_gtf_abundance_default_off(self):
        """--floor-gtf-abundance must default OFF (GTF exempt)."""
        from fin.cli import main

        opt = next(p for p in main.params if p.name == "floor_gtf_abundance")
        assert opt.default is False

    def test_config_field_defaults(self):
        """Field defaults: min_abundance 0.0, floor_gtf_abundance False."""
        c = PipelineConfig(bam_path="/tmp/x.bam")
        assert c.min_abundance == 0.0
        assert c.floor_gtf_abundance is False


class TestRunnerAbundanceFloor:
    def _runner(self, min_ab: float, floor_gtf: bool = False) -> PipelineRunner:
        cfg = PipelineConfig(
            bam_path="/tmp/x.bam",
            min_abundance=min_ab,
            floor_gtf_abundance=floor_gtf,
            # Disable every other post-EM filter so we isolate A4 abundance.
            min_max_r=0.0,
            min_novel_combined_score=0.0,
            min_isoform_fraction=0.0,
            min_fulllen_fraction=0.0,
            min_polya5p_reads=0,
        )
        runner = PipelineRunner.__new__(PipelineRunner)
        runner.config = cfg
        runner._gtf_reader = None
        return runner

    def test_default_floors_novel_only_keeps_gtf_and_fusion(self):
        """Default (floor_gtf_abundance=False): novel dropped, GTF+fusion kept."""
        runner = self._runner(min_ab=0.5)
        agg = {
            "g_low": _qr("g_low", "gtf", 0.2),
            "n_low": _qr("n_low", "novel", 0.2),
            "n_high": _qr("n_high", "novel", 5.0),
            "f_low": _qr("f_low", "fusion", 0.0),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert "n_low" not in out          # novel below floor dropped
        assert "g_low" in out              # GTF exempt by default
        assert "n_high" in out             # novel above floor survives
        assert "f_low" in out              # fusion always exempt

    def test_floor_gtf_abundance_drops_low_gtf_keeps_fusion(self):
        """floor_gtf_abundance=True: low-abundance GTF also dropped, fusion kept."""
        runner = self._runner(min_ab=0.5, floor_gtf=True)
        agg = {
            "g_low": _qr("g_low", "gtf", 0.2),
            "g_high": _qr("g_high", "gtf", 5.0),
            "n_low": _qr("n_low", "novel", 0.2),
            "f_low": _qr("f_low", "fusion", 0.0),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert "g_low" not in out          # GTF now subject to floor
        assert "n_low" not in out
        assert "g_high" in out
        assert "f_low" in out              # fusion still exempt

    def test_floor_disabled_keeps_everything(self):
        runner = self._runner(min_ab=0.0)
        agg = {
            "g_low": _qr("g_low", "gtf", 0.0),
            "n_low": _qr("n_low", "novel", 0.0),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert set(out) == {"g_low", "n_low"}
