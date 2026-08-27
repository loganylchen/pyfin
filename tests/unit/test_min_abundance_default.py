"""Tests for the baked-in A4 abundance floors (novel + GTF).

Covers:
1. The `fin` assembly CLI defaults: --min-abundance 3.0, --min-gtf-abundance 1.0,
   --floor-gtf-abundance False.
2. The runner A4 abundance filter applies per-source floors on soft EM abundance:
   NOVEL faces min_abundance, GTF faces its own lighter min_gtf_abundance, fusion
   is always exempt.
3. With floor_gtf_abundance=True the GTF floor is raised to
   max(min_gtf_abundance, min_abundance) — never lowered — reproducing the SIRV
   gffcompare T>=3 with-GTF operating point; fusion stays exempt regardless.

The dataclass field defaults stay min_abundance=0.0 / min_gtf_abundance=0.0 /
floor_gtf_abundance=False / strict_novel_abundance_floor=False, so programmatic
callers see no floor unless they opt in. Named profiles may change both value and
boundary semantics.
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
        """--floor-gtf-abundance must default OFF (GTF uses its own floor)."""
        from fin.cli import main

        opt = next(p for p in main.params if p.name == "floor_gtf_abundance")
        assert opt.default is False

    def test_cli_strict_novel_floor_default_off(self):
        from fin.cli import main

        opt = next(
            p for p in main.params if p.name == "strict_novel_abundance_floor"
        )
        assert opt.default is False

    def test_cli_min_gtf_abundance_default_is_1(self):
        """The `fin` assembly command must default --min-gtf-abundance to 1.0."""
        from fin.cli import main

        opt = next(p for p in main.params if p.name == "min_gtf_abundance")
        assert opt.default == 1.0

    def test_config_field_defaults(self):
        """Field defaults: min_abundance 0.0, min_gtf_abundance 0.0, floor off."""
        c = PipelineConfig(bam_path="/tmp/x.bam")
        assert c.min_abundance == 0.0
        assert c.min_gtf_abundance == 0.0
        assert c.strict_novel_abundance_floor is False
        assert c.floor_gtf_abundance is False


class TestRunnerAbundanceFloor:
    def _runner(
        self,
        min_ab: float,
        floor_gtf: bool = False,
        min_gtf_ab: float = 0.0,
        strict_novel: bool = False,
    ) -> PipelineRunner:
        cfg = PipelineConfig(
            bam_path="/tmp/x.bam",
            min_abundance=min_ab,
            min_gtf_abundance=min_gtf_ab,
            strict_novel_abundance_floor=strict_novel,
            floor_gtf_abundance=floor_gtf,
            # Disable every other post-EM filter so we isolate A4 abundance.
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

    def test_strict_novel_floor_excludes_exact_boundary_only(self):
        runner = self._runner(min_ab=1.0, strict_novel=True)
        agg = {
            "n_below": _qr("n_below", "novel", 0.999),
            "n_exact": _qr("n_exact", "novel", 1.0),
            "n_above": _qr("n_above", "novel", 1.000001),
            "g_exact": _qr("g_exact", "gtf", 1.0),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert set(out) == {"n_above", "g_exact"}

    def test_inclusive_novel_floor_keeps_exact_boundary(self):
        runner = self._runner(min_ab=1.0, strict_novel=False)
        agg = {
            "n_below": _qr("n_below", "novel", 0.999),
            "n_exact": _qr("n_exact", "novel", 1.0),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert set(out) == {"n_exact"}

    def test_floor_disabled_keeps_everything(self):
        runner = self._runner(min_ab=0.0)
        agg = {
            "g_low": _qr("g_low", "gtf", 0.0),
            "n_low": _qr("n_low", "novel", 0.0),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert set(out) == {"g_low", "n_low"}

    def test_gtf_floor_drops_low_gtf_independent_of_novel_floor(self):
        """min_gtf_abundance gates GTF on soft abundance, even below the novel
        floor. A GTF with soft-mass < 1 is dropped though it has hard reads; a
        GTF at/above 1 is kept; fusion stays exempt."""
        runner = self._runner(min_ab=3.0, min_gtf_ab=1.0)
        agg = {
            "g_soft": _qr("g_soft", "gtf", 0.667),   # EM-diluted below 1 -> drop
            "g_ok": _qr("g_ok", "gtf", 2.0),         # >= 1 -> keep (below novel 3)
            "n_low": _qr("n_low", "novel", 2.0),     # < 3 novel floor -> drop
            "n_high": _qr("n_high", "novel", 5.0),   # >= 3 -> keep
            "f_low": _qr("f_low", "fusion", 0.0),    # fusion exempt
        }
        out = runner._finalize_and_write(agg, None, None)
        assert "g_soft" not in out
        assert "g_ok" in out
        assert "n_low" not in out
        assert "n_high" in out
        assert "f_low" in out

    def test_gtf_floor_runs_even_when_novel_floor_disabled(self):
        """With min_abundance=0 (novel floor off) the GTF floor still applies."""
        runner = self._runner(min_ab=0.0, min_gtf_ab=1.0)
        agg = {
            "g_soft": _qr("g_soft", "gtf", 0.5),     # < 1 -> drop
            "g_ok": _qr("g_ok", "gtf", 1.0),         # >= 1 -> keep
            "n_zero": _qr("n_zero", "novel", 0.0),   # novel floor 0 -> keep
        }
        out = runner._finalize_and_write(agg, None, None)
        assert "g_soft" not in out
        assert "g_ok" in out
        assert "n_zero" in out

    def test_floor_gtf_abundance_overrides_min_gtf_abundance(self):
        """--floor-gtf-abundance raises GTF to the novel floor, overriding the
        lighter min_gtf_abundance."""
        runner = self._runner(min_ab=3.0, floor_gtf=True, min_gtf_ab=1.0)
        agg = {
            "g_mid": _qr("g_mid", "gtf", 2.0),       # >= 1 but < 3 -> dropped now
            "g_high": _qr("g_high", "gtf", 5.0),     # >= 3 -> keep
            "f_low": _qr("f_low", "fusion", 0.0),    # fusion exempt
        }
        out = runner._finalize_and_write(agg, None, None)
        assert "g_mid" not in out
        assert "g_high" in out
        assert "f_low" in out
