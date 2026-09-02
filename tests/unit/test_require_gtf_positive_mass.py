"""P0 strict-positive GTF evidence gate (require_gtf_positive_mass).

Contract under test:
1. Field/CLI default OFF everywhere (no profile enables it); default behavior
   is byte-identical to before the gate existed.
2. With the gate ON, a GTF-source model must end the A4 abundance filter with
   EM soft-mass > 0 even when every abundance floor is 0. Zero-mass catalog
   entries are unexpressed pass-through and are dropped, with a dedicated
   ``gtf_zero_mass`` outcome reason when the strict gate (not a floor) is the
   binding constraint.
3. Novel and fusion semantics are untouched; fusion stays exempt.
4. validate() fails loudly without a GTF: the gate scores GTF-source models,
   so requesting it in pure de novo mode is a caller error, never a no-op.
"""

from __future__ import annotations

import pytest

from fin.analysis.quantification import QuantResult
from fin.pipeline.config import PIPELINE_PROFILES, PipelineConfig
from fin.pipeline.runner import PipelineRunner
from fin.pipeline.selection import (
    enforce_gtf_positive_after_refit,
    select_global,
    strip_zero_mass_gtf_survivors,
)


def _qr(cid: str, source: str, abundance: float) -> QuantResult:
    return QuantResult(
        candidate_id=cid,
        abundance=abundance,
        confidence=0.0,
        num_assigned_reads=0,
        source=source,
    )


def _cfg(**overrides) -> PipelineConfig:
    base = dict(
        bam_path="/tmp/x.bam",
        min_abundance=0.0,
        min_gtf_abundance=0.0,
        floor_gtf_abundance=False,
        # Disable every other post-EM filter so A4 is isolated.
        min_isoform_fraction=0.0,
        min_fulllen_fraction=0.0,
        min_polya5p_reads=0,
    )
    base.update(overrides)
    return PipelineConfig(**base)


def _runner(**overrides) -> PipelineRunner:
    runner = PipelineRunner.__new__(PipelineRunner)
    runner.config = _cfg(**overrides)
    runner._gtf_reader = None
    return runner


class TestDefaults:
    def test_config_field_default_off(self):
        assert PipelineConfig(bam_path="/tmp/x.bam").require_gtf_positive_mass is False

    def test_cli_default_off(self):
        from fin.cli import main

        opt = next(
            p for p in main.params if p.name == "require_gtf_positive_mass"
        )
        assert opt.default is False

    def test_no_profile_enables_it(self):
        for name, values in PIPELINE_PROFILES.items():
            assert "require_gtf_positive_mass" not in values, name

    def test_default_off_keeps_zero_mass_gtf(self):
        """Backward compatibility: without the gate, floors=0 keeps GTF zeros."""
        runner = _runner()
        agg = {
            "g_zero": _qr("g_zero", "gtf", 0.0),
            "n_zero": _qr("n_zero", "novel", 0.0),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert set(out) == {"g_zero", "n_zero"}


class TestValidation:
    def test_requires_gtf(self, tmp_path):
        bam = tmp_path / "x.bam"
        bam.write_bytes(b"")
        cfg = PipelineConfig(
            bam_path=str(bam), require_gtf_positive_mass=True
        )
        with pytest.raises(ValueError, match="require_gtf_positive_mass requires --gtf"):
            cfg.validate()

    def test_passes_with_gtf(self, tmp_path):
        bam = tmp_path / "x.bam"
        gtf = tmp_path / "a.gtf"
        bam.write_bytes(b"")
        gtf.write_text("")
        cfg = PipelineConfig(
            bam_path=str(bam),
            gtf_path=str(gtf),
            require_gtf_positive_mass=True,
        )
        cfg.validate()
        assert cfg.require_gtf_positive_mass is True


class TestStrictPositiveGate:
    def test_zero_mass_gtf_dropped_when_all_floors_zero(self):
        """The gate must activate the A4 block on its own (floors all 0)."""
        runner = _runner(require_gtf_positive_mass=True)
        agg = {
            "g_zero": _qr("g_zero", "gtf", 0.0),
            "g_tiny": _qr("g_tiny", "gtf", 1e-9),
            "n_zero": _qr("n_zero", "novel", 0.0),
            "f_zero": _qr("f_zero", "fusion", 0.0),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert "g_zero" not in out         # zero-mass catalog entry dropped
        assert "g_tiny" in out             # smallest positive mass survives
        assert "n_zero" in out             # novel floor is 0 -> untouched
        assert "f_zero" in out             # fusion always exempt

    def test_negative_mass_gtf_dropped(self):
        runner = _runner(require_gtf_positive_mass=True)
        agg = {"g_neg": _qr("g_neg", "gtf", -1.0)}
        out = runner._finalize_and_write(agg, None, None)
        assert "g_neg" not in out

    def test_composes_with_gtf_floor(self):
        """Floor and strict gate apply together; floor stays binding above 0."""
        runner = _runner(
            require_gtf_positive_mass=True, min_gtf_abundance=1.0
        )
        agg = {
            "g_zero": _qr("g_zero", "gtf", 0.0),
            "g_half": _qr("g_half", "gtf", 0.5),
            "g_one": _qr("g_one", "gtf", 1.0),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert set(out) == {"g_one"}

    def test_novel_floor_semantics_unchanged(self):
        runner = _runner(require_gtf_positive_mass=True, min_abundance=3.0)
        agg = {
            "n_exact": _qr("n_exact", "novel", 3.0),
            "n_below": _qr("n_below", "novel", 2.999),
            "g_low": _qr("g_low", "gtf", 0.4),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert "n_exact" in out            # inclusive novel boundary kept
        assert "n_below" not in out
        assert "g_low" in out              # gtf floor 0, positive mass -> kept

    def test_outcome_reason_gtf_zero_mass(self):
        """When only the strict gate binds, the drop reason is gtf_zero_mass."""
        cfg = _cfg(require_gtf_positive_mass=True)
        agg = {
            "g_zero": _qr("g_zero", "gtf", 0.0),
            "g_pos": _qr("g_pos", "gtf", 0.7),
        }
        out, outcomes = select_global(cfg, agg, lambda a: a)
        assert set(out) == {"g_pos"}
        reasons = {o.candidate_id: o.reason for o in outcomes}
        assert reasons["g_zero"] == "gtf_zero_mass"

    def test_outcome_reason_floor_still_min_abundance(self):
        """When a floor also fails the model, keep the historical reason."""
        cfg = _cfg(require_gtf_positive_mass=True, min_gtf_abundance=1.0)
        agg = {
            "g_zero": _qr("g_zero", "gtf", 0.0),
            "g_one": _qr("g_one", "gtf", 1.0),
        }
        out, outcomes = select_global(cfg, agg, lambda a: a)
        assert set(out) == {"g_one"}
        reasons = {o.candidate_id: o.reason for o in outcomes}
        assert reasons["g_zero"] == "min_abundance"


class TestPostRefitSweep:
    """The refit can zero a GTF survivor that passed the pre-refit A4 gate;
    the post-refit sweep removes exactly those rows (mass/read-neutral)."""

    def test_disabled_returns_same_object(self):
        agg = {"g_zero": _qr("g_zero", "gtf", 0.0)}
        out, dropped = strip_zero_mass_gtf_survivors(agg, enabled=False)
        assert out is agg and dropped == []

    def test_removes_only_zero_mass_gtf(self):
        agg = {
            "g_zero": _qr("g_zero", "gtf", 0.0),
            "g_neg": _qr("g_neg", "gtf", -0.5),
            "g_pos": _qr("g_pos", "gtf", 1e-9),
            "n_zero": _qr("n_zero", "novel", 0.0),
            "f_zero": _qr("f_zero", "fusion", 0.0),
        }
        out, dropped = strip_zero_mass_gtf_survivors(agg, enabled=True)
        assert set(out) == {"g_pos", "n_zero", "f_zero"}
        assert dropped == ["g_neg", "g_zero"]  # sorted, deterministic

    def test_no_zero_rows_returns_same_object(self):
        agg = {"g_pos": _qr("g_pos", "gtf", 2.0)}
        out, dropped = strip_zero_mass_gtf_survivors(agg, enabled=True)
        assert out is agg and dropped == []


class TestDiagnosticsConsistency:
    """abundance_refit.json must describe the EMITTED survivor set: the sweep
    runs before the diagnostics write and records what it removed."""

    def test_annotates_diagnostics_with_drop_and_output_count(self):
        agg = {
            "g_zero": _qr("g_zero", "gtf", 0.0),
            "g_pos": _qr("g_pos", "gtf", 2.0),
            "n_pos": _qr("n_pos", "novel", 1.0),
        }
        diagnostics = {"surviving_candidates": 3}
        out, dropped = enforce_gtf_positive_after_refit(
            agg, diagnostics, enabled=True
        )
        assert dropped == ["g_zero"]
        assert diagnostics["post_refit_zero_gtf_dropped"] == 1
        assert diagnostics["post_refit_zero_gtf_ids"] == ["g_zero"]
        assert diagnostics["output_surviving_candidates"] == len(out) == 2
        assert (
            diagnostics["surviving_candidates"]
            - diagnostics["post_refit_zero_gtf_dropped"]
            == diagnostics["output_surviving_candidates"]
        )

    def test_disabled_adds_no_fields(self):
        """Default runs must stay byte-identical: no fields when disabled."""
        agg = {"g_zero": _qr("g_zero", "gtf", 0.0)}
        diagnostics = {"surviving_candidates": 1}
        out, dropped = enforce_gtf_positive_after_refit(
            agg, diagnostics, enabled=False
        )
        assert out is agg and dropped == []
        assert set(diagnostics) == {"surviving_candidates"}

    def test_none_diagnostics_tolerated(self):
        agg = {"g_zero": _qr("g_zero", "gtf", 0.0)}
        out, dropped = enforce_gtf_positive_after_refit(
            agg, None, enabled=True
        )
        assert set(out) == set() and dropped == ["g_zero"]

    def test_live_refit_zeroed_gtf_survivor_is_swept_consistently(self):
        """End-to-end over the real refit: a GTF survivor no ledger read
        references is refit to exactly zero, then swept, and the written
        diagnostics agree with the emitted set."""
        import numpy as np

        from fin.analysis.abundance_refit import (
            build_responsibility_ledger,
            refit_survivor_abundance,
        )

        ledger = build_responsibility_ledger(
            np.array([[1.0]]), ["r1"], ["n_pos"], input_read_ids=["r1"]
        )
        before = {
            "g_zero": _qr("g_zero", "gtf", 0.5),
            "n_pos": _qr("n_pos", "novel", 0.5),
        }
        out, diagnostics = refit_survivor_abundance(before, [ledger])
        assert out["g_zero"].abundance == 0.0  # the refit zeroed the GTF row
        out, dropped = enforce_gtf_positive_after_refit(
            out, diagnostics, enabled=True
        )
        assert dropped == ["g_zero"]
        assert set(out) == {"n_pos"}
        assert diagnostics["output_surviving_candidates"] == 1
        assert diagnostics["post_refit_zero_gtf_ids"] == ["g_zero"]
        # Mass conservation: the swept row carried zero mass.
        assert out["n_pos"].abundance == 1.0
