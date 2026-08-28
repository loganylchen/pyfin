"""CLI profile resolution and manifest tests."""

from __future__ import annotations

import json

from click.testing import CliRunner

from fin.cli import main


class _FakePipelineRunner:
    configs = []

    def __init__(self, config):
        self.config = config
        self.__class__.configs.append(config)

    def setup(self):
        return None

    def run(self):
        return []

    def cleanup(self):
        return None


def _required_args(tmp_path):
    paths = {}
    for name in ("input.bam", "genome.fa", "reads.fq", "signal.blow5"):
        path = tmp_path / name
        path.touch()
        paths[name] = str(path)
    output = tmp_path / "out"
    return [
        "--bam", paths["input.bam"],
        "--genome", paths["genome.fa"],
        "--fastq", paths["reads.fq"],
        "--signal", paths["signal.blow5"],
        "--output-dir", str(output),
        "--no-gpu",
    ], output


def _invoke(monkeypatch, tmp_path, *extra):
    import fin.pipeline.runner as runner_module

    _FakePipelineRunner.configs.clear()
    monkeypatch.setattr(runner_module, "PipelineRunner", _FakePipelineRunner)
    args, output = _required_args(tmp_path)
    result = CliRunner().invoke(main, [*args, *extra])
    assert result.exit_code == 0, result.output
    assert _FakePipelineRunner.configs
    return _FakePipelineRunner.configs[-1], output


def test_default_cli_profile_is_real_drna(monkeypatch, tmp_path):
    cfg, output = _invoke(monkeypatch, tmp_path)
    assert cfg.profile == "real-drna"
    assert cfg.min_abundance == 1.0
    assert cfg.strict_novel_abundance_floor is True
    assert cfg.max_soft_mass_ratio == 0.0
    assert cfg.min_fulllen_fraction == 0.0
    assert cfg.min_polya5p_reads == 0
    assert cfg.isoform_fraction_locus == "family"
    assert cfg.post_selection_refit is True
    assert cfg.post_selection_refit_effective is True
    assert cfg.post_selection_refit_disable_reason is None
    assert cfg.drop_mono_exon_novel is True
    assert cfg.min_mono_exon_reads == 5
    assert cfg.junction_snap is True
    assert cfg.junction_snap_tolerance == 6
    assert cfg.junction_snap_min_support == 2
    assert cfg.junction_snap_min_ratio == 2.0
    assert cfg.m2_metric == "summed_llr"
    assert cfg.m2_metric_route == "fixed"
    assert cfg.m2_summed_llr_margin == 1.0
    assert cfg.m2_summed_llr_flank == 8
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["profile"] == "real-drna"
    assert manifest["config"]["m2_metric"] == "summed_llr"
    assert manifest["config"]["isoform_fraction_locus"] == "family"
    assert manifest["config"]["post_selection_refit_effective"] is True
    assert manifest["source_root"].endswith("/pyfin")
    assert len(manifest["source_sha256"]) == 64
    assert manifest["git_dirty"] in (True, False, None)


def test_real_precision_profile_adds_generation_support_floor(monkeypatch, tmp_path):
    cfg, _ = _invoke(
        monkeypatch, tmp_path, "--profile", "real-drna-precision"
    )
    assert cfg.profile == "real-drna-precision"
    assert cfg.drop_mono_exon_novel is True
    assert cfg.min_mono_exon_reads == 5
    assert cfg.junction_snap is True
    assert cfg.min_novel_reads == 2
    assert cfg.m2_metric == "summed_llr"


def test_isoform_fraction_overlap_mode_is_explicit(monkeypatch, tmp_path):
    cfg, output = _invoke(
        monkeypatch, tmp_path, "--isoform-fraction-locus", "overlap"
    )
    assert cfg.isoform_fraction_locus == "overlap"
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["config"]["isoform_fraction_locus"] == "overlap"


def test_named_profile_refit_can_be_explicitly_disabled(monkeypatch, tmp_path):
    cfg, output = _invoke(
        monkeypatch, tmp_path, "--no-post-selection-refit"
    )
    assert cfg.post_selection_refit is False
    assert cfg.post_selection_refit_effective is False
    assert "post_selection_refit" in cfg.profile_overrides
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["config"]["post_selection_refit"] is False


def test_real_profile_mono_gate_can_be_explicitly_disabled(monkeypatch, tmp_path):
    cfg, _ = _invoke(
        monkeypatch, tmp_path, "--no-drop-mono-exon-novel"
    )
    assert cfg.drop_mono_exon_novel is False
    assert cfg.min_mono_exon_reads == 5
    assert "drop_mono_exon_novel" in cfg.profile_overrides


def test_real_profile_junction_snap_can_be_explicitly_disabled(
    monkeypatch, tmp_path
):
    cfg, _ = _invoke(monkeypatch, tmp_path, "--no-junction-snap")
    assert cfg.junction_snap is False
    assert "junction_snap" in cfg.profile_overrides


def test_sirv_auto_metric_routes_without_guide(monkeypatch, tmp_path):
    cfg, output = _invoke(monkeypatch, tmp_path, "--profile", "sirv")
    assert cfg.m2_metric == "summed_llr"
    assert cfg.m2_metric_route == "auto-unguided"
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["config"]["m2_metric"] == "summed_llr"
    assert manifest["config"]["m2_metric_route"] == "auto-unguided"


def test_sirv_auto_metric_routes_stub_and_valid_gtf(monkeypatch, tmp_path):
    stub = tmp_path / "stub.gtf"
    stub.write_text(
        "chr1\tt\ttranscript\t1\t10\t.\t+\t.\ttranscript_id \"a\";\n"
    )
    cfg, _ = _invoke(
        monkeypatch, tmp_path, "--profile", "sirv", "--gtf", str(stub)
    )
    assert cfg.m2_metric == "summed_llr"
    assert cfg.m2_metric_route == "auto-unguided"

    guide = tmp_path / "guide.gtf"
    guide.write_text(
        stub.read_text()
        + "chr1\tt\ttranscript\t20\t30\t.\t+\t.\ttranscript_id \"b\";\n"
    )
    cfg, _ = _invoke(
        monkeypatch, tmp_path, "--profile", "sirv", "--gtf", str(guide)
    )
    assert cfg.m2_metric == "mean"
    assert cfg.m2_metric_route == "auto-guided"


def test_sirv_profile_and_explicit_cli_override(monkeypatch, tmp_path):
    cfg, output = _invoke(
        monkeypatch,
        tmp_path,
        "--profile", "sirv",
        "--min-abundance", "7",
        "--m2-metric", "summed_llr",
        "--m2-summed-llr-margin", "3.5",
    )
    assert cfg.profile == "sirv"
    assert cfg.min_abundance == 7.0
    assert cfg.strict_novel_abundance_floor is False
    assert cfg.floor_gtf_abundance is True
    assert cfg.min_fulllen_fraction == 0.1
    assert cfg.min_polya5p_reads == 0
    assert cfg.m2_metric == "summed_llr"
    assert cfg.m2_metric_route == "fixed"
    assert cfg.m2_summed_llr_margin == 3.5
    assert cfg.profile_overrides == (
        "m2_metric",
        "m2_summed_llr_margin",
        "min_abundance",
    )
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["profile_overrides"] == list(cfg.profile_overrides)


def test_boolean_no_flag_overrides_profile(monkeypatch, tmp_path):
    cfg, _ = _invoke(
        monkeypatch,
        tmp_path,
        "--profile", "sirv",
        "--no-floor-gtf-abundance",
    )
    assert cfg.floor_gtf_abundance is False
    assert cfg.profile_overrides == ("floor_gtf_abundance",)
