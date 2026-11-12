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
        # empty placeholder files: there is nothing for the input checks to inspect
        "--no-preflight",
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
    # the endpoint channels and folded-TSS recovery are part of the default
    # operating point: the CLI must APPLY them, not just resolve them
    assert cfg.assignment_channels == "ends"
    assert cfg.folded_tss_recovery_mode == "require"
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


def test_cli_explicit_flags_override_the_promoted_defaults(monkeypatch, tmp_path):
    cfg, _ = _invoke(monkeypatch, tmp_path,
                     "--assignment-channels", "none",
                     "--folded-tss-recovery-mode", "off")
    assert cfg.assignment_channels == "none" and cfg.folded_tss_recovery_mode == "off"
    assert {"assignment_channels", "folded_tss_recovery_mode"} <= set(cfg.profile_overrides)


def test_sirv_profile_keeps_the_channels_off(monkeypatch, tmp_path):
    cfg, _ = _invoke(monkeypatch, tmp_path, "--profile", "sirv")
    assert cfg.assignment_channels == "none" and cfg.folded_tss_recovery_mode == "off"


def test_real_profile_enables_the_snap_chain_merge(monkeypatch, tmp_path):
    """The promoted value has to be APPLIED, not merely resolved: this repo has
    shipped a promoted profile setting that the CLI read and then dropped."""
    cfg, output = _invoke(monkeypatch, tmp_path)
    assert cfg.snap_merge_chain is True
    assert cfg.snap_merge_chain_window_bp == 25
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["config"]["snap_merge_chain"] is True


def test_real_profile_snap_chain_merge_can_be_explicitly_disabled(
    monkeypatch, tmp_path
):
    cfg, _ = _invoke(monkeypatch, tmp_path, "--no-snap-merge-chain")
    assert cfg.snap_merge_chain is False
    assert "snap_merge_chain" in cfg.profile_overrides


def test_the_snap_chain_merge_window_is_not_a_profile_field(monkeypatch, tmp_path):
    """The window is a knob, not part of the promoted operating point, and
    supplying it must not make the mode itself an explicit request."""
    cfg, _ = _invoke(monkeypatch, tmp_path, "--snap-merge-chain-window-bp", "40")
    assert cfg.snap_merge_chain_window_bp == 40
    assert cfg.snap_merge_chain is True
    assert "snap_merge_chain" not in cfg.profile_overrides


def test_the_precision_profile_inherits_the_snap_chain_merge(monkeypatch, tmp_path):
    """`real-drna-precision` is `real-drna` plus a higher novel-read threshold,
    so it inherits this too. The four-grid evidence was collected under
    `real-drna`; the inheritance is by construction, not by measurement."""
    cfg, _ = _invoke(monkeypatch, tmp_path, "--profile", "real-drna-precision")
    assert cfg.snap_merge_chain is True
    assert cfg.min_novel_reads == 2


# --- the support-gate attestation exemption, as a profile default ------------
#
# Promoted under docs/PRODUCTION_ATTESTATION_DEFAULT_DECISION.md. These go
# through the real Click command and assert on the CONFIG THE RUNNER IS HANDED,
# because the bug this lever actually hit was invisible to any test of the
# resolver: the profile entry existed, `resolve_profile_values` computed 3, and
# the CLI never lifted the key out of `resolved`, so the flag default of 0 was
# what ran. Nothing raised. The only symptom was a benchmark arm that came out
# byte-identical to the arm it was supposed to differ from.


def test_attestation_is_on_by_default(monkeypatch, tmp_path):
    cfg, _ = _invoke(monkeypatch, tmp_path)
    assert cfg.m2_support_gate_attest_reads == 1   # K=1 since 2026-09-15; K=3 landed first


def test_attestation_is_on_under_real_drna(monkeypatch, tmp_path):
    cfg, _ = _invoke(monkeypatch, tmp_path, "--profile", "real-drna")
    assert cfg.m2_support_gate_attest_reads == 1


def test_the_precision_derivative_inherits_it_by_design(monkeypatch, tmp_path):
    # real-drna-precision splats real-drna, so it inherits K=3. Intentional --
    # the exemption raised precision as well as recall -- but this exact
    # combination was not benchmarked, so it is pinned rather than assumed.
    cfg, _ = _invoke(monkeypatch, tmp_path, "--profile", "real-drna-precision")
    assert cfg.m2_support_gate_attest_reads == 1
    assert cfg.min_novel_reads == 2


def test_the_untested_sirv_profile_leaves_it_off(monkeypatch, tmp_path):
    cfg, _ = _invoke(monkeypatch, tmp_path, "--profile", "sirv")
    assert cfg.m2_support_gate_attest_reads == 0


def test_custom_leaves_it_off(monkeypatch, tmp_path):
    cfg, _ = _invoke(monkeypatch, tmp_path, "--profile", "custom")
    assert cfg.m2_support_gate_attest_reads == 0


def test_an_explicit_zero_disables_it_against_the_profile(monkeypatch, tmp_path):
    cfg, _ = _invoke(monkeypatch, tmp_path,
                     "--m2-support-gate-attest-reads", "0")
    assert cfg.m2_support_gate_attest_reads == 0


def test_an_explicit_value_overrides_the_profile(monkeypatch, tmp_path):
    cfg, _ = _invoke(monkeypatch, tmp_path,
                     "--m2-support-gate-attest-reads", "5")
    assert cfg.m2_support_gate_attest_reads == 5


def test_the_veto_it_was_validated_with_is_what_runs(monkeypatch, tmp_path):
    # K was frozen together with the veto; a different veto is a different
    # lever and would need its own decision run.
    cfg, _ = _invoke(monkeypatch, tmp_path)
    assert cfg.m2_support_gate_attest_veto == "novel"
    assert cfg.m2_support_gate_attest_veto_bp == 10


def test_the_manifest_records_the_operating_value(monkeypatch, tmp_path):
    # The manifest is how an arm's configuration is audited after the fact, so
    # this has to read the artifact on disk. Asserting on the in-memory config
    # instead would pass even if the manifest omitted the key entirely, which
    # is exactly the failure that makes a campaign unreconstructable.
    _, output = _invoke(monkeypatch, tmp_path)
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["config"]["m2_support_gate_attest_reads"] == 1


def test_no_profile_value_is_silently_dropped_on_the_way_to_the_runner(
        monkeypatch, tmp_path):
    """Every profile key must survive into the constructed config.

    The general form of the bug above. `m2_metric` is exempt because the CLI
    deliberately resolves the `auto` sentinel to a concrete metric before
    constructing the config (covered by the m2_metric tests); every other key
    must arrive unchanged.
    """
    from fin.pipeline.config import PIPELINE_PROFILES

    transformed = {"m2_metric"}
    for name, values in PIPELINE_PROFILES.items():
        if not values:
            continue
        cfg, _ = _invoke(monkeypatch, tmp_path, "--profile", name)
        for key, want in values.items():
            if key in transformed:
                continue
            assert getattr(cfg, key) == want, (
                f"profile {name!r} sets {key}={want!r} but the runner was "
                f"handed {getattr(cfg, key)!r} -- the value is computed and "
                f"then dropped somewhere between resolve_profile_values and "
                f"PipelineConfig"
            )
