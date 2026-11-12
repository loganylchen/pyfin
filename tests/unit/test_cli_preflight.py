"""How the command line resolves the signal format and pore model, and when it refuses to start."""
from __future__ import annotations

import pytest
from click.testing import CliRunner

import fin.pipeline.runner as runner_module
from fin import preflight
from fin.cli import main
from fin.preflight import PreflightError, PreflightReport


class _Runner:
    configs = []

    def __init__(self, config):
        _Runner.configs.append(config)

    def setup(self):
        pass

    def run(self):
        return {}

    def cleanup(self):
        pass


@pytest.fixture
def invoke(tmp_path, monkeypatch):
    _Runner.configs.clear()
    monkeypatch.setattr(runner_module, "PipelineRunner", _Runner)

    def _invoke(*extra, signal="signal.blow5", with_fastq=False):
        args = []
        for flag, name in (("--bam", "in.bam"), ("--genome", "g.fa"), ("--signal", signal)):
            (tmp_path / name).touch()
            args += [flag, str(tmp_path / name)]
        if with_fastq:
            (tmp_path / "r.fq").touch()
            args += ["--fastq", str(tmp_path / "r.fq")]
        args += ["--output-dir", str(tmp_path / "out"), "--no-gpu", *extra]
        return CliRunner().invoke(main, args)

    return _invoke


def _config():
    assert len(_Runner.configs) == 1
    return _Runner.configs[0]


def test_the_reads_are_not_a_required_input(invoke):
    result = invoke("--no-preflight")
    assert result.exit_code == 0, result.output
    assert not _config().fastq_path


def test_signal_format_defaults_to_the_file_name(invoke):
    assert invoke("--no-preflight").exit_code == 0
    assert _config().signal_format == "slow5"


def test_an_unrecognisable_signal_name_is_a_usage_error(invoke):
    result = invoke("--no-preflight", signal="signal.dat")
    assert result.exit_code == 2 and "--signal-format" in result.output
    assert _Runner.configs == []
    assert invoke("--no-preflight", "--signal-format", "slow5", signal="signal.dat").exit_code == 0


@pytest.mark.parametrize("extra", [(), ("--no-preflight",), ("--signal-format", "pod5", "--no-preflight")])
def test_pod5_is_refused_whether_or_not_the_checks_run(invoke, extra):
    result = invoke(*extra, signal="reads.pod5")
    assert result.exit_code == 2
    assert "SLOW5/BLOW5 only" in result.output and "blue-crab p2s" in result.output
    assert _Runner.configs == []


@pytest.mark.parametrize("requested, expected", [("rna002", "rna002"), ("rna004", "rna004"), (None, "rna002")])
def test_without_the_checks_the_pore_is_what_was_asked_for(invoke, requested, expected):
    extra = ("--pore", requested) if requested else ()
    result = invoke("--no-preflight", *extra)
    assert result.exit_code == 0, result.output
    assert _config().krill_pore == expected


@pytest.mark.parametrize("requested", ["auto", "rna002"])
def test_the_checks_decide_the_pore(invoke, monkeypatch, requested):
    seen = {}

    def fake_preflight(**kwargs):
        seen.update(kwargs)
        return PreflightReport(signal_format="slow5", pore="rna004")

    monkeypatch.setattr(preflight, "run_preflight", fake_preflight)
    result = invoke("--pore", requested, with_fastq=True)
    assert result.exit_code == 0, result.output
    assert _config().krill_pore == "rna004"
    assert seen["pore"] == requested and seen["signal_format"] == "slow5"
    assert seen["m2_metric"] == _config().m2_metric != "auto"
    assert seen["fastq_path"].endswith("r.fq") and seen["fusion"] is False


def test_a_failed_check_stops_before_anything_is_written(invoke, monkeypatch, tmp_path):
    def failing(**_kwargs):
        raise PreflightError("the BAM has no index.\n  Fix:   samtools index in.bam")

    monkeypatch.setattr(preflight, "run_preflight", failing)
    result = invoke()
    assert result.exit_code == 2
    assert "Error: the BAM has no index." in result.output and "samtools index" in result.output
    assert _Runner.configs == []
    assert not (tmp_path / "out").exists()   # not even the output directory


def test_placeholder_inputs_are_caught_by_the_real_checks(invoke):
    result = invoke("--m2-metric", "off")
    assert result.exit_code == 2 and "cannot open the BAM" in result.output
    assert _Runner.configs == []


def test_the_summary_is_written_after_the_fusion_analysis_and_not_if_it_fails(invoke, monkeypatch, tmp_path):
    import json

    import fin.fusion.pipeline as fusion_pipeline

    order = []
    monkeypatch.setattr(fusion_pipeline, "run_fusion_analysis",
                        lambda *_a, **_k: order.append((tmp_path / "out" / "run_summary.json").exists()))
    result = invoke("--no-preflight", "--fusion")
    assert result.exit_code == 0, result.output
    assert order == [False]            # not yet written while the fusion analysis ran
    summary = json.loads((tmp_path / "out" / "run_summary.json").read_text())
    assert summary["fusion_analysis"] is True

    # the same directory again, and this time the run fails: the summary of the
    # earlier, successful run must not be left standing
    _Runner.configs.clear()

    def failing(*_a, **_k):
        raise RuntimeError("fusion analysis failed")

    monkeypatch.setattr(fusion_pipeline, "run_fusion_analysis", failing)
    result = invoke("--no-preflight", "--fusion")
    assert result.exit_code != 0
    assert (tmp_path / "out" / "run_manifest.json").exists()
    assert not (tmp_path / "out" / "run_summary.json").exists()


def test_a_failed_assembly_in_a_reused_directory_leaves_no_summary(invoke, monkeypatch, tmp_path):
    assert invoke("--no-preflight").exit_code == 0
    assert (tmp_path / "out" / "run_summary.json").exists()

    def failing_run(_self):
        raise RuntimeError("assembly failed")

    monkeypatch.setattr(_Runner, "run", failing_run)
    assert invoke("--no-preflight").exit_code != 0
    assert not (tmp_path / "out" / "run_summary.json").exists()


def test_a_run_without_fusion_says_so_in_its_summary(invoke, tmp_path):
    import json

    assert invoke("--no-preflight").exit_code == 0
    assert json.loads((tmp_path / "out" / "run_summary.json").read_text())["fusion_analysis"] is False


def test_an_unknown_pore_is_rejected_by_the_option(invoke):
    result = invoke("--no-preflight", "--pore", "rna003")
    assert result.exit_code == 2 and "rna003" in result.output
