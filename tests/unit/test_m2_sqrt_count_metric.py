"""Experimental sqrt-count M2 metric and contrast observability counters."""
import math
from collections import defaultdict

import pytest

from fin.pipeline.assignment import _finalize_contrast, _log_contrast_stats
from fin.pipeline.config import PipelineConfig


def test_summed_values_pass_through_and_stats_recorded():
    stats = defaultdict(float)
    nlls = {0: 10.0, 1: 14.0}
    out = _finalize_contrast("summed_llr", dict(nlls), {0: 5, 1: 5}, stats)
    assert out == nlls
    assert stats["decided"] == 1
    assert stats["same_event_count"] == 1
    assert stats["diff_event_count"] == 0
    assert stats["min_ev_sum"] == 5
    assert stats["margin_sum"] == pytest.approx(4.0)


def test_unequal_event_counts_are_counted():
    stats = defaultdict(float)
    _finalize_contrast("summed_llr", {0: 10.0, 1: 14.0}, {0: 3, 1: 9}, stats)
    assert stats["diff_event_count"] == 1
    assert stats["min_ev_sum"] == 3


def test_sqrt_count_rescales_by_sqrt_min_event_count():
    stats = defaultdict(float)
    out = _finalize_contrast(
        "sqrt_count_mean_llr", {0: 1.0, 1: 1.5}, {0: 9, 1: 16}, stats)
    scale = math.sqrt(9)
    assert out[0] == pytest.approx(1.0 * scale)
    assert out[1] == pytest.approx(1.5 * scale)
    # ordering preserved; margin recorded on the scaled values
    assert stats["margin_sum"] == pytest.approx(0.5 * scale)


def test_mean_metric_is_untouched():
    stats = defaultdict(float)
    nlls = {0: 1.0, 1: 2.0}
    assert _finalize_contrast("mean", dict(nlls), {}, stats) == nlls
    assert not stats


def test_log_contrast_stats_smoke(caplog):
    import logging
    stats = defaultdict(float)
    _finalize_contrast("sqrt_count_mean_llr", {0: 1.0, 1: 2.0}, {0: 4, 1: 4}, stats)
    with caplog.at_level(logging.INFO, logger="fin.pipeline.assignment"):
        _log_contrast_stats("sqrt_count_mean_llr", stats)
    assert "M2 sqrt_count_mean_llr contrasts" in caplog.text
    assert "decided=1" in caplog.text


def test_jsonl_records_schema(tmp_path):
    """Capture mode materializes per-attempt records with full context."""
    import json
    from types import SimpleNamespace
    from fin.pipeline.assignment import (
        _log_contrast_stats,
        _record_abstention,
    )

    stats = defaultdict(float)
    stats["_records"] = []
    cand = SimpleNamespace(candidate_id="nov_a")
    _record_abstention(stats, "read1", [0], [cand], "invalid_payload")
    _finalize_contrast(
        "summed_llr", {0: 10.0, 1: 14.0}, {0: 5, 1: 4}, stats,
        read_id="read2", cand_ids={0: "nov_a", 1: "nov_b"}, coverage=True,
    )
    config = SimpleNamespace(m2_contrast_stats_jsonl=True, work_dir=str(tmp_path))
    _log_contrast_stats("summed_llr", stats, config)
    files = list((tmp_path / "m2_contrasts").glob("*.jsonl"))
    assert len(files) == 1
    lines = [json.loads(l) for l in files[0].read_text().splitlines()]
    kinds = [l["kind"] for l in lines]
    assert kinds.count("abstention") == 1
    assert kinds.count("comparison") == 1
    assert kinds.count("aggregate") == 1
    comp = next(l for l in lines if l["kind"] == "comparison")
    assert comp["read_id"] == "read2"
    assert comp["winner_id"] == "nov_a" and comp["runner_id"] == "nov_b"
    assert comp["nll_sum_best"] == 10.0
    assert comp["nll_mean_best"] == 2.0   # 10/5
    assert comp["nll_mean_runner"] == 3.5  # 14/4
    assert comp["ev_best"] == 5 and comp["ev_runner"] == 4
    assert comp["coverage"] is True
    ab = next(l for l in lines if l["kind"] == "abstention")
    assert ab["reason"] == "invalid_payload" and ab["read_id"] == "read1"
    assert ab["tie_candidates"] == ["nov_a"]


def test_all_abstain_still_writes_jsonl(tmp_path):
    """A fully-abstaining locus must still appear in the calibration data."""
    import json
    from types import SimpleNamespace
    from fin.pipeline.assignment import _log_contrast_stats, _record_abstention

    stats = defaultdict(float)
    stats["_records"] = []
    cand = SimpleNamespace(candidate_id="c1")
    for reason in ("no_window", "unequal_intron_count", "invalid_payload",
                   "lt2_scored", "no_alignment", "no_backend",
                   "batch_output_missing"):
        _record_abstention(stats, f"r_{reason}", [0], [cand], reason)
    config = SimpleNamespace(m2_contrast_stats_jsonl=True, work_dir=str(tmp_path))
    _log_contrast_stats("summed_llr", stats, config)
    files = list((tmp_path / "m2_contrasts").glob("*.jsonl"))
    assert len(files) == 1
    lines = [json.loads(x) for x in files[0].read_text().splitlines()]
    reasons = {l["reason"] for l in lines if l["kind"] == "abstention"}
    assert reasons == {"no_window", "unequal_intron_count", "invalid_payload",
                       "lt2_scored", "no_alignment", "no_backend",
                       "batch_output_missing"}
    assert sum(1 for l in lines if l["kind"] == "aggregate") == 1


@pytest.mark.parametrize(
    "reason", ["no_backend", "no_alignment", "batch_output_missing"]
)
def test_backend_and_alignment_only_abstentions_still_write_jsonl(tmp_path, reason):
    """A locus that abstained ONLY for these reasons must still be recorded."""
    import json
    from types import SimpleNamespace
    from fin.pipeline.assignment import _log_contrast_stats, _record_abstention

    stats = defaultdict(float)
    stats["_records"] = []
    _record_abstention(stats, "r1", [0], [SimpleNamespace(candidate_id="c")], reason)
    config = SimpleNamespace(m2_contrast_stats_jsonl=True, work_dir=str(tmp_path))
    _log_contrast_stats("summed_llr", stats, config)
    files = list((tmp_path / "m2_contrasts").glob("*.jsonl"))
    assert len(files) == 1, f"{reason} produced no JSONL"
    lines = [json.loads(x) for x in files[0].read_text().splitlines()]
    assert any(l["kind"] == "abstention" and l["reason"] == reason for l in lines)
    assert any(l["kind"] == "aggregate" for l in lines)


def test_comparison_record_carries_explicit_deltas(tmp_path):
    import json
    from types import SimpleNamespace
    from fin.pipeline.assignment import _log_contrast_stats

    stats = defaultdict(float)
    stats["_records"] = []
    _finalize_contrast("summed_llr", {0: 10.0, 1: 14.0}, {0: 5, 1: 4}, stats,
                       read_id="r", cand_ids={0: "a", 1: "b"}, coverage=False)
    config = SimpleNamespace(m2_contrast_stats_jsonl=True, work_dir=str(tmp_path))
    _log_contrast_stats("summed_llr", stats, config)
    f = next((tmp_path / "m2_contrasts").glob("*.jsonl"))
    comp = next(json.loads(x) for x in f.read_text().splitlines()
                if json.loads(x)["kind"] == "comparison")
    required = {"read_id", "winner_col", "winner_id", "runner_id", "margin",
                "nll_mean_best", "nll_mean_runner", "nll_sum_best",
                "nll_sum_runner", "nll_mean_delta", "nll_sum_delta",
                "ev_best", "ev_runner", "n_scored", "coverage",
                "same_intron_count", "metric"}
    assert required <= set(comp)
    assert comp["nll_sum_delta"] == 4.0        # 14 - 10
    assert comp["nll_mean_delta"] == 1.5       # 14/4 - 10/5


def test_resolve_tie_applies_same_intron_count_guard():
    """Standalone/cluster path must abstain exactly like the batch path."""
    from types import SimpleNamespace
    from fin.scoring.m2_junction_nll import m2_resolve_tie

    def cand(introns):
        return SimpleNamespace(
            candidate_id=f"c{len(introns)}",
            intron_chain=SimpleNamespace(introns=tuple(introns)),
            sequence="ACGT" * 40, start=0, end=100, strand="+", chrom="chr1",
        )

    mixed = [cand([(10, 20)]), cand([(10, 20), (30, 40)])]
    for metric in ("summed_llr", "sqrt_count_mean_llr"):
        out = m2_resolve_tie(
            "r1", "ACGT" * 20, mixed, "/nonexistent.blow5", metric=metric,
        )
        assert out == (None, 0.0), f"{metric} must abstain on unequal chains"


def test_config_accepts_sqrt_count_metric(tmp_path):
    bam = tmp_path / "x.bam"
    bam.touch()
    cfg = PipelineConfig(bam_path=str(bam), m2_metric="sqrt_count_mean_llr")
    cfg.validate()
    assert cfg.m2_metric == "sqrt_count_mean_llr"


def test_config_rejects_unknown_metric(tmp_path):
    bam = tmp_path / "x.bam"
    bam.touch()
    cfg = PipelineConfig(bam_path=str(bam), m2_metric="bogus")
    with pytest.raises(ValueError):
        cfg.validate()
