"""run_summary.json: what a run did, summed from per-interval shards."""
from __future__ import annotations

import json
import shutil

import pytest

from fin.pipeline import run_summary

ENGINE = dict(assigned_reads=9, tied_reads=4, scored_contrasts=2)


def test_shards_are_summed_and_removed(tmp_path):
    run_summary.record_interval(tmp_path, "chr1:0-100:+", reads=10, candidates=3, **ENGINE)
    run_summary.record_interval(tmp_path, "chr1:0-100:-", reads=5, candidates=1,
                                assigned_reads=5, tied_reads=0, scored_contrasts=0)
    totals = run_summary.merge_interval_shards(tmp_path)
    assert totals == {"reads": 15, "candidates": 4, "assigned_reads": 14, "tied_reads": 4,
                      "scored_contrasts": 2, "intervals_assigned": 2}
    assert not (tmp_path / run_summary.SHARD_DIR).exists()


def test_an_engine_that_records_no_signal_counters_reports_unknown_not_zero(tmp_path):
    run_summary.record_interval(tmp_path, "chr1:0-100:+", reads=10, candidates=3)
    totals = run_summary.merge_interval_shards(tmp_path)
    assert (totals["reads"], totals["candidates"], totals["intervals_assigned"]) == (10, 3, 1)
    assert totals["assigned_reads"] is totals["tied_reads"] is totals["scored_contrasts"] is None


def test_a_recorded_zero_is_a_zero(tmp_path):
    run_summary.record_interval(tmp_path, "chr1:0-100:+", reads=10, candidates=3,
                                assigned_reads=10, tied_reads=0, scored_contrasts=0)
    totals = run_summary.merge_interval_shards(tmp_path)
    assert (totals["tied_reads"], totals["scored_contrasts"]) == (0, 0)


def test_the_same_interval_recorded_twice_counts_once(tmp_path):
    # a re-assigned interval overwrites its own shard instead of adding to it
    run_summary.record_interval(tmp_path, "chr1:0-100:+", reads=10, candidates=3, **{**ENGINE, "tied_reads": 1})
    run_summary.record_interval(tmp_path, "chr1:0-100:+", reads=10, candidates=3, **{**ENGINE, "tied_reads": 7})
    totals = run_summary.merge_interval_shards(tmp_path)
    assert (totals["intervals_assigned"], totals["tied_reads"]) == (1, 7)


def test_no_shards_is_a_run_in_which_nothing_reached_an_engine(tmp_path):
    totals = run_summary.merge_interval_shards(tmp_path)
    assert totals["intervals_assigned"] == 0 and totals["reads"] == 0
    assert totals["tied_reads"] is None


def test_counter_names_are_checked(tmp_path):
    with pytest.raises(KeyError, match="refined"):
        run_summary.record_interval(tmp_path, "chr1:0-1:+", reads=1, candidates=1, refined=1)
    with pytest.raises(KeyError, match="missing \\['candidates'\\]"):
        run_summary.record_interval(tmp_path, "chr1:0-1:+", reads=1)


def test_stale_counters_are_cleared_or_the_run_stops(tmp_path, monkeypatch):
    run_summary.record_interval(tmp_path, "chr1:0-100:+", reads=10, candidates=3)
    run_summary.clear_interval_shards(tmp_path)
    assert not (tmp_path / run_summary.SHARD_DIR).exists()
    run_summary.clear_interval_shards(tmp_path)          # nothing to clear is fine

    run_summary.record_interval(tmp_path, "chr1:0-100:+", reads=10, candidates=3)
    monkeypatch.setattr(shutil, "rmtree", lambda *_a, **_k: None)   # removal that silently fails
    with pytest.raises(RuntimeError, match="could not clear the stale run-summary counters"):
        run_summary.clear_interval_shards(tmp_path)


def test_models_are_counted_by_origin():
    assert run_summary.model_counts(["gtf", "novel", "gtf", "fusion"]) == {
        "total": 4, "annotated": 2, "novel": 2}
    assert run_summary.model_counts([]) == {"total": 0, "annotated": 0, "novel": 0}


def _write(tmp_path, **over):
    kwargs = dict(intervals_total=3, models={"total": 2, "annotated": 1, "novel": 1},
                  elapsed_seconds=12.345, pore="rna004", signal_scoring="summed_llr",
                  quant_mode="m2_em", warnings=["w"])
    kwargs.update(over)
    return run_summary.write_run_summary(tmp_path, **kwargs)


def test_the_summary_file_says_whether_the_signal_was_scored(tmp_path):
    run_summary.record_interval(tmp_path, "chr1:0-100:+", reads=10, candidates=3, **ENGINE)
    written = _write(tmp_path)
    on_disk = json.loads((tmp_path / "run_summary.json").read_text())
    assert on_disk == written
    assert on_disk["intervals"] == {"total": 3, "assigned": 1}
    assert on_disk["reads"] == {"in_assigned_intervals": 10, "assigned": 9}
    assert on_disk["signal"] == {"scoring": "summed_llr", "pore": "rna004", "tied_reads": 4,
                                 "scored_contrasts": 2}
    assert on_disk["models"] == {"total": 2, "annotated": 1, "novel": 1}
    assert on_disk["quant_mode"] == "m2_em"
    assert on_disk["elapsed_seconds"] == 12.3 and on_disk["warnings"] == ["w"]
    assert on_disk["peak_rss_gib"] > 0 and on_disk["libraries"]["numpy"]
    assert not (tmp_path / run_summary.SHARD_DIR).exists()


def test_an_engine_without_signal_counters_is_summarised_as_unknown(tmp_path):
    run_summary.record_interval(tmp_path, "chr1:0-100:+", reads=10, candidates=3)
    on_disk = _write(tmp_path, quant_mode="argmax")
    assert on_disk["quant_mode"] == "argmax"
    assert on_disk["reads"] == {"in_assigned_intervals": 10, "assigned": None}
    assert on_disk["signal"]["tied_reads"] is None and on_disk["signal"]["scored_contrasts"] is None


def test_counters_handed_over_by_the_runner_are_used_as_they_are(tmp_path):
    run_summary.record_interval(tmp_path, "chr1:0-100:+", reads=10, candidates=3, **ENGINE)
    counters = run_summary.merge_interval_shards(tmp_path)      # what PipelineRunner.run() keeps
    on_disk = _write(tmp_path, counters=counters)
    assert on_disk["intervals"]["assigned"] == 1 and on_disk["signal"]["tied_reads"] == 4
