"""Post-selection abundance refit invariants."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from fin.analysis.abundance_refit import (
    ResponsibilityLedger,
    annotate_selection_metadata,
    build_responsibility_ledger,
    refit_survivor_abundance,
    write_abundance_refit_diagnostics,
)
from fin.analysis.quantification import QuantResult
from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner, _order_interval_outputs


def _result(cid, abundance=0.0, source="novel", exons=((0, 100),)):
    return QuantResult(
        candidate_id=cid,
        abundance=float(abundance),
        confidence=0.5,
        num_assigned_reads=1 if abundance else 0,
        source=source,
        chrom="chr1",
        strand="+",
        start=exons[0][0],
        end=exons[-1][1],
        exons=tuple(exons),
    )


def _candidate(cid, introns=()):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=400,
        sequence="A",
        source="novel",
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=0,
        end=400,
    )


def test_build_ledger_keeps_positive_finite_cells_and_input_population():
    ledger = build_responsibility_ledger(
        np.array([[0.75, 0.25], [1.0, 0.0]]),
        ["r1", "r2"],
        ["A", "B"],
        input_read_ids=["r0", "r1", "r2"],
    )
    assert ledger.weights == {
        "r1": {"A": 0.75, "B": 0.25},
        "r2": {"A": 1.0},
    }
    assert ledger.input_read_ids == ("r0", "r1", "r2")


def test_column_deletion_matches_beta_zero_softmax_rerun():
    distance = np.array([0.0, 1.0, 2.0])
    original = np.exp(-distance) / np.exp(-distance).sum()
    ledger = build_responsibility_ledger(
        original[None, :], ["read"], ["A", "B", "C"],
        input_read_ids=["read"],
    )
    before = {
        "A": _result("A", original[0], source="gtf"),
        "C": _result("C", original[2]),
    }
    out, diagnostics = refit_survivor_abundance(before, [ledger])

    rerun = np.exp(-distance[[0, 2]])
    rerun /= rerun.sum()
    assert out["A"].abundance == pytest.approx(rerun[0])
    assert out["C"].abundance == pytest.approx(rerun[1])
    assert out["A"].assigned_read_ids == ("read",)
    assert out["A"].exons == before["A"].exons
    assert diagnostics["mass_balance_error"] == pytest.approx(0.0)


def test_selection_orphan_is_separate_from_alignment_unassigned():
    ledger = ResponsibilityLedger(
        weights={"assigned": {"dropped": 1.0}},
        input_read_ids=("assigned", "never_entered_R"),
    )
    out, diagnostics = refit_survivor_abundance({}, [ledger])
    assert out == {}
    assert diagnostics["assignable_reads"] == 1
    assert diagnostics["alignment_unassigned_reads"] == 1
    # Scoped replacement key carries the same value as the legacy alias.
    assert diagnostics["interval_quantification_unassigned_reads"] == 1
    assert diagnostics["selection_orphaned_reads"] == 1
    assert diagnostics["unassigned_mass"] == pytest.approx(1.0)
    assert diagnostics["mass_balance_error"] == pytest.approx(0.0)


def test_forced_orphan_is_counted_in_both_overlapping_counters():
    """A forced read whose parent did not survive is forced AND orphaned."""
    ledger = ResponsibilityLedger(
        weights={"r1": {"gone": 1.0}, "r2": {"kept": 0.5, "gone": 0.5}},
        input_read_ids=("r1", "r2"),
        forced={"r1": "gone"},
    )
    results = {"kept": _result("kept", abundance=1.0)}
    _out, diagnostics = refit_survivor_abundance(results, [ledger])

    assert diagnostics["forced_reads"] == 1
    assert diagnostics["selection_orphaned_reads"] == 1
    assert diagnostics["forced_orphaned_reads"] == 1
    assert diagnostics["renormalized_reads"] == 1
    # The three counters overlap by exactly forced_orphaned_reads.
    assert (
        diagnostics["forced_reads"]
        + diagnostics["renormalized_reads"]
        + diagnostics["selection_orphaned_reads"]
        - diagnostics["forced_orphaned_reads"]
        == diagnostics["assignable_reads"]
    )
    assert diagnostics["mass_balance_error"] == pytest.approx(0.0)


def test_containment_redirect_moves_soft_column_before_renormalization():
    ledger = ResponsibilityLedger(
        weights={"r": {"parent": 0.2, "shadow": 0.3, "other": 0.5}},
        input_read_ids=("r",),
        redirects={"shadow": "parent"},
    )
    results = {
        "parent": _result("parent", 0.2, source="gtf"),
        "other": _result("other", 0.5),
    }
    out, _ = refit_survivor_abundance(results, [ledger])
    assert out["parent"].abundance == pytest.approx(0.5)
    assert out["other"].abundance == pytest.approx(0.5)
    assert out["parent"].assigned_read_ids == ("r",)  # gtf-first hard tie


def test_forced_mono_parent_wins_and_is_excluded_from_confidence_mean():
    forced_interval = ResponsibilityLedger(
        weights={"r": {"A": 0.9, "B": 0.1}},
        input_read_ids=("r",),
        interval_key="chr1:0-100:+",
        forced={"r": "B"},
    )
    uncovered_interval = ResponsibilityLedger(
        weights={"r": {"A": 1.0}},
        input_read_ids=("r",),
        interval_key="chr1:100-200:+",
    )
    out, diagnostics = refit_survivor_abundance(
        {"A": _result("A", 0.9), "B": _result("B", 0.1)},
        [uncovered_interval, forced_interval],
    )
    assert out["A"].abundance == 0.0
    assert out["B"].abundance == 1.0
    assert out["B"].assigned_read_ids == ("r",)
    assert out["B"].confidence == 0.0
    assert out["B"].max_R == 1.0
    assert diagnostics["forced_reads"] == 1


def test_forced_mapping_wins_across_intervals_and_conflicts_are_deterministic():
    ledgers = [
        ResponsibilityLedger(
            weights={"r": {"A": 0.9, "B": 0.1}},
            input_read_ids=("r",),
            forced={"r": "A"},
        ),
        ResponsibilityLedger(
            weights={"r": {"A": 0.2, "B": 0.8}},
            input_read_ids=("r",),
            forced={"r": "B"},
        ),
    ]
    results = {"A": _result("A", 2.0), "B": _result("B", 5.0)}
    out, diagnostics = refit_survivor_abundance(results, ledgers)
    assert out["B"].assigned_read_ids == ("r",)
    assert out["B"].abundance == 1.0
    assert diagnostics["forced_conflict_reads"] == 1


def test_duplicate_interval_observations_are_summed_then_normalized():
    ledgers = [
        ResponsibilityLedger(
            weights={"r": {"A": 0.8, "B": 0.2}}, input_read_ids=("r",)
        ),
        ResponsibilityLedger(
            weights={"r": {"A": 0.2, "B": 0.8}}, input_read_ids=("r",)
        ),
    ]
    results = {
        "A": _result("A", 1.0, source="gtf"),
        "B": _result("B", 1.0),
    }
    out, diagnostics = refit_survivor_abundance(results, ledgers)
    assert out["A"].abundance == pytest.approx(0.5)
    assert out["B"].abundance == pytest.approx(0.5)
    assert out["A"].assigned_read_ids == ("r",)
    assert diagnostics["refit_assigned_mass"] == pytest.approx(1.0)


def test_selection_metadata_distinguishes_containment_and_mono_semantics():
    parent = _candidate("parent", ((100, 200),))
    shadow = _candidate("shadow", ((100, 200),))
    mono = _candidate("mono")
    ledger = ResponsibilityLedger(
        weights={
            "shadow_read": {"shadow": 1.0},
            "mono_read": {"mono": 1.0},
        },
        input_read_ids=("mono_read", "shadow_read"),
    )
    surviving = [
        _result("parent", 2.0, source="gtf", exons=((0, 100), (200, 400))),
    ]
    surviving[0].assigned_read_ids = ("mono_read", "shadow_read")
    annotate_selection_metadata(
        ledger,
        candidates=[parent, shadow, mono],
        read_ids=["shadow_read", "mono_read"],
        hard_assignments=[1, 2],
        surviving_results=surviving,
        outcomes=[SimpleNamespace(
            action="fold", candidate_id="shadow", fold_into="parent"
        )],
        mono_resolve_applied=True,
    )
    assert ledger.redirects == {"shadow": "parent"}
    assert ledger.forced == {"mono_read": "parent"}


def test_refit_is_deterministic_under_parallel_ledger_order():
    ledgers = [
        ResponsibilityLedger(
            weights={"r": {"A": 0.1, "B": 0.9}},
            input_read_ids=("r",),
            interval_key="chr1:300-400:+",
        ),
        ResponsibilityLedger(
            weights={"r": {"A": 0.2, "B": 0.8}},
            input_read_ids=("r",),
            interval_key="chr1:100-200:+",
        ),
        ResponsibilityLedger(
            weights={"r": {"A": 0.7, "B": 0.3}},
            input_read_ids=("r",),
            interval_key="chr1:200-300:+",
        ),
    ]
    results = {
        "A": _result("A", 1.0, source="gtf"),
        "B": _result("B", 1.0),
    }
    forward, fd = refit_survivor_abundance(results, ledgers)
    reverse, rd = refit_survivor_abundance(results, list(reversed(ledgers)))
    assert forward == reverse
    assert fd == rd


def test_interval_output_order_is_canonical_only_when_refit_is_effective():
    late = ResponsibilityLedger(
        weights={"r2": {"A": 1.0}},
        input_read_ids=("r2",),
        interval_key="chr2:0-10:+",
    )
    early = ResponsibilityLedger(
        weights={"r1": {"A": 1.0}},
        input_read_ids=("r1",),
        interval_key="chr1:0-10:+",
    )
    outputs = [([_result("late")], late), ([_result("early")], early)]
    assert _order_interval_outputs(outputs, False) is outputs
    ordered = _order_interval_outputs(outputs, True)
    assert [ledger.interval_key for _, ledger in ordered] == [
        "chr1:0-10:+",
        "chr2:0-10:+",
    ]
    with pytest.raises(RuntimeError, match="requires a responsibility ledger"):
        _order_interval_outputs([([], None)], True)


def test_snap_redirect_is_applied_before_final_renormalization():
    ledger = ResponsibilityLedger(
        weights={"r": {"absorbed": 1.0}}, input_read_ids=("r",)
    )
    results = {"representative": _result("representative", 1.0)}
    out, diagnostics = refit_survivor_abundance(
        results, [ledger], snap_redirects={"absorbed": "representative"}
    )
    assert out["representative"].abundance == 1.0
    assert diagnostics["selection_orphaned_reads"] == 0


def test_runner_final_pass_writes_mass_diagnostics(tmp_path):
    bam = tmp_path / "input.bam"
    bam.touch()
    cfg = PipelineConfig(
        bam_path=str(bam),
        work_dir=str(tmp_path),
        post_selection_refit=True,
        enable_score_filter=False,
        junction_snap=False,
        use_gpu=False,
    )
    cfg.validate()
    runner = PipelineRunner(cfg)
    runner._gtf_reader = None
    runner._genome_fasta = None
    ledger = ResponsibilityLedger(
        weights={"r": {"kept": 0.25, "dropped": 0.75}},
        input_read_ids=("r",),
    )
    out = runner._finalize_and_write(
        {"kept": _result("kept", 0.25)},
        output_gtf=None,
        output_tsv=None,
        responsibility_ledgers=[ledger],
    )
    assert out["kept"].abundance == 1.0
    diagnostics = (tmp_path / "abundance_refit.json").read_text()
    assert '"mass_balance_error": 0.0' in diagnostics
    assert '"selection_orphaned_reads": 0' in diagnostics


def test_diagnostics_writer_is_atomic_json(tmp_path):
    path = tmp_path / "abundance_refit.json"
    write_abundance_refit_diagnostics(path, {"effective": True, "value": 1})
    assert path.read_text().endswith("\n")
    assert '"effective": true' in path.read_text()
    assert not path.with_suffix(".json.tmp").exists()
