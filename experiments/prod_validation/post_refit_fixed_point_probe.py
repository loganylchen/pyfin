#!/usr/bin/env python3
"""Simulate `once` and `fixed_point` post-refit filtering without changing the run.

v9 production is `freeze_once`: gates run on pre-refit evidence and the survivor
set is frozen before requantification. This probe answers what would happen if
the assignment-dependent gates were reapplied *after* the refit.

Two differences from the simpler gate audit matter and are honoured here:

* Gates are applied **sequentially in production order**, each seeing the
  survivors of the previous gate, exactly like ``select_global``. Computing
  gates independently and taking a union is not equivalent.
* Every shrink is followed by a real refit, so released mass is actually
  re-dealt. That is the only way to learn whether mass is absorbed by other
  candidates or becomes orphaned - a forced read whose target is deleted can
  orphan immediately rather than flow to a sibling.

Full-length support is recomputed each round from the current
``assigned_read_ids``; cached values are stale after a refit. Candidates never
return once dropped, so the shrink-only loop terminates by construction.

The audited run itself is untouched: the original freeze-once refit output is
returned, so GTF/TSV stay byte-identical.

    PYFIN_FIXED_POINT_JSON=/path/out.json python3 post_refit_fixed_point_probe.py --profile ...
"""
from __future__ import annotations

import json
import os
import time
from dataclasses import replace
from pathlib import Path

import pysam

import fin.pipeline.runner as runner
from fin.analysis.quantification import (
    compute_fulllen_frac,
    fulllen_fraction_drops,
    isoform_fraction_drops,
    mono_exon_drops,
    soft_mass_ratio_drops,
)

_ORIGINAL_REFIT = runner.refit_survivor_abundance
_ORIGINAL_FINALIZE = runner.PipelineRunner._finalize_and_write
_OUTPUT = Path(
    os.environ.get("PYFIN_FIXED_POINT_JSON", "post_refit_fixed_point.json")
)
_MAX_ROUNDS = 25
_STATE = {}


def _finalize(self, *args, **kwargs):
    _STATE["config"] = self.config
    return _ORIGINAL_FINALIZE(self, *args, **kwargs)


def _read_ends(bam_path: str) -> dict[str, tuple[int, int]]:
    ends: dict[str, tuple[int, int]] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for record in bam.fetch(until_eof=True):
            if record.is_unmapped or record.is_secondary or record.is_supplementary:
                continue
            rid = record.query_name
            if rid is None or rid in ends:
                continue
            start, end = record.reference_start, record.reference_end
            if start is None or end is None:
                continue
            ends[rid] = (int(start), int(end))
    return ends


def _abundance_floor_drops(config, results) -> set[str]:
    gtf_floor = (
        max(config.min_gtf_abundance, config.min_abundance)
        if config.floor_gtf_abundance
        else config.min_gtf_abundance
    )
    if not (config.min_abundance > 0.0 or gtf_floor > 0.0):
        return set()
    strict = bool(getattr(config, "strict_novel_abundance_floor", False))
    drops = set()
    for cid, qr in results.items():
        if qr.source == "fusion":
            continue
        if qr.source == "gtf":
            ok = qr.abundance >= gtf_floor
        else:
            ok = (
                qr.abundance > config.min_abundance
                if strict
                else qr.abundance >= config.min_abundance
            )
        if not ok:
            drops.add(cid)
    return drops


def _sequential_gates(config, results, ends) -> tuple[dict[str, list[str]], set[str]]:
    """Apply assignment-dependent gates in production order, cascading."""
    live = dict(results)
    per_gate: dict[str, list[str]] = {}

    def apply(name, drops):
        drops = {cid for cid in drops if cid in live}
        per_gate[name] = sorted(drops)
        for cid in drops:
            live.pop(cid, None)

    apply("abundance_floor", _abundance_floor_drops(config, live))

    if config.min_isoform_fraction > 0.0:
        apply(
            "isoform_fraction",
            isoform_fraction_drops(
                live,
                config.min_isoform_fraction,
                locus=getattr(config, "isoform_fraction_locus", "family"),
            ),
        )
    else:
        per_gate["isoform_fraction"] = []

    if config.max_soft_mass_ratio > 0.0:
        apply("soft_mass_ratio", soft_mass_ratio_drops(live, config.max_soft_mass_ratio))
    else:
        per_gate["soft_mass_ratio"] = []

    if getattr(config, "drop_mono_exon_novel", False) and (
        config.min_mono_exon_reads > 0 or config.min_mono_exon_length > 0
    ):
        apply(
            "mono_exon",
            mono_exon_drops(
                live,
                min_reads=config.min_mono_exon_reads,
                min_len=config.min_mono_exon_length,
            ),
        )
    else:
        per_gate["mono_exon"] = []

    if config.min_fulllen_fraction > 0.0:
        window = getattr(config, "fulllen_window_bp", 25)
        min_reads = getattr(config, "fulllen_min_reads", 4)
        shadow = {cid: replace(qr) for cid, qr in live.items()}
        for qr in shadow.values():
            if qr.source == "novel" and len(qr.exons) >= 3:
                qr.fulllen_frac = compute_fulllen_frac(qr, ends, window, min_reads)
        apply("fulllen_fraction", fulllen_fraction_drops(shadow, config.min_fulllen_fraction))
    else:
        per_gate["fulllen_fraction"] = []

    # poly(A) needs the signal-derived map; inactive in all audited profiles.
    per_gate["polya5p"] = []
    return per_gate, set(results) - set(live)


def _simulate(config, base_results, ledgers, snap_redirects, mode, ends) -> dict:
    """Run `once` or `fixed_point` shrink-only filtering."""
    started = time.monotonic()
    current = dict(base_results)
    rounds = []
    dropped_total: set[str] = set()

    for index in range(_MAX_ROUNDS):
        refitted, diag = _ORIGINAL_REFIT(
            current, ledgers, snap_redirects=snap_redirects
        )
        per_gate, drops = _sequential_gates(config, refitted, ends)
        rounds.append(
            {
                "round": index + 1,
                "survivors_before": len(refitted),
                "drops": len(drops),
                "per_gate_drops": {k: len(v) for k, v in per_gate.items()},
                "dropped_ids": sorted(drops),
                "assigned_mass": round(float(diag["refit_assigned_mass"]), 4),
                "orphaned_reads": diag["selection_orphaned_reads"],
                "unassigned_mass": round(float(diag["unassigned_mass"]), 4),
                "mass_balance_error": diag["mass_balance_error"],
            }
        )
        if not drops:
            final = refitted
            final_diag = diag
            break
        assert not (drops & dropped_total), "resurrection or double-drop detected"
        dropped_total |= drops
        current = {cid: qr for cid, qr in refitted.items() if cid not in drops}
        if mode == "once":
            final, final_diag = _ORIGINAL_REFIT(
                current, ledgers, snap_redirects=snap_redirects
            )
            rounds.append(
                {
                    "round": "final_refit",
                    "survivors_before": len(final),
                    "drops": 0,
                    "per_gate_drops": {},
                    "dropped_ids": [],
                    "assigned_mass": round(float(final_diag["refit_assigned_mass"]), 4),
                    "orphaned_reads": final_diag["selection_orphaned_reads"],
                    "unassigned_mass": round(float(final_diag["unassigned_mass"]), 4),
                    "mass_balance_error": final_diag["mass_balance_error"],
                }
            )
            break
    else:
        raise RuntimeError(f"{mode} did not converge in {_MAX_ROUNDS} rounds")

    converged = mode == "fixed_point" and rounds[-1]["drops"] == 0
    return {
        "mode": mode,
        "converged": converged if mode == "fixed_point" else None,
        "shrink_rounds": sum(1 for r in rounds if r["drops"]),
        "total_dropped": len(dropped_total),
        "dropped_ids": sorted(dropped_total),
        "final_survivors": len(final),
        "final_assigned_mass": round(float(final_diag["refit_assigned_mass"]), 4),
        "final_orphaned_reads": final_diag["selection_orphaned_reads"],
        "final_unassigned_mass": round(float(final_diag["unassigned_mass"]), 4),
        "final_mass_balance_error": final_diag["mass_balance_error"],
        "rounds": rounds,
        "seconds": round(time.monotonic() - started, 3),
    }


def _probe(results, ledgers, *, snap_redirects=None):
    refitted, diagnostics = _ORIGINAL_REFIT(
        results, ledgers, snap_redirects=snap_redirects
    )
    config = _STATE.get("config")
    report: dict = {"schema_version": 1}
    try:
        if config is None:
            raise RuntimeError("pipeline config was not captured")
        ends = (
            _read_ends(config.bam_path)
            if config.min_fulllen_fraction > 0.0
            else {}
        )
        report["profile"] = getattr(config, "profile", None)
        report["freeze_once"] = {
            "mode": "freeze_once",
            "final_survivors": len(refitted),
            "final_assigned_mass": round(
                float(diagnostics["refit_assigned_mass"]), 4
            ),
            "final_orphaned_reads": diagnostics["selection_orphaned_reads"],
            "final_unassigned_mass": round(
                float(diagnostics["unassigned_mass"]), 4
            ),
            "surviving_ids_count": len(refitted),
        }
        for mode in ("once", "fixed_point"):
            report[mode] = _simulate(
                config, results, ledgers, snap_redirects or {}, mode, ends
            )
        report["fixed_point_vs_once_same_set"] = (
            report["once"]["dropped_ids"] == report["fixed_point"]["dropped_ids"]
        )
        report["orphan_mass_delta_fixed_point"] = round(
            report["fixed_point"]["final_unassigned_mass"]
            - report["freeze_once"]["final_unassigned_mass"],
            4,
        )
    except Exception as exc:  # never break the audited run
        report["error"] = f"{type(exc).__name__}: {exc}"
    _OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    _OUTPUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(f"[fixed-point-probe] wrote {_OUTPUT}", flush=True)
    return refitted, diagnostics


runner.refit_survivor_abundance = _probe
runner.PipelineRunner._finalize_and_write = _finalize


if __name__ == "__main__":
    from fin.cli import main

    main()
