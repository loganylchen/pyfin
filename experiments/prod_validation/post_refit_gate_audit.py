#!/usr/bin/env python3
"""Offline audit: which assignment-dependent gates does the refitted set violate?

v9 freezes existence selection on pre-refit evidence, so the emitted survivor
set is correct quantification under that set but is not necessarily a fixed
point of the gates that consume abundance or assigned reads. This probe
reapplies exactly those gates to the *refitted* results and reports what would
be dropped. It never changes the run: the original refit output is returned
untouched, so the pipeline's GTF/TSV stay byte-identical.

Structural/evidence gates (canonical, junction evidence, containment, M2) do not
consume final abundance and are deliberately out of scope.

Production ``fin/`` source is not modified. Run it exactly like ``fin.cli``:

    PYFIN_GATE_AUDIT_JSON=/path/out.json python3 post_refit_gate_audit.py --profile ...
"""
from __future__ import annotations

import json
import os
from collections import Counter
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
    os.environ.get("PYFIN_GATE_AUDIT_JSON", "post_refit_gate_audit.json")
)
_CONFIG = {}


def _finalize(self, *args, **kwargs):
    _CONFIG["config"] = self.config
    return _ORIGINAL_FINALIZE(self, *args, **kwargs)


def _read_ends(bam_path: str) -> dict[str, tuple[int, int]]:
    """Genome-wide primary-alignment span map (mirrors _annotate_fulllen_frac)."""
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


def _abundance_floor_violations(config, results) -> set[str]:
    """Mirror of the select_global absolute-floor keep condition."""
    gtf_floor = (
        max(config.min_gtf_abundance, config.min_abundance)
        if config.floor_gtf_abundance
        else config.min_gtf_abundance
    )
    if not (config.min_abundance > 0.0 or gtf_floor > 0.0):
        return set()
    strict = bool(getattr(config, "strict_novel_abundance_floor", False))
    violations = set()
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
            violations.add(cid)
    return violations


def _describe(results, cids) -> list[dict]:
    out = []
    for cid in sorted(cids):
        qr = results[cid]
        out.append(
            {
                "candidate_id": cid,
                "source": qr.source,
                "exons": len(qr.exons),
                "abundance": round(float(qr.abundance), 4),
                "num_assigned_reads": int(qr.num_assigned_reads),
                "chrom": qr.chrom,
                "strand": qr.strand,
            }
        )
    return out


def _audit(config, refitted) -> dict:
    gates: dict[str, dict] = {}

    def record(name, active, cids, note=None, evaluated=True):
        gates[name] = {
            "active_in_profile": active,
            "evaluated": bool(active and evaluated),
            "violations": len(cids) if (active and evaluated) else None,
            "rate_pct": (
                round(len(cids) / max(1, len(refitted)) * 100, 4)
                if (active and evaluated)
                else None
            ),
            "note": note,
            "examples": (
                _describe(refitted, list(cids)[:10])
                if (active and evaluated)
                else []
            ),
        }

    floor = _abundance_floor_violations(config, refitted)
    record("absolute_abundance_floor", True, floor)

    frac_on = config.min_isoform_fraction > 0.0
    frac = (
        isoform_fraction_drops(
            refitted,
            config.min_isoform_fraction,
            locus=getattr(config, "isoform_fraction_locus", "family"),
        )
        if frac_on
        else set()
    )
    record("isoform_fraction", frac_on, frac)

    ratio_on = config.max_soft_mass_ratio > 0.0
    ratio = (
        soft_mass_ratio_drops(refitted, config.max_soft_mass_ratio)
        if ratio_on
        else set()
    )
    record("soft_mass_ratio", ratio_on, ratio)

    mono_on = bool(
        getattr(config, "drop_mono_exon_novel", False)
        and (config.min_mono_exon_reads > 0 or config.min_mono_exon_length > 0)
    )
    mono = (
        mono_exon_drops(
            refitted,
            min_reads=config.min_mono_exon_reads,
            min_len=config.min_mono_exon_length,
        )
        if mono_on
        else set()
    )
    record("mono_exon", mono_on, mono)

    # Full-length support must be recomputed: refit changed assigned_read_ids,
    # so the cached fulllen_frac is stale evidence.
    full_on = config.min_fulllen_fraction > 0.0
    full: set[str] = set()
    full_note = None
    if full_on:
        # Copy before recomputing: the audit must never mutate the real results.
        shadow = {cid: replace(qr) for cid, qr in refitted.items()}
        ends = _read_ends(config.bam_path)
        window = getattr(config, "fulllen_window_bp", 25)
        min_reads = getattr(config, "fulllen_min_reads", 4)
        recomputed = 0
        for qr in shadow.values():
            if qr.source == "novel" and len(qr.exons) >= 3:
                qr.fulllen_frac = compute_fulllen_frac(
                    qr, ends, window, min_reads
                )
                recomputed += 1
        full = fulllen_fraction_drops(shadow, config.min_fulllen_fraction)
        full_note = (
            f"recomputed fulllen_frac from post-refit assigned reads for "
            f"{recomputed} novel multi-exon models"
        )
    record("fulllen_fraction", full_on, full, full_note)

    # poly(A) needs the signal-derived polya_map, which is not reconstructible
    # here. Report it as unevaluated rather than as zero violations.
    polya_on = getattr(config, "min_polya5p_reads", 0) > 0
    record(
        "polya5p",
        polya_on,
        set(),
        "active but not evaluated: requires the signal-derived polyA map"
        if polya_on
        else "inactive in this profile",
        evaluated=False,
    )

    union = floor | frac | ratio | mono | full
    overlap = Counter()
    for cid in union:
        hits = tuple(
            sorted(
                name
                for name, cids in (
                    ("absolute_abundance_floor", floor),
                    ("isoform_fraction", frac),
                    ("soft_mass_ratio", ratio),
                    ("mono_exon", mono),
                    ("fulllen_fraction", full),
                )
                if cid in cids
            )
        )
        overlap["+".join(hits)] += 1

    by_source = Counter(refitted[c].source for c in union)
    by_exons = Counter(
        "mono" if len(refitted[c].exons) < 2 else "multi" for c in union
    )
    return {
        "schema_version": 1,
        "profile": getattr(config, "profile", None),
        "surviving_candidates": len(refitted),
        "gates": gates,
        "union_violations": len(union),
        "union_rate_pct": round(len(union) / max(1, len(refitted)) * 100, 4),
        "union_by_source": dict(by_source),
        "union_by_exon_class": dict(by_exons),
        "gate_overlap": dict(overlap),
        "union_candidate_ids": sorted(union),
        "released_mass_if_dropped": round(
            sum(float(refitted[c].abundance) for c in union), 4
        ),
        "total_assigned_mass": round(
            sum(float(qr.abundance) for qr in refitted.values()), 4
        ),
    }


def _probe(results, ledgers, *, snap_redirects=None):
    refitted, diagnostics = _ORIGINAL_REFIT(
        results, ledgers, snap_redirects=snap_redirects
    )
    config = _CONFIG.get("config")
    try:
        if config is None:
            raise RuntimeError("pipeline config was not captured")
        stats = _audit(config, refitted)
    except Exception as exc:  # never break the audited run
        stats = {"schema_version": 1, "error": f"{type(exc).__name__}: {exc}"}
    _OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    _OUTPUT.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(f"[gate-audit] wrote {_OUTPUT}", flush=True)
    return refitted, diagnostics


runner.refit_survivor_abundance = _probe
runner.PipelineRunner._finalize_and_write = _finalize


if __name__ == "__main__":
    from fin.cli import main

    main()
