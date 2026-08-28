#!/usr/bin/env python3
"""Instrumentation-only probe for cross-interval responsibility ledgers.

Wraps ``fin.pipeline.runner.refit_survivor_abundance`` to measure how often a
read carries responsibilities from more than one interval, and whether those
per-interval rows disagree. Production ``fin/`` source is untouched: the probe
records statistics and then delegates to the original function.

Run it exactly like ``fin.cli``:

    PYFIN_LEDGER_PROBE_JSON=/path/out.json python3 ledger_overlap_probe.py --profile ...
"""
from __future__ import annotations

import json
import os
from collections import Counter, defaultdict
from pathlib import Path

import fin.pipeline.runner as runner
from fin.analysis.abundance_refit import _resolve_redirect

_ORIGINAL = runner.refit_survivor_abundance
_OUTPUT = Path(
    os.environ.get("PYFIN_LEDGER_PROBE_JSON", "ledger_overlap_probe.json")
)
_TOLERANCE = 1e-9


def _normalized(row: dict[str, float]) -> dict[str, float]:
    total = sum(row.values())
    if total <= 0.0:
        return {}
    return {cid: weight / total for cid, weight in row.items()}


def _rows_match(left: dict[str, float], right: dict[str, float]) -> bool:
    if set(left) != set(right):
        return False
    return all(abs(left[cid] - right[cid]) <= _TOLERANCE for cid in left)


def _probe(results, ledgers, *, snap_redirects=None):
    snap_redirects = snap_redirects or {}
    survivors = set(results)

    raw_sets: dict[str, list[frozenset]] = defaultdict(list)
    effective_sets: dict[str, list[frozenset]] = defaultdict(list)
    effective_rows: dict[str, list[dict[str, float]]] = defaultdict(list)
    forced_targets: dict[str, set[str]] = defaultdict(set)

    for ledger in ledgers:
        for rid, row in ledger.weights.items():
            raw_sets[rid].append(frozenset(row))
            mapped: dict[str, float] = defaultdict(float)
            for cid, weight in row.items():
                target = _resolve_redirect(cid, ledger.redirects, snap_redirects)
                mapped[target] += float(weight)
            effective_sets[rid].append(frozenset(mapped))
            effective_rows[rid].append(_normalized(dict(mapped)))
        for rid, cid in ledger.forced.items():
            forced_targets[rid].add(
                _resolve_redirect(cid, ledger.redirects, snap_redirects)
            )

    ledger_counts = Counter(len(v) for v in raw_sets.values())
    multi = [rid for rid, v in raw_sets.items() if len(v) > 1]

    raw_differ = [rid for rid in multi if len(set(raw_sets[rid])) > 1]
    effective_differ = [rid for rid in multi if len(set(effective_sets[rid])) > 1]

    identical_rows = 0
    differing_rows = 0
    for rid in multi:
        rows = effective_rows[rid]
        if all(_rows_match(rows[0], other) for other in rows[1:]):
            identical_rows += 1
        else:
            differing_rows += 1

    # Restrict to what actually changes an abundance: candidates that survived.
    surviving_differ = []
    for rid in multi:
        seen = {frozenset(s & survivors) for s in effective_sets[rid]}
        if len(seen) > 1:
            surviving_differ.append(rid)

    stats = {
        "schema_version": 1,
        "ledgers": len(ledgers),
        "surviving_candidates": len(survivors),
        "reads_with_responsibility": len(raw_sets),
        "reads_in_one_ledger": ledger_counts.get(1, 0),
        "reads_in_multiple_ledgers": len(multi),
        "ledger_count_distribution": {
            str(k): v for k, v in sorted(ledger_counts.items())
        },
        "multi_ledger_raw_candidate_sets_differ": len(raw_differ),
        "multi_ledger_effective_candidate_sets_differ": len(effective_differ),
        "multi_ledger_surviving_candidate_sets_differ": len(surviving_differ),
        "multi_ledger_identical_normalized_rows": identical_rows,
        "multi_ledger_differing_normalized_rows": differing_rows,
        "forced_reads": len(forced_targets),
        "forced_target_conflicts": sum(
            1 for targets in forced_targets.values() if len(targets) > 1
        ),
        "examples_multi_ledger": sorted(multi)[:5],
        "examples_effective_differ": sorted(effective_differ)[:5],
    }
    _OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    _OUTPUT.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    print(f"[ledger-probe] wrote {_OUTPUT}", flush=True)
    return _ORIGINAL(results, ledgers, snap_redirects=snap_redirects)


runner.refit_survivor_abundance = _probe


if __name__ == "__main__":
    from fin.cli import main

    main()
