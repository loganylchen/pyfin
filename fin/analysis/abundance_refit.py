"""Post-selection abundance refit over the finalized transcript set.

The existence cascade deliberately runs on the original interval quantification.
This module is an additive final pass: it carries the beta=0 read/candidate
responsibilities through interval selection, redirects explicit folds, and
renormalizes each read over the candidates that survive global selection and
junction snapping.
"""
from __future__ import annotations

import json
import logging
import math
from collections import defaultdict
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple

import numpy as np

from fin.analysis.quantification import QuantResult, compute_tpm

logger = logging.getLogger(__name__)


@dataclass
class ResponsibilityLedger:
    """Picklable per-interval responsibility and selection metadata."""

    weights: Dict[str, Dict[str, float]]
    input_read_ids: Tuple[str, ...]
    interval_key: str = ""
    forced: Dict[str, str] = field(default_factory=dict)
    redirects: Dict[str, str] = field(default_factory=dict)


def build_responsibility_ledger(
    R: np.ndarray,
    read_ids: Sequence[str],
    candidate_ids: Sequence[str],
    *,
    input_read_ids: Iterable[str] = (),
    interval_key: str = "",
) -> ResponsibilityLedger:
    """Convert a dense final responsibility matrix into a sparse read map."""
    matrix = np.asarray(R, dtype=np.float64)
    if matrix.shape != (len(read_ids), len(candidate_ids)):
        raise ValueError(
            "responsibility axes do not match matrix shape: "
            f"R={matrix.shape}, reads={len(read_ids)}, candidates={len(candidate_ids)}"
        )

    weights: Dict[str, Dict[str, float]] = {}
    for i, rid in enumerate(read_ids):
        row: Dict[str, float] = {}
        for j in np.flatnonzero(np.isfinite(matrix[i]) & (matrix[i] > 0.0)):
            row[candidate_ids[int(j)]] = float(matrix[i, int(j)])
        if row:
            weights[str(rid)] = row

    return ResponsibilityLedger(
        weights=weights,
        input_read_ids=tuple(sorted({str(rid) for rid in input_read_ids})),
        interval_key=interval_key,
    )


def annotate_selection_metadata(
    ledger: ResponsibilityLedger,
    *,
    candidates: Sequence,
    read_ids: Sequence[str],
    hard_assignments: Sequence[int],
    surviving_results: Sequence[QuantResult],
    outcomes: Sequence,
    mono_resolve_applied: bool,
) -> ResponsibilityLedger:
    """Record candidate redirects and explicit mono-read assignments.

    Containment folds redirect a candidate's complete soft column. Mono
    resolution has different semantics: each resolved fragment read is pinned
    to its selected parent with mass one. The latter mapping is reconstructed
    from the pre-selection hard owner and the post-selection assigned-read sets,
    so the selection implementation itself remains unchanged.
    """
    ledger.redirects = {
        outcome.candidate_id: outcome.fold_into
        for outcome in outcomes
        if outcome.action == "fold" and outcome.fold_into
    }
    if not mono_resolve_applied:
        return ledger

    mono_ids = {
        candidate.candidate_id
        for candidate in candidates
        if not candidate.intron_chain.introns
    }
    original_owner: Dict[str, str] = {}
    candidate_ids = [candidate.candidate_id for candidate in candidates]
    for rid, col in zip(read_ids, hard_assignments):
        index = int(col)
        if 0 <= index < len(candidate_ids):
            original_owner[rid] = candidate_ids[index]

    post_owners: Dict[str, list[str]] = defaultdict(list)
    result_by_id = {result.candidate_id: result for result in surviving_results}
    for result in surviving_results:
        for rid in result.assigned_read_ids:
            post_owners[rid].append(result.candidate_id)

    for rid in sorted(original_owner):
        if original_owner[rid] not in mono_ids:
            continue
        owners = sorted(set(post_owners.get(rid, ())))
        if not owners:
            continue
        target = min(
            owners,
            key=lambda cid: (-result_by_id[cid].abundance, cid),
        )
        ledger.forced[rid] = target
    return ledger


def _resolve_redirect(candidate_id: str, *redirects: Mapping[str, str]) -> str:
    current = candidate_id
    seen = set()
    for mapping in redirects:
        while current in mapping and current not in seen:
            seen.add(current)
            current = mapping[current]
    return current


def _transcript_lengths(results: Mapping[str, QuantResult]) -> Dict[str, int]:
    return {
        cid: sum(end - start for start, end in result.exons)
        for cid, result in results.items()
    }


def refit_survivor_abundance(
    results: Dict[str, QuantResult],
    ledgers: Sequence[ResponsibilityLedger],
    *,
    snap_redirects: Mapping[str, str] | None = None,
) -> tuple[Dict[str, QuantResult], dict]:
    """Renormalize each assignable read over final surviving candidates.

    Interval observations of the same read/candidate pair are summed. For an
    ordinary read, deleting candidate columns and dividing each surviving weight
    by the surviving sum is numerically identical to rerunning a beta=0 softmax
    over those columns. A read with no surviving column contributes one unit of
    explicit selection-orphaned mass. Mono-resolved reads use their forced
    parent mapping instead of this renormalization.
    """
    snap_redirects = snap_redirects or {}
    survivors = set(results)
    merged_weights: Dict[str, Dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    forced_targets: Dict[str, list[str]] = defaultdict(list)
    input_reads = set()

    for ledger in sorted(ledgers, key=lambda item: item.interval_key):
        input_reads.update(ledger.input_read_ids)
        for rid in sorted(ledger.weights):
            for candidate_id, weight in sorted(ledger.weights[rid].items()):
                if not math.isfinite(weight) or weight <= 0.0:
                    continue
                target = _resolve_redirect(
                    candidate_id, ledger.redirects, snap_redirects
                )
                merged_weights[rid][target] += float(weight)
        for rid, candidate_id in sorted(ledger.forced.items()):
            target = _resolve_redirect(
                candidate_id, ledger.redirects, snap_redirects
            )
            forced_targets[rid].append(target)

    assignable_reads = set(merged_weights)
    abundance: Dict[str, float] = defaultdict(float)
    hard_ids: Dict[str, list[str]] = defaultdict(list)
    confidence_sum: Dict[str, float] = defaultdict(float)
    confidence_count: Dict[str, int] = defaultdict(int)
    max_responsibility: Dict[str, float] = defaultdict(float)

    forced_reads = 0
    forced_conflicts = 0
    renormalized_reads = 0
    orphaned_reads = 0
    pre_refit_survivor_mass = 0.0
    max_conservation_error = 0.0

    source_rank = {"gtf": 0, "novel": 1, "fusion": 2}

    for rid in sorted(assignable_reads):
        forced = sorted(set(forced_targets.get(rid, ())))
        probabilities: Dict[str, float] = {}
        is_forced = bool(forced)
        if forced:
            forced_reads += 1
            if len(forced) > 1:
                forced_conflicts += 1
            viable = [cid for cid in forced if cid in survivors]
            if viable:
                target = min(
                    viable,
                    key=lambda cid: (-results[cid].abundance, cid),
                )
                probabilities[target] = 1.0
                pre_refit_survivor_mass += 1.0
        else:
            weights = merged_weights[rid]
            total = sum(weights.values())
            surviving = {
                cid: weight for cid, weight in weights.items() if cid in survivors
            }
            surviving_total = sum(surviving.values())
            if total > 0.0:
                pre_refit_survivor_mass += surviving_total / total
            if surviving_total > 0.0:
                probabilities = {
                    cid: weight / surviving_total
                    for cid, weight in sorted(surviving.items())
                }
                renormalized_reads += 1

        if not probabilities:
            orphaned_reads += 1
            max_conservation_error = max(max_conservation_error, 0.0)
            continue

        assigned_mass = sum(probabilities.values())
        max_conservation_error = max(
            max_conservation_error, abs(assigned_mass - 1.0)
        )
        for candidate_id, probability in probabilities.items():
            abundance[candidate_id] += probability
            max_responsibility[candidate_id] = max(
                max_responsibility[candidate_id], probability
            )

        best_probability = max(probabilities.values())
        tied = [
            cid
            for cid, probability in probabilities.items()
            if abs(probability - best_probability) <= 1e-12
        ]
        winner = min(
            tied,
            key=lambda cid: (
                source_rank.get(results[cid].source, 99),
                cid,
            ),
        )
        hard_ids[winner].append(rid)
        if not is_forced:
            confidence_sum[winner] += probabilities[winner]
            confidence_count[winner] += 1

    if max_conservation_error > 1e-7:
        raise RuntimeError(
            "post-selection refit violated per-read mass conservation: "
            f"max_error={max_conservation_error:.3g}"
        )

    old_abundance = {cid: result.abundance for cid, result in results.items()}
    refitted: Dict[str, QuantResult] = {}
    for cid, result in results.items():
        assigned = tuple(sorted(hard_ids.get(cid, ())))
        count = len(assigned)
        confidence = (
            confidence_sum[cid] / confidence_count[cid]
            if confidence_count[cid] > 0
            else 0.0
        )
        refitted[cid] = replace(
            result,
            abundance=float(abundance.get(cid, 0.0)),
            confidence=float(confidence),
            num_assigned_reads=count,
            assigned_read_ids=assigned,
            max_R=float(max_responsibility.get(cid, 0.0)),
        )

    lengths = _transcript_lengths(results)
    old_tpm = compute_tpm(results, lengths)
    new_tpm = compute_tpm(refitted, lengths)
    abundance_shifts = {
        cid: refitted[cid].abundance - old_abundance[cid]
        for cid in results
    }
    tpm_shifts = {cid: new_tpm[cid] - old_tpm[cid] for cid in results}
    largest_abundance = max(
        abundance_shifts,
        key=lambda cid: (abs(abundance_shifts[cid]), cid),
        default=None,
    )
    top_tpm = sorted(
        results,
        key=lambda cid: (-abs(tpm_shifts[cid]), cid),
    )[:20]

    assigned_mass = sum(result.abundance for result in refitted.values())
    mass_target = float(len(assignable_reads))
    unassigned_mass = float(orphaned_reads)
    mass_balance_error = abs(mass_target - assigned_mass - unassigned_mass)
    if mass_balance_error > 1e-6 * max(1.0, mass_target):
        raise RuntimeError(
            "post-selection refit total mass did not close: "
            f"target={mass_target:.6f} assigned={assigned_mass:.6f} "
            f"unassigned={unassigned_mass:.6f}"
        )

    diagnostics = {
        "schema_version": 1,
        "mode": "post_selection_survivor_renormalization",
        "effective": True,
        "disable_reason": None,
        "selection_frozen_on_pre_refit_abundance": True,
        "confidence_excludes_forced_reads": True,
        "structural_identity": "not_checked_in_single_run",
        "input_reads": len(input_reads),
        "assignable_reads": len(assignable_reads),
        "alignment_unassigned_reads": len(input_reads - assignable_reads),
        "forced_reads": forced_reads,
        "forced_conflict_reads": forced_conflicts,
        "renormalized_reads": renormalized_reads,
        "selection_orphaned_reads": orphaned_reads,
        "unassigned_mass": unassigned_mass,
        "pre_refit_survivor_mass": pre_refit_survivor_mass,
        "rescued_mass": assigned_mass - pre_refit_survivor_mass,
        "old_survivor_abundance": sum(old_abundance.values()),
        "refit_assigned_mass": assigned_mass,
        "mass_balance_target": mass_target,
        "mass_balance_error": mass_balance_error,
        "max_per_read_conservation_error": max_conservation_error,
        "surviving_candidates": len(results),
        "candidates_with_abundance_shift": sum(
            abs(shift) > 1e-9 for shift in abundance_shifts.values()
        ),
        "abundance_l1_shift": sum(abs(x) for x in abundance_shifts.values()),
        "largest_abundance_shift": (
            {
                "candidate_id": largest_abundance,
                "before": old_abundance[largest_abundance],
                "after": refitted[largest_abundance].abundance,
                "delta": abundance_shifts[largest_abundance],
            }
            if largest_abundance is not None
            else None
        ),
        "tpm_l1_shift": sum(abs(x) for x in tpm_shifts.values()),
        "top_tpm_shifts": [
            {
                "candidate_id": cid,
                "before": old_tpm[cid],
                "after": new_tpm[cid],
                "delta": tpm_shifts[cid],
            }
            for cid in top_tpm
        ],
    }
    return refitted, diagnostics


def write_abundance_refit_diagnostics(path: str | Path, diagnostics: dict) -> None:
    """Atomically write the post-selection refit audit artifact."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_suffix(output.suffix + ".tmp")
    tmp.write_text(json.dumps(diagnostics, indent=2, sort_keys=True) + "\n")
    tmp.replace(output)
    logger.info("Wrote abundance-refit diagnostics: %s", output)
