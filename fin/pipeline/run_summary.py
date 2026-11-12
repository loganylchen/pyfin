"""What a finished run did, in one small file.

``run_manifest.json`` records what a run was asked to do. ``run_summary.json``
records what happened: how many reads reached assignment, how many of them were
in an exact alignment-score tie, how many of those ties the signal scored, what
came out, how long it took and how much memory it used. It answers "did the
signal contribute?" without reading a log.

Interval counters are written as shards by whichever process handled the
interval and summed at the end, like the other per-interval evidence. Nothing
here is read back by the pipeline.
"""
from __future__ import annotations

import hashlib
import json
import resource
import shutil
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional

SHARD_DIR = "_run_summary_shards"
# Known for every interval that reaches a quantification engine.
COUNTERS = ("reads", "candidates")
# Known only to the engine that breaks alignment-score ties with the signal
# (quant_mode m2_em, the engine of every shipped profile). Another engine does
# not record them, and the summary then says null, not zero.
ENGINE_COUNTERS = ("assigned_reads", "tied_reads", "scored_contrasts")


def clear_interval_shards(work_dir) -> None:
    """Remove counters left by an earlier run in this directory, or refuse to go on.

    A shard that survives would be summed into this run's record, so a failure
    to remove it is fatal rather than ignorable.
    """
    shard_dir = Path(work_dir) / SHARD_DIR
    shutil.rmtree(shard_dir, ignore_errors=True)
    if shard_dir.exists():
        raise RuntimeError(
            f"could not clear the stale run-summary counters at {shard_dir}; "
            "refusing to add this run's counters to another run's"
        )


def record_interval(work_dir, interval_key: str, **counters: int) -> None:
    """Write one interval's counters. Unknown counter names are a programming error."""
    unknown = set(counters) - set(COUNTERS) - set(ENGINE_COUNTERS)
    missing = set(COUNTERS) - set(counters)
    if unknown or missing:
        raise KeyError(f"run-summary counters: unknown {sorted(unknown)}, missing {sorted(missing)}")
    shard_dir = Path(work_dir) / SHARD_DIR
    shard_dir.mkdir(parents=True, exist_ok=True)
    name = hashlib.sha1(interval_key.encode()).hexdigest()[:20]
    payload = {"interval": interval_key, **{k: int(v) for k, v in counters.items()}}
    tmp = shard_dir / f"{name}.json.tmp"
    tmp.write_text(json.dumps(payload))
    tmp.replace(shard_dir / f"{name}.json")


def merge_interval_shards(work_dir) -> Dict[str, Optional[int]]:
    """Sum the shards and remove them.

    An engine counter that no interval recorded is ``None``: unknown, which is
    different from a count of zero.
    """
    shard_dir = Path(work_dir) / SHARD_DIR
    totals: Dict[str, Optional[int]] = {key: 0 for key in COUNTERS}
    totals.update({key: None for key in ENGINE_COUNTERS})
    totals["intervals_assigned"] = 0
    if shard_dir.is_dir():
        for shard in sorted(shard_dir.glob("*.json")):
            payload = json.loads(shard.read_text())
            totals["intervals_assigned"] += 1
            for key in COUNTERS:
                totals[key] += int(payload[key])
            for key in ENGINE_COUNTERS:
                if key in payload:
                    totals[key] = (totals[key] or 0) + int(payload[key])
        shutil.rmtree(shard_dir)
    return totals


def model_counts(sources: Iterable[str]) -> Dict[str, int]:
    """Output models by origin; ``sources`` is the ``source`` of each written model."""
    counts = {"total": 0, "annotated": 0, "novel": 0}
    for source in sources:
        counts["total"] += 1
        counts["annotated" if source == "gtf" else "novel"] += 1
    return counts


def peak_rss_gib() -> float:
    """Largest resident set of this process or of any finished worker, in GiB.

    It is the high-water mark of ONE process over its lifetime, not the memory
    of the run: workers run concurrently, and the kernel reports only the
    largest finished child, so the two are compared, not summed. A run with N
    workers can need up to N times this.
    """
    scale = 1 if sys.platform == "darwin" else 1024  # ru_maxrss: bytes on macOS, kB on Linux
    peak = max(
        resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,
    )
    return round(peak * scale / 1024 ** 3, 3)


def write_run_summary(
    work_dir,
    *,
    intervals_total: int,
    models: Dict[str, int],
    elapsed_seconds: float,
    pore: str,
    signal_scoring: str,
    quant_mode: str,
    fusion_analysis: bool = False,
    warnings: Optional[Iterable[str]] = None,
    counters: Optional[Dict[str, Optional[int]]] = None,
) -> Dict[str, object]:
    """Write ``run_summary.json``; return what was written.

    ``counters`` are the merged interval counters, as the runner keeps them
    after a run. Without them the shards are merged here.
    """
    from fin import __version__
    from fin.versions import library_versions

    if counters is None:
        counters = merge_interval_shards(work_dir)
    summary: Dict[str, object] = {
        "schema_version": 1,
        "pyfin_version": __version__,
        "libraries": library_versions(),
        "elapsed_seconds": round(float(elapsed_seconds), 1),
        "peak_rss_gib": peak_rss_gib(),
        "intervals": {"total": int(intervals_total), "assigned": counters["intervals_assigned"]},
        "reads": {"in_assigned_intervals": counters["reads"], "assigned": counters["assigned_reads"]},
        "quant_mode": quant_mode,
        # whether elapsed_seconds and peak_rss_gib include the optional fusion analysis
        "fusion_analysis": bool(fusion_analysis),
        "signal": {
            "scoring": signal_scoring,
            "pore": pore,
            # reads whose best alignment score is shared by two or more candidates
            "tied_reads": counters["tied_reads"],
            # ties in which at least two candidates received a valid signal score;
            # counted before the decision margin, so "scored", not "decided"
            "scored_contrasts": counters["scored_contrasts"],
        },
        "models": dict(models),
        "warnings": list(warnings or []),
    }
    path = Path(work_dir) / "run_summary.json"
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)
    return summary
