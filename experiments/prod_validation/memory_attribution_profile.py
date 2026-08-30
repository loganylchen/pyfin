#!/usr/bin/env python3
"""Attribute pipeline peak memory before optimizing anything (P6 Milestone 6).

The ~65 GB whole-job peak RSS is a measurement, not an attribution. This
probe samples per-process RSS for a live pyfin run (parent + spawn workers)
and separately sizes the main known in-memory structures on the parent side
by monkeypatching the runner, so the dominant component is identified before
any CSR/COO or buffering rework is attempted.

Usage: run any fin.cli command through this module instead of fin.cli:

    PYFIN_MEMPROF_JSON=/path/out.json python3 memory_attribution_profile.py \
        --profile real-drna-precision --bam ... --threads 8 ...

Requires psutil only for the sampler; falls back to /proc parsing without it.
"""
from __future__ import annotations

import json
import os
import sys
import threading
import time
from pathlib import Path

OUT = Path(os.environ.get("PYFIN_MEMPROF_JSON", "memory_attribution.json"))
_SAMPLES: dict = {"parent_mb": [], "workers_mb": [], "n_workers_max": 0}
_SIZES: dict = {}


def _rss_mb(pid: int) -> float:
    try:
        with open(f"/proc/{pid}/statm") as handle:
            pages = int(handle.read().split()[1])
        return pages * os.sysconf("SC_PAGE_SIZE") / 1e6
    except Exception:
        return 0.0


def _children(pid: int) -> list[int]:
    kids = []
    try:
        for entry in os.listdir("/proc"):
            if not entry.isdigit():
                continue
            try:
                with open(f"/proc/{entry}/stat") as handle:
                    parts = handle.read().split()
                if int(parts[3]) == pid:
                    kids.append(int(entry))
            except Exception:
                continue
    except Exception:
        pass
    return kids


def _sampler(stop: threading.Event) -> None:
    me = os.getpid()
    while not stop.is_set():
        parent = _rss_mb(me)
        kids = _children(me)
        workers = sum(_rss_mb(k) for k in kids)
        _SAMPLES["parent_mb"].append(round(parent, 1))
        _SAMPLES["workers_mb"].append(round(workers, 1))
        _SAMPLES["n_workers_max"] = max(_SAMPLES["n_workers_max"], len(kids))
        stop.wait(2.0)


def _deep_mb(obj) -> float:
    """Rough recursive payload size for the ledger structures (MB)."""
    import sys as _sys

    seen = set()

    def walk(x):
        oid = id(x)
        if oid in seen:
            return 0
        seen.add(oid)
        size = _sys.getsizeof(x)
        if isinstance(x, dict):
            for k, v in x.items():
                size += walk(k) + walk(v)
        elif isinstance(x, (list, tuple, set, frozenset)):
            for v in x:
                size += walk(v)
        return size

    return round(walk(obj) / 1e6, 1)


def _patch_measurements() -> None:
    import fin.pipeline.runner as runner

    original_refit = runner.refit_survivor_abundance

    def measured_refit(results, ledgers, **kwargs):
        _SIZES["n_ledgers"] = len(ledgers)
        _SIZES["ledgers_payload_mb"] = _deep_mb(
            [l.weights for l in ledgers]
        )
        _SIZES["ledger_input_ids_mb"] = _deep_mb(
            [l.input_read_ids for l in ledgers]
        )
        _SIZES["n_survivors"] = len(results)
        return original_refit(results, ledgers, **kwargs)

    runner.refit_survivor_abundance = measured_refit

    original_finalize = runner.PipelineRunner._finalize_and_write

    def measured_finalize(self, aggregated, *args, **kwargs):
        genome = getattr(self, "_genome_fasta", None) or {}
        try:
            from fin.io.lazy_genome import LazyGenomeFasta
        except Exception:
            LazyGenomeFasta = ()  # type: ignore[assignment]
        if isinstance(genome, LazyGenomeFasta):
            # Never traverse values(): that would fetch every chromosome and
            # misreport the lazy mapping as an eager whole-genome load.
            _SIZES["genome_mode"] = "lazy"
            _SIZES["genome_cached_mb"] = round(
                sum(len(v) for v in genome._cache.values()) / 1e6, 1
            )
            _SIZES["genome_n_refs"] = len(genome)
        else:
            _SIZES["genome_mode"] = "eager"
            _SIZES["genome_fasta_mb"] = round(
                sum(len(v) for v in genome.values()) / 1e6, 1
            )
        _SIZES["n_aggregated_candidates"] = len(aggregated)
        _SIZES["aggregated_read_ids_mb"] = _deep_mb(
            [getattr(q, "assigned_read_ids", ()) for q in aggregated.values()]
        )
        return original_finalize(self, aggregated, *args, **kwargs)

    runner.PipelineRunner._finalize_and_write = measured_finalize


def main() -> None:
    _patch_measurements()
    stop = threading.Event()
    thread = threading.Thread(target=_sampler, args=(stop,), daemon=True)
    thread.start()
    started = time.monotonic()
    try:
        from fin.cli import main as cli_main

        cli_main()
    finally:
        stop.set()
        thread.join(timeout=5)
        parent = _SAMPLES["parent_mb"] or [0.0]
        workers = _SAMPLES["workers_mb"] or [0.0]
        report = {
            "schema_version": 1,
            "wall_seconds": round(time.monotonic() - started, 1),
            "parent_peak_mb": max(parent),
            "parent_final_mb": parent[-1],
            "workers_peak_sum_mb": max(workers),
            "n_workers_max": _SAMPLES["n_workers_max"],
            "whole_job_peak_mb_lower_bound": max(
                p + w for p, w in zip(parent, workers)
            ),
            "structure_sizes": _SIZES,
            "note": (
                "workers_peak_sum_mb is the sampled sum across spawn workers; "
                "structure_sizes are parent-side payload estimates at the "
                "finalize/refit boundary. Attribution guide: if workers "
                "dominate, the hotspot is per-worker state (genome copies, "
                "read caches); if the parent spike at refit dominates, the "
                "ledger conversion is the target."
            ),
        }
        OUT.parent.mkdir(parents=True, exist_ok=True)
        OUT.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"[memprof] wrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
