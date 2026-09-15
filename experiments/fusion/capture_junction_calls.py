#!/usr/bin/env python3
"""Capture every align_junction_window call of a real run, then measure.

Codex's instruction before optimising anything here: measure exact-key
duplication and the reach of the provably-exact shortcuts on the REAL call
distribution, not on a guess. The capture doubles as the differential corpus
any replacement must reproduce field for field.

Nothing in fin/ is modified: the wrapper is installed on the name
fin.fusion.junctions holds, because that module imported the function directly.
"""
import argparse
import json
import pickle
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, "/SSD/logan/dev/pyfin")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fastq", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--corpus", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--realign-workers", type=int, default=8)
    a = ap.parse_args()

    import fin.fusion.junctions as junctions
    from fin.fusion.alignment import _aligned_window, align_junction_window
    from fin.fusion.analysis import FusionSettings
    from fin.fusion.pipeline import run_fusion_analysis

    calls = []

    def capturing(query, target, junction, **kwargs):
        # Time every call COLD. `align_junction_window` memoises its recurrence,
        # so timing it as-is would record a dict lookup for exactly the repeated
        # windows whose avoided work this harness exists to measure -- the
        # saving would look small precisely where it is large. Clearing first
        # cannot change any result: the core is a pure function.
        _aligned_window.cache_clear()
        started = time.perf_counter()
        result = align_junction_window(query, target, junction, **kwargs)
        calls.append({
            "query": query, "target": target, "junction": junction,
            "kwargs": dict(sorted(kwargs.items())),
            "seconds": time.perf_counter() - started,
            "result": {f: getattr(result, f) for f in result.__dataclass_fields__},
        })
        return result

    junctions.align_junction_window = capturing

    settings = FusionSettings(min_clip_bases=250, max_candidate_reads=500000)
    run_fusion_analysis(a.bam, a.reference, a.out, gtf_path=a.gtf, sequence_path=a.fastq,
                        settings=settings, threads=a.threads, realign_workers=a.realign_workers)

    corpus = Path(a.corpus)
    corpus.parent.mkdir(parents=True, exist_ok=True)
    with corpus.open("wb") as handle:
        pickle.dump(calls, handle, protocol=4)

    total_seconds = sum(call["seconds"] for call in calls)
    keys = Counter((c["query"], c["target"], c["junction"],
                    tuple(c["kwargs"].items())) for c in calls)
    duplicate_calls = sum(count - 1 for count in keys.values())
    saved = 0.0
    seen = set()
    for call in calls:
        key = (call["query"], call["target"], call["junction"], tuple(call["kwargs"].items()))
        if key in seen:
            saved += call["seconds"]
        seen.add(key)

    substring = [c for c in calls if c["query"] and c["query"].upper().replace("U", "T")
                 in c["target"].upper().replace("U", "T")]
    substring_seconds = sum(c["seconds"] for c in substring)

    sizes = sorted(len(c["query"]) * len(c["target"]) for c in calls)
    report = {
        "calls": len(calls),
        "total_seconds": round(total_seconds, 2),
        "distinct_keys": len(keys),
        "duplicate_calls": duplicate_calls,
        "duplicate_fraction": round(duplicate_calls / len(calls), 4) if calls else 0,
        "seconds_a_perfect_cache_would_save": round(saved, 2),  # uncached cost of repeats
        "cache_saving_fraction": round(saved / total_seconds, 4) if total_seconds else 0,
        "query_is_exact_substring_of_target": len(substring),
        "substring_fraction": round(len(substring) / len(calls), 4) if calls else 0,
        "substring_seconds": round(substring_seconds, 2),
        "substring_saving_fraction": round(substring_seconds / total_seconds, 4) if total_seconds else 0,
        "cells_median": sizes[len(sizes) // 2] if sizes else 0,
        "cells_p95": sizes[int(len(sizes) * 0.95)] if sizes else 0,
        "cells_max": sizes[-1] if sizes else 0,
        "status_counts": dict(Counter(c["result"]["status"] for c in calls)),
        "reason_counts": dict(Counter(c["result"]["reason"] for c in calls)),
    }
    print(json.dumps(report, indent=2))
    corpus.with_suffix(".report.json").write_text(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
