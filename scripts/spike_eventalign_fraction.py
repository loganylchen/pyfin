#!/usr/bin/env python3
"""Phase-0 spike: measure what FRACTION of m2_em wall time is krill eventalign.

No production source edits -- monkeypatches krill.align_reads_variants /
align_read_variants (accumulate wall time) and PipelineRunner.process_interval
(per-interval wall), then runs the pipeline in-process via CliRunner with
--threads 1 (serial path, no spawn, so the patches hold). Prints the eventalign
work fraction, which gates the whole batching refactor (>= ~20% to be worth it).

Usage (inside the SIF):
  python -m scripts.spike_eventalign_fraction --bam ... --genome ... --fastq ...
      --signal ... --gtf ... --out /tmp/work
Backend is whatever krill the PYTHONPATH resolves (production = CPU krill).
"""
import argparse
import time

import krill
from click.testing import CliRunner

from fin.cli import main
from fin.pipeline.runner import PipelineRunner

_EV = {"t": 0.0, "n": 0}      # eventalign wall + call count
_IV = {"t": 0.0, "n": 0}      # process_interval wall + count


def _wrap(mod, name, acc):
    orig = getattr(mod, name)

    def timed(*a, **k):
        t0 = time.perf_counter()
        try:
            return orig(*a, **k)
        finally:
            acc["t"] += time.perf_counter() - t0
            acc["n"] += 1
    setattr(mod, name, timed)


def main_cli():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--signal", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--gpu", action="store_true")
    args = ap.parse_args()

    # eventalign timers (both batch + per-read paths)
    if hasattr(krill, "align_reads_variants"):
        _wrap(krill, "align_reads_variants", _EV)
    if hasattr(krill, "align_read_variants"):
        _wrap(krill, "align_read_variants", _EV)
    # per-interval total timer
    _orig_pi = PipelineRunner.process_interval

    def _timed_pi(self, interval):
        t0 = time.perf_counter()
        try:
            return _orig_pi(self, interval)
        finally:
            _IV["t"] += time.perf_counter() - t0
            _IV["n"] += 1
    PipelineRunner.process_interval = _timed_pi

    argv = [
        "--bam", args.bam, "--genome", args.genome, "--fastq", args.fastq,
        "--signal", args.signal, "--signal-format", "slow5",
        "--output-dir", args.out, "--quant-mode", "m2_em", "--threads", "1",
        "--gtf", args.gtf,
    ]
    argv += ["--gpu"] if args.gpu else ["--no-gpu"]

    wall0 = time.perf_counter()
    res = CliRunner().invoke(main, argv, catch_exceptions=False)
    wall = time.perf_counter() - wall0

    ev, iv = _EV["t"], _IV["t"]
    print("\n==== EVENTALIGN FRACTION SPIKE ====")
    print(f"exit_code           : {res.exit_code}")
    print(f"total run wall      : {wall:8.1f}s")
    print(f"process_interval sum: {iv:8.1f}s over {_IV['n']} intervals")
    print(f"eventalign sum      : {ev:8.1f}s over {_EV['n']} calls")
    if iv > 0:
        print(f"eventalign / interval-sum : {100*ev/iv:5.1f}%  (work fraction)")
    if wall > 0:
        print(f"eventalign / total-wall   : {100*ev/wall:5.1f}%")
    print("GATE #3: refactor worth pursuing only if eventalign >= ~20% of wall.")


if __name__ == "__main__":
    main_cli()
