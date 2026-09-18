#!/usr/bin/env python3
"""How much does leaving one interval's reads out change the endpoint fit?

`_endpoint_context` refits the sample's endpoint mixtures for EVERY interval
with that interval's reads excluded -- a leave-one-interval-out design so an
interval's own reads never inform the model that judges them. On a real
library that is ~13-16k refits per cell, each over ~127k observations, and it
is the dominant cost of a genome-wide run even after the 2.9x hoist.

This measures what the exclusion buys: for a sample of real intervals, fit
with and without exclusion and report the parameter deltas and the largest
change in the per-read log-likelihood ratio the assignment would see. It
changes nothing; it is evidence for a decision the owner has to make.
"""
import argparse
import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, "/SSD/logan/dev/pyfin")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, help="_gwstage/<sample>")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--out", required=True)
    ap.add_argument("--evidence", default="", help="sample_evidence json recorded by a finished cell")
    a = ap.parse_args()

    import pysam

    from fin.analysis.sample_nuisance import fit_from_evidence
    from fin.io.interval_manager import generate_isolated_intervals
    from fin.pipeline.config import PipelineConfig, PIPELINE_PROFILES
    from fin.pipeline.runner import PipelineRunner

    st = Path(a.stage)
    cfg = PipelineConfig(bam_path=str(st / "input.bam"), fastq_path=str(st / "reads.fq.gz"),
                         signal_path=str(st / "signal.blow5"), genome_fasta_path=a.genome,
                         gtf_path=a.gtf, work_dir=a.out, use_gpu=False)
    for k, v in PIPELINE_PROFILES["real-drna"].items():
        if hasattr(cfg, k):
            setattr(cfg, k, v)
    if a.evidence:
        cfg.sample_evidence_path = a.evidence   # bypass identity addressing
    runner = PipelineRunner(cfg)
    runner.setup()
    evidence = runner._nuisance_evidence()
    if evidence is None:
        raise SystemExit("no nuisance evidence (prepass cache cold?)")
    intervals = generate_isolated_intervals(cfg.bam_path, gtf_path=cfg.gtf_path,
                                            max_gap=cfg.max_gap, max_reads=cfg.max_reads,
                                            gtf_reader=runner._gtf_reader)["intervals"]
    generator = random.Random(7)
    sample = generator.sample(intervals, min(a.n, len(intervals)))

    base = fit_from_evidence(evidence)
    b_end, b_start = base.endpoint, base.start
    tes_ids = {rid for rid, _ in (getattr(evidence, "tes", ()) or ())}
    rows = []
    with pysam.AlignmentFile(cfg.bam_path) as bam:
        for iv in sample:
            reads = {r.query_name for r in bam.fetch(iv.chrom, iv.start, iv.end)
                     if not r.is_secondary and not r.is_supplementary}
            fit = fit_from_evidence(evidence, exclude_reads=reads)
            e, s = fit.endpoint, fit.start
            # largest per-read log-likelihood change over the distances the
            # assignment would actually evaluate: the excluded interval's own
            row = dict(interval=iv.region_string, n_reads=len(reads),
                       n_in_tes=len(reads & tes_ids),
                       d_b_narrow=None if e is None or b_end is None else e.b_narrow - b_end.b_narrow,
                       d_b_broad=None if e is None or b_end is None else e.b_broad - b_end.b_broad,
                       d_eta=None if e is None or b_end is None else e.eta - b_end.eta,
                       rel_b_narrow=None if e is None or b_end is None or b_end.b_narrow == 0
                       else abs(e.b_narrow - b_end.b_narrow) / b_end.b_narrow)
            if e is not None and b_end is not None:
                worst = 0.0
                for d in (0, 1, 3, 10, 30, 100, 300):
                    worst = max(worst, abs(e.loglik(d) - b_end.loglik(d)))
                row["max_abs_dloglik_over_grid"] = worst
            rows.append(row)
    Path(a.out).mkdir(parents=True, exist_ok=True)
    (Path(a.out) / "leave_out_effect.json").write_text(json.dumps(
        {"base": None if b_end is None else dict(b_narrow=b_end.b_narrow, b_broad=b_end.b_broad, eta=b_end.eta, n_fit=b_end.n_fit),
         "rows": rows}, indent=1))
    ok = [r for r in rows if r.get("rel_b_narrow") is not None]
    if ok:
        rel = sorted(r["rel_b_narrow"] for r in ok)
        ll = sorted(r["max_abs_dloglik_over_grid"] for r in ok)
        print(f"intervals measured: {len(ok)}   base fit: b_narrow={b_end.b_narrow:.4f} b_broad={b_end.b_broad:.3f} eta={b_end.eta:.5f} n={b_end.n_fit}")
        print(f"reads excluded per interval: median {sorted(r['n_reads'] for r in ok)[len(ok)//2]}, max {max(r['n_reads'] for r in ok)}; of which in TES evidence: max {max(r['n_in_tes'] for r in ok)}")
        print(f"relative change in b_narrow: median {rel[len(rel)//2]:.2e}  p95 {rel[int(len(rel)*0.95)]:.2e}  max {rel[-1]:.2e}")
        print(f"max |d loglik| over distance grid: median {ll[len(ll)//2]:.2e}  p95 {ll[int(len(ll)*0.95)]:.2e}  max {ll[-1]:.2e}")
    runner.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
