#!/usr/bin/env python3
"""Run the independent fusion analysis on a real dRNA sample, standalone.

K562 carries BCR-ABL1 (chr22 -> chr9), so the output has a known positive
control before any new work is designed against it.
"""
import argparse, sys, time
from pathlib import Path

sys.path.insert(0, "/SSD/logan/dev/pyfin")
from fin.fusion.pipeline import run_fusion_analysis
from fin.fusion.analysis import FusionSettings


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--fastq", default=None)
    ap.add_argument("--signal", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=4, help="mappy index-build threads")
    ap.add_argument("--realign-workers", type=int, default=1,
                    help="forked workers sharing the built index; mappy holds the GIL in map() so "
                         "this is the only parallelism that helps the per-molecule genome search")
    ap.add_argument("--min-support", type=int, default=2)
    # dRNA reads carry ~80 bp of adapter/polyA soft clip on most molecules, so
    # the package default of 100 bp nominates ~8.6% of a real library (85k reads
    # on K562 rep4, which trips the 100k cap once readthrough nominations are
    # added). 250 bp is the threshold the assembly side settled on for the
    # same reason (the sc250 finding); it leaves 17k nominations here.
    ap.add_argument("--min-clip", type=int, default=250)
    ap.add_argument("--max-candidates", type=int, default=500000)
    ap.add_argument("--min-read-coverage", type=float, default=0.9,
                    help="fraction of the FULL read the stitched template must cover; real dRNA "
                         "reads carry ~80 bp of non-genomic 3' tail plus a few 5' bases, so 0.9 "
                         "rejects genuine fusions on shorter molecules")
    a = ap.parse_args()
    t0 = time.time()
    settings = FusionSettings(min_support=a.min_support, min_clip_bases=a.min_clip,
                              max_candidate_reads=a.max_candidates, min_read_coverage=a.min_read_coverage)
    analysis = run_fusion_analysis(
        a.bam, a.reference, a.out, gtf_path=a.gtf, sequence_path=a.fastq,
        settings=settings, threads=a.threads, signal_path=a.signal,
        realign_workers=a.realign_workers,
    )
    print(f"reads examined: {len(analysis.reads)}   chains: {len(analysis.chains)}   "
          f"records: {analysis.record_counts}   {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
