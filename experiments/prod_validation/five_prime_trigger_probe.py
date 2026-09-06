#!/usr/bin/env python3
"""What distinguishes a FALSE 5' contradiction from a true one?

On the uncorrupted human annotation every 5' repair is wrong by construction:
the annotation IS the truth there, so any read that makes the pipeline move a
model's TSS upstream is a false trigger. This measures the properties of the
triggering reads, to find the specificity gate the 5' test is missing (the 3'
test gets its specificity from requiring a polyA PASS call).

Per read: the upstream distance u to the nearest annotated TSS on its strand,
the 5' soft-clip length of its primary alignment (an alignment that starts
mid-read is not evidence about where the MOLECULE started), and how many other
reads start within +-tol of it (a genuine unannotated TSS has a cluster; a
run-through or degraded read is alone).
"""
import argparse, bisect, collections, json, re
from pathlib import Path

import pysam


def annotated_tss(gtf, chroms):
    out = collections.defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript" or f[0] not in chroms:
                continue
            s, e = int(f[3]) - 1, int(f[4])
            out[(f[0], f[6])].append(s if f[6] == "+" else e)
    return {k: sorted(set(v)) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--chroms", default="chr1,chr19,chr21")
    ap.add_argument("--tol", type=int, default=10)
    a = ap.parse_args()
    chroms = set(a.chroms.split(","))
    tss = annotated_tss(a.gtf, chroms)

    reads = []                      # (chrom, strand, start5, clip5, clip3)
    with pysam.AlignmentFile(a.bam, "rb") as bam:
        for r in bam.fetch():
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if r.reference_name not in chroms:
                continue
            strand = "-" if r.is_reverse else "+"
            ct = r.cigartuples or []
            lead = ct[0][1] if ct and ct[0][0] in (4, 5) else 0
            trail = ct[-1][1] if ct and ct[-1][0] in (4, 5) else 0
            # the read's 5' end in transcript orientation, and the clip there
            if strand == "+":
                five, clip5, clip3 = r.reference_start, lead, trail
            else:
                five, clip5, clip3 = r.reference_end, trail, lead
            reads.append((r.reference_name, strand, five, clip5, clip3))

    starts = collections.defaultdict(list)
    for c, s, five, _, _ in reads:
        starts[(c, s)].append(five)
    for k in starts:
        starts[k].sort()

    stats = collections.Counter()
    trig_clip, other_clip, trig_cluster, trig_u = [], [], [], []
    for c, s, five, clip5, clip3 in reads:
        pos = tss.get((c, s))
        if not pos:
            continue
        i = bisect.bisect_left(pos, five)
        near = [pos[j] for j in (i - 1, i) if 0 <= j < len(pos)]
        if not near:
            continue
        d = min((five - t for t in near), key=abs)
        u = -d if s == "+" else d          # >0 = upstream of the annotated TSS
        arr = starts[(c, s)]
        lo = bisect.bisect_left(arr, five - a.tol)
        hi = bisect.bisect_right(arr, five + a.tol)
        cluster = hi - lo
        if u > 50:                          # a trigger the kernel would not absorb
            stats["trigger"] += 1
            trig_clip.append(clip5); trig_cluster.append(cluster); trig_u.append(u)
        else:
            stats["consistent"] += 1
            other_clip.append(clip5)
    import statistics as st
    n = stats["trigger"] + stats["consistent"]
    print(f"reads with an annotated TSS nearby: {n}")
    print(f"  upstream >50bp (would trigger): {stats['trigger']} "
          f"({100*stats['trigger']/max(n,1):.1f}%)")
    for name, v in (("5' soft clip, TRIGGERING", trig_clip),
                    ("5' soft clip, consistent", other_clip)):
        if v:
            print(f"  {name}: median {st.median(v):.0f} "
                  f"mean {st.mean(v):.0f} >20bp {100*sum(1 for x in v if x>20)/len(v):.1f}% "
                  f">100bp {100*sum(1 for x in v if x>100)/len(v):.1f}%")
    if trig_cluster:
        print(f"  reads sharing the triggering start (+-{a.tol}bp): "
              f"median {st.median(trig_cluster):.0f}, "
              f"alone {100*sum(1 for x in trig_cluster if x == 1)/len(trig_cluster):.1f}%, "
              f">=3 {100*sum(1 for x in trig_cluster if x >= 3)/len(trig_cluster):.1f}%")
        print(f"  upstream distance: median {st.median(trig_u):.0f} "
              f"p90 {sorted(trig_u)[int(0.9*len(trig_u))]:.0f}")


if __name__ == "__main__":
    main()
