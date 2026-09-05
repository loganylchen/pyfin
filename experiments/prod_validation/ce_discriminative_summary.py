#!/usr/bin/env python3
"""Does the existence statistic separate TRUE from CORRUPTED annotation?

The grid's corruption conditions hand PyFIN an annotation whose intron chains
have been deliberately damaged for a known fraction of transcripts. Truth is
therefore available per candidate WITHOUT ever entering the pipeline: a supplied
annotated transcript is "corrupted" iff its chain in the input annotation differs
from its chain in the reference truth. That is the cleanest possible test of the
criterion -- it is exactly the discrimination the abundance gates cannot make.

Evaluation-only; nothing here is an input to PyFIN.
"""
import argparse
import collections
import csv
import glob
import re
import statistics as st


def chains(gtf):
    ex = collections.defaultdict(list)
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if m:
            ex[m.group(1)].append((f[0], int(f[3]), int(f[4])))
    out = {}
    for t, e in ex.items():
        e.sort(key=lambda x: x[1])
        out[t] = (e[0][0], tuple((a[2], b[1] - 1) for a, b in zip(e, e[1:])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True, help="run root with <sample>/<cond>/out/")
    ap.add_argument("--ref-root", required=True, help="benchmark _ref dir")
    ap.add_argument("--truth", required=True)
    a = ap.parse_args()
    truth = chains(a.truth)
    groups = collections.defaultdict(list)
    for f in sorted(glob.glob(f"{a.runs}/*/*/out/candidate_existence.tsv")):
        parts = f.split("/")
        sample, cond = parts[-4], parts[-3]
        try:
            anno = chains(f"{a.ref_root}/{sample}/{cond}/annotation.gtf")
        except OSError:
            continue
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r["identifiable"] != "true":
                continue
            cid = r["candidate_id"]
            sup, tru = anno.get(cid), truth.get(cid)
            if sup is None or tru is None:
                continue
            kind = "true_chain" if sup == tru else "corrupted_chain"
            groups[(cond, kind)].append(float(r["delta_predictive_ll"]))
    print(f"{'condition':20s} {'kind':17s} {'n':>4s} {'median':>9s} {'mean':>9s} "
          f"{'>0':>5s} {'<=0':>5s}")
    tot = collections.defaultdict(list)
    for (cond, kind), v in sorted(groups.items()):
        tot[kind].extend(v)
        print(f"{cond:20s} {kind:17s} {len(v):4d} {st.median(v):9.2f} {st.mean(v):9.2f} "
              f"{sum(1 for x in v if x > 1e-9):5d} {sum(1 for x in v if x <= 1e-9):5d}")
    print()
    for kind, v in sorted(tot.items()):
        if v:
            print(f"ALL {kind:17s} n={len(v):4d} median {st.median(v):8.2f} "
                  f"mean {st.mean(v):8.2f} positive {sum(1 for x in v if x > 1e-9)}/{len(v)}")
    if len(tot) == 2:
        t, c = tot["true_chain"], tot["corrupted_chain"]
        # rank-based separation (AUC): P(true > corrupted)
        wins = sum((1.0 if x > y else 0.5 if x == y else 0.0) for x in t for y in c)
        print(f"\nAUC(true > corrupted) = {wins / (len(t) * len(c)):.3f}"
              f"   (0.5 = no discrimination)")


if __name__ == "__main__":
    raise SystemExit(main())
