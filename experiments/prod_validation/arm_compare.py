#!/usr/bin/env python3
"""Paired comparison of two pyfin labels in an expression-contract scores.tsv.

Reports per-cell Sn/Pr/F1 deltas, the mean delta with a bootstrap CI, and
per-family means. Cells present for only one label, or with a non-ok status,
are LISTED as missing rather than silently dropped -- a promotion decision
needs the complete manifest. The bootstrap resamples LIBRARIES (samples) as
blocks when the grid has more than one, keeping each library's conditions
together; with a single library it resamples cells and says so, because one
library's conditions are not independent replicates.
"""
import argparse, csv, collections, statistics, random

FAMILY = {
    "full": "full",
    "p00": "ratio", "p01": "ratio", "p05": "ratio", "p10": "ratio", "p20": "ratio",
    "p50": "ratio", "p80": "ratio", "p90": "ratio", "p99": "ratio",
    "c_jitter10_10bp": "c_jitter", "c_jitter50_10bp": "c_jitter",
    "c_merge10": "c_merge", "c_merge50": "c_merge",
    "c_skip10": "c_skip", "c_skip50": "c_skip",
    "c_flip10": "c_flip", "c_flip50": "c_flip",
    "c_ir10": "c_ir", "c_ir50": "c_ir",
    "c_spurious10": "c_spurious", "c_spurious50": "c_spurious",
    "c_trunc3p10": "c_trunc3p", "c_trunc3p50": "c_trunc3p",
    "c_trunc5p10": "c_trunc5p", "c_trunc5p50": "c_trunc5p",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", required=True)
    ap.add_argument("--base", required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--sample", default=None, help="restrict to one library")
    ap.add_argument("--seed", type=int, default=20260831)
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    rows = list(csv.DictReader(open(a.scores), delimiter="\t"))
    by = collections.defaultdict(dict)
    dup = 0
    for r in rows:
        if a.sample and r["sample"] != a.sample:
            continue
        if r["tool"] not in (a.base, a.arm):
            continue
        k = (r["sample"], r["cond"])
        if r["tool"] in by[k]:
            dup += 1
        by[k][r["tool"]] = r
    if dup:
        print(f"WARNING: {dup} duplicate (sample, cond, tool) rows; last one wins")
    ok = lambda r: r is not None and str(r.get("status", "")).startswith("ok")
    paired, missing = [], []
    for k in sorted(by):
        b, x = by[k].get(a.base), by[k].get(a.arm)
        if ok(b) and ok(x):
            paired.append(k)
        else:
            missing.append((k, "-" if b is None else b["status"], "-" if x is None else x["status"]))
    f = lambda r, c: float(r[c])
    d = []
    if not a.quiet:
        print(f"{'sample':38s} {'cond':18s} {a.base + ' Sn/Pr/F1':>22s} {a.arm + ' Sn/Pr/F1':>22s} "
              f"{'dSn':>6s} {'dPr':>6s} {'dF1':>6s}")
    for k in paired:
        b, x = by[k][a.base], by[k][a.arm]
        d.append((k[0], k[1], f(x, "tx_sn") - f(b, "tx_sn"), f(x, "tx_pr") - f(b, "tx_pr"),
                  f(x, "tx_f1") - f(b, "tx_f1")))
        if not a.quiet:
            print(f"{k[0]:38s} {k[1]:18s} {f(b,'tx_sn'):6.1f}/{f(b,'tx_pr'):5.1f}/{f(b,'tx_f1'):6.2f}   "
                  f"{f(x,'tx_sn'):6.1f}/{f(x,'tx_pr'):5.1f}/{f(x,'tx_f1'):6.2f}  "
                  f"{d[-1][2]:+6.2f} {d[-1][3]:+6.2f} {d[-1][4]:+6.2f}")
    if missing:
        print(f"\nMISSING / not-ok cells ({len(missing)}): " +
              "; ".join(f"{k[0]}/{k[1]} base={bs} arm={xs}" for k, bs, xs in missing[:20]) +
              (" ..." if len(missing) > 20 else ""))
    n = len(d)
    if n == 0:
        raise SystemExit("no paired cells")
    m = [statistics.mean(x[i] for x in d) for i in (2, 3, 4)]
    libs = sorted({x[0] for x in d})
    random.seed(a.seed)
    if len(libs) > 1:
        by_lib = collections.defaultdict(list)
        for x in d:
            by_lib[x[0]].append(x[4])
        bs = []
        for _ in range(2000):
            pick = [by_lib[libs[random.randrange(len(libs))]] for _ in libs]
            bs.append(statistics.mean(v for lib in pick for v in lib))
        unit = f"library-block bootstrap over {len(libs)} libraries"
    else:
        bs = [statistics.mean(d[random.randrange(n)][4] for _ in range(n)) for _ in range(2000)]
        unit = "cell bootstrap within ONE library (conditions are not independent replicates)"
    bs.sort()
    print(f"\n{a.arm} vs {a.base}: {n} paired cells over {len(libs)} librar{'y' if len(libs)==1 else 'ies'}; "
          f"mean dSn {m[0]:+.2f} dPr {m[1]:+.2f} dF1 {m[2]:+.2f} 95% CI [{bs[50]:+.2f}, {bs[1949]:+.2f}] "
          f"({unit}); up {sum(1 for x in d if x[4] > 0)} down {sum(1 for x in d if x[4] < 0)} "
          f"same {sum(1 for x in d if x[4] == 0)}; missing {len(missing)}")
    fam = collections.defaultdict(list)
    for x in d:
        fam[FAMILY.get(x[1], "UNMAPPED:" + x[1])].append(x[4])
    for k in sorted(fam):
        print(f"  family {k:14s} n={len(fam[k]):3d} mean dF1 {statistics.mean(fam[k]):+.2f} min {min(fam[k]):+.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
