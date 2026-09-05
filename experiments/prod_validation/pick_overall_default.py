#!/usr/bin/env python3
"""Choose the production default: the configuration that is best OVERALL.

"Overall" is deliberately not a pooled average. The two grids differ in design,
depth and complexity, so a mean over 260 mixed cells is not a quantity anyone
can act on. A configuration qualifies as the default only if it is good on BOTH,
which is what makes it an algorithm rather than a fit:

  1. score every arm on the EXPRESSION contract of each dataset separately
     (the catalog contract rewards echoing the input annotation);
  2. report Sn AND Pr, never F1 alone;
  3. rank by the WORST of the two datasets' F1 (a maximin rule -- a default that
     wins on the development set by giving up the held-out one is exactly the
     failure mode we are trying to leave behind);
  4. show the paired bootstrap of each arm against the incumbent default, per
     dataset, so a tie is visible as a tie.
"""
import argparse
import collections
import csv
import random
import statistics as st
from pathlib import Path

WORK = {"sirv4": Path("experiments/prod_validation/_sirv4grid"),
        "sequin": Path("experiments/prod_validation/_sequingrid")}


def load(dataset):
    p = WORK[dataset] / "expression_contract" / "scores.tsv"
    if not p.exists():                       # SIRV keeps its historical name
        p = WORK[dataset] / "expression_score" / "scores.tsv"
    out = collections.defaultdict(dict)
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["status"].startswith("ok"):
            out[(r["sample"], r["cond"])][r["tool"]] = r
    return out


def boot(d, seed=20260831, n_boot=10000):
    rng = random.Random(seed)
    n = len(d)
    if n == 0:
        return (float("nan"),) * 2
    bs = sorted(sum(d[rng.randrange(n)] for _ in range(n)) / n
                for _ in range(n_boot))
    return bs[int(0.025 * n_boot)], bs[int(0.975 * n_boot)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", action="append", default=[],
                    help="NAME=sirv4_tool,sequin_tool (repeatable)")
    ap.add_argument("--incumbent", default="")
    a = ap.parse_args()
    data = {ds: load(ds) for ds in WORK}
    arms = {}
    for spec in a.arm:
        name, _, tools = spec.partition("=")
        s4, _, sq = tools.partition(",")
        arms[name] = {"sirv4": s4, "sequin": sq}

    stats = {}
    for name, tools in arms.items():
        row = {}
        for ds, tool in tools.items():
            cells = [k for k, v in data[ds].items() if tool in v]
            if not cells:
                row[ds] = None
                continue
            row[ds] = dict(
                n=len(cells),
                sn=st.mean(float(data[ds][k][tool]["tx_sn"]) for k in cells),
                pr=st.mean(float(data[ds][k][tool]["tx_pr"]) for k in cells),
                f1=st.mean(float(data[ds][k][tool]["tx_f1"]) for k in cells))
        stats[name] = row

    print(f"{'arm':26s} " + "".join(
        f"{ds+' Sn':>9s}{ds+' Pr':>9s}{ds+' F1':>9s}" for ds in WORK)
        + f"{'worst F1':>10s}")
    ranked = []
    for name, row in stats.items():
        vals = [row[ds]["f1"] for ds in WORK if row.get(ds)]
        worst = min(vals) if len(vals) == len(WORK) else float("nan")
        ranked.append((worst, name, row))
    for worst, name, row in sorted(ranked, reverse=True,
                                   key=lambda x: (-1e9 if x[0] != x[0] else x[0])):
        cells = ""
        for ds in WORK:
            r = row.get(ds)
            cells += (f"{r['sn']:9.2f}{r['pr']:9.2f}{r['f1']:9.2f}" if r
                      else f"{'-':>9s}{'-':>9s}{'-':>9s}")
        print(f"{name:26s} {cells}{worst:10.2f}")

    if a.incumbent and a.incumbent in arms:
        print(f"\npaired vs incumbent '{a.incumbent}' (dF1, 95% CI), per dataset:")
        for name in arms:
            if name == a.incumbent:
                continue
            out = []
            for ds in WORK:
                ti, tn = arms[a.incumbent][ds], arms[name][ds]
                cells = [k for k, v in data[ds].items() if ti in v and tn in v]
                if not cells:
                    out.append(f"{ds}: n/a")
                    continue
                d = [float(data[ds][k][tn]["tx_f1"]) - float(data[ds][k][ti]["tx_f1"])
                     for k in cells]
                lo, hi = boot(d)
                out.append(f"{ds}: {st.mean(d):+6.2f} [{lo:+.2f},{hi:+.2f}] n={len(d)}")
            print(f"  {name:24s} " + "   ".join(out))


if __name__ == "__main__":
    raise SystemExit(main())
