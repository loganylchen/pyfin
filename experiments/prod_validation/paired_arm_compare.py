"""Paired comparison of two pyfin labels on one grid's EXPRESSION-contract scores.

Prints per-condition mean deltas (Sn / Pr / F1) and the overall paired
bootstrap CI of the mean F1 delta (seed 20260831). Cells are paired by
(sample, cond); a cell missing from either label is dropped and reported.
"""
import argparse, collections, csv, random
from pathlib import Path

ROOT = Path("/SSD/logan/dev/pyfin/experiments/prod_validation")
GRID = {"sequin": "_sequingrid", "sirv4": "_sirv4grid", "heya8": "_heya8grid"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=sorted(GRID))
    ap.add_argument("--base", required=True)
    ap.add_argument("--arm", required=True)
    ap.add_argument("--boot", type=int, default=5000)
    a = ap.parse_args()
    p = ROOT / GRID[a.dataset] / "expression_contract" / "scores.tsv"
    by = collections.defaultdict(dict)
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["status"].startswith("ok") and r["tool"] in (a.base, a.arm):
            by[(r["sample"], r["cond"])][r["tool"]] = r
    pairs = {k: v for k, v in by.items() if a.base in v and a.arm in v}
    print(f"{a.dataset}: {len(pairs)} paired cells ({len(by) - len(pairs)} unpaired dropped)")
    per = collections.defaultdict(list)
    for (s, c), v in pairs.items():
        b, x = v[a.base], v[a.arm]
        per[c].append(tuple(float(x[k]) - float(b[k]) for k in ("tx_sn", "tx_pr", "tx_f1")))
    print(f"{'cond':14s} {'n':>2s} {'dSn':>7s} {'dPr':>7s} {'dF1':>7s}")
    for c in sorted(per):
        d = per[c]; n = len(d)
        print(f"{c:14s} {n:2d} {sum(x[0] for x in d)/n:+7.2f} {sum(x[1] for x in d)/n:+7.2f} {sum(x[2] for x in d)/n:+7.2f}")
    allf = [x for d in per.values() for x in d]
    n = len(allf)
    mean = [sum(x[i] for x in allf) / n for i in range(3)]
    rng = random.Random(20260831)
    boots = []
    for _ in range(a.boot):
        s = [allf[rng.randrange(n)] for _ in range(n)]
        boots.append(sum(x[2] for x in s) / n)
    boots.sort()
    lo, hi = boots[int(0.025 * a.boot)], boots[int(0.975 * a.boot) - 1]
    up = sum(1 for x in allf if x[2] > 1e-9); dn = sum(1 for x in allf if x[2] < -1e-9)
    print(f"OVERALL dSn {mean[0]:+.2f} dPr {mean[1]:+.2f} dF1 {mean[2]:+.2f} "
          f"CI95 [{lo:+.2f},{hi:+.2f}]  cells up/down/tie {up}/{dn}/{n-up-dn}")


if __name__ == "__main__":
    main()
