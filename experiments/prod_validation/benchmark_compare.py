#!/usr/bin/env python3
"""Compare PyFIN against the benchmark's competitors on either grid.

Reads the competitor gffcompare stats the benchmark tree already contains and
the PyFIN metrics TSVs written by ``sirv4_grid_runner.py``, and reports Sn AND
Pr (never F1 alone) per tool, per condition, plus how often each tool is top or
tied on a cell.

CONTRACT WARNING: this is the CATALOG contract -- truth is every annotated
transcript, expressed or not. In `full`/`p*` cells the supplied annotation IS
(a subset of) the truth, so a tool that simply echoes its input scores
perfectly. Bambu does exactly that. Read this table together with the
expression contract, which subsets truth to what is actually expressed.
"""
import argparse
import collections
import csv
import re
import statistics as st
from pathlib import Path

ROOT = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results")
GRIDS = {"sirv4": ROOT / "sirv4_full_sweep_ratio_grid_v2",
         "sequin": ROOT / "sequin_v1_full_sweep_ratio_grid_v2"}


def parse_stats(path: Path):
    if not path.exists():
        return None
    m = re.search(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)", path.read_text())
    if not m:
        return None
    sn, pr = float(m.group(1)), float(m.group(2))
    return sn, pr, (0.0 if sn + pr == 0 else 2 * sn * pr / (sn + pr))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="sequin", choices=sorted(GRIDS))
    ap.add_argument("--pyfin", action="append", default=[],
                    help="LABEL=path/to/metrics.tsv (repeatable)")
    ap.add_argument("--conds", default="", help="comma list; default all")
    a = ap.parse_args()

    gc = GRIDS[a.dataset] / "benchmark" / "gffcompare"
    cells = collections.defaultdict(dict)          # (sample, cond) -> tool -> triple
    for sdir in sorted(gc.glob("SGNex_*")):
        for cdir in sorted(p for p in sdir.glob("*") if p.is_dir()):
            for tdir in sorted(p for p in cdir.glob("*") if p.is_dir()):
                v = parse_stats(tdir / "cmp.stats")
                if v:
                    cells[(sdir.name, cdir.name)][tdir.name] = v
    for spec in a.pyfin:
        label, _, path = spec.partition("=")
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["status"] != "ok":
                continue
            cells[(r["sample"], r["cond"])][label] = (
                float(r["tx_sn"]), float(r["tx_pr"]), float(r["tx_f1"]))

    keys = sorted(cells)
    if a.conds:
        want = set(a.conds.split(","))
        keys = [k for k in keys if k[1] in want]
    tools = sorted({t for k in keys for t in cells[k]})
    # only cells where every tool is present, so the comparison is paired
    paired = [k for k in keys if all(t in cells[k] for t in tools)]
    print(f"dataset={a.dataset}  cells={len(keys)}  fully paired={len(paired)}  "
          f"tools={len(tools)}")
    print(f"\n{'tool':14s} {'cells':>5s} {'Sn':>6s} {'Pr':>6s} {'F1':>6s} {'top/tied':>9s}")
    best = {k: max(v[2] for v in cells[k].values()) for k in paired}
    rows = []
    for t in tools:
        vals = [cells[k][t] for k in paired]
        rows.append((st.mean(v[2] for v in vals), t, len(vals),
                     st.mean(v[0] for v in vals), st.mean(v[1] for v in vals),
                     sum(1 for k in paired if cells[k][t][2] >= best[k] - 1e-9)))
    for f1, t, n, sn, pr, top in sorted(rows, reverse=True):
        print(f"{t:14s} {n:5d} {sn:6.2f} {pr:6.2f} {f1:6.2f} {top:5d}/{len(paired)}")

    print("\nper-condition Sn/Pr:")
    conds = sorted({k[1] for k in paired})
    print(f"{'cond':20s} " + " ".join(f"{t[:11]:>12s}" for t in tools))
    for c in conds:
        ks = [k for k in paired if k[1] == c]
        cellsv = []
        for t in tools:
            sn = st.mean(cells[k][t][0] for k in ks)
            pr = st.mean(cells[k][t][1] for k in ks)
            cellsv.append(f"{sn:5.1f}/{pr:5.1f}")
        print(f"{c:20s} " + " ".join(f"{v:>12s}" for v in cellsv))


if __name__ == "__main__":
    raise SystemExit(main())
