#!/usr/bin/env python3
"""Human head-to-head on the GENCODE grid.

The GENCODE sweep is filling in: 541 competitor assemblies now exist where
there were 4. This scores every competitor cell that overlaps a condition we
have a PyFIN run for, on the one sample whose reads are staged locally
(H9 replicate2 run2).

Two caveats that belong with every number this prints. The competitors
assembled GENOME-WIDE while our runs use the chr1/19/21 subset BAM, so their
GTF is restricted to those chromosomes before scoring -- reads elsewhere cannot
create transcripts here, but any threshold they applied was computed over the
whole genome. And the tools with coverage on this grid are not the tools that
are strongest on the spike-ins: espresso and isoquant appear once each.
"""
import argparse, collections, csv, sys
from pathlib import Path
sys.path.insert(0, "/SSD/logan/dev/pyfin/experiments/prod_validation")
from expression_rescore import GRIDS, build_truth, expressed_ids, score

BENCH = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results"
             "/gencode_full_sweep_ratio_grid_v2")
SAMPLE = "SGNex_H9_directRNA_replicate2_run2"
KEEP = {"chr1", "chr19", "chr21"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pyfin-label", default="GCSMC")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    g = GRIDS["gencode_sub"]
    work = g["work"]
    out = work / "expression_contract"
    tmp = work / "_headtohead"
    tmp.mkdir(exist_ok=True)
    truth = out / "_truth" / f"{SAMPLE}.gtf"
    n_truth = build_truth(work / "_ref" / g["truth_name"],
                          expressed_ids(g["bench"], SAMPLE, g.get("nanocount")),
                          truth)

    pyf = collections.defaultdict(dict)
    for r in csv.DictReader(open(out / "scores.tsv"), delimiter="\t"):
        if r["sample"] == SAMPLE and r["status"].startswith("ok"):
            pyf[r["cond"]][r["tool"]] = r
    ours = {c: v[a.pyfin_label] for c, v in pyf.items() if a.pyfin_label in v}
    print(f"truth transcripts (chr1/19/21, NanoCount TPM>0): {n_truth}")
    print(f"pyfin cells available ({a.pyfin_label}): {len(ours)}\n")

    rows = []
    for cond in sorted(ours):
        for src in sorted((BENCH / SAMPLE / cond / "assembly").glob("*.gtf")):
            tool = src.stem
            dst = tmp / f"{tool}.{cond}.chr1_19_21.gtf"
            n_all = n_keep = 0
            with open(src) as fh, open(dst, "w") as o:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    n_all += 1
                    if line.split("\t", 1)[0] in KEEP:
                        n_keep += 1
                        o.write(line)
            r = score((SAMPLE, cond, tool, "competitor", str(dst), str(truth),
                       n_truth, tmp / "cells" / cond / tool, a.force))
            r["kept_frac"] = 100 * n_keep / max(1, n_all)
            rows.append(r)

    print(f"{'cond':18s} {'tool':12s} {'Sn':>6s} {'Pr':>6s} {'F1':>6s} "
          f"{'n_query':>8s}   {'dF1 vs pyfin':>12s}")
    per_tool = collections.defaultdict(list)
    for r in sorted(rows, key=lambda x: (x["cond"], x["tool"])):
        if not str(r["status"]).startswith("ok"):
            print(f"{r['cond']:18s} {r['tool']:12s} {r['status']}")
            continue
        p = ours[r["cond"]]
        d = float(p["tx_f1"]) - float(r["tx_f1"])
        per_tool[r["tool"]].append((r, p))
        print(f"{r['cond']:18s} {r['tool']:12s} {float(r['tx_sn']):6.2f} "
              f"{float(r['tx_pr']):6.2f} {float(r['tx_f1']):6.2f} "
              f"{r.get('n_query',''):>8}   {d:+12.2f}")
        print(f"{'':18s} {'pyfin':12s} {float(p['tx_sn']):6.2f} "
              f"{float(p['tx_pr']):6.2f} {float(p['tx_f1']):6.2f} "
              f"{p['n_query']:>8}")

    print(f"\n{'tool':12s} {'cells':>5s} {'their Sn/Pr/F1':>22s} {'pyfin Sn/Pr/F1':>22s} "
          f"{'mean dF1':>9s}  pyfin wins")
    for t, pairs in sorted(per_tool.items()):
        n = len(pairs)
        m = lambda f, k: sum(float(x[f][k]) for x in pairs) / n
        wins = sum(1 for r, p in pairs if float(p["tx_f1"]) > float(r["tx_f1"]))
        print(f"{t:12s} {n:5d} "
              f"{m(0,'tx_sn'):7.2f}/{m(0,'tx_pr'):6.2f}/{m(0,'tx_f1'):6.2f}   "
              f"{m(1,'tx_sn'):7.2f}/{m(1,'tx_pr'):6.2f}/{m(1,'tx_f1'):6.2f}   "
              f"{m(1,'tx_f1')-m(0,'tx_f1'):+8.2f}  {wins}/{n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
