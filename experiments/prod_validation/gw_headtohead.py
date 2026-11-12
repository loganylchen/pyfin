#!/usr/bin/env python3
"""Genome-wide head-to-head on the cells the campaign has finished.

Every tool -- PyFIN included, now -- assembled from the SAME full BAM, so no
GTF is chromosome-restricted here. Truth is the expression contract built
PER SAMPLE and genome-wide: GENCODE v44 transcripts with NanoCount TPM > 0 in
that sample. Matching is gffcompare's default `=` (intron chain). Sn and Pr
are always reported; F1 is per cell, never derived from mean Sn and mean Pr.

Cells are taken from the campaign's completion markers, not from whatever is
in the output directory, so an unscored or interrupted attempt is never
compared. Scoring lives in its own namespace (`_gwcampaign/h2h/`) so it can
never collide with the chr1/19/21 ladder.
"""
import argparse
import collections
import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, "/SSD/logan/dev/pyfin/experiments/prod_validation")
from expression_rescore import build_truth, expressed_ids, score  # noqa: E402
from gw_cells import finished_cells, scoring_lock  # noqa: E402

BENCH = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results"
             "/gencode_full_sweep_ratio_grid_v2")
CAMP = Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_gwcampaign")
MASTER = Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid/_ref/gencode.v44.gtf")
# Every tool the grid ran. flair.gtf is byte-identical to FLAIR's own
# flair_work/flair.isoforms.gtf in all 44 cells that have it.
COMPETITORS = ("stringtie3", "bambu", "isotools", "isoquant", "espresso",
               "lafite", "talon", "flair")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    # Part A (manifest.tsv / _done) and part B (manifest_b.tsv / _done_b) are
    # two pinned manifests with their own marker directories; both index the
    # same output tree. Part B indices are offset by 1000 in the `idx` column
    # so a row's origin stays readable and the two ranges cannot collide.
    # finished_cells() refuses a missing, edited or overlapping manifest: the
    # table below replaces the previous one, so it must not be a subset of it.
    #
    # One pass at a time, and never under the compactor. The markers are read
    # INSIDE the lock: a pass that listed its cells first and then waited
    # would publish that older list over the table of the pass it waited for.
    h2h = CAMP / "h2h"
    with scoring_lock(CAMP, wait=True):
        (h2h / "_truth").mkdir(parents=True, exist_ok=True)
        status = _score(finished_cells(CAMP), h2h, a.force)
        _drop_scratch(h2h)
        return status


def _drop_scratch(h2h):
    """Remove what a pass leaves behind and nothing reads back.

    The scorer copies each query GTF beside its gffcompare output and
    gffcompare writes an annotated copy: ~200 MB per cell. The cache key is
    (query identity, truth hash, image), not these files, so the cells stay
    cached. cmp.stats, cache_key.json and the tmaps are never touched. Done
    under the scoring lock, so it cannot delete another pass's working copy.
    """
    for name in ("query.gtf", "cmp.annotated.gtf"):
        for path in (h2h / "cells").glob(f"*/*/*/{name}"):
            path.unlink(missing_ok=True)


def _score(cells, h2h, force):
    rows = []
    for idx, sample, cond in cells:
        truth = h2h / "_truth" / f"{sample}.gtf"
        if not truth.exists() or force:
            n_truth = build_truth(MASTER, expressed_ids(BENCH, sample, BENCH), truth)
        else:
            n_truth = sum(1 for line in open(truth) if "\ttranscript\t" in line)
        queries = {"pyfin": CAMP / "gw" / sample / cond / "out" / "assembly.gtf"}
        for tool in COMPETITORS:
            src = BENCH / sample / cond / "assembly" / f"{tool}.gtf"
            if src.exists():
                queries[tool] = src
        for tool, query in queries.items():
            r = score((sample, cond, tool, "gw", str(query), str(truth), n_truth,
                       h2h / "cells" / sample / cond / tool, force))
            r["idx"] = idx
            rows.append(r)

    # published whole or not at all: a scorer that dies part-way leaves the
    # previous table, never a truncated one that reads as a finished pass
    out = h2h / "scores.tsv"
    tmp = out.with_name(out.name + ".partial")
    with open(tmp, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted({k for r in rows for k in r}), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    os.replace(tmp, out)

    by_cell = collections.defaultdict(dict)
    for r in rows:
        if str(r.get("status", "")).startswith("ok"):
            by_cell[(r["idx"], r["sample"], r["cond"])][r["tool"]] = r
    print(f"{'cell':44s} {'tool':11s} {'Sn':>6s} {'Pr':>6s} {'F1':>6s} {'n_query':>8s} {'truth':>6s}")
    wins = collections.defaultdict(lambda: [0, 0])
    for (idx, sample, cond), tools in sorted(by_cell.items()):
        p = tools.get("pyfin")
        label = f"{idx:2d} {sample.replace('SGNex_', '').replace('_directRNA_', ' ')}/{cond}"
        for tool in ("pyfin",) + COMPETITORS:
            r = tools.get(tool)
            if r is None:
                continue
            mark = ""
            if p is not None and tool != "pyfin":
                wins[tool][1] += 1
                if float(p["tx_f1"]) > float(r["tx_f1"]):
                    wins[tool][0] += 1
                mark = f"  dF1 {float(p['tx_f1']) - float(r['tx_f1']):+6.2f}"
            print(f"{label if tool == 'pyfin' else '':44s} {tool:11s} {float(r['tx_sn']):6.2f} "
                  f"{float(r['tx_pr']):6.2f} {float(r['tx_f1']):6.2f} {r.get('n_query', ''):>8} "
                  f"{r['truth_transcripts']:>6}{mark}")
        print()
    print(f"{'tool':11s} {'cells':>5s} {'pyfin wins (F1)':>16s}")
    for tool, (w, n) in sorted(wins.items()):
        print(f"{tool:11s} {n:5d} {f'{w}/{n}':>16s}")
    print(f"\nwritten: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
