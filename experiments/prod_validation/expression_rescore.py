#!/usr/bin/env python3
"""Score every tool on the EXPRESSION contract, on either grid.

The catalog contract (truth = every annotated transcript) rewards echoing the
input annotation: in a `full` cell the supplied annotation IS the truth, so a
tool that re-emits it scores ~100/100 without discovering anything. Bambu and
TALON do exactly that.

The expression contract subsets the structural truth to the transcripts the
sample actually expresses (NanoCount TPM > 0 on the same reads), so re-emitting
164 annotated models when 30 are expressed is paid for in precision. The
pseudo-truth is used ONLY for scoring and is never an input to any tool.

Both PyFIN run dirs and the benchmark's competitor GTFs are scored with the same
gffcompare container the tree was built with.
"""
import argparse
import collections
import csv
import re
import shutil
import statistics as st
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results")
GRIDS = {
    "sirv4": dict(bench=ROOT / "sirv4_full_sweep_ratio_grid_v2",
                  work=Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid"),
                  truth_name="SIRV_Set4.gtf"),
    "sequin": dict(bench=ROOT / "sequin_v1_full_sweep_ratio_grid_v2",
                   work=Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid"),
                   truth_name="sequin_chrIS_truth.gtf"),
}
IMAGE = "quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"
COLS = ["sample", "cond", "tool", "family", "status", "truth_transcripts",
        "n_query", "n_match", "tx_sn", "tx_pr", "tx_f1"]


def expressed_ids(bench: Path, sample: str) -> set:
    path = bench / "_ref" / sample / "expression" / "nanocount.tsv"
    out = set()
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if float(row["tpm"]) > 0.0:
                out.add(row["transcript_name"])
    return out


def build_truth(master: Path, keep: set, dest: Path) -> int:
    seen = set()
    with open(master) as src, open(dest, "w") as out:
        for line in src:
            if line.startswith("#"):
                continue
            m = re.search(r'transcript_id "([^"]+)"', line)
            if m and m.group(1) in keep:
                seen.add(m.group(1))
                out.write(line)
    return len(seen)


def parse_stats(path: Path):
    if not path.exists():
        return {}
    txt = path.read_text()
    out = {}
    m = re.search(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)", txt)
    if m:
        sn, pr = float(m.group(1)), float(m.group(2))
        out.update(tx_sn=sn, tx_pr=pr,
                   tx_f1=0.0 if sn + pr == 0 else 2 * sn * pr / (sn + pr))
    for pat, key in ((r"Query mRNAs\s*:\s*(\d+)", "n_query"),
                     (r"Matching transcripts:\s*(\d+)", "n_match")):
        m = re.search(pat, txt)
        if m:
            out[key] = int(m.group(1))
    return out


def score(job):
    sample, cond, tool, family, query, truth, n_truth, cell, force = job
    row = dict(sample=sample, cond=cond, tool=tool, family=family,
               status="?", truth_transcripts=n_truth)
    if not Path(query).exists():
        row["status"] = "missing_query"
        return row
    cell = Path(cell)
    if force and cell.exists():
        shutil.rmtree(cell)
    cell.mkdir(parents=True, exist_ok=True)
    local = cell / "query.gtf"
    if not (cell / "cmp.stats").exists():
        shutil.copyfile(query, local)
        rc = subprocess.run(
            ["docker", "run", "--rm", "-v", "/SSD:/SSD", IMAGE, "gffcompare",
             "-r", str(truth), "-o", str(cell / "cmp"), str(local)],
            capture_output=True).returncode
        if rc != 0:
            row["status"] = f"gffcompare_rc{rc}"
            return row
    stats = parse_stats(cell / "cmp.stats")
    if not stats:
        row["status"] = "no_stats"
        return row
    row.update(stats)
    if "tx_f1" not in row:
        # gffcompare omits the transcript level when the query has no usable
        # transcripts: the tool emitted nothing for this cell. Recall is 0 and
        # there is nothing correct among the predictions, so both are 0 -- and
        # the status records WHY, so it is never confused with a real score.
        row.update(tx_sn=0.0, tx_pr=0.0, tx_f1=0.0)
        row["status"] = "ok_empty_query"
        return row
    row["status"] = "ok"
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="sequin", choices=sorted(GRIDS))
    ap.add_argument("--pyfin", action="append", default=[],
                    help="LABEL=run_root (repeatable)")
    ap.add_argument("--out", default="")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--force", action="store_true",
                    help="rescore even when a cached cmp.stats exists. Cells are "
                         "cached by (sample, cond, TOOL LABEL), so reusing a label "
                         "for a different run silently returns the old numbers.")
    a = ap.parse_args()
    g = GRIDS[a.dataset]
    bench, work = g["bench"], g["work"]
    out_dir = Path(a.out) if a.out else work / "expression_contract"
    out_dir.mkdir(parents=True, exist_ok=True)
    master = work / "_ref" / g["truth_name"]
    truth_dir = out_dir / "_truth"
    truth_dir.mkdir(exist_ok=True)

    samples = sorted(p.name for p in bench.glob("SGNex_*") if p.is_dir())
    truths = {}
    for s in samples:
        dest = truth_dir / f"{s}.gtf"
        n = build_truth(master, expressed_ids(bench, s), dest)
        truths[s] = (dest, n)
    print("expressed truth per sample:",
          {s: truths[s][1] for s in samples}, flush=True)

    jobs = []
    for s in samples:
        truth, n_truth = truths[s]
        for cdir in sorted(p for p in (bench / s).glob("*") if p.is_dir()):
            cond = cdir.name
            for gtf in sorted((cdir / "assembly").glob("*.gtf")):
                tool = gtf.stem
                if tool.startswith("cmp."):
                    continue
                jobs.append((s, cond, tool, "competitor", str(gtf), str(truth),
                             n_truth, out_dir / "cells" / s / cond / tool, a.force))
        for spec in a.pyfin:
            label, _, root = spec.partition("=")
            for q in sorted(Path(root).glob(f"{s}/*/out/assembly.gtf")):
                cond = q.parents[1].name
                jobs.append((s, cond, label, "pyfin", str(q), str(truth),
                             n_truth, out_dir / "cells" / s / cond / label, a.force))
    print(f"scoring {len(jobs)} cells", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(score, j) for j in jobs]
        for i, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if i % 200 == 0:
                print(f"  {i}/{len(jobs)}", flush=True)
    rows.sort(key=lambda r: (r["sample"], r["cond"], r["family"], r["tool"]))
    with open(out_dir / "scores.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    ok = [r for r in rows if r["status"] in ("ok", "ok_empty_query")]
    n_empty = sum(1 for r in ok if r["status"] == "ok_empty_query")
    print(f"ok {len(ok)}/{len(rows)} ({n_empty} empty-query cells) "
          f"-> {out_dir/'scores.tsv'}")
    by = collections.defaultdict(dict)
    for r in ok:
        by[(r["sample"], r["cond"])][r["tool"]] = r
    tools = sorted({r["tool"] for r in ok})
    paired = [k for k in by if all(t in by[k] for t in tools)]
    print(f"\nEXPRESSION contract, {len(paired)} fully paired cells")
    print(f"{'tool':16s} {'Sn':>6s} {'Pr':>6s} {'F1':>6s} {'top/tied':>10s}")
    best = {k: max(float(by[k][t]["tx_f1"]) for t in tools) for k in paired}
    out = []
    for t in tools:
        v = [by[k][t] for k in paired]
        out.append((st.mean(float(x["tx_f1"]) for x in v), t,
                    st.mean(float(x["tx_sn"]) for x in v),
                    st.mean(float(x["tx_pr"]) for x in v),
                    sum(1 for k in paired if float(by[k][t]["tx_f1"]) >= best[k] - 1e-9)))
    for f1, t, sn, pr, top in sorted(out, reverse=True):
        print(f"{t:16s} {sn:6.2f} {pr:6.2f} {f1:6.2f} {top:6d}/{len(paired)}")


if __name__ == "__main__":
    raise SystemExit(main())
