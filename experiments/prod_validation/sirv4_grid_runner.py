#!/usr/bin/env python3
"""Run PyFIN across the SIRV4 ratio-grid benchmark and score it exactly like
the competitors were scored.

The external benchmark tree
(``.../results/sirv4_full_sweep_ratio_grid_v2``) evaluates 8 assemblers over
6 samples x 26 annotation conditions but contains **no PyFIN cells**. This
script fills that gap under the identical contract:

  * same BAM per sample,
  * same per-condition input annotation (``_ref/<sample>/<cond>/annotation.gtf``),
  * same truth (``SIRV_Set4.gtf``),
  * same scorer (gffcompare 0.12.6, the version the tree was built with).

Only PyFIN is run here; competitor numbers are read from the tree as-is.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# --- dataset registry -------------------------------------------------------
# The grid contract (per-sample BAM/fastq/blow5, per-condition input annotation,
# one fixed structural truth, gffcompare 0.12.6) is dataset-INDEPENDENT: only
# these paths change. Keeping one runner over a registry (instead of a forked
# copy per spike-in) is what makes a change testable on every dataset at once.
_ROOT = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex")
DATASETS = {
    "sirv4": dict(
        bench=_ROOT / "results" / "sirv4_full_sweep_ratio_grid_v2",
        work=Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid"),
        truth_src=_ROOT / "refs" / "sirv" / "SIRV_Set4.gtf",
        genome_src=_ROOT / "refs" / "sirv" / "SIRV_Set4.fa",
        truth_name="SIRV_Set4.gtf", genome_name="SIRV_Set4.fa",
        concentration=None,
    ),
    "sequin": dict(
        bench=_ROOT / "results" / "sequin_v1_full_sweep_ratio_grid_v2",
        work=Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid"),
        truth_src=(_ROOT / "results" / "sequin_v1_full_sweep_ratio_grid_v2"
                   / "_ref" / "_source" / "annotation.gtf"),
        genome_src=_ROOT / "refs" / "sequin" / "sequin_chrIS.fa",
        truth_name="sequin_chrIS_truth.gtf", genome_name="sequin_chrIS.fa",
        # Sequin ships an independent CONCENTRATION truth (attomol/ul), which
        # SIRV4 does not: it is the held-out quantification axis.
        concentration=_ROOT / "refs" / "concentration" / "RNAsequins_MixA_v1.tsv",
    ),
    # Dual spike-in (SIRV Set1 + Sequin) in HEYA8: 234 transcripts over chrIS and
    # SIRV1-7, 5 samples. Third design, so a rule must now hold on three grids.
    "heya8": dict(
        bench=_ROOT / "results" / "heya8_dual_full_sweep_ratio_grid_v2",
        work=Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_heya8grid"),
        truth_src=(_ROOT / "results" / "heya8_dual_full_sweep_ratio_grid_v2"
                   / "_ref" / "_source" / "annotation.gtf"),
        genome_src=_ROOT / "refs" / "spikein_sequin_sirv1.fa",
        truth_name="heya8_dual_truth.gtf", genome_name="spikein_sequin_sirv1.fa",
        concentration=None,
    ),
}
DATASET = os.environ.get("PYFIN_GRID_DATASET", "sirv4")
_D = DATASETS[DATASET]
BENCH = _D["bench"]
# gffcompare cannot read the /autofs automount from inside Docker, so the
# truth is staged on /SSD once and read from there.
WORK = _D["work"]
TRUTH = WORK / "_ref" / _D["truth_name"]
GENOME = WORK / "_ref" / _D["genome_name"]


def stage_refs() -> None:
    """Copy this dataset's truth GTF + genome FASTA onto /SSD (docker-visible)."""
    (WORK / "_ref").mkdir(parents=True, exist_ok=True)
    for src, dst in ((_D["truth_src"], TRUTH), (_D["genome_src"], GENOME)):
        if not dst.exists() or dst.stat().st_size != Path(src).stat().st_size:
            shutil.copyfile(src, dst)
    fai = Path(str(_D["genome_src"]) + ".fai")
    if fai.exists() and not Path(str(GENOME) + ".fai").exists():
        shutil.copyfile(fai, str(GENOME) + ".fai")
SIF = Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_img/"
           "pyfin_gpu_e268c9b.sif")
GFFCOMPARE_IMAGE = "quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"
# Source tree the SIF runs (PYTHONPATH). Override with PYFIN_REPO to score a
# different commit (e.g. a `git worktree` at HEAD) without touching this file.
REPO = os.environ.get("PYFIN_REPO", "/SSD/logan/dev/pyfin")
# Optional extra site-packages prepended to PYTHONPATH inside the SIF (e.g. a newer
# krill wheel installed with `pip install --target`). Set PYFIN_PYPATH_PREFIX.
PYPATH_PREFIX = os.environ.get("PYFIN_PYPATH_PREFIX", "")


def samples() -> list[str]:
    return sorted(p.name for p in BENCH.glob("SGNex_*") if p.is_dir())


def conditions() -> list[str]:
    gc = BENCH / "benchmark" / "gffcompare"
    first = next(iter(sorted(gc.glob("SGNex_*"))), None)
    if first is None:
        return []
    return sorted(p.name for p in first.glob("*") if p.is_dir())


def n_annotation_tx(gtf: Path) -> int:
    if not gtf.exists():
        return 0
    ids = set()
    with open(gtf) as fh:
        for line in fh:
            if 'transcript_id "' in line:
                ids.add(line.split('transcript_id "')[1].split('"')[0])
    return len(ids)


def parse_stats(path: Path) -> dict:
    if not path.exists():
        return {}
    txt = path.read_text()
    out: dict = {}
    for label, key in (("Transcript level", "tx"), ("Intron chain level", "ic"),
                       ("Locus level", "loc")):
        m = re.search(rf"{re.escape(label)}:\s*([\d.]+)\s*\|\s*([\d.]+)", txt)
        if m:
            sn, pr = float(m.group(1)), float(m.group(2))
            out[f"{key}_sn"], out[f"{key}_pr"] = sn, pr
            out[f"{key}_f1"] = 0.0 if sn + pr == 0 else 2 * sn * pr / (sn + pr)
    for pat, key in ((r"Query mRNAs\s*:\s*(\d+)", "n_query"),
                     (r"Matching transcripts:\s*(\d+)", "n_match")):
        m = re.search(pat, txt)
        if m:
            out[key] = int(m.group(1))
    return out


def run_cell(args) -> dict:
    sample, cond, outroot, extra, tag, extra_gtf, profile = args
    t0 = time.time()
    bam = BENCH / sample / "align" / f"{sample}.sorted.bam"
    fq = BENCH / sample / "subset" / "mapped.fq.gz"
    sig = BENCH / sample / "subset" / "mapped.blow5"
    anno = BENCH / "_ref" / sample / cond / "annotation.gtf"
    cell = Path(outroot) / sample / cond
    row = {"sample": sample, "cond": cond, "tool": f"pyfin{tag}", "status": "?"}

    for req in (bam, fq, sig):
        if not req.exists():
            row["status"] = f"missing_input:{req.name}"
            return row

    if cell.exists():
        shutil.rmtree(cell)
    cell.mkdir(parents=True, exist_ok=True)

    cmd = ["singularity", "exec", "-B", "/autofs,/SSD", str(SIF),
           "env", f"PYTHONPATH={PYPATH_PREFIX + ':' if PYPATH_PREFIX else ''}{REPO}", "/usr/bin/python3.10", "-m", "fin.cli",
           "--bam", str(bam), "--genome", str(GENOME),
           "--fastq", str(fq), "--signal", str(sig),
           "--output-dir", str(cell / "out"),
           "--profile", profile, "--no-gpu", "--threads", "2"]
    # p00 ships a zero-transcript annotation: that is the de novo cell, so no
    # --gtf is passed at all rather than an empty one.
    n_anno = n_annotation_tx(anno)
    row["n_anno"] = n_anno
    if n_anno > 0:
        cmd += ["--gtf", str(anno)]
        # Flags that are only meaningful with an annotation (e.g.
        # --catalog-mode) would fail loudly on the de novo p00 cell, which is
        # correct behaviour -- so they are simply not passed there.
        cmd += list(extra_gtf)
    cmd += list(extra)

    log = cell / "run.log"
    with open(log, "w") as fh:
        rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT).returncode
    gtf = cell / "out" / "assembly.gtf"
    if rc != 0 or not gtf.exists():
        row["status"] = f"pyfin_rc{rc}"
        row["runtime_s"] = round(time.time() - t0, 1)
        return row

    gcdir = cell / "gc"
    gcdir.mkdir(exist_ok=True)
    gcmd = ["docker", "run", "--rm", "-v", "/SSD:/SSD", GFFCOMPARE_IMAGE,
            "gffcompare", "-r", str(TRUTH), "-o", str(gcdir / "cmp"), str(gtf)]
    grc = subprocess.run(gcmd, capture_output=True).returncode
    stats = parse_stats(gcdir / "cmp.stats")
    if grc != 0 or not stats:
        row["status"] = f"gffcompare_rc{grc}"
        row["runtime_s"] = round(time.time() - t0, 1)
        return row

    row.update(stats)
    row["status"] = "ok"
    row["runtime_s"] = round(time.time() - t0, 1)
    return row


COLS = ["sample", "cond", "tool", "status", "n_anno", "n_query", "n_match",
        "tx_sn", "tx_pr", "tx_f1", "ic_sn", "ic_pr", "ic_f1",
        "loc_sn", "loc_pr", "loc_f1", "runtime_s"]


def main() -> int:
    # Rebound below BEFORE the pool forks, so workers inherit the dataset paths.
    global DATASET, _D, BENCH, WORK, TRUTH, GENOME
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="output root for run cells")
    ap.add_argument("--tsv", required=True, help="metrics table to write")
    ap.add_argument("--tag", default="", help="suffix for the tool name")
    ap.add_argument("--samples", default="", help="comma list (default: all)")
    ap.add_argument("--conds", default="", help="comma list (default: all)")
    ap.add_argument("--dataset", default=DATASET, choices=sorted(DATASETS),
                    help="which benchmark grid to run (default: $PYFIN_GRID_DATASET or sirv4)")
    ap.add_argument("--profile", default="sirv",
                    help="pyfin --profile for every cell (sirv | real-drna | "
                         "real-drna-precision | custom). The SIRV working point is a "
                         "six-H9-sample FIT: running another dataset under it measures "
                         "how much of the score is threshold transfer.")
    ap.add_argument("--workers", type=int, default=6)
    # NOTE: pass these with the `=` form (`--extra=--flag`) when the value is a
    # SINGLE token starting with `-`. argparse only accepts a leading-dash value
    # when it contains a space, so `--extra "--flag"` fails while
    # `--extra "--flag --other"` silently works -- a trap worth naming here.
    ap.add_argument("--extra", default="",
                    help="extra pyfin flags as one shell-quoted string "
                         "(use --extra=... when the value starts with '-')")
    ap.add_argument("--extra-with-gtf", default="",
                    help="extra flags applied ONLY to cells that get a --gtf "
                         "(use --extra-with-gtf=... when the value starts with '-')")
    a = ap.parse_args()
    DATASET = a.dataset
    _D = DATASETS[DATASET]
    BENCH, WORK = _D["bench"], _D["work"]
    TRUTH = WORK / "_ref" / _D["truth_name"]
    GENOME = WORK / "_ref" / _D["genome_name"]
    stage_refs()
    print(f"dataset={DATASET} bench={BENCH}", flush=True)

    smp = [s for s in (a.samples.split(",") if a.samples else samples()) if s]
    cds = [c for c in (a.conds.split(",") if a.conds else conditions()) if c]
    extra = shlex.split(a.extra)
    extra_gtf = shlex.split(a.extra_with_gtf)
    jobs = [(s, c, a.out, extra, a.tag, extra_gtf, a.profile)
            for s in smp for c in cds]
    print(f"cells={len(jobs)} samples={len(smp)} conds={len(cds)} "
          f"profile={a.profile} repo={REPO} "
          f"workers={a.workers} extra={extra} extra_with_gtf={extra_gtf}",
          flush=True)

    Path(a.tsv).parent.mkdir(parents=True, exist_ok=True)
    rows, done = [], 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(run_cell, j): j for j in jobs}
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            done += 1
            if done % 10 == 0 or r["status"] != "ok":
                print(f"[{done}/{len(jobs)}] {r['sample'][-12:]}/{r['cond']}"
                      f" {r['status']} f1={r.get('tx_f1','-')}", flush=True)
            with open(a.tsv, "w") as fh:
                w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t",
                                   extrasaction="ignore")
                w.writeheader()
                for x in sorted(rows, key=lambda z: (z["sample"], z["cond"])):
                    w.writerow(x)
    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"done: {ok}/{len(rows)} ok -> {a.tsv}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
