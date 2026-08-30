#!/usr/bin/env python3
"""Fit and evaluate the calibrated candidate ranker (offline, leakage-guarded).

Training data: candidate_evidence.tsv (inference-time observable features)
joined with gffcompare exact-match labels. GENCODE class codes and NanoCount
counts appear ONLY as labels/evaluation, never as features.

Discipline:
* fit on the tuning sample only (H9 r2r2);
* grouped cross-validation by chromosome (no locus leakage);
* the operating point is chosen on the tuning frontier under the constraint
  "simulated T1 honest F1 >= current T1" and frozen BEFORE r3r1 is scored;
* r3r1 is evaluated once with the frozen model+threshold via the real
  gffcompare/truth/scorer contract (no refit, no re-tuning).

Model: L2 logistic regression in pure NumPy on standardized features; the
frozen coefficients are emitted as JSON for the in-pipeline scorer.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parents[2]))
import profile_sweep as ps  # noqa: E402

# Single source of truth for the feature transform: the production module.
# Training and inference features cannot drift apart.
from fin.analysis.candidate_ranking import RANKER_V1, featurize  # noqa: E402

FEATURES = RANKER_V1["features"]


def load_dataset(evidence_tsv: Path, tmap: Path):
    labels = {}
    with open(tmap) as handle:
        handle.readline()
        for line in handle:
            x = line.rstrip("\n").split("\t")
            labels[x[4]] = x[2]
    X, y, ids, chroms = [], [], [], []
    with open(evidence_tsv) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            cls = labels.get(row["candidate_id"])
            if cls is None:  # pre-snap row merged away before write
                continue
            X.append(featurize(row))
            y.append(1.0 if cls == "=" else 0.0)
            ids.append(row["candidate_id"])
            chroms.append(row["chrom"])
    return np.array(X), np.array(y), ids, chroms


def standardize(X, mean=None, std=None):
    if mean is None:
        mean = X.mean(axis=0)
        std = X.std(axis=0)
        std[std == 0] = 1.0
    return (X - mean) / std, mean, std


def fit_logistic(X, y, l2=1.0, lr=0.1, iters=3000, seed=7):
    rng = np.random.default_rng(seed)
    w = rng.normal(0, 0.01, X.shape[1])
    b = 0.0
    n = len(y)
    for _ in range(iters):
        z = X @ w + b
        p = 1.0 / (1.0 + np.exp(-z))
        gw = X.T @ (p - y) / n + l2 * w / n
        gb = float(np.mean(p - y))
        w -= lr * gw
        b -= lr * gb
    return w, b


def auc(scores, y):
    order = np.argsort(scores)
    ranks = np.empty(len(scores))
    ranks[order] = np.arange(1, len(scores) + 1)
    pos = y == 1
    n1, n0 = pos.sum(), (~pos).sum()
    if n1 == 0 or n0 == 0:
        return float("nan")
    return (ranks[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def grouped_cv_auc(X, y, chroms, folds=4, l2=1.0):
    uniq = sorted(set(chroms))
    assign = {c: i % folds for i, c in enumerate(uniq)}
    fold_of = np.array([assign[c] for c in chroms])
    aucs = []
    for k in range(folds):
        tr, te = fold_of != k, fold_of == k
        Xs, mean, std = standardize(X[tr])
        w, b = fit_logistic(Xs, y[tr], l2=l2)
        Xe = (X[te] - mean) / std
        aucs.append(auc(Xe @ w + b, y[te]))
    return aucs


def run_gffcompare(gtf: Path, truth: Path, out_prefix: Path):
    subprocess.run(
        ["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}",
         "-v", "/SSD:/SSD", ps.GFFCOMPARE_IMAGE, "gffcompare",
         "-r", str(truth.resolve()), "-o", str(out_prefix), str(gtf)],
        cwd=ps.REPO, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        check=True,
    )


def frontier_point(run_dir: Path, drop_ids: set, truth_paths, workdir: Path, tag: str):
    """Filter the FINAL assembly.gtf by dropped IDs and rescore for real."""
    workdir.mkdir(parents=True, exist_ok=True)
    sim = (workdir / tag).resolve()
    sim.mkdir(exist_ok=True)
    kept_lines = []
    for line in open(run_dir / "assembly.gtf"):
        m = re.search(r'transcript_id "([^"]+)"', line)
        if m and m.group(1) in drop_ids:
            continue
        kept_lines.append(line)
    (sim / "assembly.gtf").write_text("".join(kept_lines))
    run_gffcompare(sim / "assembly.gtf", truth_paths["truth"], sim / "gc")
    metrics = ps.parse_metrics(sim, truth_paths["by_threshold"], 0.0)
    return metrics


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tune-dir", required=True, type=Path)
    ap.add_argument("--tune-sample", default="gencode_p00")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--l2", type=float, default=1.0)
    ap.add_argument("--grid", default="0,250,500,750,1000,1250,1500,1750,2000,2500,3000")
    args = ap.parse_args()

    ds = ps.real_dataset(args.tune_sample)
    truth_paths = {
        "truth": Path(ds.truth),
        "by_threshold": {
            1: ps.expressed_truth(ds.nanocount, threshold=1.0),
            3: ps.expressed_truth(ds.nanocount, threshold=3.0),
        },
    }

    X, y, ids, chroms = load_dataset(
        args.tune_dir / "candidate_evidence.tsv",
        args.tune_dir / "gc.assembly.gtf.tmap",
    )
    print(f"tuning rows={len(y)} exact={int(y.sum())} non-exact={int((1-y).sum())}")

    cv = grouped_cv_auc(X, y, chroms, l2=args.l2)
    print("grouped CV AUC:", [round(a, 4) for a in cv],
          "mean", round(float(np.nanmean(cv)), 4))

    Xs, mean, std = standardize(X)
    w, b = fit_logistic(Xs, y, l2=args.l2)
    scores = Xs @ w + b
    print("full-fit train AUC:", round(auc(scores, y), 4))

    # Frontier: drop the k lowest-scoring candidates, rescore with the real
    # contract. All rows are eligible (p00 = all novel).
    order = np.argsort(scores)
    baseline = ps.parse_metrics(args.tune_dir, truth_paths["by_threshold"], 0.0)
    t1_floor = baseline["honest_f1_t1"]
    print(f"baseline T1={baseline['honest_f1_t1']:.4f} T3={baseline['honest_f1_t3']:.4f}"
          f" (T1 floor for operating point)")
    rows = []
    workdir = args.out / "frontier"
    for k in [int(v) for v in args.grid.split(",")]:
        drop = {ids[i] for i in order[:k]}
        m = frontier_point(args.tune_dir, drop, truth_paths, workdir, f"drop{k}")
        rows.append({
            "drop": k,
            "t1_f1": m["honest_f1_t1"], "t3_f1": m["honest_f1_t3"],
            "t1_ok": m["honest_f1_t1"] >= t1_floor - 1e-9,
        })
        print(f"  drop={k:5d} T1={m['honest_f1_t1']:.4f} T3={m['honest_f1_t3']:.4f}"
              f" {'OK' if rows[-1]['t1_ok'] else 'T1_VIOLATION'}")

    feasible = [r for r in rows if r["t1_ok"]]
    best = max(feasible, key=lambda r: (r["t3_f1"], r["t1_f1"])) if feasible else rows[0]
    k = best["drop"]
    threshold = float(np.sort(scores)[k]) if 0 < k < len(scores) else float("-inf")
    print(f"chosen operating point: drop={k} threshold={threshold:.6f} "
          f"T1={best['t1_f1']:.4f} T3={best['t3_f1']:.4f}")

    model = {
        "schema_version": 1,
        "features": list(FEATURES),
        "mean": mean.tolist(),
        "std": std.tolist(),
        "weights": w.tolist(),
        "bias": float(b),
        "score_threshold": threshold,
        "chosen_drop_k": k,
        "train_sample": args.tune_sample,
        "train_rows": int(len(y)),
        "grouped_cv_auc": [float(a) for a in cv],
        "frontier": rows,
        "baseline_t1": baseline["honest_f1_t1"],
        "baseline_t3": baseline["honest_f1_t3"],
    }
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "ranker_model.json").write_text(json.dumps(model, indent=2) + "\n")
    print("wrote", args.out / "ranker_model.json")


if __name__ == "__main__":
    main()
