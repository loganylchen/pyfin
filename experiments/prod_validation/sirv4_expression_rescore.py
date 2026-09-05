#!/usr/bin/env python3
"""Expression-conditioned re-score of the SIRV4 annotation grid.

The official grid runs gffcompare against the full 176-transcript SIRV4 GTF.
That measures catalog/structure recovery: a catalog transcript is rewarded even
when its sample-specific NanoCount TPM is zero. This script evaluates the
separate question "which transcripts are expressed in this sample?" by
subsetting the same structural truth to NanoCount TPM > 0 for each sample.

Important: NanoCount was estimated from the same reads, not from an independent
mixture-concentration truth. Results are an expression-conditioned pseudo-truth
and must be reported alongside, never in place of, the official catalog score.
The pseudo-truth is used only for scoring; it is never passed to PyFIN.
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from statistics import mean

import numpy as np

ROOT = Path(
    "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/"
    "results/sirv4_full_sweep_ratio_grid_v2"
)
WORK = Path(
    "/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid"
)
MASTER_TRUTH = WORK / "_ref" / "SIRV_Set4.gtf"
IMAGE = "quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"
REPS = 10_000
SEED = 20260831

COLS = (
    "sample", "cond", "tool", "family", "status", "truth_transcripts",
    "n_query", "n_match", "tx_sn", "tx_pr", "tx_f1", "ic_sn", "ic_pr",
    "ic_f1", "loc_sn", "loc_pr", "loc_f1",
)


def transcript_ids(path: Path) -> set[str]:
    ids = set()
    pattern = re.compile(r'transcript_id "([^"]+)"')
    with path.open() as handle:
        for line in handle:
            match = pattern.search(line)
            if match:
                ids.add(match.group(1))
    return ids


def expressed_ids(sample: str) -> set[str]:
    path = ROOT / "_ref" / sample / "expression" / "nanocount.tsv"
    ids = set()
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if float(row["tpm"]) > 0.0:
                ids.add(row["transcript_name"])
    return ids


def build_truth(sample: str, out_dir: Path) -> tuple[Path, int]:
    keep = expressed_ids(sample)
    master_ids = transcript_ids(MASTER_TRUTH)
    missing = keep - master_ids
    if missing:
        raise ValueError(
            f"{sample}: {len(missing)} expressed IDs absent from master truth"
        )
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{sample}.gtf"
    pattern = re.compile(r'transcript_id "([^"]+)"')
    with MASTER_TRUTH.open() as source, path.open("w") as target:
        for line in source:
            match = pattern.search(line)
            if match and match.group(1) in keep:
                target.write(line)
    written = transcript_ids(path)
    if written != keep:
        raise ValueError(
            f"{sample}: expressed truth mismatch, expected {len(keep)}, "
            f"wrote {len(written)}"
        )
    return path, len(keep)


def parse_stats(path: Path) -> dict:
    if not path.exists():
        return {}
    text = path.read_text()
    out = {}
    for label, key in (
        ("Transcript level", "tx"),
        ("Intron chain level", "ic"),
        ("Locus level", "loc"),
    ):
        match = re.search(
            rf"{re.escape(label)}:\s*([\d.]+)\s*\|\s*([\d.]+)", text
        )
        if match:
            sensitivity, precision = map(float, match.groups())
            out[f"{key}_sn"] = sensitivity
            out[f"{key}_pr"] = precision
            out[f"{key}_f1"] = (
                0.0 if sensitivity + precision == 0
                else 2 * sensitivity * precision / (sensitivity + precision)
            )
    for pattern, key in (
        (r"Query mRNAs\s*:\s*(\d+)", "n_query"),
        (r"Matching transcripts:\s*(\d+)", "n_match"),
    ):
        match = re.search(pattern, text)
        if match:
            out[key] = int(match.group(1))
    return out


def score_job(job) -> dict:
    sample, cond, tool, family, query, truth, n_truth, out_root, force = job
    cell = Path(out_root) / sample / cond / tool
    stats_path = cell / "cmp.stats"
    if force and cell.exists():
        shutil.rmtree(cell)
    cell.mkdir(parents=True, exist_ok=True)
    local_query = cell / "query.gtf"
    if not stats_path.exists():
        shutil.copyfile(query, local_query)
        command = [
            "docker", "run", "--rm", "-v", "/SSD:/SSD", IMAGE,
            "gffcompare", "-r", str(truth), "-o", str(cell / "cmp"),
            str(local_query),
        ]
        result = subprocess.run(command, capture_output=True)
        if result.returncode != 0:
            return {
                "sample": sample, "cond": cond, "tool": tool,
                "family": family, "status": f"gffcompare_rc{result.returncode}",
                "truth_transcripts": n_truth,
            }
    metrics = parse_stats(stats_path)
    row = {
        "sample": sample, "cond": cond, "tool": tool, "family": family,
        "status": "ok" if metrics else "missing_stats",
        "truth_transcripts": n_truth,
    }
    row.update(metrics)
    return row


def bootstrap_ci(values) -> tuple[float, float]:
    array = np.asarray(values, dtype=float)
    rng = np.random.default_rng(SEED)
    indices = rng.integers(0, len(array), size=(REPS, len(array)))
    low, high = np.quantile(array[indices].mean(axis=1), (0.025, 0.975))
    return float(low), float(high)


def write_tsv(path: Path, rows: list[dict], fields) -> None:
    with path.open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict], out_dir: Path) -> None:
    valid = [row for row in rows if row["status"] == "ok"]
    competitors = defaultdict(dict)
    modes = defaultdict(dict)
    for row in valid:
        key = (row["sample"], row["cond"])
        if row["family"] == "competitor":
            competitors[key][row["tool"]] = row
        else:
            modes[row["tool"]][key] = row
    best = {
        key: max(tool_rows.values(), key=lambda row: (row["tx_f1"], row["tool"]))
        for key, tool_rows in competitors.items()
    }
    bambu_cells = {key for key, tools in competitors.items() if "bambu" in tools}

    mode_rows = []
    for mode, data in sorted(modes.items()):
        cells = sorted(set(data) & set(best))
        f1 = [data[key]["tx_f1"] for key in cells]
        low, high = bootstrap_ci(f1)
        top = sum(data[key]["tx_f1"] >= best[key]["tx_f1"] - 1e-9
                  for key in cells)
        b_cells = [key for key in cells if key in bambu_cells]
        mode_rows.append({
            "mode": mode,
            "paired_cells": len(cells),
            "top_or_tied": top,
            "bambu_present_cells": len(b_cells),
            "top_or_tied_bambu_present": sum(
                data[key]["tx_f1"] >= best[key]["tx_f1"] - 1e-9
                for key in b_cells
            ),
            "mean_truth_transcripts": mean(
                data[key]["truth_transcripts"] for key in cells
            ),
            "mean_tx_sn": mean(data[key]["tx_sn"] for key in cells),
            "mean_tx_pr": mean(data[key]["tx_pr"] for key in cells),
            "mean_tx_f1": mean(f1),
            "ci_low": low,
            "ci_high": high,
            "mean_best_competitor_f1": mean(best[key]["tx_f1"] for key in cells),
        })
    mode_fields = (
        "mode", "paired_cells", "top_or_tied", "bambu_present_cells",
        "top_or_tied_bambu_present", "mean_truth_transcripts", "mean_tx_sn",
        "mean_tx_pr", "mean_tx_f1", "ci_low", "ci_high",
        "mean_best_competitor_f1",
    )
    write_tsv(out_dir / "mode_summary.tsv", mode_rows, mode_fields)

    conditions = sorted({key[1] for key in best})
    condition_rows = []
    for mode, data in sorted(modes.items()):
        for cond in conditions:
            cells = sorted(
                key for key in set(data) & set(best) if key[1] == cond
            )
            if not cells:
                continue
            f1 = [data[key]["tx_f1"] for key in cells]
            delta = [data[key]["tx_f1"] - best[key]["tx_f1"] for key in cells]
            low, high = bootstrap_ci(f1)
            delta_low, delta_high = bootstrap_ci(delta)
            condition_rows.append({
                "mode": mode,
                "cond": cond,
                "samples": len(cells),
                "mean_truth_transcripts": mean(
                    data[key]["truth_transcripts"] for key in cells
                ),
                "mean_tx_sn": mean(data[key]["tx_sn"] for key in cells),
                "mean_tx_pr": mean(data[key]["tx_pr"] for key in cells),
                "mean_tx_f1": mean(f1),
                "ci_low": low,
                "ci_high": high,
                "mean_best_competitor_f1": mean(
                    best[key]["tx_f1"] for key in cells
                ),
                "mean_delta": mean(delta),
                "delta_ci_low": delta_low,
                "delta_ci_high": delta_high,
                "top_or_tied": sum(
                    data[key]["tx_f1"] >= best[key]["tx_f1"] - 1e-9
                    for key in cells
                ),
                "bambu_present": sum(key in bambu_cells for key in cells),
            })
    condition_fields = (
        "mode", "cond", "samples", "mean_truth_transcripts", "mean_tx_sn",
        "mean_tx_pr", "mean_tx_f1", "ci_low", "ci_high",
        "mean_best_competitor_f1", "mean_delta", "delta_ci_low",
        "delta_ci_high", "top_or_tied", "bambu_present",
    )
    write_tsv(out_dir / "condition_summary.tsv", condition_rows,
              condition_fields)

    tool_wins = Counter()
    for tools in competitors.values():
        maximum = max(row["tx_f1"] for row in tools.values())
        for tool, row in tools.items():
            if row["tx_f1"] >= maximum - 1e-9:
                tool_wins[tool] += 1
    lines = [
        "# Expression-conditioned SIRV4 score", "",
        "NanoCount TPM > 0 is a same-read expression pseudo-truth, not an "
        "independent concentration truth.", "",
        f"- Valid score rows: {len(valid)}",
        f"- Competitor cells: {sum(len(v) for v in competitors.values())}",
        f"- Bambu present: {len(bambu_cells)}/156", "",
        "| Mode | Top/tied | Top/tied (Bambu present) | Mean truth tx | "
        "Mean Sn | Mean Pr | Mean F1 (95% CI) | Best competitor |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in mode_rows:
        lines.append(
            f"| {row['mode']} | {row['top_or_tied']}/"
            f"{row['paired_cells']} | {row['top_or_tied_bambu_present']}/"
            f"{row['bambu_present_cells']} | "
            f"{row['mean_truth_transcripts']:.1f} | "
            f"{row['mean_tx_sn']:.2f} | {row['mean_tx_pr']:.2f} | "
            f"{row['mean_tx_f1']:.2f} ({row['ci_low']:.2f}, "
            f"{row['ci_high']:.2f}) | "
            f"{row['mean_best_competitor_f1']:.2f} |"
        )
    lines.extend(["", "Competitor best-or-tied counts: " + ", ".join(
        f"{tool}={count}" for tool, count in tool_wins.most_common()
    )])
    (out_dir / "SUMMARY.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--out", type=Path,
                        default=WORK / "expression_score")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    samples = sorted(path.name for path in ROOT.glob("SGNex_*") if path.is_dir())
    truth_dir = WORK / "_ref" / "expressed"
    truths = {sample: build_truth(sample, truth_dir) for sample in samples}
    print("expression truths:", ", ".join(
        f"{sample[-15:]}={truths[sample][1]}" for sample in samples
    ), flush=True)

    jobs = []
    competitor_table = WORK / "report" / "competitor_cells.tsv"
    for row in csv.DictReader(competitor_table.open(), delimiter="\t"):
        if row["status"] != "ok":
            continue
        sample, cond, tool = row["sample"], row["cond"], row["tool"]
        query = ROOT / sample / cond / "assembly" / f"{tool}.gtf"
        truth, n_truth = truths[sample]
        jobs.append((sample, cond, tool, "competitor", query, truth,
                     n_truth, args.out / "cells", args.force))

    pyfin_modes = {
        "pyfin_off": WORK / "base",
        "pyfin_complete": WORK / "catalog",
        "pyfin_adaptive": WORK / "adaptive",
        # Existing 6-sample x 8-condition ablation: preserve ordinary novel
        # filtering, but do not raise the GTF abundance floor to the SIRV
        # novel floor. Partial by design; summarized only on its 48 paired cells.
        "pyfin_keepgtf": WORK / "keepgtf",
        "pyfin_keepgtf_full": WORK / "keepgtf_full",
        "pyfin_keepgtf_positive": WORK / "keepgtf_positive",
        "pyfin_guided_suppress": WORK / "guided_suppress",
        # P0 live implementation: --require-gtf-positive-mass with
        # --min-gtf-abundance 0 --no-floor-gtf-abundance (strict >0 GTF gate
        # inside selection, not an offline row filter).
        "pyfin_p0": WORK / "p0_strict",
        # P2 complete-chain rescue on top of P0: --complete-chain-rescue-reads 2
        # (weakest-junction clean attestation by exact complete-chain reads
        # exempts a candidate from the support-gate and soft-mass drops).
        "pyfin_p2": WORK / "p2_rescue",
        # Containment decision layer validation modes (advisor protocol):
        # audit on top of P0 (must equal pyfin_p0 scores exactly), and the
        # full P0 + rescue + containment-require mode.
        "pyfin_p2caudit": WORK / "p2_caudit",
        "pyfin_p2require": WORK / "p2_require",
        # Step 5a per-read junction correction pilot (5 conditions x 6 samples):
        # audit must equal jc_off exactly; require is the behavioural arm.
        "pyfin_jc_off": WORK / "jc_off",
        "pyfin_jc_audit": WORK / "jc_audit",
        "pyfin_jc_require": WORK / "jc_req",
        "pyfin_jc2_off": WORK / "jc2_off",
        "pyfin_jc2_require": WORK / "jc2_req",
        # P3 signal-ballot junction election on top of P0+rescue.
        "pyfin_je_audit": WORK / "je_full_audit",
        "pyfin_je_require": WORK / "je_full_require",
        # Ballot ablations (screen = nomination alone, full grid; m1 /
        # shuffle = 24-cell subsets) proving the raw signal carries the
        # decision.
        "pyfin_je_screen": WORK / "je_screen",
        "pyfin_je_m1": WORK / "je_m1",
        "pyfin_je_shuffle": WORK / "je_shuffle",
        "pyfin_je_defense": WORK / "je_defense",
        # 2026-09-03 generation refactor (GTF chains as ordinary variant keys,
        # CIGAR-level collapse span_guard) -- NOT byte-inert by design. gv_off =
        # plain defaults (compare pyfin_off); gv_best = P0+rescue+defense
        # (compare pyfin_je_defense); tss_audit must equal gv_best exactly;
        # tss_req = + --folded-tss-recovery-mode require (behavioural arm).
        # gv_off_foldgtf = the rejected Codex rule "pure-GTF exact sub-chain folds
        # unconditionally" (kept as the record of WHY annotated keys are never foldees:
        # 82.55 F1 vs pyfin_off 85.58, Sn -3.8, every guided condition down).
        "pyfin_gv_off_foldgtf": WORK / "gv_off_foldgtf",
        # gv_off_clusterorder = annotated keys never foldees but candidates emitted
        # in cluster order (novel columns before annotation): Sn held (-0.02) but
        # Pr -0.93 because R.argmax ties went to novel shadows -> annotation-first
        # column order restored.
        "pyfin_gv_off_clusterorder": WORK / "gv_off_clusterorder",
        # Current tree (explicit hard tie-break + non-bridging annotation keys +
        # CIGAR span_guard), evaluated on BOTH grids by cross_dataset_eval.sbatch.
        "pyfin_xd_off": WORK / "xd_off",
        "pyfin_xd_best": WORK / "xd_best",
        "pyfin_gv_best": WORK / "gv_best",
        "pyfin_gv_tss_audit": WORK / "gv_tss_audit",
        "pyfin_gv_tss_require": WORK / "gv_tss_req",
    }
    for tool, root in pyfin_modes.items():
        for query in sorted(root.glob("SGNex_*/*/out/assembly.gtf")):
            sample = query.parents[2].name
            cond = query.parents[1].name
            truth, n_truth = truths[sample]
            jobs.append((sample, cond, tool, "pyfin", query, truth,
                         n_truth, args.out / "cells", args.force))

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(score_job, job) for job in jobs]
        for done, future in enumerate(as_completed(futures), 1):
            row = future.result()
            rows.append(row)
            if done % 100 == 0 or row["status"] != "ok":
                print(f"[{done}/{len(jobs)}] {row['status']}", flush=True)
    rows.sort(key=lambda row: (
        row["sample"], row["cond"], row["family"], row["tool"]
    ))
    write_tsv(args.out / "scores.tsv", rows, COLS)
    summarize(rows, args.out)
    ok = sum(row["status"] == "ok" for row in rows)
    print(f"done: {ok}/{len(rows)} scores -> {args.out}", flush=True)
    return 0 if ok == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
