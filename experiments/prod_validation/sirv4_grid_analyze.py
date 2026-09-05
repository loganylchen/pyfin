#!/usr/bin/env python3
"""Summarize a NanoRNATrans-style gffcompare benchmark grid.

The scanner tolerates an incomplete result tree. Each PyFIN cell is compared
with the best competitor currently available in that exact (sample,
condition) cell; the report also states how many expected competitor cells are
still missing, so provisional standings cannot be mistaken for a closed
benchmark.

Example:
    python experiments/prod_validation/sirv4_grid_analyze.py \
      --pyfin baseline=.../pyfin_base.tsv \
      --pyfin complete=.../pyfin_catalog.tsv \
      --pyfin adaptive=.../pyfin_adaptive.tsv \
      --out .../report
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from statistics import mean

import numpy as np

DEFAULT_ROOT = Path(
    "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/"
    "results/sirv4_full_sweep_ratio_grid_v2"
)
DEFAULT_OUT = Path(
    "/SSD/logan/dev/pyfin/experiments/prod_validation/"
    "_sirv4grid/report"
)

METRIC_COLS = (
    "tx_sn", "tx_pr", "tx_f1", "ic_sn", "ic_pr", "ic_f1",
    "loc_sn", "loc_pr", "loc_f1",
)
COMP_COLS = (
    "sample", "cond", "tool", "status", "n_query", "n_ref", "n_match",
    *METRIC_COLS,
)
BOOTSTRAP_REPS = 10_000
BOOTSTRAP_SEED = 20260831


def bootstrap_mean_ci(values, *, reps: int = BOOTSTRAP_REPS) -> tuple[float, float]:
    """Deterministic cell-level bootstrap CI for a mean.

    Cells, not transcript rows, are the resampling unit. The same seed is used
    for paired mode differences, preserving cell pairing.
    """
    array = np.asarray(values, dtype=float)
    if array.size == 0:
        return float("nan"), float("nan")
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    indices = rng.integers(0, array.size, size=(reps, array.size))
    boot_means = array[indices].mean(axis=1)
    low, high = np.quantile(boot_means, (0.025, 0.975))
    return float(low), float(high)


def parse_gffcompare_stats(path: Path) -> dict:
    """Parse the stable summary fields emitted by gffcompare 0.12.x."""
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
        (r"Reference mRNAs\s*:\s*(\d+)", "n_ref"),
        (r"Matching transcripts:\s*(\d+)", "n_match"),
    ):
        match = re.search(pattern, text)
        if match:
            out[key] = int(match.group(1))
    return out


def scan_competitors(root: Path) -> list[dict]:
    base = root / "benchmark" / "gffcompare"
    rows = []
    for stats in sorted(base.glob("*/*/*/cmp.stats")):
        rel = stats.relative_to(base)
        sample, cond, tool = rel.parts[:3]
        status_path = stats.parent / "status.tsv"
        status = "ok"
        if status_path.exists():
            status_rows = list(csv.DictReader(status_path.open(), delimiter="\t"))
            status = status_rows[0].get("status", "unknown") if status_rows else "empty"
        metrics = parse_gffcompare_stats(stats) if status == "ok" else {}
        row = {"sample": sample, "cond": cond, "tool": tool,
               "status": "ok" if metrics else status}
        row.update(metrics)
        rows.append(row)
    return rows


def load_pyfin(spec: str) -> tuple[str, dict]:
    if "=" not in spec:
        raise ValueError("--pyfin must be LABEL=PATH")
    label, raw_path = spec.split("=", 1)
    rows = {}
    with Path(raw_path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("status") != "ok":
                continue
            try:
                for col in METRIC_COLS:
                    if row.get(col, ""):
                        row[col] = float(row[col])
                for col in ("n_query", "n_match", "n_anno"):
                    if row.get(col, ""):
                        row[col] = int(float(row[col]))
            except ValueError:
                continue
            rows[(row["sample"], row["cond"])] = row
    return label, rows


def write_tsv(path: Path, rows: list[dict], fields) -> None:
    with path.open("w") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def summarize_mode(
    label: str, pyfin: dict, best: dict, bambu_cells: set,
) -> dict:
    cells = sorted(set(pyfin) & set(best))
    wins = sum(pyfin[k]["tx_f1"] > best[k]["tx_f1"] + 1e-9 for k in cells)
    ties = sum(abs(pyfin[k]["tx_f1"] - best[k]["tx_f1"]) <= 1e-9
               for k in cells)
    gaps = [best[k]["tx_f1"] - pyfin[k]["tx_f1"] for k in cells
            if pyfin[k]["tx_f1"] < best[k]["tx_f1"] - 1e-9]
    f1_values = [pyfin[k]["tx_f1"] for k in cells]
    ci_low, ci_high = bootstrap_mean_ci(f1_values)
    bambu_present = [k for k in cells if k in bambu_cells]
    return {
        "mode": label,
        "paired_cells": len(cells),
        "wins": wins,
        "ties": ties,
        "top_or_tied": wins + ties,
        "behind": len(cells) - wins - ties,
        "bambu_present_cells": len(bambu_present),
        "top_or_tied_bambu_present": sum(
            pyfin[k]["tx_f1"] >= best[k]["tx_f1"] - 1e-9
            for k in bambu_present
        ),
        "mean_tx_sn": mean(pyfin[k]["tx_sn"] for k in cells),
        "mean_tx_pr": mean(pyfin[k]["tx_pr"] for k in cells),
        "mean_tx_f1": mean(f1_values),
        "mean_tx_f1_ci_low": ci_low,
        "mean_tx_f1_ci_high": ci_high,
        "mean_best_competitor_f1": mean(best[k]["tx_f1"] for k in cells),
        "mean_gap_when_behind": mean(gaps) if gaps else 0.0,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--pyfin", action="append", required=True,
                        help="repeatable LABEL=PATH PyFIN result table")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    competitors = scan_competitors(args.root)
    ok_comp = [r for r in competitors if r["status"] == "ok"]
    write_tsv(args.out / "competitor_cells.tsv", ok_comp, COMP_COLS)

    by_cell = defaultdict(list)
    for row in ok_comp:
        by_cell[(row["sample"], row["cond"])].append(row)
    best = {
        key: max(rows, key=lambda row: (row["tx_f1"], row["tool"]))
        for key, rows in by_cell.items()
    }
    expected_tools = sorted({row["tool"] for row in ok_comp})
    samples = sorted({key[0] for key in by_cell})
    conds = sorted({key[1] for key in by_cell})
    expected = {(s, c, t) for s in samples for c in conds
                for t in expected_tools}
    observed = {(r["sample"], r["cond"], r["tool"]) for r in ok_comp}
    missing = [dict(zip(("sample", "cond", "tool"), key))
               for key in sorted(expected - observed)]
    write_tsv(args.out / "missing_competitor_cells.tsv", missing,
              ("sample", "cond", "tool"))

    bambu_cells = {
        (row["sample"], row["cond"])
        for row in ok_comp if row["tool"] == "bambu"
    }
    modes = dict(load_pyfin(spec) for spec in args.pyfin)
    summaries = [summarize_mode(label, rows, best, bambu_cells)
                 for label, rows in modes.items()]
    summary_fields = (
        "mode", "paired_cells", "wins", "ties", "top_or_tied", "behind",
        "bambu_present_cells", "top_or_tied_bambu_present",
        "mean_tx_sn", "mean_tx_pr", "mean_tx_f1", "mean_tx_f1_ci_low",
        "mean_tx_f1_ci_high", "mean_best_competitor_f1",
        "mean_gap_when_behind",
    )
    write_tsv(args.out / "mode_summary.tsv", summaries, summary_fields)

    difference_rows = []
    if "adaptive" in modes and "complete" in modes:
        cells = sorted(set(modes["adaptive"]) & set(modes["complete"]))
        differences = [
            modes["adaptive"][key]["tx_f1"]
            - modes["complete"][key]["tx_f1"]
            for key in cells
        ]
        ci_low, ci_high = bootstrap_mean_ci(differences)
        difference_rows.append({
            "contrast": "adaptive-complete",
            "paired_cells": len(cells),
            "mean_f1_difference": mean(differences),
            "ci_low": ci_low,
            "ci_high": ci_high,
            "bootstrap_reps": BOOTSTRAP_REPS,
            "bootstrap_seed": BOOTSTRAP_SEED,
        })
    difference_fields = (
        "contrast", "paired_cells", "mean_f1_difference", "ci_low",
        "ci_high", "bootstrap_reps", "bootstrap_seed",
    )
    write_tsv(args.out / "paired_differences.tsv", difference_rows,
              difference_fields)

    condition_rows = []
    for label, rows in modes.items():
        for cond in conds:
            cells = sorted(k for k in set(rows) & set(best) if k[1] == cond)
            if not cells:
                continue
            condition_rows.append({
                "mode": label,
                "cond": cond,
                "paired_cells": len(cells),
                "top_or_tied": sum(
                    rows[k]["tx_f1"] >= best[k]["tx_f1"] - 1e-9
                    for k in cells
                ),
                "mean_tx_sn": mean(rows[k]["tx_sn"] for k in cells),
                "mean_tx_pr": mean(rows[k]["tx_pr"] for k in cells),
                "mean_tx_f1": mean(rows[k]["tx_f1"] for k in cells),
                "mean_best_competitor_f1": mean(
                    best[k]["tx_f1"] for k in cells
                ),
            })
    condition_fields = (
        "mode", "cond", "paired_cells", "top_or_tied", "mean_tx_sn",
        "mean_tx_pr", "mean_tx_f1", "mean_best_competitor_f1",
    )
    write_tsv(args.out / "condition_summary.tsv", condition_rows,
              condition_fields)

    lines = [
        "# Grid summary", "",
        f"- Samples: {len(samples)}",
        f"- Conditions: {len(conds)}",
        f"- Competitor tools observed: {len(expected_tools)} "
        f"({', '.join(expected_tools)})",
        f"- Valid competitor cells: {len(ok_comp)}/{len(expected)}",
        f"- Missing competitor cells: {len(missing)}",
        f"- Bambu present: {len(bambu_cells)}/{len(samples) * len(conds)} "
        f"cells ({len(samples) * len(conds) - len(bambu_cells)} missing)",
        f"- Mean-F1 CIs: {BOOTSTRAP_REPS:,} cell-level bootstrap resamples, "
        f"seed {BOOTSTRAP_SEED}",
        "- Standings are provisional wherever competitor cells remain missing; "
        "the Bambu-present column quantifies current exposure to its missing cells.",
        "", "| Mode | Paired | Top/tied all | Top/tied with Bambu present | "
        "Mean Sn | Mean Pr | Mean F1 (95% CI) | Best competitor |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summaries:
        lines.append(
            f"| {row['mode']} | {row['paired_cells']} | "
            f"{row['top_or_tied']} | "
            f"{row['top_or_tied_bambu_present']}/"
            f"{row['bambu_present_cells']} | "
            f"{row['mean_tx_sn']:.2f} | {row['mean_tx_pr']:.2f} | "
            f"{row['mean_tx_f1']:.2f} "
            f"({row['mean_tx_f1_ci_low']:.2f}, "
            f"{row['mean_tx_f1_ci_high']:.2f}) | "
            f"{row['mean_best_competitor_f1']:.2f} |"
        )
    if difference_rows:
        difference = difference_rows[0]
        lines.extend([
            "",
            f"Adaptive - complete paired mean F1: "
            f"{difference['mean_f1_difference']:+.3f} "
            f"(95% bootstrap CI {difference['ci_low']:+.3f}, "
            f"{difference['ci_high']:+.3f}).",
        ])
    (args.out / "SUMMARY.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {args.out}: {len(ok_comp)} competitor cells, "
          f"{len(missing)} missing, {len(modes)} PyFIN modes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
