#!/usr/bin/env python
"""Evaluate em_subset_ablation outputs: 14 GTFs × intron-chain P/R/F1.

Reads the summary.tsv from em_subset_ablation.py, opens each cell's
predicted GTF, builds (chrom, strand, intron_chain) tuples for multi-exon
transcripts (plus (chrom, strand, start, end) for single-exon), and
compares against a ground-truth GTF.

Outputs:
    eval_results.tsv   — one row per (subset, mode) with TP/FP/FN/P/R/F1
    eval_results.md    — markdown table grouped by mode

CLI example:
    python benchmarks/em_subset_evaluate.py \\
        --summary       /tmp/em_subset/summary.tsv \\
        --ground-truth  /data/sirv/sirv.gtf \\
        --out-tsv       /tmp/em_subset/eval_results.tsv \\
        --out-md        /tmp/em_subset/eval_results.md
"""
from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Set, Tuple

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("em_subset_evaluate")

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def _parse_attrs(s: str) -> dict[str, str]:
    return {m.group(1): m.group(2) for m in _ATTR_RE.finditer(s)}


@dataclass
class _Tx:
    tid: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]] = field(default_factory=list)

    @property
    def intron_chain(self) -> Tuple[Tuple[int, int], ...]:
        if len(self.exons) < 2:
            return ()
        se = sorted(self.exons)
        return tuple((se[i][1], se[i + 1][0]) for i in range(len(se) - 1))

    @property
    def span(self) -> Tuple[int, int]:
        if not self.exons:
            return (0, 0)
        return (min(e[0] for e in self.exons), max(e[1] for e in self.exons))


def _load_gtf(path: str) -> dict[str, _Tx]:
    """Parse a GTF into {transcript_id -> _Tx} with 0-based half-open exons."""
    txs: dict[str, _Tx] = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feat = parts[2]
            if feat not in ("transcript", "exon"):
                continue
            chrom = parts[0]
            start = int(parts[3]) - 1  # GTF is 1-based inclusive → 0-based half-open
            end = int(parts[4])
            strand = parts[6]
            attrs = _parse_attrs(parts[8])
            tid = attrs.get("transcript_id", "")
            if not tid:
                continue
            tx = txs.get(tid)
            if tx is None:
                tx = _Tx(tid=tid, chrom=chrom, strand=strand)
                txs[tid] = tx
            if feat == "exon":
                tx.exons.append((start, end))
    return txs


def _build_chain_keys(
    txs: dict[str, _Tx],
) -> Tuple[Set[tuple], Set[tuple]]:
    """Return (multi_exon_chain_keys, single_exon_keys).

    multi_exon: (chrom, strand, intron_chain)
    single_exon: (chrom, strand, start, end)
    """
    multi: Set[tuple] = set()
    mono: Set[tuple] = set()
    for tx in txs.values():
        if not tx.exons:
            continue
        ic = tx.intron_chain
        if ic:
            multi.add((tx.chrom, tx.strand, ic))
        else:
            s, e = tx.span
            mono.add((tx.chrom, tx.strand, s, e))
    return multi, mono


@dataclass
class CellResult:
    subset: str
    mode: str
    tp_multi: int = 0
    fp_multi: int = 0
    fn_multi: int = 0
    tp_mono: int = 0
    fp_mono: int = 0
    fn_mono: int = 0
    n_pred_total: int = 0

    @property
    def tp(self) -> int:
        return self.tp_multi + self.tp_mono

    @property
    def fp(self) -> int:
        return self.fp_multi + self.fp_mono

    @property
    def fn(self) -> int:
        return self.fn_multi + self.fn_mono

    @property
    def precision(self) -> float:
        d = self.tp + self.fp
        return self.tp / d if d else 0.0

    @property
    def recall(self) -> float:
        d = self.tp + self.fn
        return self.tp / d if d else 0.0

    @property
    def f1(self) -> float:
        p, r = self.precision, self.recall
        return 2 * p * r / (p + r) if (p + r) > 0 else 0.0


def _evaluate_gtf(
    gtf_path: str,
    gt_multi: Set[tuple],
    gt_mono: Set[tuple],
    subset: str,
    mode: str,
) -> CellResult:
    pred_txs = _load_gtf(gtf_path)
    pred_multi, pred_mono = _build_chain_keys(pred_txs)
    cr = CellResult(subset=subset, mode=mode, n_pred_total=len(pred_txs))
    cr.tp_multi = len(pred_multi & gt_multi)
    cr.fp_multi = len(pred_multi - gt_multi)
    cr.fn_multi = len(gt_multi - pred_multi)
    cr.tp_mono = len(pred_mono & gt_mono)
    cr.fp_mono = len(pred_mono - gt_mono)
    cr.fn_mono = len(gt_mono - pred_mono)
    return cr


def _load_summary(path: str) -> list[dict]:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _write_tsv(out: Path, results: list[CellResult]) -> None:
    fields = [
        "subset", "mode",
        "tp", "fp", "fn",
        "tp_multi", "fp_multi", "fn_multi",
        "tp_mono", "fp_mono", "fn_mono",
        "precision", "recall", "f1",
        "n_predicted",
    ]
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for cr in results:
            w.writerow({
                "subset": cr.subset,
                "mode": cr.mode,
                "tp": cr.tp, "fp": cr.fp, "fn": cr.fn,
                "tp_multi": cr.tp_multi, "fp_multi": cr.fp_multi, "fn_multi": cr.fn_multi,
                "tp_mono": cr.tp_mono, "fp_mono": cr.fp_mono, "fn_mono": cr.fn_mono,
                "precision": f"{cr.precision:.4f}",
                "recall": f"{cr.recall:.4f}",
                "f1": f"{cr.f1:.4f}",
                "n_predicted": cr.n_pred_total,
            })


def _write_md(
    out: Path,
    results: list[CellResult],
    gt_multi: Set[tuple],
    gt_mono: Set[tuple],
) -> None:
    lines: list[str] = []
    lines.append("# EM matrix-subset ablation: intron-chain P/R/F1\n")
    lines.append(
        f"Ground truth: {len(gt_multi)} multi-exon chains + "
        f"{len(gt_mono)} single-exon = {len(gt_multi) + len(gt_mono)} total\n"
    )

    by_mode = defaultdict(list)
    for cr in results:
        by_mode[cr.mode].append(cr)

    for mode in sorted(by_mode):
        lines.append(f"## Mode: {mode}\n")
        lines.append(
            "| subset | TP | FP | FN | Precision | Recall | F1 | n_pred |"
        )
        lines.append(
            "|--------|----|----|----|-----------|--------|----|--------|"
        )
        # Stable subset order
        order = {s: i for i, s in enumerate(
            ("m1", "m2", "m3", "m1+m2", "m1+m3", "m2+m3", "all")
        )}
        for cr in sorted(by_mode[mode], key=lambda r: order.get(r.subset, 99)):
            lines.append(
                f"| {cr.subset} | {cr.tp} | {cr.fp} | {cr.fn} | "
                f"{cr.precision:.3f} | {cr.recall:.3f} | {cr.f1:.3f} | "
                f"{cr.n_pred_total} |"
            )
        lines.append("")

    out.write_text("\n".join(lines))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", required=True,
                    help="summary.tsv from em_subset_ablation.py")
    ap.add_argument("--ground-truth", required=True,
                    help="Reference annotation GTF")
    ap.add_argument("--out-tsv", default=None,
                    help="Per-cell results TSV (default: alongside summary)")
    ap.add_argument("--out-md", default=None,
                    help="Markdown table output (default: alongside summary)")
    args = ap.parse_args(argv)

    summary_path = Path(args.summary)
    out_tsv = Path(args.out_tsv or summary_path.with_name("eval_results.tsv"))
    out_md = Path(args.out_md or summary_path.with_name("eval_results.md"))

    logger.info("Loading ground truth from %s", args.ground_truth)
    gt_txs = _load_gtf(args.ground_truth)
    gt_multi, gt_mono = _build_chain_keys(gt_txs)
    logger.info(
        "GT: %d transcripts → %d multi-exon chains + %d single-exon",
        len(gt_txs), len(gt_multi), len(gt_mono),
    )

    rows = _load_summary(args.summary)
    logger.info("Evaluating %d cells from %s", len(rows), summary_path)

    results: list[CellResult] = []
    for r in rows:
        gtf_path = r.get("gtf_out", "")
        if not gtf_path or not Path(gtf_path).exists():
            logger.warning("Missing GTF for cell %s/%s: %s",
                           r["subset"], r["mode"], gtf_path)
            continue
        cr = _evaluate_gtf(gtf_path, gt_multi, gt_mono,
                           r["subset"], r["mode"])
        results.append(cr)
        logger.info(
            "%s/%s  TP=%d FP=%d FN=%d  P=%.3f R=%.3f F1=%.3f",
            cr.subset, cr.mode, cr.tp, cr.fp, cr.fn,
            cr.precision, cr.recall, cr.f1,
        )

    _write_tsv(out_tsv, results)
    _write_md(out_md, results, gt_multi, gt_mono)
    logger.info("Wrote %s and %s", out_tsv, out_md)

    # Echo markdown to stdout for easy copy-paste
    print(out_md.read_text())


if __name__ == "__main__":
    main()
