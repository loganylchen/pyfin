"""Quickly rebuild the class-partitioned figure from an existing .npz.

Reuses the heavy DTW matrices that viz_class_partitioned.py already computed
and saved to {out}.npz; redoes only the cheap candidate discovery + classing
+ GTF evaluation, then renders the figure with the latest plotting code.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fin.candidates.discovery import discover_candidates
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_fasta import FASTAReader
from fin.io.io_gtf import GTFReader
from fin.scoring.diff_region_dtw import (
    cluster_candidates_by_chain,
    extract_diff_regions,
)

# Reuse helpers from the main script (they're all pure plotting / lookup).
from viz_class_partitioned import (
    _cluster_order,
    _plot_assignment,
    _plot_class_structures,
    _plot_rr,
    evaluate_candidates_vs_gtf,
    remap_novel_ids,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger(__name__)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf-eval", required=True)
    ap.add_argument("--eventalign-root", required=True)
    ap.add_argument("--interval", required=True)
    ap.add_argument("--npz", required=True,
                    help="Existing .npz with M1/M2/M3/M4_via_m1/M4_via_m2.")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chain-wobble", type=int, default=6)
    ap.add_argument("--gtf-eval-tol", type=int, default=6)
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    log.info("Loading GTF + genome ...")
    gtf_eval = GTFReader(args.gtf_eval); gtf_eval.open(); gtf_eval.parse()
    fasta = FASTAReader(args.genome)
    genome = {rec.id: rec.sequence for rec in fasta.get_records()}

    log.info("Generating intervals (no-GTF) ...")
    res = generate_isolated_intervals(args.bam, gtf_path=None,
                                      max_gap=200, max_reads=None)
    want = args.interval.replace(":", "_").replace("-", "_")
    target = None
    for iv in res["intervals"]:
        if iv.region_string.replace(":", "_").replace("-", "_") == want:
            target = iv; break
    assert target is not None, f"interval {args.interval} not found"

    log.info("Discovering candidates (no-GTF) ...")
    cs = discover_candidates(
        interval=target, bam_path=args.bam, gtf_reader=None,
        genome_fasta=genome.get(target.chrom, ""),
        threshold=24, min_novel_reads=1,
    )
    cands = list(cs.candidates)
    remap_novel_ids(cands, Path(args.eventalign_root) / args.interval)

    eval_rows = evaluate_candidates_vs_gtf(
        cands, gtf_eval, target.chrom, target.start, target.end,
        tol=args.gtf_eval_tol,
    )
    n_exact = sum(1 for r in eval_rows if r["match"] == "exact")
    n_wobble = sum(1 for r in eval_rows if r["match"].startswith("wobble"))
    n_novel = sum(1 for r in eval_rows if r["match"] == "novel")

    classes = cluster_candidates_by_chain(cands, chain_wobble=args.chain_wobble)
    diff_regs = extract_diff_regions(cands)

    log.info("Loading matrices from %s ...", args.npz)
    d = np.load(args.npz, allow_pickle=True)
    M3 = d["M3"]; M4_via_m1 = d["M4_via_m1"]; M4_via_m2 = d["M4_via_m2"]
    assign_m1 = d["assign_m1"].tolist()
    assign_m2 = d["assign_m2"].tolist()
    read_ids = d["read_ids"].tolist()
    n = len(read_ids)
    off = ~np.eye(n, dtype=bool)
    m3_nan_frac = float(np.isnan(M3[off]).mean())

    # Reorder npz cand order is preserved; remap candidate_id back so the
    # plot labels match what's in the saved arrays.
    saved_cand_ids = d["cand_ids"].tolist()
    by_id = {c.candidate_id: c for c in cands}
    cands = [by_id[cid] for cid in saved_cand_ids if cid in by_id]
    # Rebuild classes against this re-ordered list.
    classes = cluster_candidates_by_chain(cands, chain_wobble=args.chain_wobble)

    # Clustering reorder driven by M4_via_m2.
    order = _cluster_order(M4_via_m2)
    if order.size == n and not np.array_equal(order, np.arange(n)):
        M3_o = M3[np.ix_(order, order)]
        M4_via_m1_o = M4_via_m1[np.ix_(order, order)]
        M4_via_m2_o = M4_via_m2[np.ix_(order, order)]
        assign_m1_o = [assign_m1[i] for i in order]
        assign_m2_o = [assign_m2[i] for i in order]
    else:
        M3_o, M4_via_m1_o, M4_via_m2_o = M3, M4_via_m1, M4_via_m2
        assign_m1_o, assign_m2_o = assign_m1, assign_m2

    log.info("Plotting ...")
    struct_h = 0.4 + 0.22 * len(cands)
    fig = plt.figure(figsize=(22, 6 + struct_h))
    gs = fig.add_gridspec(
        nrows=4, ncols=3,
        height_ratios=[0.18, struct_h, 2.0, 0.4],
        hspace=0.35, wspace=0.25,
    )

    ax_info = fig.add_subplot(gs[0, :]); ax_info.axis("off")
    eval_by_id = {r["candidate_id"]: r for r in eval_rows}
    info_lines = [
        f"{target.region_string}  |  reads={n}  cands={len(cands)}  "
        f"classes={len(classes)} (wobble=±{args.chain_wobble}bp)  "
        f"diff_regions={len(diff_regs)}",
        f"GTF eval (tol=±{args.gtf_eval_tol}bp):  "
        f"exact={n_exact}  wobble={n_wobble}  novel={n_novel}  "
        f"total={len(eval_rows)}",
    ]
    for ci, members in enumerate(classes):
        labels = []
        for m in members:
            cid = cands[m].candidate_id
            ev = eval_by_id.get(cid, {})
            tag = ev.get("match", "?")
            best = ev.get("best_gtf_tx") or "-"
            labels.append(f"{cid[:10]}[{tag}->{best}]")
        info_lines.append(f"  class {ci}: n={len(members)} | "
                          + ", ".join(labels))
    ax_info.text(0.0, 0.95, "\n".join(info_lines),
                 ha="left", va="top", family="monospace", fontsize=9,
                 transform=ax_info.transAxes)

    ax_struct = fig.add_subplot(gs[1, :])
    _plot_class_structures(ax_struct, cands, classes, eval_rows, target,
                           diff_regions=diff_regs)

    ax_a = fig.add_subplot(gs[2, 0])
    _plot_rr(ax_a, M3_o,
             f"M3 diff-region DTW (off-diag NaN={m3_nan_frac:.0%})")
    ax_b = fig.add_subplot(gs[2, 1])
    _plot_rr(ax_b, M4_via_m1_o,
             "class-partitioned M4 (assign via M1 mappy)")
    ax_c = fig.add_subplot(gs[2, 2])
    _plot_rr(ax_c, M4_via_m2_o,
             "class-partitioned M4 (assign via M2 eventalign)")

    ax_sb = fig.add_subplot(gs[3, 1])
    _plot_assignment(ax_sb, assign_m1_o, len(classes), "read class via M1")
    ax_sc = fig.add_subplot(gs[3, 2])
    _plot_assignment(ax_sc, assign_m2_o, len(classes), "read class via M2")
    fig.add_subplot(gs[3, 0]).axis("off")

    fig.suptitle(
        f"Class-partitioned M4 vs diff-region M3  —  {target.region_string}",
        fontsize=12,
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    log.info("Wrote %s", out_path)


if __name__ == "__main__":
    main()
