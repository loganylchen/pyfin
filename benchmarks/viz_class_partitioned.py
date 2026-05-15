"""Visualize class-partitioned M4 vs diff-region M3 on a SIRV interval.

Discovery is performed WITHOUT GTF (gtf_reader=None) so candidates come from
the BAM alone. The GTF is used only for downstream evaluation: each
discovered candidate is compared to overlapping GTF transcripts, and we
label each candidate as "exact" (intron chain matches a GTF tx within
±tol bp per junction) or "novel".

Usage:
    python benchmarks/viz_class_partitioned.py \
        --bam testdata/mapped.bam \
        --genome testdata/SIRV.genome.fa \
        --gtf-eval testdata/SIRV.genome.gtf \
        --signal testdata/mapped.blow5 \
        --eventalign-root testdata/out_no_gtf \
        --interval SIRVomeERCCome_34508_36900 \
        --out testdata/viz/class_partitioned.png
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from fin.candidates.discovery import discover_candidates
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_fasta import FASTAReader
from fin.io.io_gtf import GTFReader
from fin.io.io_slow5 import Slow5Reader
from fin.scoring.diff_region_dtw import (
    cluster_candidates_by_chain,
    compute_class_partitioned_m4,
    compute_diff_region_m4,
    extract_diff_regions,
)
from fin.scoring.eventalign_parser import (
    build_distance_matrix,
    parse_eventalign_tsv,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger(__name__)


def compute_m1_mappy(read_seqs: Dict[str, str], candidates):
    import mappy
    read_ids = sorted(read_seqs.keys())
    cand_ids = [c.candidate_id for c in candidates]
    aligners = []
    for c in candidates:
        if not c.sequence:
            aligners.append(None); continue
        a = mappy.Aligner(seq=c.sequence, preset="map-ont")
        aligners.append(a if a else None)
    n_r, n_c = len(read_ids), len(cand_ids)
    M = np.full((n_r, n_c), 1e6, dtype=np.float32)
    for i, rid in enumerate(read_ids):
        seq = read_seqs.get(rid, "")
        if not seq:
            continue
        for j, aln in enumerate(aligners):
            if aln is None:
                continue
            best = None
            for hit in aln.map(seq):
                s = getattr(hit, "score", None)
                if s is None:
                    s = hit.mlen
                if best is None or s > best:
                    best = s
            if best is not None:
                M[i, j] = float(-best)
    return M, read_ids, cand_ids


def compute_m2_eventalign(eventalign_root: Path, interval_name: str, candidates, read_ids):
    cand_ids = [c.candidate_id for c in candidates]
    cand_lengths = {c.candidate_id: len(c.sequence) for c in candidates}
    all_scores = []
    for cid in cand_ids:
        tsv = eventalign_root / interval_name / cid / "eventalign.tsv"
        if not tsv.exists():
            log.warning("Missing eventalign TSV: %s", tsv)
            continue
        all_scores.extend(parse_eventalign_tsv(str(tsv), cand_lengths, collect_events=True))
    M = build_distance_matrix(all_scores, read_ids, cand_ids).astype(np.float32)
    scores_by_pair = {(s.read_name, s.candidate_id): s for s in all_scores}
    return M, scores_by_pair


def remap_novel_ids(cands, ev_root: Path):
    """Mutate cand.candidate_id to match eventalign dir by read-overlap."""
    dir_read_sets: Dict[str, set] = {}
    if not ev_root.exists():
        return
    for d in sorted(ev_root.iterdir()):
        if not d.is_dir():
            continue
        tsv = d / "eventalign.tsv"
        if not tsv.exists():
            continue
        reads = set()
        with open(tsv) as fh:
            fh.readline()
            for line in fh:
                p = line.split("\t", 5)
                if len(p) > 4:
                    reads.add(p[3])
        dir_read_sets[d.name] = reads
    used = set()
    for c in cands:
        if c.source != "novel":
            continue
        best, best_ov = None, 0
        for d, rs in dir_read_sets.items():
            if d in used or not d.startswith("novel_"):
                continue
            ov = len(rs & c.supporting_read_ids)
            if ov > best_ov:
                best, best_ov = d, ov
        if best is not None and best_ov > 0:
            c.candidate_id = best
            used.add(best)


def _cluster_order(M: np.ndarray) -> np.ndarray:
    """Return a row/col permutation from hierarchical clustering of M.

    Falls back to identity order if M is too small, non-square, or has too
    many non-finite entries to compute a valid condensed distance.
    """
    n = M.shape[0]
    if M.ndim != 2 or M.shape[0] != M.shape[1] or n < 3:
        return np.arange(n)
    # Symmetrize and fill non-finite with a large finite sentinel based on
    # the observed max so non-finite rows don't dominate linkage.
    A = np.array(M, dtype=np.float64, copy=True)
    A = 0.5 * (A + A.T)
    finite = A[np.isfinite(A)]
    if finite.size == 0:
        return np.arange(n)
    fill = float(np.nanmax(finite)) + 1.0
    A = np.where(np.isfinite(A), A, fill)
    np.fill_diagonal(A, 0.0)
    try:
        cond = squareform(A, checks=False)
        Z = linkage(cond, method="average")
        order = leaves_list(Z)
        return np.asarray(order, dtype=int)
    except Exception:
        return np.arange(n)


def _vrange(M: np.ndarray):
    flat = M[np.isfinite(M)]
    flat = flat[flat < 1e5]
    if flat.size == 0:
        return 0.0, 1.0
    lo, hi = np.percentile(flat, [2, 98])
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def _plot_rr(ax, M, title, cmap="viridis_r"):
    lo, hi = _vrange(M)
    Md = np.where(np.isfinite(M), M, np.nan)
    im = ax.imshow(Md, aspect="auto", cmap=cmap, vmin=lo, vmax=hi,
                   interpolation="nearest")
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("read")
    ax.set_ylabel("read")
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)


def _plot_class_structures(ax, candidates, classes, eval_rows, interval,
                           diff_regions=None):
    """Draw transcript structures grouped by candidate class.

    Each class gets a unique color from tab10. Within a class, candidates are
    stacked vertically. Diff regions are shown as gold vertical bands.
    Each candidate is labeled with its id + GTF-eval tag (exact/wobble/novel)
    + best matching GTF transcript id.
    """
    # Order candidates by class so members of one class are contiguous.
    ordered: List[Tuple[int, int]] = []  # (class_idx, cand_idx)
    for ci, members in enumerate(classes):
        for j in members:
            ordered.append((ci, j))
    n = len(ordered)
    if n == 0:
        ax.axis("off")
        return
    x_lo = min(c.start for c in candidates)
    x_hi = max(c.end for c in candidates)
    eval_by_id = {r["candidate_id"]: r for r in eval_rows}
    cmap = plt.get_cmap("tab10")

    PAD_BP = 5
    if diff_regions:
        for k, (g_lo, g_hi) in enumerate(diff_regions):
            g_lo_v, g_hi_v = g_lo - PAD_BP, g_hi + PAD_BP
            ax.axvspan(g_lo_v, g_hi_v, color="gold", alpha=0.30, zorder=0)
            ax.text((g_lo + g_hi) / 2, n - 0.4,
                    f"D{k}\n{g_hi - g_lo}bp(+{PAD_BP}bp)",
                    ha="center", va="bottom",
                    fontsize=7, color="darkgoldenrod", fontweight="bold")

    # Class separator stripes on the left edge.
    for slot, (ci, j) in enumerate(ordered):
        y = n - 1 - slot
        c = candidates[j]
        color = cmap(ci % 10)
        introns = c.intron_chain.introns
        if not introns:
            exons = [(c.start, c.end)]
        else:
            exons = [(c.start, introns[0][0])]
            for k in range(len(introns) - 1):
                exons.append((introns[k][1], introns[k + 1][0]))
            exons.append((introns[-1][1], c.end))
        ax.hlines(y, c.start, c.end, color=color, lw=0.8, alpha=0.5)
        for ex_s, ex_e in exons:
            if ex_e > ex_s:
                ax.add_patch(plt.Rectangle(
                    (ex_s, y - 0.35), ex_e - ex_s, 0.7,
                    facecolor=color, edgecolor="black", lw=0.3,
                ))
        ev = eval_by_id.get(c.candidate_id, {})
        tag = ev.get("match", "?")
        best = ev.get("best_gtf_tx") or "-"
        label = f"c{ci} | {c.candidate_id[:12]} [{tag} -> {best}]"
        ax.text(x_hi + (x_hi - x_lo) * 0.01, y, label,
                va="center", fontsize=7, color=color)

    # Horizontal dividers between classes.
    cum = 0
    for ci, members in enumerate(classes[:-1]):
        cum += len(members)
        y_div = n - cum - 0.5
        ax.axhline(y_div, color="black", lw=0.4, alpha=0.3, ls=":")

    ax.set_xlim(x_lo - (x_hi - x_lo) * 0.02, x_hi + (x_hi - x_lo) * 0.28)
    ax.set_ylim(-0.8, n - 0.2)
    ax.set_yticks([])
    ax.set_xlabel(f"genomic position on {interval.chrom}")
    ax.set_title(
        f"candidate classes  {interval.region_string}  "
        f"(color = class; dotted line separates classes; "
        f"gold band = diff region)",
        fontsize=10,
    )
    ax.grid(axis="x", alpha=0.3, ls=":")


def _plot_assignment(ax, read_class: List[int], n_classes: int, title: str):
    """Strip showing each read's class assignment, using tab10 for class color
    (matches the candidate-class structure panel)."""
    arr = np.asarray(read_class, dtype=int).reshape(1, -1)
    # We render via an RGBA image so unassigned (-1) is a clearly distinct
    # neutral gray and classes 0..n_classes-1 map 1:1 to tab10.
    cmap = plt.get_cmap("tab10")
    n = arr.shape[1]
    rgba = np.zeros((1, n, 4), dtype=np.float32)
    for k in range(n):
        ci = int(arr[0, k])
        if ci < 0:
            rgba[0, k] = (0.85, 0.85, 0.85, 1.0)  # light gray = unassigned
        else:
            rgba[0, k] = cmap(ci % 10)
    ax.imshow(rgba, aspect="auto", interpolation="nearest")
    ax.set_yticks([])
    ax.set_xlabel("read")
    ax.set_title(title, fontsize=9)


def _exons_from_gtf_tx(tx) -> List[Tuple[int, int]]:
    """Return sorted exon list (start, end) from a GTFTranscript."""
    return sorted(tx.exons)


def _introns_from_exons(exons: List[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    if len(exons) < 2:
        return ()
    out = []
    for k in range(len(exons) - 1):
        out.append((exons[k][1], exons[k + 1][0]))
    return tuple(out)


def evaluate_candidates_vs_gtf(
    candidates,
    gtf_reader: GTFReader,
    chrom: str,
    start: int,
    end: int,
    tol: int = 6,
) -> List[Dict[str, object]]:
    """For each candidate report best-matching GTF transcript at ±tol/junction.

    Returns a list of dicts: {candidate_id, match, best_gtf_tx, jn_diff_max}.
    match ∈ {"exact", f"wobble±{tol}", "novel"}.
    """
    gtf_txs = gtf_reader.get_transcripts_in_region(chrom, start, end)
    gtf_chains = []
    for tx in gtf_txs:
        ex = _exons_from_gtf_tx(tx)
        gtf_chains.append((tx.transcript_id, _introns_from_exons(ex)))

    results = []
    for c in candidates:
        cand_chain = c.intron_chain.introns
        best_id, best_max_diff = None, None
        exact_id = None
        for tx_id, chain in gtf_chains:
            if len(chain) != len(cand_chain):
                continue
            if len(cand_chain) == 0:
                # both single-exon; accept if any chrom/strand overlap (loose).
                if best_id is None:
                    best_id = tx_id
                    best_max_diff = 0
                continue
            max_diff = 0
            for (sa, ea), (sb, eb) in zip(cand_chain, chain):
                max_diff = max(max_diff, abs(sa - sb), abs(ea - eb))
            if max_diff == 0 and exact_id is None:
                exact_id = tx_id
            if best_max_diff is None or max_diff < best_max_diff:
                best_max_diff, best_id = max_diff, tx_id
        if exact_id is not None:
            match = "exact"
            best_id = exact_id
            best_max_diff = 0
        elif best_id is not None and best_max_diff is not None and best_max_diff <= tol:
            match = f"wobble±{tol}"
        else:
            match = "novel"
        results.append({
            "candidate_id": c.candidate_id,
            "source": c.source,
            "match": match,
            "best_gtf_tx": best_id,
            "jn_diff_max": best_max_diff if best_max_diff is not None else -1,
            "introns": cand_chain,
        })
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf-eval", required=True,
                    help="GTF used ONLY for downstream evaluation, "
                         "not for discovery.")
    ap.add_argument("--signal", required=True)
    ap.add_argument("--eventalign-root", required=True,
                    help="No-GTF eventalign root (e.g. testdata/out_no_gtf).")
    ap.add_argument("--interval", required=True,
                    help="Interval region_string as named in --eventalign-root "
                         "(no-GTF interval names may differ from with-GTF).")
    ap.add_argument("--out", required=True)
    ap.add_argument("--chain-wobble", type=int, default=6)
    ap.add_argument("--gtf-eval-tol", type=int, default=6,
                    help="Per-junction wobble for GTF-based evaluation.")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    log.info("Loading GTF (evaluation only) + genome ...")
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
    if target is None:
        log.error("Interval %s not found in no-GTF interval set", args.interval)
        log.error("Available intervals near target:")
        for iv in res["intervals"]:
            if iv.chrom and "34" in iv.region_string:
                log.error("  %s", iv.region_string)
        sys.exit(2)

    log.info("Discovering candidates (no-GTF) for %s ...", target.region_string)
    cs = discover_candidates(
        interval=target, bam_path=args.bam, gtf_reader=None,
        genome_fasta=genome.get(target.chrom, ""),
        threshold=24, min_novel_reads=1,
    )
    cands = list(cs.candidates)
    remap_novel_ids(cands, Path(args.eventalign_root) / args.interval)
    log.info("Got %d candidates (all should be novel; %d gtf / %d novel)",
             len(cands),
             sum(1 for c in cands if c.source == "gtf"),
             sum(1 for c in cands if c.source == "novel"))

    # GTF-based evaluation of discovered candidates.
    eval_rows = evaluate_candidates_vs_gtf(
        cands, gtf_eval, target.chrom, target.start, target.end,
        tol=args.gtf_eval_tol,
    )
    log.info("=== GTF evaluation (discovery without GTF, tol=±%dbp) ===",
             args.gtf_eval_tol)
    for row in eval_rows:
        log.info("  %-16s match=%-12s best=%s  max_jn_diff=%s",
                 row["candidate_id"][:16],
                 row["match"],
                 row["best_gtf_tx"],
                 row["jn_diff_max"])
    n_exact = sum(1 for r in eval_rows if r["match"] == "exact")
    n_wobble = sum(1 for r in eval_rows if r["match"].startswith("wobble"))
    n_novel = sum(1 for r in eval_rows if r["match"] == "novel")
    log.info("Summary: exact=%d  wobble=%d  novel=%d  total=%d",
             n_exact, n_wobble, n_novel, len(eval_rows))

    cand_ids = [c.candidate_id for c in cands]

    log.info("Computing M1 (mappy) ...")
    M1, read_ids, _ = compute_m1_mappy(cs.read_sequences, cands)

    log.info("Computing M2 (eventalign) ...")
    M2, scores_by_pair = compute_m2_eventalign(
        Path(args.eventalign_root), args.interval, cands, read_ids,
    )

    classes = cluster_candidates_by_chain(cands, chain_wobble=args.chain_wobble)
    log.info("Clustered into %d classes (chain_wobble=%d):", len(classes), args.chain_wobble)
    for ci, members in enumerate(classes):
        ids = [cands[m].candidate_id[:14] for m in members]
        introns_repr = cands[members[0]].intron_chain.introns
        log.info("  class %d  n=%d  chain=%s  members=%s",
                 ci, len(members), introns_repr, ids)

    diff_regs = extract_diff_regions(cands)
    log.info("%d diff regions (junction-anchored)", len(diff_regs))

    log.info("Computing diff-region M3 (existing) ...")
    with Slow5Reader(args.signal) as sr:
        M3 = compute_diff_region_m4(
            read_ids=read_ids, candidates=cands,
            scores_by_pair=scores_by_pair, signal_reader=sr,
            interval_start=target.start, interval_end=target.end,
            signal_format="slow5", use_gpu=True, normalize=True,
        )
    n = len(read_ids)
    off = ~np.eye(n, dtype=bool)
    m3_nan_frac = float(np.isnan(M3[off]).mean())
    log.info("M3 off-diag NaN frac: %.2f", m3_nan_frac)

    log.info("Computing class-partitioned M4 with M1 as read->cand distance ...")
    with Slow5Reader(args.signal) as sr:
        M4_via_m1 = compute_class_partitioned_m4(
            read_ids=read_ids, candidates=cands, cand_ids=cand_ids,
            dist_read_cand=M1, signal_reader=sr,
            scores_by_pair=scores_by_pair,
            signal_format="slow5", use_gpu=True, normalize=True,
            chain_wobble=args.chain_wobble, inter_class_fill="max",
        )

    log.info("Computing class-partitioned M4 with M2 as read->cand distance ...")
    with Slow5Reader(args.signal) as sr:
        M4_via_m2 = compute_class_partitioned_m4(
            read_ids=read_ids, candidates=cands, cand_ids=cand_ids,
            dist_read_cand=M2, signal_reader=sr,
            scores_by_pair=scores_by_pair,
            signal_format="slow5", use_gpu=True, normalize=True,
            chain_wobble=args.chain_wobble, inter_class_fill="max",
        )

    # Assignments (for the strip plots)
    from fin.scoring.diff_region_dtw import _assign_reads_to_classes
    # Build cand_idx ordering matching cand_ids.
    cand_idx_to_class = {}
    for ci, members in enumerate(classes):
        for j in members:
            cand_idx_to_class[j] = ci
    id_to_idx = {c.candidate_id: i for i, c in enumerate(cands)}
    classes_in_col = [[] for _ in classes]
    for col, cid in enumerate(cand_ids):
        ci = cand_idx_to_class.get(id_to_idx.get(cid, -1))
        if ci is not None and ci >= 0:
            classes_in_col[ci].append(col)
    assign_m1 = _assign_reads_to_classes(read_ids, cand_ids, classes_in_col, M1)
    assign_m2 = _assign_reads_to_classes(read_ids, cand_ids, classes_in_col, M2)
    log.info("Class distribution via M1: %s", np.bincount(np.array(assign_m1) + 1).tolist())
    log.info("Class distribution via M2: %s", np.bincount(np.array(assign_m2) + 1).tolist())

    # Save raw matrices.
    np.savez(out_path.with_suffix(".npz"),
             M1=M1, M2=M2, M3=M3,
             M4_via_m1=M4_via_m1, M4_via_m2=M4_via_m2,
             read_ids=np.array(read_ids),
             cand_ids=np.array(cand_ids),
             assign_m1=np.array(assign_m1),
             assign_m2=np.array(assign_m2))

    # Hierarchical clustering reorder (derived from M4_via_m2; applied to all
    # heatmaps + assignment strips so panels are directly comparable).
    order = _cluster_order(M4_via_m2)
    if order.size == n and not np.array_equal(order, np.arange(n)):
        log.info("Applying clustering reorder (n=%d)", n)
        M3_o = M3[np.ix_(order, order)]
        M4_via_m1_o = M4_via_m1[np.ix_(order, order)]
        M4_via_m2_o = M4_via_m2[np.ix_(order, order)]
        assign_m1_o = [assign_m1[i] for i in order]
        assign_m2_o = [assign_m2[i] for i in order]
    else:
        M3_o, M4_via_m1_o, M4_via_m2_o = M3, M4_via_m1, M4_via_m2
        assign_m1_o, assign_m2_o = assign_m1, assign_m2

    # Plot
    log.info("Plotting ...")
    struct_h = 0.4 + 0.22 * len(cands)
    fig = plt.figure(figsize=(22, 6 + struct_h))
    gs = fig.add_gridspec(
        nrows=4, ncols=3,
        height_ratios=[0.18, struct_h, 2.0, 0.4],
        hspace=0.35, wspace=0.25,
    )

    # Row 0: title bar with class summary.
    ax_info = fig.add_subplot(gs[0, :])
    ax_info.axis("off")
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

    # Row 1: candidate class structures (transcript boxes colored by class).
    ax_struct = fig.add_subplot(gs[1, :])
    _plot_class_structures(ax_struct, cands, classes, eval_rows, target,
                           diff_regions=diff_regs)

    # Row 2: three heatmaps (clustering-reordered for direct visual comparison).
    ax_a = fig.add_subplot(gs[2, 0])
    _plot_rr(ax_a, M3_o,
             f"M3 diff-region DTW (off-diag NaN={m3_nan_frac:.0%})")
    ax_b = fig.add_subplot(gs[2, 1])
    _plot_rr(ax_b, M4_via_m1_o,
             f"class-partitioned M4 (assign via M1 mappy)")
    ax_c = fig.add_subplot(gs[2, 2])
    _plot_rr(ax_c, M4_via_m2_o,
             f"class-partitioned M4 (assign via M2 eventalign)")

    # Row 3: class-assignment strips (only for the M1/M2 panels).
    ax_sb = fig.add_subplot(gs[3, 1])
    _plot_assignment(ax_sb, assign_m1_o, len(classes), "read class via M1")
    ax_sc = fig.add_subplot(gs[3, 2])
    _plot_assignment(ax_sc, assign_m2_o, len(classes), "read class via M2")
    # Hide unused row-3 slot.
    ax_sa = fig.add_subplot(gs[3, 0])
    ax_sa.axis("off")

    fig.suptitle(
        f"Class-partitioned M4 vs diff-region M3  —  {target.region_string}",
        fontsize=12,
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    log.info("Wrote %s and %s", out_path, out_path.with_suffix(".npz"))


if __name__ == "__main__":
    main()
