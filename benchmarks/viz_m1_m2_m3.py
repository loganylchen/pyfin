"""Visualize M1, M2, M3 matrices for a single interval as heatmaps.

M1: read x candidate, mappy alignment distance (= -score, lower=better)
M2: read x candidate, eventalign mean -log-likelihood (lower=better)
M3: read x read, diff-region DTW distance (lower=more coherent)

GTF candidates are marked on the M1/M2 column axis (cyan bar).

Usage:
    python benchmarks/viz_m1_m2_m3.py \
        --bam testdata/mapped.bam \
        --genome testdata/SIRV.genome.fa \
        --gtf testdata/SIRV.genome.gtf \
        --signal testdata/mapped.blow5 \
        --eventalign-root testdata/out_with_gtf \
        --interval SIRVomeERCCome_34497_36900 \
        --out testdata/viz/m_matrices.png
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fin.candidates.discovery import discover_candidates
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_fasta import FASTAReader
from fin.io.io_gtf import GTFReader
from fin.io.io_slow5 import Slow5Reader
from fin.scoring.diff_region_dtw import compute_diff_region_m4
from fin.scoring.eventalign_parser import (
    build_distance_matrix,
    parse_eventalign_tsv,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger(__name__)


def compute_m1_mappy(
    read_seqs: Dict[str, str],
    candidates,
) -> Tuple[np.ndarray, List[str], List[str]]:
    """M1 = -mappy.score per (read, candidate). float32, missing = +1e6."""
    import mappy

    read_ids = sorted(read_seqs.keys())
    cand_ids = [c.candidate_id for c in candidates]

    aligners = []
    for c in candidates:
        if not c.sequence:
            aligners.append(None)
            continue
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


def compute_m2_eventalign(
    eventalign_root: Path,
    interval_name: str,
    candidates,
    read_ids: List[str],
):
    """M2 from eventalign TSVs. Returns (matrix, scores_by_pair)."""
    cand_ids = [c.candidate_id for c in candidates]
    cand_lengths = {c.candidate_id: len(c.sequence) for c in candidates}

    all_scores = []
    for cid in cand_ids:
        tsv = eventalign_root / interval_name / cid / "eventalign.tsv"
        if not tsv.exists():
            log.warning("Missing eventalign TSV: %s", tsv)
            continue
        scores = parse_eventalign_tsv(str(tsv), cand_lengths, collect_events=True)
        all_scores.extend(scores)

    M = build_distance_matrix(all_scores, read_ids, cand_ids)
    scores_by_pair = {(s.read_name, s.candidate_id): s for s in all_scores}
    return M.astype(np.float32), scores_by_pair


def compute_m3(
    read_ids: List[str],
    candidates,
    scores_by_pair,
    signal_reader,
    interval_start: int,
    interval_end: int,
) -> np.ndarray:
    return compute_diff_region_m4(
        read_ids=read_ids,
        candidates=list(candidates),
        scores_by_pair=scores_by_pair,
        signal_reader=signal_reader,
        interval_start=interval_start,
        interval_end=interval_end,
        signal_format="slow5",
        use_gpu=False,
        normalize=True,
    )


def _robust_vrange(M: np.ndarray):
    """Use 2/98 percentiles ignoring NaN/large sentinel for color scale."""
    flat = M[np.isfinite(M)]
    flat = flat[flat < 1e5]
    if flat.size == 0:
        return (0.0, 1.0)
    lo, hi = np.percentile(flat, [2, 98])
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def _row_zscore(M: np.ndarray) -> np.ndarray:
    """Per-row z-score; rows that are all-missing or constant become 0."""
    out = np.where(M >= 1e5, np.nan, M).astype(np.float32)
    mu = np.nanmean(out, axis=1, keepdims=True)
    sd = np.nanstd(out, axis=1, keepdims=True)
    sd = np.where(sd > 0, sd, 1.0)
    z = (out - mu) / sd
    return np.where(np.isnan(z), 0.0, z)


def _plot_rc(ax, M, cand_ids, sources, title, ylabel="read",
             cmap="viridis_r", vrange=None):
    """Plot a read x candidate heatmap with GTF column markers."""
    if vrange is None:
        lo, hi = _robust_vrange(M)
    else:
        lo, hi = vrange
    M_disp = np.where(M >= 1e5, np.nan, M) if vrange is None else M
    im = ax.imshow(M_disp, aspect="auto", cmap=cmap, vmin=lo, vmax=hi,
                   interpolation="nearest")
    ax.set_title(title)
    ax.set_xlabel("candidate")
    ax.set_ylabel(ylabel)
    ax.set_xticks(range(len(cand_ids)))
    labels = [
        f"{'★ ' if s == 'gtf' else ''}{cid[:12]}"
        for cid, s in zip(cand_ids, sources)
    ]
    ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)
    for j, s in enumerate(sources):
        if s == "gtf":
            ax.axvline(j, color="cyan", lw=0.6, alpha=0.6)
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)


def _plot_structures(ax, candidates, sources, interval, diff_regions=None):
    """Draw transcript structures: exon boxes + intron lines, one track per candidate.

    If ``diff_regions`` is provided, highlight each diff region as a vertical
    yellow band spanning all tracks.
    """
    n = len(candidates)
    x_lo = min(c.start for c in candidates)
    x_hi = max(c.end for c in candidates)
    # Highlight diff regions first (so transcript boxes draw on top of bands).
    # Each region is padded by ±5bp on both sides for visual context.
    PAD_BP = 5
    if diff_regions:
        for k, (g_lo, g_hi) in enumerate(diff_regions):
            g_lo_v, g_hi_v = g_lo - PAD_BP, g_hi + PAD_BP
            ax.axvspan(g_lo_v, g_hi_v, color="gold", alpha=0.30, zorder=0)
            ax.text((g_lo + g_hi) / 2, n - 0.4,
                    f"D{k}\n{g_hi - g_lo}bp(+{PAD_BP}bp)",
                    ha="center", va="bottom",
                    fontsize=7, color="darkgoldenrod", fontweight="bold")
    for j, c in enumerate(candidates):
        y = n - 1 - j  # plot top-down to match heatmap col order
        color = "tab:cyan" if c.source == "gtf" else "tab:orange"
        introns = c.intron_chain.introns
        if not introns:
            exons = [(c.start, c.end)]
        else:
            exons = [(c.start, introns[0][0])]
            for k in range(len(introns) - 1):
                exons.append((introns[k][1], introns[k + 1][0]))
            exons.append((introns[-1][1], c.end))
        # Thin spine across full transcript span
        ax.hlines(y, c.start, c.end, color=color, lw=0.8, alpha=0.5)
        # Thick exon boxes
        for ex_s, ex_e in exons:
            if ex_e > ex_s:
                ax.add_patch(plt.Rectangle(
                    (ex_s, y - 0.35), ex_e - ex_s, 0.7,
                    facecolor=color, edgecolor="black", lw=0.3,
                ))
        # Intron chevrons (just thin lines; exons already drawn)
        ax.text(x_hi + (x_hi - x_lo) * 0.01, y,
                f"{'★ ' if c.source == 'gtf' else ''}{c.candidate_id[:14]} ({c.strand})",
                va="center", fontsize=7, color=color)
    ax.set_xlim(x_lo - (x_hi - x_lo) * 0.02, x_hi + (x_hi - x_lo) * 0.20)
    ax.set_ylim(-0.8, n - 0.2)
    ax.set_yticks([])
    ax.set_xlabel(f"genomic position on {interval.chrom}")
    ax.set_title(
        f"transcript structures  {interval.region_string}  "
        f"(cyan=GTF, orange=novel)",
        fontsize=10,
    )
    ax.grid(axis="x", alpha=0.3, ls=":")


def _plot_rr(ax, M, n_reads, title, cmap="viridis_r"):
    lo, hi = _robust_vrange(M)
    M_disp = np.where(np.isfinite(M), M, np.nan)
    im = ax.imshow(M_disp, aspect="auto", cmap=cmap, vmin=lo, vmax=hi,
                   interpolation="nearest")
    ax.set_title(title)
    ax.set_xlabel("read")
    ax.set_ylabel("read")
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--signal", required=True, help="blow5/slow5 file")
    ap.add_argument("--eventalign-root", required=True,
                    help="dir holding {interval}/{candidate}/eventalign.tsv")
    ap.add_argument("--interval", required=True,
                    help="target interval region_string e.g. chr_start_end")
    ap.add_argument("--out", required=True, help="output PNG path")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    log.info("Loading GTF + genome ...")
    gtf = GTFReader(args.gtf); gtf.open(); gtf.parse()
    fasta = FASTAReader(args.genome)
    genome = {rec.id: rec.sequence for rec in fasta.get_records()}

    log.info("Generating intervals ...")
    result = generate_isolated_intervals(args.bam, gtf_path=args.gtf,
                                         max_gap=200, max_reads=None)
    intervals = result["intervals"]
    def _norm(s: str) -> str:
        return s.replace(":", "_").replace("-", "_")
    want = _norm(args.interval)
    target = None
    for iv in intervals:
        if _norm(iv.region_string) == want:
            target = iv; break
    if target is None:
        log.error("Interval %s not found among %d intervals",
                  args.interval, len(intervals))
        sys.exit(2)

    log.info("Discovering candidates for %s ...", target.region_string)
    cs = discover_candidates(
        interval=target, bam_path=args.bam, gtf_reader=gtf,
        genome_fasta=genome.get(target.chrom, ""),
        threshold=24, min_novel_reads=1,
    )
    cands = list(cs.candidates)
    log.info("Got %d candidates (%d GTF, %d novel)",
             len(cands),
             sum(1 for c in cands if c.source == "gtf"),
             sum(1 for c in cands if c.source == "novel"))

    # The eventalign TSVs were produced by an earlier discovery run with
    # different novel_* UUIDs. Remap our novel candidates to those stable
    # dir IDs by matching supporting_read_ids (per-novel-dir read set is the
    # set of read_names appearing in eventalign.tsv).
    ev_root = Path(args.eventalign_root) / args.interval
    dir_read_sets: Dict[str, set] = {}
    if ev_root.exists():
        for d in sorted(ev_root.iterdir()):
            if not d.is_dir():
                continue
            tsv = d / "eventalign.tsv"
            if not tsv.exists():
                continue
            reads = set()
            with open(tsv) as fh:
                fh.readline()  # header
                for line in fh:
                    parts = line.split("\t", 5)
                    if len(parts) > 4:
                        reads.add(parts[3])
            dir_read_sets[d.name] = reads

    used_dir_ids = set()
    for c in cands:
        if c.source != "novel":
            continue
        best_dir, best_ov = None, 0
        for dname, rset in dir_read_sets.items():
            if dname in used_dir_ids or not dname.startswith("novel_"):
                continue
            ov = len(rset & c.supporting_read_ids)
            if ov > best_ov:
                best_dir, best_ov = dname, ov
        if best_dir is not None and best_ov > 0:
            log.info("Remapping %s -> %s (overlap=%d reads)",
                     c.candidate_id, best_dir, best_ov)
            c.candidate_id = best_dir
            used_dir_ids.add(best_dir)

    sources = [c.source for c in cands]
    cand_ids = [c.candidate_id for c in cands]

    if not cs.read_sequences:
        log.error("CandidateSet has empty read_sequences; cannot compute M1")
        sys.exit(2)

    log.info("Computing M1 (mappy) ...")
    M1, read_ids, _ = compute_m1_mappy(cs.read_sequences, cands)
    log.info("M1 shape=%s, n_reads=%d", M1.shape, len(read_ids))

    log.info("Computing M2 (eventalign) ...")
    M2, scores_by_pair = compute_m2_eventalign(
        Path(args.eventalign_root), args.interval, cands, read_ids
    )
    log.info("M2 shape=%s, scored pairs=%d", M2.shape, len(scores_by_pair))

    log.info("Computing M3 (diff-region DTW) ...")
    # Diagnostic: report how many diff regions and how many reads have a host
    from fin.scoring.diff_region_dtw import (
        extract_diff_regions, _best_candidate_per_read,
    )
    diff_regs = extract_diff_regions(cands)
    n_hosts = sum(
        1 for rid in read_ids
        if _best_candidate_per_read(rid, cands, scores_by_pair) is not None
    )
    log.info("M3 diag: %d diff regions, %d/%d reads have a host candidate",
             len(diff_regs), n_hosts, len(read_ids))
    for i, (g_lo, g_hi) in enumerate(diff_regs[:10]):
        log.info("  diff region %d: %d-%d (%d bp)", i, g_lo, g_hi, g_hi - g_lo)
    with Slow5Reader(args.signal) as sr:
        M3 = compute_m3(read_ids, cands, scores_by_pair, sr,
                        target.start, target.end)
    log.info("M3 shape=%s", M3.shape)

    # Save raw matrices alongside the figure for downstream inspection.
    np.savez(out_path.with_suffix(".npz"),
             M1=M1, M2=M2, M3=M3,
             read_ids=np.array(read_ids),
             cand_ids=np.array(cand_ids),
             sources=np.array(sources))

    log.info("Plotting ...")
    M1z = _row_zscore(M1)
    M2z = _row_zscore(M2)
    nan_frac = float(np.isnan(M3[~np.eye(len(read_ids), dtype=bool)]).mean())
    log.info("M3 off-diagonal NaN fraction = %.2f", nan_frac)

    fig = plt.figure(figsize=(22, 12))
    gs = fig.add_gridspec(
        nrows=2, ncols=3,
        height_ratios=[0.5 + 0.18 * len(cands), 2.0],
        hspace=0.35, wspace=0.30,
    )
    ax_struct = fig.add_subplot(gs[0, :])
    _plot_structures(ax_struct, cands, sources, target,
                     diff_regions=diff_regs)

    ax_m1 = fig.add_subplot(gs[1, 0])
    _plot_rc(ax_m1, M1z, cand_ids, sources,
             "M1 mappy  (row-z, lower=preferred)",
             cmap="RdBu_r", vrange=(-2, 2))
    ax_m2 = fig.add_subplot(gs[1, 1])
    _plot_rc(ax_m2, M2z, cand_ids, sources,
             "M2 eventalign  (row-z, lower=preferred)",
             cmap="RdBu_r", vrange=(-2, 2))
    ax_m3 = fig.add_subplot(gs[1, 2])
    _plot_rr(ax_m3, M3, len(read_ids),
             f"M3 diff-region DTW  (off-diag NaN={nan_frac:.0%})")

    fig.suptitle(
        f"{target.region_string}  |  reads={len(read_ids)}  cands={len(cands)} "
        f"(★/cyan = GTF transcript, orange = novel)",
        fontsize=12,
    )
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    log.info("Wrote %s and %s", out_path, out_path.with_suffix(".npz"))


if __name__ == "__main__":
    main()
