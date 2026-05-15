"""POC: extract signal slices at a diff region using any covering candidate.

The current `compute_diff_region_m4` projects diff regions through each read's
single max-LL host, which fails when the diff region falls in the host's
intron. This script demonstrates that using any candidate whose cDNA covers
the region recovers most of the lost coverage and produces real signal slices.

Output:
  - A summary table comparing current vs proposed event coverage per region.
  - PNG overlay plot of extracted signal slices, color-coded by max-LL host.

Usage (mirrors viz_m1_m2_m3.py args):
  python benchmarks/diff_region_mapping_poc.py \
    --bam testdata/mapped.bam --genome testdata/SIRV.genome.fa \
    --gtf testdata/SIRV.genome.gtf --signal testdata/mapped.blow5 \
    --eventalign-root testdata/out_with_gtf \
    --interval SIRVomeERCCome_34497_36900 \
    --out testdata/viz/diff_region_signals.png
"""
from __future__ import annotations

import argparse
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from fin.candidates.discovery import discover_candidates
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_fasta import FASTAReader
from fin.io.io_gtf import GTFReader
from fin.io.io_slow5 import Slow5Reader
from fin.scoring.diff_region_dtw import extract_diff_regions
from fin.scoring.eventalign_parser import parse_eventalign_tsv

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
log = logging.getLogger(__name__)


def cdna_to_genomic(cand, cdna_pos: int) -> Optional[int]:
    """Inverse of genomic_region_to_cdna for a single position.

    Returns genomic coord (0-based) or None if cdna_pos is out of range.
    """
    introns = cand.intron_chain.introns
    if not introns:
        exons = [(cand.start, cand.end)]
    else:
        exons = [(cand.start, introns[0][0])]
        for k in range(len(introns) - 1):
            exons.append((introns[k][1], introns[k + 1][0]))
        exons.append((introns[-1][1], cand.end))
    spliced_len = sum(max(0, e - s) for s, e in exons)
    if cdna_pos < 0 or cdna_pos >= spliced_len:
        return None
    # '-' strand: cDNA k in rev-comp frame = (spliced_len - 1 - k) in LTR frame
    ltr_pos = (spliced_len - 1 - cdna_pos) if cand.strand == "-" else cdna_pos
    cum = 0
    for ex_s, ex_e in exons:
        ex_len = max(0, ex_e - ex_s)
        if cum + ex_len > ltr_pos:
            return ex_s + (ltr_pos - cum)
        cum += ex_len
    return None


def remap_novel_ids(cands, ev_root: Path) -> None:
    """Mutate cand.candidate_id to match eventalign dir by supporting-read overlap."""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--signal", required=True)
    ap.add_argument("--eventalign-root", required=True)
    ap.add_argument("--interval", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-reads-plot", type=int, default=15,
                    help="Cap on overlaid signal traces per region in the plot")
    args = ap.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    gtf = GTFReader(args.gtf); gtf.open(); gtf.parse()
    fasta = FASTAReader(args.genome)
    genome = {rec.id: rec.sequence for rec in fasta.get_records()}
    result = generate_isolated_intervals(args.bam, gtf_path=args.gtf,
                                         max_gap=200, max_reads=None)
    target = None
    want = args.interval.replace(":", "_").replace("-", "_")
    for iv in result["intervals"]:
        if iv.region_string.replace(":", "_").replace("-", "_") == want:
            target = iv; break
    assert target is not None, f"interval {args.interval} not found"

    cs = discover_candidates(
        interval=target, bam_path=args.bam, gtf_reader=gtf,
        genome_fasta=genome.get(target.chrom, ""),
        threshold=24, min_novel_reads=1,
    )
    cands = list(cs.candidates)
    remap_novel_ids(cands, Path(args.eventalign_root) / args.interval)
    log.info("Loaded %d candidates", len(cands))

    diff_regs = extract_diff_regions(cands)
    log.info("Diff regions: %s", diff_regs)

    # Parse all eventalign TSVs with collect_events=True
    cand_lengths = {c.candidate_id: len(c.sequence) for c in cands}
    cand_by_id = {c.candidate_id: c for c in cands}
    # events_by_read_cand[(rid, cid)] = list of (cdna_pos, s_lo, s_hi)
    events_by_read_cand: Dict[Tuple[str, str], List[Tuple[int, int, int]]] = {}
    # totalLL_by_pair[(rid, cid)] = total_log_likelihood
    totalLL_by_pair: Dict[Tuple[str, str], float] = {}
    for c in cands:
        tsv = Path(args.eventalign_root) / args.interval / c.candidate_id / "eventalign.tsv"
        if not tsv.exists():
            continue
        scores = parse_eventalign_tsv(str(tsv), cand_lengths, collect_events=True)
        for s in scores:
            events_by_read_cand[(s.read_name, s.candidate_id)] = s.events
            totalLL_by_pair[(s.read_name, s.candidate_id)] = s.total_log_likelihood

    all_reads = sorted({rid for rid, _ in events_by_read_cand.keys()})
    log.info("Reads with at least one eventalign hit: %d", len(all_reads))

    # max-LL host per read
    host_by_read: Dict[str, str] = {}
    for rid in all_reads:
        best_ll, best_c = -1e18, None
        for cid in cand_by_id.keys():
            ll = totalLL_by_pair.get((rid, cid))
            if ll is not None and ll > best_ll:
                best_ll, best_c = ll, cid
        if best_c is not None:
            host_by_read[rid] = best_c

    # For each diff region, score per-read coverage via:
    #  (A) host-only (current M3 logic): events from host that lie in the host's cDNA at this region
    #  (B) any-candidate: events from any cand whose cDNA→genomic projection lies in region
    summary_rows = []
    region_slices: Dict[int, Dict[str, Tuple[int, int, str]]] = {
        i: {} for i in range(len(diff_regs))
    }  # region_idx -> {rid: (sig_lo, sig_hi, used_cand)}

    for r_idx, (g_lo, g_hi) in enumerate(diff_regs):
        n_host_only = 0
        n_any = 0
        for rid in all_reads:
            # (A) host-only
            host = host_by_read.get(rid)
            host_hit = False
            if host is not None:
                events = events_by_read_cand.get((rid, host), [])
                host_cand = cand_by_id[host]
                for cdna_pos, _s, _e in events:
                    g = cdna_to_genomic(host_cand, cdna_pos)
                    if g is not None and g_lo <= g < g_hi:
                        host_hit = True
                        break
            if host_hit:
                n_host_only += 1

            # (B) any-candidate, pick the cand with most in-region events
            best_cand_for_region: Optional[str] = None
            best_n = 0
            best_sig: Optional[Tuple[int, int]] = None
            for cid in cand_by_id.keys():
                events = events_by_read_cand.get((rid, cid), [])
                if not events:
                    continue
                cand = cand_by_id[cid]
                sig_lo, sig_hi = None, None
                cnt = 0
                for cdna_pos, s_lo, s_hi in events:
                    g = cdna_to_genomic(cand, cdna_pos)
                    if g is not None and g_lo <= g < g_hi:
                        cnt += 1
                        if sig_lo is None or s_lo < sig_lo: sig_lo = s_lo
                        if sig_hi is None or s_hi > sig_hi: sig_hi = s_hi
                if cnt > best_n:
                    best_n = cnt
                    best_cand_for_region = cid
                    best_sig = (sig_lo, sig_hi) if sig_lo is not None else None
            if best_cand_for_region is not None and best_sig is not None:
                n_any += 1
                region_slices[r_idx][rid] = (best_sig[0], best_sig[1], best_cand_for_region)

        summary_rows.append((r_idx, g_lo, g_hi, n_host_only, n_any, len(all_reads)))

    # Print summary
    print("\n=== Diff region coverage: current (host-only) vs proposed (any-cand) ===")
    print(f"{'region':<7} {'g_lo':<8} {'g_hi':<8} {'bp':<6} {'cur':<6} {'prop':<6} "
          f"{'total':<6} {'gain':<8}")
    for r_idx, g_lo, g_hi, n_cur, n_prop, n_tot in summary_rows:
        gain = n_prop - n_cur
        gain_pct = (gain / n_tot * 100) if n_tot else 0
        print(f"D{r_idx:<6} {g_lo:<8} {g_hi:<8} {g_hi-g_lo:<6} "
              f"{n_cur:<6} {n_prop:<6} {n_tot:<6} +{gain} ({gain_pct:.0f}%)")

    # Extract signal slices and plot
    log.info("Extracting signal slices from %s ...", args.signal)
    fig, axes = plt.subplots(
        len(diff_regs), 1,
        figsize=(14, 3.0 * len(diff_regs)),
        squeeze=False,
    )
    cmap = plt.get_cmap("tab10")
    host_colors: Dict[str, str] = {}
    next_color = 0
    with Slow5Reader(args.signal) as sr:
        for r_idx, (g_lo, g_hi) in enumerate(diff_regs):
            ax = axes[r_idx, 0]
            slices = region_slices[r_idx]
            # Sort reads by host so legend groups them
            items = sorted(
                slices.items(),
                key=lambda kv: (host_by_read.get(kv[0], "zzz"), kv[0]),
            )[: args.max_reads_plot]
            for rid, (s_lo, s_hi, used_cand) in items:
                try:
                    res = sr.get_picoamp_signal(rid)
                except Exception:
                    res = None
                if res is None:
                    continue
                sig, _meta = res
                sig = np.asarray(sig, dtype=np.float32)
                s_lo = max(0, int(s_lo))
                s_hi = min(len(sig), int(s_hi))
                if s_hi <= s_lo:
                    continue
                seg = sig[s_lo:s_hi]
                # Robust z-score so traces overlay
                med = float(np.median(seg))
                mad = float(np.median(np.abs(seg - med)))
                scale = mad if mad > 0 else 1.0
                seg_z = (seg - med) / scale
                host = host_by_read.get(rid, "?")
                if host not in host_colors:
                    host_colors[host] = cmap(next_color % 10)
                    next_color += 1
                ax.plot(seg_z, color=host_colors[host], lw=0.6, alpha=0.55)
            n_cur = summary_rows[r_idx][3]
            n_prop = summary_rows[r_idx][4]
            n_tot = summary_rows[r_idx][5]
            ax.set_title(
                f"D{r_idx}  {target.chrom}:{g_lo}-{g_hi}  ({g_hi-g_lo}bp)  |  "
                f"reads with signal: host-only {n_cur}/{n_tot} → any-cand {n_prop}/{n_tot}"
            )
            ax.set_xlabel("signal sample (relative)")
            ax.set_ylabel("z-norm pA")

    # One legend (host colors)
    handles = [
        plt.Line2D([0], [0], color=c, lw=2,
                   label=f"host={h[:14]}")
        for h, c in host_colors.items()
    ]
    if handles:
        axes[0, 0].legend(handles=handles, loc="upper right", fontsize=8)

    fig.suptitle(
        f"{target.region_string}  |  diff-region signal overlays "
        f"(via any-candidate cDNA projection)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_path, dpi=150)
    log.info("Wrote %s", out_path)


if __name__ == "__main__":
    main()
