#!/usr/bin/env python3
"""Figures for the junction signal probe (spec: .omc/codex_figure_spec_junction_signal_2026-09-02.txt).

Per-read figure, one read scored under two junction hypotheses:
  A1/A2  genomic orientation, biological donor and acceptor windows, ALWAYS both,
         transcript 5'->3' left to right (so on '-' strand genomic coords decrease).
  B      the two SPLICED transcript sequences, pairwise aligned, indel bases coloured,
         splice seam marked per hypothesis, run IDs bracketed onto the k-mers.
  C      event -> k-mer assignment strip on the raw-sample axis (two rows, one per
         hypothesis), consecutive identical k-mer pairs merged into labelled runs.
  D1/2/3 observed current, then each hypothesis' model track, identical pA limits.
  E      per-event dNLL bars (no cumulative overlay).
  F      run-level dependence audit: one bar per run, dominant-run fraction and
         leave-one-run-out stability, because N events parked on one k-mer pair are
         NOT N independent observations.

Distribution figure: score distributions from a probe TSV by read class / shift size.
"""
from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pysam
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
from junction_signal_probe import build_tx, substitute, read_chain_and_indels, clean_at  # noqa: E402

C_H1, C_H2 = "#2a78d6", "#eb6834"          # hypothesis 1 (true / read-own), hypothesis 2
C_OBS, C_TXT, C_TXT2 = "#0b0b0b", "#0b0b0b", "#52514e"
C_SURF, C_GRID, C_PALE = "#fcfcfb", "#e6e6e3", "#f0efec"
FLANK_GENOMIC = 15
SEQ_CONTEXT = 16

plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": C_TXT2, "axes.labelcolor": C_TXT,
                     "xtick.color": C_TXT2, "ytick.color": C_TXT2, "figure.facecolor": C_SURF,
                     "axes.facecolor": C_SURF, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 9})


# ----------------------------------------------------------------------------- helpers
def donor_acceptor(intron, strand):
    """Biological (donor, acceptor) genomic coordinates for a genomic intron [lo, hi)."""
    lo, hi = intron
    return (lo, hi) if strand == "+" else (hi, lo)


def events_table(res):
    """Per-event arrays. sigma_eff = model_stdv * sqrt(var) reproduces krill's standardized_level."""
    d = dict(pos=np.asarray(res["position"], np.int64),
             lev=np.asarray(res["event_level_mean"], float), sd=np.asarray(res["event_stdv"], float),
             mm=np.asarray(res["model_mean"], float), ms=np.asarray(res["model_stdv"], float),
             s0=np.asarray(res["start_idx"], np.int64), s1=np.asarray(res["end_idx"], np.int64),
             kmer=np.asarray(res["model_kmer"]))
    d["lo"] = np.minimum(d["s0"], d["s1"]); d["hi"] = np.maximum(d["s0"], d["s1"])
    var = float(res.get("var", 1.0)) or 1.0
    d["sig"] = d["ms"] * math.sqrt(var)
    good = (d["sig"] > 0) & np.isfinite(d["lev"]) & np.isfinite(d["mm"])
    z = np.where(good, (d["lev"] - d["mm"]) / np.where(good, d["sig"], 1.0), np.nan)
    d["nll"] = np.where(good, 0.5 * z * z + np.log(np.where(good, d["sig"], 1.0)), np.nan)
    return d


def seam_index(g_map, intron):
    """Transcript index of the LAST base before `intron` in this hypothesis (or None)."""
    lo, hi = intron
    jumps = np.nonzero(np.abs(np.diff(g_map)) != 1)[0]
    for i in jumps:
        a, b = int(g_map[i]), int(g_map[i + 1])
        if {a, b} == {lo - 1, hi} or (a == lo - 1 and b == hi) or (a == hi and b == lo - 1):
            return int(i)
    return None


def local_columns(tx1, g1, j1, tx2, g2, j2, strand, ctx=SEQ_CONTEXT):
    """Column alignment of the two spliced transcripts around their junctions, anchored on
    GENOMIC position. The window is the UNION genomic span of both seams +- ctx exonic bases,
    so a base missing from a row means that hypothesis SPLICES IT OUT, never that the window
    was cut short. Returns (columns, seam_col_1, seam_col_2); columns = [(gpos, b1, b2)] in
    transcript 5'->3' order, b = base or None."""
    i1 = seam_index(g1, j1); i2 = seam_index(g2, j2)
    if i1 is None or i2 is None:
        return None
    # exonic genomic positions kept by each hypothesis, ctx on each side of its own seam
    def side(gm, i):
        a = gm[max(0, i - ctx + 1): i + 1]
        b = gm[i + 1: i + 1 + ctx]
        return [int(x) for x in a], [int(x) for x in b]
    up1, dn1 = side(g1, i1); up2, dn2 = side(g2, i2)
    keep = set(up1) | set(dn1) | set(up2) | set(dn2)
    lo, hi = min(keep), max(keep)
    m1 = {int(g): tx1[k] for k, g in enumerate(g1) if lo <= int(g) <= hi}
    m2 = {int(g): tx2[k] for k, g in enumerate(g2) if lo <= int(g) <= hi}
    gs = sorted(set(m1) | set(m2), reverse=(strand == "-"))
    cols = [(g, m1.get(g), m2.get(g)) for g in gs]
    s1 = next((k for k, (g, _, _) in enumerate(cols) if g == int(g1[i1])), None)
    s2 = next((k for k, (g, _, _) in enumerate(cols) if g == int(g2[i2])), None)
    return cols, s1, s2


def runs_of(kmer1, kmer2, dnll, n_samp):
    """Maximal consecutive events sharing the same ordered k-mer pair."""
    out = []
    i = 0
    while i < len(kmer1):
        j = i
        while j + 1 < len(kmer1) and kmer1[j + 1] == kmer1[i] and kmer2[j + 1] == kmer2[i]:
            j += 1
        idx = np.arange(i, j + 1)
        out.append(dict(i0=i, i1=j, k1=kmer1[i], k2=kmer2[i], n=j - i + 1,
                        differs=kmer1[i] != kmer2[i], raw=float(np.nansum(dnll[idx])),
                        nsamp=float(np.nansum(n_samp[idx])), idx=idx))
        i = j + 1
    return out


# ----------------------------------------------------------------------------- per-read figure
def plot_read(args, rid, own_j, alt_j, own_is_true, tag, pair_type="artifact"):
    import krill
    import pyslow5

    with pysam.AlignmentFile(args.bam) as bam:
        rd = next(r for r in bam.fetch(until_eof=True)
                  if r.query_name == rid and not r.is_secondary and not r.is_supplementary)
    chain, ind = read_chain_and_indels(rd)
    strand = "-" if rd.is_reverse else "+"
    genome = pysam.FastaFile(args.genome)
    gseq = genome.fetch(rd.reference_name)
    alt_chain = substitute(chain, own_j, alt_j)
    tx_own, g_own = build_tx(gseq, rd.reference_start, rd.reference_end, chain, strand)
    tx_alt, g_alt = build_tx(gseq, rd.reference_start, rd.reference_end, alt_chain, strand)

    s5 = pyslow5.Open(args.blow5, "r"); rec = s5.get_read(rid, pA=False); s5.close()
    raw = rec["signal"].astype("float32")
    pA = (raw + float(rec["offset"])) * float(rec["range"]) / float(rec["digitisation"])
    meta = dict(digitisation=float(rec["digitisation"]), offset=float(rec["offset"]),
                range=float(rec["range"]), sample_rate=float(rec["sampling_rate"]))
    al = krill.Aligner(pore=args.pore, use_gpu=bool(args.gpu), batch_size=2, hmm_confidence=False)
    full = al.align([dict(read_id=rid, sequence=tx_own, start=0, signal=raw, **meta)])[0]
    assert full["status"] == 0, "full-read alignment failed"
    isc, ish = float(full["init_scale"]), float(full["init_shift"])
    ro, ra = al.align([dict(read_id=rid, sequence=tx_own, start=0, signal=raw,
                            init_scale=isc, init_shift=ish, **meta),
                       dict(read_id=rid, sequence=tx_alt, start=0, signal=raw,
                            init_scale=isc, init_shift=ish, **meta)])
    assert ro["status"] == 0 and ra["status"] == 0, "hypothesis alignment failed"

    # H1 = the hypothesis we score positively toward (TRUE for artifact pairs, read-own for controls)
    if pair_type == "artifact":
        (H1, g1, j1, lab1), (H2, g2, j2, lab2) = (
            ((ro, g_own, own_j, "TRUE (= read's own)") if own_is_true else (ra, g_alt, alt_j, "TRUE (annotation)")),
            ((ra, g_alt, alt_j, "FAKE (artifact)") if own_is_true else (ro, g_own, own_j, "FAKE (= read's own)")))
    else:
        (H1, g1, j1, lab1) = (ro, g_own, own_j, "H1: read's own genuine junction")
        (H2, g2, j2, lab2) = (ra, g_alt, alt_j, "H2: nearby genuine junction")
    e1, e2 = events_table(H1), events_table(H2)

    # ---- scored window: genomic +-FLANK around EVERY splice coordinate of either hypothesis
    coords = sorted({j1[0], j1[1], j2[0], j2[1]})
    wins = [(c - 8, c + 8) for c in coords]

    def in_any_win(e, gmap):
        idx = np.clip(e["pos"], 0, len(gmap) - 1)
        gp = np.where((e["pos"] >= 0) & (e["pos"] < len(gmap)), gmap[idx], -1)
        m = np.zeros(len(gp), bool)
        for lo, hi in wins:
            m |= (gp >= lo) & (gp <= hi)
        return m

    m1, m2 = in_any_win(e1, g1), in_any_win(e2, g2)
    if not m1.any() or not m2.any():
        print(f"skip {rid[:8]}: no events in window"); return
    lo_s = int(min(e1["lo"][m1].min(), e2["lo"][m2].min()))
    hi_s = int(max(e1["hi"][m1].max(), e2["hi"][m2].max()))

    sel1 = (e1["hi"] > lo_s) & (e1["lo"] < hi_s)
    sel2 = (e2["hi"] > lo_s) & (e2["lo"] < hi_s)
    key1 = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(e1["lo"], e1["hi"])) if sel1[i]}
    pairs = [(key1[(int(a), int(b))], i) for i, (a, b) in enumerate(zip(e2["lo"], e2["hi"]))
             if sel2[i] and (int(a), int(b)) in key1]
    pairs.sort(key=lambda p: e1["lo"][p[0]])
    i1 = np.array([p[0] for p in pairs]); i2 = np.array([p[1] for p in pairs])
    n_unmatched_1 = int(sel1.sum()) - len(i1); n_unmatched_2 = int(sel2.sum()) - len(i2)
    dnll = e2["nll"][i2] - e1["nll"][i1]                    # > 0 favours H1
    nsamp = (e1["hi"][i1] - e1["lo"][i1]).astype(float)
    tot1 = float(np.nansum(e1["nll"][i1])); tot2 = float(np.nansum(e2["nll"][i2]))
    score = tot2 - tot1
    R = runs_of(e1["kmer"][i1], e2["kmer"][i2], dnll, nsamp)
    diff_runs = [r for r in R if r["differs"]]
    absum = sum(abs(r["raw"]) for r in R) or 1.0
    dom = max((abs(r["raw"]) for r in R), default=0.0) / absum
    loo = [score - r["raw"] for r in R]
    loo_lo, loo_hi = (min(loo), max(loo)) if loo else (score, score)
    robust = (len(diff_runs) >= 2 and dom <= 0.5
              and (min(loo) > 0) == (max(loo) > 0) and (min(loo) > 0) == (score > 0))
    verdict = ("H1" if score > 0 else "H2") if robust else "ABSTAIN (run-dominated)"

    # ---- sequence alignment around the junction (genomic-column anchored)
    LC = local_columns(tx_own if H1 is ro else tx_alt, g1, j1,
                       tx_own if H2 is ro else tx_alt, g2, j2, strand)
    seq_ok = LC is not None
    only1 = only2 = ""
    if seq_ok:
        cols, seam1, seam2 = LC
        only1 = "".join(b1 for _, b1, b2 in cols if b1 is not None and b2 is None)
        only2 = "".join(b2 for _, b1, b2 in cols if b2 is not None and b1 is None)
        cdna1 = "".join(b1 for _, b1, _ in cols if b1 is not None)
        cdna2 = "".join(b2 for _, _, b2 in cols if b2 is not None)

    # ---------------------------------------------------------------- figure
    fig = plt.figure(figsize=(15.5, 14.2))
    gs = fig.add_gridspec(7, 5, height_ratios=[1.05, 1.30, 1.15, 1.5, 1.1, 1.1, 1.35],
                          hspace=0.66, wspace=0.55, left=0.055, right=0.985, top=0.90, bottom=0.045)

    d1, a1 = donor_acceptor(j1, strand); d2, a2 = donor_acceptor(j2, strand)
    edit_txt = "n/a"
    if seq_ok:
        if cdna1 == cdna2:
            edit_txt = ("IDENTICAL cDNA in this window (intron slide) -- "
                        "NO sequence-identifiable contrast")
        elif only1 and not only2:
            edit_txt = f"H1 has {len(only1)} exonic base(s) H2 splices out: {only1}"
        elif only2 and not only1:
            edit_txt = f"H2 has {len(only2)} exonic base(s) H1 splices out: {only2}"
        else:
            edit_txt = f"H1-only bases {only1 or '-'} ; H2-only bases {only2 or '-'}"

    head = (
        f"SIRV_Set4 spike-in  |  {args.sample}  |  {rd.reference_name}  |  read {rid[:8]}   "
        f"pair={pair_type}   strand {strand}   "
        f"read CIGAR carries: {'H1' if (own_is_true and pair_type=='artifact') or pair_type!='artifact' else 'H2'}   "
        f"junction indel-clean: {'yes' if clean_at(ind, *own_j) else 'no'}   {len(chain)} introns\n"
        f"H1 genomic intron [{j1[0]}, {j1[1]})   H2 genomic intron [{j2[0]}, {j2[1]})   |   "
        f"donor shift {d2-d1:+d} bp, acceptor shift {a2-a1:+d} bp   |   {edit_txt}\n"
        f"krill {args.krill_ver} pore={args.pore} HMM=off   fixed init_scale={isc:.5f} init_shift={ish:.3f} "
        f"(both hypotheses)   window [{lo_s}, {hi_s}) = {hi_s-lo_s} samples\n"
        f"paired events {len(i1)} (unmatched H1 {n_unmatched_1}, H2 {n_unmatched_2})   "
        f"differing-k-mer events {int(np.sum(e1['kmer'][i1] != e2['kmer'][i2]))} in {len(diff_runs)} distinct runs\n"
        f"per-event NLL sum:  H1 = {tot1:.1f}   H2 = {tot2:.1f}   "
        f"score = NLL(H2) - NLL(H1) = {score:+.1f}  (>0 favours H1)   "
        f"dominant-run fraction {dom:.0%}   leave-one-run-out [{loo_lo:+.1f}, {loo_hi:+.1f}]   "
        f"verdict: {verdict}")
    fig.text(0.06, 0.995, head, va="top", ha="left", fontsize=8.6, family="monospace", color=C_TXT,
             linespacing=1.55)

    # ---- Panel A1 / A2: genomic donor and acceptor windows, transcript 5'->3'
    for col, (name, c1, c2) in enumerate((("donor", d1, d2), ("acceptor", a1, a2))):
        ax = fig.add_subplot(gs[0, 0:2] if col == 0 else gs[0, 3:5])
        lo = min(c1, c2) - FLANK_GENOMIC; hi = max(c1, c2) + FLANK_GENOMIC
        # exon side: for donor the exon is upstream (transcript-wise), for acceptor downstream
        for y, (c, col_h, lab) in enumerate(((c2, C_H2, lab2), (c1, C_H1, lab1))):
            if name == "donor":
                ex = (lo, c) if strand == "+" else (c, hi)
                intr = (c, hi) if strand == "+" else (lo, c)
            else:
                ex = (c, hi) if strand == "+" else (lo, c)
                intr = (lo, c) if strand == "+" else (c, hi)
            ax.add_patch(Rectangle((ex[0], y - 0.17), ex[1] - ex[0], 0.34, color=col_h, alpha=0.85, lw=0))
            ax.plot(intr, [y, y], color=col_h, lw=1.0, ls="--")
            ax.plot([c, c], [y - 0.32, y + 0.32], color=col_h, lw=1.6)
            ax.plot([c], [y], marker="o" if col_h == C_H1 else "D", ms=6,
                    mfc=col_h if col_h == C_H1 else "none", mec=col_h, mew=1.4)
        ax.set_ylim(-0.75, 1.75); ax.set_yticks([0, 1]); ax.set_yticklabels(["H2", "H1"], fontsize=8)
        if strand == "-":
            ax.set_xlim(hi, lo)
        else:
            ax.set_xlim(lo, hi)
        same = c1 == c2
        ax.set_title(f"A{col+1}. biological {name}: H1 {c1}  H2 {c2}   "
                     + ("SAME coordinate (not part of the scored difference)" if same else f"shift {c2-c1:+d} bp"),
                     loc="left", color=C_TXT2 if same else C_TXT)
        ax.set_xlabel(f"{rd.reference_name} genomic position   (left -> right = transcript 5' -> 3')", fontsize=8)
        for lo_w, hi_w in wins:
            if lo_w < hi and hi_w > lo:
                ax.axvspan(lo_w, hi_w, color=C_GRID, zorder=0)
    fig.text(0.5, 0.905, f"intron  {abs(j1[1]-j1[0])} bp  (not to scale)", ha="center", va="bottom",
             fontsize=8, color=C_TXT2, style="italic")

    # ---- Panel B: aligned spliced-transcript sequences
    axB = fig.add_subplot(gs[1, :]); axB.axis("off")
    axB.set_title("B. spliced transcript sequences around the junction, pairwise aligned "
                  "(transcript 5'->3'; '|' = splice seam; coloured = present in only one hypothesis)",
                  loc="left", color=C_TXT)
    if not seq_ok:
        axB.text(0.0, 0.5, "sequence context unavailable (junction not found in transcript map)",
                 fontsize=9, color=C_TXT2)
    else:
        n = len(cols)
        x0, dx = 0.135, min(0.0135, 0.855 / max(n, 1))
        fs = 8.5 if dx > 0.011 else 7.0
        for row, (ytop, col_h, lab, seam) in enumerate(((0.80, C_H1, lab1, seam1),
                                                        (0.34, C_H2, lab2, seam2))):
            axB.text(0.0, ytop, lab, fontsize=8.5, color=col_h, va="center", family="monospace")
            axB.text(x0 - 0.020, ytop, "5'", fontsize=7.5, color=C_TXT2, va="center", family="monospace")
            for k, (gp, b1, b2) in enumerate(cols):
                b = b1 if row == 0 else b2
                other = b2 if row == 0 else b1
                x = x0 + k * dx
                if b is None:
                    axB.add_patch(Rectangle((x - dx * 0.45, ytop - 0.115), dx * 0.9, 0.23,
                                            color=C_PALE, lw=0, transform=axB.transAxes))
                    axB.text(x, ytop, "-", fontsize=fs, color=C_TXT2, ha="center", va="center",
                             family="monospace")
                elif other is None:
                    axB.add_patch(Rectangle((x - dx * 0.45, ytop - 0.115), dx * 0.9, 0.23,
                                            color=col_h, lw=0, transform=axB.transAxes))
                    axB.text(x, ytop, b, fontsize=fs, color="white", ha="center", va="center",
                             family="monospace", fontweight="bold")
                else:
                    axB.text(x, ytop, b, fontsize=fs, color=C_TXT, ha="center", va="center",
                             family="monospace")
                if seam is not None and k == seam:
                    axB.plot([x + dx * 0.5] * 2, [ytop - 0.165, ytop + 0.165], color=col_h, lw=2.0,
                             transform=axB.transAxes)
            axB.text(x0 + n * dx + 0.004, ytop, "3'", fontsize=7.5, color=C_TXT2, va="center",
                     family="monospace")
        for k, (gp, b1, b2) in enumerate(cols):
            axB.text(x0 + k * dx, 0.57, "|" if (b1 is not None and b2 is not None and b1 == b2) else " ",
                     fontsize=fs, color=C_TXT2, ha="center", va="center", family="monospace")
        # genomic coordinate ticks every 5 columns
        for k in range(0, n, 5):
            axB.text(x0 + k * dx, 0.10, str(cols[k][0]), fontsize=6.2, color=C_TXT2, ha="center",
                     va="center", rotation=90, family="monospace")
        axB.text(0.0, 0.10, "genomic pos", fontsize=7, color=C_TXT2, va="center")
        axB.text(0.0, -0.06, edit_txt + ("   |   shown as the reverse complement of the genome (strand -)"
                                         if strand == "-" else "")
                 + "   |   thick bar = splice seam of that hypothesis",
                 fontsize=8.5, color=C_TXT, va="center")
        axB.set_xlim(0, 1); axB.set_ylim(-0.12, 1)

    # ---- shared x-limits for C / D / E
    x0s, x1s = lo_s - 60, hi_s + 60

    # ---- Panel C: event -> k-mer assignment strip with runs
    axC = fig.add_subplot(gs[2, :])
    for y, (e, idx, col_h, lab) in enumerate(((e2, i2, C_H2, lab2), (e1, i1, C_H1, lab1))):
        for k, i in enumerate(idx):
            differs = e1["kmer"][i1[k]] != e2["kmer"][i2[k]]
            axC.add_patch(Rectangle((e["lo"][i], y - 0.30), e["hi"][i] - e["lo"][i], 0.60,
                                    facecolor=col_h if differs else "#efefec",
                                    edgecolor=col_h if differs else "#d8d7d2",
                                    alpha=0.9 if differs else 1.0, lw=0.6))
    axC.set_ylim(-0.75, 3.5); axC.set_yticks([0, 1]); axC.set_yticklabels(["H2", "H1"], fontsize=8)
    axC.set_xlim(x0s, x1s)
    top_runs = sorted(diff_runs, key=lambda r: -abs(r["raw"]))[:3]
    for r in diff_runs:
        xa = float(e1["lo"][i1[r["i0"]]]); xb = float(e1["hi"][i1[r["i1"]]])
        axC.add_patch(Rectangle((xa, -0.55), max(xb - xa, 4), 2.05, fill=False,
                                edgecolor=C_TXT if r in top_runs else "#b9b8b3",
                                lw=1.1 if r in top_runs else 0.7, ls=":"))
        if r in top_runs:
            k = top_runs.index(r)
            axC.annotate(f"R{diff_runs.index(r)+1}  {r['k1']}/{r['k2']}  {r['n']} ev, "
                         f"{int(r['nsamp'])} smp, raw {r['raw']:+.1f}",
                         ((xa + xb) / 2, 1.55), xytext=(0, 4 + 14 * k), textcoords="offset points",
                         ha="center", fontsize=7.0, color=C_TXT, family="monospace",
                         arrowprops=dict(arrowstyle="-", color="#b9b8b3", lw=0.6))
    axC.text(0.0, 1.20, "C. event -> assigned k-mer (same event boundaries under both hypotheses; "
             "coloured = hypotheses assign DIFFERENT k-mers; dotted box = one run of the same k-mer pair)",
             transform=axC.transAxes, ha="left", va="bottom", fontsize=9, color=C_TXT)
    axC.set_xticklabels([])

    # ---- Panel D1/D2/D3: observed then per-hypothesis model tracks
    ylim = (min(np.nanmin(e1["lev"][i1]), np.nanmin(e1["mm"][i1]), np.nanmin(e2["mm"][i2])) - 8,
            max(np.nanmax(e1["lev"][i1]), np.nanmax(e1["mm"][i1]), np.nanmax(e2["mm"][i2])) + 8)
    xs = np.arange(max(0, x0s), min(len(pA), x1s))
    axD1 = fig.add_subplot(gs[3, :])
    axD1.plot(xs, pA[xs], color="#c2c1bc", lw=0.55)
    for k, i in enumerate(i1):
        axD1.plot([e1["lo"][i], e1["hi"][i]], [e1["lev"][i]] * 2, color=C_OBS, lw=1.8)
    axD1.set_ylabel("D1 observed\ncurrent (pA)", fontsize=8)
    axD1.set_title("D. same raw samples under each hypothesis (identical pA limits; black = observed event mean)",
                   loc="left", color=C_TXT)
    axes_D = [axD1]
    for row, (e, idx, col_h, lab) in enumerate(((e1, i1, C_H1, lab1), (e2, i2, C_H2, lab2))):
        ax = fig.add_subplot(gs[4 + row, :])
        for k, i in enumerate(idx):
            ax.fill_between([e["lo"][i], e["hi"][i]], e["mm"][i] - e["sig"][i], e["mm"][i] + e["sig"][i],
                            color=col_h, alpha=0.16, lw=0)
            ax.plot([e["lo"][i], e["hi"][i]], [e["mm"][i]] * 2, color=col_h, lw=1.8)
            ax.plot([e1["lo"][i1[k]], e1["hi"][i1[k]]], [e1["lev"][i1[k]]] * 2, color=C_OBS, lw=0.9, alpha=0.75)
        ax.set_ylabel(f"D{row+2} model\n{lab.split('(')[0].strip()} (pA)", fontsize=8, color=col_h)
        axes_D.append(ax)
    for ax in axes_D:
        ax.set_xlim(x0s, x1s); ax.set_ylim(*ylim); ax.set_xticklabels([])
        ax.axvspan(lo_s, hi_s, color=C_PALE, zorder=0)

    # ---- Panel E: per-event dNLL
    axE = fig.add_subplot(gs[6, 0:3])
    xe = [(e1["lo"][i] + e1["hi"][i]) / 2 for i in i1]
    we = [max(e1["hi"][i] - e1["lo"][i] - 1, 1) for i in i1]
    differs = e1["kmer"][i1] != e2["kmer"][i2]
    cols = [(C_H1 if v > 0 else C_H2) if d else "#dcdbd6" for v, d in zip(dnll, differs)]
    axE.bar(xe, np.nan_to_num(dnll), width=we, color=cols, align="center", lw=0)
    axE.axhline(0, color=C_TXT2, lw=0.8)
    axE.set_xlim(x0s, x1s)
    axE.set_ylabel("E  per-event\nNLL(H2) - NLL(H1)", fontsize=8)
    axE.set_xlabel("raw signal sample index", fontsize=8)
    axE.set_title(f"E. per-event evidence (metric: per-event mean-current NLL sum). "
                  f"blue = favours H1, orange = favours H2, grey = same k-mer",
                  loc="left", color=C_TXT, fontsize=8.5)

    # ---- Panel F: run-level dependence audit
    axF = fig.add_subplot(gs[6, 3:5])
    show = sorted(R, key=lambda r: -abs(r["raw"]))[:8]
    lbl = []
    for r in show:
        tag_r = (f"R{diff_runs.index(r)+1} " if r["differs"] else "= ")
        lbl.append(f"{tag_r}{r['k1']}/{r['k2']}  n={r['n']}")
    y = np.arange(len(show))[::-1]
    axF.barh(y, [r["raw"] for r in show], color=[C_H1 if r["raw"] > 0 else C_H2 for r in show], height=0.62)
    axF.set_yticks(y); axF.set_yticklabels(lbl, fontsize=7, family="monospace")
    axF.axvline(0, color=C_TXT2, lw=0.8)
    axF.set_xlabel("raw summed contribution of the run", fontsize=8)
    axF.set_title(f"F. dependence audit: {len(diff_runs)} informative run(s), dominant run = {dom:.0%} of |score|,\n"
                  f"leave-one-run-out score in [{loo_lo:+.1f}, {loo_hi:+.1f}]  ->  {verdict}",
                  loc="left", color=C_TXT if robust else "#d03b3b", fontsize=8.5)

    out = Path(args.out) / f"fig_read_{tag}_{rid[:8]}.png"
    fig.savefig(out, dpi=145, bbox_inches="tight"); plt.close(fig)
    print(f"wrote {out}  score={score:+.1f} dom={dom:.2f} runs={len(diff_runs)} verdict={verdict}")


# ----------------------------------------------------------------------------- distribution figure
def plot_distribution(tsv, out):
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))

    def f(r, k):
        return float(r[k]) if r.get(k) not in (None, "") else None

    bins_shift = ((1, 3, "shift 1-3 bp"), (4, 6, "shift 4-6 bp"), (7, 99, "shift >= 7 bp"))
    fig, axes = plt.subplots(2, 3, figsize=(13, 6.4), sharex=True)
    edges = np.linspace(-40, 80, 61)
    clip = lambda v: np.clip(np.asarray(v, float), -40, 80)
    for j, (lo, hi, lab) in enumerate(bins_shift):
        ax = axes[0, j]
        fk = [f(r, "ev_true") for r in rows if r["pair_type"] == "artifact" and r["side"] == "FAKE"
              and r["ev_true"] and lo <= int(r["shift_max"]) <= hi]
        tr = [f(r, "ev_true") for r in rows if r["pair_type"] == "artifact" and r["side"] == "TRUE"
              and r["ev_true"] and lo <= int(r["shift_max"]) <= hi]
        ax.hist(clip(fk), bins=edges, color=C_H2, alpha=0.75,
                label=f"CIGAR says FAKE  (n={len(fk)}, {100*np.mean(np.asarray(fk)>0):.0f}% favour TRUE)")
        ax.hist(clip(tr), bins=edges, color=C_H1, alpha=0.55,
                label=f"CIGAR says TRUE  (n={len(tr)}, {100*np.mean(np.asarray(tr)>0):.0f}% favour TRUE)")
        ax.axvline(0, color=C_TXT2, lw=0.8)
        ax.set_title(f"artifact pair, {lab}", loc="left", fontsize=9, color=C_TXT)
        ax.legend(fontsize=7.5, frameon=False, loc="upper right")
        ax = axes[1, j]
        gg = [f(r, "ev_contrast") for r in rows if r["pair_type"] == "genuine" and r["ev_contrast"]
              and lo <= int(r["shift_max"]) <= hi]
        if gg:
            ax.hist(clip(gg), bins=edges, color="#1baf7a", alpha=0.85,
                    label=f"genuine pair  (n={len(gg)}, {100*np.mean(np.asarray(gg)>0):.0f}% keep own)")
            ax.legend(fontsize=7.5, frameon=False, loc="upper right")
        ax.axvline(0, color=C_TXT2, lw=0.8)
        ax.set_title(f"control: two genuine junctions, {lab}", loc="left", fontsize=9, color=C_TXT)
    axes[0, 0].set_ylabel("reads"); axes[1, 0].set_ylabel("reads")
    fig.supxlabel("signal score  =  NLL(the other junction)  -  NLL(this read's scored junction)        "
                  "positive  ->  the raw signal fits the TRUE / own junction better",
                  fontsize=9.5, color=C_TXT2, y=-0.01)
    fig.suptitle("SIRV_Set4 spike-in reads (H9 dRNA): one read, two junction hypotheses, "
                 "same raw signal, same initial scaling (krill 0.1.6)",
                 fontsize=11, color=C_TXT, y=1.0)
    fig.tight_layout()
    o = Path(out) / f"fig_distribution_{Path(tsv).stem}.png"
    fig.savefig(o, dpi=150, bbox_inches="tight"); plt.close(fig)
    print("wrote", o)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", default=""); ap.add_argument("--bam"); ap.add_argument("--blow5")
    ap.add_argument("--genome"); ap.add_argument("--pore", default="rna002")
    ap.add_argument("--krill-ver", default="0.1.6+cu122")
    ap.add_argument("--gpu", type=int, default=1); ap.add_argument("--out", required=True)
    ap.add_argument("--read", action="append", default=[],
                    help="rid:own_d-own_a:alt_d-alt_a:own_is_true(0/1):tag[:pair_type]")
    ap.add_argument("--tsv", action="append", default=[])
    a = ap.parse_args()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    for spec in a.read:
        parts = spec.split(":")
        rid, own, alt, oit, tag = parts[:5]
        ptype = parts[5] if len(parts) > 5 else "artifact"
        od, oa = map(int, own.split("-")); ad, aa = map(int, alt.split("-"))
        plot_read(a, rid, (od, oa), (ad, aa), oit == "1", tag, ptype)
    for t in a.tsv:
        plot_distribution(t, a.out)


if __name__ == "__main__":
    main()
