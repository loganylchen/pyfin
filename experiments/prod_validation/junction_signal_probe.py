#!/usr/bin/env python3
"""Do krill eventalign events separate TRUE from FALSE junctions on SIRV?

Read-only, evaluation-only probe (truth GTF is used ONLY to label). SIRV_Set4
is synthetic, so a read junction absent from the catalog is an alignment
artifact; a catalog junction within +-NEAR bp of it is the truth it was
displaced from.

For every (artifact junction J_fake -> nearest true junction J_true) pair and,
as a control, every genuine (J_true1, J_true2) nearby pair, take reads whose
own CIGAR chain carries one side, build TWO full transcript hypotheses from
the genome -- the read's own chain, and the same chain with that single
junction replaced by the partner -- and eventalign the read's raw signal
against both with krill 0.1.6:

  run 1: own hypothesis alone -> init_scale / init_shift of the whole read
  run 2: [own, alt] with that SAME fixed initial scaling

Over the raw-signal window covering the junction neighbourhood (genomic
+-FLANK bp around every coordinate that differs, mapped through each
hypothesis; the union of both hypotheses' sample ranges, so both score the
same physical signal) we compute per hypothesis:

  ev_sum  = sum_e [0.5 z^2 + log sigma]                 (old per-EVENT metric)
  ps_sum  = sum_e n_e [log sigma + ((m-mu)^2 + s^2)/(2 sigma^2)]   (per-SAMPLE)

with m = event_level_mean, s = event_stdv (pA), mu = model_mean (pA), sigma =
model_stdv * sqrt(var) (verified: reproduces krill's standardized_level),
n_e = end_idx - start_idx (samples in the event).

score = cost(alt) - cost(own): > 0 means the signal favours the read's OWN
junction. For artifact pairs the correct answer is always J_true, so we also
report score_true = cost(fake hyp) - cost(true hyp) (> 0 = signal is right).

Outputs per-read TSV + summary JSON per sample and pooled.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
from junction_clean_probe import truth_junctions  # noqa: E402

NEAR = 15
FLANK = 8
K_RNA002 = 5
CENTER = (K_RNA002 - 1) // 2
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_chain_and_indels(rd):
    pos = rd.reference_start
    chain, ind = [], []
    for op, ln in rd.cigartuples:
        if op in (0, 7, 8):
            pos += ln
        elif op == 2:
            ind.append((pos, pos + ln))
            pos += ln
        elif op == 1:
            ind.append((pos, pos))
        elif op == 3:
            chain.append((pos, pos + ln))
            pos += ln
    return tuple(chain), tuple(ind)


def clean_at(ind, d, a, flank=10):
    return all(not ((lo - flank <= d <= hi + flank) or (lo - flank <= a <= hi + flank))
               for lo, hi in ind)


def exon_blocks(start, end, chain):
    blocks, cur = [], start
    for d, a in chain:
        blocks.append((cur, d))
        cur = a
    blocks.append((cur, end))
    return blocks


def valid_chain(start, end, chain):
    prev = start
    for d, a in chain:
        if not (prev < d < a < end):
            return False
        prev = a
    return True


def build_tx(genome_seq: str, start: int, end: int, chain, strand: str):
    """Spliced transcript sequence + tx_index -> genomic position array (transcript 5'->3')."""
    blocks = exon_blocks(start, end, chain)
    seq = "".join(genome_seq[s:e] for s, e in blocks)
    g = np.concatenate([np.arange(s, e) for s, e in blocks]).astype(np.int64)
    if strand == "-":
        seq = revcomp(seq)
        g = g[::-1].copy()
    return seq, g


def substitute(chain, j_old, j_new):
    return tuple(j_new if j == j_old else j for j in chain)


def hkey(s: str) -> str:
    return hashlib.sha1(s.encode()).hexdigest()


def event_costs(res, g_map, gwin_lo, gwin_hi, use_hmm):
    """Per-event arrays restricted to events whose central base maps inside [gwin_lo, gwin_hi].
    Returns dict with arrays: raw_lo, raw_hi, ev (per-event nll), ps (per-sample nll), n."""
    if res is None or res.get("status", -1) != 0 or len(res["position"]) == 0:
        return None
    pos = np.asarray(res["position"], np.int64)
    keep = np.ones(len(pos), bool)
    if use_hmm and "hmm_state" in res:
        hs = np.asarray(res["hmm_state"])
        keep &= (hs != "K") & (hs != "B")
    mm = np.asarray(res["model_mean"], float)
    ms = np.asarray(res["model_stdv"], float)
    lev = np.asarray(res["event_level_mean"], float)
    sd = np.asarray(res["event_stdv"], float)
    s0 = np.asarray(res["start_idx"], np.int64)
    s1 = np.asarray(res["end_idx"], np.int64)
    keep &= np.isfinite(mm) & np.isfinite(ms) & (ms > 0) & np.isfinite(lev) & np.isfinite(sd)
    idx = np.clip(pos, 0, len(g_map) - 1)
    gpos = g_map[idx]
    keep &= (pos >= 0) & (pos < len(g_map)) & (gpos >= gwin_lo) & (gpos <= gwin_hi)
    if not keep.any():
        return None
    var = float(res.get("var", 1.0)) or 1.0
    sig = ms[keep] * math.sqrt(max(var, 1e-12))
    m, mu, s = lev[keep], mm[keep], sd[keep]
    n = np.abs(s1[keep] - s0[keep]).astype(float)
    z = (m - mu) / sig
    ev = 0.5 * z * z + np.log(sig)
    ps = n * (np.log(sig) + ((m - mu) ** 2 + s * s) / (2.0 * sig * sig))
    return dict(raw_lo=np.minimum(s0[keep], s1[keep]), raw_hi=np.maximum(s0[keep], s1[keep]),
                ev=ev, ps=ps, n=n)


def window_costs(res, g_map, gwin_lo, gwin_hi, use_hmm):
    """All events of `res` (not only in-window) as arrays + the in-window sample range."""
    c = event_costs(res, g_map, gwin_lo, gwin_hi, use_hmm)
    if c is None:
        return None, None
    return c, (int(c["raw_lo"].min()), int(c["raw_hi"].max()))


def score_over_samples(res, g_map, lo_s, hi_s, use_hmm):
    """Sum ev / ps over ALL events of `res` overlapping raw sample window [lo_s, hi_s]."""
    if res is None or res.get("status", -1) != 0 or len(res["position"]) == 0:
        return None
    mm = np.asarray(res["model_mean"], float)
    ms = np.asarray(res["model_stdv"], float)
    lev = np.asarray(res["event_level_mean"], float)
    sd = np.asarray(res["event_stdv"], float)
    s0 = np.asarray(res["start_idx"], np.int64)
    s1 = np.asarray(res["end_idx"], np.int64)
    lo = np.minimum(s0, s1)
    hi = np.maximum(s0, s1)
    keep = (hi > lo_s) & (lo < hi_s)
    if use_hmm and "hmm_state" in res:
        hs = np.asarray(res["hmm_state"])
        keep &= (hs != "K") & (hs != "B")
    keep &= np.isfinite(mm) & np.isfinite(ms) & (ms > 0) & np.isfinite(lev) & np.isfinite(sd)
    if not keep.any():
        return None
    var = float(res.get("var", 1.0)) or 1.0
    sig = ms[keep] * math.sqrt(max(var, 1e-12))
    m, mu, s = lev[keep], mm[keep], sd[keep]
    n = (hi[keep] - lo[keep]).astype(float)
    z = (m - mu) / sig
    ev = float(np.sum(0.5 * z * z + np.log(sig)))
    ps = float(np.sum(n * (np.log(sig) + ((m - mu) ** 2 + s * s) / (2.0 * sig * sig))))
    return dict(ev=ev, ps=ps, n_ev=int(keep.sum()), n_samp=float(n.sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fq", required=True)
    ap.add_argument("--blow5", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pore", default="rna002")
    ap.add_argument("--hmm", type=int, default=0)
    ap.add_argument("--cap", type=int, default=30, help="max reads per (pair, side)")
    ap.add_argument("--gpu", type=int, default=1)
    ap.add_argument("--batch", type=int, default=64)
    args = ap.parse_args()

    import krill
    import pyslow5

    truth = truth_junctions(Path(args.truth))
    truth_by_chrom = collections.defaultdict(list)
    for (c, d, a) in truth:
        truth_by_chrom[c].append((d, a))
    genome = pysam.FastaFile(args.genome)
    gcache = {}

    def gseq(chrom):
        if chrom not in gcache:
            gcache[chrom] = genome.fetch(chrom)
        return gcache[chrom]

    # ---- pass 1: reads, chains, junction classes -------------------------------------------
    reads = {}      # rid -> dict(chrom, strand, start, end, chain, ind)
    carriers = collections.defaultdict(list)   # (chrom, strand, d, a) -> [rid]
    with pysam.AlignmentFile(args.bam) as bam:
        for rd in bam.fetch(until_eof=True):
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary or rd.mapping_quality < 1:
                continue
            if not rd.cigartuples:
                continue
            chain, ind = read_chain_and_indels(rd)
            if not chain:
                continue
            rid = rd.query_name
            strand = "-" if rd.is_reverse else "+"
            reads[rid] = dict(chrom=rd.reference_name, strand=strand, start=rd.reference_start,
                              end=rd.reference_end, chain=chain, ind=ind)
            for j in chain:
                carriers[(rd.reference_name, strand, j[0], j[1])].append(rid)

    def cls_of(chrom, d, a):
        if (chrom, d, a) in truth:
            return "TRUE", None
        best = None
        for td, ta in truth_by_chrom.get(chrom, ()):
            if abs(td - d) <= NEAR and abs(ta - a) <= NEAR:
                dist = abs(td - d) + abs(ta - a)
                if best is None or dist < best[0]:
                    best = (dist, td, ta)
        return ("NEAR_FAKE", (best[1], best[2])) if best else ("FAR_FAKE", None)

    # ---- pairs -------------------------------------------------------------------------------
    pairs = []   # dict(type, chrom, strand, j1, j2)  j1 = own side A, j2 = partner
    seen_true_pairs = set()
    for (chrom, strand, d, a), rids in sorted(carriers.items()):
        cls, near = cls_of(chrom, d, a)
        if cls == "NEAR_FAKE":
            pairs.append(dict(type="artifact", chrom=chrom, strand=strand,
                              fake=(d, a), true=near))
        elif cls == "TRUE":
            for td, ta in truth_by_chrom[chrom]:
                if (td, ta) == (d, a):
                    continue
                if abs(td - d) <= NEAR and abs(ta - a) <= NEAR and (chrom, strand, td, ta) in carriers:
                    key = (chrom, strand) + tuple(sorted([(d, a), (td, ta)]))
                    if key in seen_true_pairs:
                        continue
                    seen_true_pairs.add(key)
                    pairs.append(dict(type="genuine", chrom=chrom, strand=strand,
                                      fake=None, true=None, t1=(d, a), t2=(td, ta)))
    print(f"[{args.sample}] reads={len(reads)} artifact_pairs={sum(p['type']=='artifact' for p in pairs)} "
          f"genuine_pairs={sum(p['type']=='genuine' for p in pairs)}", file=sys.stderr)

    # ---- tasks: (pair_id, side_label, rid, j_own, j_alt) -----------------------------------
    tasks = []
    for pid, p in enumerate(pairs):
        if p["type"] == "artifact":
            sides = [("FAKE", p["fake"], p["true"]), ("TRUE", p["true"], p["fake"])]
        else:
            sides = [("T1", p["t1"], p["t2"]), ("T2", p["t2"], p["t1"])]
        for label, j_own, j_alt in sides:
            rids = carriers.get((p["chrom"], p["strand"]) + j_own, [])
            rids = sorted(set(rids), key=hkey)[: args.cap]
            for rid in rids:
                r = reads[rid]
                alt_chain = substitute(r["chain"], j_own, j_alt)
                if alt_chain == r["chain"] or not valid_chain(r["start"], r["end"], alt_chain):
                    continue
                tasks.append(dict(pid=pid, side=label, rid=rid, j_own=j_own, j_alt=j_alt))
    print(f"[{args.sample}] tasks={len(tasks)}", file=sys.stderr)

    # ---- read sequences + signal ---------------------------------------------------------------
    need = {t["rid"] for t in tasks}
    seqs = {}
    with pysam.FastxFile(args.fq) as fh:
        for rec in fh:
            if rec.name in need:
                seqs[rec.name] = rec.sequence
    s5 = pyslow5.Open(args.blow5, "r")
    sig = {}
    for rid in need:
        rec = s5.get_read(rid, pA=False)
        if rec is None:
            continue
        sig[rid] = dict(signal=rec["signal"].astype("float32"), digitisation=float(rec["digitisation"]),
                        offset=float(rec["offset"]), range=float(rec["range"]),
                        sample_rate=float(rec["sampling_rate"]))
    s5.close()
    tasks = [t for t in tasks if t["rid"] in seqs and t["rid"] in sig]
    print(f"[{args.sample}] tasks with seq+signal={len(tasks)}", file=sys.stderr)

    aligner = krill.Aligner(pore=args.pore, use_gpu=bool(args.gpu), batch_size=args.batch,
                            hmm_confidence=bool(args.hmm))

    # ---- run 1: own hypothesis -> init scaling (batched, one per distinct read) --------------
    own_init = {}
    rid_list = sorted({t["rid"] for t in tasks})
    for i in range(0, len(rid_list), args.batch):
        chunk = rid_list[i:i + args.batch]
        items = []
        for rid in chunk:
            r = reads[rid]
            tx, _ = build_tx(gseq(r["chrom"]), r["start"], r["end"], r["chain"], r["strand"])
            items.append(dict(read_id=rid, sequence=tx, start=0, **sig[rid]))
        for rid, res in zip(chunk, aligner.align(items)):
            if res.get("status", -1) == 0:
                own_init[rid] = (float(res["init_scale"]), float(res["init_shift"]))
    print(f"[{args.sample}] run1 ok={len(own_init)}/{len(rid_list)}", file=sys.stderr)

    # ---- run 2: [own, alt] with fixed init scaling -------------------------------------------
    rows = []
    tasks = [t for t in tasks if t["rid"] in own_init]
    for i in range(0, len(tasks), max(1, args.batch // 2)):
        chunk = tasks[i:i + max(1, args.batch // 2)]
        items, meta = [], []
        for t in chunk:
            r = reads[t["rid"]]
            g = gseq(r["chrom"])
            tx_own, g_own = build_tx(g, r["start"], r["end"], r["chain"], r["strand"])
            alt_chain = substitute(r["chain"], t["j_own"], t["j_alt"])
            tx_alt, g_alt = build_tx(g, r["start"], r["end"], alt_chain, r["strand"])
            isc, ish = own_init[t["rid"]]
            base = dict(read_id=t["rid"], start=0, init_scale=isc, init_shift=ish, **sig[t["rid"]])
            items.append(dict(base, sequence=tx_own))
            items.append(dict(base, sequence=tx_alt))
            meta.append((t, g_own, g_alt))
        results = aligner.align(items)
        for k, (t, g_own, g_alt) in enumerate(meta):
            res_own, res_alt = results[2 * k], results[2 * k + 1]
            r = reads[t["rid"]]
            d1, a1 = t["j_own"]
            d2, a2 = t["j_alt"]
            # genomic window: every differing coordinate +- FLANK (donor side and/or acceptor side)
            wins = []
            if d1 != d2:
                wins.append((min(d1, d2) - FLANK, max(d1, d2) + FLANK))
            if a1 != a2:
                wins.append((min(a1, a2) - FLANK, max(a1, a2) + FLANK))
            row = dict(sample=args.sample, pair=t["pid"], pair_type=pairs[t["pid"]]["type"],
                       side=t["side"], rid=t["rid"], chrom=r["chrom"], strand=r["strand"],
                       own=f"{d1}-{a1}", alt=f"{d2}-{a2}", shift_d=d2 - d1, shift_a=a2 - a1,
                       shift_max=max(abs(d2 - d1), abs(a2 - a1)),
                       own_clean=int(clean_at(r["ind"], d1, a1)),
                       n_introns=len(r["chain"]), status_own=res_own.get("status", -1),
                       status_alt=res_alt.get("status", -1),
                       init_scale=own_init[t["rid"]][0], init_shift=own_init[t["rid"]][1])
            ok = row["status_own"] == 0 and row["status_alt"] == 0
            ev_c = ps_c = None
            n_ev_o = n_ev_a = 0
            n_samp = 0.0
            if ok:
                # sample window = union over both hypotheses of in-window event sample ranges
                lo_s, hi_s = None, None
                for res, gm in ((res_own, g_own), (res_alt, g_alt)):
                    for wl, wh in wins:
                        c, rng = window_costs(res, gm, wl, wh, bool(args.hmm))
                        if rng is None:
                            continue
                        lo_s = rng[0] if lo_s is None else min(lo_s, rng[0])
                        hi_s = rng[1] if hi_s is None else max(hi_s, rng[1])
                if lo_s is not None and hi_s is not None and hi_s > lo_s:
                    so = score_over_samples(res_own, g_own, lo_s, hi_s, bool(args.hmm))
                    sa = score_over_samples(res_alt, g_alt, lo_s, hi_s, bool(args.hmm))
                    if so and sa:
                        ev_c = sa["ev"] - so["ev"]
                        ps_c = sa["ps"] - so["ps"]
                        n_ev_o, n_ev_a = so["n_ev"], sa["n_ev"]
                        n_samp = float(hi_s - lo_s)
            row.update(n_ev_own=n_ev_o, n_ev_alt=n_ev_a, n_samp=n_samp,
                       ev_contrast=ev_c, ps_contrast=ps_c)
            # truth-oriented score: > 0 means the signal favours the TRUE junction
            if pairs[t["pid"]]["type"] == "artifact" and ev_c is not None:
                sgn = 1.0 if t["side"] == "TRUE" else -1.0
                row["ev_true"] = sgn * ev_c
                row["ps_true"] = sgn * ps_c
            else:
                row["ev_true"] = row["ps_true"] = None
            rows.append(row)

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    cols = ["sample", "pair", "pair_type", "side", "rid", "chrom", "strand", "own", "alt", "shift_d",
            "shift_a", "shift_max", "own_clean", "n_introns", "status_own", "status_alt", "init_scale",
            "init_shift", "n_ev_own", "n_ev_alt", "n_samp", "ev_contrast", "ps_contrast", "ev_true",
            "ps_true"]
    tag = f"{args.sample}.hmm{args.hmm}"
    with open(outdir / f"{tag}.reads.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in rows:
            fh.write("\t".join("" if row.get(c) is None else str(row[c]) for c in cols) + "\n")

    # ---- summary -----------------------------------------------------------------------------
    def frac_pos(vals):
        v = [x for x in vals if x is not None]
        return (round(sum(1 for x in v if x > 0) / len(v), 3), len(v)) if v else (None, 0)

    summ = {"sample": args.sample, "hmm": args.hmm, "n_rows": len(rows),
            "n_scored": sum(1 for r in rows if r["ev_contrast"] is not None)}
    art = [r for r in rows if r["pair_type"] == "artifact" and r["ev_contrast"] is not None]
    gen = [r for r in rows if r["pair_type"] == "genuine" and r["ev_contrast"] is not None]
    for name, key in (("per_event", "ev_true"), ("per_sample", "ps_true")):
        blk = {}
        for side in ("FAKE", "TRUE"):
            rr = [r for r in art if r["side"] == side]
            blk[f"{side}_reads_favour_true"] = frac_pos([r[key] for r in rr])
            for lo, hi, lab in ((1, 3, "1-3bp"), (4, 6, "4-6bp"), (7, 99, "7bp+")):
                blk[f"{side}_favour_true_{lab}"] = frac_pos(
                    [r[key] for r in rr if lo <= r["shift_max"] <= hi])
            if side == "FAKE":
                blk["FAKE_unclean_favour_true"] = frac_pos([r[key] for r in rr if not r["own_clean"]])
                blk["FAKE_clean_favour_true"] = frac_pos([r[key] for r in rr if r["own_clean"]])
        # decisive-call analysis on FAKE reads: |score| > margin
        rr = [r for r in art if r["side"] == "FAKE" and r[key] is not None]
        for mg in (1.0, 3.0, 10.0):
            dec = [r for r in rr if abs(r[key]) > mg]
            blk[f"FAKE_margin>{mg}_coverage"] = round(len(dec) / len(rr), 3) if rr else None
            blk[f"FAKE_margin>{mg}_accuracy"] = frac_pos([r[key] for r in dec])[0]
        # genuine pairs: does the signal favour the read's own (both true) junction?
        gk = "ev_contrast" if key == "ev_true" else "ps_contrast"
        blk["GENUINE_reads_favour_own"] = frac_pos([r[gk] for r in gen])
        for lo, hi, lab in ((1, 3, "1-3bp"), (4, 6, "4-6bp"), (7, 99, "7bp+")):
            blk[f"GENUINE_favour_own_{lab}"] = frac_pos(
                [r[gk] for r in gen if lo <= r["shift_max"] <= hi])
        summ[name] = blk
    (outdir / f"{tag}.summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
