#!/usr/bin/env python3
"""Read-only probe: are reads supporting FAKE (non-annotated, near-truth)
junctions indel-unclean, and do TRUE junctions ever look unclean?

SIRV_Set4 is synthetic, so every junction absent from the catalog is an
alignment/assembly artifact. For every primary MAPQ>=1 read we walk the
CIGAR exactly like fin/candidates/discovery.py (N -> junction, I/D ->
indel) and record, per junction observation, whether an indel lies within
FLANK bp of the donor or acceptor (the discovery "indel-clean" rule), the
distance to the nearest indel and its size.

Classes (per (chrom, donor, acceptor), 0-based, strand from read flag):
  TRUE       exact catalog junction
  NEAR_FAKE  not in catalog, but a catalog junction within +-NEAR bp on
             BOTH coordinates (wobble artifact)
  FAR_FAKE   everything else

Outputs (per sample + pooled): junction table, read-level clean fractions,
indel-size-vs-shift signature, and the election-screen confusion
(accused / not) for TRUE vs NEAR_FAKE junctions under (a) the current
clean<=1 rule and (b) a clean<=1 AND raw<=RAW_CAP rule.

Evaluation-only: uses the truth GTF. Never a PyFIN input.
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import pysam

FLANK = 10
NEAR = 15
SCREEN_CLEAN_MAX = 1
SCREEN_COMP_MIN = 3
SCREEN_RATIO = 3


def truth_junctions(gtf_path: Path):
    exons = collections.defaultdict(list)   # tid -> [(chrom, strand, start1, end1)]
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            tid = None
            for kv in f[8].split(";"):
                kv = kv.strip()
                if kv.startswith("transcript_id"):
                    tid = kv.split('"')[1]
                    break
            if tid is None:
                continue
            exons[tid].append((f[0], f[6], int(f[3]), int(f[4])))
    truth = {}   # (chrom, d, a) -> set(strands)
    for tid, ex in exons.items():
        ex.sort(key=lambda e: e[2])
        for (c, s, _s1, e1), (_c2, _s2, s2, _e2) in zip(ex, ex[1:]):
            # 1-based inclusive exon end == 0-based first intronic base;
            # next 1-based exon start - 1 == 0-based first exonic base.
            truth.setdefault((c, e1, s2 - 1), set()).add(s)
    return truth


def scan_bam(bam_path: Path, truth):
    """Return per-junction stats and per-read observation lists."""
    J = {}  # key (chrom, strand, d, a) -> dict
    n_reads = 0
    with pysam.AlignmentFile(str(bam_path)) as bam:
        for rd in bam.fetch(until_eof=True):
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                continue
            if rd.mapping_quality < 1:
                continue
            cig = rd.cigartuples
            if not cig:
                continue
            n_reads += 1
            strand = "-" if rd.is_reverse else "+"
            chrom = rd.reference_name
            pos = rd.reference_start
            ind = []   # (lo, hi, size, kind)
            jxs = []
            for op, ln in cig:
                if op in (0, 7, 8):
                    pos += ln
                elif op == 2:
                    ind.append((pos, pos + ln, ln, "D"))
                    pos += ln
                elif op == 1:
                    ind.append((pos, pos, ln, "I"))
                elif op == 3:
                    jxs.append((pos, pos + ln))
                    pos += ln
            for d, a in jxs:
                clean = True
                near_dist = None
                near_size = 0
                near_kind = ""
                for lo, hi, size, kind in ind:
                    dist = min(abs(d - lo), abs(d - hi), abs(a - lo), abs(a - hi))
                    if lo <= d <= hi or lo <= a <= hi:
                        dist = 0
                    if near_dist is None or dist < near_dist:
                        near_dist, near_size, near_kind = dist, size, kind
                    if (lo - FLANK <= d <= hi + FLANK) or (lo - FLANK <= a <= hi + FLANK):
                        clean = False
                key = (chrom, strand, d, a)
                ent = J.get(key)
                if ent is None:
                    ent = J[key] = {"raw": 0, "clean": 0, "dists": [], "sizes": []}
                ent["raw"] += 1
                if clean:
                    ent["clean"] += 1
                else:
                    ent["dists"].append(near_dist)
                    ent["sizes"].append(near_size if near_kind == "D" else -near_size)
    return J, n_reads


def classify(J, truth):
    # index truth by chrom for near search
    by_chrom = collections.defaultdict(list)
    for (c, d, a) in truth:
        by_chrom[c].append((d, a))
    rows = []
    for (c, s, d, a), ent in J.items():
        if (c, d, a) in truth:
            cls, shift_d, shift_a, near_t = "TRUE", 0, 0, (d, a)
        else:
            best = None
            for td, ta in by_chrom.get(c, ()):
                if abs(td - d) <= NEAR and abs(ta - a) <= NEAR:
                    dist = abs(td - d) + abs(ta - a)
                    if best is None or dist < best[0]:
                        best = (dist, td, ta)
            if best:
                cls, shift_d, shift_a, near_t = "NEAR_FAKE", d - best[1], a - best[2], (best[1], best[2])
            else:
                cls, shift_d, shift_a, near_t = "FAR_FAKE", None, None, None
        rows.append({
            "chrom": c, "strand": s, "donor": d, "acceptor": a, "cls": cls,
            "raw": ent["raw"], "clean": ent["clean"],
            "shift_d": shift_d, "shift_a": shift_a, "near_true": near_t,
            "unclean_dists": ent["dists"], "unclean_sizes": ent["sizes"],
        })
    return rows


def election_screen(rows):
    """Apply the election nomination screen interval-free (whole chrom+strand):
    accused iff clean <= 1 and a competitor within +-NEAR bp (both coords)
    has clean >= max(3, 3*max(clean,1))."""
    by_cs = collections.defaultdict(list)
    for r in rows:
        by_cs[(r["chrom"], r["strand"])].append(r)
    for r in rows:
        r["accused_clean_rule"] = False
        r["accused_raw_rule"] = False
        r["best_comp_clean"] = 0
        cj = r["clean"]
        bar = max(SCREEN_COMP_MIN, SCREEN_RATIO * max(cj, 1))
        comps = [
            k for k in by_cs[(r["chrom"], r["strand"])]
            if k is not r and abs(k["donor"] - r["donor"]) <= NEAR
            and abs(k["acceptor"] - r["acceptor"]) <= NEAR
        ]
        if comps:
            r["best_comp_clean"] = max(k["clean"] for k in comps)
            r["best_comp_raw"] = max(k["raw"] for k in comps)
        if cj <= SCREEN_CLEAN_MAX and comps and r["best_comp_clean"] >= bar:
            r["accused_clean_rule"] = True
            # variant: additionally require raw support to be dominated too
            if r["raw"] * SCREEN_RATIO <= r["best_comp_raw"]:
                r["accused_raw_rule"] = True


def summarize(rows, label):
    out = {"label": label}
    for cls in ("TRUE", "NEAR_FAKE", "FAR_FAKE"):
        rr = [r for r in rows if r["cls"] == cls]
        raw = sum(r["raw"] for r in rr)
        clean = sum(r["clean"] for r in rr)
        n_unclean_obs = raw - clean
        dists = [d for r in rr for d in r["unclean_dists"]]
        # junction-level clean-count buckets
        buckets = collections.Counter()
        for r in rr:
            c = r["clean"]
            buckets["0" if c == 0 else "1" if c == 1 else "2-4" if c <= 4 else "5+"] += 1
        info = {
            "n_junctions": len(rr),
            "reads_supporting": raw,
            "read_level_unclean_frac": round(n_unclean_obs / raw, 4) if raw else None,
            "junction_clean_count_buckets": dict(buckets),
            "junctions_raw>=2_clean<=1": sum(1 for r in rr if r["raw"] >= 2 and r["clean"] <= 1),
            "junctions_raw>=5_clean<=1": sum(1 for r in rr if r["raw"] >= 5 and r["clean"] <= 1),
            "accused_clean_rule": sum(1 for r in rr if r["accused_clean_rule"]),
            "accused_raw_rule": sum(1 for r in rr if r["accused_raw_rule"]),
            "reads_in_accused_clean_rule": sum(r["raw"] for r in rr if r["accused_clean_rule"]),
            "reads_in_accused_raw_rule": sum(r["raw"] for r in rr if r["accused_raw_rule"]),
        }
        if dists:
            dists.sort()
            info["unclean_dist_to_indel_quantiles"] = {
                "p10": dists[len(dists) // 10], "p50": dists[len(dists) // 2],
                "p90": dists[(len(dists) * 9) // 10],
            }
        if cls == "NEAR_FAKE":
            # indel signature: nearest indel size (D>0, I<0) vs shift
            sig_hit = sig_tot = 0
            shift_hist = collections.Counter()
            for r in rr:
                shift = abs(r["shift_d"]) if r["shift_d"] else abs(r["shift_a"])
                shift_hist[f"d{r['shift_d']:+d}/a{r['shift_a']:+d}"] += r["raw"]
                for s in r["unclean_sizes"]:
                    sig_tot += 1
                    if abs(s) == shift:
                        sig_hit += 1
            info["unclean_reads_with_indel_size==shift_frac"] = (
                round(sig_hit / sig_tot, 4) if sig_tot else None)
            info["top_shifts_by_reads"] = shift_hist.most_common(12)
            info["clean>=2_escapes"] = sum(1 for r in rr if r["clean"] >= 2)
            info["reads_in_clean>=2_escapes"] = sum(r["raw"] for r in rr if r["clean"] >= 2)
        out[cls] = info
    # TRUE junctions that the screen would wrongly accuse: list them
    out["TRUE_accused_examples"] = [
        {k: r[k] for k in ("chrom", "strand", "donor", "acceptor", "raw", "clean", "best_comp_clean")}
        for r in rows if r["cls"] == "TRUE" and r["accused_clean_rule"]
    ][:40]
    out["TRUE_raw>=5_clean<=1_examples"] = [
        {k: r[k] for k in ("chrom", "strand", "donor", "acceptor", "raw", "clean")}
        | {"p50_dist": sorted(r["unclean_dists"])[len(r["unclean_dists"]) // 2]}
        for r in rows if r["cls"] == "TRUE" and r["raw"] >= 5 and r["clean"] <= 1
    ][:40]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", action="append", required=True, help="sample=path")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--flank", type=int, default=10)
    args = ap.parse_args()
    global FLANK
    FLANK = args.flank
    truth = truth_junctions(Path(args.truth))
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    pooled_rows = {}
    summaries = []
    for spec in args.bam:
        sample, path = spec.split("=", 1)
        J, n_reads = scan_bam(Path(path), truth)
        rows = classify(J, truth)
        election_screen(rows)
        with open(outdir / f"{sample}.junctions.tsv", "w") as fh:
            fh.write("chrom\tstrand\tdonor\tacceptor\tcls\traw\tclean\tshift_d\tshift_a"
                     "\tbest_comp_clean\taccused_clean_rule\taccused_raw_rule\n")
            for r in sorted(rows, key=lambda r: (r["chrom"], r["strand"], r["donor"], r["acceptor"])):
                fh.write("\t".join(str(r[k]) for k in (
                    "chrom", "strand", "donor", "acceptor", "cls", "raw", "clean",
                    "shift_d", "shift_a", "best_comp_clean",
                    "accused_clean_rule", "accused_raw_rule")) + "\n")
        s = summarize(rows, sample)
        s["n_reads"] = n_reads
        summaries.append(s)
        print(json.dumps(s, indent=1), file=sys.stderr)
        for key, ent in J.items():
            p = pooled_rows.setdefault(key, {"raw": 0, "clean": 0, "dists": [], "sizes": []})
            p["raw"] += ent["raw"]; p["clean"] += ent["clean"]
            p["dists"] += ent["dists"]; p["sizes"] += ent["sizes"]
    rows = classify(pooled_rows, truth)
    election_screen(rows)
    summaries.append(summarize(rows, "POOLED"))
    (outdir / "summary.json").write_text(json.dumps(summaries, indent=1))
    print(json.dumps(summaries[-1], indent=1))


if __name__ == "__main__":
    main()
