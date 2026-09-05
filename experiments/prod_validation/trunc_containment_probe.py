#!/usr/bin/env python3
"""Offline probe: can the production containment tests settle truncation?

User hypothesis: a truncated (shorter) model contained in a longer supported
structure should be kept only if a statistical test says a genuine short
transcript exists inside the long one's read pile - e.g. an excess spike over
the expected 5' coverage-decay trend.

This is exactly what fin/analysis/tss_evidence.py implements for dRNA:

- rung `tss_only`  (5'-truncated entry, shares chain suffix + TES with the
  long form): conditional termination-hazard test `evaluate_internal_tss` -
  degradation-null vs an extra read population starting at the claimed TSS;
- rung `own_tes`   (3'-truncated entry, shares 5' prefix, own premature TES):
  `evaluate_tes_support` - dRNA 3' ends are reliable, so a distinct 3'
  cluster at the claimed TES is positive evidence, its absence is negative.

The probe runs both production functions offline over emitted P0 GTF models
in `c_trunc5p50` / `c_trunc3p50`, using the full SIRV truth structure as a
stand-in for the read-derived long competitor that P2 will supply at
inference time. Expression labels (tmap) are joined afterwards and never feed
the tests. Positive controls come from `full` cells (genuine nested SIRV
isoforms). Read-only: no fin/ code changed.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import pysam

sys.path.insert(0, "/SSD/logan/dev/pyfin")
from fin.analysis.tss_evidence import (  # noqa: E402
    build_hazard_profile,
    evaluate_internal_tss,
    evaluate_tes_support,
    genomic_to_offset,
    pooled_background_hazard,
    read_five_prime_offset,
    spliced_length,
)

BENCH = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/"
             "results/sirv4_full_sweep_ratio_grid_v2")
W = Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid")
TRUTH = W / "_ref" / "SIRV_Set4.gtf"
SAMPLES = ["SGNex_H9_directRNA_replicate3_run1",
           "SGNex_H9_directRNA_replicate4_run1"]
CONDS = ["c_trunc5p50", "c_trunc3p50", "full"]
TOL = 6


def parse_gtf(path, want_source=False):
    models = {}
    for line in Path(path).open():
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        m = models.setdefault(tid, {"chrom": f[0], "strand": f[6],
                                    "exons": [], "source": None})
        m["exons"].append((int(f[3]) - 1, int(f[4])))  # 0-based half-open
    if want_source:
        for line in Path(path).open():
            if "\ttranscript\t" not in line:
                continue
            tid = re.search(r'transcript_id "([^"]+)"', line).group(1)
            src = re.search(r'transcript_source "([^"]+)"', line)
            if tid in models and src:
                models[tid]["source"] = src.group(1)
    for m in models.values():
        m["exons"] = sorted(m["exons"])
        ex = m["exons"]
        m["chain"] = tuple((ex[i][1], ex[i + 1][0])
                           for i in range(len(ex) - 1))
    return models


def read_index(bam_path):
    per_chrom = defaultdict(list)
    with pysam.AlignmentFile(str(bam_path)) as bam:
        for aln in bam.fetch():
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue
            pos = aln.reference_start
            chain = []
            for op, ln in aln.cigartuples:
                if op in (0, 7, 8, 2):
                    pos += ln
                elif op == 3:
                    chain.append((pos, pos + ln))
                    pos += ln
            per_chrom[aln.reference_name].append(
                (tuple(chain), aln.reference_start, aln.reference_end))
    return per_chrom


def close(a, b, tol=TOL):
    return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol


def chain_eq(a, b, tol=TOL):
    return len(a) == len(b) and all(close(x, y, tol) for x, y in zip(a, b))


def is_sub_at(read, model, offset, tol=TOL):
    return all(close(r, model[offset + i], tol) for i, r in enumerate(read))


def containment(short, long, tol=TOL):
    """Strand-aware relation of S's chain inside L's chain."""
    cs, cl = short["chain"], long["chain"]
    if len(cs) >= len(cl) or not cs:
        return None
    n, m = len(cs), len(cl)
    strand = long["strand"]
    if is_sub_at(cs, cl, m - n, tol):
        genomic = "suffix"
    elif is_sub_at(cs, cl, 0, tol):
        genomic = "prefix"
    else:
        return None
    if strand == "+":
        return "tss_only" if genomic == "suffix" else "own_tes"
    return "tss_only" if genomic == "prefix" else "own_tes"


def pool_reads(per_chrom, L):
    """Reads attributable to the L-family: any contiguous subchain of L's
    chain at any offset (full reads, 5'-degraded tails, short-isoform reads
    that stop before L's later junctions), plus unspliced fragments lying
    inside L's genomic span (degradation/internal-termination products)."""
    cl = L["chain"]
    ex = L["exons"]
    lo, hi = ex[0][0], ex[-1][1]
    out = []
    for chain, start, end in per_chrom.get(L["chrom"], []):
        keep = False
        if chain:
            n, m = len(chain), len(cl)
            if n <= m:
                for off in range(m - n + 1):
                    if is_sub_at(chain, cl, off):
                        keep = True
                        break
        elif lo - TOL <= start and end <= hi + TOL:
            keep = True
        if keep:
            out.append((start, end))
    return out


def main() -> int:
    truth = parse_gtf(TRUTH)
    rows = []
    for sample in SAMPLES:
        per_chrom = read_index(
            BENCH / sample / "align" / f"{sample}.sorted.bam")
        # pooled degradation null from truth transcripts with no containment
        # partner (standalone), mirroring production's "no contained sibling"
        contained_ids = set()
        for a in truth.values():
            for bid, b in truth.items():
                if a is b or a["chrom"] != b["chrom"]:
                    continue
                if containment(b, a):
                    contained_ids.add(bid)
        profiles = []
        for tid, L in truth.items():
            if tid in contained_ids:
                continue
            pool = pool_reads(per_chrom, L)
            offs = [read_five_prime_offset((s, e), L["exons"], L["strand"])
                    for s, e in pool]
            offs = [o for o in offs if o is not None]
            if len(offs) >= 20:
                profiles.append(build_hazard_profile(
                    offs, spliced_length(L["exons"])))
        bg = pooled_background_hazard(profiles)
        for cond in CONDS:
            asm = parse_gtf(W / "p0_strict" / sample / cond /
                            "out" / "assembly.gtf", want_source=True)
            tmap = (W / "expression_score" / "cells" / sample / cond /
                    "pyfin_p0" / "cmp.query.gtf.tmap")
            lab = {r["qry_id"]: ("TP" if r["class_code"] == "=" else "FP")
                   for r in csv.DictReader(tmap.open(), delimiter="\t")}
            for tid, S in asm.items():
                if S["source"] != "gtf" or not S["chain"]:
                    continue
                best = None
                for lid, L in truth.items():
                    if L["chrom"] != S["chrom"]:
                        continue
                    rung = containment(S, L)
                    if rung and (best is None
                                 or len(L["chain"]) > len(best[1]["chain"])):
                        best = (lid, L, rung)
                if not best:
                    continue
                lid, L, rung = best
                pool = pool_reads(per_chrom, L)
                strand = L["strand"]
                offs5 = [read_five_prime_offset((s, e), L["exons"], strand)
                         for s, e in pool]
                offs5 = [o for o in offs5 if o is not None]
                if rung == "tss_only":
                    g5 = (S["exons"][0][0] if strand == "+"
                          else S["exons"][-1][1] - 1)
                    d0 = genomic_to_offset(g5, L["exons"], strand)
                    if d0 is None:
                        continue
                    ev = evaluate_internal_tss(
                        candidate_id=tid, parent_id=lid, offsets=offs5,
                        tss_offset=d0, background_hazard=bg,
                        identifiability="tss_only", n_eligible_bins=1)
                else:
                    offs3 = []
                    for s, e in pool:
                        g3 = e - 1 if strand == "+" else s
                        o = genomic_to_offset(g3, L["exons"], strand)
                        if o is not None:
                            offs3.append(o)
                    g3c = (S["exons"][-1][1] - 1 if strand == "+"
                           else S["exons"][0][0])
                    dc = genomic_to_offset(g3c, L["exons"], strand)
                    if dc is None:
                        continue
                    ev = evaluate_tes_support(
                        candidate_id=tid, parent_id=lid,
                        read_three_prime_offsets=offs3,
                        candidate_tes_offset=dc,
                        parent_tes_offset=spliced_length(L["exons"]) - 1)
                rows.append({
                    "sample": sample, "cond": cond, "tid": tid,
                    "parent": lid, "rung": rung,
                    "label": lab.get(tid, "NA"), "verdict": ev.verdict,
                    "reason": ev.reason, "n_peak": ev.n_peak,
                    "n_at_risk": ev.n_at_risk,
                    "peak_fraction": round(ev.peak_fraction, 4),
                    "p_value": (round(ev.p_value, 6)
                                if ev.p_value is not None else ""),
                })
    out = W / "trunc_containment_probe.tsv"
    with out.open("w") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(rows)
    print(f"wrote {out} rows={len(rows)}")
    agg = defaultdict(int)
    for r in rows:
        agg[(r["cond"], r["rung"], r["label"], r["verdict"])] += 1
    for k in sorted(agg):
        print("%-12s %-8s %-2s %-15s %d" % (*k, agg[k]))
    print("\nreasons:")
    rag = defaultdict(int)
    for r in rows:
        rag[(r["cond"], r["rung"], r["label"], r["verdict"],
             r["reason"].split("(")[0])] += 1
    for k in sorted(rag):
        print("%-12s %-8s %-2s %-13s %-38s %d" % (*k, rag[k]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
