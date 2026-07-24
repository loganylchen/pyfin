#!/usr/bin/env python3
"""The 30 `exact_read` cases: a genome read HAS the exact truth chain, yet pyfin did not emit a '='.
Where does each die? Recompute the exact_read set, then for each report BOTH:
  (1) gffcompare class_code assigned to that truth (from cf refmap) + the query id list, and
  (2) direct comparison of the truth chain vs pyfin's EMITTED chains (cf.gtf) at the locus:
      exact_emitted   : pyfin emitted the exact chain (=> a measurement/gffcompare quirk, not a loss)
      wobble_emitted  : emitted a <=TOLbp-off chain (generation/selection snapped to the wrong variant)
      subchain_emit   : emitted a truncated sub-chain (fold/selection kept a shorter isoform)
      super_emitted   : truth is a sub-chain of a longer emitted chain (folded INTO a longer isoform)
      none_emitted    : nothing multi-exon emitted at the locus (dropped in generation/clustering)

Args: <fl3.txt> <cf_refmap> <truth.gtf> <genome.bam> <pyfin_cf.gtf>
"""
import sys
from collections import Counter, defaultdict

import pysam

FL3, REFMAP, TRUTH, BAM, CFGTF = sys.argv[1:6]
TOL = 6


def rchain(aln):
    intr, pos = [], aln.reference_start
    for op, ln in aln.cigartuples or ():
        if op in (0, 7, 8, 2):
            pos += ln
        elif op == 3:
            intr.append((pos, pos + ln)); pos += ln
    return tuple(intr)


def close(c1, c2):
    return (len(c1) == len(c2) and c1 and
            all(abs(a - c) <= TOL and abs(b - d) <= TOL for (a, b), (c, d) in zip(c1, c2)))


def is_sub(short, long_):
    n, m = len(long_), len(short)
    if m == 0 or m > n:
        return False
    return any(all(abs(short[k][0] - long_[o + k][0]) <= TOL and
                   abs(short[k][1] - long_[o + k][1]) <= TOL for k in range(m))
               for o in range(n - m + 1))


fl3 = {ln.strip() for ln in open(FL3) if ln.strip()}
matched, rm_class, rm_qry = set(), {}, {}
with open(REFMAP) as f:
    next(f, None)
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 4:
            rm_class[c[1]] = c[2]; rm_qry[c[1]] = c[3]
            if c[2] == "=":
                matched.add(c[1])
missed = fl3 - matched

# truth chains for missed
tx = {}
with open(TRUTH) as fh:
    for line in fh:
        if line[0] == "#":
            continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "exon":
            continue
        a = p[8]; k = a.find('transcript_id "')
        if k < 0:
            continue
        s = k + 15; tid = a[s:a.find('"', s)]
        if tid not in missed:
            continue
        tx.setdefault(tid, {"chrom": p[0], "strand": p[6], "exons": []})["exons"].append(
            (int(p[3]) - 1, int(p[4])))
gt = {}
for tid, t in tx.items():
    ex = sorted(t["exons"])
    gt[tid] = {"chrom": t["chrom"], "strand": t["strand"],
               "chain": tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)),
               "start": ex[0][0], "end": ex[-1][1]}

# find exact_read set via genome bam
gb = pysam.AlignmentFile(BAM, "rb")
exact_ids = []
for t in sorted(gt):
    g = gt[t]; tc = g["chain"]
    if not tc:
        continue
    try:
        it = gb.fetch(g["chrom"], max(g["start"], 0), g["end"])
    except (ValueError, KeyError):
        continue
    n = 0
    for aln in it:
        n += 1
        if n > 8000:
            break
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if ("-" if aln.is_reverse else "+") != g["strand"]:
            continue
        ct = aln.cigartuples
        if not ct or (ct[0][0] == 4 and ct[0][1] >= 250) or (ct[-1][0] == 4 and ct[-1][1] >= 250):
            continue
        if rchain(aln) == tc:
            exact_ids.append(t); break
print(f"exact_read cases: {len(exact_ids)}", flush=True)

# parse pyfin emitted cf.gtf into chains by chrom
emit = defaultdict(list)  # chrom -> list of (strand, chain, start, end)
cur = {}
exons = defaultdict(list)
meta = {}
with open(CFGTF) as fh:
    for line in fh:
        if line[0] == "#":
            continue
        p = line.split("\t")
        if len(p) < 9:
            continue
        a = p[8]; k = a.find('transcript_id "')
        if k < 0:
            continue
        s = k + 15; tid = a[s:a.find('"', s)]
        if p[2] == "exon":
            exons[tid].append((int(p[3]) - 1, int(p[4])))
            meta[tid] = (p[0], p[6])
for tid, ex in exons.items():
    ex.sort(); chrom, strand = meta[tid]
    ch = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
    emit[chrom].append((strand, ch, ex[0][0], ex[-1][1], tid))

fate = Counter()
cc = Counter()
print(f"\n{'transcript':<20}{'nex':>4} {'class':>6} {'emit_fate':<14} qry")
for t in exact_ids:
    g = gt[t]; tc = g["chain"]; nex = len(tc) + 1
    cc[rm_class.get(t, "(absent)")] += 1
    cands = [e for e in emit.get(g["chrom"], [])
             if e[0] == g["strand"] and not (e[3] < g["start"] or e[2] > g["end"])]
    f = "none_emitted"
    for strand, ch, s0, e0, tid in cands:
        if ch == tc:
            f = "exact_emitted"; break
        if close(ch, tc):
            f = "wobble_emitted"
        elif is_sub(tc, ch) and f not in ("wobble_emitted",):
            f = "super_emitted"       # truth is subchain of a longer emitted
        elif is_sub(ch, tc) and f == "none_emitted":
            f = "subchain_emit"       # emitted is a truncated subchain of truth
    fate[f] += 1
    print(f"{t:<20}{nex:>4} {rm_class.get(t,'-'):>6} {f:<14} {rm_qry.get(t,'-')[:40]}")

print(f"\n=== gffcompare class_code of the {len(exact_ids)} ===")
for k, v in cc.most_common():
    print(f"  {k:<10}{v}")
print(f"\n=== pyfin emitted-chain fate ===")
for k in ["exact_emitted", "wobble_emitted", "super_emitted", "subchain_emit", "none_emitted"]:
    if fate.get(k):
        print(f"  {k:<16}{fate[k]}")
