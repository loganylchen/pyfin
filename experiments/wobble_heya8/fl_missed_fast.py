#!/usr/bin/env python3
"""FAST version: why do FULL-LENGTH-supported truth transcripts produce NO candidate?

Compute FL support (transcriptome BAM) + matched set (cf refmap) FIRST, so the set of
missed FL>=3 transcript ids is known BEFORE parsing the huge truth GTF -- then parse the GTF
keeping ONLY those ids' exons (cheap set lookup + string ops, no full-genome regex).

Args: <transcriptome.bam> <genome.bam> <truth.gtf> <cf_refmap>
"""
import sys
from collections import Counter

import pysam

TBAM, GBAM, TRUTH, CF_RM = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
FL_TOL = 50
TOL = 6
MAXREADS = 6000


def enst(x):
    i = x.find("ENST")
    if i < 0:
        return x
    j = i
    while j < len(x) and (x[j].isalnum() or x[j] == "."):
        j += 1
    return x[i:j]


# 1. FL support from transcriptome BAM
tb = pysam.AlignmentFile(TBAM, "rb")
fl = Counter()
for aln in tb.fetch(until_eof=True):
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        continue
    L = tb.lengths[aln.reference_id]
    if aln.reference_start <= FL_TOL and aln.reference_end >= L - FL_TOL:
        fl[enst(tb.get_reference_name(aln.reference_id))] += 1
fl3 = {t for t, v in fl.items() if v >= 3}

# 2. matched (cf refmap class '=')
matched = set()
with open(CF_RM) as f:
    next(f, None)
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 3 and c[2] == "=":
            matched.add(c[1])

missed_ids = fl3 - matched
print(f"FL>=3={len(fl3)}  matched(of FL>=3)={len(fl3 & matched)}  MISSED={len(missed_ids)}", flush=True)

# 3. parse truth GTF for ONLY the missed ids (fast: string ops + set lookup)
tx = {}
with open(TRUTH) as fh:
    for line in fh:
        if line[0] == "#":
            continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "exon":
            continue
        a = p[8]
        k = a.find('transcript_id "')
        if k < 0:
            continue
        s = k + 15
        tid = a[s:a.find('"', s)]
        if tid not in missed_ids:
            continue
        t = tx.setdefault(tid, {"chrom": p[0], "strand": p[6], "exons": []})
        t["exons"].append((int(p[3]) - 1, int(p[4])))
gt = {}
for tid, t in tx.items():
    ex = sorted(t["exons"])
    intr = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
    gt[tid] = {"chrom": t["chrom"], "strand": t["strand"], "chain": intr,
               "start": ex[0][0], "end": ex[-1][1], "nexon": len(ex)}
print(f"truth chains built for {len(gt)} missed transcripts", flush=True)


def rchain(aln):
    intr = []
    pos = aln.reference_start
    for op, ln in aln.cigartuples:
        if op in (0, 7, 8, 2):
            pos += ln
        elif op == 3:
            intr.append((pos, pos + ln))
            pos += ln
    return tuple(intr)


def close(c1, c2):
    if len(c1) != len(c2) or not c1:
        return False
    return all(abs(a - c) <= TOL and abs(b - d) <= TOL for (a, b), (c, d) in zip(c1, c2))


def is_sub(short, long_):
    n, m = len(long_), len(short)
    if m == 0 or m >= n:
        return False
    return any(all(abs(short[k][0] - long_[o + k][0]) <= TOL and
                   abs(short[k][1] - long_[o + k][1]) <= TOL for k in range(m))
               for o in range(n - m + 1))


gb = pysam.AlignmentFile(GBAM, "rb")
cat = Counter()
exon_by_cat = Counter()
for t in sorted(gt):
    g = gt[t]
    tc = g["chain"]
    if not tc:
        cat["mono_truth"] += 1
        continue
    try:
        it = gb.fetch(g["chrom"], max(g["start"], 0), g["end"])
    except (ValueError, KeyError):
        cat["E_no_genome_reads"] += 1
        continue
    exact = wob = sub = nreads = scanned = 0
    for aln in it:
        scanned += 1
        if scanned > MAXREADS:
            break
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if ("-" if aln.is_reverse else "+") != g["strand"]:
            continue
        ct = aln.cigartuples
        if not ct:
            continue
        if (ct[0][0] == 4 and ct[0][1] >= 250) or (ct[-1][0] == 4 and ct[-1][1] >= 250):
            continue
        rc = rchain(aln)
        if not rc:
            continue
        nreads += 1
        if rc == tc:
            exact += 1
            break
        elif close(rc, tc):
            wob += 1
        elif is_sub(rc, tc):
            sub += 1
    key = ("A_has_EXACT_read" if exact else "B_wobble_only" if wob else
           "C_subchain_only" if sub else "D_other_structure" if nreads else "E_no_genome_reads")
    cat[key] += 1
    b = "2-3ex" if g["nexon"] <= 3 else ("4-8ex" if g["nexon"] <= 8 else ">8ex")
    exon_by_cat[(key, b)] += 1

print(f"\n{'category':<22}{'count':>7}{'%':>7}")
for k in ["A_has_EXACT_read", "B_wobble_only", "C_subchain_only",
          "D_other_structure", "E_no_genome_reads", "mono_truth"]:
    v = cat.get(k, 0)
    print(f"{k:<22}{v:>7}{100 * v / max(len(missed_ids), 1):>6.1f}%")
print("\nby exon-count:")
for k in ["B_wobble_only", "C_subchain_only", "D_other_structure"]:
    print(f"  {k:<20} " + " ".join(f"{b}={exon_by_cat.get((k, b), 0)}"
                                   for b in ["2-3ex", "4-8ex", ">8ex"]))
