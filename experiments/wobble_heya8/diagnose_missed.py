#!/usr/bin/env python3
"""Diagnose WHY est_count>=3 truth transcripts were NOT recovered by cluster_families.

For each missed truth, fetch its locus reads (primary, strand-matched, fusion-excluded) and
classify the truth by the BEST read evidence available:
  A_has_EXACT_read  : >=1 read whose intron chain EXACTLY equals the truth chain
                      -> RECOVERABLE (clustering/emission/interval bug: the structure exists)
  B_wobble_only     : >=1 read wobble-close (<=6bp/junction) but NONE exact
                      -> no-snap ceiling (would need consensus/snap to hit '=')
  C_subchain_only   : >=1 read is a truncated contiguous sub-chain (no exact/wobble full chain)
                      -> path-completion territory (only partial reads seen)
  D_other_structure : reads present but none exact/wobble/sub (different structure / novel isoform)
  E_no_reads        : no usable primary reads at the locus (est_count is transcriptome-based)
  mono_truth        : single-exon truth (handled by the mono finalizer, separate track)

Args: <truth.gtf> <nanocount.tsv> <cf_refmap> <bam>
"""
import re
import sys
from collections import Counter

import pysam

TRUTH, NC, REFMAP, BAM = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
TOL = 6


def chains(path):
    tx = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if not m:
            continue
        t = tx.setdefault(m.group(1), {"chrom": p[0], "strand": p[6], "exons": []})
        t["exons"].append((int(p[3]) - 1, int(p[4])))
    out = {}
    for tid, t in tx.items():
        ex = sorted(t["exons"])
        intr = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        out[tid] = {"chrom": t["chrom"], "strand": t["strand"], "chain": intr,
                    "start": ex[0][0], "end": ex[-1][1], "nexon": len(ex)}
    return out


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


def close(c1, c2, tol=TOL):
    if len(c1) != len(c2) or len(c1) == 0:
        return False
    return all(abs(a - c) <= tol and abs(b - d) <= tol for (a, b), (c, d) in zip(c1, c2))


def is_sub(short, long_, tol=TOL):
    n, m = len(long_), len(short)
    if m == 0 or m >= n:
        return False
    for off in range(n - m + 1):
        if all(abs(short[k][0] - long_[off + k][0]) <= tol and
               abs(short[k][1] - long_[off + k][1]) <= tol for k in range(m)):
            return True
    return False


gt = chains(TRUTH)
est = {t: 0.0 for t in gt}
fh = open(NC)
next(fh, None)
for ln in fh:
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 3 and p[0] in est:
        try:
            est[p[0]] = float(p[2])
        except ValueError:
            pass
est3 = {t for t, v in est.items() if v >= 3}

matched = set()
with open(REFMAP) as f:
    next(f, None)
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 3 and c[2] == "=":
            matched.add(c[1])

missed = est3 - matched
print(f"est>=3={len(est3)}  matched={len(est3 & matched)}  MISSED={len(missed)}\n")

bam = pysam.AlignmentFile(BAM, "rb")
cat = Counter()
exon_hist = Counter()   # missed truth by exon-count bucket
for t in sorted(missed):
    g = gt[t]
    tc = g["chain"]
    if not tc:
        cat["mono_truth"] += 1
        continue
    try:
        it = bam.fetch(g["chrom"], max(g["start"], 0), g["end"])
    except (ValueError, KeyError):
        cat["E_no_reads"] += 1
        continue
    exact = wob = sub = nreads = 0
    for aln in it:
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
        elif close(rc, tc):
            wob += 1
        elif is_sub(rc, tc):
            sub += 1
    if exact > 0:
        key = "A_has_EXACT_read"
    elif wob > 0:
        key = "B_wobble_only"
    elif sub > 0:
        key = "C_subchain_only"
    elif nreads > 0:
        key = "D_other_structure"
    else:
        key = "E_no_reads"
    cat[key] += 1
    b = "2-3ex" if g["nexon"] <= 3 else ("4-8ex" if g["nexon"] <= 8 else ">8ex")
    exon_hist[(key, b)] += 1

print(f"{'category':<22}{'count':>7}{'%':>7}")
for k in ["A_has_EXACT_read", "B_wobble_only", "C_subchain_only",
          "D_other_structure", "E_no_reads", "mono_truth"]:
    v = cat.get(k, 0)
    print(f"{k:<22}{v:>7}{100 * v / max(len(missed), 1):>6.1f}%")
print("\nby exon-count (multi-exon categories):")
for k in ["A_has_EXACT_read", "B_wobble_only", "C_subchain_only", "D_other_structure"]:
    row = " ".join(f"{b}={exon_hist.get((k, b), 0)}" for b in ["2-3ex", "4-8ex", ">8ex"])
    print(f"  {k:<20} {row}")
