#!/usr/bin/env python3
"""Trace the ACTUAL full-length reads of missed transcripts by READ NAME.

For each missed FL>=3 truth transcript, take the QNAMEs of its full-length reads (from the
TRANSCRIPTOME BAM), then look up those SAME QNAMEs in the GENOME BAM and see what happened to each
read at the transcript's genome locus. This is a read-name join, not a locus scan -- it distinguishes
"the read is truncated/wobbled in the genome spliced alignment" from "the read's primary genome
alignment is at another locus" from "the read was fusion-excluded" etc.

Per (transcript, full-length read) fate:
  genome_exact         primary at locus, genome chain == truth chain (should have been recalled)
  genome_wobble        primary at locus, same #introns, every junction <=6bp off
  genome_subchain      primary at locus, genome chain is a truncated contiguous sub-chain of truth
  genome_other         primary at locus, different structure
  fusion_excluded      primary at locus but >=250bp terminal soft-clip (generation drops it)
  only_sec_suppl       no primary at locus, but a secondary/supplementary alignment is there
  primary_elsewhere    the read's primary genome alignment is at a different locus/strand (multi-map)
  not_in_genome        the QNAME has no genome alignment at all

Args: <transcriptome.bam> <genome.bam> <truth.gtf> <cf_refmap>
"""
import sys
from collections import Counter, defaultdict

import pysam

TBAM, GBAM, TRUTH, CF_RM = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
FL_TOL = 50
TOL = 6


def enst(x):
    i = x.find("ENST")
    if i < 0:
        return x
    j = i
    while j < len(x) and (x[j].isalnum() or x[j] == "."):
        j += 1
    return x[i:j]


def chain_of(aln):
    intr = []
    pos = aln.reference_start
    for op, ln in aln.cigartuples or ():
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


# 1. transcriptome BAM: FL support + FL read QNAMEs per transcript
import time
t0=time.time()
print("scanning transcriptome BAM...", flush=True)
tb = pysam.AlignmentFile(TBAM, "rb")
fl_reads = defaultdict(set)   # transcript -> set(qname) of full-length reads
n_scan = 0
for aln in tb.fetch(until_eof=True):
    n_scan += 1
    if n_scan % 20000000 == 0:
        print(f'  scanned {n_scan//1000000}M in {time.time()-t0:.0f}s', flush=True)
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        continue
    L = tb.lengths[aln.reference_id]
    if aln.reference_start <= FL_TOL and aln.reference_end >= L - FL_TOL:
        fl_reads[enst(tb.get_reference_name(aln.reference_id))].add(aln.query_name)
fl3 = {t for t, q in fl_reads.items() if len(q) >= 3}
print(f"transcriptome scan done in {time.time()-t0:.0f}s", flush=True)

# 2. matched (cf refmap)
matched = set()
with open(CF_RM) as f:
    next(f, None)
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 3 and c[2] == "=":
            matched.add(c[1])
missed = fl3 - matched
print(f"FL>=3={len(fl3)}  MISSED={len(missed)}", flush=True)

# 3. truth chains for missed ids (filtered fast parse)
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

# target QNAMEs = FL reads of missed transcripts, mapped per transcript
target = set()
for t in missed:
    target |= fl_reads[t]
print(f"missed truth chains={len(gt)}  target FL-read QNAMEs={len(target)}", flush=True)

# 4. genome BAM: per missed-transcript LOCUS fetch (fast, indexed) -> qname -> alns AT THAT LOCUS
gb = pysam.AlignmentFile(GBAM, "rb")

def locus_alns(g, keep):
    out = {}
    try:
        it = gb.fetch(g["chrom"], max(g["start"], 0), g["end"])
    except (ValueError, KeyError):
        return out
    scanned = 0
    for aln in it:
        scanned += 1
        if scanned > 60000:
            break
        if aln.is_unmapped:
            continue
        q = aln.query_name
        if q not in keep:
            continue
        ct = aln.cigartuples or []
        sc0 = ct[0][1] if ct and ct[0][0] == 4 else 0
        scN = ct[-1][1] if ct and ct[-1][0] == 4 else 0
        out.setdefault(q, []).append({
            "strand": "-" if aln.is_reverse else "+", "sec": aln.is_secondary,
            "suppl": aln.is_supplementary, "chain": chain_of(aln), "sc0": sc0, "scN": scN,
        })
    return out

def fate(g, q, la):
    at_locus = [a for a in la.get(q, []) if a["strand"] == g["strand"]]
    if not at_locus:
        return "not_at_locus"          # primary elsewhere OR unmapped (read not at truth locus)
    prim = [a for a in at_locus if not a["sec"] and not a["suppl"]]
    if not prim:
        return "only_sec_suppl"
    a = prim[0]
    if a["sc0"] >= 250 or a["scN"] >= 250:
        return "fusion_excluded"
    rc = a["chain"]
    tc = g["chain"]
    if rc == tc:
        return "genome_exact"
    if close(rc, tc):
        return "genome_wobble"
    if is_sub(rc, tc):
        return "genome_subchain"
    return "genome_other"


per_read = Counter()
per_tx_best = Counter()
BEST = ["genome_exact", "genome_wobble", "genome_subchain", "genome_other",
        "fusion_excluded", "only_sec_suppl", "not_at_locus"]
done = 0
for t in sorted(missed):
    done += 1
    if done % 100 == 0:
        print(f"  ...{done}/{len(missed)}", flush=True)
    la = locus_alns(gt[t], fl_reads[t])
    fates = [fate(gt[t], q, la) for q in fl_reads[t]]
    for f in fates:
        per_read[f] += 1
    best = min((BEST.index(f) for f in fates), default=len(BEST) - 1)
    per_tx_best[BEST[best]] += 1

tot_reads = sum(per_read.values())
print(f"\n=== per FULL-LENGTH READ fate (n={tot_reads}) ===")
for k in BEST:
    v = per_read.get(k, 0)
    print(f"  {k:<18}{v:>7}{100 * v / max(tot_reads, 1):>6.1f}%")
print(f"\n=== per MISSED TRANSCRIPT, best read fate (n={len(missed)}) ===")
for k in BEST:
    v = per_tx_best.get(k, 0)
    print(f"  {k:<18}{v:>7}{100 * v / max(len(missed), 1):>6.1f}%")
