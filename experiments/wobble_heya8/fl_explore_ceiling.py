#!/usr/bin/env python3
"""Of the missed FL>=3 truth transcripts, how many can the explore step (exact read-phased
intron-graph tiling) actually reconstruct, vs are fundamentally unrecoverable from CIGAR alone?

For each missed truth locus, collect the exact intron chains of the genome reads (primary,
strand-matched, fusion-excluded) and build the exact read-phased edge set. Classify the truth:

  exact_read     : some read already has the exact truth chain (no assembly needed) -> selection bug
  explore_exact  : NOT on any single read, but every truth junction is observed on SOME read AND
                   every consecutive truth junction-pair is co-observed on SOME read -> explore
                   can tile it exactly (the recoverable path-completion ceiling)
  missing_edge   : some consecutive truth junction-pair is on NO read (both junctions may exist,
                   never adjacent) -> global-recombination risk / needs phasing beyond CIGAR
  missing_junction: some truth junction is on NO read at all -> correlated truncation; unrecoverable
                   from CIGAR, needs GTF or raw-signal evidence
  mono / no_reads : single-exon / no usable reads

Args: <fl3.txt> <cf_refmap> <truth.gtf> <genome.bam>
"""
import sys
from collections import Counter

import pysam

FL3, REFMAP, TRUTH, BAM = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
MAXREADS = 8000


def rchain(aln):
    intr = []
    pos = aln.reference_start
    for op, ln in aln.cigartuples or ():
        if op in (0, 7, 8, 2):
            pos += ln
        elif op == 3:
            intr.append((pos, pos + ln))
            pos += ln
    return tuple(intr)


fl3 = {ln.strip() for ln in open(FL3) if ln.strip()}
matched = set()
with open(REFMAP) as f:
    next(f, None)
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 3 and c[2] == "=":
            matched.add(c[1])
missed = fl3 - matched

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

gb = pysam.AlignmentFile(BAM, "rb")
cat = Counter()
for t in sorted(gt):
    g = gt[t]
    tc = g["chain"]
    if not tc:
        cat["mono"] += 1
        continue
    try:
        it = gb.fetch(g["chrom"], max(g["start"], 0), g["end"])
    except (ValueError, KeyError):
        cat["no_reads"] += 1
        continue
    chains = []
    scanned = 0
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
        if rc:
            chains.append(rc)
    if not chains:
        cat["no_reads"] += 1
        continue
    if any(rc == tc for rc in chains):
        cat["exact_read"] += 1
        continue
    introns = set()
    edges = set()
    for rc in chains:
        introns.update(rc)
        edges.update(zip(rc, rc[1:]))
    if any(j not in introns for j in tc):
        cat["missing_junction"] += 1
    elif any((tc[i], tc[i + 1]) not in edges for i in range(len(tc) - 1)):
        cat["missing_edge"] += 1
    else:
        cat["explore_exact"] += 1

tot = len(missed)
print(f"MISSED FL>=3 = {tot}\n")
print(f"{'class':<18}{'count':>7}{'%':>8}   interpretation")
rows = [
    ("exact_read", "exact chain IS on a read -> selection/emit bug"),
    ("explore_exact", "tileable across reads -> explore CAN recover"),
    ("missing_edge", "junctions seen but never adjacent -> phasing gap"),
    ("missing_junction", "a junction on NO read -> unrecoverable from CIGAR"),
    ("mono", "single-exon truth"),
    ("no_reads", "no usable primary reads"),
]
for k, desc in rows:
    v = cat.get(k, 0)
    print(f"{k:<18}{v:>7}{100*v/max(tot,1):>7.1f}%   {desc}")
absent = tot - sum(cat.values())
if absent:
    print(f"{'(not in GTF)':<18}{absent:>7}{100*absent/max(tot,1):>7.1f}%")
rec = cat.get("exact_read", 0) + cat.get("explore_exact", 0)
print(f"\nCIGAR-recoverable (exact_read + explore_exact) = {rec}/{tot} = {100*rec/max(tot,1):.1f}%")
