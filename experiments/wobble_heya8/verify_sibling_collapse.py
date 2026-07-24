#!/usr/bin/env python3
"""Are the 30 `exact_read` "misses" real losses, or gffcompare identical-chain sibling collapse?

For each of the 30 (pyfin emitted the exact chain but the truth ID is absent from the refmap),
check whether ANOTHER gencode transcript with the byte-identical intron chain DID get '=' in the
refmap. If so, pyfin structurally recovered the isoform and the '=' was merely attributed to a
sibling ID => per-tid recall undercount, NOT a selection bug.

Args: <fl3.txt> <cf_refmap> <truth.gtf> <genome.bam>
"""
import sys
from collections import defaultdict

import pysam

FL3, REFMAP, TRUTH, BAM = sys.argv[1:5]


def rchain(aln):
    intr, pos = [], aln.reference_start
    for op, ln in aln.cigartuples or ():
        if op in (0, 7, 8, 2):
            pos += ln
        elif op == 3:
            intr.append((pos, pos + ln)); pos += ln
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

# parse FULL gencode: per-transcript exons -> chain (multi-exon), and chain -> set(tid)
tx = defaultdict(lambda: {"chrom": "", "strand": "", "exons": []})
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
        t = tx[tid]; t["chrom"] = p[0]; t["strand"] = p[6]
        t["exons"].append((int(p[3]) - 1, int(p[4])))

chain_of = {}
chain2tids = defaultdict(set)
for tid, t in tx.items():
    ex = sorted(t["exons"])
    if len(ex) < 2:
        continue
    ch = (t["chrom"], t["strand"],
          tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)))
    chain_of[tid] = (ch, ex[0][0], ex[-1][1])
    chain2tids[ch].add(tid)

# recompute the 30 exact_read ids
gb = pysam.AlignmentFile(BAM, "rb")
exact_ids = []
for tid in sorted(missed):
    if tid not in chain_of:
        continue
    (chrom, strand, tc), start, end = chain_of[tid]
    if not tc:
        continue
    try:
        it = gb.fetch(chrom, max(start, 0), end)
    except (ValueError, KeyError):
        continue
    n = 0
    for aln in it:
        n += 1
        if n > 8000:
            break
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if ("-" if aln.is_reverse else "+") != strand:
            continue
        ct = aln.cigartuples
        if not ct or (ct[0][0] == 4 and ct[0][1] >= 250) or (ct[-1][0] == 4 and ct[-1][1] >= 250):
            continue
        if rchain(aln) == tc:
            exact_ids.append(tid); break

sib_recovered = real_loss = 0
print(f"{'transcript':<20}{'nex':>4}  verdict")
for tid in exact_ids:
    (ch, _, _), = (chain_of[tid],)
    ch = chain_of[tid][0]
    sibs = chain2tids[ch] - {tid}
    sib_in_matched = sibs & matched
    sib_in_fl3 = sibs & fl3
    if sib_in_matched:
        sib_recovered += 1
        v = f"SIBLING-COLLAPSE  '=' went to {sorted(sib_in_matched)[0]}  (#sibs={len(sibs)})"
    else:
        real_loss += 1
        v = f"REAL-LOSS  no identical-chain sibling matched (#sibs={len(sibs)}, sibs_in_fl3={len(sib_in_fl3)})"
    print(f"{tid:<20}{len(ch[2])+1:>4}  {v}")

print(f"\n=== of {len(exact_ids)} exact_read cases ===")
print(f"  sibling-collapse (structurally recovered, per-tid undercount): {sib_recovered}")
print(f"  real loss (pyfin emitted exact chain but NO sibling got '='):  {real_loss}")
