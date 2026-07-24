#!/usr/bin/env python3
"""Codex's decisive experiment: of the EXPRESSED truth transcripts that isoquant
matches ('=') at p00 but pyfin MISSES, how many have >=2 raw reads in pyfin's own
no-snap BAM that trace the FULL intron chain with both ends in place? That count
decides: >=350 => pyfin has a candidate-GENERATION bug (recoverable +3 recall);
far below => the recall is not in the data under no-snap (bank the corruption win).
"""
import os, re
import pysam
from collections import defaultdict

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
S = "SGNex_H9_directRNA_replicate2_run2"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
BAM = (f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam")
NC = f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
ISO_TRK = f"{HERE}/gc_frontier/{S}__p00__isoquant.tracking"
PYF_TRK = f"{HERE}/gc_frontier/{S}__p00__pyfin.tracking"
ENST = re.compile(r"(ENST\d+)")
TOL = 6
ENDTOL = 50   # read terminus must be within this of the transcript end


def matched_eq(trk):
    s = set()
    for ln in open(trk):
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 4 and c[3] == "=":
            m = ENST.search(c[2])
            if m: s.add(m.group(1))
    return s


est = {}
for ln in open(NC):
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 3:
        m = ENST.search(p[0])
        if m:
            try: est[m.group(1)] = float(p[2])
            except ValueError: pass
expr3 = {k for k, v in est.items() if v >= 3}

iso = matched_eq(ISO_TRK)
pyf = matched_eq(PYF_TRK)
iso_only = (iso - pyf) & expr3
print(f"isoquant '=' expressed: {len(iso&expr3)}  pyfin '=' expressed: {len(pyf&expr3)}")
print(f"isoquant-only expressed (iso matched, pyfin missed): {len(iso_only)}")

# truth chains + spans for iso_only
tx = {}
for ln in open(TRUTH):
    if ln.startswith("#"): continue
    f = ln.split("\t")
    if len(f) < 9 or f[2] != "exon": continue
    m = ENST.search(f[8])
    if not m or m.group(1) not in iso_only: continue
    tx.setdefault(m.group(1), [f[0], f[6], []])[2].append((int(f[3]) - 1, int(f[4])))
TX = {}
for e, (chrom, strand, ex) in tx.items():
    ex.sort()
    introns = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
    TX[e] = dict(chrom=chrom, strand=strand, introns=introns,
                 lo=ex[0][0], hi=ex[-1][1])


def read_introns(r):
    out = []
    pos = r.reference_start
    for op, ln in r.cigartuples:
        if op in (0, 2, 7, 8):
            pos += ln
        elif op == 3:
            out.append((pos, pos + ln)); pos += ln
    return out


def chain_contains(read_ch, t_ch):
    """t_ch is a contiguous subchain of read_ch within TOL on every junction."""
    n, m = len(read_ch), len(t_ch)
    if m == 0 or m > n:
        return False
    for off in range(n - m + 1):
        if all(abs(read_ch[off + k][0] - t_ch[k][0]) <= TOL and
               abs(read_ch[off + k][1] - t_ch[k][1]) <= TOL for k in range(m)):
            return True
    return False


def chain_equal(read_ch, t_ch):
    """read chain EXACTLY equals t_ch (same length) within TOL -> read is a
    full-length copy of THIS transcript, not a longer isoform containing it."""
    if len(read_ch) != len(t_ch) or not t_ch:
        return False
    return all(abs(read_ch[k][0] - t_ch[k][0]) <= TOL and
               abs(read_ch[k][1] - t_ch[k][1]) <= TOL for k in range(len(t_ch)))


bam = pysam.AlignmentFile(BAM, "rb")
full_reads = {}     # enst -> count of full-chain reads
full_reads_ends = {}  # enst -> count of full-chain reads WITH both ends in place
full_reads_exact = {}
for e, t in TX.items():
    if not t["introns"]:
        continue
    nfull = 0; nfull_ends = 0; nexact = 0
    for r in bam.fetch(t["chrom"], max(0, t["lo"]), t["hi"]):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        rch = read_introns(r)
        if chain_contains(rch, t["introns"]):
            nfull += 1
            if (r.reference_start <= t["lo"] + ENDTOL and
                    r.reference_end >= t["hi"] - ENDTOL):
                nfull_ends += 1
        # STRICT: read chain == transcript chain, both ends snug (two-sided)
        if (chain_equal(rch, t["introns"]) and
                abs(r.reference_start - t["lo"]) <= ENDTOL and
                abs(r.reference_end - t["hi"]) <= ENDTOL):
            nexact += 1
    full_reads[e] = nfull
    full_reads_ends[e] = nfull_ends
    full_reads_exact[e] = nexact
bam.close()

ge2_span = sum(1 for e in full_reads if full_reads[e] >= 2)
ge2_full = sum(1 for e in full_reads_ends if full_reads_ends[e] >= 2)
ge1_full = sum(1 for e in full_reads_ends if full_reads_ends[e] >= 1)
ge2_exact = sum(1 for e in full_reads_exact if full_reads_exact[e] >= 2)
ge1_exact = sum(1 for e in full_reads_exact if full_reads_exact[e] >= 1)
print(f"\nOf {len(TX)} multi-exon iso-only expressed transcripts:")
print(f"  >=2 reads SPAN the full intron chain:              {ge2_span}")
print(f"  >=2 reads span full chain AND both ends in place:  {ge2_full}")
print(f"  >=2 reads with EXACT chain == transcript, snug:    {ge2_exact}  <-- strict (unambiguous bug)")
print(f"  >=1 exact full-length read:                        {ge1_exact}")
print(f"\nDECISION (Codex): >=350 => candidate-generation bug, recoverable +3 recall.")
print(f"                  far below 350 => recall not in the data under no-snap.")
