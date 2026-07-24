#!/usr/bin/env python3
"""Ceiling of lever A (5'-TSS-guarded truncation DROP), no signal needed.

For pyfin p00 output: how many 'c' (truncation) candidates are a 3'-suffix sub-chain
of a LONGER pyfin candidate at the same locus (so their reads can reassign to it)?
That longer candidate is the reassignment target. Split by whether the longer target
is itself CORRECT ('='): dropping a 'c' whose target is '=' is a clean precision win
with reads flowing to the right transcript. Also count how many longer targets exist
per 'c'. This bounds lever A before any 5'-TSS filtering (the 5' guard only REMOVES
some of these to protect real short isoforms, so this is the upper bound).
"""
import os, re, sys
from collections import defaultdict

HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
S = "SGNex_H9_directRNA_replicate2_run2"
PYGTF = f"{HERE}/prodfull/p00/pyfin.gtf"
TRACK = f"{HERE}/gc_full/{S}__p00__pyfin.tracking"


def introns_by_tx(path):
    tid = re.compile(r'transcript_id "([^"]+)"')
    exons = defaultdict(list); meta = {}
    for ln in open(path):
        if ln.startswith("#"): continue
        f = ln.split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = tid.search(f[8])
        if not m: continue
        exons[m.group(1)].append((int(f[3]) - 1, int(f[4])))
        meta[m.group(1)] = (f[0], f[6])
    chains = {}
    for t, ex in exons.items():
        ex.sort()
        chains[t] = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
    return chains, meta


chains, meta = introns_by_tx(PYGTF)
code = {}
for ln in open(TRACK):
    c = ln.rstrip("\n").split("\t")
    if len(c) < 5: continue
    for q in c[4:]:
        m = re.search(r"\|([^|]+)\|", q)
        if m: code[m.group(1)] = c[3]

# index candidates by chrom/strand for containment search
by_cs = defaultdict(list)
for t, ch in chains.items():
    if ch:
        by_cs[meta[t]].append((t, ch))


def is_3p_suffix(short, long_, strand):
    if len(short) >= len(long_): return False
    return (long_[-len(short):] == short) if strand == "+" else (long_[:len(short)] == short)


c_ids = [t for t in chains if code.get(t) == "c" and chains[t]]
n_contained = 0; n_target_correct = 0; n_target_any_correct = 0
for t in c_ids:
    ch = chains[t]; chrom, strand = meta[t]
    targets = [(lt, lch) for lt, lch in by_cs[(chrom, strand)]
               if lt != t and is_3p_suffix(ch, lch, strand)]
    if targets:
        n_contained += 1
        if any(code.get(lt) == "=" for lt, _ in targets):
            n_target_any_correct += 1

tot_c = len(c_ids)
print(f"pyfin p00 multi-exon 'c' truncations: {tot_c}")
print(f"  contained (3'-suffix) in a LONGER pyfin candidate: {n_contained} "
      f"({100*n_contained/tot_c if tot_c else 0:.0f}%)")
print(f"  of those, a CORRECT '=' longer target exists (clean drop+reassign): "
      f"{n_target_any_correct} ({100*n_target_any_correct/tot_c if tot_c else 0:.0f}%)")
print(f"  NOT contained (no reassignment target; dropping = pure recall loss): "
      f"{tot_c - n_contained}")
print(f"\nUPPER BOUND lever A droppable (before 5'-TSS guard removes real short "
      f"isoforms): {n_contained} of {tot_c} 'c'.")
print(f"Interpretation: if all {n_target_any_correct} with a '=' target are dropped and "
      f"reads reassign correctly, tx {5102}->{5102-n_target_any_correct}, struct-prec "
      f"3484/{5102-n_target_any_correct} = {100*3484/(5102-n_target_any_correct):.1f}% "
      f"(from 68.3%).")
