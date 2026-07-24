#!/usr/bin/env python3
"""Oracle upper-bound study: how high can p00 precision / honest-F1 go if we could
PERFECTLY remove wrong candidates? Builds pruned pyfin GTFs under several oracle drop
policies and writes them to prodfull_orc<N>/p00/pyfin.gtf for eval_honest scoring.

Policies (all keep every '=' correct candidate):
  orc1  drop non-'=' multi-exon candidates CONTAINED (3'-suffix) in a longer candidate
        (the realistic terminal/containment lever ceiling; reads reassign to the longer)
  orc2  drop ALL 'c' truncations (mono+multi)                     [Codex ~80.9% prec]
  orc3  drop ALL 'c' + all fake-junction candidates              [Codex ~84% prec]
  orc4  drop EVERY non-'=' candidate (perfect oracle: precision->100%, recall unchanged)
Then run eval_honest.py on each dir to get real honF1 (recall+orphan effects included).
"""
import os, re, sys
from collections import defaultdict

HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
S = "SGNex_H9_directRNA_replicate2_run2"
PYGTF = f"{HERE}/prodfull/p00/pyfin.gtf"
TRACK = f"{HERE}/gc_full/{S}__p00__pyfin.tracking"
TOL = 5


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


truth_chains, truth_meta = introns_by_tx(TRUTH)
truth_junc = defaultdict(set)
for t, ch in truth_chains.items():
    for j in ch:
        truth_junc[truth_meta[t][0]].add(j)


def is_real(chrom, j):
    d, a = j
    if (d, a) in truth_junc[chrom]:
        return True
    for rd in range(d - TOL, d + TOL + 1):
        for ra in range(a - TOL, a + TOL + 1):
            if (rd, ra) in truth_junc[chrom]:
                return True
    return False


chains, meta = introns_by_tx(PYGTF)
code = {}
for ln in open(TRACK):
    c = ln.rstrip("\n").split("\t")
    if len(c) < 5: continue
    for q in c[4:]:
        m = re.search(r"\|([^|]+)\|", q)
        if m: code[m.group(1)] = c[3]

by_cs = defaultdict(list)
for t, ch in chains.items():
    if ch:
        by_cs[meta[t]].append((t, ch))


def contained_in_longer(t):
    ch = chains[t]; chrom, strand = meta[t]
    for lt, lch in by_cs[(chrom, strand)]:
        if lt == t or len(lch) <= len(ch):
            continue
        if (lch[-len(ch):] == ch) if strand == "+" else (lch[:len(ch)] == ch):
            return True
    return False


def has_fake(t):
    return any(not is_real(meta[t][0], j) for j in chains[t])


drop = {"orc1": set(), "orc2": set(), "orc3": set(), "orc4": set()}
for t, ch in chains.items():
    cc = code.get(t, "?")
    if cc == "=":
        continue
    # orc4: drop every non-'='
    drop["orc4"].add(t)
    # orc3: all 'c' + fake-junction
    if cc == "c" or (ch and has_fake(t)):
        drop["orc3"].add(t)
    # orc2: all 'c'
    if cc == "c":
        drop["orc2"].add(t)
    # orc1: non-'=' multi-exon contained in a longer candidate
    if ch and contained_in_longer(t):
        drop["orc1"].add(t)

lines = open(PYGTF).readlines()
tid = re.compile(r'transcript_id "([^"]+)"')
for name, dropset in drop.items():
    outdir = f"{HERE}/prodfull_{name}/p00"
    os.makedirs(outdir, exist_ok=True)
    with open(f"{outdir}/pyfin.gtf", "w") as o:
        for ln in lines:
            if ln.startswith("#"):
                o.write(ln); continue
            m = tid.search(ln)
            if m and m.group(1) in dropset:
                continue
            o.write(ln)
    kept = sum(1 for l in open(f"{outdir}/pyfin.gtf") if "\ttranscript\t" in l)
    print(f"{name}: dropped {len(dropset):5}  kept {kept:5} tx  -> prodfull_{name}/p00/pyfin.gtf")
print("\nNow score: for N in orc1 orc2 orc3 orc4; do python3 eval_honest.py prodfull_$N gc_$N; done")
