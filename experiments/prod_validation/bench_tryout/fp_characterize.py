#!/usr/bin/env python3
"""Characterize FALSE POSITIVES in the pre-gate (max-ablation, 175k) p00 candidate set,
so we can DESIGN a targeted gate instead of toggling SIRV-tuned ones.

Label each candidate by gffcompare class ('=' = true structure; else FP). Join per-
candidate features (scores.unfiltered.tsv) + structural features (from GTF). Report:
  - FP class-code breakdown (what KIND of junk dominates)
  - TP vs FP distribution for every feature (median / quartiles)
  - for each feature, the separability: at a threshold that keeps 95% of TP, how much
    FP does it remove? (a good gate keeps TP, kills FP)
"""
import re, statistics
from collections import defaultdict, Counter

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S = "SGNex_H9_directRNA_replicate2_run2"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
NC = f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
SCORES = f"{HERE}/prodfull_ablate/p00/work/scores.unfiltered.tsv"
GTF = f"{HERE}/prodfull_ablate/p00/pyfin.gtf"
TRK = f"{HERE}/gc_ablate/{S}__p00__pyfin.tracking"
ENST = re.compile(r"(ENST\d+)")

est = {}
for ln in open(NC):
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 3:
        m = ENST.search(p[0])
        if m:
            try: est[m.group(1)] = float(p[2])
            except ValueError: pass
expr3 = {k for k, v in est.items() if v >= 3}

# candidate_id -> class code (+ ref)
cls = {}
for ln in open(TRK):
    c = ln.rstrip("\n").split("\t")
    if len(c) < 5: continue
    code = c[3]; ref = ENST.search(c[2]); ref = ref.group(1) if ref else None
    for q in c[4:]:
        m = re.search(r":([^|]+)\|", q)
        if m: cls[m.group(1)] = (code, ref)

# structural features from GTF: candidate_id -> (n_introns, tx_len, chrom, strand, chain)
tid = re.compile(r'transcript_id "([^"]+)"')
exons = defaultdict(list); meta = {}
for ln in open(GTF):
    if ln.startswith("#"): continue
    f = ln.split("\t")
    if len(f) < 9 or f[2] != "exon": continue
    m = tid.search(f[8])
    if not m: continue
    exons[m.group(1)].append((int(f[3]) - 1, int(f[4])))
    meta[m.group(1)] = (f[0], f[6])
struct = {}
chain_by = {}
for t, ex in exons.items():
    ex.sort()
    tx_len = sum(b - a for a, b in ex)
    ch = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
    struct[t] = (len(ch), tx_len)
    chain_by.setdefault(meta[t], []).append((t, ch))

# is_contained: chain is a strict 3'-suffix subchain of a longer candidate at same locus
contained = set()
for (chrom, strand), lst in chain_by.items():
    lst_sorted = sorted(lst, key=lambda x: -len(x[1]))
    chains = [c for _, c in lst_sorted]
    for t, ch in lst:
        if not ch: continue
        for oc in chains:
            if len(oc) <= len(ch): continue
            if (oc[-len(ch):] == ch) if strand == "+" else (oc[:len(ch)] == ch):
                contained.add(t); break

# features from scores
feat = {}
with open(SCORES) as fh:
    h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
    for ln in fh:
        c = ln.rstrip("\n").split("\t")
        cid = c[ci["candidate_id"]]
        feat[cid] = dict(num_reads=int(float(c[ci["num_reads"]])),
                         abundance=float(c[ci["abundance"]]),
                         tpm=float(c[ci["tpm"]]),
                         coherence=float(c[ci["coherence_score"]]),
                         discrimination=float(c[ci["discrimination_score"]]),
                         combined=float(c[ci["combined_score"]]),
                         max_R=float(c[ci["max_R"]]))

ids = list(feat)
def is_tp(cid): return cls.get(cid, ("", None))[0] == "="
TP = [c for c in ids if is_tp(c)]
FP = [c for c in ids if not is_tp(c)]
print(f"total candidates: {len(ids)}   TP(structural '='): {len(TP)}   FP: {len(FP)}")

# FP class-code breakdown
codes = Counter(cls.get(c, ("u", None))[0] for c in FP)
print("\nFP class-code breakdown:")
for k, v in codes.most_common():
    print(f"  {k:3} {v:7} ({100*v/len(FP):.0f}%)")

# structural: mono vs multi, contained fraction
def frac(sub, whole): return 100 * sum(sub) / len(whole) if whole else 0
print(f"\nmono-exon:   TP {frac([struct.get(c,(1,0))[0]==0 for c in TP],TP):.0f}%   "
      f"FP {frac([struct.get(c,(1,0))[0]==0 for c in FP],FP):.0f}%")
print(f"contained:   TP {frac([c in contained for c in TP],TP):.0f}%   "
      f"FP {frac([c in contained for c in FP],FP):.0f}%")

def q(vals):
    vals = sorted(vals)
    n = len(vals)
    return (vals[n//4], statistics.median(vals), vals[3*n//4]) if n else (0,0,0)

print(f"\n{'feature':14}{'TP q25/med/q75':>26}{'FP q25/med/q75':>26}"
      f"{'keep95%TP->FP kill':>20}")
for fname in ["num_reads", "abundance", "tpm", "combined", "coherence",
              "discrimination", "max_R"]:
    tv = [feat[c][fname] for c in TP]; fv = [feat[c][fname] for c in FP]
    tq = q(tv); fq = q(fv)
    # threshold that keeps 95% of TP (5th percentile of TP), then FP removed below it
    thr = sorted(tv)[max(0, int(0.05 * len(tv)) - 1)]
    fp_removed = 100 * sum(1 for v in fv if v < thr) / len(fv)
    print(f"{fname:14}{tq[0]:7.2f}/{tq[1]:6.2f}/{tq[2]:7.2f}"
          f"{fq[0]:9.2f}/{fq[1]:6.2f}/{fq[2]:7.2f}   thr>={thr:7.2f} kills {fp_removed:4.0f}%FP")

# also add structural n_introns + tx_len
print("\nstructural:")
for fname, getter in [("n_introns", lambda c: struct.get(c,(0,0))[0]),
                      ("tx_len", lambda c: struct.get(c,(0,0))[1])]:
    tv = [getter(c) for c in TP]; fv = [getter(c) for c in FP]
    tq = q(tv); fq = q(fv)
    print(f"  {fname:12} TP {tq}   FP {fq}")
