#!/usr/bin/env python3
"""Offline: from the max-ablation p00 candidate set (175k, scores.unfiltered.tsv with
per-candidate abundance + num_reads) and its gffcompare tracking, sweep a simple
min-num-reads / min-abundance filter and show the precision/recall tradeoff WITHOUT
re-running pyfin. Finds whether a simple read/abundance floor recovers a good
operating point (high corrRec + high precision) on the broadly-generated set.

CAVEAT: in the ablation run reads are EM-split across 175k candidates, so a correct
candidate's num_reads is a LOWER bound vs a cleaner run -> this understates recall at
a given threshold. Directional, free, decisive for the shape.
"""
import re
from collections import defaultdict

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S = "SGNex_H9_directRNA_replicate2_run2"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
NC = f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
SCORES = f"{HERE}/prodfull_ablate/p00/work/scores.unfiltered.tsv"
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
DEN = len(expr3)

# candidate_id -> (class_code, ref_enst)
cls = {}
for ln in open(TRK):
    c = ln.rstrip("\n").split("\t")
    if len(c) < 5:
        continue
    code = c[3]
    ref = ENST.search(c[2]); ref = ref.group(1) if ref else None
    for q in c[4:]:
        m = re.search(r":([^|]+)\|", q)   # q1:candidate_id|...
        if m:
            cls[m.group(1)] = (code, ref)

# candidate_id -> (num_reads, abundance)
feat = {}
with open(SCORES) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ci = {c: i for i, c in enumerate(h)}
    for ln in fh:
        c = ln.rstrip("\n").split("\t")
        feat[c[ci["candidate_id"]]] = (int(float(c[ci["num_reads"]])),
                                       float(c[ci["abundance"]]))

def sweep(metric_idx, name, thresholds):
    print(f"\n=== filter by {name} >= T (on 175k ablated candidates) ===")
    print(f"{'T':>5}{'tx':>9}{'=cands':>8}{'structPr':>9}{'exprMatched':>12}{'corrRec':>8}")
    for T in thresholds:
        kept = [cid for cid, ft in feat.items() if ft[metric_idx] >= T]
        tot = len(kept)
        eq = sum(1 for cid in kept if cls.get(cid, ("", None))[0] == "=")
        expr_refs = set()
        for cid in kept:
            code, ref = cls.get(cid, ("", None))
            if code == "=" and ref in expr3:
                expr_refs.add(ref)
        sp = 100 * eq / tot if tot else 0
        rec = 100 * len(expr_refs) / DEN
        print(f"{T:>5}{tot:>9}{eq:>8}{sp:>8.1f}%{len(expr_refs):>12}{rec:>8.1f}")

sweep(0, "num_reads", [1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30])
sweep(1, "abundance", [0, 1, 2, 3, 4, 5, 6, 8, 10])
print(f"\n(OFF production: tx5102 structPr~68 corrRec22.8 | isoquant structPr91 corrRec34.4 |"
      f" ablate-all tx175531 corrRec53.2)")
