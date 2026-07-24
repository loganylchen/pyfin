#!/usr/bin/env python3
"""Per-mode precision/recall frontier from a broad (min-ab 0) run's unfiltered scores
+ its gffcompare tracking. Sweep an abundance floor; report corrRec at several
structPr levels so M1/M2/M3 modes can be compared: the mode whose frontier gives the
HIGHEST corrRec at a given structPr concentrates read mass best (kills FP without
killing TP). Usage: frontier_of.py <pyfin_dir> <gc_dir> <label>
"""
import re, sys
from collections import defaultdict

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S = "SGNex_H9_directRNA_replicate2_run2"
NC = f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
ENST = re.compile(r"(ENST\d+)")
pyf, gc, label = sys.argv[1], sys.argv[2], sys.argv[3]
SCORES = f"{pyf}/p00/work/scores.unfiltered.tsv"
TRK = f"{gc}/{S}__p00__pyfin.tracking"

est = {}
for ln in open(NC):
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 3:
        m = ENST.search(p[0])
        if m:
            try: est[m.group(1)] = float(p[2])
            except ValueError: pass
expr3 = {k for k, v in est.items() if v >= 3}; DEN = len(expr3)
cls = {}
for ln in open(TRK):
    c = ln.rstrip("\n").split("\t")
    if len(c) < 5: continue
    code = c[3]; ref = ENST.search(c[2]); ref = ref.group(1) if ref else None
    for q in c[4:]:
        m = re.search(r":([^|]+)\|", q)
        if m: cls[m.group(1)] = (code, ref)
feat = {}
with open(SCORES) as fh:
    h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
    for ln in fh:
        c = ln.rstrip("\n").split("\t")
        feat[c[ci["candidate_id"]]] = (float(c[ci["abundance"]]),
                                       int(float(c[ci["num_reads"]])))
ids = list(feat)
def evalset(keep):
    tot = len(keep)
    if not tot: return (0, 0)
    eq = [c for c in keep if cls.get(c, ("", None))[0] == "="]
    refs = {cls[c][1] for c in eq if cls[c][1] in expr3}
    return 100 * len(eq) / tot, 100 * len(refs) / DEN
print(f"### {label}: frontier (min-abundance floor sweep on broad set)")
print(f"{'floor':>7}{'tx':>9}{'structPr':>9}{'corrRec':>8}")
for T in [0, 0.5, 1, 2, 3, 4, 5, 6, 8, 10]:
    keep = [c for c in ids if feat[c][0] >= T]
    sp, rec = evalset(keep)
    print(f"{T:>7}{len(keep):>9}{sp:>8.1f}%{rec:>8.1f}")
