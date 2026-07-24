#!/usr/bin/env python3
"""RECALL oracle: the precision oracle showed honF1 caps at 32.3 because corrRec is
stuck at 22.8. corrRec only rises if we recover MORE expressed-truth transcripts. A
truncated 'c' candidate is a fragment of a real truth transcript; if that transcript
is EXPRESSED (est>=3) and not already matched '=', then PERFECT assembly of the 'c' to
full length would add it to corrRec. This measures that recall ceiling.

Also 'k' (ref contained in query) and other partials point at truth transcripts we
partially see. We count distinct expressed truth transcripts that are *pointed at* by a
non-'=' candidate but not yet '=' matched -> the assembly-recoverable recall pool.
"""
import os, re
from collections import defaultdict

HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S = "SGNex_H9_directRNA_replicate2_run2"
NC = f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
TRACK = f"{HERE}/gc_full/{S}__p00__pyfin.tracking"
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
DENOM = len(expr3)

# per class code, the truth transcript each candidate points at
matched_eq = set()
pointed = defaultdict(set)   # code -> set(ENST it points at)
for ln in open(TRACK):
    c = ln.rstrip("\n").split("\t")
    if len(c) < 4: continue
    m = ENST.search(c[2])
    if not m: continue
    e = m.group(1); code = c[3]
    if code == "=":
        matched_eq.add(e)
    else:
        pointed[code].add(e)

cur_expr = matched_eq & expr3
print(f"expressed-truth denom: {DENOM}")
print(f"current '=' matched (all): {len(matched_eq)}  expressed: {len(cur_expr)}  "
      f"corrRec={100*len(cur_expr)/DENOM:.1f}%")

# assembly-recoverable: expressed truth pointed at by 'c' (truncation) but not '=' matched
for codes, label in [(["c"], "c only"),
                     (["c", "k"], "c+k"),
                     (["c", "k", "j", "m", "n", "e", "o", "i"], "all partial/near")]:
    pool = set()
    for cc in codes:
        pool |= pointed.get(cc, set())
    recover = (pool & expr3) - cur_expr
    new_rec = 100 * (len(cur_expr) + len(recover)) / DENOM
    # honF1 at this recall, assuming precision oracle (honPr up to ~55 as orc4)
    for hp in (46.4, 55.1):
        h = 0 if (hp + new_rec) == 0 else 2 * hp * new_rec / (hp + new_rec)
        print(f"  [{label:16}] +{len(recover):4} expressed recovered -> corrRec "
              f"{100*len(cur_expr)/DENOM:.1f}->{new_rec:.1f}   honF1@honPr{hp:.0f} = {h:.1f}")
print(f"\nisoquant p00 reference: honF1=42.1")
