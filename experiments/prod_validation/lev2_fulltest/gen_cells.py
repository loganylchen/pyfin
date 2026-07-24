#!/usr/bin/env python3
"""Generate the full Lever-2 test matrix cell list.

cells.tsv columns: ds  sample  ratio  gtf(path|-)  off_src(path|-)  gpu(0/1)  run_off(0/1)
- SIRV4/gencode: OFF baselines are e268c9b (byte-identical) -> reuse (run_off=0).
- heya8: out_bp20 OFF is a pre-e268c9b build -> run OFF live (run_off=1).
ratio 'nogtf' => no --gtf (true de novo); else --gtf _ref/<ratio>/annotation.gtf.
"""
import os, glob
REPO = "/SSD/logan/dev/pyfin"
EXP = f"{REPO}/experiments"
rows = []

# --- SIRV4: reuse OFF (prod_validation/sirv4/<S>/<ratio>/pyfin.gtf) ---
S4 = f"{EXP}/prod_validation/sirv4"
for off in glob.glob(f"{S4}/SGNex_*/*/pyfin.gtf"):
    parts = off.split("/"); samp, ratio = parts[-3], parts[-2]
    gtf = f"{S4}/_ref/{ratio}/annotation.gtf"
    if not os.path.exists(gtf):
        gtf = "-"
    rows.append(("sirv", samp, ratio, gtf, off, "0", "0"))

# --- gencode: reuse OFF (prod_validation/gencode/<S>/<ratio>/pyfin.gtf) ---
GC = f"{EXP}/prod_validation/gencode"
for off in glob.glob(f"{GC}/SGNex_*/*/pyfin.gtf"):
    parts = off.split("/"); samp, ratio = parts[-3], parts[-2]
    if ratio == "stage":
        continue
    gtf = "-" if ratio == "nogtf" else f"{GC}/_ref/{ratio}/annotation.gtf"
    if gtf != "-" and not os.path.exists(gtf):
        continue
    rows.append(("gencode", samp, ratio, gtf, off, "1", "0"))

# --- heya8: run BOTH off+on live (provenance) ---
HM = f"{EXP}/wobble_heya8/matrix"
ratios = [os.path.basename(os.path.dirname(p)) for p in glob.glob(f"{HM}/_ref/*/annotation.gtf")]
samples = sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(f"{HM}/SGNex_HEYA8*/stage"))
for samp in samples:
    for ratio in ratios:
        gtf = f"{HM}/_ref/{ratio}/annotation.gtf"
        rows.append(("heya8", samp, ratio, gtf, "-", "0", "1"))
    rows.append(("heya8", samp, "nogtf", "-", "-", "0", "1"))  # true de novo

out = f"{EXP}/prod_validation/lev2_fulltest/cells.tsv"
with open(out, "w") as fh:
    for r in rows:
        fh.write("\t".join(r) + "\n")
from collections import Counter
c = Counter(r[0] for r in rows)
print(f"wrote {len(rows)} cells -> {out}")
print("by dataset:", dict(c))
print("runs (on + off-live):", len(rows) + sum(1 for r in rows if r[6] == "1"))
