#!/usr/bin/env python3
"""TPM-stratified structural Sn/Pr for a tool GTF vs truth, using NanoCount TPM as the
expression axis. Also prints the pure NanoCount stratification (how many truth transcripts
sit in each TPM band) so recall can be read against a real denominator.

Sn@T = |matched & {truth with tpm>=T}| / |{truth with tpm>=T}|   (recall over expressed truth)
Pr@T = |matched & {truth with tpm>=T}| / n_output                (honest precision at floor T)
Pr(any) = |matched (any truth)| / n_output

Args: <refmap> <truth.gtf> <nanocount.tsv> <n_output>
"""
import re
import sys

REFMAP, TRUTH, NC, NOUT = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
THRS = [0.0, 0.1, 1.0, 3.0, 10.0, 30.0, 100.0]

# truth transcript ids
truth = set()
for ln in open(TRUTH):
    if ln.startswith("#"):
        continue
    p = ln.split("\t")
    if len(p) < 9 or p[2] != "exon":
        continue
    m = re.search(r'transcript_id "([^"]+)"', p[8])
    if m:
        truth.add(m.group(1))

# nanocount tpm (col 3), est_count (col 2)
tpm = {t: 0.0 for t in truth}
est = {t: 0.0 for t in truth}
fh = open(NC)
next(fh, None)
for ln in fh:
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 4 and p[0] in tpm:
        try:
            est[p[0]] = float(p[2])
            tpm[p[0]] = float(p[3])
        except ValueError:
            pass

# matched (gffcompare class '=' in refmap col1=ref_id col2=class)
matched = set()
with open(REFMAP) as f:
    next(f, None)
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 3 and c[2] == "=":
            matched.add(c[1])

n_expr_any = sum(1 for v in tpm.values() if v > 0)
print(f"truth_total={len(truth)}  truth_with_tpm>0={n_expr_any}  "
      f"output={NOUT}  matched(any_truth)={len(matched)}  "
      f"Pr(any)={100 * len(matched) / max(NOUT, 1):.1f}%")
print(f"{'TPM>=':>7}{'#truth':>9}{'#matched':>9}{'Sn%(recall)':>13}{'Pr%':>7}")
prev = None
for T in THRS:
    expr = truth if T == 0.0 else {t for t, v in tpm.items() if v >= T}
    m = len(matched & expr)
    sn = 100 * m / max(len(expr), 1)
    pr = 100 * m / max(NOUT, 1)
    band = ""
    if prev is not None:
        band = f"  (band {prev}-{T}: {len(prevset - expr)} truth)"
    print(f"{T:>7}{len(expr):>9}{m:>9}{sn:>13.1f}{pr:>7.1f}{band}")
    prev, prevset = T, expr

# est_count>=3 reference point (the standard honest-metrics denominator)
e3 = {t for t, v in est.items() if v >= 3}
m3 = len(matched & e3)
print(f"\n[ref] est_count>=3: #truth={len(e3)}  #matched={m3}  "
      f"Sn={100 * m3 / max(len(e3), 1):.1f}%  Pr={100 * m3 / max(NOUT, 1):.1f}%")
