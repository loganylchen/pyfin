#!/usr/bin/env python3
"""Re-score recall/precision using FULL-LENGTH transcriptome support as the expressed-truth
denominator (instead of NanoCount est_count, which counts partial reads too).

A transcript is "full-length supported" if the transcriptome BAM has >=k reads whose alignment
reaches BOTH ends of the transcript (reference_start <= FL_TOL and reference_end >= L - FL_TOL).
This is the fair de-novo-detectable denominator: an assembler can only recover a transcript that
has full-length read evidence.

Args: <transcriptome.bam> <truth.gtf> <cf_refmap> <prod_refmap> <cf_nout> <prod_nout>
"""
import re
import sys
from collections import Counter

import pysam

TBAM, TRUTH, CF_RM, PROD_RM = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
CF_NOUT, PROD_NOUT = int(sys.argv[5]), int(sys.argv[6])
FL_TOL = 50
ENST = re.compile(r'(ENST\d+\.\d+)')


def norm(x):
    m = ENST.search(x)
    return m.group(1) if m else x


b = pysam.AlignmentFile(TBAM, "rb")
reflen = b.lengths
fl = Counter()
tot = Counter()
n = 0
for aln in b.fetch(until_eof=True):
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        continue
    L = reflen[aln.reference_id]
    tid = norm(b.get_reference_name(aln.reference_id))
    tot[tid] += 1
    if aln.reference_start <= FL_TOL and aln.reference_end >= L - FL_TOL:
        fl[tid] += 1
    n += 1
print(f"primary reads={n}  transcripts_with_any_read={len(tot)}  "
      f"with_FL>=1={sum(1 for v in fl.values() if v >= 1)}  "
      f"with_FL>=3={sum(1 for v in fl.values() if v >= 3)}", flush=True)

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


def matched(rm):
    s = set()
    with open(rm) as f:
        next(f, None)
        for ln in f:
            c = ln.rstrip("\n").split("\t")
            if len(c) >= 3 and c[2] == "=":
                s.add(c[1])
    return s


cf = matched(CF_RM)
prod = matched(PROD_RM)
fl_in_truth = {t: fl.get(t, 0) for t in truth}

print(f"\ndenominator = TRUTH transcripts with >=k FULL-LENGTH reads (FL_TOL={FL_TOL}bp both ends)")
print(f"cf output={CF_NOUT}  prod output={PROD_NOUT}\n")
print(f"{'FL>=':>5}{'#truth':>9} | {'cf_Sn%':>7}{'cf_Pr%':>7} | {'prod_Sn%':>9}{'prod_Pr%':>9}")
print("-" * 60)
for k in [1, 2, 3, 5, 10]:
    denom = {t for t, v in fl_in_truth.items() if v >= k}
    d = max(len(denom), 1)
    cm, pm = len(cf & denom), len(prod & denom)
    print(f"{k:>5}{len(denom):>9} | {100 * cm / d:>7.1f}{100 * cm / max(CF_NOUT, 1):>7.1f} | "
          f"{100 * pm / d:>9.1f}{100 * pm / max(PROD_NOUT, 1):>9.1f}")

# distribution of FL support across truth
print("\nFL-support distribution over truth transcripts:")
for lo, hi in [(0, 1), (1, 2), (2, 3), (3, 5), (5, 10), (10, 10 ** 9)]:
    c = sum(1 for v in fl_in_truth.values() if lo <= v < hi)
    print(f"  FL in [{lo},{hi if hi < 10**9 else 'inf'}): {c}")
