#!/usr/bin/env python3
"""Why do FULL-LENGTH-supported truth transcripts still produce NO matching candidate?

Denominator = truth with >=3 full-length reads in the TRANSCRIPTOME BAM. Of those NOT matched
('=') by cluster_families, probe the GENOME BAM at each locus and classify by best read evidence
there -- to see whether the transcriptome-full-length reads simply don't yield a clean exact
GENOME intron chain (genome-alignment wobble / truncation / mapped elsewhere).

  A_has_EXACT_read : >=1 genome read whose chain EXACTLY equals the truth chain  -> emission bug
  B_wobble_only    : genome reads wobble-close (<=6bp) but none exact             -> genome-align wobble
  C_subchain_only  : only truncated sub-chain genome reads                        -> genome truncation/soft-clip
  D_other_structure: genome reads present, different structure                    -> multi-map / novel
  E_no_genome_reads: no usable primary genome reads at the locus                  -> primary-mapped elsewhere
  mono_truth       : single-exon truth

Args: <transcriptome.bam> <genome.bam> <truth.gtf> <cf_refmap>
"""
import re
import sys
from collections import Counter

import pysam

TBAM, GBAM, TRUTH, CF_RM = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
FL_TOL = 50
TOL = 6
ENST = re.compile(r'(ENST\d+\.\d+)')


def norm(x):
    m = ENST.search(x)
    return m.group(1) if m else x


def chains(path):
    tx = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if not m:
            continue
        t = tx.setdefault(m.group(1), {"chrom": p[0], "strand": p[6], "exons": []})
        t["exons"].append((int(p[3]) - 1, int(p[4])))
    out = {}
    for tid, t in tx.items():
        ex = sorted(t["exons"])
        intr = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        out[tid] = {"chrom": t["chrom"], "strand": t["strand"], "chain": intr,
                    "start": ex[0][0], "end": ex[-1][1], "nexon": len(ex)}
    return out


def rchain(aln):
    intr = []
    pos = aln.reference_start
    for op, ln in aln.cigartuples:
        if op in (0, 7, 8, 2):
            pos += ln
        elif op == 3:
            intr.append((pos, pos + ln))
            pos += ln
    return tuple(intr)


def close(c1, c2, tol=TOL):
    if len(c1) != len(c2) or len(c1) == 0:
        return False
    return all(abs(a - c) <= tol and abs(b - d) <= tol for (a, b), (c, d) in zip(c1, c2))


def is_sub(short, long_, tol=TOL):
    n, m = len(long_), len(short)
    if m == 0 or m >= n:
        return False
    return any(all(abs(short[k][0] - long_[off + k][0]) <= tol and
                   abs(short[k][1] - long_[off + k][1]) <= tol for k in range(m))
               for off in range(n - m + 1))


# FL support per transcript from transcriptome BAM
tb = pysam.AlignmentFile(TBAM, "rb")
fl = Counter()
for aln in tb.fetch(until_eof=True):
    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
        continue
    L = tb.lengths[aln.reference_id]
    if aln.reference_start <= FL_TOL and aln.reference_end >= L - FL_TOL:
        fl[norm(tb.get_reference_name(aln.reference_id))] += 1

gt = chains(TRUTH)
matched = set()
with open(CF_RM) as f:
    next(f, None)
    for ln in f:
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 3 and c[2] == "=":
            matched.add(c[1])

fl3 = {t for t in gt if fl.get(t, 0) >= 3}
missed = fl3 - matched
print(f"FL>=3 truth={len(fl3)}  matched={len(fl3 & matched)}  MISSED={len(missed)}\n", flush=True)

gb = pysam.AlignmentFile(GBAM, "rb")
MAXREADS = 6000   # cap per-locus scan (dense genes) -- enough to see exact/wobble/sub evidence
cat = Counter()
done = 0
for t in sorted(missed):
    done += 1
    if done % 100 == 0:
        print(f"  ...{done}/{len(missed)}", flush=True)
    g = gt[t]
    tc = g["chain"]
    if not tc:
        cat["mono_truth"] += 1
        continue
    try:
        it = gb.fetch(g["chrom"], max(g["start"], 0), g["end"])
    except (ValueError, KeyError):
        cat["E_no_genome_reads"] += 1
        continue
    exact = wob = sub = nreads = scanned = 0
    for aln in it:
        scanned += 1
        if scanned > MAXREADS:
            break
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if ("-" if aln.is_reverse else "+") != g["strand"]:
            continue
        ct = aln.cigartuples
        if not ct:
            continue
        if (ct[0][0] == 4 and ct[0][1] >= 250) or (ct[-1][0] == 4 and ct[-1][1] >= 250):
            continue
        rc = rchain(aln)
        if not rc:
            continue
        nreads += 1
        if rc == tc:
            exact += 1
            break                       # exact found -> category A, stop scanning
        elif close(rc, tc):
            wob += 1
        elif is_sub(rc, tc):
            sub += 1
    if exact:
        cat["A_has_EXACT_read"] += 1
    elif wob:
        cat["B_wobble_only"] += 1
    elif sub:
        cat["C_subchain_only"] += 1
    elif nreads:
        cat["D_other_structure"] += 1
    else:
        cat["E_no_genome_reads"] += 1

print(f"{'category':<22}{'count':>7}{'%':>7}")
for k in ["A_has_EXACT_read", "B_wobble_only", "C_subchain_only",
          "D_other_structure", "E_no_genome_reads", "mono_truth"]:
    v = cat.get(k, 0)
    print(f"{k:<22}{v:>7}{100 * v / max(len(missed), 1):>6.1f}%")
