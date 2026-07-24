#!/usr/bin/env python3
"""Root-cause split for the 1064 recoverable missed transcripts: did pyfin EMIT a
near-miss candidate (wobbled junctions -> fails '='; fix = read-consensus snap) or
emit NOTHING at that locus (dropped by a gate/grouping; fix = the gate)?

For each iso-only expressed transcript with >=2 exact full-length reads, check
whether pyfin's p00 GTF has a candidate whose intron chain equals the truth chain
within an increasing bp tolerance. tol=0 already '=' (shouldn't be, they're 'missed');
matched only at tol>0 => WOBBLE (case W). Never matched at any tol => DROPPED (case D).
"""
import os, re
import pysam
from collections import defaultdict

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
S = "SGNex_H9_directRNA_replicate2_run2"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
BAM = f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
NC = f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
PYGTF = f"{HERE}/prodfull/p00/pyfin.gtf"
ISO_TRK = f"{HERE}/gc_frontier/{S}__p00__isoquant.tracking"
PYF_TRK = f"{HERE}/gc_frontier/{S}__p00__pyfin.tracking"
ENST = re.compile(r"(ENST\d+)")
RTOL = 6


def matched_eq(trk):
    s = set()
    for ln in open(trk):
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 4 and c[3] == "=":
            m = ENST.search(c[2])
            if m: s.add(m.group(1))
    return s


def introns_by_tx(path, keep=None):
    tid = re.compile(r'transcript_id "([^"]+)"')
    exons = defaultdict(list); meta = {}
    for ln in open(path):
        if ln.startswith("#"): continue
        f = ln.split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = (ENST.search(f[8]) if keep is not None else tid.search(f[8]))
        if not m: continue
        t = m.group(1)
        if keep is not None and t not in keep: continue
        exons[t].append((int(f[3]) - 1, int(f[4])))
        meta[t] = (f[0], f[6])
    chains = {}
    for t, ex in exons.items():
        ex.sort()
        chains[t] = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
    return chains, meta


est = {}
for ln in open(NC):
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 3:
        m = ENST.search(p[0])
        if m:
            try: est[m.group(1)] = float(p[2])
            except ValueError: pass
expr3 = {k for k, v in est.items() if v >= 3}
iso_only = (matched_eq(ISO_TRK) - matched_eq(PYF_TRK)) & expr3

truth_ch, truth_meta = introns_by_tx(TRUTH, keep=iso_only)
truth_ch = {t: c for t, c in truth_ch.items() if c}   # multi-exon only

# pyfin candidates indexed by chrom for wobble search
py_ch, py_meta = introns_by_tx(PYGTF)
py_by_chrom = defaultdict(list)
for t, c in py_ch.items():
    if c:
        py_by_chrom[py_meta[t][0]].append(c)


def chain_match(a, b, tol):
    if len(a) != len(b):
        return False
    return all(abs(a[k][0] - b[k][0]) <= tol and abs(a[k][1] - b[k][1]) <= tol
               for k in range(len(a)))


# restrict to those with >=2 exact full-length reads (the 1064)
bam = pysam.AlignmentFile(BAM, "rb")
def read_introns(r):
    out = []; pos = r.reference_start
    for op, ln in r.cigartuples:
        if op in (0, 2, 7, 8): pos += ln
        elif op == 3: out.append((pos, pos + ln)); pos += ln
    return out

recoverable = []
for t, tch in truth_ch.items():
    chrom, strand = truth_meta[t]
    lo = min(j[0] for j in tch); hi = max(j[1] for j in tch)
    n = 0
    for r in bam.fetch(chrom, max(0, lo - 100), hi + 100):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        if chain_match(tuple(read_introns(r)), tch, RTOL):
            n += 1
    if n >= 2:
        recoverable.append(t)
bam.close()

# classify: does a pyfin candidate wobble-match the truth chain?
buckets = defaultdict(int)
for t in recoverable:
    tch = truth_ch[t]; chrom = truth_meta[t][0]
    cands = py_by_chrom.get(chrom, [])
    best = None
    for tol in (0, 2, 5, 10, 20):
        if any(chain_match(pc, tch, tol) for pc in cands):
            best = tol; break
    buckets["W_tol%s" % best if best is not None else "D_dropped"] += 1

print(f"iso-only expressed multi-exon: {len(truth_ch)}")
print(f"  with >=2 exact full-length reads (recoverable): {len(recoverable)}")
print("\nRoot-cause split of the recoverable set:")
tot = len(recoverable)
wob = sum(v for k, v in buckets.items() if k.startswith("W"))
drp = buckets.get("D_dropped", 0)
for k in sorted(buckets):
    print(f"  {k:12} {buckets[k]:5} ({100*buckets[k]/tot:.0f}%)")
print(f"\n  WOBBLE (pyfin emitted a near-miss cand, off by <=20bp) = {wob} ({100*wob/tot:.0f}%)"
      f"  -> fix = read-consensus junction snap (no annotation)")
print(f"  DROPPED (no pyfin candidate within 20bp)              = {drp} ({100*drp/tot:.0f}%)"
      f"  -> fix = candidate-generation / gate")
