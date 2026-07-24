#!/usr/bin/env python3
"""Premise check for a krill-signal junction-validation gate (NO signal needed).

A per-junction signal-NLL gate can only drop candidates that carry a FAKE junction
(a splice site not in any real transcript -> minimap wobble / misalignment). It
CANNOT drop a truncation (gffcompare 'c') or a recombination of REAL junctions,
because those junctions have genuine signal support.

So: classify every pyfin p00 novel candidate by gffcompare class code, and for each
count how many of its junctions are absent from the truth junction set (a "fake"
junction) within a small tolerance. If wrong candidates ('c','j','i','o','x','?')
are mostly built from REAL junctions, a signal gate has little to catch.
"""
import os, re, sys
from collections import defaultdict

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
S = "SGNex_H9_directRNA_replicate2_run2"
PYGTF = f"{HERE}/prodfull/p00/pyfin.gtf"
TRACK = f"{HERE}/gc_full/{S}__p00__pyfin.tracking"
TOL = int(sys.argv[1]) if len(sys.argv) > 1 else 5   # bp tol for "real junction"


def introns_by_tx(path, want_tid_re=re.compile(r'transcript_id "([^"]+)"')):
    exons = defaultdict(list)
    meta = {}
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = want_tid_re.search(f[8])
        if not m:
            continue
        t = m.group(1)
        exons[t].append((int(f[3]) - 1, int(f[4])))
        meta[t] = (f[0], f[6])
    chains = {}
    for t, ex in exons.items():
        ex.sort()
        chains[t] = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
    return chains, meta


# truth junction set, keyed by chrom for tolerance matching
truth_chains, truth_meta = introns_by_tx(TRUTH)
truth_donors = defaultdict(set)
truth_accept = defaultdict(set)
truth_junc = defaultdict(set)
for t, ch in truth_chains.items():
    chrom = truth_meta[t][0]
    for d, a in ch:
        truth_junc[chrom].add((d, a))
        truth_donors[chrom].add(d)
        truth_accept[chrom].add(a)


def is_real(chrom, j):
    d, a = j
    if (d, a) in truth_junc[chrom]:
        return True
    # within tol on both sides of SOME real junction
    for rd in range(d - TOL, d + TOL + 1):
        if rd in truth_donors[chrom]:
            for ra in range(a - TOL, a + TOL + 1):
                if ra in truth_accept[chrom] and (rd, ra) in truth_junc[chrom]:
                    return True
    return False


# pyfin candidate chains
py_chains, py_meta = introns_by_tx(PYGTF)

# map pyfin transcript_id -> class code from tracking qData
code_by_tid = {}
for ln in open(TRACK):
    c = ln.rstrip("\n").split("\t")
    if len(c) < 5:
        continue
    code = c[3]
    for q in c[4:]:
        m = re.search(r"\|([^|]+)\|", q)  # q1:gene|transcript|...
        if m:
            code_by_tid[m.group(1)] = code

# aggregate: per class code, how many candidates have >=1 fake junction
stats = defaultdict(lambda: [0, 0, 0, 0])  # [n_cands, n_multiexon, n_with_fake_junc, total_fake_junc]
for t, ch in py_chains.items():
    code = code_by_tid.get(t, "?")
    chrom = py_meta[t][0]
    st = stats[code]
    st[0] += 1
    if not ch:
        continue
    st[1] += 1
    fakes = sum(0 if is_real(chrom, j) else 1 for j in ch)
    if fakes:
        st[2] += 1
    st[3] += fakes

order = ["=", "c", "j", "k", "m", "n", "e", "o", "i", "x", "p", "s", "?"]
print(f"pyfin p00 novel candidates — fake-junction analysis (tol=+/-{TOL}bp)")
print(f"{'code':5}{'n_cand':>8}{'multiexon':>10}{'w/fake_j':>9}{'%fake':>7}{'fakeJ':>7}")
for code in order + [c for c in stats if c not in order]:
    if code not in stats:
        continue
    n, me, wf, tf = stats[code]
    pct = 100 * wf / me if me else 0
    print(f"{code:5}{n:>8}{me:>10}{wf:>9}{pct:>6.0f}%{tf:>7}")

# the decisive numbers
wrong = ["c", "j", "k", "m", "n", "e", "o", "i", "x", "p", "s", "?"]
wr_me = sum(stats[c][1] for c in wrong if c in stats)
wr_wf = sum(stats[c][2] for c in wrong if c in stats)
cor_me = stats["="][1] if "=" in stats else 0
cor_wf = stats["="][2] if "=" in stats else 0
print(f"\nCORRECT (=)   multiexon={cor_me}  with>=1 fake junction={cor_wf} "
      f"({100*cor_wf/cor_me if cor_me else 0:.0f}%)")
print(f"WRONG (rest)  multiexon={wr_me}  with>=1 fake junction={wr_wf} "
      f"({100*wr_wf/wr_me if wr_me else 0:.0f}%)")
print(f"\n=> Signal gate can only catch WRONG candidates that carry a fake junction: "
      f"{wr_wf} of {wr_me} wrong multi-exon.")
print(f"   Truncations 'c' with ALL-real junctions (signal-invisible): "
      f"{stats['c'][1]-stats['c'][2] if 'c' in stats else 0} of {stats['c'][1] if 'c' in stats else 0}")
