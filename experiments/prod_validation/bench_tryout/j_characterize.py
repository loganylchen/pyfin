#!/usr/bin/env python3
"""What ARE the 'j' FPs, and what co-occurs with them? For the broad m2_em set:
(1) classify each 'j' by how its intron chain differs from the ref transcript it is
    'j' to: subset (exon-skip/partial), superset (extra intron), or shares-some;
(2) for each 'j', check whether pyfin ALSO produced a '=' candidate for the SAME ref
    at the same locus -> if so the 'j' is a redundant shadow (safe to drop, correct
    answer already present); if not, dropping it loses that region.
Also report read support (num_reads) of the 'j' vs the co-located '='.
"""
import re
from collections import defaultdict, Counter

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
S = "SGNex_H9_directRNA_replicate2_run2"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
GTF = f"{HERE}/prodfull_ladder_m2_em/p00/pyfin.gtf"
TRK = f"{HERE}/gc_ladder_m2_em/{S}__p00__pyfin.tracking"
SC = f"{HERE}/prodfull_ladder_m2_em/p00/work/scores.unfiltered.tsv"
ENST = re.compile(r"(ENST\d+)")

# candidate_id -> (class, ref_enst)
cls = {}
for ln in open(TRK):
    c = ln.rstrip("\n").split("\t")
    if len(c) < 5: continue
    code = c[3]; ref = ENST.search(c[2]); ref = ref.group(1) if ref else None
    for q in c[4:]:
        m = re.search(r":([^|]+)\|", q)
        if m: cls[m.group(1)] = (code, ref)

# pyfin candidate chains
tid = re.compile(r'transcript_id "([^"]+)"')
pex = defaultdict(list)
for ln in open(GTF):
    if ln.startswith("#"): continue
    f = ln.split("\t")
    if len(f) < 9 or f[2] != "exon": continue
    m = tid.search(f[8])
    if m: pex[m.group(1)].append((int(f[3]) - 1, int(f[4])))
pchain = {}
for t, ex in pex.items():
    ex.sort(); pchain[t] = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))

# truth chains for refs we need
need = {ref for code, ref in cls.values() if ref}
tex = defaultdict(list)
for ln in open(TRUTH):
    if ln.startswith("#"): continue
    f = ln.split("\t")
    if len(f) < 9 or f[2] != "exon": continue
    m = ENST.search(f[8])
    if not m or m.group(1) not in need: continue
    tex[m.group(1)].append((int(f[3]) - 1, int(f[4])))
tchain = {}
for t, ex in tex.items():
    ex.sort(); tchain[t] = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))

# num_reads per candidate
nr = {}
with open(SC) as fh:
    h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
    for ln in fh:
        c = ln.rstrip("\n").split("\t")
        nr[c[ci["candidate_id"]]] = int(float(c[ci["num_reads"]]))

# ref -> set of pyfin class codes pointing at it (co-occurrence)
ref_codes = defaultdict(set)
for cid, (code, ref) in cls.items():
    if ref: ref_codes[ref].add(code)

def jset(ch):  # junctions as a set (coordinate-exact)
    return set(ch)

diff_kinds = Counter()
cooc = Counter()
j_nr = []; eq_nr = []
n_j = 0
for cid, (code, ref) in cls.items():
    if code != "j" or not ref or ref not in tchain or cid not in pchain:
        continue
    n_j += 1
    jc = jset(pchain[cid]); rc = jset(tchain[ref])
    shared = jc & rc
    j_only = jc - rc; r_only = rc - jc
    # classify structural difference
    if jc < rc:
        diff_kinds["subset_of_ref (exon-skip/partial: fewer junctions, all real)"] += 1
    elif jc > rc:
        diff_kinds["superset_of_ref (extra junction vs ref)"] += 1
    elif not j_only:
        diff_kinds["equal-set?? (shouldn't be j)"] += 1
    else:
        # shares some, has some novel junctions not in this ref
        if len(j_only) == 1:
            diff_kinds["shares-most, 1 novel junction"] += 1
        else:
            diff_kinds[f"shares {len(shared)}, {len(j_only)} novel junctions"] += 1
    # co-occurrence with a '=' for the same ref
    has_eq = "=" in ref_codes.get(ref, set())
    cooc["ref ALSO has a '=' pyfin candidate (redundant j)" if has_eq
         else "ref has NO '=' (j is the only shot at this region)"] += 1
    j_nr.append(nr.get(cid, 0))
    if has_eq:
        # read support of the co-located '=' for this ref
        for cid2, (c2, r2) in cls.items():
            if c2 == "=" and r2 == ref:
                eq_nr.append(nr.get(cid2, 0)); break

import statistics
print(f"total 'j' FPs analyzed: {n_j}\n")
print("STRUCTURAL: how each 'j' differs from its ref transcript's true chain:")
for k, v in diff_kinds.most_common():
    print(f"  {v:6} ({100*v/n_j:4.0f}%)  {k}")
print("\nCO-OCCURRENCE: is the correct answer ('=') already present at the same ref?")
for k, v in cooc.most_common():
    print(f"  {v:6} ({100*v/n_j:4.0f}%)  {k}")
print(f"\nREAD SUPPORT (num_reads):")
print(f"  'j' FPs:            median {statistics.median(j_nr):.0f}  mean {statistics.mean(j_nr):.1f}")
if eq_nr:
    print(f"  co-located '=' TP: median {statistics.median(eq_nr):.0f}  mean {statistics.mean(eq_nr):.1f}")
    print(f"  -> if 'j' has far fewer reads than the co-located '=', a relative "
          f"read-support drop is safe.")
