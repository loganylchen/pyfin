#!/usr/bin/env python3
"""Break down OLD's '=' advantage over NEW:
  (1) how much is pure REDUNDANCY (multiple '=' candidates matching one truth), and
  (2) of the distinct truth OLD matches but NEW misses, WHY did NEW miss each one
      (exact-subset-collapsed away / present as a wobbled non-'=' member / truly absent).
"""
import re
from collections import defaultdict

HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
OUT = f"{HERE}/cluster_cmp"
B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
ENST = re.compile(r"(ENST\d+)")


def eq_refs_and_count(trk):
    """distinct '=' refs and raw '=' candidate count."""
    refs, n = set(), 0
    for ln in open(trk):
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 4 and c[3] == "=":
            n += 1
            m = ENST.search(c[2])
            if m:
                refs.add(m.group(1))
    return refs, n


old_refs, old_n = eq_refs_and_count(f"{OUT}/gc_OLD.tracking")
new_refs, new_n = eq_refs_and_count(f"{OUT}/gc_NEW.tracking")
new_misses = old_refs - new_refs

print("=== OLD '=' advantage breakdown ===")
print(f"OLD '=' candidates: {old_n}  -> distinct truth: {len(old_refs)}  "
      f"(redundant duplicates: {old_n - len(old_refs)})")
print(f"NEW '=' candidates: {new_n}  -> distinct truth: {len(new_refs)}  "
      f"(redundant duplicates: {new_n - len(new_refs)})")
print(f"OLD raw '=' lead: {old_n - new_n}  =  "
      f"{old_n - len(old_refs)} OLD-duplicates + {len(new_misses)} distinct NEW-misses "
      f"- {new_n - len(new_refs)} NEW-duplicates")
print(f"\ndistinct truth OLD matches but NEW misses: {len(new_misses)}")

# --- classify WHY NEW missed each of the new_misses ---
# truth chains for the missed refs
tex = defaultdict(list)
tid = re.compile(r'transcript_id "([^"]+)"')
for ln in open(TRUTH):
    if ln.startswith("#"):
        continue
    f = ln.split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    m = ENST.search(f[8])
    if not m or m.group(1) not in new_misses:
        continue
    tex[m.group(1)].append((f[0], int(f[3]) - 1, int(f[4])))
truth_chain = {}
for e, ex in tex.items():
    ex.sort(key=lambda x: x[1])
    chrom = ex[0][0]
    introns = tuple((ex[i][2], ex[i + 1][1]) for i in range(len(ex) - 1))
    truth_chain[e] = (chrom, introns)

# NEW candidate chains, indexed by chrom
new_ex = defaultdict(list)
for ln in open(f"{OUT}/NEW.gtf"):
    f = ln.split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    m = tid.search(f[8])
    if m:
        new_ex[m.group(1)].append((f[0], int(f[3]) - 1, int(f[4])))
new_by_chrom = defaultdict(list)
for t, ex in new_ex.items():
    ex.sort(key=lambda x: x[1])
    chrom = ex[0][0]
    introns = tuple((ex[i][2], ex[i + 1][1]) for i in range(len(ex) - 1))
    if introns:
        new_by_chrom[chrom].append(introns)


def exact_sub(short, long_):
    n, mm = len(long_), len(short)
    if mm == 0 or mm >= n:
        return False
    return any(long_[o:o + mm] == short for o in range(n - mm + 1))


def wobble(a, b, bp=6):
    if len(a) != len(b):
        return False
    return all(abs(s1 - s2) <= bp and abs(e1 - e2) <= bp
               for (s1, e1), (s2, e2) in zip(a, b))


def tol_contains(short, long_, bp=6):
    n, mm = len(long_), len(short)
    if mm == 0 or mm >= n:
        return False
    for o in range(n - mm + 1):
        if all(abs(short[k][0] - long_[o + k][0]) <= bp and
               abs(short[k][1] - long_[o + k][1]) <= bp for k in range(mm)):
            return True
    return False


cat = defaultdict(int)
for e, (chrom, tch) in truth_chain.items():
    cands = new_by_chrom.get(chrom, [])
    if not tch:
        cat["mono (no introns)"] += 1
        continue
    if any(exact_sub(tch, nc) for nc in cands):
        cat["A. exact-subset-collapsed into a NEW candidate"] += 1
    elif any(wobble(tch, nc) for nc in cands):
        cat["B. present as a WOBBLED member (coords off -> not '=')"] += 1
    elif any(tol_contains(tch, nc) for nc in cands):
        cat["C. present as a tolerant-contained member (not '=')"] += 1
    else:
        cat["D. truly absent from NEW (no nearby candidate)"] += 1

tot = sum(cat.values())
print(f"\n=== WHY NEW missed the {tot} distinct truth (with chains) ===")
for k in sorted(cat):
    print(f"  {cat[k]:5} ({100*cat[k]/tot:4.0f}%)  {k}")
print("\nA = recoverable by relaxing exact-subset collapse; B/C = recoverable by "
      "snapping/keeping the true-coord chain; D = genuinely not generated.")
