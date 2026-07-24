#!/usr/bin/env python3
"""Scan every (sample, ratio, config) cell for cassette-skip candidate pairs.

A cassette-skip pair (A, B) is:
  - |#introns_A - #introns_B| = 1
  - one chain (A) has K introns, the other (B) has K-1
  - aligned chain-wise: all positions match within MAX_WOBBLE bp except a
    single position where A has 2 introns flanking a small exon (size < MAX_EXON_BP)
    and B has 1 long intron spanning the same span (within wobble)
  - BOTH have hard read support (num_reads >= 1)

Output: tables/cassette_pairs.tsv with per-cell counts + total / hard-supported.
"""
import os, re, glob, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
MATRIX = os.path.join(ROOT, "matrix")
MAX_WOBBLE = 20
MAX_EXON_BP = 70

SAMPLES = sorted(d for d in os.listdir(MATRIX)
                 if d.startswith("SGNex_HEYA8") and os.path.isdir(os.path.join(MATRIX, d)))
RATIOS  = ["full", "p00", "c_jitter10bp", "c_skip10", "c_spurious5",
           "c_merge5", "c_flip5", "c_ir10"]
CONFIGS = ["bp10", "bp20", "bp20cas"]

def parse_chains(path):
    tx = {}
    if not os.path.exists(path): return {}
    for line in open(path):
        if line.startswith("#"): continue
        p = line.rstrip().split("\t")
        if len(p) < 9 or p[2] not in ("exon", "transcript"): continue
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if not m: continue
        t = m.group(1)
        d = tx.setdefault(t, {"c": p[0], "s": p[6], "ex": [], "nr": 0, "ab": 0.0, "src": ""})
        if p[2] == "transcript":
            for k in ("transcript_source", "abundance", "num_reads"):
                mm = re.search(fr'{k} "([^"]+)"', p[8])
                if mm:
                    val = mm.group(1)
                    if k == "transcript_source":   d["src"] = val
                    elif k == "abundance":         d["ab"]  = float(val or 0)
                    elif k == "num_reads":         d["nr"]  = int(float(val or 0))
        else:
            d["ex"].append((int(p[3])-1, int(p[4])))
    out = {}
    for t, d in tx.items():
        ex = sorted(d["ex"])
        intr = tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))
        out[t] = {"chrom": d["c"], "strand": d["s"], "introns": intr,
                  "nr": d["nr"], "ab": d["ab"], "src": d["src"]}
    return out

def cassette_diff(A, B, max_exon=MAX_EXON_BP, max_wobble=MAX_WOBBLE):
    """A has K introns, B has K-1. Returns (i, cassette_size) or None."""
    if len(A) != len(B) + 1: return None
    for i in range(len(A) - 1):
        if any(abs(A[k][0]-B[k][0])>max_wobble or abs(A[k][1]-B[k][1])>max_wobble
               for k in range(i)): continue
        if abs(A[i][0]-B[i][0])>max_wobble:    continue
        if abs(A[i+1][1]-B[i][1])>max_wobble: continue
        rest = len(B) - i - 1
        if any(abs(A[i+2+k][0]-B[i+1+k][0])>max_wobble or abs(A[i+2+k][1]-B[i+1+k][1])>max_wobble
               for k in range(rest)): continue
        cas = A[i+1][0] - A[i][1]
        if 0 < cas < max_exon:
            return (i, cas)
    return None

def scan_cell(gtf_path):
    cands = parse_chains(gtf_path)
    if not cands:
        return {"n_tx": 0, "pairs_total": 0, "pairs_with_hard": 0,
                "pairs_both_novel": 0, "pairs_both_novel_with_hard": 0}
    by_loc = defaultdict(list)
    for t, d in cands.items():
        by_loc[(d["chrom"], d["strand"])].append((t, d))
    n_total = 0; n_hard = 0
    n_bn = 0;    n_bn_hard = 0
    for (chrom, strand), items in by_loc.items():
        by_n = defaultdict(list)
        for t, d in items: by_n[len(d["introns"])].append((t, d))
        for n, A_items in by_n.items():
            for ta, da in A_items:
                for tb, db in by_n.get(n-1, []):
                    r = cassette_diff(da["introns"], db["introns"])
                    if r is None: continue
                    n_total += 1
                    both_hard = da["nr"] >= 1 and db["nr"] >= 1
                    both_novel = da["src"] == "novel" and db["src"] == "novel"
                    if both_hard:  n_hard += 1
                    if both_novel: n_bn += 1
                    if both_hard and both_novel: n_bn_hard += 1
    return {"n_tx": len(cands), "pairs_total": n_total, "pairs_with_hard": n_hard,
            "pairs_both_novel": n_bn, "pairs_both_novel_with_hard": n_bn_hard}

out = os.path.join(ROOT, "tables", "cassette_pairs.tsv")
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as fh:
    fh.write("sample\tratio\tconfig\tn_tx\tcassette_pairs_total\twith_hard_reads\tboth_novel\tboth_novel_with_hard\n")
    for s in SAMPLES:
        for r in RATIOS:
            for c in CONFIGS:
                gtf = os.path.join(MATRIX, s, r, f"out_{c}", "assembly.gtf")
                stats = scan_cell(gtf)
                fh.write(f"{s}\t{r}\t{c}\t{stats['n_tx']}\t{stats['pairs_total']}"
                         f"\t{stats['pairs_with_hard']}\t{stats['pairs_both_novel']}"
                         f"\t{stats['pairs_both_novel_with_hard']}\n")
print(f"wrote {out}")

# Quick markdown summary highlighting cells where the rule WOULD matter
md = os.path.join(ROOT, "tables", "cassette_summary.md")
with open(md, "w") as fh:
    fh.write("# Cassette-pair scan (per cell)\n\n")
    fh.write("A cassette pair = two candidates differing by ONE small exon (<70bp), all other junctions within 20bp.\n")
    fh.write("`both_novel_with_hard` is the key column: if non-zero, the cassette-cluster rule WOULD remove FP that current filters miss.\n\n")
    fh.write("| sample | ratio | config | nout | total | hard | both_novel | both_novel_hard |\n")
    fh.write("|---|---|---|---|---|---|---|---|\n")
    # re-read TSV to sort/print
    with open(out) as src:
        next(src)
        for ln in src:
            f = ln.rstrip().split("\t")
            mark = " ⚠️" if int(f[7]) > 0 else ""
            fh.write(f"| {f[0][-9:]} | {f[1]} | {f[2]} | {f[3]} | {f[4]} | {f[5]} | {f[6]} | **{f[7]}**{mark} |\n")
print(f"wrote {md}")
