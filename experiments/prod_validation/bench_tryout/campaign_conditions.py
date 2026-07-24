#!/usr/bin/env python3
"""Per-condition pyfin-vs-competitors matrices (F1/Sn/Pr) across all 12 conditions:
6 completeness ratios (full..p00) and 6 corruption operators (+dF1 vs full).
Complements campaign_rank.py (which gives the single composite number).

Usage: campaign_conditions.py [PYFIN_GTF_DIR] [SCORE_CACHE_DIR]
"""
import os, re, sys, subprocess
B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
METRICS = f"{B}/analysis/gencode_full_sweep/scoring/gffcompare_metrics.tsv"
IMG = "quay.io/biocontainers/gffcompare:0.12.10--h9948957_0"
S = "SGNex_H9_directRNA_replicate2_run2"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
RATIOS = ["full","p99","p90","p50","p10","p00"]
CORRUPT = ["c_skip20","c_jitter20_10bp","c_spurious20","c_merge20","c_flip20","c_ir20"]
COMP = ["bambu","talon","stringtie3","isoquant","lafite","isotools"]
ORDER = ["pyfin"] + COMP
pyf_dir = sys.argv[1] if len(sys.argv) > 1 else f"{HERE}/prodfull"
scdir = sys.argv[2] if len(sys.argv) > 2 else f"{HERE}/gc_full"
os.makedirs(scdir, exist_ok=True)

def f1(sn, sp): return 0.0 if not sn or not sp or (sn+sp)==0 else 2*sn*sp/(sn+sp)
def parse(stats):
    if not os.path.exists(stats): return (None, None)
    t = open(stats).read(); m = re.search(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)", t)
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)
def gff(gtf, pre):
    if not os.path.exists(pre+".stats"):
        subprocess.run(["docker","run","--rm","-u",f"{os.getuid()}:{os.getgid()}",
            "-v","/autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD","-v","/SSD:/SSD",
            IMG,"gffcompare","-r",TRUTH,"-o",pre,gtf], capture_output=True)
    return pre+".stats"

cm = {}
with open(METRICS) as fh:
    h = fh.readline().rstrip("\n").split("\t"); ci = {c:i for i,c in enumerate(h)}
    for ln in fh:
        c = ln.rstrip("\n").split("\t")
        if c[ci["sample"]] != S: continue
        cm[(c[ci["ratio"]], c[ci["tool"]])] = (float(c[ci["sn_tx"]]), float(c[ci["sp_tx"]]), float(c[ci["f1_tx"]]))

def triple(tool, r):
    if tool == "pyfin":
        g = f"{pyf_dir}/{r}/pyfin.gtf"
        if not os.path.exists(g): return None
        sn, sp = parse(gff(g, f"{scdir}/{S}__{r}__pyfin"))
        return None if sn is None else (sn, sp, f1(sn, sp))
    return cm.get((r, tool))

def matrix(title, cols, which, delta=False):
    print(f"\n### {title}")
    hdr = f"{'tool':11}" + "".join(f"{c.replace('c_','').replace('20_10bp','jit').replace('20',''):>9}" for c in cols)
    print(hdr)
    for t in ORDER:
        full = triple(t, "full")
        cells = []
        for c in cols:
            v = triple(t, c)
            if v is None: cells.append(f"{'-':>9}"); continue
            val = v[which]
            if delta and full: val = v[2] - full[2]
            cells.append(f"{val:>9.1f}")
        mark = " <=PYFIN" if t == "pyfin" else ""
        print(f"{t:11}" + "".join(cells) + mark)

print(f"# PER-CONDITION — {S} | pyfin_dir={os.path.basename(pyf_dir)}  (which=2:F1)")
print("## COMPLETENESS RATIOS (keep X% of annotation)")
matrix("F1", RATIOS, 2)
matrix("Sn (recall)", RATIOS, 0)
matrix("Pr (precision)", RATIOS, 1)
print("\n## CORRUPTIONS (20% of annotation damaged)")
matrix("F1", CORRUPT, 2)
matrix("dF1 vs full (robustness; 0=immune)", CORRUPT, 2, delta=True)
matrix("Pr (precision) — jitter is the pyfin weak spot", CORRUPT, 1)
