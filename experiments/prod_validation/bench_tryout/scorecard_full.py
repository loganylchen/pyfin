#!/usr/bin/env python3
"""Single-sample tryout: score pyfin (production Lever 2) with the benchmark's OWN
procedures and put it head-to-head with the 6 competitors on the same sample.

- gffcompare vs the pristine GENCODE v44 truth (same image 0.12.10, same -r).
- orphan% @full and corrected_recall>=3 @full via the benchmark's nanocount, using
  the exact matched-ENST + est_count logic from run_orphan.py.
- corruption worst ΔF1 = min over c_* of (F1@c - F1@full).
Competitor numbers come from the benchmark's own gffcompare_metrics.tsv (per-cell)
and their existing .tracking files; pyfin's are computed here from its on.gtf.
"""
import os, re, subprocess, glob
B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
CGC = f"{B}/analysis/gencode_full_sweep/scoring/gffcompare"       # competitor tracking
METRICS = f"{B}/analysis/gencode_full_sweep/scoring/gffcompare_metrics.tsv"
IMG = "quay.io/biocontainers/gffcompare:0.12.10--h9948957_0"
PYF = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/prodfull"
SC = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/gc_full"
os.makedirs(SC, exist_ok=True)

S = "SGNex_H9_directRNA_replicate2_run2"
RATIOS = ["full","p99","p90","p50","p10","p00",
          "c_skip20","c_jitter20_10bp","c_spurious20","c_merge20","c_flip20","c_ir20"]
CORR = ["c_skip20","c_jitter20_10bp","c_spurious20","c_merge20","c_flip20","c_ir20"]
COMP = ["bambu","isoquant","isotools","lafite","stringtie3","talon"]
ENST_RE = re.compile(r"(ENST\d+\.\d+)")

def parse_stats(p):
    if not os.path.exists(p): return None
    t = open(p).read()
    def g(l):
        m = re.search(rf"{l} level:\s*([\d.]+)\s*\|\s*([\d.]+)", t)
        return (float(m.group(1)), float(m.group(2))) if m else (None, None)
    sn, sp = g("Transcript")
    return sn, sp

def f1(sn, sp):
    return 0.0 if not sn or not sp or (sn+sp)==0 else round(2*sn*sp/(sn+sp),2)

def pyfin_gtf(ratio):
    p = f"{PYF}/{ratio}/pyfin.gtf"
    return p if os.path.exists(p) else None

def gff_pyfin(ratio):
    """Run gffcompare for pyfin at a ratio; return (sn,sp,tracking_path)."""
    gtf = pyfin_gtf(ratio)
    if not gtf: return (None, None, None)
    pre = f"{SC}/{S}__{ratio}__pyfin"
    stats = pre + ".stats"; trk = pre + ".tracking"
    if not os.path.exists(stats):
        subprocess.run(["docker","run","--rm","-u",f"{os.getuid()}:{os.getgid()}",
            "-v","/autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD","-v","/SSD:/SSD",
            IMG,"gffcompare","-r",TRUTH,"-o",pre,gtf], capture_output=True)
    sn, sp = (parse_stats(stats) or (None, None))
    return sn, sp, (trk if os.path.exists(trk) else None)

# --- competitor per-cell F1 from their metrics table ---
comp_f1 = {}  # (ratio,tool)->(sn,sp,f1)
with open(METRICS) as fh:
    h = fh.readline().rstrip("\n").split("\t"); ci = {c:i for i,c in enumerate(h)}
    for ln in fh:
        c = ln.rstrip("\n").split("\t")
        if c[ci["sample"]] != S: continue
        comp_f1[(c[ci["ratio"]], c[ci["tool"]])] = (
            float(c[ci["sn_tx"]]), float(c[ci["sp_tx"]]), float(c[ci["f1_tx"]]))

# --- nanocount (expressed truth) ---
counts = {}
nc = f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
with open(nc) as fh:
    next(fh, None)
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 3: continue
        try: counts[p[0].split("|")[0].split(".")[0]] = float(p[2])
        except ValueError: pass
expr3 = {k for k,v in counts.items() if v >= 3}

def matched_base(trk):
    s = set()
    if not trk or not os.path.exists(trk): return s
    for ln in open(trk):
        c = ln.rstrip("\n").split("\t")
        if len(c) < 4 or c[3] != "=": continue
        m = ENST_RE.search(c[2])
        if m: s.add(m.group(1).split(".")[0])
    return s

def orphan_and_recall(trk):
    mb = matched_base(trk); n = len(mb)
    if n == 0: return (None, None)
    c0 = sum(1 for e in mb if counts.get(e, 0.0) < 1)
    orphan = round(100.0*c0/n, 1)
    rec = round(100.0*len(mb & expr3)/max(len(expr3),1), 1)
    return orphan, rec

# --- build the scorecard ---
rows = []  # tool, F1@full, Pr@full, orphan%, corr_recall3, worst_dF1
for tool in ["pyfin"] + COMP:
    fulls = {}
    for r in RATIOS:
        if tool == "pyfin":
            sn, sp, trk = gff_pyfin(r); fulls[r] = (sn, sp, f1(sn, sp), trk)
        else:
            v = comp_f1.get((r, tool))
            trk = f"{CGC}/{S}__{r}__{tool}.tracking"
            fulls[r] = (v[0], v[1], v[2], trk) if v else (None, None, None, trk)
    fsn, fsp, ff1, ftrk = fulls["full"]
    orphan, rec = orphan_and_recall(ftrk)
    # only count corruptions whose gffcompare actually ran (Sn present) — missing != F1 0
    dfs = [fulls[c][2]-ff1 for c in CORR if fulls[c][0] is not None and ff1 is not None]
    worst = round(min(dfs), 1) if dfs else None
    rows.append((tool, ff1, fsp, orphan, rec, worst))

print(f"# Single-sample tryout: {S}  (pyfin = production Lever 2, vs 6 competitors)")
print(f"# expressed truth (nanocount est>=3): {len(expr3)} transcripts\n")
print(f"{'tool':12s} {'F1@full':>8} {'Pr@full':>8} {'orphan%':>8} {'corrRec>=3':>11} {'worstΔF1':>9}")
for t, ff1, fsp, orph, rec, worst in sorted(rows, key=lambda x:-(x[4] or 0)):
    print(f"{t:12s} {str(ff1):>8} {str(fsp):>8} {str(orph):>8} {str(rec):>11} {str(worst):>9}")
