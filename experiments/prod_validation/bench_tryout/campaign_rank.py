#!/usr/bin/env python3
"""Single-sample 6-axis composite ranking (mirrors the benchmark's derive_tables
composite, but computed on ONE sample so pyfin is apples-to-apples with the 6
competitors). Axes/weights: F1@full .25, corrRec>=3 .20, precision(sp_tx@full)
.15, (1-orphan%) .15, sweep_auc .15, worst-corruption-dF1 .10. Min-max normalized
across the 7 tools, weighted x100, ranked.

Usage: campaign_rank.py [PYFIN_GTF_DIR] [SCORE_CACHE_DIR]
  PYFIN_GTF_DIR/<ratio>/pyfin.gtf  (default: prodfull)   -- the lever's output
Reuses competitor per-sample F1/Sn/Sp from the benchmark metrics TSV + their
tracking files for orphan/corrRec.
"""
import os, re, sys, subprocess
B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
METRICS = f"{B}/analysis/gencode_full_sweep/scoring/gffcompare_metrics.tsv"
CGC = f"{B}/analysis/gencode_full_sweep/scoring/gffcompare"
NC = f"{B}/results/gencode_full_sweep/SGNex_H9_directRNA_replicate2_run2/full/assembly/nanocount.tsv"
IMG = "quay.io/biocontainers/gffcompare:0.12.10--h9948957_0"
S = "SGNex_H9_directRNA_replicate2_run2"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
COMPLETE = ["full", "p99", "p90", "p50", "p10", "p00"]
KEEP = {"full":100,"p99":99,"p90":90,"p50":50,"p10":10,"p00":0}
CORRUPT = ["c_skip20","c_jitter20_10bp","c_spurious20","c_merge20","c_flip20","c_ir20"]
COMP = ["bambu","talon","stringtie3","isoquant","lafite","isotools"]
ENST = re.compile(r"(ENST\d+)")

pyf_dir = sys.argv[1] if len(sys.argv) > 1 else f"{HERE}/prodfull"
scdir = sys.argv[2] if len(sys.argv) > 2 else f"{HERE}/gc_rank"
os.makedirs(scdir, exist_ok=True)

def f1(sn, sp): return 0.0 if not sn or not sp or (sn+sp)==0 else 2*sn*sp/(sn+sp)

# expressed truth
est = {}
with open(NC) as fh:
    next(fh, None)
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if len(p) < 3: continue
        m = ENST.search(p[0])
        if m:
            try: est[m.group(1)] = float(p[2])
            except ValueError: pass
expr = {k for k,v in est.items() if v >= 3}

def gff(gtf, pre):
    if not os.path.exists(pre+".stats"):
        subprocess.run(["docker","run","--rm","-u",f"{os.getuid()}:{os.getgid()}",
            "-v","/autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD","-v","/SSD:/SSD",
            IMG,"gffcompare","-r",TRUTH,"-o",pre,gtf], capture_output=True)
    return pre+".stats", pre+".tracking"

def parse(stats):
    if not os.path.exists(stats): return (None, None)
    t = open(stats).read(); m = re.search(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)", t)
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)

def orphan_recall(trk):
    if not trk or not os.path.exists(trk): return (None, None)
    mb = set()
    for ln in open(trk):
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 4 and c[3] == "=":
            m = ENST.search(c[2])
            if m: mb.add(m.group(1))
    if not mb: return (None, None)
    orphan = 100.0*sum(1 for e in mb if est.get(e,0.0) < 1)/len(mb)
    rec = 100.0*len(mb & expr)/len(expr)
    return orphan, rec

# competitor per-(ratio,tool) F1/Sn/Sp for THIS sample
cm = {}
with open(METRICS) as fh:
    h = fh.readline().rstrip("\n").split("\t"); ci = {c:i for i,c in enumerate(h)}
    for ln in fh:
        c = ln.rstrip("\n").split("\t")
        if c[ci["sample"]] != S: continue
        cm[(c[ci["ratio"]], c[ci["tool"]])] = (float(c[ci["sn_tx"]]), float(c[ci["sp_tx"]]), float(c[ci["f1_tx"]]))

def f1_of(tool, r):
    if tool == "pyfin":
        g = f"{pyf_dir}/{r}/pyfin.gtf"
        if not os.path.exists(g): return None
        st, _ = gff(g, f"{scdir}/{S}__{r}__pyfin"); sn, sp = parse(st)
        return None if sn is None else (sn, sp, f1(sn, sp))
    return cm.get((r, tool))

def axes_for(tool):
    full = f1_of(tool, "full")
    if not full: return None
    f1full = full[2]; prec = full[1]
    # sweep auc
    pts = sorted([(KEEP[r], v[2]) for r in COMPLETE if (v := f1_of(tool, r))])
    auc = float("nan")
    if len(pts) >= 2:
        area = sum((pts[i][0]-pts[i-1][0])*(pts[i][1]+pts[i-1][1])/2 for i in range(1,len(pts)))
        auc = area/(pts[-1][0]-pts[0][0])
    # worst corruption dF1
    ds = [v[2]-f1full for c in CORRUPT if (v := f1_of(tool, c))]
    worst = min(ds) if ds else float("nan")
    # orphan + corrRec
    if tool == "pyfin":
        _, trk = gff(f"{pyf_dir}/full/pyfin.gtf", f"{scdir}/{S}__full__pyfin")
    else:
        trk = f"{CGC}/{S}__full__{tool}.tracking"
    orph, corr3 = orphan_recall(trk)
    clean = None if orph is None else 100.0 - orph
    return {"f1_full":f1full, "corr3":corr3, "prec":prec, "clean":clean, "auc":auc, "robust":worst}

tools = ["pyfin"] + COMP
raw = {t: a for t in tools if (a := axes_for(t))}
AX = ["f1_full","corr3","prec","clean","auc","robust"]
W = {"f1_full":.25,"corr3":.20,"prec":.15,"clean":.15,"auc":.15,"robust":.10}
norm = {}
for a in AX:
    vals = [raw[t][a] for t in raw if raw[t][a] is not None and raw[t][a]==raw[t][a]]
    lo, hi = (min(vals), max(vals)) if vals else (0.0, 1.0)
    for t in raw:
        v = raw[t][a]
        norm.setdefault(t, {})[a] = 0.0 if (v is None or v!=v or hi==lo) else (v-lo)/(hi-lo)
rows = []
for t in raw:
    sc = 100.0*sum(W[a]*norm[t][a] for a in AX)
    r = raw[t]
    rows.append((t, sc, r["f1_full"], r["corr3"], r["prec"], r["clean"], r["auc"], r["robust"]))
rows.sort(key=lambda x: -x[1])

print(f"# CAMPAIGN composite rank — {S} | pyfin_dir={os.path.basename(pyf_dir)}")
print(f"{'rank tool':16} {'composite':>9} {'F1@full':>8} {'corrRec':>8} {'prec':>6} {'clean':>6} {'sweepAUC':>9} {'worstdF1':>9}")
for i,(t,sc,ff,cr,pr,cl,au,ro) in enumerate(rows,1):
    star = " <== PYFIN" if t=="pyfin" else ""
    def s(x): return "  NA" if x is None or x!=x else f"{x:.1f}"
    print(f"{i:>2} {t:12} {sc:>9.1f} {s(ff):>8} {s(cr):>8} {s(pr):>6} {s(cl):>6} {s(au):>9} {s(ro):>9}{star}")
