#!/usr/bin/env python3
"""Per-condition honest metrics for ONE pyfin output dir, all 12 conditions:
standard F1, orphan%(est<3), honest precision (Pr*(1-orphan)), corrRec(est>=3),
honest F1(est>=3). Run on two dirs and diff to get OFF-vs-ON across every ratio
and corruption.

Usage: eval_honest.py [PYFIN_GTF_DIR] [SCORE_CACHE_DIR]
"""
import re, os, sys, subprocess
B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
IMG = "quay.io/biocontainers/gffcompare:0.12.10--h9948957_0"
S = "SGNex_H9_directRNA_replicate2_run2"
NC = f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
ENST = re.compile(r"(ENST\d+)")
CONDS = ["full","p99","p90","p50","p10","p00",
         "c_skip20","c_jitter20_10bp","c_spurious20","c_merge20","c_flip20","c_ir20"]
pyf = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else f"{HERE}/prodfull"
sc = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else f"{HERE}/gc_full"
os.makedirs(sc, exist_ok=True)

est = {}
for ln in open(NC):
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 3:
        m = ENST.search(p[0])
        if m:
            try: est[m.group(1)] = float(p[2])
            except ValueError: pass
expr3 = {k for k,v in est.items() if v >= 3}
def gff(g, pre):
    if not os.path.exists(pre+".stats"):
        subprocess.run(["docker","run","--rm","-u",f"{os.getuid()}:{os.getgid()}",
            "-v","/autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD","-v","/SSD:/SSD",
            IMG,"gffcompare","-r",TRUTH,"-o",pre,g], capture_output=True)
    return pre+".stats", pre+".tracking"
def parse(st):
    m = re.search(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)", open(st).read())
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)
def matched(trk):
    s = set()
    for ln in open(trk):
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 4 and c[3] == "=":
            m = ENST.search(c[2])
            if m: s.add(m.group(1))
    return s
def f1(a, b): return 0.0 if a+b == 0 else 2*a*b/(a+b)

print(f"# eval_honest — dir={os.path.basename(pyf)}  {S}")
print(f"{'cond':16} {'tx':>7} {'F1std':>7} {'orphan':>7} {'honPr':>7} {'corrRec':>8} {'honF1':>7}")
for c in CONDS:
    g = f"{pyf}/{c}/pyfin.gtf"
    if not os.path.exists(g):
        print(f"{c:16} {'MISSING':>7}"); continue
    st, trk = gff(g, f"{sc}/{S}__{c}__pyfin")
    sn, sp = parse(st); mb = matched(trk)
    tx = sum(1 for l in open(g) if "\ttranscript\t" in l)
    if not mb:
        print(f"{c:16} {tx:>7} {f1(sp or 0,sn or 0):>7.1f} {'-':>7} {'-':>7} {'-':>8} {'-':>7}"); continue
    me = mb & expr3
    orph = 100*sum(1 for e in mb if est.get(e,0) < 3)/len(mb)
    hpr = sp*(len(me)/len(mb)); hrec = 100*len(me)/len(expr3)
    print(f"{c:16} {tx:>7} {f1(sp,sn):>7.1f} {orph:>7.1f} {hpr:>7.1f} {hrec:>8.1f} {f1(hpr,hrec):>7.1f}")
