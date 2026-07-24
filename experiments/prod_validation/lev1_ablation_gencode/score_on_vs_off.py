#!/usr/bin/env python3
"""Gencode Lever-1 ablation: ON vs OFF, transcript-level Sn/Pr + class 'c'% delta.

No nanocount for gencode -> raw gffcompare transcript-level vs full GENCODE
(diluted Sn, but the OFF-vs-ON DELTA is meaningful since the denominator is
fixed). Also reports the class-code 'c' (5'-truncated/contained) query fraction,
which is exactly what Lever 1 targets. gffcompare via biocontainer.
"""
import os, re, glob, subprocess
REPO = "/SSD/logan/dev/pyfin"
GC = f"{REPO}/experiments/prod_validation/gencode"
OUT = f"{REPO}/experiments/prod_validation/lev1_ablation_gencode"
TRUTH = os.path.realpath(f"{GC}/_ref/full/annotation.gtf")
IMG = "quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"

def nout(p):
    n = 0
    if not os.path.exists(p): return 0
    for ln in open(p):
        if not ln.startswith("#"):
            c = ln.split("\t")
            if len(c) > 2 and c[2] == "transcript": n += 1
    return n

def gff(query, sd):
    os.makedirs(sd, exist_ok=True)
    tm = glob.glob(f"{sd}/gc.*.tmap")
    if not tm:
        if not (os.path.exists(query) and os.path.getsize(query) > 0): return (None, None, None, None)
        import shutil; shutil.copyfile(query, f"{sd}/in.gtf")
        subprocess.run(["docker","run","--rm","-u",f"{os.getuid()}:{os.getgid()}",
            "-v",f"{REPO}:{REPO}","-v","/autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD:ro",
            "-w",sd,IMG,"gffcompare","-r",TRUTH,"-o","gc","in.gtf"],
            capture_output=True)
        tm = glob.glob(f"{sd}/gc.*.tmap")
    if not tm: return (None, None, None, None)
    # parse stats + class codes
    sn = pr = None
    st = f"{sd}/gc.stats"
    if os.path.exists(st):
        for ln in open(st):
            if ln.strip().startswith("Transcript level"):
                nums = re.findall(r"[0-9]+\.[0-9]+", ln)
                if len(nums) >= 2: sn, pr = float(nums[0]), float(nums[1])
    cc = {}
    fh = open(tm[0]); h = next(fh).rstrip("\n").split("\t"); ci = {x:i for i,x in enumerate(h)}
    tot = 0
    for ln in fh:
        f = ln.rstrip("\n").split("\t"); cc[f[ci["class_code"]]] = cc.get(f[ci["class_code"]],0)+1; tot += 1
    cpct = 100*cc.get("c",0)/tot if tot else 0.0
    eqpct = 100*cc.get("=",0)/tot if tot else 0.0
    return (sn, pr, cpct, eqpct)

samples = sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(f"{GC}/SGNex_*/stage"))
rows = []
for S in samples:
    for mode in ("full","nogtf"):
        d = f"{OUT}/{S}/{mode}"
        for arm in ("off","on"):
            q = f"{d}/{arm}.gtf"
            if not os.path.exists(q): continue
            r = gff(q, f"{d}/score_{arm}")
            if r[0] is None and r[2] is None and r[3] is None: continue
            rows.append((S,mode,arm,nout(q),r[0],r[1],r[2],r[3]))

md = f"{OUT}/SUMMARY.md"
with open(md,"w") as fh:
    fh.write("# Lever-1 ablation on HUMAN GENCODE (raw gffcompare vs full annotation)\n\n")
    fh.write("OFF = e268c9b SIF output (== live OFF, byte-identical). ON = live --containment-collapse.\n")
    fh.write("Sn/Pr = transcript-level (diluted by full-annotation denominator; DELTA is the signal). c%/=% = query class-code share.\n\n")
    fh.write("| sample | mode | arm | nout | Tx_Sn | Tx_Pr | c% | =% |\n|"+"---|"*8+"\n")
    for S,mode,arm,n,sn,pr,cp,eq in rows:
        ss = S.replace("SGNex_","").replace("_directRNA","")
        fh.write(f"| {ss} | {mode} | {arm} | {n} | {sn if sn is not None else 'NA'} | {pr if pr is not None else 'NA'} | {cp:.1f} | {eq:.1f} |\n")
    fh.write("\n## Delta (on - off)\n\n| sample | mode | dNout | dTx_Pr | dc% | d=% |\n|"+"---|"*6+"\n")
    by = {(r[0],r[1],r[2]):r for r in rows}
    for S in samples:
        for mode in ("full","nogtf"):
            o=by.get((S,mode,"off")); n=by.get((S,mode,"on"))
            if not o or not n: continue
            dpr = (n[5]-o[5]) if (o[5] is not None and n[5] is not None) else None
            ss = S.replace("SGNex_","").replace("_directRNA","")
            fh.write(f"| {ss} | {mode} | {n[3]-o[3]} | {('%+.1f'%dpr) if dpr is not None else 'NA'} | {n[6]-o[6]:+.1f} | {n[7]-o[7]:+.1f} |\n")
print(open(md).read())
