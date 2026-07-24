#!/usr/bin/env python3
"""Score heya8 Lever-1 ablation: honest F1@1/@3, off vs on, per (sample,mode)."""
import os, re, glob, subprocess
REPO = "/SSD/logan/dev/pyfin"
DS = f"{REPO}/experiments/wobble_heya8/matrix"
ABL = f"{REPO}/experiments/prod_validation/lev1_ablation_heya8"
TRUTH = f"{DS}/_ref/full/annotation.gtf"
IMG = "quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"

def parse_tx(p):
    s=set()
    if not os.path.exists(p): return s
    for ln in open(p):
        if ln.startswith("#"): continue
        c=ln.rstrip("\n").split("\t")
        if len(c)>=9 and c[2]=="transcript":
            m=re.search(r'transcript_id "([^"]+)"',c[8])
            if m: s.add(m.group(1))
    return s
def nanocount(p):
    e={}
    if not os.path.exists(p): return e
    fh=open(p); next(fh,None)
    for ln in fh:
        c=ln.rstrip("\n").split("\t")
        if len(c)>=3:
            try: e[c[0]]=float(c[2])
            except ValueError: pass
    return e
def gffcompare(query, sd):
    os.makedirs(sd, exist_ok=True)
    tm=glob.glob(f"{sd}/gc.*.tmap")
    if tm: return tm[0]
    if not (os.path.exists(query) and os.path.getsize(query)>0): return None
    import shutil; shutil.copyfile(query, f"{sd}/in.gtf")
    subprocess.run(["docker","run","--rm","-u",f"{os.getuid()}:{os.getgid()}",
        "-v",f"{REPO}:{REPO}","-w",sd,IMG,"gffcompare","-r",TRUTH,"-o","gc","in.gtf"],
        capture_output=True)
    tm=glob.glob(f"{sd}/gc.*.tmap"); return tm[0] if tm else None
def matched(tm):
    out=set()
    if not tm or not os.path.exists(tm): return out
    fh=open(tm); h=next(fh).rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for ln in fh:
        f=ln.rstrip("\n").split("\t")
        if f[ci["class_code"]]=="=": out.add(f[ci["ref_id"]])
    return out
def f1(a,b): return 2*a*b/(a+b) if (a+b) else 0.0
truth=parse_tx(TRUTH)
samples=sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(f"{DS}/SGNex_HEYA8*/stage"))
rows=[]
for S in samples:
    est=nanocount(f"{DS}/{S}/stage/nanocount.tsv")
    e1={t for t,v in est.items() if v>=1 and t in truth}
    e3={t for t,v in est.items() if v>=3 and t in truth}
    for mode in ("denovo","guided"):
        for arm in ("off","on"):
            q=f"{ABL}/{S}/{mode}/{arm}/pyfin.gtf"
            if not os.path.exists(q): continue
            tm=gffcompare(q, f"{ABL}/{S}/{mode}/{arm}/score")
            nout=len(parse_tx(q)); mref=matched(tm)
            m1,m3=len(mref&e1),len(mref&e3)
            sn3=100*m3/max(len(e3),1); pr3=100*m3/max(nout,1)
            sn1=100*m1/max(len(e1),1); pr1=100*m1/max(nout,1)
            rows.append((S,mode,arm,nout,len(e3),sn1,pr1,f1(sn1,pr1),sn3,pr3,f1(sn3,pr3)))
md=f"{ABL}/SUMMARY.md"
with open(md,"w") as fh:
    fh.write("# Lever-1 containment ablation (heya8 dense-locus, honest metrics, CPU)\n\n")
    fh.write("| sample | mode | arm | nout | expr3 | Sn@3 | Pr@3 | F1@3 | F1@1 |\n|"+"---|"*9+"\n")
    for r in rows:
        S,mode,arm,nout,e3n,sn1,pr1,f11,sn3,pr3,f13=r
        fh.write(f"| {S.replace('SGNex_','').replace('_directRNA','')} | {mode} | {arm} | {nout} | {e3n} | {sn3:.1f} | {pr3:.1f} | {f13:.1f} | {f11:.1f} |\n")
    fh.write("\n## Delta (on - off)\n\n| sample | mode | dNout | dSn@3 | dPr@3 | dF1@3 |\n|"+"---|"*6+"\n")
    by={(r[0],r[1],r[2]):r for r in rows}
    for S in samples:
        for mode in ("denovo","guided"):
            o=by.get((S,mode,"off")); n=by.get((S,mode,"on"))
            if not o or not n: continue
            fh.write(f"| {S.replace('SGNex_','').replace('_directRNA','')} | {mode} | {n[3]-o[3]} | {n[8]-o[8]:+.1f} | {n[9]-o[9]:+.1f} | {n[10]-o[10]:+.1f} |\n")
    fh.write("\n## Mean over samples\n\n| mode | arm | mean nout | mean Sn@3 | mean Pr@3 | mean F1@3 |\n|"+"---|"*6+"\n")
    for mode in ("denovo","guided"):
        for arm in ("off","on"):
            g=[r for r in rows if r[1]==mode and r[2]==arm]
            if not g: continue
            k=len(g)
            fh.write(f"| {mode} | {arm} | {sum(r[3] for r in g)/k:.0f} | {sum(r[8] for r in g)/k:.1f} | {sum(r[9] for r in g)/k:.1f} | {sum(r[10] for r in g)/k:.1f} |\n")
print(open(md).read())
