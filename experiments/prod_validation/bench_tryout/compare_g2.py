#!/usr/bin/env python3
"""OFF (prodfull) vs ON (prodfull_g2, --guided-junction-min-reads 2) ablation.
Scores ON with the benchmark's own gffcompare + orphan% + corrected-recall,
prints per-condition Sn/Pr/F1 and the key verdict (jitter precision, orphan, recall)."""
import os, re, subprocess
B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
IMG = "quay.io/biocontainers/gffcompare:0.12.10--h9948957_0"
S = "SGNex_H9_directRNA_replicate2_run2"
RATIOS = ["full","p99","p90","p50","p10","p00","c_skip20","c_jitter20_10bp",
          "c_spurious20","c_merge20","c_flip20","c_ir20"]
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
ENST_RE = re.compile(r"(ENST\d+\.\d+)")

def gff(gtf, pre):
    if not os.path.exists(pre+".stats"):
        subprocess.run(["docker","run","--rm","-u",f"{os.getuid()}:{os.getgid()}",
            "-v","/autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD","-v","/SSD:/SSD",
            IMG,"gffcompare","-r",TRUTH,"-o",pre,gtf], capture_output=True)
    return pre+".stats", pre+".tracking"

def parse(stats):
    if not os.path.exists(stats): return (None,None)
    t=open(stats).read(); m=re.search(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)",t)
    return (float(m.group(1)),float(m.group(2))) if m else (None,None)

def f1(sn,sp): return 0.0 if not sn or not sp or sn+sp==0 else round(2*sn*sp/(sn+sp),2)

# expressed truth (nanocount est>=3)
counts={}
with open(f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv") as fh:
    next(fh,None)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        if len(p)<3: continue
        try: counts[p[0].split("|")[0].split(".")[0]]=float(p[2])
        except ValueError: pass
expr3={k for k,v in counts.items() if v>=3}

def orphan_recall(trk):
    if not trk or not os.path.exists(trk): return (None,None)
    mb=set()
    for ln in open(trk):
        c=ln.rstrip("\n").split("\t")
        if len(c)<4 or c[3]!="=": continue
        m=ENST_RE.search(c[2])
        if m: mb.add(m.group(1).split(".")[0])
    if not mb: return (None,None)
    orphan=round(100.0*sum(1 for e in mb if counts.get(e,0.0)<1)/len(mb),1)
    rec=round(100.0*len(mb&expr3)/max(len(expr3),1),1)
    return orphan,rec

print(f"# ABLATION  guided-junction-min-reads: OFF(prodfull) vs ON=2(prodfull_g2)")
print(f"# {S} | expressed truth(est>=3)={len(expr3)}\n")
print(f"{'condition':16} | {'OFF Sn/Pr/F1':>20} | {'ON  Sn/Pr/F1':>20} | {'dF1':>6} {'dPr':>6}")
off_full=on_full=None
for r in RATIOS:
    offs=f"{HERE}/gc_full/{S}__{r}__pyfin"
    ons=f"{HERE}/gc_g2/{S}__{r}__pyfin"
    o_sn,o_sp=parse(offs+".stats")
    on_gtf=f"{HERE}/prodfull_g2/{r}/pyfin.gtf"
    n_sn=n_sp=None
    if os.path.exists(on_gtf):
        st,_=gff(on_gtf,ons); n_sn,n_sp=parse(st)
    of1=f1(o_sn,o_sp); nf1=f1(n_sn,n_sp)
    d1 = round(nf1-of1,2) if (o_sn and n_sn) else None
    dpr= round((n_sp or 0)-(o_sp or 0),1) if (o_sp and n_sp) else None
    mark=" <== JITTER" if r=="c_jitter20_10bp" else ""
    print(f"{r:16} | {f'{o_sn}/{o_sp}/{of1}':>20} | {f'{n_sn}/{n_sp}/{nf1}':>20} | {str(d1):>6} {str(dpr):>6}{mark}")

# orphan + corrRec at full
o_orph,o_rec=orphan_recall(f"{HERE}/gc_full/{S}__full__pyfin.tracking")
n_orph,n_rec=orphan_recall(f"{HERE}/gc_g2/{S}__full__pyfin.tracking")
print(f"\n@full  orphan%: OFF {o_orph} -> ON {n_orph}   |   corrRec>=3: OFF {o_rec} -> ON {n_rec}")
