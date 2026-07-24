#!/usr/bin/env python3
"""Score a prod_validation dataset: pyfin (production) + competitor tools vs full truth.

Usage: score_dataset.py <dataset_dir> <nfs_sweep_dir>
  <dataset_dir>   = .../prod_validation/sirv4  (holds <sample>/<ratio>/pyfin.gtf,
                    _ref/<ratio>/annotation.gtf, _ref/full/annotation.gtf,
                    <sample>/stage/nanocount.tsv)
  <nfs_sweep_dir> = .../benchmark/sgnex/results/sirv4_full_sweep (competitor GTFs:
                    <sample>/<ratio>/assembly/<tool>.gtf)

Honest metrics vs full truth: expressed-truth = nanocount est_count >= K.
Runs gffcompare via the biocontainer (idempotent: skips existing tmap).
Writes tables/<sample>.tsv and tables/SUMMARY.md under <dataset_dir>.
"""
import sys, os, re, glob, subprocess
from collections import defaultdict

DS, NFS = sys.argv[1], sys.argv[2]
TRUTH = os.path.join(DS, "_ref", "full", "annotation.gtf")
GFFC_IMG = "quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"
COMPETITORS = ["bambu", "espresso", "flair", "isoquant", "isotools", "lafite", "stringtie3", "talon"]

samples = sorted(os.path.basename(os.path.dirname(p))
                 for p in glob.glob(os.path.join(DS, "*", "stage")))
ratios = sorted(os.path.basename(os.path.dirname(p))
                for p in glob.glob(os.path.join(DS, "_ref", "*", "annotation.gtf")))
# canonical ratio order
order = ["full","p99","p90","p50","p10","p00","c_skip10","c_skip20","c_jitter10bp",
         "c_jitter20_10bp","c_spurious5","c_spurious20","c_merge5","c_merge20",
         "c_flip5","c_flip20","c_ir10","c_ir20"]
ratios = [r for r in order if r in ratios] + [r for r in ratios if r not in order]

def parse_tx(path):
    s=set()
    if not path or not os.path.exists(path): return s
    for ln in open(path):
        if ln.startswith("#"): continue
        p=ln.rstrip("\n").split("\t")
        if len(p)>=9 and p[2]=="transcript":
            m=re.search(r'transcript_id "([^"]+)"',p[8])
            if m: s.add(m.group(1))
    return s

def nanocount(path):
    e={}
    if not os.path.exists(path): return e
    fh=open(path); next(fh,None)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        if len(p)>=3:
            try: e[p[0]]=float(p[2])
            except ValueError: pass
    return e

MNT="/SSD/logan/dev/pyfin"
import shutil
def gffcompare(query, prefix, scoredir=None):
    """run gffcompare; tmap lands next to the (writable) query as <prefix>.<qbase>.tmap.

    gffcompare writes tmap/refmap next to the input, so a read-only NFS query
    (competitor GTFs under /autofs) must be copied into a writable local dir first.
    """
    if scoredir:
        os.makedirs(scoredir, exist_ok=True)
        tmaps=glob.glob(os.path.join(scoredir, f"{prefix}.*.tmap"))
        if tmaps: return tmaps[0]
    else:
        tmaps=glob.glob(os.path.join(os.path.dirname(query), f"{prefix}.*.tmap"))
        if tmaps: return tmaps[0]
    if not (os.path.exists(query) and os.path.getsize(query)>0): return None
    # Always run inside the writable scoredir (gffcompare writes tmap next to the
    # query; NFS competitor GTFs are read-only). Copy the query in with a tool-
    # specific name so multiple tools' tmaps coexist.
    if scoredir:
        local_q=os.path.join(scoredir, f"{prefix}.input.gtf")
        shutil.copyfile(query, local_q); query=local_q
    qd=os.path.dirname(os.path.realpath(query))
    subprocess.run(["docker","run","--rm","-u",f"{os.getuid()}:{os.getgid()}",
        "-v",f"{MNT}:{MNT}","-v","/autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD:ro",
        "-w",qd,GFFC_IMG,"gffcompare","-r",os.path.realpath(TRUTH),"-o",prefix,
        os.path.realpath(query)],capture_output=True)
    g=glob.glob(os.path.join(qd,f"{prefix}.*.tmap"))
    return g[0] if g else None

def matched(tmap):
    out=set()
    if not tmap or not os.path.exists(tmap): return out
    with open(tmap) as fh:
        h=next(fh).rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
        for ln in fh:
            f=ln.rstrip("\n").split("\t")
            if f[ci["class_code"]]=="=": out.add(f[ci["ref_id"]])
    return out

truth_tx=parse_tx(TRUTH)
def f1(s,p): return 2*s*p/(s+p) if (s+p) else 0.0
os.makedirs(os.path.join(DS,"tables"),exist_ok=True)
allrows=[]
for s in samples:
    est=nanocount(os.path.join(DS,s,"stage","nanocount.tsv"))
    e1={t for t,v in est.items() if v>=1 and t in truth_tx}
    e3={t for t,v in est.items() if v>=3 and t in truth_tx}
    with open(os.path.join(DS,"tables",f"{s}.tsv"),"w") as fh:
        fh.write("ratio\ttool\tnout\tSn1\tPr1\tF1@1\tSn3\tPr3\tF1@3\n")
        for r in ratios:
            cells=[("pyfin_prod", os.path.join(DS,s,r,"pyfin.gtf"))]
            for t in COMPETITORS:
                cells.append((t, os.path.join(NFS,s,r,"assembly",f"{t}.gtf")))
            for tool,q in cells:
                tm=gffcompare(q, f"gc_{tool}", scoredir=os.path.join(DS,s,r,"scoring"))
                if not tm:
                    fh.write(f"{r}\t{tool}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"); continue
                nout=len(parse_tx(q)); mref=matched(tm)
                m1=len(mref&e1); m3=len(mref&e3)
                sn1=100*m1/max(len(e1),1); pr1=100*m1/max(nout,1)
                sn3=100*m3/max(len(e3),1); pr3=100*m3/max(nout,1)
                fh.write(f"{r}\t{tool}\t{nout}\t{sn1:.1f}\t{pr1:.1f}\t{f1(sn1,pr1):.1f}"
                         f"\t{sn3:.1f}\t{pr3:.1f}\t{f1(sn3,pr3):.1f}\n")
                allrows.append((s,r,tool,nout,sn1,pr1,f1(sn1,pr1),sn3,pr3,f1(sn3,pr3)))
    print(f"scored {s} (expr1={len(e1)} expr3={len(e3)})")

# cross-sample mean F1@3 per (ratio, tool)
md=os.path.join(DS,"tables","SUMMARY.md")
tools=["pyfin_prod"]+COMPETITORS
with open(md,"w") as fh:
    fh.write(f"# {os.path.basename(DS)} — production pyfin vs competitors (honest F1@3, mean over samples)\n\n")
    fh.write("| ratio | " + " | ".join(tools) + " |\n")
    fh.write("|" + "---|"*(len(tools)+1) + "\n")
    for r in ratios:
        row=[r]
        for t in tools:
            vals=[x[9] for x in allrows if x[1]==r and x[2]==t]
            row.append(f"{sum(vals)/len(vals):.1f}" if vals else "NA")
        fh.write("| "+" | ".join(row)+" |\n")
    # global means
    fh.write("\n**Global mean F1@3 (all ratios × samples):**\n\n")
    fh.write("| tool | F1@3 | Sn@3 | Pr@3 |\n|---|---|---|---|\n")
    for t in tools:
        f=[x[9] for x in allrows if x[2]==t]; sn=[x[7] for x in allrows if x[2]==t]; pr=[x[8] for x in allrows if x[2]==t]
        if f: fh.write(f"| {t} | {sum(f)/len(f):.1f} | {sum(sn)/len(sn):.1f} | {sum(pr)/len(pr):.1f} |\n")
print(f"wrote {md}")
