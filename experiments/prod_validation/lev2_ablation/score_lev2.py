#!/usr/bin/env python3
"""Lever-3 (mono gate) ablation scorer: off vs on per cell.

Key metric = mono% (single-exon share of output) — should DROP on. Plus:
- sirv/heya8: honest F1@3/Sn@3 (nanocount expressed-truth), to check recall.
- gencode: raw transcript-level Sn/Pr vs full GENCODE.
gffcompare via biocontainer (/autofs mounted for gencode truth symlink).
"""
import os, re, glob, subprocess
REPO = "/SSD/logan/dev/pyfin"
OUT = f"{REPO}/experiments/prod_validation/lev2_ablation"
IMG = "quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"
DSROOT = {
    "sirv":   f"{REPO}/experiments/prod_validation/sirv4",
    "heya8":  f"{REPO}/experiments/wobble_heya8/matrix",
    "gencode":f"{REPO}/experiments/prod_validation/gencode",
}
TRUTH = {k: os.path.realpath(f"{v}/_ref/full/annotation.gtf") for k,v in DSROOT.items()}

def tx_and_mono(p):
    """(nout, mono_count): count transcripts and how many are single-exon."""
    if not os.path.exists(p): return 0, 0
    exons = {}
    tx = set()
    for ln in open(p):
        if ln.startswith("#"): continue
        c = ln.rstrip("\n").split("\t")
        if len(c) < 9: continue
        m = re.search(r'transcript_id "([^"]+)"', c[8])
        if not m: continue
        t = m.group(1)
        if c[2] == "transcript": tx.add(t)
        elif c[2] == "exon": exons[t] = exons.get(t,0)+1
    mono = sum(1 for t in tx if exons.get(t,0) == 1)
    return len(tx), mono

def nano(p):
    e={}
    if not os.path.exists(p): return e
    fh=open(p); next(fh,None)
    for ln in fh:
        c=ln.rstrip("\n").split("\t")
        if len(c)>=3:
            try: e[c[0]]=float(c[2])
            except ValueError: pass
    return e

def gff_matched(query, sd, truth):
    os.makedirs(sd, exist_ok=True)
    tm=glob.glob(f"{sd}/gc.*.tmap")
    if not tm:
        if not (os.path.exists(query) and os.path.getsize(query)>0): return None,None,None
        import shutil; shutil.copyfile(query,f"{sd}/in.gtf")
        subprocess.run(["docker","run","--rm","-u",f"{os.getuid()}:{os.getgid()}",
            "-v",f"{REPO}:{REPO}","-v","/autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD:ro",
            "-w",sd,IMG,"gffcompare","-r",truth,"-o","gc","in.gtf"],capture_output=True)
        tm=glob.glob(f"{sd}/gc.*.tmap")
    if not tm: return None,None,None
    matched=set();
    fh=open(tm[0]); h=next(fh).rstrip("\n").split("\t"); ci={x:i for i,x in enumerate(h)}
    for ln in fh:
        f=ln.rstrip("\n").split("\t")
        if f[ci["class_code"]]=="=": matched.add(f[ci["ref_id"]])
    sn=pr=None
    st=f"{sd}/gc.stats"
    if os.path.exists(st):
        for ln in open(st):
            if ln.strip().startswith("Transcript level"):
                n=re.findall(r"[0-9]+\.[0-9]+",ln)
                if len(n)>=2: sn,pr=float(n[0]),float(n[1])
    return matched, sn, pr

def truth_tx(p):
    s=set()
    for ln in open(p):
        if ln.startswith("#"): continue
        c=ln.rstrip("\n").split("\t")
        if len(c)>=9 and c[2]=="transcript":
            m=re.search(r'transcript_id "([^"]+)"',c[8])
            if m: s.add(m.group(1))
    return s
def f1(a,b): return 2*a*b/(a+b) if (a+b) else 0.0

rows=[]
for ds in ("sirv","heya8","gencode"):
    truth=TRUTH[ds]; ttx=truth_tx(truth)
    for S in sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(f"{OUT}/{ds}/SGNex_*/")):
        for mode in ("denovo","guided","nogtf"):
            d=f"{OUT}/{ds}/{S}/{mode}"
            if not os.path.isdir(d): continue
            e3=None
            if ds in ("sirv","heya8"):
                est=nano(f"{DSROOT[ds]}/{S}/stage/nanocount.tsv")
                e3={t for t,v in est.items() if v>=3 and t in ttx}
            for arm in ("off","on"):
                q=f"{d}/{arm}.gtf"
                if not os.path.exists(q): continue
                nout,mono=tx_and_mono(q)
                mref,sn,pr=gff_matched(q,f"{d}/score_{arm}",truth)
                monop=100*mono/nout if nout else 0.0
                if e3 is not None and mref is not None:
                    m3=len(mref&e3); sn3=100*m3/max(len(e3),1); pr3=100*m3/max(nout,1); f13=f1(sn3,pr3)
                    rows.append((ds,S,mode,arm,nout,monop,sn3,pr3,f13))
                else:
                    rows.append((ds,S,mode,arm,nout,monop,sn if sn else 0.0,pr if pr else 0.0,0.0))

md=f"{OUT}/SUMMARY.md"
with open(md,"w") as fh:
    fh.write("# Lever-3 (mono gate) ablation — off vs on\n\n")
    fh.write("ON = --drop-mono-exon-novel --min-mono-exon-reads 3 --min-mono-exon-length 200.\n")
    fh.write("mono% = single-exon share of output (Lever-3 target, should drop). sirv/heya8: Sn3/Pr3/F1@3 honest. gencode: Tx_Sn/Pr raw.\n\n")
    fh.write("| ds | sample | mode | arm | nout | mono% | Sn/Sn3 | Pr/Pr3 | F1@3 |\n|"+"---|"*9+"\n")
    for ds,S,mode,arm,nout,mp,sn,pr,f13 in rows:
        ss=S.replace("SGNex_","").replace("_directRNA","")
        fh.write(f"| {ds} | {ss} | {mode} | {arm} | {nout} | {mp:.1f} | {sn:.1f} | {pr:.1f} | {f13:.1f} |\n")
    fh.write("\n## Delta (on - off)\n\n| ds | sample | mode | dNout | dmono% | dSn | dPr | dF1@3 |\n|"+"---|"*8+"\n")
    by={(r[0],r[1],r[2],r[3]):r for r in rows}
    seen=set()
    for ds,S,mode,arm,*_ in rows:
        k=(ds,S,mode)
        if k in seen: continue
        seen.add(k)
        o=by.get((ds,S,mode,"off")); n=by.get((ds,S,mode,"on"))
        if not o or not n: continue
        ss=S.replace("SGNex_","").replace("_directRNA","")
        fh.write(f"| {ds} | {ss} | {mode} | {n[4]-o[4]} | {n[5]-o[5]:+.1f} | {n[6]-o[6]:+.1f} | {n[7]-o[7]:+.1f} | {n[8]-o[8]:+.1f} |\n")
    fh.write("\n## Mean delta by (ds, mode)\n\n| ds | mode | dNout | dmono% | dSn | dPr | dF1@3 |\n|"+"---|"*7+"\n")
    groups={}
    for ds,S,mode,arm,*rest in rows: groups.setdefault((ds,mode),set()).add(S)
    for (ds,mode),ss in sorted(groups.items()):
        ds_=[]
        for S in ss:
            o=by.get((ds,S,mode,"off")); n=by.get((ds,S,mode,"on"))
            if o and n: ds_.append((n[4]-o[4],n[5]-o[5],n[6]-o[6],n[7]-o[7],n[8]-o[8]))
        if not ds_: continue
        k=len(ds_)
        fh.write(f"| {ds} | {mode} | {sum(x[0] for x in ds_)/k:.0f} | {sum(x[1] for x in ds_)/k:+.1f} | {sum(x[2] for x in ds_)/k:+.1f} | {sum(x[3] for x in ds_)/k:+.1f} | {sum(x[4] for x in ds_)/k:+.1f} |\n")
print(open(md).read())
