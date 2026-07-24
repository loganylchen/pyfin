#!/usr/bin/env python3
"""Thorough OFF-vs-ON evaluation of the guided junction-support gate on the two
decisive conditions (full = don't-break-baseline, c_jitter = the fix target).
For each: tx count, standard Sn/Pr/F1, orphan%, honest precision (Pr*(1-orphan)),
corrRec(est>=3), and honest F1 at est>=3 and est>=1 — OFF (prodfull) vs ON
(prodfull_g2), with deltas. This is the precision-lever verdict.
"""
import re, os, subprocess
B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
IMG = "quay.io/biocontainers/gffcompare:0.12.10--h9948957_0"
S = "SGNex_H9_directRNA_replicate2_run2"
NC = f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
ENST = re.compile(r"(ENST\d+)")
CONDS = ["full", "c_jitter20_10bp"]

est = {}
for ln in open(NC):
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 3:
        m = ENST.search(p[0])
        if m:
            try: est[m.group(1)] = float(p[2])
            except ValueError: pass
def exprset(thr): return {k for k,v in est.items() if v >= thr}

def gff(gtf, pre):
    if not os.path.exists(pre+".stats"):
        subprocess.run(["docker","run","--rm","-u",f"{os.getuid()}:{os.getgid()}",
            "-v","/autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD","-v","/SSD:/SSD",
            IMG,"gffcompare","-r",TRUTH,"-o",pre,gtf], capture_output=True)
    return pre+".stats", pre+".tracking"
def parse(st):
    t = open(st).read()
    m = re.search(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)", t)
    return (float(m.group(1)), float(m.group(2))) if m else (None, None)
def matched(trk):
    s = set()
    for ln in open(trk):
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 4 and c[3] == "=":
            m = ENST.search(c[2])
            if m: s.add(m.group(1))
    return s
def f1(p, r): return 0.0 if p+r == 0 else 2*p*r/(p+r)

def metrics(gtf, pre):
    if not os.path.exists(gtf): return None
    st, trk = gff(gtf, pre)
    sn, sp = parse(st)
    mb = matched(trk)
    tx = sum(1 for l in open(gtf) if "\ttranscript\t" in l)
    out = {"tx": tx, "sn": sn, "pr": sp, "f1": f1(sp, sn)}
    for thr in (1, 3):
        ex = exprset(thr); me = mb & ex
        orph = 100*sum(1 for e in mb if est.get(e,0) < thr)/len(mb) if mb else 0
        hpr = sp * (len(me)/len(mb)) if mb else 0
        hrec = 100*len(me)/len(ex)
        out[f"orphan{thr}"] = orph
        out[f"hpr{thr}"] = hpr
        out[f"hrec{thr}"] = hrec
        out[f"hf1{thr}"] = f1(hpr, hrec)
    return out

print(f"# GUIDED GATE eval — OFF(prodfull) vs ON(prodfull_g2, guided_min_reads=2)  {S}")
for c in CONDS:
    off = metrics(f"{HERE}/prodfull/{c}/pyfin.gtf", f"{HERE}/gc_full/{S}__{c}__pyfin")
    on = metrics(f"{HERE}/prodfull_g2/{c}/pyfin.gtf", f"{HERE}/gc_g2/{S}__{c}__pyfin")
    print(f"\n## {c}")
    if not off or not on:
        print(f"   OFF={'ok' if off else 'MISSING'} ON={'ok' if on else 'MISSING'}"); continue
    def row(name, k, fmt="{:.1f}"):
        o, n = off[k], on[k]
        d = n - o
        print(f"   {name:22} OFF {fmt.format(o):>8}   ON {fmt.format(n):>8}   Δ {d:+.1f}")
    row("tx count", "tx", "{:.0f}")
    row("Sn (recall)", "sn"); row("Pr (precision)", "pr"); row("F1 (standard)", "f1")
    row("orphan% (est<3)", "orphan3"); row("honest Pr (est>=3)", "hpr3")
    row("corrRec (est>=3)", "hrec3"); row("HONEST F1 (est>=3)", "hf13")
    row("HONEST F1 (est>=1)", "hf11")
