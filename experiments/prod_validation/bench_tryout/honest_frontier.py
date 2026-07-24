#!/usr/bin/env python3
"""p00 HONEST two-axis frontier: pyfin vs all competitors. For each tool compute
standard gffcompare Sn/Pr AND honest precision (Pr*(1-orphan)), honest recall
(corrRec = matched expressed-truth / expressed denom), honest-F1. Answers: at p00,
who is Pareto-best on BOTH recall and precision, and where does pyfin sit?
"""
import os, re, subprocess
B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
IMG = "quay.io/biocontainers/gffcompare:0.12.10--h9948957_0"
S = "SGNex_H9_directRNA_replicate2_run2"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
NC = f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
ASM = f"{B}/results/gencode_full_sweep/{S}/p00/assembly"
ENST = re.compile(r"(ENST\d+)")
sc = f"{HERE}/gc_frontier"
os.makedirs(sc, exist_ok=True)

est = {}
for ln in open(NC):
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 3:
        m = ENST.search(p[0])
        if m:
            try: est[m.group(1)] = float(p[2])
            except ValueError: pass
expr3 = {k for k, v in est.items() if v >= 3}
DEN = len(expr3)

GTFS = {
    "pyfin": f"{HERE}/prodfull/p00/pyfin.gtf",
    "isoquant": f"{ASM}/isoquant.gtf", "bambu": f"{ASM}/bambu.gtf",
    "talon": f"{ASM}/talon.gtf", "stringtie3": f"{ASM}/stringtie3.gtf",
    "lafite": f"{ASM}/lafite.gtf", "isotools": f"{ASM}/isotools.gtf",
}


def gff(g, pre):
    if not os.path.exists(pre + ".stats"):
        subprocess.run(["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}",
            "-v", "/autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD", "-v", "/SSD:/SSD",
            IMG, "gffcompare", "-r", TRUTH, "-o", pre, g], capture_output=True)
    return pre + ".stats", pre + ".tracking"


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


def f1(a, b): return 0.0 if a + b == 0 else 2 * a * b / (a + b)

print(f"p00 HONEST frontier  {S}  (expressed denom={DEN})")
print(f"{'tool':11}{'tx':>7}{'stdSn':>7}{'stdPr':>7}{'orphan':>7}{'honPr':>7}{'corrRec':>8}{'honF1':>7}")
rows = []
for tool, g in GTFS.items():
    if not os.path.exists(g):
        print(f"{tool:11}{'MISSING':>7}"); continue
    st, trk = gff(g, f"{sc}/{S}__p00__{tool}")
    sn, sp = parse(st); mb = matched(trk)
    tx = sum(1 for l in open(g) if "\ttranscript\t" in l)
    if not mb:
        print(f"{tool:11}{tx:>7}{sn or 0:>7.1f}{sp or 0:>7.1f}{'-':>7}"); continue
    me = mb & expr3
    orph = 100 * sum(1 for e in mb if est.get(e, 0) < 3) / len(mb)
    hpr = sp * (len(me) / len(mb)); hrec = 100 * len(me) / DEN
    rows.append((tool, tx, sn, sp, orph, hpr, hrec, f1(hpr, hrec)))
    print(f"{tool:11}{tx:>7}{sn:>7.1f}{sp:>7.1f}{orph:>7.1f}{hpr:>7.1f}{hrec:>8.1f}{f1(hpr,hrec):>7.1f}")

print("\n-- honest Pareto frontier (who dominates on BOTH honPr and corrRec) --")
for t, tx, sn, sp, o, hp, hr, hf in rows:
    dominators = [u[0] for u in rows if u[0] != t and u[5] >= hp and u[6] >= hr
                  and (u[5] > hp or u[6] > hr)]
    tag = "PARETO-BEST" if not dominators else "dominated by " + ",".join(dominators)
    print(f"  {t:11} honPr={hp:5.1f} corrRec={hr:5.1f}  -> {tag}")
