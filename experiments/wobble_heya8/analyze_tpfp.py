#!/usr/bin/env python3
"""Characterize missed-TP and surviving-FP on heya8 for a pyfin run.

Usage: analyze_tpfp.py <run_dir> <truth.gtf> <nanocount.tsv> <gffcompare_prefix>
  <run_dir>            holds assembly/pyfin_work/{assembly.gtf,scores.tsv} OR pyfin.gtf
  <gffcompare_prefix>  the -o prefix used; expects <prefix>.<gtf>.tmap and .refmap

Emits:
  - class-code distribution (overall + chrIS vs SIRV)
  - missed expressed-TP (nanocount est_count>=3 truth not matched '=') with coverage
  - surviving FP (query class != '=','c','k') with pyfin attributes
"""
import sys, re, glob, os
from collections import defaultdict, Counter

RUN, TRUTH, NANO, GCPREF = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

def find_gtf(run):
    for c in (os.path.join(run, "pyfin.gtf"),
              os.path.join(run, "assembly", "pyfin_work", "assembly.gtf"),
              os.path.join(run, "pyfin_work", "assembly.gtf"),
              os.path.join(run, "assembly.gtf")):
        if os.path.exists(c):
            return c
    g = glob.glob(os.path.join(run, "**", "assembly.gtf"), recursive=True)
    return g[0] if g else None

GTF = find_gtf(RUN)
print(f"# pyfin GTF: {GTF}")

def parse_gtf_tx(path):
    """transcript_id -> dict(chrom, attrs, region)"""
    out = {}
    for line in open(path):
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "transcript": continue
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if not m: continue
        tid = m.group(1)
        attrs = dict(re.findall(r'(\w+) "([^"]*)"', p[8]))
        region = "chrIS" if p[0] == "chrIS" else "SIRV"
        out[tid] = {"chrom": p[0], "region": region, "attrs": attrs}
    return out

pyfin = parse_gtf_tx(GTF)

# nanocount expressed truth
est = {}
fh = open(NANO); next(fh, None)
for ln in fh:
    pp = ln.rstrip("\n").split("\t")
    if len(pp) >= 3:
        try: est[pp[0]] = float(pp[2])
        except ValueError: pass
expr3 = {t for t, v in est.items() if v >= 3}
expr1 = {t for t, v in est.items() if v >= 1}

def region_of(tid):
    return "chrIS" if tid.startswith("R") or tid.startswith("chrIS") else "SIRV"

# gffcompare tmap: query -> class, ref ; refmap: ref -> class, querylist
tmap = glob.glob(f"{GCPREF}.*.tmap")
refmap = glob.glob(f"{GCPREF}.*.refmap")
print(f"# tmap={tmap} refmap={refmap}")

q_class = {}   # query tid -> (class_code, ref_id)
matched_ref = set()  # ref tids with '=' match
if tmap:
    with open(tmap[0]) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        ci = {h: i for i, h in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            cc = f[ci["class_code"]]
            qid = f[ci["qry_id"]]
            rid = f[ci.get("ref_id", 2)]
            q_class[qid] = (cc, rid)
            if cc == "=":
                matched_ref.add(rid)

# class-code distribution by region
print("\n## Class-code distribution (pyfin query transcripts)")
by_region = defaultdict(Counter)
for tid, info in pyfin.items():
    cc = q_class.get(tid, ("?", ""))[0]
    by_region[info["region"]][cc] += 1
allc = Counter()
for r in by_region: allc.update(by_region[r])
print(f"  overall: {dict(allc)}")
for r in sorted(by_region):
    print(f"  {r}: {dict(by_region[r])}  (total {sum(by_region[r].values())})")

# truth regions
truth_tx = parse_gtf_tx(TRUTH)
truth_region = {t: ("chrIS" if truth_tx[t]["chrom"] == "chrIS" else "SIRV") for t in truth_tx}

# Missed expressed TP
print(f"\n## Expressed truth: >=1read={len(expr1)}  >=3read={len(expr3)}")
missed3 = sorted(expr3 - matched_ref, key=lambda t: -est.get(t, 0))
print(f"## MISSED expressed-TP (est_count>=3, not '=' matched): {len(missed3)}")
print(f"   {'truth_id':<14}{'est_cnt':>9}  region")
for t in missed3:
    print(f"   {t:<14}{est.get(t,0):>9.1f}  {truth_region.get(t,'?')}")

# Surviving FP (query class not in good set)
GOOD = {"=", "c", "k", "m", "n"}
fp = [(tid, q_class.get(tid, ("?", ""))) for tid in pyfin
      if q_class.get(tid, ("?", ""))[0] not in GOOD]
print(f"\n## Surviving FP (class not in {GOOD}): {len(fp)}")
fp_region = Counter(pyfin[t]["region"] for t, _ in fp)
fp_class = Counter(c[0] for _, c in fp)
print(f"   by region: {dict(fp_region)}   by class: {dict(fp_class)}")
print(f"   {'pyfin_id':<14}{'cc':>3}{'src':>6}{'abund':>9}{'nreads':>7}{'conf':>7}  region")
def keyf(x):
    a = pyfin[x[0]]["attrs"]
    return -float(a.get("num_reads", 0))
for tid, (cc, rid) in sorted(fp, key=keyf):
    a = pyfin[tid]["attrs"]
    print(f"   {tid:<14}{cc:>3}{a.get('transcript_source','?'):>6}"
          f"{float(a.get('abundance',0)):>9.1f}{int(float(a.get('num_reads',0))):>7}"
          f"{float(a.get('confidence',0)):>7.3f}  {pyfin[tid]['region']}")
