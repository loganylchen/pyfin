#!/usr/bin/env python3
"""Attribute prod misses/FP by diffing prod vs gates-off candidate sets via intron chain.

Usage: drop_attrib.py <prod.gtf> <nofilter.gtf> <truth.gtf> <nanocount.tsv> \
                       <prod_tmap> <nofilter_tmap>

Reports:
  - TP recovered in gates-off but DROPPED by a gate in prod (recoverable recall) + region
  - FP present in BOTH prod and gates-off (survive every gate) -> the "can't delete" set
  - FP in gates-off only (a gate already removes them)
"""
import sys, re
from collections import Counter, defaultdict

PROD, NOFILT, TRUTH, NANO, PTMAP, NTMAP = sys.argv[1:7]

def chains(path):
    """transcript_id -> (chrom, strand, intron-chain tuple, attrs)"""
    tx = {}
    for line in open(path):
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] not in ("exon", "transcript"): continue
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if not m: continue
        t = m.group(1)
        d = tx.setdefault(t, {"chrom": p[0], "strand": p[6], "exons": [], "attrs": {}})
        if p[2] == "transcript":
            d["attrs"] = dict(re.findall(r'(\w+) "([^"]*)"', p[8]))
        else:
            d["exons"].append((int(p[3]) - 1, int(p[4])))
    out = {}
    for t, d in tx.items():
        ex = sorted(d["exons"])
        intr = tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))
        out[t] = {"chrom": d["chrom"], "strand": d["strand"], "chain": intr,
                  "key": (d["chrom"], d["strand"], intr), "attrs": d["attrs"],
                  "region": "chrIS" if d["chrom"] == "chrIS" else "SIRV"}
    return out

def load_tmap(path):
    q = {}
    with open(path) as fh:
        hdr = next(fh).rstrip("\n").split("\t"); ci = {h: i for i, h in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            q[f[ci["qry_id"]]] = (f[ci["class_code"]], f[ci["ref_id"]])
    return q

prod = chains(PROD); nofilt = chains(NOFILT)
pt = load_tmap(PTMAP); nt = load_tmap(NTMAP)

# Map chain-key -> class for each run
prod_key_cc = {prod[t]["key"]: pt.get(t, ("?", ""))[0] for t in prod}
nofilt_key_cc = {nofilt[t]["key"]: nt.get(t, ("?", ""))[0] for t in nofilt}

prod_keys = set(prod_key_cc); nofilt_keys = set(nofilt_key_cc)

# Refs matched '=' in each
prod_match_ref = {pt[t][1] for t in prod if pt.get(t, ("",""))[0] == "="}
nofilt_match_ref = {nt[t][1] for t in nofilt if nt.get(t, ("",""))[0] == "="}

est = {}
fh = open(NANO); next(fh, None)
for ln in fh:
    pp = ln.rstrip("\n").split("\t")
    if len(pp) >= 3:
        try: est[pp[0]] = float(pp[2])
        except ValueError: pass
expr3 = {t for t, v in est.items() if v >= 3}

print(f"prod transcripts={len(prod)}  gates-off transcripts={len(nofilt)}")
print(f"prod '=' refs={len(prod_match_ref)}  gates-off '=' refs={len(nofilt_match_ref)}")

# Recoverable TP: expressed truth matched in gates-off but NOT in prod
recoverable = sorted((nofilt_match_ref & expr3) - prod_match_ref, key=lambda t: -est.get(t,0))
print(f"\n## Recoverable expressed-TP (matched gates-off, dropped by a gate in prod): {len(recoverable)}")
truth = chains(TRUTH)
for t in recoverable:
    reg = "chrIS" if truth.get(t,{}).get("chrom")=="chrIS" else "SIRV"
    print(f"   {t:<14} est={est.get(t,0):>8.1f}  {reg}")

# Truly missed (not even recovered with gates off) -> discovery/assembly/EM miss
truly_missed = sorted(expr3 - nofilt_match_ref, key=lambda t: -est.get(t,0))
print(f"\n## Truly missed expressed-TP (absent even gates-off; discovery/EM miss): {len(truly_missed)}")
for t in truly_missed:
    reg = "chrIS" if truth.get(t,{}).get("chrom")=="chrIS" else "SIRV"
    print(f"   {t:<14} est={est.get(t,0):>8.1f}  {reg}")

# FP analysis
GOOD = {"=", "c", "k", "m", "n"}
prod_fp_keys = {k for k, cc in prod_key_cc.items() if cc not in GOOD}
nofilt_fp_keys = {k for k, cc in nofilt_key_cc.items() if cc not in GOOD}
survive = prod_fp_keys  # everything still in prod survived all gates by definition
print(f"\n## Surviving FP in prod: {len(prod_fp_keys)} (gates-off FP: {len(nofilt_fp_keys)})")
# attach attrs from prod
key2tid = {}
for t in prod: key2tid.setdefault(prod[t]["key"], t)
rows = []
for k in survive:
    t = key2tid[k]; a = prod[t]["attrs"]
    rows.append((t, prod_key_cc[k], a.get("transcript_source","?"),
                 float(a.get("abundance",0)), int(float(a.get("num_reads",0))),
                 float(a.get("confidence",0)), prod[t]["region"]))
print(f"   by region: {dict(Counter(r[6] for r in rows))}  by class: {dict(Counter(r[1] for r in rows))}")
print(f"   {'id':<14}{'cc':>3}{'src':>6}{'abund':>9}{'nreads':>7}{'conf':>7}  region")
for r in sorted(rows, key=lambda x: -x[4]):
    print(f"   {r[0]:<14}{r[1]:>3}{r[2]:>6}{r[3]:>9.1f}{r[4]:>7}{r[5]:>7.3f}  {r[6]}")
