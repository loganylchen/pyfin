#!/usr/bin/env python3
"""Split pyfin's MISSED expressed transcripts into
  (a) generated-then-filtered  = present in scores.unfiltered.tsv (reached EM,
      dropped by a POST-EM filter), and
  (b) pre-EM / never-generated = absent from unfiltered too,
and attribute (a) to the likely filter via abundance / num_reads.

Reuses gc_full/<S>__full__pyfin.tracking for the FINAL matched set (the unfiltered
run's assembly.gtf is production-identical). Reads the new run's unfiltered.tsv.
"""
import re, collections
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S="SGNex_H9_directRNA_replicate2_run2"
NC=f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
TRK=f"gc_full/{S}__full__pyfin.tracking"
UNF="prodfull_unfilt/work/scores.unfiltered.tsv"
ENST=re.compile(r"(ENST\d+)")

est={}
with open(NC) as fh:
    next(fh,None)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        if len(p)<3: continue
        m=ENST.search(p[0])
        if m:
            try: est[m.group(1)]=float(p[2])
            except ValueError: pass
expr={k for k,v in est.items() if v>=3}

matched=set()
for ln in open(TRK):
    c=ln.rstrip("\n").split("\t")
    if len(c)>=4 and c[3]=="=":
        m=ENST.search(c[2])
        if m: matched.add(m.group(1))

# unfiltered (post-EM, pre-filter): base ENST -> (abundance, num_reads) for GTF source
unf={}
with open(UNF) as fh:
    hdr=fh.readline().rstrip("\n").split("\t"); ix={c:i for i,c in enumerate(hdr)}
    for ln in fh:
        c=ln.rstrip("\n").split("\t")
        if c[ix["source"]]!="gtf": continue
        m=ENST.search(c[ix["candidate_id"]])
        if not m: continue
        try: unf[m.group(1)]=(float(c[ix["abundance"]]), float(c[ix["num_reads"]]))
        except (ValueError, KeyError): unf[m.group(1)]=(0.0,0.0)

missed=expr-matched
reached_em=[e for e in missed if e in unf]           # (a) filtered out post-EM
pre_em=[e for e in missed if e not in unf]            # (b) pre-EM / not generated

print(f"# recall FN split — {S} | expressed={len(expr)} matched={len(matched&expr)} missed={len(missed)}")
print(f"# unfiltered.tsv GTF-source entries: {len(unf)}")
print(f"\n(a) reached EM, dropped by POST-EM filter : {len(reached_em)} ({100*len(reached_em)/len(missed):.1f}% of missed)")
print(f"(b) pre-EM / never generated / EM-zeroed  : {len(pre_em)} ({100*len(pre_em)/len(missed):.1f}% of missed)")

# attribute (a): the gtf abundance floor is 1.0 (min_gtf_abundance). abundance<1 -> killed by floor.
lt1=sum(1 for e in reached_em if unf[e][0] < 1.0)
ge1=len(reached_em)-lt1
print(f"\n## within (a) filtered-out:")
print(f"   abundance < 1  (killed by min_gtf_abundance floor=1): {lt1} ({100*lt1/max(len(reached_em),1):.1f}%)")
print(f"   abundance >= 1 (dropped by a LATER filter: isoform-frac/soft-mass/fulllen/polyA): {ge1}")

def bx(v):
    return "3-5" if v<5 else "5-10" if v<10 else "10-50" if v<50 else "50+"
print(f"\n## split by expression bin (missed only):")
print(f"   {'bin':6} {'missed':>7} {'filtered(a)':>11} {'preEM(b)':>9}")
tot=collections.Counter(); fa=collections.Counter(); fb=collections.Counter()
for e in missed:
    b=bx(est[e]); tot[b]+=1
    if e in unf: fa[b]+=1
    else: fb[b]+=1
for b in ["3-5","5-10","10-50","50+"]:
    print(f"   {b:6} {tot[b]:>7} {fa[b]:>11} {fb[b]:>9}")
