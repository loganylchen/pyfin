import sys, json, pysam
import numpy as np
BAM=sys.argv[1]
cases=json.load(open("missed340.json")); b=pysam.AlignmentFile(BAM,"rb")
def maxsc(r):
    c=r.cigartuples
    return max(c[0][1] if c[0][0]==4 else 0, c[-1][1] if c[-1][0]==4 else 0)
# 1. of the 126 fusion-caused missed: max softclip across their full-length reads
miss_sc=[]
for case in cases:
    ch=case["chrom"]; s=case["start"]; e=case["end"]
    truth_chain=tuple((a,bb) for a,bb in case["chain"])
    fl=[]
    for r in b.fetch(ch,max(0,s-200),e+200):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        introns=[]; pos=r.reference_start
        for op,ln in r.cigartuples:
            if op in (0,2,7,8): pos+=ln
            elif op==3: introns.append((pos,pos+ln)); pos+=ln
        if tuple(introns)==truth_chain and r.reference_start<=s+20 and r.reference_end>=e-20:
            fl.append(maxsc(r))
    if not fl: continue
    nonfus50=sum(1 for x in fl if x<50)
    if nonfus50<2:  # fusion-caused
        miss_sc.append(min(fl))  # the SMALLEST softclip among its full-len reads (threshold must exceed this to recover >=2? actually need >=2 with sc<thr)
        # better: the 2nd-smallest, since we need >=2 reads under threshold
b.close()
sc=np.array(sorted(miss_sc))
print(f"126 fusion-caused missed: min-softclip per tx (need threshold > this to un-flag its full-len read)")
if len(sc):
    for thr in (100,150,200,250,300):
        print(f"  threshold {thr}: recovers {100*np.mean(sc<thr):.0f}% of them (their min full-len softclip < {thr})")
    print(f"  min-softclip distribution: median={np.median(sc):.0f} p90={np.percentile(sc,90):.0f} max={sc.max()}")
