import re,sys,statistics as st
from collections import defaultdict
import pysam
BAM=sys.argv[1]
GTF='p00val/pyfin_m2off_sc250.gtf'; TMAP='p00val/gcL_m2off_sc250.pyfin_m2off_sc250.gtf.tmap'
ex=defaultdict(list); chrom={}
for ln in open(GTF):
    if ln.startswith('#'): continue
    f=ln.split('\t')
    if len(f)<9 or f[2]!='exon': continue
    m=re.search('transcript_id "([^"]+)"',f[8])
    if not m: continue
    ex[m.group(1)].append((int(f[3]),int(f[4]))); chrom[m.group(1)]=f[0]
cls={}
for i,ln in enumerate(open(TMAP)):
    if i==0: continue
    f=ln.rstrip("\n").split("\t"); cls[f[4]]=f[2]
# candidate adjacent junction pairs (0-based half-open junctions)
cand_pairs={}; want=defaultdict(set)
for t,e in ex.items():
    e=sorted(e); js=[(e[i][1],e[i+1][0]-1) for i in range(len(e)-1)]
    if len(js)<2: continue
    pairs=set()
    for i in range(len(js)-1):
        p=(js[i],js[i+1]); pairs.add(p); want[chrom[t]].add(p)
    cand_pairs[t]=(pairs)
# BAM pass: count reads carrying each adjacent junction pair
pcount=defaultdict(int)
b=pysam.AlignmentFile(BAM,'rb')
for c in want:
    ws=want[c]
    for r in b.fetch(c):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        pos=r.reference_start; js=[]
        for op,ln in r.cigartuples or []:
            if op==3: js.append((pos,pos+ln)); pos+=ln
            elif op in (0,2,7,8): pos+=ln
        for i in range(len(js)-1):
            p=(js[i],js[i+1])
            if p in ws: pcount[(c,p)]+=1
b.close()
rows=defaultdict(list)
for t,pairs in cand_pairs.items():
    supp=[pcount.get((chrom[t],p),0) for p in pairs]
    rows[cls.get(t,'?')].append(min(supp))
print(f"{'class':6}{'n(>=2junc)':>11}{'medMinLink':>12}{'%hasBREAK(link=0)':>18}{'%weaklink<2':>13}")
for cl in ['=','c','j','k','m','n','u','i']:
    if cl not in rows: continue
    v=rows[cl]
    print(f"{cl:6}{len(v):11}{st.median(v):12.0f}{100*sum(1 for x in v if x==0)/len(v):17.0f}%{100*sum(1 for x in v if x<2)/len(v):12.0f}%")
tp=rows.get('=',[]); j=rows.get('j',[])
print(f"\nTP: %break={100*sum(1 for x in tp if x==0)/len(tp):.0f}  |  j-recombinant: %break={100*sum(1 for x in j if x==0)/len(j):.0f}  %weaklink<2={100*sum(1 for x in j if x<2)/len(j):.0f}")
