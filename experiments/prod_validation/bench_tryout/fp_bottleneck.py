import re,sys,statistics as st
from collections import defaultdict
import pysam
BAM=sys.argv[1]
GTF='p00val/pyfin_m2off_sc250.gtf'
TMAP='p00val/gcL_m2off_sc250.pyfin_m2off_sc250.gtf.tmap'

# candidate junctions (0-based half-open) + class + num_reads
ex=defaultdict(list); chrom={}; nr={}
for ln in open(GTF):
    if ln.startswith('#'): continue
    f=ln.split('\t')
    if len(f)<9: continue
    m=re.search('transcript_id "([^"]+)"',f[8])
    if not m: continue
    t=m.group(1)
    if f[2]=='exon': ex[t].append((int(f[3]),int(f[4]))); chrom[t]=f[0]
    elif f[2]=='transcript':
        rm=re.search('num_reads "([0-9]+)"',f[8]); nr[t]=int(rm.group(1)) if rm else 0
cls={}
for i,ln in enumerate(open(TMAP)):
    if i==0: continue
    f=ln.rstrip("\n").split("\t"); cls[f[4]]=f[2]

cand_junc={}   # t -> set of (chrom,a0,b0)
want=defaultdict(set)  # chrom -> set of (a0,b0) we care about
for t,e in ex.items():
    e=sorted(e)
    js=set()
    for i in range(len(e)-1):
        a0=e[i][1]; b0=e[i+1][0]-1   # 0-based half-open [a0,b0)
        js.add((chrom[t],a0,b0)); want[chrom[t]].add((a0,b0))
    cand_junc[t]=js

# one BAM pass: count read support per wanted junction
jcount=defaultdict(int)
b=pysam.AlignmentFile(BAM,'rb')
for c in want:
    ws=want[c]
    for r in b.fetch(c):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        pos=r.reference_start
        for op,ln in r.cigartuples or []:
            if op==3:
                key=(pos,pos+ln)
                if key in ws: jcount[(c,pos,pos+ln)]+=1
                pos+=ln
            elif op in (0,2,7,8): pos+=ln
b.close()

# per candidate: bottleneck = min junction support; ratio = bottleneck/num_reads
rows=defaultdict(list)
for t,js in cand_junc.items():
    if not js: continue
    supp=[jcount.get(j,0) for j in js]
    bottleneck=min(supp)
    rows[cls.get(t,'?')].append((bottleneck, nr.get(t,0)))

print(f"{'class':6}{'n':>6}{'medBottleneck':>14}{'%bottleneck<2':>14}{'%bottleneck<=1':>15}")
for cl in ['=','c','j','k','m','o','n','u','i']:
    if cl not in rows: continue
    v=rows[cl]; bn=[x[0] for x in v]
    print(f"{cl:6}{len(v):6}{st.median(bn):14.0f}{100*sum(1 for x in bn if x<2)/len(bn):13.0f}%{100*sum(1 for x in bn if x<=1)/len(bn):14.0f}%")
# TP vs FP
tp=[x[0] for x in rows.get('=',[])]; fp=[x[0] for cl in rows if cl!='=' for x in rows[cl]]
print(f"\nTP bottleneck: median={st.median(tp):.0f} %(<2)={100*sum(1 for x in tp if x<2)/len(tp):.0f}")
print(f"FP bottleneck: median={st.median(fp):.0f} %(<2)={100*sum(1 for x in fp if x<2)/len(fp):.0f}")
