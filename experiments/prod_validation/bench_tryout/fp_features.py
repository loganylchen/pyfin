import re,sys,statistics as st
from collections import defaultdict
import pysam
GEN=sys.argv[1]
GTF='p00val/pyfin_m2off_sc250.gtf'
TMAP='p00val/gcL_m2off_sc250.pyfin_m2off_sc250.gtf.tmap'
fa=pysam.FastaFile(GEN)
chrset=set(fa.references)

# parse candidates: exons, strand, num_reads
ex=defaultdict(list); strand={}; nr={}; chrom={}
for ln in open(GTF):
    if ln.startswith('#'): continue
    f=ln.split('\t')
    if len(f)<9: continue
    m=re.search('transcript_id "([^"]+)"',f[8])
    if not m: continue
    t=m.group(1)
    if f[2]=='exon': ex[t].append((int(f[3]),int(f[4]))); strand[t]=f[6]; chrom[t]=f[0]
    elif f[2]=='transcript':
        rm=re.search('num_reads "([0-9]+)"',f[8]); nr[t]=int(rm.group(1)) if rm else 0

cls={}
for i,ln in enumerate(open(TMAP)):
    if i==0: continue
    f=ln.rstrip("\n").split("\t"); cls[f[4]]=f[2]

def rc(s):
    return s.translate(str.maketrans('ACGTacgt','TGCAtgca'))[::-1]
def canon_frac(t):
    e=sorted(ex[t]); c=chrom[t]
    if len(e)<2 or c not in chrset: return None
    good=0; tot=0
    for i in range(len(e)-1):
        istart=e[i][1]+1; iend=e[i+1][0]-1   # 1-based inclusive intron
        if iend-istart+1<4: continue
        d=fa.fetch(c,istart-1,istart+1).upper()      # first 2 bases
        a=fa.fetch(c,iend-2,iend).upper()            # last 2 bases
        if strand[t]=='-': d,a=rc(a),rc(d)
        tot+=1
        if (d,a) in (('GT','AG'),('GC','AG'),('AT','AC')): good+=1
    return good/tot if tot else None

# aggregate by class
bycls_frac=defaultdict(list); bycls_allcanon=defaultdict(int); bycls_n=defaultdict(int)
for t,cl in cls.items():
    if t not in ex: continue
    cf=canon_frac(t)
    if cf is None: continue
    bycls_n[cl]+=1
    bycls_frac[cl].append(cf)
    if cf>=0.999: bycls_allcanon[cl]+=1

print(f"{'class':6}{'n':>6}{'medianCanonFrac':>16}{'%allCanonical':>15}{'medReads':>9}")
order=['=','c','j','k','m','o','n','e','u','i','p','x']
for cl in order:
    if cl not in bycls_n: continue
    fr=bycls_frac[cl]; rds=[nr.get(t,0) for t in cls if cls[t]==cl and t in ex]
    print(f"{cl:6}{bycls_n[cl]:6}{st.median(fr):16.2f}{100*bycls_allcanon[cl]/bycls_n[cl]:14.0f}%{st.median(rds):9.0f}")

# TP vs FP overall
tp=[cf for t,cl in cls.items() if cl=='=' and t in ex and (cf:=canon_frac(t)) is not None]
fp=[cf for t,cl in cls.items() if cl!='=' and t in ex and (cf:=canon_frac(t)) is not None]
print(f"\nTP all-canonical: {100*sum(1 for x in tp if x>=0.999)/len(tp):.0f}%   FP all-canonical: {100*sum(1 for x in fp if x>=0.999)/len(fp):.0f}%")
