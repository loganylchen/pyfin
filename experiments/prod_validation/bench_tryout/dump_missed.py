import sys, pysam, re, json
from collections import defaultdict
BAM=sys.argv[1]; TRUTH="p00val/truth_full.gtf"; TMAP="p00val/gcL_ceil.pyfin_ceiling.gtf.tmap"; SLACK=20
tex=defaultdict(list); tmeta={}
for ln in open(TRUTH):
    if '\ttranscript\t' in ln:
        p=ln.split('\t'); m=re.search(r'transcript_id "([^"]+)"',p[8]); tmeta[m.group(1)]=(p[0],int(p[3])-1,int(p[4]))
    elif '\texon\t' in ln:
        p=ln.split('\t'); m=re.search(r'transcript_id "([^"]+)"',p[8]); tex[m.group(1)].append((int(p[3])-1,int(p[4])))
truth_by_chain={}; chain_of={}
for tid,ex in tex.items():
    ex=sorted(ex); chain=tuple((ex[i][1],ex[i+1][0]) for i in range(len(ex)-1)); ch,s,e=tmeta[tid]
    truth_by_chain[(ch,chain)]=(tid,s,e); chain_of[tid]=(ch,chain)
b=pysam.AlignmentFile(BAM,"rb"); fl_reads=defaultdict(list)
for r in b.fetch():
    if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
    introns=[]; pos=r.reference_start
    for op,ln in r.cigartuples:
        if op in (0,2,7,8): pos+=ln
        elif op==3: introns.append((pos,pos+ln)); pos+=ln
    t=truth_by_chain.get((r.reference_name,tuple(introns)))
    if t is None: continue
    tid,s,e=t
    if r.reference_start<=s+SLACK and r.reference_end>=e-SLACK: fl_reads[tid].append(r.query_name)
b.close()
findable={tid for tid,rl in fl_reads.items() if len(rl)>=2}
eq=set()
with open(TMAP) as fh:
    h=fh.readline().split('\t'); ci=h.index('class_code'); ri=h.index('ref_id')
    for ln in fh:
        p=ln.rstrip('\n').split('\t')
        if p[ci]=='=' and p[ri] and p[ri]!='-': eq.add(p[ri])
miss=findable-eq
by_chrom=defaultdict(list)
for (ch,chain),(tid,s,e) in truth_by_chain.items(): by_chrom[ch].append((tid,chain))
def is_sub(a,bb):
    n,m=len(bb),len(a)
    return m>0 and m<n and any(bb[o:o+m]==a for o in range(n-m+1))
out=[]
for tid in miss:
    ch,chain=chain_of[tid]
    subtruth=any(is_sub(chain,c2) for (t2,c2) in by_chrom[ch] if t2!=tid)
    if subtruth: continue   # keep only the non-subchain-of-truth ones (the 182)
    tsp,s,e=None,tmeta[tid][1],tmeta[tid][2]
    out.append({"tid":tid,"chrom":ch,"start":s,"end":e,"chain":list(chain),"reads":fl_reads[tid][:8],"nreads":len(fl_reads[tid])})
json.dump(out,open("missed182.json","w"))
print(f"dumped {len(out)} non-subchain-of-truth missed cases")
