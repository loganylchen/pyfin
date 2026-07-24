import sys, pysam, re
from collections import defaultdict
BAM=sys.argv[1]; TRUTH="p00val/truth_full.gtf"; TMAP="p00val/gcL_ceil.pyfin_ceiling.gtf.tmap"
SLACK=20
# 1. truth chains: (chrom, introns) -> (tid, start, end)
tex=defaultdict(list); tmeta={}
for ln in open(TRUTH):
    if '\ttranscript\t' in ln:
        p=ln.split('\t'); m=re.search(r'transcript_id "([^"]+)"',p[8])
        tmeta[m.group(1)]=(p[0],int(p[3])-1,int(p[4]))
    elif '\texon\t' in ln:
        p=ln.split('\t'); m=re.search(r'transcript_id "([^"]+)"',p[8])
        tex[m.group(1)].append((int(p[3])-1,int(p[4])))
truth_by_chain={}
for tid,ex in tex.items():
    ex=sorted(ex); chain=tuple((ex[i][1],ex[i+1][0]) for i in range(len(ex)-1))
    ch,s,e=tmeta[tid]
    truth_by_chain[(ch,chain)]=(tid,s,e)
sys.stderr.write(f"truth transcripts: {len(truth_by_chain)}\n")
# 2. scan BAM: reads with exact chain matching a truth tx AND spanning its full extent
b=pysam.AlignmentFile(BAM,"rb")
fl_count=defaultdict(int)  # tid -> #full-length reads
for r in b.fetch():
    if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
    introns=[]; pos=r.reference_start
    for op,ln in r.cigartuples:
        if op in (0,2,7,8): pos+=ln
        elif op==3: introns.append((pos,pos+ln)); pos+=ln
    key=(r.reference_name,tuple(introns))
    t=truth_by_chain.get(key)
    if t is None: continue
    tid,s,e=t
    if r.reference_start<=s+SLACK and r.reference_end>=e-SLACK:
        fl_count[tid]+=1
b.close()
findable={tid for tid,c in fl_count.items() if c>=2}
sys.stderr.write(f"truth tx with >=2 full-length reads: {len(findable)}\n")
# 3. ceiling '=' refs
eq=set()
with open(TMAP) as fh:
    h=fh.readline().split('\t'); ci=h.index('class_code'); ri=h.index('ref_id')
    for ln in fh:
        p=ln.rstrip('\n').split('\t')
        if p[ci]=='=' and p[ri] and p[ri]!='-': eq.add(p[ri])
found=findable & eq; miss=findable - eq
print(f"\n=== transcripts with >=2 FULL-LENGTH reads (definitely present in data) ===")
print(f"  total findable: {len(findable)}")
print(f"  GENERATED as '=' (found): {len(found)} ({100*len(found)/len(findable):.1f}%)")
print(f"  NOT generated (MISSED):   {len(miss)} ({100*len(miss)/len(findable):.1f}%)")
# also >=3, >=5
for k in (3,5,10):
    f={tid for tid,c in fl_count.items() if c>=k}
    if f: print(f"  [>= {k} full-len reads] findable={len(f)} found={len(f&eq)} ({100*len(f&eq)/len(f):.0f}%) missed={len(f-eq)}")

# 4. of the MISSED, how many are an EXACT sub-chain of a LONGER truth transcript
#    (i.e. folded away)? group truth by chrom -> list of (tid, chain)
by_chrom=defaultdict(list)
for (ch,chain),(tid,s,e) in truth_by_chain.items():
    by_chrom[ch].append((tid,chain))
missed_tids=miss
# build chrom -> set of chains present (any truth)
def is_exact_subchain(short,long_):
    n,m=len(long_),len(short)
    if m==0 or m>=n: return False
    return any(long_[o:o+m]==short for o in range(n-m+1))
miss_meta={}
for (ch,chain),(tid,s,e) in truth_by_chain.items():
    if tid in missed_tids: miss_meta[tid]=(ch,chain)
sub_of_longer=0; sub_of_found=0
for tid,(ch,chain) in miss_meta.items():
    hit=False; hitfound=False
    for (tid2,chain2) in by_chrom[ch]:
        if tid2==tid: continue
        if is_exact_subchain(chain,chain2):
            hit=True
            if tid2 in found: hitfound=True
    if hit: sub_of_longer+=1
    if hitfound: sub_of_found+=1
print(f"\n=== root cause of the {len(missed_tids)} missed (>=2 full-len) ===")
print(f"  exact sub-chain of a LONGER truth tx (folded-away shape): {sub_of_longer} ({100*sub_of_longer/len(missed_tids):.0f}%)")
print(f"    ...of which the longer tx WAS generated as '=': {sub_of_found} ({100*sub_of_found/len(missed_tids):.0f}%)")
print(f"  NOT a sub-chain (missed for other reasons): {len(missed_tids)-sub_of_longer} ({100*(len(missed_tids)-sub_of_longer)/len(missed_tids):.0f}%)")
