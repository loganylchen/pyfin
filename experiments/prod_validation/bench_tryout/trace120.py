import sys, json, pysam
sys.path.insert(0,"/SSD/logan/dev/pyfin")
from fin.io.interval_manager import GenomicInterval
from fin.candidates.discovery import discover_candidates
from fin.candidates.chain_cluster import cluster_read_chains
from fin.candidates.intron_chains import extract_intron_chain
from collections import Counter
BAM=sys.argv[1]; GENOME=f"{sys.argv[2]}/genome.fa"
fa=pysam.FastaFile(GENOME); _cc={}
def cs(c):
    if c not in _cc:
        try:_cc[c]=fa.fetch(c)
        except Exception:_cc[c]=""
    return _cc[c]
cases=json.load(open("missed340.json"))
b=pysam.AlignmentFile(BAM,"rb")
# reselect the member-local ones
V=Counter(); samples=[]
for case in cases:
    tid=case["tid"]; ch=case["chrom"]; s=case["start"]; e=case["end"]
    truth_chain=tuple((a,bb) for a,bb in case["chain"])
    rc=[]
    for r in b.fetch(ch,max(0,s-200),e+200):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        rc.append(({"query_name":r.query_name}, extract_intron_chain(r.cigartuples,r.reference_start)))
    cl=cluster_read_chains(rc)
    members=[tuple(m.chain.introns) for c2 in cl for m in c2.members]
    if truth_chain not in members: continue     # only the member-local ones
    # run real discover_candidates on this interval
    strand="+"  # infer from any read
    iv=GenomicInterval(chrom=ch,start=max(0,s-200),end=e+200,strand=strand)
    try:
        cset=discover_candidates(iv,BAM,None,cs(ch),threshold=24,min_novel_reads=1,
                                 chain_cluster=True,canonical_search_bp=0)
    except Exception as ex:
        V[f"discover_error:{type(ex).__name__}"]+=1; continue
    cand_chains=[tuple(c.intron_chain.introns) for c in cset.candidates]
    if truth_chain in cand_chains:
        V["candidate_PRESENT(so lost later, not generation)"]+=1
    else:
        V["candidate_ABSENT(dropped in _chain_cluster_candidates)"]+=1
        if len(samples)<3:
            # why? check if it's a member but no seq / span
            samples.append((tid,ch,s,e,len(truth_chain)))
b.close()
print("=== the member-local cases: does real discover_candidates emit the chain? ===")
for k,v in V.most_common(): print(f"  {v:4}  {k}")
print("samples(absent):",samples)
