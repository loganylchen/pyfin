import sys, json, pysam
sys.path.insert(0,"/SSD/logan/dev/pyfin")
from fin.candidates.chain_cluster import cluster_read_chains, _exact_subchain, _related
from fin.candidates.intron_chains import extract_intron_chain
BAM=sys.argv[1]
cases=json.load(open("missed182.json"))
b=pysam.AlignmentFile(BAM,"rb")
from collections import Counter
verdict=Counter()
for case in cases:
    tid=case["tid"]; ch=case["chrom"]; s=case["start"]; e=case["end"]
    truth_chain=tuple((a,b_) for a,b_ in case["chain"])
    # fetch reads overlapping the transcript region (generous flank)
    rc=[]
    for r in b.fetch(ch,max(0,s-200),e+200):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        rc.append(({"query_name":r.query_name}, extract_intron_chain(r.cigartuples,r.reference_start)))
    clusters=cluster_read_chains(rc)
    # is truth_chain a member of any cluster?
    members=[tuple(m.chain.introns) for cl in clusters for m in cl.members]
    if truth_chain in members:
        verdict["member_generated(should be =)"]+=1
        continue
    # is it a folded shadow of some member?
    folded=[tuple(f.chain.introns) for cl in clusters for m in cl.members for f in m.folded]
    if truth_chain in folded:
        # find container
        cont=None
        for cl in clusters:
            for m in cl.members:
                if any(tuple(f.chain.introns)==truth_chain for f in m.folded):
                    cont=tuple(m.chain.introns)
        exact = cont is not None and _exact_subchain(truth_chain, cont)
        verdict["folded_shadow"]+=1
        continue
    # present in reads at all?
    read_chains=[tuple(c.introns) for _,c in rc]
    if truth_chain in read_chains:
        verdict["in_reads_but_not_member_not_folded(??)"]+=1
    else:
        verdict["truth_chain_NOT_in_any_read_here"]+=1
b.close()
print("=== what happened to the 182 truth chains in cluster_read_chains ===")
for k,v in verdict.most_common(): print(f"  {v:4}  {k}")
