import sys, json, pysam
sys.path.insert(0,"/SSD/logan/dev/pyfin")
from fin.io.interval_manager import GenomicInterval, is_fusion_read
from fin.candidates.discovery import discover_candidates, _build_spliced_sequence, _exons_from_chain
from fin.candidates.chain_cluster import cluster_read_chains
from fin.candidates.intron_chains import extract_intron_chain
from fin.candidates.dataclasses import IntronChain
from fin.io.io_bam import get_reads_in_region
from collections import Counter
BAM=sys.argv[1]; GENOME=f"{sys.argv[2]}/genome.fa"
fa=pysam.FastaFile(GENOME); _cc={}
def cs(c):
    if c not in _cc:
        try:_cc[c]=fa.fetch(c)
        except Exception:_cc[c]=""
    return _cc[c]
cases=json.load(open("missed340.json")); b=pysam.AlignmentFile(BAM,"rb")
V=Counter()
for case in cases:
    tid=case["tid"]; ch=case["chrom"]; s=case["start"]; e=case["end"]
    truth_chain=tuple((a,bb) for a,bb in case["chain"])
    rc=[]
    for r in b.fetch(ch,max(0,s-200),e+200):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        rc.append(({"query_name":r.query_name}, extract_intron_chain(r.cigartuples,r.reference_start)))
    cl=cluster_read_chains(rc)
    if truth_chain not in [tuple(m.chain.introns) for c2 in cl for m in c2.members]: continue
    iv=GenomicInterval(chrom=ch,start=max(0,s-200),end=e+200,strand="+")
    cset=discover_candidates(iv,BAM,None,cs(ch),threshold=24,min_novel_reads=1,chain_cluster=True,canonical_search_bp=0)
    if truth_chain in [tuple(c.intron_chain.introns) for c in cset.candidates]: continue  # present -> skip
    # ABSENT. why? re-run the internal path: fetch non_fusion reads, cluster, check member -> span/seq
    reads=get_reads_in_region(BAM, ch, max(0,s-200), e+200)
    n_fusion=sum(1 for rd in reads if is_fusion_read(rd))
    # build spans as _chain_cluster_candidates does (non-fusion only)
    spans={}; nfrc=[]
    for rd in reads:
        if is_fusion_read(rd): continue
        cig=rd.get("cigartuples"); rs=rd.get("reference_start")
        if cig is None or rs is None: continue
        q=rd.get("query_name")
        if q: spans[q]=(rs, rd.get("reference_end",rs))
        nfrc.append((rd, extract_intron_chain(cig,rs)))
    cl2=cluster_read_chains(nfrc)
    mem=None
    for c2 in cl2:
        for m in c2.members:
            if tuple(m.chain.introns)==truth_chain: mem=m
    if mem is None:
        V["chain_gone_after_fusion_filter(reads were fusion)"]+=1; continue
    starts=[spans[r][0] for r in mem.read_ids if r in spans]
    if not starts:
        V["no_span(reads not in spans)"]+=1; continue
    st,en=min(starts),max(spans[r][1] for r in mem.read_ids if r in spans)
    seq=_build_spliced_sequence(cs(ch), st, en, IntronChain(introns=truth_chain), "+")
    if not seq: V["seq_build_failed"]+=1
    else: V["seq_ok_but_still_absent(??)"]+=1
b.close()
print("=== why the 75 candidate-absent are dropped ===")
for k,v in V.most_common(): print(f"  {v:4}  {k}")
