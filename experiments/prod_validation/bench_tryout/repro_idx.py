import sys
sys.path.insert(0,"/SSD/logan/dev/pyfin")
import pysam
from fin.io.interval_manager import generate_isolated_intervals
from fin.candidates.discovery import discover_candidates
BAM="smoke/smoke.bam"
GENOME="/SSD/logan/dev/pyfin/experiments/prod_validation/gencode/SGNex_H9_directRNA_replicate2_run2/stage/genome.fa"
fa=pysam.FastaFile(GENOME); _cc={}
def cs(c):
    if c not in _cc:
        try:_cc[c]=fa.fetch(c)
        except Exception:_cc[c]=""
    return _cc[c]
ivs=generate_isolated_intervals(BAM,gtf_path=None,max_gap=0)["intervals"]
bad=0
for iv in ivs:
    cset=discover_candidates(iv,BAM,None,cs(iv.chrom),threshold=24,min_novel_reads=1,
                             chain_cluster=True,canonical_search_bp=0)
    n=len(cset.candidates)
    cl=cset.clusters
    if cl is None: continue
    for gi,idxs in enumerate(cl):
        for i in idxs:
            if i>=n or i<0:
                bad+=1
                print(f"BAD {iv.region_string}: cluster {gi} idx {i} >= ncand {n} (nclusters={len(cl)})")
                if bad<=5:
                    print("   cluster idxs:",idxs)
if bad==0: print("all cluster indices in range")
else: print(f"total bad: {bad}")
