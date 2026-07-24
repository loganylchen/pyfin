import os,sys
import mappy,pysam
REPO="/SSD/logan/dev/pyfin"; sys.path.insert(0,REPO)
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_gtf import GTFReader
from fin.candidates.discovery import discover_candidates
from fin.candidates.chain_cluster import _related,_exact_subchain,_cassette,_contains,_wobble
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"; S="SGNex_H9_directRNA_replicate2_run2"
BAM=f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
GENOME=f"{B}/refs/host/GRCh38.primary_assembly.genome.fa"
GTF=f"{B}/results/gencode_full_sweep/_ref/p00/annotation.gtf"
RID="c37c1e06-7a20-4582-924d-2ed5e920dd48"
fa=pysam.FastaFile(GENOME)
gtf=GTFReader(GTF); gtf.open(); gtf.parse()
intervals=generate_isolated_intervals(BAM,gtf_path=GTF,max_gap=0)["intervals"]
preset=get_m1_preset()
target=[iv for iv in intervals if iv.chrom=="chr1" and iv.start<=7985285<=iv.end]
for iv in target:
    cs=fa.fetch(iv.chrom)
    cset=discover_candidates(iv,BAM,gtf,cs,threshold=24,min_novel_reads=1,chain_cluster=True,canonical_search_bp=0)
    cl=list(cset.candidates)
    rs=getattr(cset,"read_sequences",{}) or {}
    if RID not in rs: continue
    seq=rs[RID]
    aln=[mappy.Aligner(seq=c.sequence,preset=preset) if c.sequence else None for c in cl]
    best=-1e18; row=[-1e18]*len(cl)
    for j,a in enumerate(aln):
        if a is None: continue
        b=None
        for h in a.map(seq):
            v=score_hit(h)
            if v is not None and (b is None or v>b): b=v
        if b is not None:
            row[j]=b
            if b>best: best=b
    tie=[j for j in range(len(cl)) if row[j]>0 and row[j]>=best-1e-9]
    print(f"interval {iv.chrom}:{iv.start}-{iv.end}  #candidates={len(cl)}  read tie size={len(tie)}")
    # per candidate: #reads it was built from
    for rank,j in enumerate(tie):
        c=cl[j]; ch=tuple(c.intron_chain.introns)
        nr=len(getattr(c,"read_ids",[]) or getattr(c,"support_reads",[]) or [])
        print(f"  tie[{rank}] AS={row[j]:.0f} nIntrons={len(ch)} span={c.start}-{c.end} src={c.source} nreads~{nr}")
    # pairwise: same chain_cluster?
    print("  --- pairwise chain_cluster relation among tied candidates ---")
    for a in range(len(tie)):
        for b in range(a+1,len(tie)):
            ca=tuple(cl[tie[a]].intron_chain.introns); cb=tuple(cl[tie[b]].intron_chain.introns)
            rel=_related(ca,cb,6,70)
            wob=_wobble(ca,cb,6)
            lo,hi=(ca,cb) if len(ca)<len(cb) else (cb,ca)
            cas=_cassette(hi,lo,6,70) if len(lo)<len(hi) else False
            con=_contains(lo,hi,6) if len(lo)<len(hi) else False
            exsub=_exact_subchain(lo,hi) if len(lo)<len(hi) else False
            print(f"   tie[{a}] vs tie[{b}]: related={rel} (wobble={wob} cassette={cas} contains={con} exact_sub={exsub}) |dIntrons|={abs(len(ca)-len(cb))}")
    break
