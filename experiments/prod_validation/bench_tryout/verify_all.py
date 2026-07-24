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
READS={"961e2771-c34d-4766-9dfd-0c53eb08f33f":6820250,"c37c1e06-7a20-4582-924d-2ed5e920dd48":7985285,
 "0caf98c4-621b-451b-bc34-befb7724ec88":26429735,"c9355cd8-d3e1-4de9-9193-c4fbfd6ef692":26922492,
 "4f32fcb6-2f3a-4aea-918b-8339cf44d959":28506084}
fa=pysam.FastaFile(GENOME)
gtf=GTFReader(GTF); gtf.open(); gtf.parse()
intervals=generate_isolated_intervals(BAM,gtf_path=GTF,max_gap=0)["intervals"]
preset=get_m1_preset()
for rid,anchor in READS.items():
    iv=next((iv for iv in intervals if iv.chrom=="chr1" and iv.start<=anchor<=iv.end),None)
    if iv is None: print(f"{rid[:8]}: no interval"); continue
    cs=fa.fetch(iv.chrom)
    cset=discover_candidates(iv,BAM,gtf,cs,threshold=24,min_novel_reads=1,chain_cluster=True,canonical_search_bp=0)
    cl=list(cset.candidates); rs=getattr(cset,"read_sequences",{}) or {}
    if rid not in rs: print(f"{rid[:8]}: read not in cset"); continue
    seq=rs[rid]
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
    print(f"\n{rid[:8]}  tie={len(tie)}  (interval {iv.start}-{iv.end}, {len(cl)} cands)")
    for j in tie:
        c=cl[j]; print(f"   nIntr={len(c.intron_chain.introns):2d} span={c.start}-{c.end}")
    allrel=True
    for a in range(len(tie)):
        for b in range(a+1,len(tie)):
            ca=tuple(cl[tie[a]].intron_chain.introns); cb=tuple(cl[tie[b]].intron_chain.introns)
            rel=_related(ca,cb,6,70)
            lo,hi=(ca,cb) if len(ca)<len(cb) else (cb,ca)
            ex=_exact_subchain(lo,hi) if len(lo)<len(hi) else False
            co=_contains(lo,hi,6) if len(lo)<len(hi) else False
            if not rel: allrel=False
            print(f"     [{a}]x[{b}] related={rel} exact_sub={ex} contains={co}")
    print(f"   => {'ALL pairs same-cluster (containment/wobble)' if allrel else 'has CROSS-cluster pair (M1 spurious tie)'}")
