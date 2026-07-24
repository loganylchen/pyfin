#!/usr/bin/env python3
"""Within-cluster M1: for each read's TOP-2 members, record the AS gap (best-2nd) and
whether the top-2 are a wobble pair (same intron count, every junction within 6bp).
Tells us how big an m1_tie_margin captures wobble near-ties as ties without merging
genuinely-different candidates. Signal-free, CPU."""
import os, sys, logging
logging.disable(logging.INFO)
import numpy as np
import mappy, pysam
REPO="/SSD/logan/dev/pyfin"; sys.path.insert(0,REPO)
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_gtf import GTFReader
from fin.candidates.discovery import discover_candidates
from fin.candidates.chain_cluster import _wobble, _contains, _cassette
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"; S="SGNex_H9_directRNA_replicate2_run2"
BAM=os.environ.get("BAM_OVR", f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam")
GENOME=f"{B}/refs/host/GRCh38.primary_assembly.genome.fa"
GTF=f"{B}/results/gencode_full_sweep/_ref/p00/annotation.gtf"
LIMIT=int(os.environ.get("LIMIT","3000"))
fa=pysam.FastaFile(GENOME); _cc={}
def cs(c):
    if c not in _cc:
        try:_cc[c]=fa.fetch(c)
        except Exception:_cc[c]=""
    return _cc[c]
gtf=GTFReader(GTF); gtf.open(); gtf.parse()
ivs=generate_isolated_intervals(BAM,gtf_path=GTF,max_gap=0)["intervals"]
if LIMIT: ivs=ivs[:LIMIT]
sys.stderr.write(f"intervals:{len(ivs)}\n")
preset=get_m1_preset()
gap_wob=[]; gap_other=[]; n_reads=0
for n,iv in enumerate(ivs):
    cset=discover_candidates(iv,BAM,gtf,cs(iv.chrom),threshold=24,min_novel_reads=1,
                             chain_cluster=True,canonical_search_bp=0)
    cl=list(cset.candidates); rs=getattr(cset,"read_sequences",{}) or {}
    clusters=cset.clusters
    if not clusters: continue
    id2i={c.candidate_id:i for i,c in enumerate(cl)}
    for ids in clusters:
        idx=[id2i[c] for c in ids if c in id2i]
        if len(idx)<2: continue
        cands=[cl[i] for i in idx]
        aln=[mappy.Aligner(seq=c.sequence,preset=preset) if c.sequence else None for c in cands]
        chains=[tuple(c.intron_chain.introns) for c in cands]
        crids=set()
        for c in cands: crids|=c.supporting_read_ids
        for rid in crids:
            seq=rs.get(rid)
            if not seq: continue
            row=[]
            for j,a in enumerate(aln):
                if a is None: continue
                b=None
                for h in a.map(seq):
                    v=score_hit(h)
                    if v is not None and (b is None or v>b): b=v
                if b is not None and b>0: row.append((b,j))
            if len(row)<2: continue
            row.sort(reverse=True)
            (b1,j1),(b2,j2)=row[0],row[1]
            gap=b1-b2; n_reads+=1
            c1,c2=chains[j1],chains[j2]
            if _wobble(c1,c2,6): gap_wob.append(gap)
            else: gap_other.append(gap)
    if (n+1)%500==0: sys.stderr.write(f"  {n+1}/{len(ivs)} reads={n_reads} wob_pairs={len(gap_wob)}\n")
def rep(name,v):
    if not v: print(f"{name}: none"); return
    v=np.array(v)
    ps=[np.percentile(v,p) for p in (50,75,90,95,99)]
    print(f"{name}: n={len(v)} median={ps[0]:.0f} p75={ps[1]:.0f} p90={ps[2]:.0f} p95={ps[3]:.0f} p99={ps[4]:.0f} "
          f"| <=8:{100*np.mean(v<=8):.0f}% <=12:{100*np.mean(v<=12):.0f}% <=16:{100*np.mean(v<=16):.0f}% <=20:{100*np.mean(v<=20):.0f}% <=30:{100*np.mean(v<=30):.0f}%")
print(f"\n=== TOP-2 AS gap (best-2nd), within-cluster, {n_reads} reads ===")
rep("wobble top-2 pair", gap_wob)
rep("non-wobble top-2 ", gap_other)
