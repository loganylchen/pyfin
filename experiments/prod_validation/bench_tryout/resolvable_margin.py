#!/usr/bin/env python3
"""Run M2 eventalign ONLY on the structurally-resolvable tie reads (those whose
aligned span crosses a junction that DIFFERS among the tie candidates) and report
their NLL margin distribution. Answers: for the reads that CAN in principle be
distinguished, does the signal actually distinguish them, or is NLL still tied?
"""
import os, sys
import numpy as np
import mappy, pysam

REPO="/SSD/logan/dev/pyfin"; sys.path.insert(0,REPO)
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_gtf import GTFReader
from fin.candidates.discovery import discover_candidates
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit
from fin.scoring.m2_junction_nll import class_junction_window_set, read_cand_mean_nll
from fin.scoring.krill_aligner import make_krill_aligner, krill_thread_count
import krill

B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S="SGNex_H9_directRNA_replicate2_run2"
BAM=f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
GENOME=f"{B}/refs/host/GRCh38.primary_assembly.genome.fa"
GTF=f"{B}/results/gencode_full_sweep/_ref/p00/annotation.gtf"
SIGNAL="/autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/blow5/nanopore.blow5"
LIMIT=int(os.environ.get("LIMIT","1500"))
USE_GPU=os.environ.get("USE_GPU","1")=="1"

fa=pysam.FastaFile(GENOME); _cc={}
def cseq(c):
    if c not in _cc:
        try:_cc[c]=fa.fetch(c)
        except Exception:_cc[c]=""
    return _cc[c]
bam=pysam.AlignmentFile(BAM,"rb")
gtf=GTFReader(GTF); gtf.open(); gtf.parse()
intervals=generate_isolated_intervals(BAM,gtf_path=GTF,max_gap=0)["intervals"]
if LIMIT: intervals=intervals[:LIMIT]
sys.stderr.write(f"intervals:{len(intervals)} gpu={USE_GPU}\n")
preset=get_m1_preset()
nthr=krill_thread_count()
ka,eff=make_krill_aligner(krill,"rna002",USE_GPU,hmm_confidence=False,num_thread=nthr)

def spanmap(chrom,start,end):
    m={}
    for r in bam.fetch(chrom,max(0,start-10),end+10):
        if r.is_supplementary or r.is_secondary or r.is_unmapped: continue
        m.setdefault(r.query_name,(r.reference_start,r.reference_end))
    return m

margins=[]; n_res=0; n_scored=0
for n,iv in enumerate(intervals):
    cs=cseq(iv.chrom)
    cset=discover_candidates(iv,BAM,gtf,cs,threshold=24,min_novel_reads=1,
                             chain_cluster=True,canonical_search_bp=0)
    cl=list(cset.candidates)
    if len(cl)<2: continue
    rs=getattr(cset,"read_sequences",{}) or {}
    aln=[mappy.Aligner(seq=c.sequence,preset=preset) if c.sequence else None for c in cl]
    spans=spanmap(iv.chrom,iv.start,iv.end)
    for rid,seq in rs.items():
        if not seq or rid not in spans: continue
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
        if best<=0: continue
        tie=[j for j in range(len(cl)) if row[j]>0 and row[j]>=best-1e-9]
        if len(tie)<2: continue
        jsets=[set(cl[j].intron_chain.introns) for j in tie]
        disc=set().union(*jsets)-set(jsets[0]).intersection(*jsets)
        rst,ren=spans[rid]
        if not any(rst<=s and e<=ren for (s,e) in disc):
            continue                      # containment -> skip (measured already)
        n_res+=1
        tc=[cl[j] for j in tie]; ta=[aln[j] for j in tie]
        gset=class_junction_window_set(tc,flank=2,k=10)
        if not gset: continue
        sc=[]
        for idx,c in enumerate(tc):
            nll,nev=read_cand_mean_nll(rid,seq,c,[],ka,ta[idx],SIGNAL,"rna002",
                                       gset=gset,use_gpu=eff,num_thread=nthr)
            if nev>0 and np.isfinite(nll): sc.append(nll)
        if len(sc)<2: continue
        sc.sort(); n_scored+=1; margins.append(sc[1]-sc[0])
    if (n+1)%200==0 and margins:
        mm=np.array(margins)
        sys.stderr.write(f"  [{n+1}] res={n_res} scored={n_scored} | "
          f"median={np.median(mm):.3f} >=0.1:{100*np.mean(mm>=0.1):.0f}% "
          f">=0.5:{100*np.mean(mm>=0.5):.0f}% >=1:{100*np.mean(mm>=1):.0f}%\n"); sys.stderr.flush()

m=np.array(margins) if margins else np.array([0.0])
print(f"\n=== M2 margin on RESOLVABLE ties only (read spans a differing junction) ===")
print(f"resolvable reads: {n_res}   scored(>=2 NLL): {n_scored}   margins: {len(margins)}")
print(f"  median={np.median(m):.3f} mean={np.mean(m):.3f} p25={np.percentile(m,25):.3f} p75={np.percentile(m,75):.3f} max={m.max():.2f}")
for t in (0.05,0.1,0.2,0.5,1.0,2.0):
    print(f"  margin >= {t:<4}: {100*np.mean(m>=t):.0f}%")
bam.close()
