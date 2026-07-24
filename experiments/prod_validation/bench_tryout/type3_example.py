#!/usr/bin/env python3
"""Isolate 'type-3' blind spots: a tie read whose TOP-2 lowest-NLL candidates are
STRUCTURALLY DIFFERENT within the read's covered span (not a 2-3bp wobble: differing
count, or a junction >10bp off / present-vs-absent) yet the NLL margin ~ 0. Prints
each candidate's local structure, mean NLL, and #events scored (to expose whether the
blind spot is a degenerate low-event window). Early-exits after N examples.
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
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"; S="SGNex_H9_directRNA_replicate2_run2"
BAM=f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
GENOME=f"{B}/refs/host/GRCh38.primary_assembly.genome.fa"
GTF=f"{B}/results/gencode_full_sweep/_ref/p00/annotation.gtf"
SIGNAL="/autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/blow5/nanopore.blow5"
LIMIT=int(os.environ.get("LIMIT","1200")); USE_GPU=os.environ.get("USE_GPU","1")=="1"
N_WANT=int(os.environ.get("N_WANT","8")); MAXM=float(os.environ.get("MAXM","0.02")); WOB=6
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
preset=get_m1_preset(); nthr=krill_thread_count()
ka,eff=make_krill_aligner(krill,"rna002",USE_GPU,hmm_confidence=False,num_thread=nthr)
def spanmap(chrom,s,e):
    m={}
    for r in bam.fetch(chrom,max(0,s-10),e+10):
        if r.is_supplementary or r.is_secondary or r.is_unmapped: continue
        m.setdefault(r.query_name,(r.reference_start,r.reference_end))
    return m
def local(chain,lo,hi):
    return [(s,e) for s,e in chain if not(e<lo or s>hi)]
def wobble_equal(a,b):
    if len(a)!=len(b): return False
    return all(abs(s1-s2)<=WOB and abs(e1-e2)<=WOB for (s1,e1),(s2,e2) in zip(a,b))
got=0
for n,iv in enumerate(intervals):
    if got>=N_WANT: break
    cs=cseq(iv.chrom)
    cset=discover_candidates(iv,BAM,gtf,cs,threshold=24,min_novel_reads=1,chain_cluster=True,canonical_search_bp=0)
    cl=list(cset.candidates)
    if len(cl)<2: continue
    rs=getattr(cset,"read_sequences",{}) or {}
    aln=[mappy.Aligner(seq=c.sequence,preset=preset) if c.sequence else None for c in cl]
    spans=spanmap(iv.chrom,iv.start,iv.end)
    for rid,seq in rs.items():
        if got>=N_WANT: break
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
        tc=[cl[j] for j in tie]
        jsets=[set(c.intron_chain.introns) for c in tc]
        disc=set().union(*jsets)-set(jsets[0]).intersection(*jsets)
        rst,ren=spans[rid]
        spanned=[(s,e) for (s,e) in disc if rst<=s and e<=ren]
        if not spanned: continue
        gset=class_junction_window_set(tc,flank=2,k=10)
        if not gset: continue
        ta=[aln[j] for j in tie]; sc=[]
        for idx,c in enumerate(tc):
            nll,nev=read_cand_mean_nll(rid,seq,c,[],ka,ta[idx],SIGNAL,"rna002",gset=gset,use_gpu=eff,num_thread=nthr)
            sc.append((idx,nll,nev))
        fin=[(i,v,ne) for i,v,ne in sc if ne>0 and np.isfinite(v)]
        if len(fin)<2: continue
        fin.sort(key=lambda t:t[1]); margin=fin[1][1]-fin[0][1]
        if margin>MAXM: continue
        i0,i1=fin[0][0],fin[1][0]
        lo=min(s for s,e in spanned)-40; hi=max(e for s,e in spanned)+40
        l0=local(tc[i0].intron_chain.introns,lo,hi); l1=local(tc[i1].intron_chain.introns,lo,hi)
        if wobble_equal(l0,l1): continue    # top-2 identical/wobble in coverage -> type-2, skip
        got+=1
        print(f"\n===== type3 #{got}  read {rid}  margin={margin:.4f} =====")
        print(f"  read span chr {rst}-{ren}  spanned diff junction(s): "+" ".join(f"({s},{e})" for s,e in sorted(spanned)))
        for rank,(i,v,ne) in enumerate(fin):
            tag="BEST" if rank==0 else f"#{rank+1}"
            loc=local(tc[i].intron_chain.introns,lo,hi)
            print(f"    {tag:4} NLL={v:.4f} nev={ne:3d}  local@[{lo}-{hi}]: "+(" ".join(f"({s},{e})" for s,e in loc) or "(no intron / exonic)"))
    if (n+1)%100==0: sys.stderr.write(f"  {n+1}/{len(intervals)} got={got}\n")
sys.stderr.write(f"DONE got={got}\n"); bam.close()
