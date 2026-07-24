#!/usr/bin/env python3
"""Compare 3 ways to combine per-event NLL for resolvable tie reads:
 (a) MEAN over the wide +-10bp junction window   (current pyfin M2)
 (b) SUM  over the wide window                    (no divide)
 (c) SUM  over a TIGHT window (+-FLANK bp around the exact differing junction coords only)
Reports the best-vs-2nd margin distribution for each. Tests whether summed-LLR over the
differing k-mers turns the ~0.08 mean-margin into a decisive (~2-6) margin.
"""
import os, sys, math
import numpy as np
import mappy, pysam
REPO="/SSD/logan/dev/pyfin"; sys.path.insert(0,REPO)
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_gtf import GTFReader
from fin.candidates.discovery import discover_candidates
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit
from fin.scoring.m2_junction_nll import class_junction_window_set, _zrecords, tx2genome_array, _tx2genome
from fin.scoring.krill_aligner import make_krill_aligner, krill_thread_count
import krill
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"; S="SGNex_H9_directRNA_replicate2_run2"
BAM=f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
GENOME=f"{B}/refs/host/GRCh38.primary_assembly.genome.fa"
GTF=f"{B}/results/gencode_full_sweep/_ref/p00/annotation.gtf"
SIGNAL="/autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/blow5/nanopore.blow5"
LIMIT=int(os.environ.get("LIMIT","1200")); USE_GPU=os.environ.get("USE_GPU","1")=="1"
FLANK=int(os.environ.get("FLANK","6"))
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
sys.stderr.write(f"intervals:{len(intervals)} gpu={USE_GPU} flank={FLANK}\n")
preset=get_m1_preset(); nthr=krill_thread_count()
ka,eff=make_krill_aligner(krill,"rna002",USE_GPU,hmm_confidence=False,num_thread=nthr)
def spanmap(chrom,s,e):
    m={}
    for r in bam.fetch(chrom,max(0,s-10),e+10):
        if r.is_supplementary or r.is_secondary or r.is_unmapped: continue
        m.setdefault(r.query_name,(r.reference_start,r.reference_end))
    return m

def per_event(rid,seq,cand,aln):
    """Return list of (genomic_pos, nll) for the read's eventalign vs this candidate."""
    bh=None; ba=None
    for h in aln.map(seq):
        v=score_hit(h)
        if v is not None and (ba is None or v>ba): ba,bh=v,h
    if bh is None: return []
    try:
        recs=krill.align_read_variants(SIGNAL,rid,{cand.candidate_id:cand.sequence[bh.r_st:bh.r_en]},
              pore="rna002",use_gpu=eff,num_thread=nthr,aligner=ka,start=bh.r_st)
    except Exception: return []
    res={x.get("variant_label"):x for x in recs}.get(cand.candidate_id)
    if res is None or res.get("status",-1)!=0: return []
    pos=res["position"]; z,sd=_zrecords(res)
    if len(pos)==0: return []
    gen=tx2genome_array(cand,np.asarray(pos,dtype=np.int64))
    out=[]
    for i in range(len(pos)):
        if not np.isfinite(z[i]): continue
        g=gen[i]
        if g<0: continue
        out.append((int(g),0.5*z[i]*z[i]+math.log(sd[i])))
    return out

def score(ev,gposset):
    s=0.0; n=0
    for g,nll in ev:
        if g in gposset: s+=nll; n+=1
    return (s,n)

mm={"mean":[], "sum_wide":[], "sum_tight":[]}
n_res=0; n_used=0
for n,iv in enumerate(intervals):
    cs=cseq(iv.chrom)
    cset=discover_candidates(iv,BAM,gtf,cs,threshold=24,min_novel_reads=1,chain_cluster=True,canonical_search_bp=0)
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
        tc=[cl[j] for j in tie]
        jsets=[set(c.intron_chain.introns) for c in tc]
        disc=set().union(*jsets)-set(jsets[0]).intersection(*jsets)
        rst,ren=spans[rid]
        spanned=[(s,e) for (s,e) in disc if rst<=s and e<=ren]
        if not spanned: continue
        n_res+=1
        wide=class_junction_window_set(tc,flank=2,k=10)
        if not wide: continue
        tight=set()
        for (s,e) in spanned:
            for c0 in (s,e):
                for g in range(c0-FLANK,c0+FLANK+1): tight.add(g)
        ta=[aln[j] for j in tie]
        rows={"mean":[], "sum_wide":[], "sum_tight":[]}
        ok=True
        for idx,c in enumerate(tc):
            ev=per_event(rid,seq,c,ta[idx])
            if not ev: ok=False; break
            sw,nw=score(ev,wide); st,nt=score(ev,tight)
            if nw<1: ok=False; break
            rows["mean"].append(sw/nw)
            rows["sum_wide"].append(sw)
            rows["sum_tight"].append(st)
        if not ok or len(rows["mean"])<2: continue
        n_used+=1
        for k in mm:
            v=sorted(rows[k]); mm[k].append(v[1]-v[0])
    if (n+1)%200==0 and mm["mean"]:
        a=np.array(mm["mean"]); b=np.array(mm["sum_wide"]); c=np.array(mm["sum_tight"])
        sys.stderr.write(f"  [{n+1}] used={n_used} | mean med={np.median(a):.3f} | "
          f"sum_wide med={np.median(b):.2f} | sum_tight med={np.median(c):.2f} "
          f"(tight>=2:{100*np.mean(c>=2):.0f}%)\n"); sys.stderr.flush()

print(f"\n=== combine-method comparison on {n_used} resolvable reads (FLANK={FLANK}) ===")
for k in ("mean","sum_wide","sum_tight"):
    a=np.array(mm[k]) if mm[k] else np.array([0.0])
    print(f"{k:10} median={np.median(a):.3f} mean={np.mean(a):.3f} p75={np.percentile(a,75):.3f} "
          f"max={a.max():.2f} | >=0.5:{100*np.mean(a>=0.5):.0f}% >=2:{100*np.mean(a>=2):.0f}% >=5:{100*np.mean(a>=5):.0f}%")
bam.close()
