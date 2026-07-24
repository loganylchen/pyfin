#!/usr/bin/env python3
"""Pure-structural M2 skip gate (NO krill). For each M1-tie read, decide WITHOUT
signal whether it can distinguish the tie: does the read's aligned genomic span
actually cross any junction that DIFFERS among the tie candidates? If not, it's a
containment/no-info tie -> M2 would return identical NLL -> SKIP. Reports the split
so we know M2's real workload and how many tie reads carry discriminating info.
"""
import os, sys
from collections import Counter
import mappy, pysam

REPO="/SSD/logan/dev/pyfin"; sys.path.insert(0,REPO)
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_gtf import GTFReader
from fin.candidates.discovery import discover_candidates
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit

B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S="SGNex_H9_directRNA_replicate2_run2"
BAM=f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
GENOME=f"{B}/refs/host/GRCh38.primary_assembly.genome.fa"
GTF=f"{B}/results/gencode_full_sweep/_ref/p00/annotation.gtf"
LIMIT=int(os.environ.get("LIMIT","3000"))

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
sys.stderr.write(f"intervals:{len(intervals)}\n")
preset=get_m1_preset()

def read_span_map(chrom,start,end):
    m={}
    for r in bam.fetch(chrom,max(0,start-10),end+10):
        if r.is_supplementary or r.is_secondary or r.is_unmapped: continue
        m.setdefault(r.query_name,(r.reference_start,r.reference_end))
    return m

n_tie=0; n_resolvable=0; n_skip=0; tie_hist=Counter()
for n,iv in enumerate(intervals):
    cs=cseq(iv.chrom)
    cset=discover_candidates(iv,BAM,gtf,cs,threshold=24,min_novel_reads=1,
                             chain_cluster=True,canonical_search_bp=0)
    cl=list(cset.candidates)
    if len(cl)<2: continue
    rs=getattr(cset,"read_sequences",{}) or {}
    aln=[mappy.Aligner(seq=c.sequence,preset=preset) if c.sequence else None for c in cl]
    spans=read_span_map(iv.chrom,iv.start,iv.end)
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
        n_tie+=1; tie_hist[len(tie)]+=1
        # discriminating junctions = present in some tie cand, not all
        jsets=[set(cl[j].intron_chain.introns) for j in tie]
        allj=set().union(*jsets); commonj=set(jsets[0]).intersection(*jsets)
        disc=allj-commonj
        rst,ren=spans[rid]
        # read resolves tie iff it spans (covers) any discriminating intron
        resolves=any(rst<=s and e<=ren for (s,e) in disc)
        if resolves: n_resolvable+=1
        else: n_skip+=1
    if (n+1)%500==0:
        sys.stderr.write(f"  {n+1}/{len(intervals)} tie={n_tie} M2-worth={n_resolvable} skip={n_skip}\n")

print(f"\n=== M2 skip gate (pure structural, no signal) ===")
print(f"tie reads: {n_tie}")
print(f"  RESOLVABLE (spans a differing junction -> run M2): {n_resolvable} ({100*n_resolvable/n_tie if n_tie else 0:.1f}%)")
print(f"  SKIP (no covered difference -> containment, go structural): {n_skip} ({100*n_skip/n_tie if n_tie else 0:.1f}%)")
print("tie-size hist:", dict(sorted(tie_hist.items())))
bam.close()
