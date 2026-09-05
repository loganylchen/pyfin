import sys, gzip, numpy as np, pysam, pyslow5, krill, importlib.metadata as m
from fin.io.io_gtf import GTFReader
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/sirv4_full_sweep_ratio_grid_v2/SGNex_H9_directRNA_replicate4_run1"
bam=f"{B}/align/SGNex_H9_directRNA_replicate4_run1.sorted.bam"; blow5=f"{B}/subset/mapped.blow5"; fq=f"{B}/subset/mapped.fq.gz"
gtf=GTFReader("experiments/prod_validation/_sirv4grid/_ref/SIRV_Set4.gtf"); gtf.open(); gtf.parse()
genome=pysam.FastaFile("experiments/prod_validation/_sirv4grid/_ref/SIRV_Set4.fa")
# pick a spliced read on SIRV3 - strand near the 8966/8975 alternative donors
rid=None
for rd in pysam.AlignmentFile(bam).fetch("SIRV3", 8900, 9300):
    if rd.is_secondary or rd.is_supplementary or rd.mapping_quality<1 or not rd.is_reverse: continue
    if any(op==3 for op,_ in rd.cigartuples) and rd.reference_start<8950 and rd.reference_end>9200:
        rid=rd.query_name; break
print("read", rid)
seq=None
with pysam.FastxFile(fq) as fh:
    for r in fh:
        if r.name==rid: seq=r.sequence; break
print("read len", len(seq))
# candidate = SIRV310 (donor 8966) and SIRV309 (donor 8975) spliced sequences
txs={t.transcript_id:t for t in gtf.get_transcripts_in_region("SIRV3",8000,10000)}
gfa={"SIRV3": genome.fetch("SIRV3")}
c310=txs["SIRV310"].get_spliced_sequence(gfa["SIRV3"]); c309=txs["SIRV309"].get_spliced_sequence(gfa["SIRV3"])
print("SIRV310 len", len(c310), "SIRV309 len", len(c309))
s5=pyslow5.Open(blow5,"r"); rec=s5.get_read(rid, pA=False); s5.close()
common=dict(signal=rec["signal"].astype("float32"), digitisation=float(rec["digitisation"]), offset=float(rec["offset"]), range=float(rec["range"]), sample_rate=float(rec["sampling_rate"]))
al=krill.Aligner(pore="rna002", use_gpu=krill.built_with_cuda, batch_size=4)
full=al.align([dict(read_id=rid, sequence=c310, start=0, **common), dict(read_id=rid, sequence=c309, start=0, **common)])
for r,name in zip(full,("SIRV310","SIRV309")):
    print(name, "status", r["status"], "n_events", len(r["position"]), "init_scale/shift", r.get("init_scale"), r.get("init_shift"), "post scale/shift", r.get("scale"), r.get("shift"), "keys+", [k for k in r if k.startswith("init")])
# slice run: junction region in transcript coords of SIRV310, reuse init scaling
r=full[0]; pos=np.asarray(r["position"]); 
# transcript coordinate of donor 8966 on '-' strand: find events whose tx pos maps near it -> use middle slice around median position for smoke
K=5; center=(K-1)//2
plo,phi=np.percentile(pos,25),np.percentile(pos,75); mask=(pos>=plo)&(pos<=phi)
st=np.asarray(r["start_idx"])[mask]; en=np.asarray(r["end_idx"])[mask]
lo=int(min(st.min(),en.min())); hi=int(max(st.max(),en.max()))
k0=int(pos[mask].min())-center; k1=int(pos[mask].max())-center
print("slice: kmers", k0, k1, "seq len", len(c310[k0:k1+K]), "samples", hi-lo, "n_full_events_in_region", int(mask.sum()))
for gpu in (True, False):
    al2=krill.Aligner(pore="rna002", use_gpu=gpu, batch_size=2)
    part=al2.align([dict(read_id=rid, sequence=c310[k0:k1+K], signal=common["signal"][lo:hi], start=k0, init_scale=float(r["init_scale"]), init_shift=float(r["init_shift"]), **{k:v for k,v in common.items() if k!="signal"})])[0]
    print("gpu" if gpu else "cpu", "slice status", part["status"], "n_events", len(part["position"]), "echo init", part.get("init_scale"), part.get("init_shift"))
    if part["status"]==0: break
# interior agreement: event position->kmer identical between full and slice for shared positions (drop 3 edge events)
fp={int(p):k for p,k in zip(r["position"], r["model_kmer"])}
sp=[(int(p),k) for p,k in zip(part["position"], part["model_kmer"])][3:-3]
agree=sum(1 for p,k in sp if fp.get(p)==k); print("interior agreement", agree, "/", len(sp), f"= {agree/max(len(sp),1):.4f}")
print("krill", m.version("krill"), "cuda", krill.built_with_cuda)
