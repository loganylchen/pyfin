import pysam
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S="SGNex_H9_directRNA_replicate2_run2"
BAM=f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
RID="f9f6a19e-8083-4dc6-ba69-d5830b563d12"
bam=pysam.AlignmentFile(BAM,"rb")
for r in bam.fetch("chr1",1020000,1056000):
    if r.query_name!=RID: continue
    if r.is_supplementary or r.is_secondary: 
        tag="SUPP" if r.is_supplementary else "SEC"
    else: tag="PRIMARY"
    # genomic introns from CIGAR N
    introns=[]; pos=r.reference_start; 
    for op,ln in r.cigartuples:
        if op in (0,2,7,8): pos+=ln
        elif op==3: introns.append((pos,pos+ln)); pos+=ln
    print(f"[{tag}] {r.reference_start}-{r.reference_end}  mapq={r.mapping_quality}  n_introns={len(introns)}  clip5={r.cigartuples[0]} clip3={r.cigartuples[-1]}")
    print("  read introns:", " ".join(f"({s},{e})" for s,e in introns))
bam.close()
