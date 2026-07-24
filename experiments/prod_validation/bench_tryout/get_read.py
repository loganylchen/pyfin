import pysam
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"; S="SGNex_H9_directRNA_replicate2_run2"
BAM=f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
bam=pysam.AlignmentFile(BAM,"rb")
targets={"c37c1e06-7a20-4582-924d-2ed5e920dd48":("chr1",7961000,8007000),
         "961e2771-c34d-4766-9dfd-0c53eb08f33f":("chr1",6806000,6873000)}
for rid,(c,s,e) in targets.items():
    for r in bam.fetch(c,s,e):
        if r.query_name!=rid or r.is_supplementary or r.is_secondary: continue
        pos=r.reference_start; exons=[]; introns=[]; cur=[pos]
        for op,ln in r.cigartuples:
            if op in (0,2,7,8): pos+=ln
            elif op==3:
                exons.append((cur[0],pos)); introns.append((pos,pos+ln)); pos+=ln; cur=[pos]
        exons.append((cur[0],pos))
        c5=r.cigartuples[0]; c3=r.cigartuples[-1]
        print(f"\nread {rid}")
        print(f"  span {r.reference_start}-{r.reference_end} mapq={r.mapping_quality} softclip5={c5} softclip3={c3}")
        print(f"  EXON blocks: "+" ".join(f"[{a}-{b}]" for a,b in exons))
        print(f"  INTRONS(N): "+(" ".join(f"({a},{b})" for a,b in introns) or "NONE"))
bam.close()
