One SIRV locus with its raw signal: the end-to-end fixture.

Every alignment record on SIRV1 of SG-NEx H9 direct RNA replicate4 run2 (91
reads, all primary), aligned to the SIRV Set 4 reference with
`minimap2 -ax splice -uf -k14`; the raw signal of exactly those 91 reads
(RNA002, 3012 Hz); the SIRV1 contig; and its eight annotated transcripts.
Nothing was edited -- coordinates, CIGARs, tags, read sequences, signal samples
and the signal header are the originals. The BAM's @PG line, which carried the
aligner command with site-local paths, was dropped.

| file | content |
|---|---|
| `reads.bam`, `.bai` | 91 alignment records on SIRV1 |
| `reads.blow5`, `.idx` | raw signal of the same 91 reads |
| `reference.fa`, `.fai` | the SIRV1 contig (12,643 bp) |
| `annotation.gtf` | 8 SIRV1 transcripts of SIRV Set 4 |

Under the shipped `real-drna` profile the locus forms two intervals, both with
exact alignment-score ties that the signal has to break (17 tied reads, 5
scored contrasts), which is the path a sequence-only fixture cannot reach.

Built by `experiments/prod_validation/build_sirv1_fixture.py`.
