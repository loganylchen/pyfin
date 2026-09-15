# The fusion analysis on a whole real library

Date: 2026-09-15. The locus fixture (`FUSION_REAL_DATA_2026-09-15.md`) proved
the regression on 95 reads over two 300-kb slices. This is the acceptance run
the fixture explicitly did not cover: the **entire** K562 direct-RNA library
(SG-NEx replicate4 run1, 1,074,197 alignment records) against the **whole**
GRCh38 primary assembly with GENCODE v44 — nothing told it where to look.

## What it found

69 fused-transcript models from 595 chains over 23,403 candidate molecules.
The classification does the separating:

| class | models |
| --- | ---: |
| `readthrough_or_splicing_candidate` | 48 |
| `intrachromosomal_candidate` | 19 |
| `interchromosomal_candidate` | **2** |

**Four models have three arms** — the multi-gene case — and all four are
classified readthrough rather than presented as fusion calls.

Three of those four have an *unannotated* middle arm (`chr16:66853404-66853685`,
`chr6:27952297-27952456`, `chr4:83455946-83456107`). That is exactly the
`{A} -> {} -> {B}` shape Codex identified as a false negative in the
adjacent-pair nomination rule, and fixing it is what produced them: on the
previous rule this same library yielded 66 models with **one** multi-gene case.
The fourth, `ABRAXAS1:e8 > chr4:… > ENSG00000163312.11`, names a gene on each
outer arm.

Both interchromosomal candidates in the whole library:

| support | junction | arms |
| ---: | --- | --- |
| 4 | **`BCR:e14>ABL1:e2`** | `arm_fsm` / `arm_fsm` |
| 3 | `ZFAND3:e1>chr5:51015098-51015563:?` | `single_exon_compatible` / `unannotated_locus` |

**BCR-ABL1 is the top interchromosomal call by support**, found genome-wide:

```
BCR:e14>ABL1:e2   chr22:23290413:+ > chr9:130854063:+
arm_fsm / arm_fsm    7172 bp    4 molecules    single_molecule_spanned
BCR:[ENST00000305877.13, ENST00000359540.7] > ABL1:[ENST00000318560.6, ENST00000372348.9]
```

That is K562's e14a2 (b3a2) junction, at the fused-transcript level: which
isoform of each gene, the conventional exon-rank junction name, and a sequence.
The `fusion_model_exons.tsv` rows for it are 24 exons across chr22 and chr9
summing to exactly 7172 bases — the declared length and the FASTA record length.

## The precision signal

Six models have **both arms annotated** (`arm_fsm`/`arm_ism_*` on each side).
Every one is either the known driver or a known readthrough/paralog locus:

| support | class | junction |
| ---: | --- | --- |
| 6 | readthrough | `DHRS4:e3>DHRS4L2:e4/7` |
| 5 | readthrough | `GAGE2A:e3>GAGE1:e4` |
| 4 | intrachromosomal | `BAG6:e2>SLC44A4:e15` |
| 4 | **interchromosomal** | **`BCR:e14>ABL1:e2`** |
| 2 | intrachromosomal | `BAG6:e2>SLC44A4:e15` |
| 1 | readthrough | `DHRS4:e5>DHRS4L2:e4/5/7` |

DHRS4/DHRS4L2 and GAGE2A/GAGE1 are tandem paralog clusters that produce
chimeric transcripts; BAG6-SLC44A4 is an MHC-region readthrough. No junk
reaches this tier.

## Two things a downstream user must know

**1. chrM is 12 of the models, including the two highest-support models in the whole
run (390 and 91 molecules).** The mitochondrial genome is circular, so a read
spanning its origin presents as a chr M-to-chr M "fusion". These are classified
intrachromosomal/readthrough, never interchromosomal, but any support-ranked
list puts them first. Filter `chrM` before reading the table.

**2. One event can occupy more than one row.** BAG6-SLC44A4 appears three
times; two of those differ by 6 bp at the 3' breakpoint and list *each other's*
molecules in `nearby_read_ids`. This is the merge rule working as designed:
models merge only when some molecule is an exact member of both — coordinate
proximity is not evidence, because a corrector that decisively settled on a
2-bp-shifted placement is making a claim, not a rounding error. The
cross-listing is what makes the grouping recoverable: join rows on
`nearby_read_ids` to collapse them yourself.

## Runtime, and one thing I could not explain

The run takes about 1h50: ~5 minutes to scan 1,074,197 records and build the
whole-genome mappy index, ~3 minutes of forked realignment, then ~6,277 s in
`analyze_fusion_input`, then ~350 s for the transcript layer.

What is NOT the bottleneck, measured rather than assumed: chain construction.
Its `O(clusters x eligible)` subpath sweep -- which I suspected first -- runs in
**1 second** on 1,874 eligible paths and 595 clusters. So the cost is the
per-molecule adjudication, where the cheap exits are "no eligible path" and
"one arm" and the expensive path is a stitched-template alignment per candidate.

I tried to attribute that cost per molecule and the probe did not survive
scrutiny: its elapsed values tracked wall-clock rather than per-molecule time,
and it emitted more records than the loop had iterations left. Two independent
inconsistencies mean it was measuring something other than what it claimed, so
its numbers are not reported and it was not kept. A batch counter in the same
loop reports passing 22,000 of 23,403 in 8 seconds, which cannot be reconciled
with a 6,277 s loop either. **One of those two probes is wrong and I have not
determined which**; the phase boundaries above come from the pipeline's own
per-phase logs, which are independent of both. Profiling this properly is the
next operational task -- it is not a correctness question.

Node contention is real and large: the early phases ran 33x faster on an idle
node than while a full pytest suite and several Codex reviews shared it.

## The failure that only a named diagnostic could have found

Two whole-library runs died at the very end -- after every artifact was
computed -- on `fusion input changed before publication`. The message named
nothing, so the first death was undiagnosable and I lost the run.

With the inputs and differing fields in the message, the second death answered
itself: `st_dev 144->138` on the `/autofs` fastq, with inode, size, mtime and
ctime all unchanged. The automounter had expired the share and remounted it
during the two-hour run. `st_dev` names the mount, not the file, and is now out
of `_snapshot`; every remaining field speaks to content. The work directory is
preserved on publication failure, which is why the 69 models above survived the
crash that produced them.

## Nomination, and why the library is tractable at all

`min_clip_bases` 250 and the gene-transition rule bring 1,074,197 records to
~23,400 candidate molecules. Both bounds were measured on this library rather
than assumed:

- at `min_clip_bases=100` the direct-RNA adapter/poly(A) tail alone nominates
  8.6% of primaries and trips the 100k candidate cap;
- the untightened gene-transition rule nominated 6.09% of chr9+chr22 records,
  because `labels()` attaches a gene on one base of exonic overlap. Requiring a
  genuine transition — each side carrying a gene the other lacks, across any
  pair of arms rather than adjacent ones — takes that to 0.19%.

## What this still does not claim

`supported_rna_fusion` is RNA-level chain evidence, not a demonstrated genomic
rearrangement. The adapter-chimera negative control is now in the suite
(`tests/unit/test_fusion_adapter_chimera.py`) and passes, but it is synthetic:
it shows the rules exclude an adapter between two arms, not that no real
artefact class survives them. Nothing here is a benchmark against a fusion
caller — that comparison has not been run.
