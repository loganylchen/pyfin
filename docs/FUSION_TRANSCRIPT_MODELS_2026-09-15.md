# Fused-transcript models

Date: 2026-09-15. Landed as `25e5448` on top of the independent fusion
analysis (`docs/FUSION_IMPLEMENTATION_CONTRACT_2026-09-12.md`,
`docs/FUSION_VALIDATION_2026-09-12.md`, `docs/FUSION_SIGNAL_VALIDATION.md`).
Design consulted with and reviewed by Codex (gpt-6-astra) in three rounds.

## Why a second layer

The chain layer answers "which ordered gene path and breakpoint tuple do
complete molecules support?". It already separates different internal splice
patterns at one breakpoint (its keys include every internal donor/acceptor),
but it does not say which annotated transcript of each gene took part, it
cannot extend a 5'-degraded direct-RNA read to a transcript start, its
identity depends on observed coordinate modes, and its FASTA records are
opaque hashes. Two molecules that lost different amounts of the first exon to
degradation could even split into two chains. A downstream user asking "is
this e13a2 or e14a2, and what is the full transcript sequence?" has no answer
at that level.

The transcript layer reads the finished `FusionAnalysis` and never feeds back
into chain support, correction or signal. The six chain-level artifacts stay
byte-identical; the record-order and assembly-isolation tests are unchanged.

## What a model is

`FusedTranscriptModel`: ordered arms, each `(chrom, strand, exons ASCENDING,
resolution)`. Identity is the ordered exon structure -- never support counts,
transcript ids or chain ids. `FusionArm.blocks` are in transcript order
(descending on minus); `genomic_exons()` is the single conversion point and
everything downstream works on ascending coordinates.

### Per-arm resolution (`resolve_arm`)

Against a transcript-preserving index (the gene index deliberately merges
isoforms and cannot resolve them). A transcript is compatible with an arm iff:

1. every observed intron is an intron of the transcript;
2. no annotated intron overlaps an observed exon (the molecule aligned
   continuously across a splice site of that isoform -- incompatible);
3. the breakpoint-facing boundary is exonic in the transcript (the retained
   slice from TSS/entry to exit/TES reaches it exactly);
4. every observed exon lies inside the retained slice.

Missing annotated introns are tolerated only outside the observed coverage.
Classes: `arm_fsm`, `arm_ism_5p`, `arm_ism_3p`, `arm_ism_both`,
`single_exon_compatible`; without a compatible transcript
`novel_junction`, `novel_combination_of_known_junctions`,
`single_exon_unresolved`, `unannotated_locus`. **All** compatible transcripts
are reported; the class is the tightest, for display, and is never used as
evidence against a longer compatible isoform. A completion is offered only
when every compatible transcript retains the same exons.

Exon ranks of the exit-facing and entry-receiving exons give the conventional
junction name: `BCR:e13>ABL1:e2`.

### Endpoint completeness is not splice class

A molecule that lost 30 bp of exon 1 but saw every intron is `arm_fsm` with
`five_prime_status = annotated_extension_available`. Statuses:
`annotated_boundary_observed` (the model's outer exon equals the annotated
one), `annotated_extension_available`, `observed_boundary_5p_completeness_
unknown` / `observed_end_completeness_unknown` (no agreeing completion).
`five_prime_observed_read_ids` / `three_prime_observed_read_ids` record which
molecules actually reach the model's outer extents.

## How molecules are assigned

* **Structures** = splice keys (everything observed except the outer 5'/3'
  ends), collected from every chain's full-support molecules. No filter on
  chain status is applied -- readthrough, circular-origin and
  unresolved-coordinate chains all enter model construction, and several
  individually insufficient chains can jointly clear the model's own
  distinct-molecule threshold. (This paragraph previously claimed only
  *supported* chains may become models. That never described the code;
  `build_fused_transcript_models` reads every chain.)
* **Envelope**: a model's outer exons are the maximal observed extent among
  molecules of that exact structure (min start / max end) -- an aggregate
  observed envelope, deterministic, and every base of it was observed by some
  molecule. Resolution runs on the final envelope, so a completion always
  contains the observed exons (the writer refuses a negative extension).
* **Collapse**: a structure is collapsed when *every one* of its observations
  fits inside some eligible longer envelope -- a 5' truncation of the first
  arm judged on that molecule's own blocks; a molecule may hold several
  selected paths under one splice key and each is judged separately. Maximality
  is decided among eligible structures only, so an unsupported longer structure
  never erases a supported shorter one. Collapsed structures are written to
  `fusion_model_collapsed.jsonl` with the models they could be: direct RNA
  cannot distinguish a degraded observation from a genuine downstream-start
  isoform, and the shorter hypothesis must stay visible.
* **Votes**: a molecule is UNIQUE to a reported model when it fits that model
  alone among reported models ("unique among reported models", not "observed
  the whole model"); otherwise COMPATIBLE on each it fits. One molecule, at most
  one unique vote. Partial molecules (`partial_read_ids`, from the chain's
  partial support) never acquire full-model support. A-B plus B-C never makes
  A-B-C; a conflicting middle-arm splicing is a different structure, never
  majority-merged.
* **Microhomology**: same-sequence placements of *one* junction -- identical
  structure away from the junction-facing boundaries, breakpoints within
  20 bp, identical projected sequence -- merge into one record
  (`coordinate_status = unresolved_microhomology_same_sequence`) that lists the
  other placements in `coordinate_alternatives`; the record's placement is the
  one annotation agrees with on the most arms. Unrelated events, and same-gene
  structures with different internal exons, are never merged on sequence.

## Outputs (all beside the chain artifacts)

| file | content |
| --- | --- |
| `fusion_models.tsv` | one row per model: genes, junction names, breakpoints, compatible transcripts, classes, support_unique / compatible / partial, read ids, endpoint statuses, lengths, completion id |
| `fusion_models.fa` | evidence-constrained sequence, reference-projected (read substitutions/insertions are NOT reconstructed); ids `FTM_<digest>.obs` |
| `fusion_models.annotated.fa` | annotated completion of the outer arms; ids `FTM_<digest>.ann`; header carries `extension_bases_5p/3p` and `extension_observed=false` |
| `fusion_model_exons.tsv` | multi-chromosome exon table: `model_id, sequence_id, record_kind, arm_index, exon_index_in_arm, exon_index_in_model, chrom, strand, start0, end0, tx_start0, tx_end0, gene_ids, compatible_transcript_ids, arm_match_class, provenance` (observed / annotated_extension / observed_partial_annotated_extension) |
| `fusion_model_evidence.jsonl` | per model: read ids by class, arms with resolution details |
| `fusion_model_collapsed.jsonl` | folded structures with their reads and possible models |

Not a cross-chromosome GTF: the contract forbids a fictitious single-
chromosome transcript, and a GTF whose exons change `seqname` within one
`transcript_id` breaks most parsers. The exon table reconstructs either FASTA
exactly (tested).

FASTA header schema:

```
>FTM_<digest>.obs genes=gA>gB junction=gA:e3>gB:e1 tx=gA:[A1]>gB:[B1] bp=A:600:+>B:1000:+ arms=2 model=FTM_<digest> parent_chain=fusion_chain_... support_unique=2 support_compatible=1 support_partial=0 end5=annotated_extension_available end3=annotated_boundary_observed coordinates=sequence_resolved coordinate_alternatives=0 class=intergenic_fusion_candidate kind=evidence_constrained length=570 sequence_source=reference_projected
```

Values are tokenised (spaces -> `_`, `=` -> `:`, `;` -> `,`); full lists live
in the tables.

## Acceptance evidence

`tests/integration/test_fusion_transcript_models.py` (12 cases through the
real pipeline with real mappy; the fixture was agreed with Codex before the
layer was written; truth sequences are built from exon slices, never with the
implementation's stitcher):

* two isoform combinations at one breakpoint -> two models (A1-B1 570 bp
  observed / 600 bp completion with `extension_bases_5p=30`; A2-B1 500 bp),
  unique support 2/2; the A-550 molecule compatible with both, unique to
  neither, audited as collapsed;
* A-B-C 700 bp with unique 2; the partial B-C molecule contributes no full-
  model vote; middle-arm exons are exactly the retained ones; independent
  A-B and B-C molecules never yield A-B-C; a conflicting middle-arm molecule
  is not majority-merged;
* an accidental 3-bp microhomology at the B->C edge exercised the merge: one
  record, two alternatives, unique 2;
* minus-strand last arm (exact sequence, descending exon table with
  contiguous offsets, `M:1700:-`), minus-strand first arm with 40-bp 5'
  degradation (same model, `five_prime_observed_reads=2`, completion 40 bp);
* both encounter orders of a fitting and a misfitting short observation
  sharing a splice key; supported-shorter / unsupported-longer retention;
  identical-sequence distinct loci; writer failure -> `status="failed"`.

`tests/unit/test_fusion_transcript_models_unit.py` (7 hand-built cases): one
molecule with two paths under one splice key judged per path; same-gene
structures with different internal exons not merged on sequence; distinct
loci with identical sequence; collapsed-audit ids surviving a merge; the
crossing-TSS envelope in both orders (no negative extension); minus-strand
conversion; the unavailable-annotated-sequence branch raising.

## What this does not claim

Sequences are reference projections of molecule alignments, not sample
consensus. An annotated completion is a hypothesis carried by the annotation,
not an observation, and is emitted as a separate product for that reason.
Support is RNA-level chain evidence; nothing here establishes a genomic
rearrangement. Real-sample performance requires independent fusion-positive
and fusion-negative validation; the K562 BCR-ABL1 run is the first positive
control and is recorded separately once complete.
