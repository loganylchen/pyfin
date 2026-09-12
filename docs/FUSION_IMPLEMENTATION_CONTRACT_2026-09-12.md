# Fusion completion: research and implementation contract

Extension: `FUSION_SIGNAL_VALIDATION.md` describes the subsequently implemented
signal and completed-assembly junction evidence; earlier sequence-only limits
below describe the original implementation stage.

Date: 2026-09-12. Status: independent implementation connected and acceptance
evidence recorded in `FUSION_VALIDATION_2026-09-12.md`; no independent fusion
accuracy or clinical-validation claim.
Baseline: `945fdae04a5a57750b47f60655ddeec34705cf5a` plus the existing worktree.

## Non-negotiable user requirements

1. Complete fusion discovery, not just the existing pairwise prototype.
2. Do not change ordinary assembly decisions, abundance, or outputs.
3. Correct junctions using read evidence; do not equate canonical sequence or annotation with expression.
4. Represent and adjudicate fusions involving more than two genes.
5. Research the design before changing the production entry point.

The implementation must keep unsupported/ambiguous cases explicit. Resolving a
fusion means determining what the evidence supports, not forcing a unique call.
An RNA chimera is not by itself proof of a genomic rearrangement or a cancer driver.

## Targeted literature retrieval

This is a targeted methods lookup, not an exhaustive systematic review.
Europe PMC search, accessed 2026-09-12:

```text
endpoint: https://www.ebi.ac.uk/europepmc/webservices/rest/search
query: TITLE:JAFFAL OR TITLE:Genion OR TITLE:FusionSeeker OR TITLE:GFSeeker
       OR TITLE:"Accurate fusion transcript identification"
       OR TITLE:"A comprehensive benchmark of transcriptome-wide fusion detection"
format=json, pageSize=30, resultType=core
hitCount=9, retrieved=9, one page
```

Raw responses and parsed full text are in `/tmp/pyfin-fusion-research-20260912/`.
JATS was parsed using the paper-lookup skill's body-validating parser. Full text
was retrieved for JAFFAL, Genion, CTAT-LR-Fusion and GFSeeker. Europe PMC returned
404 for FusionSeeker full text; its metadata is not treated as full-text evidence.
The August 2026 benchmark record is a preprint, not a peer-reviewed result.

### Primary sources and design consequences

| Source | Verified identifier | Consequence for this implementation |
| --- | --- | --- |
| JAFFAL, Genome Biology, 2022 | DOI 10.1186/s13059-021-02588-5; PMC8739696 | Long reads can phase three-gene fusions. Preserve each molecule's ordered path and actual breakpoint-spanning support. |
| Genion, BMC Genomics, 2022 | DOI 10.1186/s12864-022-08339-5; PMC8842519 | Use supplementary alignments and gene annotation, but account for homology and chance chimeric mappings instead of counting split records as independent reads. |
| CTAT-LR-Fusion, Genome Research, 2025 | DOI 10.1101/gr.279200.124; PMC12047241 | Separate nomination from validation against fusion contigs; validate anchors and competing normal/paralogous explanations. Nearby exon boundaries are proposals, not unconditional corrections. |
| GFSeeker, Briefings in Bioinformatics, 2026 | DOI 10.1093/bib/bbaf702; PMC12777712 | Splice paths, consensus and independent reference checks can refine junctions. Do not transplant a new graph aligner into ordinary assembly. |
| DeepChopper, Nature Communications, 2026 | DOI 10.1038/s41467-026-70593-2 | Direct-RNA technical chimeras require an explicit artifact policy. An internal unmapped gap is a warning signal, not a molecular identification of an adapter. |
| SAM specification | SAMv1, samtools/hts-specs | Primary/supplementary status does not define transcript order. Normalize reverse-strand and hard-clipped records into full forward-query coordinates before chaining. |

No published software's accuracy claim is a performance claim for PyFIN. In
particular, the work below requires its own end-to-end and false-positive tests.

## Baseline prototype gaps confirmed before integration

- The runner schedules ordinary read intervals and does not independently
  consume the returned `fusion_reads`; fusion-only loci can never be examined.
- `collect_chimeric_reads` chooses only the longest soft clip and one partner,
  so a middle primary alignment and two flanking partners are not represented.
- `ChimericRead`, `FusionPairCluster` and the old BEDPE writer encode two arms.
- Pairwise clustering can bridge distant breakpoints through intermediate reads.
- Annotated-arm Cartesian products can survive without reads supporting the
  complete combination. This must not be the multi-gene support model.
- Fusion candidates currently compete with ordinary candidates inside assembly
  assignment. That violates the requested isolation when fusion is enabled.

## Architecture decision: an independent fusion analysis

The public `--fusion` path will run normal assembly with the same effective
assembly configuration as without the flag, then run a separate fusion analysis.
Fusion must not insert candidates into assembly EM, withdraw reads from assembly,
modify shared selection thresholds, or write to assembly's evidence/cache shards.

Only CLI/orchestration glue may change outside the fusion package and its tests.
Shared helpers may be called without modifying them: canonical motif enumeration,
sequence orientation/splicing, lazy reference access and krill backend creation.
The legacy internal fusion APIs remain compatibility code during migration; the
new feature is not complete while the public CLI still routes to the old path.

Fusion uses a separate configuration, manifest and output namespace. It scans
the BAM independently, including supplementary records and SA tags, without
requiring an assembly interval or changing the BAM/index. Missing secondary
alignments, annotation or signal must be reported, not invented as negative evidence.

## Ordered, molecule-supported multi-gene model

`AlignmentSegment -> ReadPath -> FusionChain` is the data flow.

- An alignment segment stores full forward-query coordinates, genomic strand,
  ordered exon blocks, quality, and record provenance. SA tags and BAM records
  describing the same alignment are deduplicated before support is counted.
- A read path is an ordered tuple of segments, not a fixed pair. Query-coordinate
  order is authoritative; arbitrary SAM primary choice is not.
- A fusion chain retains all arms and all N-1 fusion junctions. Gene labels are
  derived from exon overlap and strand-aware annotation. Ambiguous overlapping
  genes retain a set of possible labels; absent annotation yields loci, not
  fabricated gene identities. Same-chromosome gene transitions also need handling.
- A-B and B-C observed on different reads do **not** establish A-B-C. A multi-gene
  call requires direct reads spanning its entire ordered chain and every edge.
  Partial reads can add separately reported edge support, never full-chain support.
- Distinct complete A-B-C and A-B-D paths may coexist. A single ambiguous read
  must not count as independent decisive evidence for both paths.
- Candidate/path search bounds are explicit. Exceeding a bound produces a
  search-limit/unresolved audit, not an apparently unique call from a truncated search.

Normal splicing, readthrough candidates, intragenic rearrangements, intergenic
fusion candidates and ambiguous mappings must have distinct reported classes.
Default calls require independent read IDs, not multiple records for one molecule.
Confidence in the gene path and confidence in exact junction placement are separate.

## Junction correction and competing explanations

Two levels need correction: introns inside each arm and the fusion edges between
arms. Ordinary assembly's same-chromosome junction data structures cannot encode
an interchromosomal fusion edge by pretending it is an intron.

1. Retain the observed path and raw coordinates as the baseline hypothesis.
2. Nominate nearby alternatives from distinct-read support, annotation boundaries
   and canonical motifs; motif/annotation alone does not confer support.
3. Keep exons and query anchors valid after each edit. Corrected templates must
   preserve the ordered gene path and must not create an unsupported third arm.
4. Compare alternatives on the same read evidence window, including an unmodified
   baseline. Penalize indels and test the full stitched junction with anchors on
   both sides; do not snap individual coordinates to unrelated medians.
5. Use signal only where a valid competing contrast can be built. Reuse the
   backend, not the falsified unpenalized raw-DTW minimum-SSE rule. Missing signal,
   incomplete event coverage, or incomparable calibration requires abstention.
6. Microhomology or identical alternative templates may leave exact coordinates
   unresolved even when the gene chain is supported; expose alternatives/ranges.
7. Cluster complete breakpoint tuples with a bounded diameter. Transitive
   single-linkage must not merge two distant events through one bridge read.

The implementation must also compare plausible non-fusion/paralog explanations,
flag large internal gaps and low-quality/multimapping arms, and retain sufficient
audit evidence to distinguish technical failure from biological rejection.

## Outputs

- `fusions.tsv`: one complete chain per row, status, class, ordered genes/loci,
  complete-read support, per-edge support and coordinate-resolution status.
- `fusions.bedpe`: one row per adjacent fusion edge, with chain ID and edge index
  encoded in the name; never silently flatten A-B-C into one A-C pair.
- Fusion transcript FASTA plus an exon/arm table or JSON that can represent
  multiple chromosomes. Do not emit a fictitious single-chromosome GTF transcript.
- Per-read/path and junction-correction audits, plus `fusion_manifest.json`.
- Ordinary `assembly.gtf`, `scores.tsv` and `abundance_refit.json` remain unchanged.

New writes must use a run-local lifecycle and atomic publication. A failed rerun
must not publish an old fusion table as the new result.

## Frozen acceptance checklist

### Assembly isolation

- Feature-off results equal the pre-change baseline byte-for-byte on a controlled
  end-to-end fixture: GTF, scoring TSV and refit diagnostics.
- Feature-on ordinary assembly outputs equal feature-off outputs byte-for-byte.
- No ordinary discovery/assignment/selection/scoring implementation changes.
- Existing unit/integration regressions pass subject to explicitly documented
  environment skips; a missing backend is not silently called a passing test.

### Fusion discovery and chain support

- Fusion-only locus is discovered without ordinary-read coverage or GTF.
- Plus/minus orientations, hard/soft clipping, primary reassignment and duplicate
  SA/BAM representations give the same biological ordered path/support.
- A middle primary alignment can yield A-B-C; three and four arms are covered.
- Independent A-B plus B-C reads never generate a supported A-B-C chain.
- Multiple true chains can coexist; ambiguous alternative partner mappings are
  unresolved rather than broken by input order.
- One molecule cannot supply multiple support votes through supplementary records.
- Same-chromosome transitions, readthrough and ordinary splicing are distinguished.

### Junctions, negatives and observability

- Known shifted introns and fusion edges are corrected on positive controls,
  while already-correct junctions remain unchanged on negative controls.
- Read-supported noncanonical edges are not forcibly canonicalized.
- Equal/microhomologous alternatives preserve ambiguity; invalid or unanchored
  edits cannot win. Signal failure is explicit abstention.
- Normal, paralogous, ambiguous, high-gap and malformed-alignment controls are
  not presented as confident fusion calls.
- Limits, missing annotation/backend and partial scans are visible in the audit.
- Chain/event IDs and outputs are deterministic under record/read order changes.

### Evidence required before calling the feature complete

An executable CLI path, the isolation checks, synthetic truth-based end-to-end
positive/negative tests, accurate multi-gene output, and junction correction must
all be present. A literature document, parser, or mocked unit tests alone are
not completion. Real-sample performance claims require additional independent
validation and will not be inferred from assembly's F1 or from these unit tests.

## Earlier implementation checkpoint — 2026-09-12

At this checkpoint the production entry point was not switched yet. The work below
is in new fusion modules and their tests; ordinary assembly source remains
unchanged. This is progress toward the acceptance checklist, not completion.

- `records.py`, `paths.py`, `annotation.py`, `chains.py`: full forward-query
  coordinates, N-arm paths, explicit gene ambiguity, complete-molecule support,
  bounded-diameter clustering and deterministic duplicate handling.
- `input.py`: independent two-pass, index-free BAM scan, including hard-clipped
  primaries without SA; optional FASTQ/FASTA recovery validates full length and
  sequence consistency. Invalid alignments cannot be rescued by replacing their
  sequence. Incomplete or changing input files fail the analysis.
- `realignment.py`: an independent mappy splice index searches the entire
  molecule, existing arm windows and recursively clipped tails. It keeps
  secondary/zero-MAPQ alternatives, including normal full-length and paralogous
  mappings. Hit/window/segment caps return unresolved, with no partial evidence
  presented as a successful search. This remains heuristic genome alignment,
  not proof that every possible genomic explanation has been enumerated.
- `alignment.py`, `junctions.py`: fully consumed fixed-query affine contrasts,
  both-side anchor checks, canonical/read/annotation nominations, explicit
  coordinate ties and unavailable sequence/reference handling. Signal is
  explicitly `not_scored`; no raw-SSE rule or signal-validation claim is added.
- `correction.py`: proposals for every internal splice and inter-arm edge on
  an immutable complete path. All decisions are made on observed geometry;
  combined edits are then rechecked for anchored sequence support. One
  unresolved edge withholds the whole edit transaction. These are per-molecule
  proposals, not sample-level supported fusion calls.

Validation in `/tmp/pyfin-review-20260912-env/bin/python`:

```text
python -m pytest tests/unit/test_fusion*.py -q
151 passed
```

The tests include real mappy searches (not mocked) for three/four-arm synthetic
reads, reverse middle arms, full-length normal alternatives and duplicated
paralog loci; strand-aware fusion/intron correction, four-gene joint correction,
noncanonical retention, microhomology abstention and FASTQ conflicts. A randomly
generated minus/minus positive fixture exposed an equal-score canonical
microhomology alternative; it is retained as a separate ambiguity regression,
not bypassed by weakening the decision margin.

Work still required at that checkpoint (subsequently connected and tested; see
`FUSION_VALIDATION_2026-09-12.md` for current evidence):

1. Connect discovery, reference realignment, path alternatives and correction in
   one sample-level analysis; compare competing explanations without filtering
   away low-MAPQ competitors. Do not count duplicated alternative paths as votes.
2. Determine supported corrected whole-chain coordinates across distinct read
   IDs, keeping full-chain, partial-edge and coordinate-ambiguity evidence
   separate. Recheck complete stitched transcripts, not only local windows.
3. Publish run-local atomic chain/edge/FASTA/arm/audit/manifest outputs and wire
   CLI fusion orchestration outside assembly EM.
4. Execute truth-based end-to-end positive/negative cases and byte comparisons
   for assembly GTF, scoring TSV and refit diagnostics against the pre-change
   baseline and between fusion-off/on. No such end-to-end isolation claim has
   been established by the 151 unit tests.
