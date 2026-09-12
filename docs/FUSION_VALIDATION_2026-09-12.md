# Independent fusion implementation: validation and scope

Historical sequence-stage report. The signal/final-assembly extension and its
additional validation are documented in `FUSION_SIGNAL_VALIDATION.md`.

Date: 2026-09-12. Baseline: `945fdae04a5a57750b47f60655ddeec34705cf5a`.
The changes are uncommitted. This report concerns implementation correctness
and assembly isolation, not comparative fusion accuracy on an independent cohort.

## Implemented behavior

The public `--fusion` command runs ordinary assembly with
`PipelineConfig.fusion_enabled=False` and no assembly BEDPE writer. After
assembly cleanup, it calls the independent `fin.fusion.pipeline.run_fusion_analysis`.
Only CLI orchestration changes in the pre-existing production Python files;
ordinary discovery, assignment, scoring, selection, refit and runner source are
unchanged. The standalone Python entry point never constructs an assembly runner.

Independent fusion analysis now connects:

- complete index-free BAM scanning, SA/supplementary deduplication, forward-read
  coordinates and optional full-sequence recovery from FASTQ/FASTA;
- full-molecule, existing-arm and clipped-tail genome searches, including
  normal and paralogous alternatives and explicit search limits;
- annotation-derived gene labels and observed N-arm complete-read paths;
- read/motif/annotation nominations for internal introns and fusion edges,
  fixed-query affine contrasts, joint coordinate edits and full stitched-
  transcript sequence validation;
- complete-molecule support, separately counted partial-edge support, bounded
  cluster diameter and explicit coordinate alternatives;
- deterministic chain tables, N-1 BEDPE edges, alternative FASTA templates and
  arm/exon paths, per-read and junction audits, provenance and atomic publication.

Local corrections that yield no publishable proposal are still written to the
junction audit. Canonical motifs or annotation do not supply expression votes.
Low-MAPQ competing mappings are not silently filtered before adjudication.
One read never supplies two independent full-chain votes.

## Requirement-to-evidence audit

| Requirement | Evidence inspected |
| --- | --- |
| Research precedes integration | `FUSION_IMPLEMENTATION_CONTRACT_2026-09-12.md`: targeted primary-source retrieval and design consequences, recorded before CLI integration |
| Assembly decisions and abundance unchanged | Three-way byte comparisons below, plus production diff restricted to `fin/cli.py` orchestration |
| Fusion-only locus, no assembly interval dependency | `test_fusion_input.py::test_fusion_only_locus_is_scanned_without_index_annotation_or_intervals`; four-gene BAM-to-output integration without ordinary read coverage |
| Missing annotation does not invent gene identities | `test_annotation_free_file_pipeline_reports_loci_without_inventing_gene_names`: three loci recovered with empty gene labels and an explicit manifest status |
| Primary/supplementary/strand/hard-clip invariance | Record/path tests and `test_reverse_middle_primary_and_hard_clip_recovery_preserve_complete_chain`: first primary, reverse middle primary and hard-only representations give the same chain table |
| Three/four genes and complete-read support | Real mappy three/four-arm tests, sample-level three-gene tests, four-gene file-to-output test checking every BEDPE edge |
| No A-B + B-C transitive invention | `test_independent_pairs_never_manufacture_a_three_gene_fusion` and complete/partial support unit cases |
| Internal and inter-gene junction correction | `test_shifted_fusion_junction_is_corrected_in_file_to_output_path`; `test_internal_splice_and_fusion_edges_are_jointly_corrected_end_to_end` |
| Normal, paralogous, ambiguous and high-gap negatives | Actual mappy normal/paralog competition tests; whole-template wrong-middle-arm/adapter tests; input/path ambiguity and gap tests |
| Noncanonical retention and microhomology abstention | `test_fusion_junction_decisions.py`: annotation cannot override exact read sequence, equal-template and canonical microhomology ties remain unresolved |
| Complete rather than separately aligned edge evidence | `test_fusion_transcript_validation.py`: one forward alignment must span every stitched boundary; a partial read cannot validate a full chain |
| Determinism and bounded clusters | `test_record_order_does_not_change_published_biological_artifacts`; duplicate-terminal-extent regression; complete-link diameter bridge test |
| Limits, failures and missing evidence visible | Input/search limit tests; malformed SA/FASTQ tests; `test_failed_local_corrections_remain_in_published_audit` |
| No stale-success publication | Failure-rerun, concurrent-writer, invalid-setting and unmanaged-directory tests in `test_fusion_independent_pipeline.py` |

All test names in the table refer to the current tests under `tests/unit/` or
`tests/integration/`, not to proposed tests. Real mappy alignment is used in the
sample and file-level integration cases. The CLI isolation case runs the actual
assembly runner, not a mock producing sentinel files.

## Assembly byte isolation: synthetic positive fixture

The fixture has two three-arm fusion reads plus four ordinary two-exon reads.
The resulting assembly contains two transcripts and nonempty scoring/refit
outputs. The fixture's small SLOW5 record only supplies a valid input file:
**M2 and endpoint channels are explicitly disabled in this test**. It is not a
signal-scoring validation.

Three independently executed CLI runs:

1. Baseline `fin/` extracted using `git archive` from the commit above.
2. Current source, fusion off.
3. Current source, fusion on, with an independently recovered three-gene chain.

All three runs have the following SHA-256 values:

| Artifact | SHA-256 |
| --- | --- |
| `assembly.gtf` | `b3731bf89bd91d817b9813461cd3a544c2378383c5634cde2c295c843222b3f3` |
| `scores.tsv` | `bcffdbd5c5399898fe73333cc2993acfa908224b24370506d6084ceec980f872` |
| `abundance_refit.json` | `1cfc99f3b6538d5acd44a5de55bd087fef29191f446fe1b5c69b753d4f4b6a0e` |

The baseline hashes and a reproducible input generator are now exercised by
`test_real_cli_assembly_bytes_match_prechange_baseline_with_fusion_off_and_on`.
Initial live run artifacts: `/tmp/pyfin-fusion-cli-blyaq7b0/`.

## Assembly byte isolation: real sequence and signal

Sample: `SGNex_K562_directRNA_replicate4_run1`, Sequin `chrIS`, 575 BAM records.
The BAM was restricted to `chrIS` using indexed fetch; the original input files
were not rewritten. The full-annotation GTF, FASTQ, BLOW5 and Sequin reference
come from the existing Sequin grid. This is one real-sample slice, not a claim
of validation across all chromosomes or the complete benchmark grid.

Execution environment:

- `experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif`;
- Python 3.10, CPU, one worker;
- krill package from `experiments/prod_validation/_img/krill016_cpu`;
- mappy 2.31;
- ordinary `real-drna` defaults: `m2_metric=summed_llr`,
  `assignment_channels=ends`, `folded_tss_recovery_mode=require`.

The logs contain executed M2 `summed_llr` contrasts and nonzero refined-read
counts. This was not merely an M2-off or missing-backend fallback. Baseline,
current-off and current-on were run into fresh output directories and all
three primary outputs match byte-for-byte:

| Artifact | SHA-256 for all three runs |
| --- | --- |
| `assembly.gtf` | `2a393d30f83189d39d772c5d672c042bd83eb47d1e672d910a1143c2082c6a9f` |
| `scores.tsv` | `99d8fe803007ae488f30907a11ca872a9b44a563dc31bca55162350c54fa637b` |
| `abundance_refit.json` | `5e6b5465231195590f0f719eec614404c17f021b81fc99682f0495b378d444a0` |

Current-off/on `run_manifest.json` source hash:
`ec9b024f967c306e8e100a9a005b736cc7b7e2ccd42b08ab6c3abbb9ba8c73bd`.
Baseline source hash:
`691a0ecba47acb30a15f4739c6d10e006dc5e86933d85f3fdf1706593b1fd68d`.
Their `source_root` fields independently identify the current and archived
source trees. Ordinary configuration has fusion disabled in all three runs.

Raw evidence directory:
`/tmp/pyfin-fusion-cli-blyaq7b0/sequin/` contains `inputs.json`, the extracted
`chrIS.bam`, `final_*.command.json`, `final_*.log`, three `final_*` output
directories and `final_comparison.json`. Each command JSON is the exact argv
passed to `singularity exec`; the input paths originate from
`experiments/prod_validation/_sequingrid/smc/SGNex_K562_directRNA_replicate4_run1/full/out/run_manifest.json`.
The independent scan examined 103 nominated molecules; its own fusion manifest
explicitly says that fusion signal scoring was not performed.

## Regression results

Native Python 3.12 environment:
`/tmp/pyfin-review-20260912-env/bin/python`.

Focused fusion and CLI validation:

```text
python -m pytest tests/unit/test_fusion*.py tests/integration/test_fusion_independent_pipeline.py tests/unit/test_cli.py -q
196 passed
```

This includes eleven file/CLI-level integration tests. Log:
`/tmp/pyfin-fusion-focused-20260912.log`.

Full regression:

```text
python -m pytest tests/unit tests/integration -q
2333 passed, 27 skipped, 2 failed
```

The two failures are the existing
`TestTieNllKrillUnavailable.test_gate_on_populates_ties_from_raw` and
`test_gate_off_leaves_ties_empty`. Both require an import of krill that this
native environment does not provide. The exact two tests were also executed
against the separately extracted baseline source and failed at the same
`fin/pipeline/assignment.py:369` import. They were not fixed or relabeled as
passing. Logs: `/tmp/pyfin-fusion-regression-20260912.log` and
`/tmp/pyfin-fusion-baseline-known-failures.log`.

The relevant source compiles. Whitespace checks pass on this change's tracked
and new files; the unrelated pre-existing `HANDOFF.md` EOF whitespace warning
was left untouched. No files were staged and no commit was created.

## Interpretation and remaining scientific limits

- The implementation uses **sequence-only fusion correction**. Assembly retains
  its normal krill signal scoring, but fusion records state `not_scored` rather
  than claiming signal support.
- `supported_rna_fusion` denotes validated RNA-level chain evidence with resolved
  sequence coordinates, not a demonstrated genomic rearrangement or driver.
  `supported_chain_unresolved_coordinates` preserves gene/locus-chain support
  without claiming a unique breakpoint. Insufficient-support and readthrough
  rows are not accepted fusion calls.
- Search completion is completion of a bounded heuristic mapping procedure, not
  proof that every possible genomic explanation was enumerated. Limits and
  competing mappings can lead to abstention.
- No external fusion accuracy, sensitivity, clinical validity or runtime ranking
  is inferred from these tests or from assembly's SIRV/Sequin F1. Such claims
  require independent fusion-positive and fusion-negative biological validation.
- The legacy interval-coupled helpers remain unchanged for compatibility. The
  isolated interfaces are the public CLI path and `run_fusion_analysis`; directly
  enabling the old runner's `PipelineConfig.fusion_enabled` is not this API.
