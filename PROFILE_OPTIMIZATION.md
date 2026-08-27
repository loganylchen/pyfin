# PyFIN Profile Optimization and Production Integration

Date: 2026-08-26 to 2026-08-27

Final live-source SHA-256:
`e4b2764b12dd5e966b85ba914421fe3f9e2e0d8c4e40aa90922d0d13988e5256`

This report records the experiments used to turn the former single,
SIRV-oriented CLI operating point into separate SIRV and real-dRNA profiles.
It distinguishes metrics selected for production from experiments that remain
available but are not defaults.

## 1. Objective and acceptance rule

The optimization had three required outputs:

1. A reproducible SIRV operating point.
2. A separate, recall-safe real-dRNA operating point.
3. Production wiring in `fin/`, not experiment-only flags.

There is no single threshold that simultaneously maximizes all-truth recall,
expressed-truth precision, and low-abundance recall. The production objective is
therefore explicit:

- Primary metric: NanoCount expressed-truth transcript F1 at `est_count >= 3`
  (T3), matching the repository's honest-F1 convention.
- Real-dRNA recall guard: also report T1 and reject apparent improvements that
  merely erase low-abundance truth.
- Compatibility metric: standard gffcompare transcript sensitivity, precision,
  and F1 against the full reference.
- Tie-break: when biological metrics are effectively tied, prefer the faster
  configuration.
- Every final run must carry a resolved `run_manifest.json`.

## Validation scope

SIRV began with two tuning replicates and the untouched `replicate3_run2`
holdout. The later guide-routing decision used all six complete H9 SIRV samples
in a paired mean/SUM analysis and separately reports p00 and full-guide cells.

Real-dRNA now has two independent sample checks beyond tuning
`SGNex_H9_directRNA_replicate2_run2` (`gencode_p00`):

- `SGNex_H9_directRNA_replicate3_run1`
- `SGNex_H9_directRNA_replicate4_run2` (blind threshold arbiter)

Their old `subset/mapped.fq.gz` and `mapped.blow5` links are broken, but the
original sample-sheet FASTQ and BLOW5 files remain available. The harness uses
those full raw stores with each sample's own mapped BAM and NanoCount table;
extra non-BAM reads are inert because PyFIN retrieves sequence/signal by BAM read
ID. Thus the samples are independent biological checks, although GENCODE remains
a pseudo-truth rather than absolute isoform truth.

Competitor tools were not re-executed. Their existing p00 GTFs (June/July 2026)
were freshly rescored in this work with the same full GENCODE truth,
sample-specific NanoCount T1/T3 sets, gffcompare image, and `parse_metrics()`.
Historical per-tool runtimes are not hardware-normalized and are not used for
competitive accuracy claims.

## 2. Evaluation data

### SIRV tuning

Two independent direct-RNA replicates, each evaluated with p00 and full guide
annotation:

- `SGNex_H9_directRNA_replicate4_run2`
- `SGNex_H9_directRNA_replicate2_run2`

Truth: `experiments/prod_validation/sirv4/_ref/full/annotation.gtf`.
Expressed truth: each sample's `stage/nanocount.tsv`.

### SIRV holdout

`SGNex_H9_directRNA_replicate3_run2` was not used to select parameters. It was
run once after the tuning matrix on p00 and full guide conditions.

### Real-dRNA optimization

The complete staged GENCODE p00 dataset under
`experiments/prod_validation/gencode/_p00val/stage`:

- 15,190 genomic intervals.
- 21 GB BLOW5 signal.
- Full GENCODE structural truth.
- NanoCount truth at T1 and T3.

The two independent H9 samples use their original raw FASTQ/BLOW5 stores after
verifying the mapped-subset links were the only missing artifacts. Existing
HEYA8/SIRV structural sweeps remain supporting evidence for the already-landed
recheck/containment defaults.

## 3. SIRV results

### Tuning aggregate (2 samples x 2 guide conditions)

| Configuration | Mean honest F1 T3 | Mean corrected recall T3 | Mean honest precision T3 | Mean standard F1 |
| --- | ---: | ---: | ---: | ---: |
| Legacy CLI profile | 83.437 | 84.602 | 83.368 | 56.927 |
| PolyA gate off | 85.035 | 88.965 | 81.703 | 59.068 |
| GTF floor on | 84.549 | 84.260 | 85.594 | 55.643 |
| Final combination: polyA off + GTF floor on + mean M2 | **86.147** | **88.623** | **83.929** | **57.784** |
| Best summed-LLR combination | 86.077 | 88.281 | 84.120 | 57.684 |
| M2 off on final base | 85.229 | 87.956 | 82.777 | 57.270 |

The initial static SUM grid did not beat all-mean M2 on the four-cell tuning
aggregate. A later six-sample paired analysis exposed a consistent guide
interaction using true p00 (the one-transcript stub was omitted):

| Guide condition | Mean M2 F1 T3 | SUM (1,4) F1 T3 | SUM delta | Wins/ties/losses |
| --- | ---: | ---: | ---: | ---: |
| p00 / no usable GTF | 80.370 | **81.091** | **+0.722** | 4 / 1 / 1 |
| Full guide | **92.866** | 92.581 | -0.284 | 2 / 0 / 4 |

The SIRV profile therefore uses `m2_metric=auto`: SUM margin 1/flank 4 when the
GTF is absent or a one-transcript stub, and mean when at least two guide
transcripts exist. A manifest records the concrete metric and
`auto-unguided`/`auto-guided` route.

### Independent holdout

| Configuration | Mean honest F1 T3 | Mean corrected recall T3 | Mean honest precision T3 | Mean standard F1 |
| --- | ---: | ---: | ---: | ---: |
| Legacy profile | 82.892 | 86.842 | 79.906 | **58.505** |
| All-mean optimized profile | 85.513 | **90.132** | 81.366 | 58.034 |
| All-SUM candidate | 85.458 | 89.474 | 81.810 | 57.872 |
| Guide-aware auto profile | **85.778** | **90.132** | **81.861** | 58.144 |

The auto profile improves holdout honest F1 by 2.89 and corrected recall by
3.29 points over legacy. The legacy profile retains a 0.36 all-truth
standard-F1 advantage because it preserves unexpressed GTF truth. Users
benchmarking that alternative objective can explicitly pass
`--no-floor-gtf-abundance`.

## 4. Real-dRNA one-factor results

The first table uses the pre-optimization real profile (`min_abundance=1`,
wide mean M2, soft-mass gate on) as baseline.

| Variant | T1 F1 | T1 recall | T3 F1 | T3 recall | T3 precision | Runtime (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Baseline mean M2 | 32.477 | 34.159 | 30.827 | 45.578 | 23.289 | 1391 |
| Tight summed LLR (2,6) | 32.566 | 34.146 | 30.916 | 45.530 | 23.404 | 611 |
| M2 off | 32.475 | 34.139 | 30.847 | 45.578 | 23.313 | **388** |
| `min_abundance=0` | 31.934 | 34.159 | 30.179 | 45.578 | 22.558 | 1390 |
| `min_abundance=2` | **37.361** | 33.172 | 37.729 | 44.816 | 32.577 | 1378 |
| `min_abundance=3` | 37.096 | 29.803 | **41.179** | 42.745 | **39.724** | 1377 |
| Full-length gate on | 25.906 | 25.171 | 24.415 | 32.627 | 19.506 | 1381 |
| Junction minimum 1 | 27.121 | **40.234** | 21.944 | **49.244** | 14.118 | 1380 |
| Containment cluster off | 32.185 | 34.609 | 30.247 | 45.983 | 22.535 | 1377 |
| Structural recheck off | 31.623 | 34.454 | 29.753 | 45.983 | 21.992 | 1376 |
| Canonical gate off | 32.430 | 34.199 | 30.772 | 45.649 | 23.209 | 1603 |
| Soft/hard-mass gate off | 32.452 | 34.944 | 30.687 | 46.733 | 22.844 | 1334 |
| Isoform-fraction gate off | 32.480 | 34.293 | 30.782 | 45.733 | 23.199 | 1320 |

Interpretation:

- The SIRV full-length gate is decisively invalid for real dRNA.
- `novel_junction_min_reads=2`, containment cluster, and structural recheck are
  real precision levers; turning them off gains recall but loses much more F1.
- Canonical and isoform-fraction gates are near-neutral, with the enabled state
  slightly ahead on F1.
- `max_soft_mass_ratio=2` costs more than one T3 recall point for almost no F1;
  real-dRNA disables it.
- Floors 2 and 3 maximize precision/T3 F1 but impose material T1 recall loss.

### Abundance Pareto refinement

Filtering the baseline final candidate set offline and rescoring with the same
gffcompare image showed the large transition is `abundance == 1` versus
`abundance > 1`, not a special biological meaning of 1.1:

| Floor | T1 F1 | T1 recall | T3 F1 | T3 recall |
| ---: | ---: | ---: | ---: | ---: |
| 1.00 | 32.477 | 34.159 | 30.827 | 45.578 |
| >1.0 (strict boundary) | 36.522 | 33.669 | 36.212 | 45.197 |
| 1.50 | 36.878 | 33.414 | 36.858 | 44.971 |
| 2.00 | 37.361 | 33.172 | 37.729 | 44.816 |
| 3.00 | 37.096 | 29.803 | 41.179 | 42.745 |

The intended setting is a strict boundary (`abundance > 1.0`), not a tuned
`1.01` magnitude. During exploratory output naming, `1.000001` formatted to the
same `t1` cache directory as 1.0, producing stale metrics despite the same
12,575-candidate filtered set as 1.001. That row was discarded. Production now
encodes strict/inclusive semantics explicitly, so EM precision changes cannot
move a biological-looking magic epsilon.

Combining `abundance > 1` with soft-mass gate off recovers the floor's recall
loss while retaining most of its precision gain. Under mean M2 the combination
already dominated the original baseline:

- T1 F1 36.383, recall 34.454.
- T3 F1 35.857, recall 46.352.

## 5. Summed-LLR optimization

Nine real-data combinations were tested:

- Margins: 1, 2, 3.
- Tight-window flanks: 4, 6, 8 bp.

The winner was margin 1, flank 8:

| Metric | Legacy wide mean | Summed LLR (1,8) |
| --- | ---: | ---: |
| T1 F1 | 32.477 | **32.570** |
| T1 recall | **34.159** | 34.152 |
| T3 F1 | 30.827 | **30.929** |
| T3 recall | **45.578** | 45.554 |
| Runtime | 1391 s | **599 s** |

This is a controlled metric comparison: both rows use the one-factor baseline
profile (`min_abundance=1` inclusive and soft-mass gate on); only M2 metric and
its window/margin differ. The intermediate strict-floor profile ran in 503
seconds; the final balanced profile adds mono selection and one global junction
support scan and runs in about 542 seconds. Neither total is attributed solely
to M2.

This confirms the recorded design rationale: on real data, tight sum avoids the
wide mean's dilution and cuts runtime by avoiding wide class windows. It is not
used for containment/cassette contrasts; those abstain because their
candidate-private event populations are asymmetric.

## 6. Final named-profile validation

The final live `--profile real-drna` run used no scientific CLI overrides and
wrote a manifest proving these resolved values:

- `m2_metric=summed_llr`, margin 1, flank 8
- strict novel abundance `>1.0`
- novel mono hard-read floor 5
- finalized junction consensus: tolerance 6, support 2, ratio 2
- soft-mass, full-length, and polyA gates off
- `novel_junction_min_reads=2`
- canonical, containment cluster, and structural recheck enabled
- stable structural novel IDs

Final balanced result:

| Metric | Old real baseline | Final balanced profile | Delta |
| --- | ---: | ---: | ---: |
| T1 honest F1 | 32.477 | **38.691** | +6.214 |
| T1 corrected recall | **34.159** | 33.978 | -0.181 |
| T1 honest precision | 30.953 | **44.922** | +13.969 |
| T3 honest F1 | 30.827 | **39.422** | +8.595 |
| T3 corrected recall | 45.578 | **46.149** | +0.571 |
| T3 honest precision | 23.289 | **34.406** | +11.117 |
| Runtime | 1391 s | **542 s** | 2.57x faster |
| Output transcripts | 16,436 | 11,281 | -5,155 artifacts/low-support models |

The balanced combination preserves the real low-expression guard: T1 recall
changes by -0.18 pp while T1/T3 F1, precision, T3 recall, and runtime improve.
The separate precision profile makes the larger ~1 pp recall trade-off explicit
rather than embedding it in the default:

| Metric | Balanced | Precision |
| --- | ---: | ---: |
| T1 F1 / recall / precision | 38.699 / **33.984** / 44.931 | **38.854** / 32.749 / **47.755** |
| T3 F1 / recall / precision | 39.422 / **46.149** / 34.406 | **40.510** / 44.888 / **36.910** |
| Output transcripts | 11,281 | 10,217 |
| Runtime | 542 s | **369 s** |

Runtime values remain the representative pre-family timings; paired v8
family/overlap runs differed by only 0.13 seconds on tuning and 1.46 seconds on
the larger holdout, so family propagation adds negligible measured overhead.
Two pre-family precision repeats produced byte-identical GTF/TSV files,
confirming that stable structural IDs removed the earlier UUID tie variation.

### Cross-sample real existence and junction improvements

The former final profile was then evaluated against two independent H9 samples.
The balanced profile adds a novel-mono hard-read floor of 5 and finalized-model
junction consensus; the precision profile additionally requires two discovery
reads for every novel candidate. The replicate4_run2 balanced/precision rows
apply the validated exact-key offline transform to its live baseline/mono5_r2
outputs; that transform reproduced the live structure set exactly on tuning and
within 0.007 F1 on replicate3_run1. Tuning and replicate3 rows are live snapping
runs.

| Sample | Previous T3 F1 | Balanced T3 F1 / recall | Precision T3 F1 / recall |
| --- | ---: | ---: | ---: |
| replicate2_run2 (tuning p00) | 36.386 | **39.422 / 46.149** | 40.510 / 44.888 |
| replicate3_run1 | 35.370 | **38.601 / 44.261** | 39.731 / 42.920 |
| replicate4_run2 (blind) | 35.344 | **38.626 / 45.372** | 39.788 / 44.342 |

The mono>=5 gate contributes most of the balanced gain (+2.7 to +3.1 T3 F1)
while losing at most 0.17 T3 recall points. Junction consensus adds +0.15 to
+0.57 F1 and increases recall in every sample. Raising generation support to 2
adds about +1.1 F1 but costs roughly 1 pp T1/T3 recall, so it is exposed as a
separate precision operating point. Raising the mono floor from 5 to 9-11 was
rejected for the balanced default because the blind sample exceeded the
pre-registered 0.6 pp T1 recall-loss guard.

Fresh same-truth rescoring of pre-existing competitor GTFs gives this three-
sample arithmetic mean (accuracy only; competitor tools were not rerun):

| Tool/profile | Mean T1 F1 | Mean T3 F1 | Mean T3 recall | Mean T3 precision |
| --- | ---: | ---: | ---: | ---: |
| IsoQuant | 33.449 | **41.652** | 33.711 | **54.571** |
| PyFIN precision (pre-family holdouts) | **38.491** | 40.020 | 44.050 | 36.676 |
| PyFIN balanced | 38.315 | 38.883 | **45.261** | 34.085 |
| StringTie3 | 38.046 | 38.506 | 44.637 | 33.918 |
| Bambu | 36.139 | 38.483 | 38.706 | 38.460 |

PyFIN balanced therefore moves from fourth on the original tuning sample to
second by three-sample mean T3 F1, while retaining the highest recall. IsoQuant
still leads T3 F1 by 2.77 points; the remaining gap is precision, not discovery.
The precision three-sample row retains the pre-family holdout artifacts; v8
paired tuning changed precision-profile T1/T3 F1 by -0.014/-0.031, both within
the 0.10 acceptance guard. Precision holdouts were not rerun for this bugfix.

The v7 audit exposed a separate ownership bug: discovery families were not
propagated to global isoform-fraction selection. v8 assigns stable family IDs to
read/GTF chains and carries them through assignment, aggregation, and junction
merges. Paired live validation compares the new family denominator with the
historical overlap denominator under identical source:

| Sample | Family vs overlap T1 F1 | Family vs overlap T3 F1 | T3 recall / precision delta | Output delta |
| --- | ---: | ---: | ---: | ---: |
| tuning p00 | **+0.0076** | 0.0000 | 0.0000 / 0.0000 | 0 |
| replicate3_run1 | **+0.0218** | **+0.0059** | +0.0248 / -0.0056 | +19 |
| SIRV p00/full | 0.0000 | 0.0000 | 0.0000 / 0.0000 | 0 |

Both real samples passed the pre-registered 0.10-point maximum F1-loss guard, so family is
the default and `--isoform-fraction-locus overlap` remains the explicit legacy
fallback. The original 87/95 audit reconstructed families only after earlier
gates; persistent discovery families can retain connections through bridge
variants later removed by those gates, explaining why the live effect is much
smaller than that upper bound.

## 7. Resolved profile values

| Setting | `sirv` | `real-drna` | `real-drna-precision` |
| --- | ---: | ---: | ---: |
| `m2_metric` | auto: SUM p00, mean guided | `summed_llr` | `summed_llr` |
| sum margin/flank | 1 / 4 | 1 / 8 | 1 / 8 |
| novel abundance boundary | inclusive >=3.0 | strict >1.0 | strict >1.0 |
| novel mono hard-read floor | off | 5 | 5 |
| finalized junction consensus | off | 6 bp / 2 reads / 2x | 6 bp / 2 reads / 2x |
| generation `min_novel_reads` | 1 | 1 | 2 |
| `min_gtf_abundance` | 1.0 | 1.0 | 1.0 |
| `floor_gtf_abundance` | on | off | off |
| `min_isoform_fraction` | 0.01 | 0.01 | 0.01 |
| isoform-fraction locus | family | family | family |
| `max_soft_mass_ratio` | 2.0 | off | off |
| `min_fulllen_fraction` | 0.1 | off | off |
| `min_polya5p_reads` | off | off | off |
| `canonical_gate` | on | on | on |
| `novel_junction_min_reads` | 2 | 2 | 2 |
| `containment_cluster` | on | on | on |
| `m2_cluster_recheck` | on | on | on |

Every listed setting remains explicitly overridable. `custom` applies no profile
overlay; reproducing pre-family runs also requires
`--isoform-fraction-locus overlap`.

## 8. Production changes

Integrated into `fin/`:

1. Named profile registry plus `PipelineConfig.from_profile`.
2. CLI `--profile real-drna|real-drna-precision|sirv|custom`, defaulting to balanced real dRNA.
3. Explicit override precedence for numeric and boolean options.
4. `m2_metric=off|mean|summed_llr|auto`; SIRV auto resolves by GTF usability.
5. Tight summed windows, separate sum margin/flank, and same-intron-count scope.
6. Two-score abstention plus per-read invalid-candidate-sequence abstention; one
   ambiguous base can no longer poison a locus batch.
7. Cross-sample validated real mono support and finalized junction consensus,
   with a separate precision support floor instead of a hidden recall trade-off.
8. Mass-preserving merges after junction correction.
9. Stable structural BLAKE2b novel IDs instead of run-random UUIDs.
10. `PipelineConfig.validate()` on the actual CLI path.
11. Atomic `run_manifest.json` with profile, route, overrides, config, source
    SHA-256, commit when available, and result-changing environment.
12. Secondary BAM alignments excluded from candidate generation.
13. Benchmark harness recovers original raw FASTQ/BLOW5 when mapped-subset links
    are unavailable and handles nested automount binds.
14. Benchmark smoke tests use a discovered/current Python executable.
15. M3 read-by-read DTW removed from production config, CLI, assignment, and
    active ablations; its prototype remains under `experiments/m3_coherence/`.
16. Stable splice-family IDs propagated through discovery/quantification and
    selected by default for isoform-fraction denominators, with overlap fallback.

Post-removal v7/v5 validation produced byte-identical `assembly.gtf` and
`scores.tsv` files for real balanced, real precision, SIRV p00, and SIRV full.
Manifests intentionally differ because the three retired config keys are absent
and the source hash changed. These v7/v5 jobs are correctness reruns, not new
runtime estimates: v7 balanced includes the family audit and scheduler/runtime
variance affected v7 precision. The table retains representative v6 timing for
the byte-identical outputs.

Kept optional or rejected as defaults:

- Static all-SUM SIRV: rejected; SUM is used only for unguided runs, mean for guided runs.
- Full-length/polyA real gates: disabled due severe recall loss/cost.
- `min_abundance=2/3` and mono floors 9-11: available precision points, not the balanced default.
- Generation `min_novel_reads=2`: promoted only to `real-drna-precision` because
  three samples show about 1 pp lower T1/T3 recall.
- Junction minimum 1: available recall-heavy point, not a balanced default.
- M3: removed from production; the generic EM engine remains at `beta=0` and
  the prototype is retained only for a future masked redesign.
- M4, candidate-align existence gates, family-path exploration, and TSS/TES/APA
  recovery: no winning end-to-end evidence for the final mode.
- Cluster quantification: lost to `m2_em` in the existing experiment record.

## 9. Reproducible artifacts

- Driver: `experiments/prod_validation/profile_sweep.py`
- SIRV one-factor results:
  `experiments/prod_validation/sirv4/_goal_opt/profile_sweep/results.tsv`
- SIRV combination results:
  `experiments/prod_validation/sirv4/_goal_opt/combo_sweep/results.tsv`
- SIRV holdout:
  `experiments/prod_validation/sirv4/_goal_opt/holdout/results.tsv`
- Six-sample true-p00 SUM/mean matrix:
  `experiments/prod_validation/sirv4/_goal_opt/sum_guide_split_true_p00/results.tsv`
- Final SIRV auto manifests:
  `experiments/prod_validation/sirv4/_goal_opt/final_auto_v6/baseline/SGNex_H9_directRNA_replicate3_run2/p00/run_manifest.json`
  and the sibling `full/run_manifest.json`
- Real one-factor aggregate:
  `experiments/prod_validation/gencode/_goal_opt/profile_sweep_retry/real_one_factor_results.tsv`
- Real sum grid:
  `experiments/prod_validation/gencode/_goal_opt/real_sum_grid/real_sum_grid_results.tsv`
- Independent real holdouts:
  `experiments/prod_validation/gencode/_goal_opt/real_holdout_r3r1_v2/results.tsv`
  and `real_holdout_r4r2/results.tsv`
- Cross-sample competitor rescoring:
  `experiments/prod_validation/gencode/_goal_opt/competitor_p00_score*/results*.tsv`
- Family/overlap paired acceptance results:
  `experiments/prod_validation/gencode/_goal_opt/final_profiles_v8/balanced/results.tsv`
- Original family/isoform-fraction audit:
  `experiments/prod_validation/gencode/_goal_opt/final_profiles_v7/balanced/family_fraction_audit.json`
- Final balanced manifests:
  `experiments/prod_validation/gencode/_goal_opt/final_profiles_v8/balanced/baseline/`
- Final precision result/manifests:
  `experiments/prod_validation/gencode/_goal_opt/final_profiles_v8/precision/`
- Deterministic precision-repeat evidence (pre-canonical-consistency source):
  `experiments/prod_validation/gencode/_goal_opt/final_profiles_v5/precision_a/`
  and `precision_b/`
- Frozen SLURM provenance source:
  `experiments/prod_validation/gencode/_goal_opt/profile_code_20260826/`
  (`SOURCE_SHA256.txt` =
  `3c63af34bb954df6d47b5f958ce351e3a95a77a81a9859710f5712532109c0be`)

The frozen source is not authoritative production code; it is retained only to
make the parallel one-factor/sum-grid runs byte-reproducible. Live source and
final manifests remain authoritative.

The initial 16-worker/120-GB real sweep attempt was OOM-killed because every
worker loads the 3-GB genome. The completed sweep used 8 workers and 200 GB per
job. This is benchmark scheduling evidence, not a scientific profile setting.
