# PyFIN Performance Strategy: Winning Every Benchmark State

Goal: make PyFIN's results superior to every fairly-comparable tool in every
evaluated state, not merely on average. This document defines what "every
state" means, records the authoritative current standings, explains — at the
mechanism level, from live artifact analysis — exactly where PyFIN wins and
loses today, and lays out a prioritized, evidence-gated roadmap.

Baseline identity for every number attributed to "current PyFIN" unless noted:
commit `04fb802`, source SHA-256
`90c8117a958a2b0653ecef1235d6521ebc0a614315f5ea6e8cc55f1909e63eeb`
(byte-identical GTF/TSV to the validated `05e303ab` v9 release).

The theory batch documented below (BAM completeness, M2 observability +
`sqrt_count_mean_llr`, candidate evidence layer, evidence ranker, EndpointRefine,
lazy genome) closes at working-tree source SHA-256
`9ea0b8858dda7456123b8a56f66ddaf03b2669c5fa453edb5a5ace70fe51fdc4`
(uncommitted). Closure evidence, all at that hash and sharing one manifest
hash across seven runs
(`experiments/prod_validation/gencode/_goal_opt/goal_closure/`):
ranking-off and SIRV outputs byte-identical to the v9/v7 releases (the whole
batch is numerically inert with every new mode off), two independent
`--ranking-mode filter` runs byte-identical, two independent
`--endpoint-refine` runs byte-identical across all four artifacts with mass
error 1.46e-10, and peak RSS 40,659 MB versus 67,110 MB eager.

---

## 1. What "all states" means (the verifiable definition)

A claim of superiority is only meaningful over a pre-registered matrix. The
matrix has seven axes:

| Axis | Values |
| --- | --- |
| Data domain | SIRV spike-in (absolute truth) / human real dRNA (GENCODE pseudo-truth) |
| Guide state | `p00` (de novo), `p10`, `p50`, `p90`, `p99`, `full` |
| Guide corruption | `skip`, `jitter`, `spurious`, `merge`, `flip`, `IR` at 5-20% |
| Biological holdout | H9 r2/r3/r4, HEYA8, K562, MCF7, HCT116 (SG-NEx) |
| Expression operating point | T1 (>=1 est. read), T3 (>=3 est. reads) |
| Runtime state | CPU/GPU, thread counts, read depth |
| Competitors | IsoQuant, StringTie3, Bambu, ESPRESSO (fair); TALON/IsoTools excluded from p00 verdicts until their wrappers run true de novo |

**Winning definition — two distinct gates:**

1. *Release safety (non-inferiority):* in each benchmark cell, the same
   single `auto` invocation of PyFIN scores an F1 not below the best
   fairly-run competitor in that cell.
2. *Superiority claim:* the paired per-cell ΔF1 is positive with a bootstrap
   confidence interval excluding zero. A delta whose CI includes zero is
   reported as a statistical tie, never as a win.

Averages never excuse a losing cell. SIRV (absolute truth) and GENCODE
(pseudo-truth) are reported in separate tables and never mixed into one mean.

**What cannot be promised:** mathematical superiority on arbitrary unseen
distributions. What can be promised and verified: cell-by-cell superiority on
the pre-registered matrix, with independent holdouts, bootstrap confidence
intervals, and a fixed truth/scorer contract to bound extrapolation risk.

Scoring contract (already enforced, must never drift): same BAM, same truth
GTF, same NanoCount expressed-truth sets, same gffcompare version and
`profile_sweep.py` scorer, competitor GTFs scored byte-for-byte as produced by
their own pipelines.

---

## 2. Authoritative current standings

### 2.1 SIRV domain (absolute truth) — PyFIN leads the historical grid

Historical full-grid (`experiments/prod_validation/sirv4/tables/`, mean
honest F1@3 across samples; `pyfin_prod` there predates the current auto
profile):

| Tool | Global mean F1@3 | Sn@3 | Pr@3 |
| --- | ---: | ---: | ---: |
| **PyFIN** | **85.0** | 92.2 | 78.9 |
| ESPRESSO | 78.7 | 79.2 | 81.9 |
| StringTie3 | 76.6 | 78.7 | 75.6 |
| IsoQuant | 72.5 | 66.9 | 82.0 |
| Bambu | 50.4 | 77.2 | 42.3 |

In that historical grid PyFIN leads every ratio and corruption row except
`full`/`p99` (-0.1 vs ESPRESSO 90.6), `c_merge5` (-3.0), `c_merge20` (-4.6).
Comparability caveats: those rows are **means across samples**, so a won mean
does not prove every underlying sample cell was won; and the current auto
profile's **92.567 (full)** / **78.988 (p00)** are single-holdout results at
the current source, not directly comparable to the historical means. The
whole SIRV grid (ratios, corruption rows, per-sample cells) must be re-run at
the current source in Milestone 0 before any SIRV victory claim; `c_merge*`
is the only SIRV region where PyFIN has ever trailed and gets priority
scrutiny there.

### 2.2 Real dRNA p00 (GENCODE pseudo-truth) — the contested domain

Honest F1 under the shared contract (competitor GTFs rescored on the same
BAM/truth/scorer; `_goal_opt/competitor_p00_score*/`):

**Tuning sample (H9 replicate2_run2):**

| Tool | T1 F1 | T3 F1 | T3 recall | T3 honest precision |
| --- | ---: | ---: | ---: | ---: |
| IsoQuant | 33.654 | **42.134** | 34.377 | 54.411 |
| **PyFIN precision** | **38.854** | 40.510 | **44.888** | 36.910 |
| Bambu | 36.848 | 39.439 | 40.043 | 38.852 |
| **PyFIN balanced** | 38.699 | 39.422 | 46.149 | 34.406 |
| StringTie3 | 38.257 | 38.681 | 45.280 | 33.761 |
| ESPRESSO | 14.369 | 20.818 | 12.903 | 53.850 |

**Holdout r3r1 (H9 replicate3_run1):**

| Tool | T1 F1 | T3 F1 |
| --- | ---: | ---: |
| IsoQuant | 33.994 | **41.429** |
| StringTie3 | **38.355** | 39.579 |
| Bambu | 35.514 | 38.645 |
| **PyFIN balanced** | 37.863 | 38.602 |
| **PyFIN precision** | 38.044 | 39.732 |
| ESPRESSO | 16.792 | 22.791 |

**Holdout r4r2 (H9 replicate4_run2):**

| Tool | T1 F1 | T3 F1 |
| --- | ---: | ---: |
| IsoQuant | 32.700 | **41.394** |
| **PyFIN precision** | **38.580** | 39.788 |
| **PyFIN balanced** | 38.415 | 38.622 |
| StringTie3 | 37.527 | 37.258 |
| Bambu | 36.057 | 37.365 |

(PyFIN r3r1/r4r2 rows: fresh runs at source `90c8117a`, manifests confirm
refit on; `_goal_opt/plan_baseline/`.)

**Current verdict, real p00 (precision profile, the stronger of the two):**

| Deficit cell | Gap | Consistency |
| --- | ---: | --- |
| T3 vs IsoQuant — tuning | -1.624 | uniform ~1.6-1.7 across all three samples |
| T3 vs IsoQuant — r3r1 | -1.697 | |
| T3 vs IsoQuant — r4r2 | -1.606 | |
| T1 vs StringTie3 — r3r1 only | -0.311 | tuning +0.60 and r4r2 +1.05 ahead |

Every other real-p00 cell is already won: precision beats StringTie3/Bambu at
T3 on all three samples (r3r1 +0.15/+1.09, r4r2 +2.53/+2.42) and everyone but
StringTie3-on-r3r1 at T1. The T3-vs-IsoQuant gap is remarkably stable
(-1.6 to -1.7), consistent with a single systematic cause — the
selection-precision mechanism in Section 3 — rather than sample-specific
noise.

### 2.3 Known matrix holes (cannot claim victory before they are filled)

1. Real guided ratios `p10..full`: PyFIN guided profiles exist but have no
   validated numbers; competitor GTFs existed at all ratios in the
   NanoRNATrans sweep. (Historical guided-real spot data —
   `comparison.tsv` — shows guided PyFIN ahead of guided IsoQuant per-sample,
   but under an older scorer; not verdict-grade.)
   **Data availability update (2026-08-29):** the external
   `gencode_full_sweep` results tree (competitor GTFs, NanoCount oracle
   tables, r3r1/r4r2 BAM/fastq/blow5 targets) was deleted outside this
   repository while the user prepares new data. Every number in this
   document was computed and archived on /SSD before the deletion
   (evidence tables, filtered GTFs, gffcompare outputs, audits). Until the
   new data arrives, live re-validation is only possible on the local tuning
   sample and SIRV; r3r1/r4r2 live pipeline-filter validation, guided
   ratios, non-H9 lines, and competitor re-scoring are all data-gated.
2. SIRV corruption rows with the current source (only p00/full re-validated).
3. Non-H9 cell lines (HEYA8/K562/MCF7/HCT116): competitor GTFs exist; PyFIN
   has no current-source p00 numbers.
4. TALON/IsoTools p00 wrappers consume a one-transcript stub, so their p00
   numbers are unfairly low; excluded from verdicts until fixed.
5. Runtime: competitor runtimes are historical and used different
   hardware/threads; no fair runtime table exists yet.

---

## 3. Mechanism-level gap analysis (from live artifacts, not intuition)

All numbers below are from tmap-level analysis of the released precision run
on the tuning sample against the same-contract IsoQuant rescore
(`gc.assembly.gtf.tmap`, normalized reference IDs).

### 3.1 The retained output has the recall to win; the gap is selection precision

- PyFIN precision matches **3,771** unique T3 truth IDs; IsoQuant matches
  **2,888**.
- Overlap decomposition: both 2,762; **PyFIN-only 1,009**; IsoQuant-only 126.
- Of the 126 IsoQuant-only truths, PyFIN is entirely absent at 104 and has a
  near-miss (`j`/`c`/`m`) at 22.

PyFIN's retained output covers ~35% more expressed true structures than the
T3 leader, so the **primary cause of the current relative T3 gap is
selection/ranking precision**: PyFIN emits 10,217 rows of which 2,352 (23.0%)
are reference-non-exact, while IsoQuant emits 5,309 with 463 (8.7%) —
gffcompare structural precision 72.0 vs 91.3. Two scoping limits: this
final-output analysis cannot see candidates removed earlier in
discovery/selection (upstream stages may hold additional recoverable recall),
and discovery remains a hard recall ceiling for the 104 IsoQuant-only truths
PyFIN never emits and for truths both tools miss (T3 truth is 8,401; both
tools together match only 3,897).

### 3.2 Reference-non-exact taxonomy (the "error budget")

Terminology note: a non-`=` gffcompare row is *reference-non-exact*, not
necessarily a biological false positive — GENCODE is a pseudo-truth and some
`j` rows are real novel isoforms. Under the agreed metric they still cost
precision, so they are the budget to reduce; deciding which are truly wrong
needs orthogonal evidence (SIRV analogues, cross-sample recurrence,
independent long-read support).

PyFIN precision tuning output, 2,352 non-exact rows:

| Class | n | Mono/multi | Interpretation |
| --- | ---: | --- | --- |
| `j` (chain mismatch) | 939 | all multi (932 standalone, 7 sibling-subchain) | hypothesis: wrong chains and/or real novel isoforms |
| `c` (contained) | 762 | 377 mono / 385 multi | hypothesis: 5' degradation fragments kept standalone |
| `k` (contains ref) | 418 | 58 mono / 360 multi | hypothesis: extensions/readthrough beyond truth chain |
| `o/e/i/n/u/x/p/m` | 233 | mostly mono | hypothesis: mono overlap/intronic/intergenic artifacts |

Counting note: this table and the distribution statistics below are
**query-row-level** (a reference transcript can be matched by more than one
query row; e.g. 4,096 exact-T3 rows collapse to 3,771 unique T3 truth IDs).
The authoritative unique-truth-ID counts remain those of Section 3.1.
Relatedly, the exact-row ratio 7,865/10,217 = 77.0% is not the official
gffcompare transcript precision (72.0) — gffcompare de-duplicates and counts
differently; official honest-F1 numbers always come from `gc.stats` + the
scorer, never from row ratios. Every number in this section is reproducible
via `experiments/prod_validation/t3_gap_decomposition.py`, which fixes the ID
normalization, truth thresholds, and row-vs-ID semantics; the exact
invocation backing this section is:

```
python3 experiments/prod_validation/t3_gap_decomposition.py \
  --pyfin-dir experiments/prod_validation/gencode/_goal_opt/\
final_profiles_v9_release/precision/baseline/gencode_p00/p00 \
  --competitor-tmap experiments/prod_validation/gencode/_goal_opt/\
competitor_p00_score/isoquant/gc.assembly.gtf.tmap
```

Evidence boundary: GENCODE class codes label *annotation agreement*, and
NanoCount labels *oracle expression* — neither is biological truth. Before a
non-exact family is treated as biologically wrong (rather than
metric-penalized), it needs orthogonal evidence: SIRV absolute-truth
analogues, cross-sample recurrence, and independent junction/end support. A
ranker trained on these labels optimizes the agreed metric, and that
limitation is accepted explicitly.

Aggregate mono non-exact rows: **595**. IsoQuant emits **zero** mono
non-exact rows (it hard-drops novel mono at p00); PyFIN keeps novel mono at a
hard-read floor of 5. This single policy difference is a large, quantified
chunk of the raw-precision gap.

### 3.3 Negative results that discipline the roadmap (measured, not assumed)

1. **Abundance thresholding cannot close the gap.** Ranking all outputs by
   abundance and cutting the lowest 1,000 removes only 430 non-exact rows at
   the cost of 570 exact rows (~0.75 true per false removed). Distributions
   overlap: median abundance TP-T3 13.0 vs `j` 5.2 vs sub-T3 exact 7.0.
   Threshold tuning is exhausted; this is why v5-v9 plateaued.
2. **Observable sibling-containment cannot be hard-dropped.** Only 35/762 `c`
   rows are subchains of a surviving sibling (the parent was usually already
   removed); a hard drop of sibling-contained multi rows removes more exact
   rows (85) than non-exact (56).
3. **`confidence`/`max_R` are saturated** (median 1.0 in every category
   except `c` at 0.721) and cannot rank. The one usable signal found:
   contained-fragment rows have visibly depressed EM confidence.
4. **Post-refit fixed-point filtering is settled** — enforcing gate
   consistency after the abundance refit deletes 2-5 real transcripts,
   converts their reads to orphan mass, improves nothing (T3 F1 -0.0027 on
   r3r1), and is kept off (`V9_REVIEW_NOTES.md`).
5. **EndpointRefine upper bound is small right now**: endpoint-equivalent
   splitting covers only 62/4,477 (1.38%) of missed T3 multi-exon truths. It
   stays behind junction-chain precision in priority.
6. **Exact-but-unexpressed matches dilute honest precision at both
   thresholds.** Exact query-row strata on the tuning run
   (`t3_gap_decomposition.py`): 4,096 rows / 3,771 unique IDs are
   T3-expressed; 1,142 rows / 1,108 IDs are T1-but-not-T3; **2,627 rows /
   2,477 IDs are below T1** — real annotated structures the NanoCount oracle
   deems unexpressed (<1 est. read) that PyFIN nonetheless emits (median 8
   assigned reads: an assignment disagreement with the oracle, not a
   structural error). The below-T1 stratum dilutes honest precision at BOTH
   T1 and T3 while contributing recall at neither. Its median abundance
   (8.0) is *higher* than the T1-only stratum's (5.0), so no abundance floor
   can isolate it — expression calibration against the oracle's assignment
   model, not a threshold, is required, and its removal is expected (not
   guaranteed: honest precision couples raw gffcompare precision with the
   expressed-match ratio, so only a rescored frontier can confirm) to raise
   honest precision at both thresholds. The T1-only stratum genuinely
   carries T1 recall and must be protected: PyFIN's T1 standing vs the best fair competitor (StringTie3)
   is competitive, not safe — +0.597 (tuning), **-0.311 (r3r1)**, +1.053
   (r4r2).

### 3.4 The winning geometry

To beat IsoQuant's tuning T3 F1 of 42.134 at PyFIN's current recall (44.888),
honest precision must rise from 36.910 to >=39.7. A back-of-envelope
heuristic — removing roughly 900-1,100 of the 2,352 non-exact rows at
near-zero exact-row loss — sizes the job, but the real requirement is a
rescored precision-recall frontier, since honest precision couples raw
gffcompare precision with the expressed-match ratio and cannot be predicted
exactly from row deletions. PyFIN's T3 recall surplus over IsoQuant
(+10.24 to +10.51 points: tuning 44.888 vs 34.377, r3r1 42.936 vs 32.657,
r4r2 44.342 vs 34.100) is a war chest: even sacrificing 1-2 recall points for 5+ precision points wins T3
while T1 stays ahead — **provided** the removals concentrate on true
structural errors, which Sections 3.2-3.3 show requires new evidence, not
thresholds. The same lever is expected to address the one T1 deficit cell
(StringTie3 on r3r1, -0.31): removing non-exact rows and below-T1 exact rows
should raise honest precision at T1 as well while leaving T1-expressed
matches untouched — subject to the same rescored-frontier confirmation. The
uniformity of the IsoQuant gap
(-1.6 to -1.7 on all three samples) predicts that a lever validated on
tuning + r3r1 will transfer to r4r2; the acceptance contract still requires
proving it.

---

### 3.5 Repository architecture and workstream ownership

Complete module map: `REPO_MAP.md` (71 modules, 121 import edges, 89 formal test
files); algorithm formalization: `ALGORITHM.md`; independent review with
resolved/open findings: `ALGORITHM_REVIEW.md`; profile tuning evidence:
`PROFILE_OPTIMIZATION.md`. The pipeline stages and where each roadmap
workstream lands:

| Stage | Modules | Owning workstream |
| --- | --- | --- |
| CLI / profiles / config | `fin/cli.py`, `fin/pipeline/config.py` (138 fields, 113 options, profile overlays, `validate()`) | P5 adaptive auto |
| Interval generation | `fin/io/interval_manager.py` (read clustering, strand-split intervals, fusion detection) | P1 (region-fetch completeness) |
| Discovery / candidates | `fin/candidates/` (intron chains, `cluster_families` single-linkage families, chain clustering, canonical gate, stable BLAKE2b IDs) | P3 (`c`/`k`/mono family models) |
| Assignment / EM | `fin/pipeline/assignment.py`, `fin/analysis/assignments.py` (`em_with_coherence`, production `beta=0` one-softmax), `fin/scoring/m2_junction_nll.py` (M2 junction NLL, mono resolve, containment collapse), `fin/scoring/em_inputs.py` | P1 (M2 calibration), P2 (ranker features) |
| Selection | `fin/pipeline/selection.py` (interval cascade + `select_global`: abundance floors, family isoform fraction, soft-mass, mono, full-length, polyA) | P2 (replace hard gates with evidence ranking) |
| Finalization | `fin/pipeline/junction_snap.py` (consensus snapping + merge redirects), `fin/analysis/abundance_refit.py` (v9 survivor refit, sparse ledgers, mass conservation) | P4 (EndpointRefine builds on refit), P6 (ledger memory) |
| Writers / scoring | `fin/io/io_gtf.py`, `fin/io/io_tsv.py`, `fin/pipeline/finalize.py` (no filtering in writers — architectural guarantee) | P2 (ranker score column; probability calibration future) |
| Parallelism / GPU | `fin/pipeline/parallel.py` (spawn pool, canonical aggregation order), CuPy EM (GPU SIRV / CPU real) | P6 (RSS ~65 GB, CSR ledgers, batching) |
| Validation harness | `experiments/prod_validation/` (`profile_sweep.py` runner+scorer, refit/gate/fixed-point/ledger probes, SBATCH acceptance scripts) | P0 (matrix + diagnostic tables) |
| Tests | `tests/unit` (1,032), `tests/integration` (15), retired-M3 prototype (5) | every promotion gate |

Key dependencies: pysam (BAM), mappy (SIF-only), NumPy/CuPy (EM), krill
(signal/eventalign backend), Singularity image
`pyfin_gpu_e268c9b.sif` as the reproducible runtime; gffcompare (Docker) +
NanoCount tables as the scoring contract.

---

## 4. Roadmap (prioritized, evidence-gated)

Ordering principle: measurement infrastructure first, then evidence quality,
then an evidence ranker, then per-error-family structure models, then
endpoints/adaptivity, then speed. No milestone may promote a default without
passing the acceptance contract in Section 5.

### P0 — Complete the fair benchmark matrix and error taxonomy (Milestones 0 and 1)

- Run current-source PyFIN (`auto`/balanced/precision) on: r4r2 + non-H9 cell
  lines at p00; guided ratios p10..full on at least tuning + one holdout;
  SIRV corruption grid (`c_skip/jitter/spurious/merge/flip/ir` at 10/20%).
- Rescore TALON/IsoTools only after their wrappers run true de novo.
- **M0 deliverable:** every matrix cell filled at the current source with
  manifests.
- **M1 deliverable:** one per-candidate diagnostic table per run (class code,
  exon count, annotation status, expression stratum, family, junction
  support, M2 coverage/abstention, canonicality, end support) plus the
  orthogonal-evidence split of each non-exact family into "truly wrong" vs
  "novel-but-real" (SIRV analogues, cross-sample recurrence, independent
  junction/end evidence) — the substrate for P1-P3.
- Freeze sample roles now, honestly accounting for what this planning work
  already opened: tuning = H9 r2r2; threshold-fitting validation = r3r1;
  secondary validation/audit = r4r2 (its PyFIN and competitor numbers were
  inspected while writing this plan, so it can no longer serve as an
  untouched holdout); **untouched final holdouts = K562 replicate5_run1 +
  MCF7 replicate4_run1** (never inspected in this work). Final holdouts are
  spent once, at the end, for the promotion decision only.
- Record competitor versions, commands, input/truth/scorer hashes, and
  matched hardware/threads for the runtime table.

### P1 — Make M2/junction evidence trustworthy (Milestone 2) — partially DONE

Status 2026-08-29: the BAM completeness fix shipped (`RegionReadBatch`;
partial region scans now abstain instead of counting as zero support — review
Medium 9 resolved); inert single-softmax EM knobs now warn on non-default
values (Medium 7 resolved); M2 tight-window contrasts log per-interval
observability counters (event counts, same/diff-count, margins, abstention;
optional JSONL via `m2_contrast_stats_jsonl`), and the experimental
`sqrt_count_mean_llr` metric (mean NLL x sqrt(min event count)) is selectable
but never auto-routed. First observation point (single sample, NOT promotion
evidence): SIRV p00 T3 F1 79.464 vs summed 78.988 (+0.48), and live counters
confirm unequal event counts occur even inside the same-intron-count niche
(e.g. one interval decided=17 with diff_events=5) — the calibration concern
is real and the observation phase has begun. Promotion still requires the
full SIRV pair grid plus real-sample contrasts.

- Fix partial BAM region fetch to abstain instead of emitting false
  hard-negative junction evidence (open High-severity correctness item).
- Calibrate summed-LLR across unequal hypothesis event counts; report
  coverage, abstention, pair accuracy, and reliability curves per sample.
- Remove or reject inert `m2_em` iterative-EM/prior options so configuration
  space equals behavior space.
- Acceptance: `j`-class non-exact count drops on tuning AND r3r1 with T1
  recall unchanged; direction consistent across samples.

### P2 — Evidence-ranked candidate selection (Milestone 3, the centerpiece) — v1 DONE (experimental mode)

Status 2026-08-29 (all numbers real-gffcompare, frozen-model discipline):

* Feature layer `fin/analysis/candidate_evidence.py` + one-pass whole-BAM
  collector; `--candidate-evidence` writes the per-survivor observable table.
* Frozen L2-logistic ranker `fin/analysis/candidate_ranking.py`
  (provenance: `experiments/prod_validation/models/candidate_ranker_v1.json`,
  bit-identical constants enforced by unit test). Trained on tuning only;
  chromosome-grouped CV AUC 0.809; threshold frozen on the tuning frontier
  under the T1-not-lower constraint before any validation sample was scored.
* Results at the frozen operating point (offline filter of final GTFs, real
  scorer, while the oracle tables were still available):

| Sample | T1 F1 | T3 F1 | vs IsoQuant T3 | vs best-competitor T1 |
| --- | --- | --- | --- | --- |
| tuning | 38.854 -> 39.173 | 40.510 -> 42.466 | **+0.33** (was -1.62) | +0.92 |
| r3r1 (frozen validation) | 38.044 -> 38.600 | 39.732 -> 42.077 | **+0.65** (was -1.70) | +0.25 (flips the last T1 deficit) |
| r4r2 (one-shot audit) | 38.580 -> 38.256 | 39.788 -> 40.758 | -0.64 (was -1.61) | +0.73 |

* Live pipeline validation on tuning: `--ranking-mode filter` runs
  filter->snap->refit; output differs from the offline simulation by exactly
  one transcript (8,217 vs 8,216; a snap-representative changed because
  filtering precedes snapping), mass conservation error 2.9e-11, SIRV output
  byte-identical with ranking off, manifests record the mode.
* Per the pre-registered promotion gate (r4r2 T1 -0.32 below its own
  baseline; T3 gap to IsoQuant halved but not closed), the ranker ships as an
  EXPERIMENTAL opt-in (`--ranking-mode filter`), not a profile default;
  default-on awaits the reserved unopened holdouts (K562/MCF7) and the
  user's incoming data.

Original plan follows.

Replace saturated confidence with a monotonic, calibrated score built ONLY
from inference-time observables:

- weakest-junction support; unique/perfect junction reads;
- M2 margin, effective event count, abstention state;
- family-relative abundance share and assignment entropy;
- 5'/3' read-end mode agreement; poly(A)-supported TES;
- canonical/noncanonical junction evidence;
- containment/extension geometry vs family siblings;
- primary/supplementary/chimeric read-support composition (NOT in ranker v1:
  the one-pass collector counts them globally, but no per-candidate feature
  is computed yet — a v2 candidate feature).

Discipline: train/evaluate with sample-level splits (never transcript-random
splits — locus leakage), labels from gffcompare class + expression stratum on
tuning/validation samples only; GENCODE class codes and NanoCount counts are
forbidden as inference features. Cross-sample recurrence is a cohort-only
feature: it defaults to missing/neutral in ordinary single-sample runs and is
never computed against the holdout being scored. Features, labels, operating
points, and `auto` routing are frozen before any final holdout is opened;
competitor results on final holdouts are evaluation-only and never tune
anything. Output: a raw-logit ranker score (probability calibration is future work) and a movable operating
point per profile.
Acceptance: at fixed T1 (>= current per-sample values), T3 honest precision
rises enough to clear IsoQuant per-sample T3 F1 on tuning AND validation,
with bootstrap CIs excluding zero.

### P3 — Per-error-family structural models (Milestone 4) — evidence DONE, dedicated models data-gated

Status 2026-08-29: every error family named below is now represented as
inference-time evidence inside the ranker feature layer rather than as a new
hard gate — `j` by weakest/median junction support and
`n_junctions_below3`, `c` by `is_subchain_of_sibling` plus 5'/3' end-mode
agreement, `k` by `is_superchain_of_sibling` and length/exon geometry, mono
by `is_mono` with end-support fractions (fitted weights: junction support
+0.563, canonical +0.361, mono -0.371). That is the deliberate ordering the
plan asked for (features first, then dedicated models only after each family
improves in the same direction across samples). The dedicated per-family
structural models below remain **data-gated**: proving a family-specific
lever needs per-family cross-sample direction evidence, and the external
competitor/oracle tree required for it was deleted (see the availability
note in Section 2.3).

- `j` rows: junction-chain consensus + per-junction evidence calibration
  (extends the P1 machinery to chain-level decisions).
- `c` rows: degradation-aware 5' fragment model — decide fragment-vs-isoform
  from end-mode distributions and family context, not suffix collapse alone.
  The 350 standalone multi `c` rows whose containing structure was already
  removed are the main target; the 35 with a surviving sibling need a
  separate audit of why containment collapse retained them — hard-dropping
  them is measured harmful (Section 3.3.2), so the answer is per-case
  evidence, not a blanket rule.
- `k` rows: extension/readthrough control with adaptor/chimera decoys.
- mono rows: policy informed by the 595-row budget — distinguish genuine
  single-exon genes (histones, ncRNA) from intronic/internal-priming
  artifacts via poly(A) support, end-mode sharpness, cross-sample recurrence;
  evaluate stricter default floors against the T1 cost explicitly rather
  than adopting IsoQuant's hard drop blindly.
- Noncanonical junctions: replace canonical-only exclusion with
  evidence-aware high-support exceptions.
- Fusion stays experimental and off by default.

### P4 — EndpointRefine (Milestone 5a) — experimental implementation DONE

Status 2026-08-29: implemented off-by-default (`--endpoint-refine`,
`fin/analysis/endpoint_refine.py`): strand-aware end-mode clustering,
supported (TSS, TES) pairs (>=3 reads and >=15% of end-mapped reads),
interior-TSS degradation guard (2x support), split cap, stable BLAKE2b
endpoint IDs, and mandatory post-split requantification — the split plan
emits per-read routes consumed by `refit_survivor_abundance`
(`split_routes`/`split_primary`), and `validate()` rejects the mode without
the effective refit. Poly(A)-supported TES uses the same krill whole-read
polyA pass as the polyA+5' gate when a signal file is configured
(poly(A)-confident reads strengthen non-primary TES votes); without signal
it degrades — documented and logged — to end-sharpness + guards only.
Per-run `endpoint_refine.json` records modes, pair support, poly(A) support,
and routes. Enabled in no profile: promotion is data-gated on the incoming
holdouts (measured upper bound of the addressed error class: 1.38% of
missed T3 multi-exon truths).

Original plan follows.

### P5 — State-adaptive `auto` (Milestone 5b) — mechanism exists, extension data-gated

Status 2026-08-29: the routing mechanism is in place and already carries the
one routing decision that has paired live evidence — `m2_metric=auto` selects
mean with a usable guide and summed_llr when unguided (`resolve_m2_metric`,
validated on six SIRV pairs). Extending `auto` to further axes (read
depth/error profile, junction-evidence coverage) is **not a theory gap but a
data gap**: each new routing rule needs the benchmark cell that proves which
setting wins in which state, and that matrix cannot be completed while the
external competitor/oracle tree is deleted. Adding unvalidated routing rules
would violate the promotion contract in Section 5.

One `auto` profile must serve the whole matrix, choosing behavior only from
observable state: guide usability (already routes SIRV SUM/mean), read
depth/error profile, junction-evidence coverage, signal availability. Never
from benchmark names or truth files. Acceptance: the single `auto` invocation
reproduces or beats the per-profile numbers in every matrix cell.

### P6 — Runtime and memory (Milestone 6) — attribution DONE, hypothesis refuted

Measured 2026-08-29 with `experiments/prod_validation/memory_attribution_profile.py`
(live tuning run, 8 threads, sampled per-process RSS + parent structure
sizes; artifact `_goal_opt/ranking_live_v1/memory_attribution.json`):

| Component | Size |
| --- | ---: |
| Whole-job peak (lower bound) | 67.1 GB |
| Spawn workers combined (max 9 procs) | **59.7 GB (~7.5 GB each)** |
| Parent peak | 7.5 GB |
| Genome FASTA dict (parent; every worker holds its own copy) | 3.1 GB |
| ResponsibilityLedger payload + input IDs | **0.14 GB** |
| Aggregated read-ID tuples | 0.02 GB |

**The CSR/COO ledger hypothesis is refuted** — ledgers are 0.2% of the peak.
The dominant cost is per-worker duplicated state, led by the full genome
dict (~3.1 GB x 9 processes ≈ 28 GB) plus per-worker BAM/signal/aligner
caches.

**Fix implemented and measured (`fin/io/lazy_genome.py`,
`lazy_genome=True` default):** the genome now opens as an indexed lazy
mapping holding at most `genome_cache_chroms` (2) chromosomes per process,
falling back to the eager dict when the FASTA cannot be indexed. Re-profiled
on the same tuning run: **whole-job peak 67.1 GB -> 40.9 GB (-39.1%)**,
workers 59.7 -> 36.5 GB, parent 7.5 -> 4.7 GB, wall time unchanged
(286 -> 282 s); real and SIRV GTF/TSV outputs stayed byte-identical. The
reduction was re-measured on the final source `e2856440` after the profiler
itself was corrected to report the lazy mapping's cache payload instead of
forcing a whole-genome traversal (`goal_final_v2/memory_final.json`:
40,893.6 MB peak, `genome_mode=lazy`, `genome_cached_mb=0.0`,
`ledgers_payload_mb=111.3`), and that run's `assembly.gtf` is byte-identical
to the v9 release. The remaining ~4.5 GB/worker (krill/aligner/BAM/signal
state) is the next profiling target.
- Then profile-guided: BAM/signal pass consolidation, worker reference
  duplication, GPU batching.
- Preserve canonical aggregation order and byte determinism (two independent
  runs byte-identical remains a release gate).
- Fair runtime table vs competitors on identical hardware/threads; historical
  runtimes are not verdict evidence.

---

## 5. Acceptance contract (applies to every promotion)

1. Paired live validation (change vs baseline, same everything else);
   structure, abundance, runtime reported separately.
2. Tuning + at least two independent real holdouts move in the same
   direction; SIRV p00/full/corruption rows show no regression.
3. T1 and T3 must each stay >= the current profile per sample; per-cell F1
   must be >= the best fair competitor — averages never excuse a losing cell.
4. GENCODE and SIRV tables stay separate; bootstrap CIs reported with the
   resampling unit being loci/genes (paired between tools), never individual
   transcripts — transcript-level resampling ignores within-locus correlation
   and understates variance; deltas inside measurement noise are not
   victories.
5. Source/input/truth/scorer/container hashes recorded in manifests; two
   independent runs byte-identical; unit/integration/M3/static suites green.
6. The final untouched holdout is evaluated once, only for the final
   promotion decision.

## 6. Decision gates and stop conditions

| Milestone | Gate to proceed | Stop/fallback |
| --- | --- | --- |
| M0 matrix | all cells filled, diagnostic tables generated | if guided-real or corruption cells show unknown regressions, fix before any algorithm work |
| M1 taxonomy | error families quantified per sample with orthogonal-evidence split | if `j` rows are mostly recurrent cross-sample (likely real novel), re-weight roadmap toward expression calibration instead of chain suppression |
| M2 evidence | `j` reduction at fixed T1 on 2 samples | if M2 calibration cannot separate, fall back to per-junction perfect-read gates only |
| M3 ranker | beats IsoQuant T3 on tuning+validation at fixed T1 | if no feature set separates (AUC ~0.5), stop; revisit after P3 adds features |
| M4 structure | per-family budget reductions replicate on holdout | any T1 loss > 0.10 F1 reverts that family's lever |
| M5 endpoints/auto | single `auto` >= per-profile everywhere | if endpoint splits do not requantify cleanly, keep families unsplit |
| M6 speed | RSS/runtime improved, outputs byte-identical | any accuracy or determinism change reverts |

## 7. Ceilings and honesty

- GENCODE is a pseudo-truth: some `j`/`u` rows are real biology; the metric
  still counts them against precision. Winning under this contract means
  matching annotation-shaped output, and the plan accepts that explicitly.
- NanoCount is the expression oracle; PyFIN's own abundance disagreeing with
  it is scored as dilution at T3 even when reads genuinely support the model.
- SIRV is absolute truth but synthetic: 5'-truncation tails and mono policies
  tuned there do not transfer to real tissue without revalidation (already
  institutionalized in the two-profile system).
- TALON/IsoTools p00 comparisons stay excluded until their wrappers are fixed;
  claiming victory over a crippled baseline is not victory.

## 8. Immediate next actions (in order)

1. ~~Finish and record the three in-flight baseline runs~~ Done: Section 2.2
   is complete at source `90c8117a`; deficit set = 3 uniform IsoQuant-T3
   cells + 1 StringTie3-T1 cell (r3r1).
2. M0: SIRV corruption grid + guided real ratios + one non-H9 line with
   current source; per-candidate diagnostic tables.
3. M1: orthogonal-evidence split of the 2,352-row budget (cross-sample
   recurrence + SIRV analogues) to size the truly-wrong fraction per family.
4. P1 BAM-completeness fix (also the open correctness item) and M2
   calibration.
5. P2 ranker prototype on tuning+r3r1 with sample-split evaluation.
