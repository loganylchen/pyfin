# PyFIN Optimization — Session Handoff

Last updated: 2026-09-13. This document is the authoritative handoff. Read it
TOGETHER with `SIRV4_IMPROVEMENT_PLAN.md` (evidence + acceptance gates per
lever) and `SIRV4_GRID_STUDY.md` (benchmark semantics). Treat live source and
generated artifacts as authoritative over any narrative.

**START HERE (2026-09-13).** Sections 1-3 and the ladder in section 2 are the
HISTORICAL SIRV-profile campaign (kept for provenance). The shipped production
default is the `real-drna` profile chosen by cross-dataset maximin (4d), with
the endpoint channels (4g), snap-merge (4l) and -- since 2026-09-11 -- the
support-gate attestation exemption at K=3 (4t). Current standing on real human
GENCODE vs stringtie3, the one real peer: ahead on 6 of 7 conditions, mean
+1.4 Sn / +5.8 Pr / +2.54 F1; the single open structural dimension is de novo
(`p00`, -0.30). Where things stand per dimension: 4u and
`docs/PRODUCTION_ALL_DIMENSION_AUDIT_2026-09-11.md`. Git state: section 6.
Next lever, chosen by measurement: section 7 (the K operating point; HEYA8
K=1 arm is SLURM 211963).

Parallel algorithm notebook (2026-09-09; not a SIRV benchmark result):
`docs/PRODUCTION_RAW_DTW_CASEBOOK.md` records the original chr19 seven-read
raw-DTW case plus five additional junction pairs / 57 reads, reproducible
inputs and paths, calibration flips, and a possible longer-template preference.
No production scoring defaults changed; no staging or commits performed.
Follow-up: `docs/PRODUCTION_PRIMARY_EVENT_DTW.md` records primary-only observed
event means vs uncalibrated pore means (64 reads; 2x2 event/raw and calibration
controls; no new eventalign). Treat numerical preferences as exploratory, not truth.
User-confirmed next-direction contract (2026-09-09):
`docs/PRODUCTION_DTW_SHARED_CALIBRATION_DECISION.md` — candidate A/B model means
must share that read's final primary scale/shift and exact krill built-in model;
observed raw/events remain unchanged. Raw vs event DTW is not yet decided.
This is a recorded future direction, not an implemented production migration.
Executed exploratory follow-up: `docs/PRODUCTION_RAW_DTW_EXPANSION.md` records
8 new pairs / 105 reads with exact built-in means and shared final calibration;
105 primary calls, zero candidate eventalign, 420 paths independently checked.
Experimental scripts support this mode; production remains unchanged.
Full frozen spike-in validation (2026-09-09):
`docs/PRODUCTION_JUNCTION_ACCURACY_BENCHMARK.md` — all 6 SIRV4 + 4 Sequin
libraries, 10,521 reads / 42,315 read-junctions, annotation-blind prediction.
The minimum raw-DTW SSE rule FAILS the primary metric: SIRV4 exact catalog-valid
read-junction fraction 93.144% -> 71.377% (-21.767 pp); Sequin 93.775% -> 93.498%
(-0.277 pp, CI includes zero). SIRV4 raw fixes 641 but damages 8,988 observations;
95.438% of damaged observations move to longer spliced templates. This is an
experimental junction re-selection test, NOT production transcript F1. 68 unit
tests and full numerical/input audit pass; production is unchanged. Do not
promote this unpenalized minimum-SSE rule based on the earlier human examples.
Follow-up cached ablation (2026-09-09):
`docs/PRODUCTION_DTW_NORMALIZATION_ABLATION.md` — all ten libraries, 16 frozen
exploratory settings, zero new eventalign. Dividing SSE by k-mer count worsens
SIRV4/Sequin to 63.447%/92.515%. Minimum 5 raw samples per model state improves
unconstrained DTW to 81.521%/94.103%, but remains below SIRV4 BAM; Sequin gain
over BAM is +0.328 pp with CI [-0.254,+0.928]. Strong length penalties reverse
the bias: alpha 1 gives Sequin 84.274%, with 412/415 damaged observations moving
to shorter templates. No new abstentions; 84 tests and full score/path audit
pass. These already-inspected libraries are development data, not held-out
validation. Keep production unchanged.

## 1. Mission and scoring contracts

Benchmark: `/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/sirv4_full_sweep_ratio_grid_v2`
— 6 H9 dRNA samples x 26 annotation conditions (`p00`..`full` richness bands +
corruption classes `c_trunc*/c_merge*/c_ir*/c_skip*/c_flip*/c_jitter*/c_spurious*`)
= 156 cells. Eight competitors (ESPRESSO, StringTie3, Bambu, IsoQuant, FLAIR,
IsoTools, LAFITE, TALON) ran the same inputs; external coverage 1224/1248
cells (18 LAFITE + 6 TALON outstanding).

Two contracts — never pool them:
- **Expression contract (the optimization target)**: truth = per-sample
  NanoCount TPM>0 subset (`_ref/expressed/<sample>.gtf` staged under
  `experiments/prod_validation/_sirv4grid/_ref/expressed/`). Emitting an
  unexpressed catalog entry is a FP. NanoCount is same-read pseudo-truth:
  scoring only, NEVER a PyFIN input/feature.
- **Official catalog contract**: truth = all 176 SIRV_Set4 transcripts.
  Rewards catalog pass-through. `--catalog-mode complete|adaptive` are
  benchmark adapters for this contract only (adaptive: 111/156 top).

"Top/tied N/156" = cells where PyFIN F1 >= best available competitor F1.

## 2. Result ladder of the SIRV-profile campaign (HISTORICAL; expression contract, 156-cell runs)

> Superseded on 2026-09-05 by 4d: the `sirv` profile and the P0 / rescue /
> election levers are benchmark working points, NOT the shipped default. They
> win on SIRV4 and lose on the held-out Sequin set. The numbers below remain
> the reproducible record of that campaign.

| Stage | Flags added | Sn | Pr | F1 | Top/tied |
|---|---|---:|---:|---:|---:|
| baseline `off` | (sirv profile) | 83.32 | 88.16 | 85.58 | 65/156 |
| P0 | `--min-gtf-abundance 0 --no-floor-gtf-abundance --require-gtf-positive-mass` | 90.21 | 88.22 | 89.02 | 86/156 |
| + rescue | `--complete-chain-rescue-reads 2` | 90.73 | 88.29 | 89.29 | 87/156 |
| + election (prosecution) | `--junction-election-mode require` | 90.72 | 90.39 | 90.38 | 105/156 |
| + election (**defense**) | `--junction-election-burden defense` | 90.02 | 94.03 | **91.87** | **112/156** |

- defense vs off: +6.29 (CI +5.54..+7.08), 146 up / **0 down**.
- defense vs rescue: +2.584 (CI +2.199..+2.997), 132 up / 8 down (worst −0.06).
- Official contract moved 58.26 -> 62.91 alongside.
- Remaining 44 loss cells: ESPRESSO 30 (mean gap 2.39), StringTie3 14 (3.43),
  concentrated in `c_trunc*` / `c_merge10` RECALL (missing true chains), not
  precision.
- `--containment-evidence-mode` (hazard/TES containment layer) is implemented
  and audit-validated but net-neutral-to-slightly-negative on SIRV
  (no 5'-degradation isoform tail); keep `audit` for real dRNA, do not use
  `require` on SIRV-like data.

Recommended SIRV command (all opt-in flags):
```bash
--profile sirv --no-gpu \
  --min-gtf-abundance 0 --no-floor-gtf-abundance \
  --require-gtf-positive-mass \
  --complete-chain-rescue-reads 2 \
  --junction-election-mode require \
  --junction-election-burden defense   # aggressive; prosecution = conservative
```

## 3. What each lever is (plain language) and where it lives

Pipeline: discovery -> assignment/EM -> interval selection gates -> global
selection gates -> refit -> writers.

1. **P0 strict-positive GTF gate** (`require_gtf_positive_mass`): a GTF
   transcript survives iff EM soft-mass > 0 (zero-mass catalog entries are
   unexpressed pass-through). Enforced at the A4 floor
   (`fin/pipeline/selection.py::select_global`) and again post-refit
   (`strip_zero_mass_gtf_survivors` / `enforce_gtf_positive_after_refit`,
   called from `fin/pipeline/runner.py`; diagnostics fields in
   `abundance_refit.json`). Root cause: the SIRV profile raised the GTF floor
   to 3 while 46/47 `full` FNs had est_count < 3.
2. **Complete-chain rescue** (`complete_chain_rescue_reads`, K=2): a
   candidate whose EVERY junction is attested by >= K exact complete-chain
   reads with indel-clean flanks (per-junction 10 bp; metric =
   `complete_chain_weakest_indel_clean_depth`, per-junction array
   `complete_chain_junction_depths`, computed in
   `fin/candidates/discovery.py::junction_clean_depths` only when evidence
   collection is on) is EXEMPT from the m2 support gate and soft-mass filter
   (policy: `fin/pipeline/selection.py::complete_chain_rescue_ok` — gtf >= K,
   novel multi-junction >= K+1, novel single-junction >= 10). Root cause:
   endpoint-truncated but chain-correct GTF entries lose every per-read
   alignment-score contest; corrupt merged annotation inflates true chains'
   soft/hard ratio.
3. **Containment evidence layer** (`containment_evidence_mode off|audit|require`,
   `fin/analysis/containment_evidence.py`): for short-in-long candidate pairs
   (sub-chain / equal-chain interior-start / mono-span rungs), run the
   degradation-aware ladder from `fin/analysis/tss_evidence.py`
   (TES cluster test = reliable dRNA end; conditional termination-hazard test
   for internal TSS; parent must be independently supported). Tri-state
   verdicts; hazard actions only after the global grouped-FDR pass
   (`finalize_containment` in runner). Artifact `containment_evidence.tsv`.
   SIRV verdict: works, conservative, not a SIRV leaderboard move.
4. **Signal-ballot junction election** (`junction_election_mode off|audit|require`
   + `junction_election_burden prosecution|defense` +
   `junction_election_ballot signal|m1|screen|shuffle` (ablation-only),
   `fin/analysis/junction_election.py`): K-way vote over nearby junction
   variants. Nomination = cigar screen (accused clean <= 1, competitor within
   ±15 bp with clean >= max(3, 3x); junction support map collected in
   discovery: `CandidateSet.junction_support`). Ballots = per-read krill
   eventalign mean NLL against each nominee's representative over merged
   ±8 bp windows (reuses `read_cand_mean_nll` from
   `fin/scoring/m2_junction_nll.py`; window purity + >= 3 events/hypothesis +
   `m2_tiebreak_margin`). Verdicts: quorum 5, supermajority 0.75.
   Burden: prosecution = drop only on conviction; defense = nominated
   junctions drop unless EXONERATED (upheld). Measured: decided elections
   97.6% correct vs truth, 45% abstention; screen-only ablation 91.80 F1,
   signal-prosecution 90.38, defense 91.87 (73 genuine junctions exonerated
   beyond the screen); m1 ≈ signal on SIRV; shuffle ≈ baseline. Artifact
   `junction_election.tsv` (worker shards merged deterministically).

All four are wired through `fin/pipeline/config.py` (fields + loud
validation), `fin/cli.py` (flags), and covered by unit tests
(`tests/unit/test_require_gtf_positive_mass.py`,
`test_complete_chain_rescue.py`, `test_containment_evidence.py`,
`test_junction_election.py`). Suite: **1269 passed, 4 skipped** + M3
regression 5 passed.

## 4. Benchmark tooling (reusable)

- Grid runner: `experiments/prod_validation/sirv4_grid_runner.py`
  (`--samples/--conds/--extra/--extra-with-gtf/--tag`; ABSOLUTE `--out/--tsv`
  paths required — docker gffcompare mounts /SSD only; runs the SIF at
  `experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif` with
  `PYTHONPATH=/SSD/logan/dev/pyfin`, i.e. live source).
- Expression rescore: `experiments/prod_validation/sirv4_expression_rescore.py`
  (builds TPM>0 truth GTFs, rescoress every mode in its `pyfin_modes` dict —
  ADD NEW RUN DIRS THERE; writes
  `experiments/prod_validation/_sirv4grid/expression_score/{scores.tsv,SUMMARY.md,cells/}`).
- Official analyzer: `experiments/prod_validation/sirv4_grid_analyze.py`
  (`--pyfin LABEL=TSV` repeatable; refreshes competitor coverage; writes
  `_sirv4grid/report/`). Rerun both when the 24 missing competitor cells
  arrive.
- Probes (offline, read-only): `fsm_gate_probe.py`, `trunc_containment_probe.py`,
  `discovery_diff_probe.py` under `experiments/prod_validation/`.
- Key run dirs under `experiments/prod_validation/_sirv4grid/` (untracked,
  never commit): `base` (original off), `p0_strict`, `p2_rescue`,
  `p2_caudit`, `p2_require`, `je_full_audit`, `je_full_require`, `je_screen`,
  `je_m1`, `je_shuffle`, `je_defense` (+ `*_det*` determinism reruns).
- Tests must run inside the SIF **with `--nv`** (cupy GPU tests):
  `singularity exec --nv -B /autofs,/SSD experiments/prod_validation/_img/pyfin_gpu_e268c9b.sif python3 -m pytest tests/unit tests/integration -q`.

## 4b. SECOND DATASET (2026-09-04): Sequin

`.../results/sequin_v1_full_sweep_ratio_grid_v2` — 4 samples (Hct116 rep3 run1,
K562 rep4 run1, K562 rep5 run1, MCF7 rep4 run1) x the SAME 26 conditions = 104
cells; 8 competitors already scored per cell. Truth `_ref/_source/annotation.gtf`
(165 tx, chrIS only), genome `refs/sequin/sequin_chrIS.fa`, and — unlike SIRV4 —
an INDEPENDENT concentration ladder `refs/concentration/RNAsequins_MixA_v1.tsv`
(attomol/ul, 6.9 decades). Depth is LOW: 190-575 primary reads on chrIS per
sample; only 28-42 of 165 transcripts are expressed (NanoCount TPM>0) and only
24 have est_count >= 3. A cell runs in ~26 s (the whole grid in ~8 min).

Tooling (one runner, dataset registry — do NOT fork per spike-in):
- `sirv4_grid_runner.py --dataset sirv4|sequin` (`PYFIN_GRID_DATASET` env too),
  `--profile <name>` and `PYFIN_REPO` to score another commit. `stage_refs()`
  copies truth+genome onto /SSD for the dockerised gffcompare.
- `sequin_quant_score.py` — the concentration axis. ALWAYS reports both
  `matched_*` (quality conditional on recovery) and `complete_*` (unrecovered
  truth entered at 0); `matched_*` alone is gameable by emitting only the easy
  high-abundance transcripts.
- `cross_dataset_eval.sbatch` — one code state, both grids, both axes.

**Headline finding (2026-09-04).** Same code, only `--profile` differs, 104 cells:

| operating point | Sn | Pr | matched rho | complete rho | complete NRMSE | recovered/164 |
|---|---|---|---|---|---|---|
| `sirv` (six-H9-sample fit) | 8.85 | 91.33 | 0.727 | 0.465 | 0.327 | 14.6 |
| `real-drna` | 16.05 | 92.35 | 0.718 | 0.585 | 0.278 | 26.5 |

`matched_rho` is UNCHANGED: the SIRV abundance floor never improved the
quantification quality of what was kept — it only deleted transcripts, and
precision did not benefit. The SIRV working point is a dataset fit, not an
existence criterion. Do NOT look for a better constant for Sequin; the
constants must exit the default existence decision — see
`CANDIDATE_EXISTENCE_DESIGN.md`.

## 4b2. Where PyFIN stands against the 8 competitors (2026-09-04)

Tooling: `benchmark_compare.py` (catalog contract, reads the tree's own
gffcompare stats) and `expression_rescore.py --dataset sirv4|sequin` (expression
contract for BOTH grids, scores competitor GTFs from
`<sample>/<cond>/assembly/<tool>.gtf` and any PyFIN run root with the same
gffcompare container).

**Read the expression contract, not the catalog one.** Under the catalog
contract truth is every annotated transcript, so in `full`/`p*` cells the
supplied annotation IS the truth and a tool that merely re-emits it scores
~100/100 without discovering anything — Bambu (Sn 86.8/Pr 95.9 on Sequin) and
TALON do exactly that. The expression contract subsets truth to what the sample
expresses (NanoCount TPM>0: 28-42 of 165 on Sequin), so echoing is paid for in
precision.

SEQUIN, expression contract, 100 fully paired cells:

| tool | Sn | Pr | F1 | top/tied |
|---|---|---|---|---|
| **pyfin (real-drna profile)** | **75.30** | **90.16** | **81.34** | **78/100** |
| isoquant | 61.06 | 85.64 | 69.83 | 11/100 |
| stringtie3 | 64.73 | 77.11 | 69.68 | 2/100 |
| espresso | 58.02 | 89.31 | 69.08 | 4/100 |
| lafite | 53.71 | 83.08 | 63.34 | 6/100 |
| pyfin (sirv profile) | 41.05 | 90.73 | 55.76 | 0/100 |
| flair | 36.28 | 66.21 | 46.38 | 0/100 |
| isotools | 71.06 | 24.02 | 35.76 | 0/100 |
| bambu | 74.79 | 23.28 | 31.91 | 1/100 |
| talon | 88.53 | 16.92 | 27.70 | 0/100 |

The same binary decides first place: PyFIN is #1 by a wide margin at the
real-dRNA operating point and **6th** under the SIRV-fitted one, on identical
code. Largest margins are the junction-corruption cells (c_jitter50: PyFIN
51.4/73.8 vs isoquant 6.3/7.4, stringtie3 13.0/18.2) — the signal-based
junction machinery is doing real work. Weakest cells are the sparse-annotation
ones (p01-p20), where PyFIN leads on F1 but isoquant holds higher precision.

SIRV4, expression contract, 138 fully paired cells (for reference):
pyfin+defense 91.84/94.58 (94/138 top), pyfin_off 84.25/88.39, espresso
80.47/94.44 (34/138), stringtie3 79.73/89.33, lafite 66.84/79.49, isoquant
59.83/89.58, flair 57.51/80.74, bambu 78.51/46.52, talon 96.69/19.40,
isotools 81.27/18.95.

### 4b3. Refreshed against the promoted default (2026-09-07)

Same contract and scorer, PyFIN = the current production default (`CLUST` /
`GCCLUST`), cells where all 9 tools scored:

| tool | SIRV4 (150) Sn/Pr/F1 | Sequin (100) | HEYA8 (125) |
|---|---|---|---|
| **pyfin** | **90.5 / 87.8 / 89.0** | **74.6 / 91.2 / 81.4** | **78.7 / 87.1 / 82.4** |
| espresso | 75.0 / 94.6 / 80.7 | 58.0 / 89.3 / 69.1 | 64.1 / 93.6 / 73.8 |
| stringtie3 | 75.5 / 86.7 / 80.0 | 64.7 / 77.1 / 69.7 | 60.2 / 77.7 / 67.1 |
| isoquant | 55.7 / 89.6 / 66.4 | 61.1 / 85.6 / 69.8 | 57.0 / 85.4 / 66.1 |
| lafite | 61.5 / 73.1 / 65.7 | 53.7 / 83.1 / 63.3 | 56.2 / 91.4 / 66.9 |
| flair | 53.7 / 79.8 / 62.1 | 36.3 / 66.2 / 46.4 | 42.5 / 63.3 / 49.5 |
| bambu | 74.8 / 48.5 / 55.6 | 74.8 / 23.3 / 31.9 | 72.9 / 39.2 / 48.0 |
| talon | 95.6 / 19.6 / 32.4 | 88.5 / 16.9 / 27.7 | 90.6 / 9.0 / 16.3 |
| isotools | 81.2 / 18.9 / 30.4 | 71.1 / 24.0 / 35.8 | 76.3 / 7.7 / 13.9 |

**The whole margin is RECALL, and it is bought where the annotation is poor.**
espresso holds 95-99% precision nearly everywhere against PyFIN's 80-92%, and
loses because its recall collapses exactly where a discovery tool is needed
(SIRV4 p01: Sn 5.1 vs 70.8; c_spurious50: 60.3 vs 96.5; c_flip50: 41.6 vs 81.9).
Read by family (mean F1):

| family | SIRV4 pyfin/espresso | Sequin | HEYA8 |
|---|---|---|---|
| ratio p00..p99 | 84.7 / 65.2 | 73.4 / 57.4 | 77.0 / 61.5 |
| corrupted c_* | 90.8 / 87.4 | 84.6 / 73.7 | 84.5 / 78.8 |
| full | 94.9 / **98.4** | 93.5 / 88.2 | 92.6 / 91.1 |

espresso still wins 68 of 150 SIRV4 cells (2/100 Sequin, 19/125 HEYA8) -- it takes
the clean and near-clean ones narrowly (full -3.5, p99 -3.4, c_merge50 -6.5,
c_merge10 -5.1) while PyFIN takes the degraded ones by 20-60 points. **So the
remaining gap to the best rival is PyFIN's PRECISION at high annotation
completeness (~91 vs ~98 on SIRV4 full/p99)**, which is the wobble/extent-shadow
FP population already characterised, not a recall problem.

### 4b4. The human head-to-head, refreshed (2026-09-08)

The GENCODE sweep has been filling in: **541 competitor assemblies** now exist
where there were 4 (stringtie3 182, bambu 174, isotools 173, espresso 8,
isoquant 4, spread over 15 samples). On the one sample whose reads and signal are
staged locally, **6 of our 7 conditions now have at least one competitor**, giving
13 head-to-head cells instead of 2. Same caveats as before: the competitors
assembled genome-wide and their GTF is restricted to chr1/19/21 post hoc, and
the pyfin column is the current production default (`GCSMC`).

Sn / Pr / F1 per cell:

| cond | pyfin | stringtie3 | espresso | isoquant | bambu | isotools |
| --- | --- | --- | --- | --- | --- | --- |
| full | 31.2/67.9/**42.8** | - | - | - | 100.0/18.3/30.9 | 27.3/14.2/18.7 |
| p50 | 24.4/62.6/**35.1** | 22.9/59.1/33.0 | - | - | 55.6/19.6/29.0 | 27.3/14.2/18.7 |
| p00 | 15.8/56.2/24.7 | 16.0/56.8/**25.0** | - | - | 13.8/62.2/22.6 | - |
| c_merge50 | 27.1/61.3/37.6 | 28.0/57.6/**37.7** | - | - | - | - |
| c_trunc3p50 | 30.9/66.5/**42.2** | - | - | 15.4/74.2/25.5 | - | 27.3/14.2/18.7 |
| c_trunc5p50 | 30.5/66.8/41.9 | - | 29.4/74.2/**42.1** | - | 99.8/18.2/30.8 | - |
| c_jitter50_10bp | 16.3/41.9/**23.5** | - | - | - | - | 27.1/14.1/18.5 |

Paired means over the cells each tool actually covers (**different subsets, so the
tools are not comparable to each other here — only each one to pyfin**):

| tool | n | tool Sn/Pr/F1 | pyfin Sn/Pr/F1 | dF1 | pyfin wins |
| --- | --- | --- | --- | --- | --- |
| stringtie3 | 3 | 22.3/57.8/31.9 | 22.4/60.0/32.5 | **+0.57** | **1/3** |
| espresso | 1 | 29.4/74.2/42.1 | 30.5/66.8/41.9 | **-0.23** | **0/1** |
| isoquant | 1 | 15.4/74.2/25.5 | 30.9/66.5/42.2 | +16.69 | 1/1 |
| bambu | 4 | 67.3/29.6/28.3 | 25.5/63.4/36.1 | +7.78 | 4/4 |
| isotools | 4 | 27.2/14.2/18.6 | 25.7/59.7/35.9 | +17.23 | 4/4 |

**The headline is that the spike-in dominance does not transfer, and now there
are two peers rather than one.** stringtie3 is level (it takes p00 by 0.30 and
c_merge50 by 0.10, we take p50 by 2.10) and espresso is level. The large margins
are all against tools that fail in a way the expression contract is designed to
expose: bambu re-emits the annotation (`full`: Sn 100.0, Pr 18.3, 40209 query
models against our 3377), and isotools emits ~14k models on three chromosomes
almost regardless of the input annotation -- its GTF barely moves across
conditions (661124 / 661300 / 661135 / 663365 lines) so its score is near
constant at 27.3/14.2/18.7. Neither margin says anything about the axis that
matters.

Against the two peers the shape is the same as on the spike-ins: **espresso buys
precision (74.2 vs our 66.8) and pays recall; we are the opposite.** stringtie3
sits between us on both. Our remaining work on real human data is precision, and
the competitive claim to make is the one about degraded annotations, not a
blanket lead.

*(Harness: `experiments/prod_validation/human_headtohead.py`, run on the HOST --
the scorer shells out to docker.)*

### 4b5. The human head-to-head, complete for the three sweep tools (2026-09-09)

The sweep finished for stringtie3 / bambu / isotools: **1184 assemblies, 390
each = 15 samples x 26 conditions** (espresso 10, isoquant 4). Every one of our
7 production conditions on H9 rep2 run2 now has all three, so the head-to-head
is **23 cells** (same harness, same caveats as 4b4: genome-wide assembly
restricted to chr1/19/21 post hoc; pyfin = `GCSMC`).

| cond | pyfin | stringtie3 | espresso | isoquant | bambu | isotools |
| --- | --- | --- | --- | --- | --- | --- |
| full | 31.2/67.9/42.8 | 33.9/65.5/**44.7** | - | - | 100.0/18.3/30.9 | 27.3/14.2/18.7 |
| p50 | 24.4/62.6/**35.1** | 22.9/59.1/33.0 | - | - | 55.6/19.6/29.0 | 27.3/14.2/18.7 |
| p00 | 15.8/56.2/24.7 | 16.0/56.8/**25.0** | - | - | 13.8/62.2/22.6 | 27.2/14.2/18.7 |
| c_merge50 | 27.1/61.3/37.6 | 28.0/57.6/**37.7** | - | - | 72.0/13.6/22.9 | 27.3/14.2/18.7 |
| c_trunc3p50 | 30.9/66.5/42.2 | 32.7/63.1/**43.1** | - | 15.4/74.2/25.5 | 99.8/18.2/30.8 | 27.3/14.2/18.7 |
| c_trunc5p50 | 30.5/66.8/41.9 | 33.5/62.3/**43.6** | 29.4/74.2/42.1 | - | 99.8/18.2/30.8 | 27.3/14.2/18.7 |
| c_jitter50_10bp | 16.3/41.9/**23.5** | 9.5/20.9/13.1 | - | - | 5.4/1.0/1.7 | 27.1/14.1/18.5 |

| tool | n | tool Sn/Pr/F1 | pyfin Sn/Pr/F1 | dF1 | pyfin wins |
| --- | --- | --- | --- | --- | --- |
| stringtie3 | 7 | 25.2/55.0/34.3 | 25.2/60.5/35.4 | +1.09 | **2/7** |
| espresso | 1 | 29.4/74.2/42.1 | 30.5/66.8/41.9 | -0.23 | 0/1 |
| isoquant | 1 | 15.4/74.2/25.5 | 30.9/66.5/42.2 | +16.69 | 1/1 |
| bambu | 7 | 63.8/21.6/24.1 | 25.2/60.5/35.4 | +11.29 | 7/7 |
| isotools | 7 | 27.3/14.2/18.7 | 25.2/60.5/35.4 | +16.72 | 7/7 |

**What the four new stringtie3 cells change:** the +1.09 mean is carried
entirely by `c_jitter50_10bp` (+10.41, stringtie3 collapses to 9.5/20.9 when the
annotation is jittered) and `p50` (+2.10). On the five conditions where the
annotation is intact or merely truncated/merged, **stringtie3 is ahead on every
one, by 0.10 to 1.92 F1**, and the shape is constant: it recalls ~2-3 pp more
(33.9 vs 31.2 on `full`) while we are ~2.5-4.5 pp more precise (67.9 vs 65.5).
Its extra recall is not bought with junk -- its query count is only ~12% above
ours (3804 vs 3377 on `full`). So on real human data with a usable annotation
the honest statement is "level with stringtie3, slightly behind on F1, ahead on
precision, behind on recall"; the clear lead is confined to corrupted (`jitter`)
and half-removed (`p50`) annotations, and disappears at `p00` (-0.30).

bambu and isotools are unchanged in kind: bambu echoes the annotation (`full`
Sn 100.0 / Pr 18.3; at `p00` it drops to 13.8/62.2 and we beat it by 2.08) and
isotools is annotation-blind (27.3/14.2/18.7 in all 7 cells, ~14.1k models).
These 14 "wins" say nothing about the axis that matters.

**Next thing this points at:** the 19 conditions we have not run (p01..p99,
the 10% corruptions, flip/ir/skip/spurious at 50%) are now fully covered by
stringtie3. Running `GCSMC` on them would make the human comparison a full
26-cell curve against the one real peer, and would show where between `p50`
(+2.10) and `full` (-1.92) the crossover sits. Each run is one SLURM job on the
staged chr1/19/21 BAM.

**Real human, the older two-cell snapshot (superseded by 4b4 above):**
Only two overlap conditions we ran, and the competitor assembled genome-wide, so
its GTF is restricted to chr1/19/21 before scoring (same truth, same scorer; the
restriction is post-hoc, so treat these two cells as indicative, not as a
benchmark):

| cell | tool | Sn | Pr | F1 |
|---|---|---|---|---|
| c_trunc3p50 | isoquant | 15.40 | 74.20 | 25.51 |
| | **pyfin** | **30.90** | 66.40 | **42.17** |
| | pyfin +gtf-exempt | 31.50 | 66.50 | 42.75 |
| c_trunc5p50 | espresso | 29.40 | **74.20** | 42.11 |
| | pyfin | 30.50 | 66.60 | 41.84 |
| | pyfin +gtf-exempt | **31.30** | 67.00 | **42.67** |

On a real transcriptome espresso is a PEER, not 8 points behind: it trades 1.1
points of recall for 7.6 of precision and ties on F1. The spike-in dominance does
not transfer, and the same precision gap is what decides it. (Harness:
`/tmp/human_headtohead.py` pattern -- filter to the run's chromosomes, then reuse
`expression_rescore.score`; run it on the HOST, the scorer shells out to docker.)

## 4c. The main line now: a unified candidate EXISTENCE test

`CANDIDATE_EXISTENCE_DESIGN.md` is the frozen pre-registration;
`fin/analysis/candidate_existence.py` is the audit-only v1
(`--candidate-existence-mode audit`, default off). Codex designs/reviews:
`.omc/codex_consult_principled_sequin_2026-09-04.md`,
`.omc/codex_review_existence_design_2026-09-04.md`.

Criterion: a candidate exists iff adding it gives a REPRODUCIBLE predictive
improvement, on reads HELD OUT of the fit, for evidence its conflict-component
alternatives cannot explain. Three-state (`supported` / `unresolved` /
`contradicted`); `p > alpha` is NOT evidence of absence, so v1 never drops
anything. v1 tests only externally specified annotated multi-exon candidates
(selection-honest by construction — discovery used every read, so novel
candidates need fold-specific re-discovery or an honest split) and only the
structural channel (M1 tie mask + M2 NLL + diff-region coverage).

First audit (2026-09-04), both datasets. The first pass found TWO real bugs in
the audit itself, which is precisely what audit-first is for:
1. it scored ``d_tx``, but with the diff-coverage gate on (the default) that is a
   DECISION matrix -- a distinguishing read is one-hot and a fuzzy read carries the
   component-wide covered_vote ratio -- so the strongest evidence reads were being
   SKIPPED. Fixed by adding ``QuantOutput.structural_d_tx_raw`` (raw per-read
   junction NLLs over the raw AS tie set, built before the gate) and reading that.
2. a global hash-mod-5 fold left most low-depth components with empty folds.
   Fixed by component-balanced adaptive folds, ``K = min(5, n_component_reads)``,
   degenerating to leave-one-out; ``K < 2`` abstains.
Sequin identifiable rows went 3 -> 28 as a result.

- byte-inert: Sequin 37/37 and SIRV4 54/54 outputs IDENTICAL, 0 DIFF; determinism
  rerun of ``candidate_existence.tsv`` byte-identical.
- discriminative power at 50% corruption (truth = does the SUPPLIED annotation's
  chain match the reference chain; evaluation-only, never an input):

| dataset | true_chain | corrupted_chain | AUC(true>corrupted) |
|---|---|---|---|
| SIRV4 (36 cells) | n=282, median 4.47, 261 positive | n=19, median 0.36, 13 positive | 0.625 |
| Sequin (24 cells) | n=43, median 16.46, 32 positive | n=2, median -0.00, 0 positive | 0.767 |

- The informative breakdown is per condition, not the pooled AUC. On SIRV,
  `c_jitter50` corrupted chains get median 0.00 (3/9 positive) -- the structural
  channel DOES reject junction-level corruption. But `c_merge50` corrupted chains
  get median 14.93 with **10/10 positive**: a merged annotation genuinely explains
  the junction evidence of both merged transcripts, so the structural channel
  alone cannot reject it. That is a channel gap, not a tuning problem -- merge/
  trunc3p corruption needs the 3'/TES+polyA channel, which is exactly why the
  channels are being added one at a time.

Gates this is meant to replace (only after the audit proves calibration,
determinism and selection honesty): m2 support gate + complete-chain rescue,
abundance floor, isoform fraction, soft/hard ratio, fulllen fraction.

## 4d. The production DEFAULT is now chosen by cross-dataset maximin (2026-09-05)

The user's requirement: the shipped default must be the best configuration for
the data, with no per-dataset tuning. That rules out "pick the profile that suits
your spike-in". The selection rule (`pick_overall_default.py`) is therefore
MAXIMIN on the expression contract: rank by the WORSE of the two datasets, never
by a pooled average over 260 cells whose designs and depths differ.

**Profile.** The CLI default is `real-drna`; the `sirv` profile is a benchmark
working point and must not be recommended. On the expression contract, levers
off: SIRV 86.98 vs 85.67, Sequin 80.37 vs 55.38. `real-drna` wins on BOTH.

**Levers (P0 / complete-chain rescue / junction election) stay OPT-IN.** They are
a large precision win on SIRV (Pr 85.79 -> 92.64 at nearly no recall cost) and a
net loss on the held-out set (Sequin F1 -1.86 to -4.76). The 2x2:

| arm | SIRV Sn/Pr/F1 | Sequin Sn/Pr/F1 | worst F1 |
|---|---|---|---|
| off + real-drna | 88.57/85.79/86.98 | 74.04/89.90/80.37 | **80.37** |
| best + real-drna | 88.17/92.64/90.20 | 69.91/91.27/78.50 | 78.50 |
| best + sirv | 90.40/93.73/91.91 | 68.66/87.62/75.61 | 75.61 |
| off + sirv | 83.20/88.47/85.67 | 40.72/90.36/55.38 | 55.38 |

Decomposition showed the levers' SIRV gain is largely a REPAIR of the `sirv`
profile's own GTF floor (sirv+levers 15.33 Sn vs sirv alone 8.85 on Sequin), and
that their Sequin cost is depth-driven (`require_gtf_positive_mass` and the
election quorum are statements about sampling, not existence). Both prosecution
and defense burdens behave the same at Sequin depth, so the burden is not the
cause. Removing the read-count floor entirely is worse on both (Sn AND Pr drop:
the relative gates' denominators change).

**One default DID change: annotation as a fold container, gated by
contradiction.** A zero-read annotated chain may absorb a read sub-chain only
when the reads do not contradict it (no different read junction within the
alignment tolerance of either of its splice sites). Measured, expression
contract, paired bootstrap:

| dataset | dSn | dPr | dF1 | 95% CI |
|---|---|---|---|---|
| SIRV4 (development) | -0.08 | -0.12 | -0.10 | [-0.16, -0.05] |
| Sequin (held out) | +2.46 | -0.89 | +1.11 | [+0.71, +1.53] |

maximin worst-dataset F1 79.62 -> 80.78. The ungated version scores higher still
on Sequin (81.40) but pays for it exactly where the annotation is corrupted
(c_merge50 -1.94 F1, Pr -7.75); the gate turns that into +0.58 while keeping the
trustworthy-annotation gain (c_trunc3p50 +3.38). No threshold was introduced: the
window is the project's existing junction-equivalence tolerance.

METHOD NOTE (cost me several wrong conclusions): `expression_rescore.py` caches
cells by (sample, cond, TOOL LABEL). Reusing a label for a different run returns
the OLD numbers silently. Use unique labels or `--force`.

## 4e. Codex review of the merged batch (2026-09-05)

`.omc/codex_review_production_batch_2026-09-05.md` — verdict "mergeable after
fixes", 3 P1 + 3 P2, all fixed and pinned by tests before the commit:

- P1 the contradiction registry used the MAPQ>=1 SUPPORT counts while read chains
  enter clustering unconditionally, so a low-MAPQ read could shape the topology
  yet count as "the reads are silent". Now built from every read chain that
  enters clustering (`discovery.py`, `observed_junctions=`).
- P1 `tss_shadow_recovery` mutated an existing candidate's support in the
  existing-id branch, so `audit` was not output-inert. The merge is now returned
  in `ShadowRecoveryResult.support_merges` and applied by the runner in `require`
  mode only.
- P1 the `m1_em` path only used `hard_assign` when `krill_tiebreak` was on, so
  the default still broke exact ties by column index. Now unconditional.
- P2 the election's mappy LRU was keyed by candidate id alone (ids recur across
  intervals) -> keyed by id + sequence digest.
- P2 `candidate_existence_mode` / `folded_tss_recovery_mode` silently no-oped
  outside `m2_em` -> now a loud validation error.
- P2 `junction_support.tsv` was not written for an empty registry, so "no
  junctions" was indistinguishable from "the step never ran" -> header-only file.

Final validation of the exact committed tree: 1370 tests pass; determinism
32/32 and 48/48 byte-identical across different worker counts; mass balance
16/16 within 1e-9; both grids rerun (below).

| dataset | HEAD (ca82319) Sn/Pr/F1 | committed Sn/Pr/F1 |
|---|---|---|
| Sequin (held out, 100 paired cells) | 71.42/91.44/79.62 | **73.98/90.52/80.78** |
| SIRV4 (development, 150 paired cells) | 89.29/85.96/87.44 | 89.21/85.84/87.33 |

## 4f. Endpoint channels + read-extent hypothesis (2026-09-05, branch `existence-algorithm`, UNCOMMITTED)

What: `--assignment-channels tes|ends` (requires `--signal` and m2_em). Per-sample nuisance
kernels (3' EndpointModel on polyA-PASS ends, one-sided 5' StartModel) enter the m2_em
posterior as 3'/5' terms, and an annotated model whose ends the reads contradict (a
polyA-PASS read past the claimed TES / a read starting upstream of the claimed TSS, not
explained by any other annotated model on the same or a wobble chain) gets a READ-EXTENT
SIBLING `<gtf_id>.ext` (same intron chain, `source="gtf"`, `extent_of=gtf_id`, only the
contradicted end re-estimated from the contradicting reads, outward). Selection keeps one of
the pair by posterior mass (ties -> annotation), folding the other; a pair split across
intervals is unified before aggregation (`unify_extent_pairs`) and its fold is a refit
redirect. A kernel is usable only when calibrated (n_fit >= 50, narrow scale <= 25 bp);
otherwise every endpoint term abstains. Design + measured story: CANDIDATE_EXISTENCE_DESIGN.md
SS9. Root cause it fixes: trunc-condition wobble FPs (HEYA8 SIRV107 reproduction in
`experiments/prod_validation/_heya8grid/_dbg_sirv107/`).

Measured (expression contract, paired vs the maximin default, snapshot `_snap/endsarm5`):
Sequin F1 +0.64 (CI [+0.27,+1.12], Sn +0.66, Pr +0.64, 17 up / 0 down); SIRV F1 +1.49
(CI [+0.98,+2.05], Sn +1.12, Pr +1.76, 95/26/35); HEYA8 pending. `folded_tss_recovery_mode=
require`: SIRV +0.05, Sequin 0, HEYA8 +0.03 (CI [+0.01,+0.07]).

Codex (gpt-6-astra, effort max) review rounds 1-9: every finding fixed except one it could
not clear on a scope argument: the aggregator's cross-interval dedup ratio (unique hard
reads / hard observations) is an approximation ONLY for reads that reach two intervals'
assignment axes, which in this pipeline are exactly the fusion-classified primaries
(>= 250 bp soft clip) excluded from interval construction but not from the read axis --
the PRE-EXISTING "fusion quant-axis leak" (discovery.py ~856/866; deferred in memory since
the SIRV4 era). This batch neither introduces nor worsens it; it is tracked as its own
default-path change requiring its own grid validation.

## 4g. The production DEFAULT now carries the endpoint channels (2026-09-06)

`PIPELINE_PROFILES["real-drna"]` (the CLI default; `real-drna-precision` inherits)
resolves `assignment_channels="ends"` and `folded_tss_recovery_mode="require"`.
The `sirv` benchmark profile and the raw dataclass defaults are unchanged, and an
explicit `--assignment-channels none` / `--folded-tss-recovery-mode off` still wins
(recorded in `profile_overrides`). Without `--signal` the 3'/polyA term abstains;
the 5' term still applies when a `--gtf` is supplied (it scores read starts against
the annotated TSS), so an unguided run without signal has no endpoint evidence and
says so. Under another quantifier a PROFILE-IMPLIED setting is disabled with a
warning (the `post_selection_refit` precedent); an explicit one is an error.

Measured (expression contract, paired vs the previous default, seed 20260831,
each arm one `source_sha256`, snapshot verified through `run_manifest.source_root`):

| grid | previous default | new default | dF1 (95% CI) | up/down |
| --- | --- | --- | --- | --- |
| Sequin (104) | 73.98/90.52/80.78 | 74.67/91.19/**81.44** | +0.64 [+0.27,+1.12] | 17/0 |
| SIRV4 (156) | 89.21/85.84/87.33 | 90.47/87.80/**88.99** | +1.59 [+1.09,+2.15] | 94/19 |
| HEYA8 (130) | 77.39/85.18/80.76 | 78.75/87.08/**82.43** | +1.63 [+1.04,+2.29] | 69/4 |

Sn AND Pr both rise on every grid; best competitor is isoquant 69.83 (Sequin),
espresso 80.74 (SIRV4), espresso 73.78 (HEYA8). The promoted profile reproduces
the separately measured `--assignment-channels ends` arm (Sequin 104/104 cells
identical; SIRV4 +0.05 and HEYA8 +0.04, the folded-TSS recovery the ends arm did
not carry), and the sample-evidence rework that followed the Codex rounds is
result-neutral (Sequin rerun on the final code: 104/104 identical).

**Sample-evidence artifact (the endpoint channels' input).** One krill polyA pass
plus the strand-aware read ends per sample, in the work dir. Contract, after 11
Codex rounds: the file is ADDRESSED BY ITS INPUT IDENTITY (path/size/mtime/ctime +
head-tail digest of BAM and signal, pore, requested GPU) and is IMMUTABLE -- every
write is a create (`os.link`), never a replace, so no run can change an artifact
another run is reading. A result that is incomplete (`partial`, some calls
unusable) or absent (`unavailable`, the caller failed) is written to an
attempt-unique run-local name, never shared, and removed at `cleanup()`. Payloads
carry a completeness manifest (`n_ends`/`n_polya`) so a truncated-but-parseable
file is rejected rather than read as a sample with fewer reads, and every record is
type/range validated on the way in. Inputs are re-fingerprinted after collection
(and the annotation after parsing) so a file replaced mid-run raises instead of
publishing a mixed snapshot.

**Known bound (class b, accepted):** the fast fingerprint samples the file's edges,
so an interior edit that also preserves size, both timestamps and the sampled bytes
is invisible. `PYFIN_EVIDENCE_STRICT_HASH=1` digests the whole file for setups where
inputs may be rewritten in place behind the pipeline.

## 4h. Human GENCODE evaluation of the promoted default (2026-09-06, IN PROGRESS)

The fourth grid, `results/gencode_full_sweep_ratio_grid_v2`, is the real human
transcriptome: GENCODE v44 on GRCh38, 15 samples x 26 conditions (the same
corruption ladder), NanoCount truth per sample (47667 of 252894 transcripts
expressed). Its competitor runs DO NOT EXIST -- the workflow could not resolve
its container engine (`missing-containers.md`: `apptainer` not found), so only
4 assembler GTFs exist across 390 cells. A PRIOR human run on a different grid
(rates 20, 6 tools, `analysis/gencode_full_sweep/scoring/`) is scored on the
CATALOG contract and shows what that contract rewards: bambu 99.7 F1 and talon
85.4 with **91.6% of their matched transcripts carrying zero reads**, and F1
collapsing to 5.6 / 9.3 without an annotation. PyFIN was excluded from it.

Only 3 of the 15 samples have a genome alignment; one (H9 rep2 run2) is fully
staged locally. Runs use a chr1+chr19+chr21 slice (49195 of 296508 reads, 40163
of 252894 annotated transcripts, 7398 expressed): both arms see identical
intervals, so the comparison is paired; only the absolute scale differs.

**Cost of the promoted default on real data:** the krill polyA prepass is ~51
min for 296508 reads (CPU, 6 threads), once per sample and then cached.

**Expression contract, 7 conditions, paired:**

| arm | Sn | Pr | F1 |
| --- | --- | --- | --- |
| pre-promotion default | 25.10 | 60.13 | 35.26 |
| `tes` (3' only) | 25.10 | 59.54 | 35.15 |
| promoted default (`ends`) | 25.11 | 59.97 | 35.25 |
| ECHO (emit the input annotation) | 60.89 | 12.50 | 20.53 |

The endpoint channels are **structurally neutral on human data** (mean dF1
-0.01, 3 up / 4 down) -- and they must be, because gffcompare matches INTRON
CHAINS: two models with the same chain and 500 bp different ends are both `=`.
The structural F1 used everywhere in this benchmark is blind to the extent,
which is the only thing these channels change. The echo baseline is the honest
reference for the human grid: emitting the annotation gets 60.9 recall at 12.5
precision.

**The axis the channels actually target** (`endpoint_accuracy.py`: for every
gffcompare `=` match, the strand-aware distance from the emitted ends to the
truth's ends):

| arm | TES within 20 bp | TSS within 20 bp | TES within 20 bp on c_trunc3p50 |
| --- | --- | --- | --- |
| pre-promotion default | 86.8% | 77.6% | 61.7% |
| `tes` (3' only) | 85.6% | 77.8% | **70.8%** |
| promoted default (`ends`) | 85.7% | **68.6%** | 70.9% |

The 3' repair does exactly what it claims: where the annotation's 3' ends were
truncated it recovers +9 points of TES agreement, on real human data. The 5'
start term is a REGRESSION: -9 points of TSS agreement overall (96.4% -> 79.1%
on the uncorrupted annotation), because on real data reads starting upstream of
a CORRECT annotated TSS are common (unannotated upstream starts, first-intron
retention, run-through) and the 5' contradiction test has no analogue of the
polyA-PASS requirement that makes the 3' test specific.

**Resolved.** `tes` (3' only) is worse than `ends` on every spike-in grid
(+0.36 / +1.04 / +1.01 vs +0.64 / +1.59 / +1.63 over the pre-promotion
default), so the 5' term earns about half the gain and must not be removed. Its
specificity comes from a START CLUSTER instead: a 5' start may contradict an
annotated TSS only when >= `tss_contradiction_min_reads` (3) molecules end
within 10 bp of it (`sample_nuisance.TSS_CLUSTER_TOL_BP`). Human: TSS-within-20bp
68.6% -> 75.8% overall, 79.1% -> 90.8% on the uncorrupted annotation, the 3'
repair intact (trunc3p TES 70.9%) and the genuine 5' repair IMPROVED (trunc5p
TSS 60.0% -> 62.5%). Spike-ins: SIRV cell-identical to plain `ends` (156/156,
+1.59 kept), Sequin +0.59 vs +0.64. Landed as the default.

**The raw signal cannot substitute for it.** krill's polyA sample boundary plus
`read_rate` give a molecule length independent of the alignment; the reads that
trigger falsely have a MEDIAN of 1 base of signal beyond their basecall (probe:
`signal_completeness_probe.py`). Those molecules really do end where the read
ends -- 5' truncation is chemical and happens before sequencing, so unlike the
3' end (where a polyA PASS is a positive per-molecule statement) there is no
signal-level statement to be had about a 5' end.

## 4i. What actually costs human recall: the support gate (2026-09-07)

Human recall decomposition (`human_recall_decomposition.py`, H9 chr1/19/21,
7398 expressed): 63.8% of expressed transcripts have NO read carrying their
full intron chain, and 25.4% are supported only by 3' SUFFIX fragments. We emit
30.7%, ABOVE the 29.2% full-chain ceiling, and only 3.4% is "attested by >= 2
full-chain reads but not emitted". Of the fragment-supported 1881, only 14.6%
ever become a candidate -- but 92% of those survive, so the bottleneck was the
interval gates, not selection.

The gate responsible is `m2_support_gate`: it asks a read to DISCRIMINATE a
candidate (its sole best-AS, or the M2 winner), which a 5'-degraded fragment
shared by several annotated isoforms of one gene cannot do. Turning it off on
human data is +7.7 F1 on the uncorrupted annotation with BOTH axes up (Sn
31.2 -> 39.5, Pr 67.8 -> 69.6) and costs ~1 F1 under jitter -- the corruption it
was built for. `--m2-support-gate-gtf-exempt` keeps both: an ANNOTATED
multi-exon candidate is exempt when the reads CONTRADICT none of its junctions,
reusing clustering's own predicate (a jittered/merged junction sits beside an
observed one; an unobserved junction is silence, not disagreement).

Human, expression contract, 7 conditions (pre-promotion default -> full stack):

| cond | old default | + cluster gate | + gtf exemption | dF1 |
| --- | --- | --- | --- | --- |
| full | 43.06 | 42.73 | **46.24** | +3.18 |
| p50 | 34.73 | 35.05 | **37.07** | +2.33 |
| c_trunc3p50 | 41.48 | 42.17 | **44.86** | +3.38 |
| c_trunc5p50 | 41.10 | 41.83 | **45.39** | +4.29 |
| c_merge50 | 38.02 | 37.56 | **40.17** | +2.15 |
| c_jitter50_10bp | 23.77 | 23.40 | 23.34 | -0.42 |
| p00 | 24.59 | 24.58 | 24.58 | -0.01 |
| **overall** | 25.10/60.13/35.26 | 25.17/60.16/35.34 | **27.33/60.37/37.38** | **+2.13** |

**The exemption stays DEFAULT-OFF: it fails the cross-dataset rule.** Three
forms were measured against the cluster-gated default, each stricter than the
last (paired, expression contract):

| form | human | Sequin | SIRV4 | HEYA8 |
| --- | --- | --- | --- | --- |
| any uncontradicted annotation | +2.04 | -0.27 | -0.23 | - |
| + containment-only (no wobble twin) | +1.48 | -0.27 | -0.07 | - |
| + EM soft mass >= 2 | **+0.49** (6 up / 0 down) | -0.02 (ns) | +0.03 (ns) | **-0.20** |

The mass requirement is what restored the corruption protection -- with it, the
human jitter condition turns POSITIVE (+0.09) where the unrestricted exemption
cost 0.06 -- and it removes the Sequin/SIRV4 cost entirely. HEYA8 still loses
0.20 F1 (CI [-0.32,-0.09], 20 cells down), so `--m2-support-gate-gtf-exempt` is
an option, not a default: worth +0.49 F1 on a real human transcriptome and
0.20 on the dual spike-in.

What the difference between the domains actually is, since it is not a defect to
be fixed by a better predicate: the exemption admits annotated isoforms that no
read can single out, so its value is the PRIOR that such an isoform is really
expressed. In a human locus with reads that prior is high; in a spike-in mix,
whose annotation is a defined catalog of which only a subset is present, it is
low. No junction- or signal-level test separates those cases, because the
evidence genuinely does not exist -- which is the same conclusion the existence
design reached for wobble twins (CANDIDATE_EXISTENCE_DESIGN section 8).

## 4j. Where the DRS signal is NOT being used, and one REJECTED attempt to extend it (2026-09-07)

**The measurement (stands).** On the human `full` cell the junction-signal layer
attempts 25686 contrasts and DECIDES 664 of them -- **2.6%**; only 2% of the
interval ties are refined by M2. The abstentions are structural:

* **51.7% `unequal_intron_count`** -- the tied candidates are NESTED (a
  5'-degraded fragment is compatible with both a 3-intron and a 5-intron isoform
  of one gene). `_m2_score_region` restricts the summed-LLR metric to the
  same-intron-count wobble niche.
* **40.8% `no_window`** -- same intron count and no differing junction: the tie
  is not about a junction at all. Those are EXTENT disagreements, which the
  endpoint channels decide, and the signal is correctly silent.

**The attempt, and why it was REMOVED.** `--m2-containment-contrast` scored the
nested contrast and charged each candidate, per event it was missing relative to
the best-covered candidate, the mean per-event NLL of the candidates that did
explain the window. It measured +1.25 F1 on the human uncorrupted annotation and
-0.16 / -0.06 on SIRV4 / Sequin, so it failed the cross-dataset rule -- and the
Codex review (gpt-6-astra) then showed the rationale was wrong, which is the
reason it is gone rather than merely default-off:

* the case it was written for is UNREACHABLE. `_mean_nll_in_window` returns
  `(nan, 0)` for a candidate with no events in the window, and both scoring
  paths keep only candidates with `n_ev > 0` and a finite NLL. A zero-event
  candidate never reaches the correction, so "a candidate that explains nothing
  scores a free 0 and wins" was never what the code did;
* what it actually computed, for two covered candidates, is
  `S'_A - S'_B = [2 n_A n_B / (n_A + n_B)] (mean_A - mean_B)` -- a ranking of
  candidate-specific MEAN NLLs with a count-dependent margin, over event sets
  that need not correspond to the same raw samples;
* so the human gain is not attributable to better signal discrimination. A
  plausible mechanism is that forcing decisions concentrates EM mass and lets
  more candidates clear the abundance floor, which pays where the admitted
  candidates are usually real (human) and costs where they are not (spike-ins) --
  but that mechanism is NOT established either, and `tie_ok=True` rescuing both
  candidates is an equally live explanation.

**What a defensible extension would look like** (Codex, and I agree): compare
likelihoods over a COMMON set of raw signal samples, marginalising over each
candidate's alignment paths, and treat a missing eventalign record as TECHNICAL
missingness -- these candidates are M1-tied, i.e. the read aligns about equally
well to both, so absent events are an anchoring artefact, not evidence of an
inability to generate the sequence. That inability belongs to the
sequence/alignment model. Abstain when comparable signal cannot be recovered for
both hypotheses. Do not re-attempt the shortcut above.

**What the signal DOES buy, and it is the DRS-specific part of this pipeline:**
the 3'/polyA channel, where a PASS is a genuine per-molecule statement (this
molecule was polyadenylated here). That is what the promoted default rests on
(+1.5 F1 on three grids) and what repairs a truncated annotation's 3' end on
real human data (+9 points of TES agreement). The 5' end has no such statement
-- measured, not assumed (section 4h) -- which is why its specificity comes from
a read CLUSTER instead.

## 4k. Widening the DRS channel: the polyA QC classes (2026-09-07, MEASURED -- both stay OFF)

Section 4j ended with the conclusion that the 3'/polyA channel is where the raw
signal actually pays. The measurement that follows from it: **the channel
abstains on 57% of a human sample's reads**, because it admits only
`polya_qc == PASS` (43.2% -- READ_FAILED 33%, ADAPTER 20%, NOREGION 7% of the
49k reads on H9 chr1/19/21). The abstentions are not a statement about the 3'
END; they are a statement about the tail LENGTH, which this channel never reads.

Fraction of reads whose 3' end is within 20 bp of an annotated TES:

| dataset | PASS | ADAPTER | READ_FAILED | NOREGION |
| --- | --- | --- | --- | --- |
| human H9 | 82.8 | 79.9 | 74.2 | **10.1** |
| HEYA8 | 95.3 | 92.1 | 91.3 | 74.5 (75.1 on the largest sample) |
| Sequin | 95.6 | 89.0 | 88.2 | 26.3 |
| SIRV4 | 95.1 | 95.6 | 92.6 | 34.5 |

NOREGION -- krill found no polyA region at all -- is genuinely bad on human
(median 334 bp away) and genuinely fine on HEYA8. **Which classes are worth
trusting is a property of the run**, so a whitelist is the wrong instrument.

**The design under test (`--endpoint-qc-conditional`).** Admit every read krill
returned a verdict for and fit the two-Laplace endpoint mixture SEPARATELY per QC
class; score each read under its own class's kernel; a class whose kernel fails
the existing `available` rule (n >= 50, b_narrow <= 25 bp) abstains. The broad
weight eta is then a per-sample, per-class estimate of "how often this class's
end is not a transcript end". Offline fits, all four grids
(`qc_conditional_kernel_probe.py`):

| dataset | eta PASS | eta ADAPTER | eta READ_FAILED | NOREGION |
| --- | --- | --- | --- | --- |
| human | 0.209 | 0.237 | 0.294 | b_narrow 4407 -> **abstains** |
| HEYA8 | 0.245 | 0.127 | 0.104 | 0.256 (usable, and it is 75.1% good) |
| Sequin | 0.045 | 0.126 | 0.202 | n=19 -> **abstains** |
| SIRV4 | 0.190 | 0.095 | 0.153 | eta pinned at the cap |

The eta ordering reproduces the measured accuracy ordering without being told
it. **Abstention is decided by `available` -- the narrow SCALE and the count, not
eta:** human NOREGION fits b_narrow 4407 bp and Sequin NOREGION has n=19, so both
drop out under either cap. What the eta cap does is separate: the EM capped eta
at 0.5, which is right for a POOLED fit (a sample cannot be more than half
read-through) and wrong per class, because a class can fit a genuinely narrow
component and still be mostly wrong -- SIRV4 NOREGION fits b_narrow 1.8 bp with
only 34.5% of its ends within 20 bp of a TES. Released to 0.99 for per-class fits
only, that class fits eta 0.695 (0.305 narrow weight per read) instead of the
capped 0.5, while **every informative class fits an identical eta to 3 decimals
on all four grids**. So the release calibrates weak classes; it does not create
abstention.

**Scope, forced by a Codex finding.** The conditional mode changes the
ASSIGNMENT channel only. The other consumer of the endpoint kernel -- the
discovery contradiction hook that decides whether a read's 3' end disagrees with
an annotated TES -- tests `narrow_posterior(d) < 0.5`, an ABSOLUTE threshold that
is not invariant to eta: a class fitted as LESS reliable has a lower narrow
posterior at every distance, so it would raise MORE contradictions and be
protected by FEWER annotated alternatives. Less confidence producing more
aggression is backwards, so discovery keeps the pooled kernel and the
verdict-name predicate, and the pooled kernel is fitted from a second, name-based
evidence derivation (`NuisanceEvidence.pooled`) so it is bit-for-bit the kernel
it was before the widening -- verified on the human sample
(b_narrow 1.8686806129384694, eta 0.20867754234594615, identical). Extending the
widening to discovery needs a contradiction rule that does not read eta, and that
is a separate change with its own validation.

`--polya-end-accept-adapter` is the simpler ablation: admit PASS+ADAPTER by name
(43% -> 63% of human reads), one pooled kernel, reaching BOTH consumers.

**THE RESULT: neither variant is a default. The 3' channel is SATURATED at PASS
coverage.** Paired against the promoted default (`clust` / `gcclust`), expression
contract, bootstrap seed 20260831:

| arm | human (7 cells) | Sequin (104) | SIRV4 (156) | HEYA8 (130) |
| --- | --- | --- | --- | --- |
| `--polya-end-accept-adapter` | +0.07 (6 up/1 dn) | +0.01 [+0.00,+0.04] | **-0.11 [-0.20,-0.03]** | -0.01 [-0.08,+0.07] |
| `--endpoint-qc-conditional` | +0.02 (3/2/2 tie) | 0.00 (0 cells moved) | **-0.08 [-0.15,-0.01]** | **-0.14 [-0.21,-0.08]** |

Human Sn/Pr: 25.17/60.16 -> 25.26/60.09 (adapter), -> 25.20/60.11 (conditional).
Reading 2.2x as many molecules moves human F1 by two hundredths of a point: in a
tie between two candidates what matters is whether ANY read separates them, and
22k PASS reads on three chromosomes already decide the decidable ties. The extra
molecules are mostly redundant.

The spike-in losses are NOT noise and NOT the discovery hook (the conditional arm
never touches it): both are dPr with dSn flat-to-positive -- SIRV4 dSn +0.11 /
dPr -0.25, HEYA8 dSn +0.08 / dPr -0.43. On SIRV4 the per-condition split is
clean:

* the corruption conditions the endpoint channel exists for GAIN --
  c_jitter10/50 +0.24/+0.31, c_merge50 +0.48, c_flip50 +0.60;
* the clean or near-clean annotations LOSE, all of it precision at flat recall --
  full -0.39 (dPr -0.72), p99 -0.40, c_trunc5p10 -0.41, c_ir10 -0.40,
  c_skip10 -0.40.

So the extra 3' evidence redistributes EM mass and lets more candidates clear the
abundance floor: it buys nothing in recall there and costs output. HEYA8 does not
reproduce the split (its jitter/merge conditions lose too, -0.48 to -0.52, with
dPr as large as -1.50 and p00 -0.68), so "help under corruption, hurt when clean"
is a SIRV4 pattern, not a general one -- do not build the targeting rule on it
without measuring HEYA8 first. The general statement the two grids DO agree on is
that the extra molecules cost precision at flat recall.

The lever is mis-targeted rather than broken: it belongs where the annotation's
ends are untrustworthy, which is measurable without a condition name (the channel
already counts how often reads contradict an annotated TES). That is a concrete
next experiment, not something to enable now.

Both flags are therefore committed DEFAULT-OFF as options (commit c29275e).
Default inertness is verified twice: the pooled `EndpointModel` is bit-identical
to the pre-change `_snap/clust` code on the human sample on BOTH fit paths
(annotated: b_narrow 1.8686806129384694, eta 0.20867754234594615; annotation-free
split-half: 1.051438057833809 / 0.18053688247926988), and a grid arm re-runs
three human conditions with the flags off and `cmp`s the GTFs against `gcclust`.

**Two side measurements, both negative, both worth not repeating:**

* *polyA tail LENGTH as an identity term* (`tail_length_icc_probe.py`). It is the
  one quantity only dRNA can measure and the pipeline uses it for nothing. Its
  intraclass correlation across annotated TES sites (172 sites, >= 20 PASS reads
  each) is **0.056**: a single read's tail favours its own site over a random
  other site 58.4% of the time. Real but far too weak to discriminate the
  candidates that are actually tied, which differ by tens of bp at the 3' end.
* *a signal-derived 3' end* would not help either: the pooled narrow kernel is
  already **b_narrow ~1.9 bp** on human and 0.5-2.5 bp on the spike-ins. There is
  no precision left to buy at the 3' end -- only COVERAGE, which is what this
  section is about.

*(Open, not pursued: `internal_priming_probe.py` finds that reads ending on a
genomic A-run of >= 8 nt are 0.87% of reads but sit far from any annotated TES
(A-run 12-20: 2.8% within 20 bp, 67.3% more than 200 bp away) and CLUSTER hard --
374 of 427 fall in 15 ten-bp bins, one holding 170 reads. That is the internal-
priming signature, but the measured tail length does NOT separate them (median
91-123 nt in every A-run bucket), so the discriminator is genomic sequence, not
signal, and whether those sites are artifacts or unannotated APA is unsettled.)*

## 4l. The redundant-model finding, and the first change that is a PURE precision win (2026-09-08, PROMOTED)

**The measurement.** gffcompare's transcript precision has DISTINCT reference
matches in the numerator and ALL query models in the denominator, so a query
model that duplicates another's reference match is a precision loss at zero
recall cost. Verified with the production scorer and flags on a three-reference
fixture: a query carrying one duplicate scores Sn 100.0 / Pr 75.0, the same query
without it Sn 100.0 / Pr 100.0. (`_scorer_fixture/`.)

The production arm carries 4-6.5% such models per grid -- an oracle deletion
estimate of +6.03 Pr on SIRV4, +3.79 on HEYA8, +2.78 on human, +0.00 on Sequin.

**The cause, reproduced end to end.** Human p00, interval chr1:1020099-1056116,
two emitted models with 36 identical exons: `novel_9df5871e6811ef52` (start
1020104, 2 reads) and `novel_72388b83671ff792` (start 1020129, 12 reads). Rerunning
`discover_candidates` on that interval with the run's parameters reproduces both
ids, and at discovery their chains DIFFER by a 3 bp donor wobble --
(1043454,1043537) against (1043457,1043537). `junction_snap` then corrects the
minority junction onto the consensus, at which point they are the same
transcript; but `snap_quant_results` keys its merge on the FULL EXON TUPLE, which
still differs at the first exon. That cell logs "snapped 249 junctions and merged
1 duplicate models" while carrying 141 redundant matches.

*Two things I got wrong on the way here, both caught by Codex and both worth not
repeating: the 24 bp `three_prime_threshold` merge gates (`_collapse_novel`,
`_exact_collapse_novel`, `_three_prime_within_threshold`) are NOT on the
production path -- `discover_candidates` returns from `_chain_cluster_candidates`
before them; and `max(subset) <= max(all)`, so asking whether a gated maximum
moves an end further DOWNSTREAM is a question that can only answer "never".*

**The rule (`--snap-merge-chain`, now DEFAULT ON in `real-drna`).** A second merge
pass over the first pass's survivors. It is deliberately narrow, because one
corrected intron chain is NOT one transcript -- APA and alternative-TSS isoforms
share a chain, and the chain-based scorer cannot see the difference, so a
permissive same-chain collapse would improve the number by destroying the
distinction. A merge requires ALL of: every member novel and multi-exon, one
chromosome, a known strand; the same post-snap chain; the collision CREATED by
the snap (the members' PER-SURVIVOR pre-snap chain sets non-empty and pairwise
DISJOINT, at least one corrected -- pooling them lets a wobble model license the
collapse of two models that already shared a chain); ONE non-null discovery
family across every ORIGINAL member (`_merge_group` relabels a mixed group with
its lexicographically first family, so the merged survivor's own `family_id`
cannot be trusted here); identical strand-aware 3' ends; and a 5' spread within
25 bp as a group DIAMETER, not connected components. The survivor keeps its own
geometry -- no endpoint is re-estimated. `fulllen_frac` is marked unavailable
(the refit does not recompute it). The mode requires `junction_snap` and the
post-selection refit, and is mutually exclusive with `endpoint_refine`.

**The result -- the cleanest in the project so far.** Paired bootstrap, seed
20260831, expression contract, against the previous default:

| grid | dSn | dPr | dF1 | CI95 | up / down / tie |
| --- | --- | --- | --- | --- | --- |
| SIRV4 (156) | **+0.00** | +1.45 | +0.70 | [+0.59,+0.83] | 98 / **0** / 58 |
| HEYA8 (130) | **+0.00** | +0.82 | +0.32 | [+0.25,+0.39] | 62 / **0** / 68 |
| human (7) | **+0.00** | +0.30 | +0.043 | - | 7 / **0** / 0 |
| Sequin (104) | +0.00 | +0.00 | +0.00 | [0,0] | 0 / 0 / 104 |

**Recall is unchanged to two decimals on every condition of every grid** (human
matched references 1849.29 -> 1849.29) and **no cell degrades anywhere**. Largest
gains where duplicates are commonest: SIRV4 p05 +1.65, p00 +1.62; HEYA8 p01 +0.79.

Promotion verified BOTH ways: with `snap_merge_chain` in the profile, three human
conditions under `--profile real-drna` with no extra flag reproduce the
benchmarked arm byte-for-byte (6/6 SAME) AND differ from the pre-promotion
default (3/3 CHANGED) -- the check that once caught a promoted profile value
being resolved and then never applied. Before promotion the flag-off path was
byte-identical to the baseline (6/6). The provenance guard, the family guard and
the CLI application step are each pinned by mutation testing: reverting any one
makes exactly its own tests fail.

**What the benchmark cannot see** (Codex, and it is the honest caveat): the
scorer is blind to terminal coordinates, so "no cell degrades" is necessary and
not sufficient. A genuine nearby alternative start whose models satisfy the
snap-collision rule does disappear; the surviving id and reported 5' coordinate
depend on read support and can move between samples; and the redirect+refit
changes quantification, which structural F1 does not validate. The rule's
narrowness -- exact 3' agreement, bounded 5' diameter, one family, snap-created
provenance -- is what makes that risk acceptable, not the F1 table.

**What is still on the table.** The narrow rule takes 42% / 47% / 16% of the
redundant groups on SIRV4 / HEYA8 / human; the rest are pairs whose 3' ends
genuinely differ (human especially: only 29% have identical 3' ends) or that
involve an annotated model, which is exempt by design. Those need endpoint
reconciliation, not identity reconciliation -- see 4k for why widening the polyA
channel is NOT the way to get it.

## 4m. Where the ceiling is, and what "fully use the signal" should mean (Codex consult, 2026-09-08)

The owner asked two questions at production level: what can still improve, and
how to FULLY utilise the raw signal. A full evidence dossier went to Codex
(gpt-6-astra); its answer, with the parts I checked, is the current plan.

**The verdict.** *Stop treating "consume more raw signal" as the objective; the
objective is additional IDENTIFIABLE evidence, and most proposed channels cannot
distinguish the hypotheses that currently survive.* Keep the default, run one
narrow precision experiment, and move real effort to the axes this benchmark
does not measure.

**A. The ceiling is not a number we have.** The "63.8% of expressed transcripts
have no full-chain read / we are above the 29.2% ceiling" headline conflates
three quantities (scorer recall with no precision constraint, which is ~100% by
echoing; recall supported by single complete-chain reads; recall identifiable
from ALL molecules at a stated error rate -- the one that matters, unmeasured).
Codex also found two definitional problems in `human_recall_decomposition.py`:
the "no read" bucket is really "no exact-chain read AND not emitted", and the
mono-exon read key `(chrom, strand, ())` carries no coordinates, so a mono read
anywhere on the chromosome counts as support. Checked: mono-exon is 2.1% of the
truth (157 / 7398), so the distortion is bounded and in the direction of
UNDERSTATING fragmentation; the qualitative finding stands, the percentages are
not a ceiling. The distinction to measure instead: **uniquely reconstructible /
statistically distinguishable / observationally indistinguishable** -- junction
coverage is not junction phasing, and assembling the Cartesian product of
unlinked choices invents identities.

**A, on the GTF exemption.** Do not override maximin on "it is right on human":
it agrees better with NanoCount-derived labels, and both can favour spreading
support over annotated isoforms the molecules do not distinguish. A human
operating point is defensible only as a PROSPECTIVE product decision validated
on independent human samples (not seven perturbations of one), with absent-
isoform controls and rescued identities labelled prior-supported.

**B. The one precision mechanism to build next:** extend the snap-collision
merge to a pair of ONE annotated + ONE novel model, under every existing guard
(snap-created, exact 3', bounded 5' diameter, original-family provenance). It is
representation cleanup of a collision the pipeline itself created, not a new
biological classifier. Do NOT start with different-3'-end pairs: that is an
APA-vs-artifact inference, and neither "few reads" nor "A-rich and clustered"
decides it -- report endpoint uncertainty rather than delete a biological
possibility to please an intron-chain scorer. Expected: 0 to +0.5 SIRV4, ~0
Sequin, 0 to +0.2 HEYA8, 0 to +0.1 human. It cannot close the 7-point espresso
gap; the pre-merge redundancy oracle was only +2.78 Pr on human.

**C. The raw signal, ranked by expected value (structural benchmark):**
reliability calibration > narrow common-sample contrast > tail-boundary audit >
terminal-signal audit > modification/structure identity.

1. *Common-sample, alignment-marginalised nested contrast* -- principled, but
   NOT a broad rescue of the 51.7%: when the read observes only the identical
   suffix the likelihood ratio is neutral under any sequence-conditioned model.
   Cheap kill test first: count contrasts whose candidates differ INSIDE the
   recoverable common signal; if negligible, stop before HMM engineering.
   `unequal_intron_count` is a coarse tie-set label, not proof of containment.
2. *Dwell/current as modification or structure identity* -- dead end for now.
   The missing link is prediction of UNOBSERVED structure conditional on shared
   sequence; an isoform-modification association found after sequence-based
   assignment runs the inference backwards. Credible test: establish identity
   from uniquely informative molecules, hide the distinguishing sequence,
   predict from shared-region signal, on independent runs.
3. *PolyA tail signal beyond length* -- low structural value; the honest
   question is boundary validity (adapter / A-rich / body segmentation), and a
   genomic A and an added A do not encode their origin. Never train an artifact
   classifier with "far from annotated TES" as the label.
4. *Signal-level TSS* -- our probe measured "no large hidden segment after the
   first aligned event", NOT terminal waveform / cap / motor release. One
   bounded audit with capped-vs-fragmented controls is allowed; without labels
   it is not a production project. A broken molecule also has an end.
5. *Read reliability as a covariate* -- the most practical untested lever.
   Local, candidate-NEUTRAL reliability measured in shared flanks tempers a
   contrast's influence; never down-weight a read for fitting the preferred
   candidate poorly; separate identity evidence from molecule count; treat many
   samples from one k-mer run as ONE vote. Expect calibration gains, ~0 F1.
6. *The DRS-native advantage we are not using:* joint molecular
   characterisation with honest assignment uncertainty -- chain-linked tail
   length and (independently validated) modification distributions, reported at
   the identifiable equivalence-group level when transcript identity is
   ambiguous. Low tail ICC makes tail length a weak IDENTITY feature; it does not
   make tail regulation an uninteresting OUTPUT.

**D. Stop / re-baseline:** the P2 sidecar's "44 loss cells" and +5 F1 gates
belong to an older operating point -- re-baseline before continuing; defense
election on p00: deprioritise; tail-length identity, step-5a correction, M3 DTW,
broad same-chain collapse, more opaque current/dwell features, hunting a subset
where the failed polyA widening wins: stop. Raise real-human revalidation on
INDEPENDENT samples to the top -- perturbation cells of one library are not
biological replication, and repeatedly inspected held-out grids are now
regression panels, not renewable confirmation.

**E. The plan:** (1) narrow annotated/novel collision reconciliation, default-off,
frozen gates; (2) local reliability calibration, expect ~0 F1, judged on held-out
contrast calibration; (3) freeze assembly and evaluate the OTHER axis --
log-ratio error against Sequin molar inputs (not Spearman alone),
abundance-stratified and replicate-consistent; mass in non-identifiable groups
and sensitivity to reassignment (softmax confidence is not calibrated
uncertainty); split-library reproducibility of chain-linked tail distributions;
cold/warm wall time and peak memory.

## 4n. Should the junction snap use raw signal? (2026-09-08, audit DONE -- no)

The owner's premise: raw signal is intrinsically better information and we are
under-using it; concretely, should `junction_snap` (count-based, production
default since 2026-08-27) take a signal vote?

**Three facts that bound the answer, all from our own records.**

1. **The election ablation, already in `SIRV4_IMPROVEMENT_PLAN.md`:** on the same
   24 cells with the same nomination screen, `screen` alone 92.44 F1, `m1`
   (base-space ballot) 89.15, `signal` ballot 89.13, `shuffle` 87.33, baseline
   87.29. The signal vote equals the base-space vote and BOTH are below the
   screen acting alone. The 97.6% "decisive ballot accuracy" is real, but the
   incremental value of the raw current over the basecalled alignment on that
   population was zero.
2. **Per-read signal correction is harmful** (SIRV4, margin 10): read moves
   own_TP->FP 612 vs own_FP->TP 448; new candidates 12 TP / 116 FP. Single-read
   scores are dominated by one repeated k-mer run.
3. **Snap's failure modes split by what signal can touch.** (a) ratio-abstained
   wobble pairs (both survive, precision loss) -- signal could help; (b) a
   genuine tandem site "corrected" away (recall loss; common in human, rare in
   SIRV) -- signal cannot help and a forced vote makes it WORSE. Snap's exact
   rule (Codex correction): support is a fresh primary-BAM scan, not the
   registry, and the target needs `max(min_support, floor(cur*ratio)+1)` --
   strictly MORE than 2x, so 2-vs-4 abstains and 2-vs-5 moves.

**Codex's verdict: do not put the current signal scores into production snap;
run the audit as an OPPORTUNITY assessment.** "Signal prefers A" is not
sufficient evidence to merge B into A -- the missing test is whether B has its
own credible molecular support. On the premise itself: raw signal is a richer
MEASUREMENT, not automatically a better DECISION; the basecaller already
integrates the current with a contextual model, so going back to it competes
with that inference rather than starting from untouched evidence. It recovers
discarded information only when the alternatives predict distinguishable
sequence inside the observed signal AND both align comparably; when candidates
predict the same observed sequence it is re-reading the same evidence with a
weaker scorer.

**If the audit justifies anything, the shape is:** abstention-only, aggregate,
three outcomes (collapse to A / collapse to B / preserve both), never overriding
a count decision, with a MINORITY VETO -- >= 3 technically qualified molecules
robustly supporting the loser => preserve both. The decisive statistic is not
the pooled majority but *among reads whose own chain carried the loser, do their
signal contrasts defect to the winner?* (a wobble's reads defect; a genuine
alternative keeps a loyal subset). A count x signal posterior is NOT the next
step: placement and current derive from the same underlying signal, so
multiplying them double-counts the observation.

**What the audit CANNOT deliver:** the per-read `junction_correction.tsv`
records total NLL and event counts only -- no per-run contributions, so the
"one k-mer run = one vote" robust ballot is not reconstructible from it. Totals
are a diagnostic baseline; a decision-grade vote needs retained event-level
alignments.

**Pre-registered gates (fixed before the numbers):** build a robust-scoring
prototype only if >= 100 actionable one-catalogued abstention loci across >= 50
genes AND a model-level oracle >= 0.2 F1; the vote rule if built: >= 10
qualified molecules, >= 90% winner, >= 3 from each origin group, >= 80% of
loser-origin molecules defect, preserve both if >= 3 loyal; zero merges on the
both-catalogued preservation panel; must beat the best base-space comparator
with a gene-clustered paired CI excluding zero; >= 0.1 F1 and no recall loss on
>= 3 independent human libraries (`full`/`p00` are one library).

The audit (`sigaudit.sbatch`, byte-inert: election + correction in audit mode on
human `full` + `p00`, krill 0.1.6 overlay) and its analysis
(`snap_signal_audit.py`, which replays snap's exact rule from the same BAM scan
and computes the loser-origin defection) were launched before the gates were
written down.

**RESULT (human H9 rep2 run2, chr1/19/21; both cells byte-identical to
production; `full` and `p00` give IDENTICAL numbers because the per-read audit
runs on raw reads before any annotation is consulted -- leakage-free by
construction).** 23803 registry junctions, 2031 pairs within 6 bp at both
coordinates, 146267 usable per-read signal ballots over 3107 pairs.

| snap replay | both annotated | one annotated | neither |
| --- | --- | --- | --- |
| moved | **95** | 1385 | 94 |
| abstained | 65 | **202** | 190 |

* **(a) the population signal could address -- one-annotated abstentions: 202.**
  The signal reaches quorum + supermajority on **18** of them (9%), and those 18
  agree with the annotation 100%. The pre-registered gate was >= 100 actionable
  loci across >= 50 genes: **fails by 5x.**
* **Where count and signal both act, the signal is not better:** on the 1385
  one-annotated pairs the count moved, count is right **99.9%**, signal **98.9%**
  (1067 called). Same finding as the election ablation, on human.
* **The loser-origin test** (do the count loser's own reads defect?) behaves as
  Codex predicted: on both-annotated abstentions the loser's reads stay loyal
  (median defect 0.17; 16 of 21 testable pairs keep >= 3 loyal reads), on
  one-annotated abstentions it is a coin flip (median 0.50). So a minority veto
  would correctly preserve tandem sites -- but there are only 65 such
  abstentions, and a forced signal vote would collapse 4 of them.
* **The flip -- signal loyalty as a VETO on count moves -- is net negative:** it
  would protect 17 of the 95 collapsed real tandem sites (28 testable) but block
  60 of 214 testable CORRECT moves.

**Verdict, first pass: no signal in snap** -- ~18 addressable loci, signal
slightly less accurate than count where both act, neither form survives its
gate.

**Then the owner objected, correctly, and the objection was tested
(`snap_signal_independence.py`).** The principle: two independent ~90% features
should combine to beat either. On the 1085 one-annotated pairs where both
features decide:

| | signal right | signal wrong |
| --- | --- | --- |
| count right | 1068 | 12 |
| count wrong | **5** | **0** |

Count 99.54%, signal 98.89%. The errors ARE independent: the signal is right on
all 5 count errors, joint-wrong observed 0 against 0.055 expected under
independence. So my "same source, so correlated errors" argument was WRONG on
this data. What still makes a hard combination negative is the combining rule:
letting the signal override fixes 5 and breaks 12 (net -7); with two features
and a hard disagreement you follow the more accurate one, and the pair equals
its better member. The principle needs a third voter or per-case calibration.

**Calibration exists, and it is the COUNT's margin, not the signal's.** The
signal's own supermajority is 0.75-1.00 on both the fixes and the breaks --
uninformative. The count ratio separates them almost perfectly: the 5 fixes sit
at ratio 1.3-3.2, the 12 breaks at 3.0-28. Rule "defer to the signal only when
the count ratio < 2.5": fixes 4, breaks 0, **net +4** -- and all four fall in
snap's ABSTENTION region, i.e. exactly the abstention-only tie-breaker Codex
specified. Extrapolated genome-wide that is ~15-25 junctions; the F1 effect
depends on whether each corrected model is its chain's only carrier, and is
0.0x at best.

**Corrected verdict:** the owner's principle holds and the data confirm
independence; the abstention-only signal tie-breaker is a small, real, safe
positive (+4 corrected pairs on chr1/19/21, 0 broken), buildable as opt-in under
the gates above. It is not the answer to "we are under-using the signal": the
count is already at 99.5% on this decision, so a second opinion can only touch
the residual 0.5%. The signal pays where the basecaller makes NO statement
(polyA), not where it already makes a 99.5%-correct one.

**The secondary finding is base-space and worth a follow-up:** snap's count rule
collapses **95 pairs where both coordinates are annotated** -- failure mode (b),
more than 5x the size of (a). Whether that costs recall depends on whether the
collapsed model was the only carrier of its chain; the original three-sample
validation showed recall UP under snap, so the net is positive, but a count-only
minority guard (do not move when the loser has >= N indel-clean reads of its own)
is a cheap base-space experiment with a clear falsifier.

## 4o. The human recall gap has ONE owner, and an exemption that closes most of it (2026-09-09/10; LANDED opt-in 82f0151, PROMOTED to default 17fec28 -- see 4t)

Full record (Chinese): `docs/PRODUCTION_SUPPORT_GATE_ATTESTATION_2026-09-10.md`.
Literature review commissioned for this campaign: `docs/LITERATURE_REVIEW_2026-09-09.md`.
(Status line at the time of writing: everything below was uncommitted. It has
since landed; the K=0 default described here is now K=3 under `real-drna`.)

**Where the stringtie3-only matches die.** With the complete human curve in
hand (4b5: 13/26, crossover at p80, every loss is recall), the `full` cell was
dissected: 554 stringtie3-only `=` (NanoCount median 1 read; 272 of them have a
BAM read carrying their exact chain). A new opt-in `--selection-audit` (per-
interval shards + a global file recording every SelectionOutcome from
`discovered` to `final`; collector `selection_audit_collect.py`; GTF verified
structurally identical with the flag on) attributes **475 of the 554 to
`m2_support_gate`**, and **4566 of all 7398 expressed truth transcripts (62%)**
to the same gate. None of the 554 is lost in discovery (all are GTF candidates;
none folded). The gate drops 20,294 of the 31,418 candidates entering EM on
that cell. This is 4i's finding with an exact attribution instead of an
inference.

**What the dropped annotated candidates look like** (`gate_drop_precision.py`,
truth = NanoCount TPM>0, "suffix" = read chain is a strand-aware contiguous 3'
suffix of the candidate chain, count includes exact reads):

| tier | n | truth |
| --- | --- | --- |
| exact reads >= 3 | 122 | 75% |
| exact == 2 / == 1 | 71 / 347 | 65% / 45% |
| no exact, suffix >= 3 | 1367 | **77.5%** |
| no exact, suffix 1-2 | 2160 | 59% |
| no chain read at all | 15929 | 12% |

Exact-or-suffix >= 3: 1595 candidates at 77% truth, above the cell's 67.9%
precision -- admitting them raises Sn AND Pr. Codex's caveat stands: that 77%
is agreement with a compatibility-based pseudo-truth, not independent evidence
of expression.

**The lever: `--m2-support-gate-attest-reads K` (flag default 0 = byte-identical; `real-drna` profile default 3 since 17fec28).**
An annotated multi-exon candidate is exempt from the support gate iff
`compatible_read_depth >= K` (distinct MAPQ>=1 primary reads whose chain equals
the candidate's or is its strand-aware 3' suffix, not running across the
upstream model intron, 3' end within 50 bp past the annotated TES; computed
always in `discovery.chain_compatible_read_depths`) AND no unobserved junction
of it is a JITTER of an observed one (both boundaries within
`--m2-support-gate-attest-veto-bp` 10 of a junction in the MAPQ>=1 registry;
`m2_junction_nll.unjittered_gtf_columns`, fails closed). No EM-mass bar. Three
things were tried and rejected on data, each with Codex: requiring every
junction of the read to be indel-clean (v1: a read-length filter on dRNA,
1079 -> 210 candidates for 95% vs 90% precision); the OR-shaped
`uncontradicted_gtf_columns` veto (vetoed 129 attested candidates, 128 of them
expressed -- alternative splice sites are not contradictions); and an
"exclusive" rule where an exact-chain read attests only its own chain (v3:
Codex -- an exact chain match does not make a molecule complete; inert on the
corruption ladder, -0.3 on human).

**FINAL (v6 = the committed code; snapshot `_snap/attest6`, verified
byte-identical to the pre-commit tree `_snap/attest8` on two human cells):**

| dataset | cells/libs | dSn | dPr | dF1 [95% CI] | up/down/same | worst family |
| --- | --- | --- | --- | --- | --- | --- |
| human | 26/1 | +1.65 | +0.42 | **+1.63 [+1.25, +2.02]** | **24/0/2** | c_jitter +0.34 |
| SIRV4 | 156/6 | +0.36 | -0.15 | +0.09 [-0.11, +0.34] | 37/16/103 | c_trunc3p +0.07 |
| Sequin | 104/4 | +0.00 | -0.03 | -0.01 [-0.03, +0.00] | 0/1/103 | c_ir -0.15 |
| HEYA8 | 130/5 | +0.10 | -0.02 | +0.05 [+0.01, +0.12] | 12/2/116 | c_jitter +0.05 |

Human vs stringtie3 on the end-blind contract: production 13/26 (mean +1.41)
-> **25/26 (mean +3.04)**; the only loss left is p00 (-0.30), where the lever
is inert because there is no annotation. K=0 control on human `full` is
BYTE-identical to production. **Under an end-aware contract the human picture
is different and we do not lead -- see 4p.**

**These are NOT a confirmation.** Rule 2's shape was chosen after inspecting
the v2 failures, so the panels are post-selection characterization; the
bootstrap intervals do not account for that selection; one human library gives
a within-library bound, not a population one; precision falls slightly on three
panels. The flag ships DEFAULT-OFF and no profile enables it.

**A bug the verification caught, worth remembering:** `collect_compatible_depths`
was added to the outer and inner signatures but the FORWARD line was lost in a
script exception, so the gated build silently never fired the exemption (output
fell back to production's 3377 models). The unit test only checked that the
runner passed the flag. Asserting equivalence is not verifying it.

**Sentinels (one library each, paired vs the production default):**

| cell | v2 |
| --- | --- |
| human `full` | Sn 31.2 -> **34.8**, Pr 67.9 -> **69.3**, F1 42.75 -> **46.33** (+3.58; stringtie3 44.68) |
| human c_jitter10 / c_jitter50 / p00 | +0.83 / +0.06 / 0.00 |
| SIRV4 one library, 26 cells | mean -0.12 [-0.28, 0.00]; c_ir50 -0.9, c_merge50 -0.9, c_skip50 -1.3; 23/26 identical |
| HEYA8 one library, 26 cells | mean -0.03; c_skip50 -0.68; 24/26 identical |
| K=0 same-snapshot control | model multiset identical to production |

The spike-in FPs are the corrupted entries themselves: c_merge50's SIRV409 is
SIRV409's 5' part merged onto SIRV410's chain with SIRV410's own entry REMOVED,
so SIRV410's 35 reads (equal to no annotated chain) attest the chimera as a
suffix; c_skip50/c_ir50 are the same shape. **A corrupted junction the reads
never cover is silence, and chain compatibility cannot tell it from a human
degradation fragment.** A TSS-peak block (reads with piled-up 5' ends are
complete molecules, not fragments) was measured offline (costs 181 of 874 human
truths at M>=5) and NOT adopted -- Codex: a 5' pile-up implies neither
completeness nor zero contribution from longer isoforms.

**Codex's ruling on the contract:** the pre-frozen stop rule (no corruption
family mean below -0.5) was breached by c_skip on the SIRV4 sentinel, so the
panels now running (human 26 + SIRV4 156 + Sequin 100 + HEYA8 125 on the
repaired v2, snapshot `_snap/attest4`, arms `gcatt4` / `att4`) are
CHARACTERIZATION, not the continuation of an unbroken promotion test. A
repaired arm needs a newly frozen protocol: human mean dF1 >= +1 with paired
LCB > 0; every other dataset LCB > -0.25; every corruption family mean >= -0.5;
missing cells block; no cross-dataset offsetting; library-block bootstrap.
`arm_compare.py` now lists missing cells, maps families explicitly and
bootstraps libraries as blocks.

**What this lever IS (Codex's wording, adopted):** an annotation-prior rescue
that exempts multi-exon models from the support gate when >= 3 distinct
eligible primary reads carry an exact or 3'-suffix chain and pass the span
checks, provided no unobserved junction has an observed neighbour within 10 bp
on both boundaries. It raised F1 3.58 on the inspected human full-annotation
cell with both axes up and beat stringtie3 there under the NanoCount contract;
it can also rescue incorrect annotated structures whose distinguishing
upstream region is unobserved (merge/skip/IR sentinel FPs, 1-2 models per
50%-corruption cell).

## 4p. Two dimensions nobody had measured, and what they say (2026-09-10)

Full record (Chinese, both parts): `docs/PRODUCTION_QUANT_ACCURACY_2026-09-10.md`.
Harnesses: `quant_accuracy.py`, `strict_match_rescore.py`. No production code
changed for either. Both conclusions were REVISED after Codex review -- the
first versions of both harnesses flattered us.

### QUANTIFICATION, against real concentrations

Every result in this file scores STRUCTURE against a pseudo-truth (NanoCount
TPM>0). The Sequin grid has actual quantitative truth --
`refs/concentration/RNAsequins_MixA_v1.tsv`, attomol/uL, 164/165 ids matching
our Sequin truth GTF (the grid is `sequin_v1_...`, so **v1, not V2**). 104
paired cells over 4 libraries, pyfin production default vs stringtie3:

| metric | denominator | mean d (pyfin - st3) | pyfin better |
| --- | --- | --- | --- |
| rho | each tool's OWN detections | -0.004 | 42/104 |
| rho | **COMMON detections** | **-0.022** | 44/103 |
| MARD (lower better) | COMMON detections | **+0.056** | 32/103 |
| rho_all | expressed set, misses = 0 | **+0.126** | 75/104 |
| MARD_all (lower better) | expressed set, misses = 0 | **-0.211** | 80/104 |
| unmatched TPM fraction | whole output | **0.044 vs 0.145** | - |

**On the transcripts both tools found, stringtie3's abundances are slightly
BETTER than ours.** Our advantage in the joint metrics is entirely detection.
The one abundance-side win is that only 4.4% of our TPM sits on models matching
no truth transcript, against stringtie3's 14.5%. We are also behind the
NanoCount oracle (which is handed the true annotation): 0.756/0.661/0.631 vs
our 0.725/0.543/0.707 unconditional. Never claim a better quantifier.

Codex fixed two things that had inflated us: comparing each tool on its OWN
detected set compares different problems, and the first version renormalised
away the abundance a tool put on wrong structures (no penalty for mass
misallocation). An even earlier version used all 164 mix transcripts as the
denominator, which saturates the symmetric MARD at its 2.0 ceiling for every
tool including the oracle and made us look better than it.

### END-AWARE STRUCTURE (`gffcompare --strict-match`)

gffcompare's default `=` is intron-chain only -- its own docs say matching
transcripts may differ in the length of the first and last exons -- so **every
precision number in this file is TSS/TES-blind**.

**Convention trap, measured on 0.12.6:** with `--strict-match` the `.tmap`
class code demands EXACT terminal coordinates (all 2167 survivors on a real
human cell have max |offset| = 0; a 19 bp shift is `~` even at `-e 25`) while
the `.stats` transcript level honours `-e` (same run: 2557 blind, 2468 at
e=100, 2333 at e=25, byte-identical tmap codes). The F1 numbers below are the
`.stats`, i.e. the tolerance convention; attribution must be recomputed from
coordinates, not read off the tmap.

SIRV4, 4 conditions x 6 libraries, mean F1 blind -> e=100 -> e=25:
pyfin **85.67 -> 85.24 -> 85.10** (-0.57); stringtie3 69.18 -> 64.28 ->
**59.64** (-9.54, Pr 79.4 -> 66.9); espresso 78.91 flat; isoquant 56.00 ->
54.80. **On the spike-ins strictness widens our lead, +16.5 -> +25.5.**

Human, 26 conditions, 1 library, mean dF1 vs stringtie3 (cells won):
blind +1.41 (13/26) production, **+3.60 (25/26)** with attestation; e=100
-0.26 / +2.03; **e=25 -1.52 / +0.70 (12/26)**. Per condition: with a complete
annotation stringtie3 inherits correct ends from `-G` and beats us at e=25
(full 44.54 vs 42.31); with a degraded or corrupted annotation we win under
every contract (c_jitter50 12.62 vs 8.00, c_flip50 29.19 vs 23.97).

**Cause, at the reference level** (own coordinate recount, reconciled with the
stats to within 4 of 2554): at 100 bp we lose 90 references, 87 of them matched
ONLY by an `.ext`; at **25 bp we lose 222 (8.7%), 219 of them (99%) ext-only**;
at exact ends 394 (389 ext-only). Of the 397 `.ext` models only 169 have both
boundaries within 25 bp; their 5' ends mostly extend UPSTREAM of the annotation
(median 15 bp) and their 3' ends mostly fall short. That is
`read_extent_contradicted` working as designed -- it fires only when reads
OVERRUN an annotated end, at the fitted narrow kernel (b_narrow ~1.9 bp on
human) gated by 3 clustered reads -- so three reads 15 bp upstream overturn an
annotated TSS.

**Whether that is an error, this benchmark cannot say**: truth IS the GENCODE
end, but LRGASP reports high reference-TSS agreement with low CAGE support.
Codex's ranking of the decisive experiments: (1) stage a second human library
and test whether the same `.ext` TSS reproduces (strongest, but replication
shows reproducible structure, not cap-defined initiation); (2) check whether
the `.ext` 5' end coincides with an annotated TSS of ANOTHER isoform of the
gene (diagnoses mis-assignment, not expression); (3) act on the 3' side only
(an ablation, not evidence about the 5' end).

**The lever is NOT "emit both models."** `resolve_extent_pairs` deliberately
keeps one (more EM mass wins, ties to the annotation) and `snap_merge_chain`
never merges the pair (it is novel-only). The lever is a MAGNITUDE FLOOR on the
contradiction test. It is not free -- the `.ext` competes in EM, so removing it
moves read assignment -- and Codex's condition is that any floor be calibrated
on independent endpoint evidence, **not set at 25 bp because that rescues this
benchmark**.

### The honest one-sentence claim (Codex's wording, adopted)

> Across the evaluated annotation-perturbation panels, PyFIN shows strong
> intron-chain recovery and favorable joint detection-abundance performance on
> Sequin; its relative endpoint accuracy is dataset-dependent, with robust
> SIRV4 agreement but weaker agreement with complete human reference
> annotations at a 25-bp tolerance.

## 4q. The quantification gap was a REPORTING-CONVENTION bug (2026-09-10, LANDED opt-in 8782013)

Full record: `docs/PRODUCTION_QUANT_ACCURACY_2026-09-10.md` part 3.

`scores.tsv`'s `tpm` column is `abundance / length_kb` -- the SHORT-READ
convention, which exists because fragmentation makes a longer transcript yield
proportionally more reads. Direct RNA sequences whole molecules (one molecule,
at most one read, regardless of length), so the count does not carry that
factor and dividing by length gives an INVERSE-length-weighted quantity (not a
"mass-like" one -- mass is molarity TIMES length; Codex corrected that).
NanoCount and `salmon --ont` omit the correction for the same reason.

Measured against the Sequin spike-in's known INPUT MOLARITY -- the only real
quantitative truth this project has -- 104 evaluated cells, 103 eligible
common-set comparisons across 4 libraries, library-block bootstrap:

| comparison | length-normalised | count |
| --- | --- | --- |
| paired change WITHIN pyfin, Spearman | - | **+0.074 [+0.036, +0.112]**, 83/103 |
| paired change WITHIN pyfin, MARD | - | **-0.104 [-0.134, -0.068]**, 94/103 |
| vs stringtie3, common set, Spearman | -0.022 [-0.099, +0.053], 44/103 | **+0.052 [+0.013, +0.091]**, 82/103 |
| vs NanoCount (reference-assisted), joint Spearman / MARD | 0.725 / 0.707 | **0.788 / 0.611** (its 0.756 / 0.631) |

Every library's mean Spearman improves (0.678->0.694, 0.830->0.889,
0.588->0.683, 0.694->0.823). **This is what 4p's "we are slightly behind
stringtie3 per transcript" actually was: a reporting convention, not an
estimator gap.**

**It is not uniformly better, and my first explanation of why was wrong in the
direction.** Among 124 matched transcript-cell observations with length >= 2 kb
the count has the HIGHER median symmetric error (0.911 vs 0.646); signed, it
OVER-estimates them (median log2 estimate/truth +1.42 vs -0.45) -- so length
division is not "compensating for long-molecule loss", it is pulling an
over-estimate down. Recorded, not explained. Bin the errors after normalising
ONCE per cell; normalising inside each bin destroys the between-bin effect the
length question lives in (I got that wrong first too).

`--no-tpm-length-normalize`, default unchanged. `compute_tpm`'s only consumers
are the scores.tsv column and the refit diagnostics (which take the same flag,
so "tpm" means one thing per run); verified on two real cells -- GTF
byte-identical, scores.tsv differs only in `tpm`, new `tpm` == 1e6 *
abundance / sum on every row. Codex favours making it the `real-drna` default
in a small explicit follow-up that also exercises `scores.unfiltered.tsv` and
the refit JSON numerically and adds an output-semantics note.

**The endpoint floor shipped in the same commit is a NEGATIVE result kept at 0**
(`--tss-contradiction-min-bp` / `--tes-...`): motivated by an independent
measurement (`annotation_end_precision.py`: the contradiction threshold uses the
narrow component, ~2 bp on human, while dRNA read starts disperse 20-50 bp
around an annotated TSS with the mode 20 bp INSIDE it; the 3' end is measured an
order of magnitude better, median 8 bp), but at 5' = 20 bp it suppressed 145 of
515 read-extent siblings and moved 126 transcripts onto annotated geometry for
+0.06 end-aware F1 -- because the suppressed ones were already inside the 25 bp
tolerance (median offset 12) while the failures sit at exactly 25. A floor big
enough to catch those exceeds the measured resolvable scale. **The end-aware
human gap is not fixable by a calibrated contradiction floor.**

## 4r. The endpoint dimension, answered by replication (2026-09-10)

Full record: `docs/PRODUCTION_QUANT_ACCURACY_2026-09-10.md` 8.2-8.2e.

4p left one open algorithmic question: under an end-aware contract our human
result is carried by `.ext` read-extent siblings whose 5' end differs from the
annotation, and the benchmark cannot say who is right. **It is now answered,
by replication rather than by the benchmark**, and the test needed no pipeline
run -- only a second library's BAM, and all 15 human libraries have one.

**Pre-registered** (`ext_tss_replication.py` docstring, written before any
result): population = the 133 `.ext` models whose 5' end is >25 bp from the
annotated TSS; "reproduces" in a library iff >= 3 chain-compatible MAPQ>=1
reads start within 10 bp of the claimed start AND more than start within 10 bp
of the annotated TSS; NULL = the mirror position t+(t-s) under the identical
rule.

**Result: 19/133 (14%) reproduce in >= 2 of 4 libraries, against an 8% mirror
null** (per library 12/7%, 21/8%, 13/8%, 7/6%). By the pre-registered rule that
is "~ null": **86% of the `.ext` 5' claims are library-specific read-start
dispersion, not reproducible structure.** And **no single-library statistic
finds the reproducible minority** (base rate 15%: "claimed start is the modal
start" 17% precision, ">= 30/50/70% read share" 14-16%; the EXISTING >= 3-read
rule is the best at 33% precision / 95% recall).

**So turning the 5' contradiction off is supported by independent evidence --
but it is a DOMAIN trade, not a fix:**

| | end-blind F1 | end-aware e=25 | `.ext` emitted |
| --- | --- | --- | --- |
| attestation only | 45.51 | 41.38 (2272 refs) | 515 |
| + 5' contradiction OFF | 45.40 (-0.11) | **43.29 (+1.91, 2381)** | 212 |

Human gains +1.91 end-aware F1 for -0.11 end-blind (the e=25 gap to stringtie3
on `full` narrows from -3.16 to -1.25). **SIRV4 pays -0.45 [-0.51, -0.40] over
156 cells, with `c_trunc5p` at -6.10 (worst -11.15)** and every other family
within +-0.05 -- because c_trunc5p IS the truncated-5'-annotation corruption
and SIRV reads are full-length, so the extension is RIGHT there. That is 4g's
mechanism doing its job.

**The discriminator between the two regimes is measurable per sample, and it is
the TAIL, not the median** (`annotation_end_precision.py`; SIRV 55 chains vs
human 1111): 5' |offset| median 14 vs 22 bp -- nearly the same, that is the
motor/degradation constant -- but p90 **20 vs 253** and **100% vs 56% within
25 bp**. On SIRV every locus's read-start mode tracks the annotated TSS; on
human barely half do. Where the starts track the annotation an outlier cluster
is evidence; where they do not it is noise.

**The gate is now implemented** (`--tss-contradiction-auto-gate`, default off).
`sample_nuisance.annotated_tss_tracking` computes the statistic from data the
run has ALREADY loaded (`starts_by_locus` + `_annotated_tss()`), so it costs no
new BAM pass: the fraction of annotated TSSs whose local read-start MODE lands
within `--tss-contradiction-track-bp` (25). Below `--tss-contradiction-track-frac`
the 5' contradiction is disabled for that run and the measurement is logged; the
3' contradiction is never gated; an unscorable sample leaves it ENABLED (unknown
is not evidence).

**The threshold is the midpoint of an empty band, not a tuned value.** Measured
on 10 libraries across 4 datasets: human dRNA 0.491 / 0.496 / 0.498 / 0.500 /
0.514, SIRV4 0.797 / 0.800 / 0.814, Sequin 0.917 / 0.950 -- nothing between
0.52 and 0.79. Default 0.65.

**The first version of the gate was falsified by its own test, and the failure
is the interesting part.** It measured only whether read starts LAND on
annotated TSSs -- which conflates two opposite causes of not landing: the reads
are imprecise (human, 0.50) or the ANNOTATION is wrong (`c_trunc5p50`, 0.553),
and the second is exactly what the repair exists for. The gate therefore switched
the repair off precisely where it was needed: **-4.88 F1 on the c_trunc5p family
(worst -11.15)** while every other family moved 0.00.

**The corrected statistic adds one clause and is immune to it**: a locus counts
as RESOLVED when its modal read start lands on the annotated TSS **or** holds a
majority of the local starts -- a tight cluster away from the annotation is
evidence against the ANNOTATION, not against the reads. Measured, SIRV4 goes
0.957 (intact annotation) -> 0.959 (5'-truncated) where the tracking-only
version collapsed 0.797 -> 0.438; human sits at 0.784-0.834. Default threshold
0.89, between them.

**Verified (`_snap/autogate2`, arms `gcag2` / `ag2`):**

| grid | corrected gate vs attestation-only | gate decision |
| --- | --- | --- |
| SIRV4 156 cells / 6 libs | **0.00, 156/156 identical** | enabled (0.968) |
| Sequin 104 / 4 | **0.00, 104/104 identical** | enabled |
| human `full` | **byte-identical to the manual 5'-off arm** | disabled (0.834) |

Human `full` under the gate: end-blind F1 45.40 (attestation-only 45.51,
production 42.75) and **end-aware e=25 F1 43.29 (attestation-only 41.38,
production 38.61; stringtie3 44.54)**. One setting, measured per sample,
produces the right behaviour in both domains. HEYA8 and the other 13 human
conditions are still running; the flag stays default-OFF until they land.

## 4s. The M2 signal tie-break: its value IS the signal, and it survives real human data (2026-09-10/11, LANDED d4e63fc + docs)

Record: `docs/PRODUCTION_M2_GUIDED_PATH_2026-09-10.md`.

**What production actually runs** (from the manifests, not the profile
sentinel): `real-drna` PINS `m2_metric="summed_llr"` (route `fixed`), so the
tie-break scores the TIGHT differing-junction window (`m2_summed_llr_flank=8`)
with a REAL abstention threshold (`m2_summed_llr_margin=1.0`). The wide-window
`auto -> mean` belongs to the `sirv` profile. `m2_tiebreak_margin=1e-9` is the
`argmax` mode's threshold and is not on this path.

**Three-arm decomposition.** `--m2-metric off` removes the signal AND restores
M1's 1/K split, so it cannot separate "signal is right" from "concentrating a
tie onto any structure helps". The new `--m2-scramble-seed` (default 0 = off,
never a production setting) keeps the concentration and removes only the
signal: same tie population, same abstentions, same thresholded margin, only
the winning STRUCTURE replaced by a blake2b hash.

| contrast | SIRV4 dF1 [CI] | HEYA8 dF1 [CI] | human dF1 |
|---|---|---|---|
| off - base (signal + concentration) | -0.93 [-1.46,-0.49] | -0.76 [-0.96,-0.56] | -0.27 |
| scram - base (signal only) | -0.75 [-1.05,-0.45] | -0.59 [-0.90,-0.32] | -0.14 |
| scram - off (concentration only) | +0.18 [-0.27,+0.60] | +0.17 [-0.23,+0.50] | +0.13 |

Concentration alone buys nothing measurable (CI straddles zero on both
spike-ins). Recall is untouched everywhere; the whole effect is precision, and
the signal carries 81% / 78% / 67% of it. On human the mechanism holds but the
magnitude is ~4x smaller than on spike-ins (spike-in density of wobble-shadow
contests). Per corruption family, the signal is worth most on the truncation
families (-0.96 vs -0.68 elsewhere; grids agree, r=+0.63), i.e. where reads
are most ambiguous -- the direction the NanoCount/degradation literature
predicts. Per-candidate-base normalisation (`sqrt_count_mean_llr`) is a NULL:
SIRV4 +0.08 looked directional, HEYA8 -0.02 killed it.

Building the control was wrong twice, both caught by Codex: permute
STRUCTURES not candidates (the gate reduces to per-structure minima; `.ext`
siblings are one hypothesis), and anchor reconstruction on the assigned
minimum (`old + (new-old)` cost a 1-ulp margin and flipped the bare `>=`
threshold). Three CLI guards refuse a partial control. Literature bound: no
published method does isoform assignment in squiggle space at all.

**Scoring trap (cost a full pass):** `sirv4_grid_runner.py --dataset
gencode_sub` scores against GENOME-WIDE gencode v44 while assemblies cover
chr1/19/21 -> Sn diluted ~30x, every F1 delta reads +-0.00. Rescore with
`experiments/prod_validation/rescore_human_expression_contract.sh` against
`_gencodegrid/expression_contract/_truth/<sample>.gtf`; the check is `full`
landing at 42.75.

## 4t. The attestation exemption is now the real-drna DEFAULT at K=3 (2026-09-11, LANDED 17fec28; decision record fc14e15)

Pre-registered rule and outcome: `docs/PRODUCTION_ATTESTATION_DEFAULT_DECISION.md`.

**The reframing that unblocked it (Codex).** Adding n predictions that recover
t unmatched truths improves F1 whenever `t/n > F1_old/2` -- at human `full`
that is **21.4%** marginal precision, not the current 67.9%. The gate-dropped
suffix>=3 tier (1367 annotated candidates, 77.5% truth) was far above the bar.

**Decision run** (one snapshot `_snap/att_dec` = 735cbf8, one lever, default
re-run on the same snapshot, K=3 + `novel` veto frozen before the arms):

| grid | dSn | dPr | dF1 |
|---|---|---|---|
| human GENCODE (7 conds) | +1.47 | +0.36 | **+1.45** |
| SIRV4 (156 / 6 libs) | +0.18 | -0.15 | +0.12 [-0.05,+0.32] |
| HEYA8 (130 / 5) | +0.04 | -0.03 | +0.04 [+0.01,+0.09] |

`full` 31.2/67.9/42.75 -> **34.0/68.8/45.51**; cells beating stringtie3
**2/7 -> 6/7**; `p00` exactly inert. Precision ROSE: 194 distinct newly
matched truths from 232 new models (84%). The pre-registered failure mode
(one shared 3' suffix credits every compatible annotated extension) is an
ASSOCIATION, not a demonstrated cost: shared-suffix buckets yield 8-15 pp less
(alone 90 / pairs 75 / groups>=3 82) and still clear break-even by 3.5-4x.
My first audit said the opposite and was wrong three ways (query-side `=`
labels, a containment predicate that missed `[A,C]`/`[B,C]` sharing `[C]`, and
non-merging grouping); `attest_sharing_audit.py` is the corrected version.

**Silent configuration trap, fixed and guarded:** a profile key nobody reads
back out of `resolved` in `fin/cli.py` is a NO-OP -- the profile computed 3,
click's default 0 ran, nothing raised, and the only symptom was a run coming
out byte-identical to the arm it should have differed from. Guard is now
behavioural (`tests/unit/test_profiles.py`, through the real click command,
incl. one test iterating every profile key). `real-drna-precision` inherits
K=3 by design, unvalidated in combination; dataclass default stays 0; `sirv`
untouched.

**The K operating curve (2026-09-13, human arms scored today):**

| K | human mean Sn/Pr/F1 | vs K=3 | beats stringtie3 |
|---|---|---|---|
| 0 | 25.17/60.46/35.38 | -1.45 | 2/7 |
| 3 (landed) | 26.64/60.81/36.83 | -- | 6/7 |
| 2 | 27.14/61.04/37.34 | +0.50 (5 up / 1 down) | 6/7 |
| 1 | 27.63/61.16/**37.80** | **+0.97** (5 up / 1 down) | 6/7 |

Monotone in both Sn AND Pr; `full` 34.0/68.8 -> 35.7/69.4 at K=1; the K=3->1
step is 85% marginal precision on `full`. The one cell down is
`c_jitter50_10bp` at -0.05. On SIRV4, K=3/2/1 produce BYTE-IDENTICAL
assemblies (the depth threshold never binds on a high-depth spike-in; only the
veto does), so lowering K carries zero SIRV4 risk. HEYA8 at K=1 is the last
grid the frozen rule requires: SLURM 211963 -> `_heya8grid/attK1.tsv`,
compare vs `_heya8grid/adef.tsv` with `arm_compare.py`. If it is >= 0 in the
mean with no family <= -1.0, K=1 is promotable by the same rule.

## 4u. All-dimension audit, and where the remaining loss lives (2026-09-11, 945fdae)

`docs/PRODUCTION_ALL_DIMENSION_AUDIT_2026-09-11.md`. Every axis, vs the best
alternative: structure recall/precision on annotated human ahead or level
(+1.4 / +5.8); robustness to corrupt (+10.5 jitter) and sparse (+3.6 p50)
annotation far ahead; **de novo `p00` is the one open structural dimension
(-0.30)**, and it is exactly where an annotation-conditioned exemption is
inert; transcript ends: matched models place both ends exactly (median 0,
100% within 25 bp; `gffcompare --strict-match` returns ZERO on real dRNA at
any `-e` -- wrong instrument, not a 0% score); quantification: downstream of
detection (4q); signal: genuinely discriminative (4s).

**Remaining recall loss under the new default** (`remaining_loss_attribution.py`
on the `--selection-audit` run `_gencodegrid/_audit_after_attest`): of 4823
missed expressed multi-exon truths, `m2_support_gate` still owns **81.8%**,
`cluster_recheck` 7.2%, survived-then-lost 9.1%, generation **1.1%**. So on
annotated human data generation is NOT the bottleneck and the next question
is the gate's operating point (4t), not a new mechanism.

**A literature premise tested and mostly not confirmed:** minimap2's default
`--splice-flank=yes` prior (base after GT is A/G) holds on human (93.9%) and is
violated on SIRV4 (62.3%) / HEYA8 (50.9%). But the per-junction placement
penalty is ~3 pp everywhere and SMALLEST on HEYA8 (+0.3), so grids that break
the prior do not pay for it per junction -- only the violating fraction
differs. ~0.2 pp of wobble on human, ~1.1 on SIRV4. Not worth realigning;
use `--splice-flank=no` on any future SIRV alignment. No grid uses
`--junc-bed` (no annotation leak through the aligner).
`experiments/prod_validation/splice_flank_artifact.py`.

**Broad literature review (2026-09-11, agent-run, in the audit doc's sources):**
top implementable ideas, none yet tested -- (1) IsoQuant's LOCUS-RELATIVE floor
(`>=N FSM reads AND >=2% of max locus coverage`) plus its non-polyA
terminal-attachment rule, aimed squarely at `p00`; (2) StringTie3's verbatim
poly(A/T) terminal rules; (3) a degradation-ramp-vs-sharp-peak statistic on
read 5' starts as a SELECTION feature (scATS alpha/beta; Telos
`5'_degradation_score`); (4) publish operating curves and a spike-in-derived
noise floor rather than single points; (5) TUSCO as a fourth, endogenous grid.
Also: dRNA chimeras are adapter-bridged, so an internal-adapter scan
(DeepChopper) is the mechanistically right replacement for the `sc250`
soft-clip heuristic.

## 4v. Fusion: isolated analysis (2d1644e), fused-transcript models (25e5448), and the real-data repair (9a200b7)

A separate line from the assembly campaign above; the only contact is the
pre-existing quant-axis leak noted in 4f. Contracts and records:
`docs/FUSION_IMPLEMENTATION_CONTRACT_2026-09-12.md`,
`docs/FUSION_VALIDATION_2026-09-12.md`, `docs/FUSION_SIGNAL_VALIDATION.md`,
`docs/FUSION_TRANSCRIPT_MODELS_2026-09-15.md`,
`docs/FUSION_REAL_DATA_2026-09-15.md`. READ THOSE before touching fusion.

**The owner's independent analysis (2d1644e).** `--fusion` runs ordinary
assembly unchanged, then a separate `fin/fusion/` pipeline with its own
BAM scan, N-arm read paths, chain support, junction correction, optional
krill signal, and its own atomically-published output namespace. Assembly
byte-identity is asserted by test.

**Fused-transcript models (25e5448 + 992d8be).** The chain layer says which
ordered gene path and breakpoints molecules support; it does not say WHICH
isoform of each gene, cannot extend a 5'-degraded dRNA read to a transcript
start, and emits opaque templates. The transcript layer adds, per model:
per-arm isoform resolution on observed introns (all compatible transcripts
reported; class is display-only), the conventional exon-rank junction name
(`BCR:e14>ABL1:e2`), endpoint completeness kept SEPARATE from splice class,
two FASTAs that never share an id (`fusion_models.fa` evidence-constrained,
`fusion_models.annotated.fa` annotated completion with
`extension_observed=false`), a multi-chromosome exon table with per-exon
provenance, and molecule assignment (unique / compatible / partial / nearby).
Three Codex rounds; the two rejected versions are recorded in the doc
(permuting candidates instead of structures; reconstructing a group minimum
through a shift instead of anchoring on it).

**The real-data repair (9a200b7) -- the important part.** Validated on
synthetic error-free reads, the analysis found NOTHING on the canonical
positive control (BCR-ABL1 in K562 dRNA): 0 chains at defaults; at relaxed
coverage, one wrong chain with the ABL1 arm unlabelled. Four real defects:
an SA CIGAR is a SUMMARY (minimap2 writes one M plus one I/D; the spliced
BCR arm appears as `175M723D` and that fiction beat the truth) -> SA entries
are locators, superseded by any detailed alignment of the same placement; a
12/12 deletion-free anchor rejects nanopore reads with a 2-bp deletion near
a real junction -> windowed anchors with an identity floor proportional to
what aligned; a template carrying 700 bp of intron VALIDATED because the
read never aligned that part -> every template block must be represented by
the same forward alignment; and dRNA's ~80 bp non-genomic tail fails a 0.9
full-read coverage gate -> a per-molecule coverage domain excluding only
terminal bases no detailed alignment could place. Nomination clip 100 -> 250
(at 100 it nominates 8.6% of a real library and trips the cap).

The transcript layer had to follow: the chain layer keys on exact topology,
so five molecules of one fusion that lost different numbers of 5' exons are
five chains of support 1. Models are now aggregate ENVELOPES with
`envelope_status` (`single_molecule_spanned` only when ONE observation
reaches both ends), EXACT membership (a molecule that decisively settled
2 bp away is `nearby`, reported and NOT counted), and merging only when a
molecule is an exact member of both models AND the group's diameter is
bounded pairwise.

Result, frozen as `tests/data/fusion_k562_bcr_abl1/` (95 real records incl.
the `175M723D` SA, 1.3 MB, seconds): ONE model `BCR:e14>ABL1:e2` (K562's
canonical b3a2), `arm_fsm`/`arm_fsm`, 7172 bp, the right BCR and ABL1
transcript ids, 4 molecules of unique support. Five Codex rounds; every
finding has a regression that fails on the previous code. 2417 tests pass.

**Throughput:** mappy holds the GIL in `map()` (`n_threads` only builds the
index), so a full-library run sat at 17k molecules for 7 h.
`FusionGenomeAligner.map_many` forks a pool sharing one index; serial and
multi-worker results are asserted equal per read. Limits are in its
docstring (not reentrant, fork/Linux only, no crash watchdog).

**WHOLE-LIBRARY ACCEPTANCE — DONE.** `docs/FUSION_WHOLE_LIBRARY_2026-09-15.md`.
The entire K562 direct-RNA library (1,074,197 records) against the whole
GRCh38 + GENCODE v44, nothing told where to look: 23,403 candidate molecules →
595 chains → **69 fused-transcript models, 4 of them three-arm (multi-gene)**,
split 48 readthrough / 19 intrachromosomal / **2 interchromosomal**. The top
interchromosomal call by support is **`BCR:e14>ABL1:e2`** — K562's e14a2
driver, `arm_fsm`/`arm_fsm`, 7172 bp, 4 molecules, correct transcript ids, and
its 24 exons across chr22+chr9 sum to exactly the FASTA record length. All six
both-arms-annotated models are the driver or known readthrough/paralog loci
(DHRS4-DHRS4L2, GAGE2A-GAGE1, BAG6-SLC44A4). Run ~1h50 on an idle node.

Two things to tell a downstream user: **chrM is 12 of the 69 and holds the two
highest-support models** (390, 91) because the mitochondrial genome is circular
— filter it; and **one event can occupy several rows** (BAG6-SLC44A4 is 3, two
differing by 6 bp), joinable via the mutual `nearby_read_ids`, because models
merge only on exact shared membership and never on proximity.

Three commits since: the gene-transition nomination had to be tightened
(6.09%→0.19% of records) to make the library tractable, then Codex found the
false negative in the tightening — `{A}->{}->{B}` and `{A}->{A,B}->{B}` cross
an intermediate arm, so the test is over ALL pairs, not adjacent ones. That fix
is what produced 3 of the 4 multi-gene models (66/1 → 69/4). And `_snapshot`
dropped `st_dev`: two runs died at publication on an autofs REMOUNT (`st_dev
144->138`, content identical) — the work dir survives publication failure,
which is the only reason the 69 models above exist.

**PERFORMANCE, 2026-09-16 (`docs/FUSION_DP_AND_CACHE_2026-09-16.md`).** The
whole-library run is **6,631 s -> 1,643 s (4.0x), all 13 artifacts
byte-identical**. Two causes, both found by profiler after two wrong structural
guesses (chain construction is 1 s, measured). (1) `LazyGenomeFasta` returns a
WHOLE chromosome and its LRU held two -- 499 misses, **497 evictions**, 66.6 Gb
fetched from a 3.1 Gb genome; now bounded by SIZE, since a count bound cannot be
chosen without knowing the reference. This is also why the two-locus fixture was
always fast: a two-chromosome reference cannot miss. (2) 29-38% of
`align_junction_window` calls are exact repeats, so memoising the pure
recurrence is an identity. Phase: 667 -> 370 (cache) -> 213 s (memo).

**Do NOT try to replace that DP with mappy.** Codex's kill case: the function
must return an optimal alignment for EVERY valid window, including one with no
matching seeds (`"A"*75` vs `"C"*75` scores -154.0, 75 insertions, reason
large_indel) -- a seed-chain-extend mapper produces no hit at all. And score
equality is not a sufficient gate for any replacement: two EQUALLY optimal
intervals carry anchor counts 0/2 and 1/1, so one rejects and one passes. The
contract is this recurrence plus this canonical traceback, not "affine
semi-global alignment". Going further needs a compiled DP = numba or parasail =
a new dependency and an SIF rebuild; that is an owner decision. Amdahl caps the
prize well below 2x on the phase.

**Open:** adapter-bearing negative controls are now in the suite
(`tests/unit/test_fusion_adapter_chimera.py`) but are synthetic; the two-locus
fixture still proves the regression, not genome-wide mapping uniqueness; no
benchmark against a fusion caller has been run. One loose end: the per-molecule
adjudication is ~6,277 s of the ~6,600 s run (chain construction is 1 s — it
was NOT the bottleneck, measured), but two progress probes disagree about where
inside it the time goes and I did not resolve which is wrong. Profiling that is
operational, not correctness.

## 5. Validation protocol for ANY future change (frozen invariants)

1. Opt-in flag, default off, no profile enables it; explicit misuse fails
   loudly (see `PipelineConfig.validate()` precedents).
2. Default-path byte identity: rerun `full`+`p50` (r2r1) with default flags,
   `cmp` against `_sirv4grid/base/` outputs (assembly.gtf, scores.tsv,
   abundance_refit.json) — must be identical. `p00` must stay identical too.
3. Audit modes must be byte-inert (max per-cell dF1 = 0.000000 on both
   contracts).
4. Determinism: two independent reruns of >= 2 cells byte-identical including
   the feature's artifact TSV/JSON.
5. Mass conservation (`abundance_refit.json.mass_balance_error` ~1e-12), zero
   zero-abundance transcript rows.
6. Fresh full 156-cell run + rescore for any behavioral claim; paired
   cell-level bootstrap (10k, seed 20260831); report ups/downs, per-condition
   means, top/tied; compare against the CURRENT best mode, not just `off`.
7. Full test suite + M3 regression must pass.
8. Leakage red lines: NanoCount/gffcompare/truth GTFs/condition names are
   evaluation-only; never inputs or features. No condition-name branching.
9. All thresholds are six-H9-sample development fits; bootstrap CIs after
   tuning are descriptive.
10. CROSS-DATASET (2026-09-04, supersedes the old "data-gated" note): SIRV4 is
   the DEVELOPMENT set, Sequin is the HELD-OUT generalisation set. Nothing may
   be tuned on Sequin. Bootstrap each dataset SEPARATELY (never pooled). A
   change "transfers" iff the direction agrees on both and, on Sequin, neither
   Sn nor Pr degrades (paired bootstrap 95% CI lower bound > -0.5) and complete
   Spearman does not drop.

## 6. Git state (IMPORTANT; refreshed 2026-09-13)

Branch `dev` is AHEAD of `origin/dev`. Since 2d1644e, all Codex-reviewed
(gpt-6-astra) before each commit: 25e5448 + 992d8be (fused-transcript models),
c0c83a6 (attestation K=1), 9a200b7 + 951ee6c (real-data fusion repair),
d10c9c3 (nomination tightening), 9486814 (Codex's cross-arm transition fix),
d6adabf (adjudication progress log), a11e55a (adapter-chimera control),
6da527e (autofs remount / `SNAPSHOT_FIELDS`), 87689bd (genome cache bounded by
size), 040d8a0 (junction-DP memoisation), 6c95ea6 (clock starts before the GTF
parse). Suite on this tree: **2460 passed, 4 skipped**, and the suite itself
dropped from 273 s to 66 s because the memo hits in tests too.

Historical: `2d1644e` was (owner's fusion commit on top of the
2026-09-10/11 campaign). Everything in 4s-4u is COMMITTED and Codex-reviewed
(gpt-6-astra) before each commit: d4e63fc (scramble control), a15fdb4 /
b8fe84c / e389faf / 735cbf8 (signal decomposition records), 90ca80a
(attestation pre-registration), 17fec28 (attestation default), fc14e15
(decision record), 945fdae (all-dimension audit). Test suite on that tree:
2203 passed, 4 skipped (`--ignore=tests/unit/test_junction_dtw_variants.py`,
which needs g++). Run tests inside the SIF with `--nv`.

Tracked-modified: only `HANDOFF.md` (this file; pure additions -- the owner's
raw-DTW notebook block at the top plus 4b5/4o-4v). Untracked and DELIBERATELY
NOT STAGED (owner's instruction, do not commit without explicit approval):
the `docs/PRODUCTION_*DTW*/NEIGHBOR*/HMM*/DISCORDANT*/ALL_READS*/QUESTIONS*/
SIGNAL_WORKED*/REPO_REVIEW*.md` batch and its scripts, `SIRV4_GRID_STUDY.md`,
`SIRV4_IMPROVEMENT_PLAN.md`, `CLAUDE.md`/`AGENTS.md` and the other agent
config files, plus ~600 run artefacts under `experiments/prod_validation/_*`
(exclude forever: grid outputs, SIF images, FASTA/FASTQ, `_img/`).

Source snapshots used by arms (`experiments/prod_validation/_snap/`):
`m2off`=fecd197, `m2scram`=fecd197+control, `att_dec`=735cbf8,
`att_land`=17fec28-equivalent (used by the K-curve and `_audit_after_attest`).
Always verify `run_manifest.json.source_root` and one `source_sha256` per arm.

## 7. Next work, in priority order (refreshed 2026-09-13)

0. **OWNER GOAL (2026-09-13) -- DELIVERED AND ACCEPTED ON A WHOLE REAL LIBRARY;
   see 4v and `docs/FUSION_WHOLE_LIBRARY_2026-09-15.md`.** Fused-transcript
   level multi-gene output plus a final fusion-transcript FASTA:
   `fusion_models.tsv` / `.fa` / `.annotated.fa` / `_exons.tsv` /
   `_evidence.jsonl` / `_collapsed.jsonl`, landed 25e5448, repaired for real
   direct RNA in 9a200b7, and validated on the ENTIRE K562 library against the
   WHOLE genome: 69 models, 4 multi-gene (three-arm), and `BCR:e14>ABL1:e2` as
   the top interchromosomal call. Adapter-bearing negative controls are in the
   suite (a11e55a). REMAINING on this axis: the negative controls are
   synthetic, and no head-to-head against a fusion caller has been run --
   neither blocks downstream use of the FASTA.
1. **Attestation operating point: DONE (c0c83a6).** HEYA8 K=1 came in at
   +0.07 F1 vs default, CI [+0.02,+0.09], worst family -0.01, so `real-drna`
   moved K=3 -> K=1 under the same frozen rule. Next: re-run
   `remaining_loss_attribution.py` under K=1 to see what the gate still owns.
2. **De novo (`p00`), the one open structural dimension (4u):** first
   candidate is IsoQuant's locus-relative floor (`>=N FSM reads AND >=2% of
   max locus coverage`) + non-polyA terminal-attachment rule, as an opt-in
   with the frozen protocol; earlier global-floor attempts gutted recall.
3. **Cheap literature levers not yet tested (4u):** StringTie3 poly(A/T)
   terminal rules; degradation-ramp-vs-peak 5' statistic as a selection
   feature; operating curves + spike-in noise floor in the benchmark report;
   internal-adapter scan to replace `sc250`.
4. **Real-data replication of 4s on a second human library** (the `.ext`
   sibling density makes structure-level grouping matter more on human).
5. Historical items (sidecar trunc/merge recall, polyA wiring, election on
   p00) remain in `SIRV4_IMPROVEMENT_PLAN.md`; they were SIRV-profile levers
   and are superseded by the real-drna default for production purposes.

Cluster-scoped raw-DTW correction, the owner's proposal (2026-09-09):
`docs/PRODUCTION_CLUSTER_DTW_REPLAY_2026-09-09.md` -- a masked replay of the
frozen benchmark restricting alternatives to candidate-level junctions from the
read's interval (offline production discovery + canonical gate), zero new
eventalign, protocol frozen by Codex before scoring. Kill test: 89.3% of the
8,988 SIRV4 raw damages land on junctions that ARE candidate-level, so the
restriction can remove at most 9.5% of the harm. Restricted raw SSE: SIRV4
73.834% (+651/-8,056), -19.311 pp vs BAM -- the frozen stop rule fires and the
proposal is falsified in its stated form; no signal arm passes the gate (best,
dwell5 + 5% abstention, is -3.912 on SIRV4 with all six libraries negative).
The only positive arm is the BAM-support-only count rule (+3.390 / +1.739 pp,
zero catalog->off-catalog), and an unrestricted-registry count arm does as
well or better (96.644%), so candidate restriction contributes nothing there
either; that rule is production's existing junction_snap applied at read level
and shows 439 catalog->catalog moves that the metric cannot judge. Codex
independently recomputed every count. Production unchanged; nothing staged.

