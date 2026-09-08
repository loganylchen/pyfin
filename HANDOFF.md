# PyFIN SIRV4 Optimization — Session Handoff

Last updated: 2026-09-02. This document is the authoritative handoff for
continuing the SIRV4 expression-performance work. Read it TOGETHER with
`SIRV4_IMPROVEMENT_PLAN.md` (full evidence + acceptance gates per lever) and
`SIRV4_GRID_STUDY.md` (benchmark semantics). Treat live source and generated
artifacts as authoritative over any narrative.

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

## 2. Current result ladder (expression contract, all live 156-cell runs)

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

## 4n. Should the junction snap use raw signal? (2026-09-08, audit RUNNING)

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

## 6. Git state (IMPORTANT)

Branch `dev`; origin/dev = f628f58. On 2026-09-02 the batch was split by
a Codex review (gpt-5.6-sol; full output in
`.omc/codex_review_sirv4_batch_2026-09-02.txt`):

- **Committed as `ca82319`** (Codex-clean, 1214 tests + M3 regression
  passed on exactly that tree): P0 strict-positive GTF gate, complete-chain
  rescue, catalog completion, and their config/CLI/README wiring + tests.
- **Still uncommitted (Codex found 3 issues, fixes pending)**: the
  junction-election and containment-evidence levers — tracked modified
  `README.md`, `fin/cli.py`, `fin/pipeline/config.py`,
  `fin/pipeline/runner.py`, `fin/pipeline/selection.py`,
  `fin/analysis/quantification.py`, `fin/candidates/dataclasses.py`,
  `fin/candidates/discovery.py` (election/containment hunks only) and new
  `fin/analysis/containment_evidence.py`, `fin/analysis/junction_election.py`,
  `tests/unit/test_containment_evidence.py`,
  `tests/unit/test_junction_election.py`. Findings: (P1) election
  electorate built from `reps[j]` only, omitting reads of other carriers of
  the same junction (defense burden then drops well-supported junctions on
  `no_quorum`); (P2) `finalize_containment` merges multi-row verdicts
  last-row-wins instead of by precedence, losing `supported` exemptions;
  (P2) `default_signal_scorer` rebuilds a krill aligner per contested
  interval. Must be fixed and re-reviewed before commit.
- **Uncommitted 2026-09-03 batch (NOT yet Codex-reviewed; user must run
  `/codex:review` before any commit)** — three generation/quant changes plus
  the additive junction correction:
  1. GTF multi-exon chains are ORDINARY zero-read variant keys in
     `cluster_families` (default, no flag): `fin/candidates/chain_cluster.py`
     (`gtf_variants=`, `_attach_gtf_variants` removed; annotated keys are
     legal containers but NEVER foldees — Codex's "pure-GTF foldee folds
     unconditionally" was implemented first and measured at -3.03 F1 /
     Sn -3.8 on the 156-cell plain default, run kept as
     `_sirv4grid/gv_off_foldgtf`), `fin/candidates/discovery.py`
     (`gtf_models` built before clustering, every GTF model of a surviving key
     emitted with annotated geometry, `emitted_gtf`), tests rewritten in
     `tests/unit/test_cluster_families.py::TestGtfVariantKeys`,
     `test_collapse.py`, `test_discovery_families.py`. **Guided-condition
     outputs change by design** — `_sirv4grid/base` is no longer the
     byte-identity reference for full/p50; p00 unaffected by this item.
  2. CIGAR-level collapse span_guard: `chain_cluster.read_alignment_from_cigar`
     + `_cigar_exonically_covers` (W=[d-10,a+10): footprint covers W, no N,
     D<=3 each / <=6 total, >=8/10 matched per flank, body matched >= len-6);
     `collapse(read_alignments=)` replaces `read_spans=`; discovery records
     alignments in its CIGAR walk. **Changes default outputs in every
     condition** (fold decisions). Figure:
     `_sirv4grid/fig/span_guard_explainer.png` (`span_guard_figure.py`).
  3. `--folded-tss-recovery-mode off|audit|require` (default off):
     `fin/analysis/tss_shadow_recovery.py` (post-EM shadow test via
     `evaluate_internal_tss`, `tss_shadow_recovery.tsv` shards merged at
     finalize), `_quant_m2_em` two-pass wiring, and the frozen read-5'
     positional likelihood `assignment.apply_start_likelihood`
     (`Assigner.assign(start_likelihood_context=)`, second pass only,
     exact {parent, recovered} ties; mirrored into `nlls_by_read` for the
     support gate). `CandidateSet.indels_by_read` added.
     Tests: `tests/unit/test_tss_shadow_recovery.py` (9). Designs:
     `.omc/codex_consult_spanguard_gtf_lost_2026-09-03.txt`,
     `.omc/codex_consult_tss_5p_likelihood_2026-09-03.md`.
  4. Additive per-read junction correction (`fin/analysis/junction_correction.py`,
     `junction_support.py`, config/cli `--junction-correction-*`, runner
     ballots, always-on `junction_support.tsv`; tests
     `test_junction_correction.py` (20), `test_junction_support_record.py`).
     Candidate-level probe grid: `_sirv4grid/jc_additive_m{M}_r{R}_k{K}/`
     (`jc_additive_grid*.sbatch`, `junction_correction_candidate_probe.py`).
  Additive junction-correction candidate probe (96 cells, done): lost
  chains 0 everywhere, but new candidates are ~90% FP (m10/K2: 12 TP /
  116 FP; m20/K2: 0 / 64; read moves own_TP->FP 612 vs own_FP->TP 448)
  -> the per-read K-way ballot is NOT a viable completion mechanism
  (`_sirv4grid/jc_additive_probe_summary.txt`).
  Validation running as SLURM 211628 (`experiments/prod_validation/gv_grid.sbatch`):
  156-cell `gv_off` (plain defaults, vs `base`),
  `gv_best` (P0+rescue+defense, vs `je_defense` 91.87), `gv_tss_audit` (must
  cmp-equal `gv_best`), `gv_tss_req`, determinism rerun on 3 cells, then
  `sirv4_expression_rescore.py` (modes registered). Full suite: 1319 passed,
  4 skipped (SLURM 211622, before the never-foldee change; rerun 211629).
- Also uncommitted, not yet reviewed: `SIRV4_GRID_STUDY.md`,
  `SIRV4_IMPROVEMENT_PLAN.md`, `HANDOFF.md`, and the
  `experiments/prod_validation/sirv4_*.py` + probe scripts. Exclude
  forever: all `_sirv4grid/` outputs, SIF images, FASTA/FASTQ, `_img/`.
  `.gitignore` has a pre-existing user-owned modification — leave it.

## 7. Next work, in priority order

1. **Trunc/merge recall (the remaining 44 loss cells)** — P2 sidecar:
   preserve annotation-independent read-derived candidates through joint EM
   (do not let GTF absorb/outcompete them before evidence), then contrastive
   reconciliation + read-supported endpoint repair. Evidence: 42-64% of
   residual trunc/merge FN chains already exist in the same sample's p00
   output; remaining FNs die at MULTIPLE stages (no single gate ablation
   recovers them). Gates frozen in the plan: +5 F1 per `c_trunc*50`, no
   p00-p50 regression, recover >= 50% of p00-recoverable FN chains.
2. **PolyA wiring (P4)**: polyA-qualified TES evidence into the containment
   own_tes rung + a mono-exon span route (mono corrupt entries: 0/8 have
   polyA at claimed TES vs reads reading through; genuine mono TPs often lack
   polyA too — tri-state only, never absence=rejection). Machinery exists
   (`fin/scoring/polya.py` krill pass, `evaluate_tes_support(polya_supported=)`).
3. **Junction election on `p00`/de novo cells**: the mechanism is
   annotation-free but the grid applied it via `--extra-with-gtf` only;
   extending to unannotated runs is unexplored upside.
4. **Real-data re-validation**: defense burden and all screen thresholds on
   real dRNA (genuine nearby splice sites, lower coverage) — prosecution is
   the safe default until then. Also rescan external tree for the 24 missing
   competitor cells and refresh standings.
5. Optional consolidation: retire superseded probe scripts, promote the
   election/containment audits into `candidate_evidence.tsv` schema.
