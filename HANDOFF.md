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

Spike-in validation of the exemption is the gate for defaulting it.

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
