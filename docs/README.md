# PyFIN documentation index

The repository carries several documents that each read as authoritative
because each *was*, on the day it was written. This index says which one to
open for which question, and which are historical.

Four destinations:

## 1. Using the tool

| Question | Read |
|---|---|
| Install, including the signal backend | root `README.md`, Installation |
| Run it, and what the profiles mean | root `README.md`, Quick Start |
| Prepare inputs, diagnose a run, read the outputs | [`USAGE.md`](USAGE.md) |
| Every CLI flag | `pyfin --help`, and the parameter table in root `README.md` |

Two things the parameter table does not yet do: it does not separate raw
defaults from profile overrides, and it omits some settings the shipped
`real-drna` profile turns on. `pyfin_out/run_manifest.json` records the
**resolved** configuration for a given run and is the authority when the two
disagree.

## 2. What the algorithm currently does

[`CURRENT_ALGORITHM.md`](CURRENT_ALGORITHM.md) -- the shipped `real-drna` path
stage by stage, plus the fusion acceptance rules and what is off by default.
It contains no benchmark results by design.

For the reasoning behind particular defaults, the dated production notes:

- [`PRODUCTION_ATTESTATION_DEFAULT_DECISION.md`](PRODUCTION_ATTESTATION_DEFAULT_DECISION.md) -- why the support gate exempts attested annotated candidates by default
- [`PRODUCTION_QUANT_ACCURACY_2026-09-10.md`](PRODUCTION_QUANT_ACCURACY_2026-09-10.md) -- the quantification axis and the TPM reporting convention
- [`PRODUCTION_M2_GUIDED_PATH_2026-09-10.md`](PRODUCTION_M2_GUIDED_PATH_2026-09-10.md) -- what the M2 tie-break contributes
- [`PRODUCTION_ALL_DIMENSION_AUDIT_2026-09-11.md`](PRODUCTION_ALL_DIMENSION_AUDIT_2026-09-11.md) -- the audit that set the current priorities

## 3. Fusion

Fusion is an independent analysis run after an unchanged assembly (`--fusion`).
It is **opt-in** and has not been benchmarked against any fusion caller.

| Question | Read |
|---|---|
| How the whole thing works, at length, in Chinese | [`FUSION_EXPLAINED_zh_2026-09-18.md`](FUSION_EXPLAINED_zh_2026-09-18.md) |
| The contract the implementation must satisfy | [`FUSION_IMPLEMENTATION_CONTRACT_2026-09-12.md`](FUSION_IMPLEMENTATION_CONTRACT_2026-09-12.md) |
| The fused-transcript model layer and its outputs | [`FUSION_TRANSCRIPT_MODELS_2026-09-15.md`](FUSION_TRANSCRIPT_MODELS_2026-09-15.md) |
| What a whole real library produced | [`FUSION_WHOLE_LIBRARY_2026-09-15.md`](FUSION_WHOLE_LIBRARY_2026-09-15.md) |
| Signal validation of fusion junctions | [`FUSION_SIGNAL_VALIDATION.md`](FUSION_SIGNAL_VALIDATION.md) |

## 4. Historical, kept for the reasoning

These describe a state the code has moved past. They carry banners saying so.
Do not quote a threshold from them as current behaviour.

| Document | What it was |
|---|---|
| `../ALGORITHM.md` | full algorithm description, snapshot 2026-08-26 |
| `../PIPELINE_AS_IMPLEMENTED.md` | step-by-step walkthrough, 2026-09-03, **SIRV** profile as its "default" |
| `PIPELINE_WALKTHROUGH.md` | production walkthrough, 2026-09-08, Chinese; superseded by `CURRENT_ALGORITHM.md` but its worked examples remain useful |
| `../REPO_MAP.md` | repository layout, snapshot at `97029eae`; some entries are stale |
| `../ALGORITHM_REVIEW.md`, `../CANDIDATE_EXISTENCE_DESIGN.md`, `../TSS_CONTAINMENT_STUDY.md`, `../PROFILE_OPTIMIZATION.md` | design and study notes from particular investigations |
| `../HANDOFF.md` | a working session handoff, not a user document |
| `experiments/**/SUMMARY.md` | per-experiment results; scientific provenance for specific claims |

## Known documentation gaps

Recorded rather than implied, so a reader is not misled by their absence:

- **No clean-checkout reproduction recipe.** The runners under
  `experiments/prod_validation/` read dataset trees configured for the
  authors' machine. Reproducing the published comparison needs an input
  manifest with checksums, the annotation-condition generation step,
  competitor commands and versions, and the table-generation step. None of
  that is written yet.
- **No end-to-end example that exercises signal.** `tests/testdata/` is
  referenced by the root README and `tests/README.md` but no files under it
  are tracked. The tracked K562 fusion fixture
  (`tests/data/fusion_k562_bcr_abl1/`) is sequence-only.
- **The figures in `diagrams/README.md` are not in the tree** -- that
  directory holds the generation prompts, not the PNGs it lists.
- **No effective-configuration reference** separating raw defaults, profile
  overrides, inert mode-specific options and experimental switches.
  `CURRENT_ALGORITHM.md` §0 and §8 cover the production profile and the
  disabled options, which is most of what such a reference would say, but it
  is not exhaustive over every flag.
