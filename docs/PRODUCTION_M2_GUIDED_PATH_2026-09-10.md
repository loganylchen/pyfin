# The signal tie-break on the guided path: is it earning its place?

Pre-registered before either arm finished. Written because the literature
review's top signal item (NanoSplicer S1.7 — normalise the junction score per
candidate base, turn it into a posterior over the closed candidate set, and
gate it with an explicit abstention threshold) turned out to describe machinery
this codebase already has, and the interesting question became which of that
machinery production actually runs.

## What production runs

**Correction, made after the arms had already run.** This section first claimed
that every production number came from the wide-window per-event `mean` with no
abstention, on the reasoning that `resolve_m2_metric` maps the sentinel `auto`
to `mean` whenever a usable GTF is present, and all four grids are guided. That
inference was wrong, and the run manifests say so: the `real-drna` profile — the
one every grid arm uses — **pins `m2_metric = "summed_llr"` explicitly**, so the
route is `fixed`, not `auto-guided`, and `auto` never runs. The wide-window
`mean` describes the older `sirv` profile, which is where the `mean` manifests I
first sampled came from.

What production actually runs is therefore:

| | production (`real-drna`) | `sirv` profile |
| --- | --- | --- |
| metric | `summed_llr` (route `fixed`) | `auto` → `mean` |
| window | tight, `±m2_summed_llr_flank` (8) around the DIFFERING junction boundaries | wide, `±junction_k` (10 tx-bp) |
| reduction | undivided sum | per-event mean |
| threshold | `m2_summed_llr_margin` = **1.0** | `m2_diff_cover_margin` = 0.5 |

That matters for the literature question, and in production's favour: the
tight discriminating window and a real abstention threshold are two of the
three things NanoSplicer S1.7 prescribes. The one production does *not* do is
normalise per candidate base — `summed_llr` is an undivided sum, which is the
form with the known length bias. `sqrt_count_mean_llr` is the normalised
variant and exists, but is never auto-routed.

Note also that `m2_tiebreak_margin = 1e-09`, which I originally read as "no
abstention", is **not the threshold on this path at all**. Under
`quant_mode=m2_em` the diff-cover gate picks `m2_summed_llr_margin` for the
tight metrics and `m2_diff_cover_margin` otherwise; `m2_tiebreak_margin`
belongs to the `argmax` mode's separate resolver
(`_quant_argmax_keep` → `m2_resolve_tie`), which production does not use.

One gap survives the correction: **there is still no observability on either
path's decisions beyond a per-interval count.** `_log_contrast_stats` does fire
for `summed_llr`, but the per-comparison JSONL is behind
`m2_contrast_stats_jsonl`, which is off in every arm, so the margin
distribution production actually operates at has not been recorded.

## Result: the tie-break is load-bearing, and what it buys is precision

`m2off` vs `m2base`, 156 paired cells over 6 libraries, library-block bootstrap:

| | mean Δ |
| --- | --- |
| Sn | **−0.16** |
| Pr | **−3.04** |
| F1 | −0.93, 95% CI [−1.46, −0.49] |

13 cells up, 136 down, 7 unchanged. Every corruption family is negative, worst
`c_trunc3p` (−1.50) and `c_trunc5p` (−1.39); `c_jitter` is the mildest (−0.33).

The shape is the finding. Recall barely moves because the tie-break removes no
candidate; the whole effect is precision. The mechanism is concentration: M2
gives a tied read's full mass to one structure, and the spurious wobble
siblings that would otherwise have received a share of it fall below the
abundance floor. Per the decision rule below, this is the "clearly worse"
branch — the tie-break is earning its place, and the follow-up is calibration,
not removal.

**What this does NOT establish.** It does not show that the *signal* is right.
`off` changes two things at once: it removes the signal and it restores M1's
fractional 1/K split. An arbitrary-but-decisive tie-break would concentrate
mass just as well. Separating those is what the control arm below is for.

## The tie-break is live, so the test is not vacuous

Aggregating the per-interval `ties=/refined=` log line over ~60 guided SIRV4
cells of an existing run:

| quantity | value |
| --- | --- |
| intervals | 1447 |
| M1 ties presented to M2 | 14002 |
| ties the signal actually decided | 3937 (28.1%) |
| intervals with at least one refinement | 260 (18.0%) |

So roughly one in four ties is resolved by signal on the guided path. Whether
those 3937 decisions are better than the M1 fallback split is the open
question.

## The arms

Both run the full 156-cell SIRV4 grid from **one** source snapshot
(`_snap/m2off`, `fecd197`) under the `real-drna` profile, differing in a single
lever, per the one-lever rule.

| arm | lever | job |
| --- | --- | --- |
| `m2base` | none (production `real-drna`: `summed_llr`, margin 1.0) | 211934 |
| `m2off` | `--m2-metric off` | 211935 |
| `m2scram` | `--m2-scramble-seed 20260831` (signal-free control) | see below |

`off` does not merely change the score — it returns the raw AS tie sets and
skips signal entirely, so the read keeps M1's default split. The contrast is
therefore "signal decides the tie" versus "signal does not participate".

### The signal-free control, and the version of it that was wrong

`m2scram` keeps the concentration and removes only the signal, so
`m2base − m2scram` isolates the signal content that `m2base − m2off` cannot.

The first implementation permuted the finalized scores across the **candidates**
of a contrast, preserving the candidate-level value multiset. Codex rejected it
as a confounded control, correctly, and the mechanism is worth recording because
it is not obvious: the diff-cover gate does not read candidate values. It first
reduces them to **per-structure minima** (`by_group[k] = min(...)`, keyed by
intron chain) and thresholds the margin between those minima. Two candidates
that share a chain — an annotated model and its `.ext` read-extent sibling — are
one hypothesis to it. So with values `A=1, A_extent=1, B=3` (structure margin
2.0, a hard assignment), the candidate-level permutation that sends `B→1` and
`A→3` leaves the multiset `[1,1,3]` intact while collapsing the structure margin
to 0.0 — flipping the read from hard-assigned to abstaining, and changing the
covered-vote prior that redistributes *other* reads. The "control" would have
been a second lever.

The corrected control permutes **structures**, using the gate's own
`structure_key` (lifted from a closure to a module-level function so the two
cannot drift). Each member keeps its offset from its own structure's minimum, so
the scored column set, the structure partition, the multiset of structure
minima, the thresholded margin, and every abstention are all preserved exactly;
only which structure wins changes. Tests assert the gate's margin and its
hard/abstain classification over 40 seeds on the duplicate-structure shape,
including Codex's counter-example seed.

Two further review points are handled rather than argued: `--m2-scramble-seed`
now **fails loudly** outside `quant_mode=m2_em`, because the `argmax` mode
reaches signal through `m2_resolve_tie`, which this control does not touch and
would silently half-scramble; and the "byte-identical when off" claim is scoped
to the assembly, verified by rerunning SIRV4 H9 rep4 run2 `full` on the new code
with the default and matching the stored arm's `assembly.gtf` sha256 exactly
(`run_manifest.json` does change — it serialises the config, which now has one
more field).

## Result: the benefit is the signal, not the concentration

All three arms, 156 paired cells over 6 libraries, one source snapshot,
library-block bootstrap (seed 20260831):

| contrast | what it removes | dSn | dPr | dF1 | 95% CI | cells up/down |
| --- | --- | ---: | ---: | ---: | --- | --- |
| `off` − `base` | signal **and** concentration | −0.16 | −3.04 | **−0.93** | [−1.46, −0.49] | 13 / 136 |
| `scram` − `base` | signal only | −0.07 | −2.83 | **−0.75** | [−1.05, −0.45] | 6 / 142 |
| `scram` − `off` | *adds* concentration, no signal | +0.08 | +0.21 | **+0.18** | [−0.27, +0.60] | 72 / 48 |

The third row is the answer. Concentrating a tied read onto an **arbitrary**
structure, with the signal carrying no information, buys **+0.18 F1 with a
confidence interval that straddles zero** — nothing measurable, and the
cell-level split (72 up, 48 down, 36 unchanged) is close to a coin flip.
Meanwhile removing only the signal, while keeping every read concentrated onto
exactly one structure, costs −0.75 of the −0.93 the whole tie-break is worth.

So the mechanism I proposed after the two-arm result — "M2 helps by
concentrating mass" — is **wrong**. It predicted that an arbitrary tie-break
would recover most of the benefit, and it recovers none of it. The tie-break's
value is the signal picking the *right* structure. Concentration is the
delivery mechanism, not the source.

That is the first direct evidence in this project that the nanopore signal is
doing real discriminative work on the production path, as opposed to being a
tie-breaking ritual whose benefit could have come from any deterministic rule.
It is also the reason the earlier raw-DTW failures did not sink the M2 channel:
the channel survives because it scores *inside* an M1 tie, on the tight
differing-junction window, with an abstention threshold — and within that niche
its content is what matters.

### It replicates on an independent grid

The same three arms on HEYA8 (130 cells, 5 libraries, a different spike-in mix,
different libraries, same snapshot and protocol):

| contrast | SIRV4 dF1 [95% CI] | HEYA8 dF1 [95% CI] |
| --- | --- | --- |
| `off` − `base` | −0.93 [−1.46, −0.49] | −0.76 [−0.96, −0.56] |
| `scram` − `base` | −0.75 [−1.05, −0.45] | −0.59 [−0.90, −0.32] |
| `scram` − `off` | **+0.18 [−0.27, +0.60]** | **+0.17 [−0.23, +0.50]** |

The concentration-only effect is +0.18 on one grid and +0.17 on the other, with
the interval straddling zero both times. The share of the tie-break's value that
disappears when only the signal is removed is 81% on SIRV4 and 78% on HEYA8.
HEYA8 makes the precision mechanism even plainer: `off` costs 5.14 points of
precision there for **+0.11** recall.

Caveats that still bound the claim. Both grids are spike-ins. `.ext`
read-extent siblings are far more common on human than here, so the
structure-level grouping carries more weight there and the human behaviour is
not implied — that is the replication still worth running. And the control
measures the **total** effect of structural junction preference inside the
production pipeline, endpoint-channel interactions included, not an
endpoint-free junction-only effect.

## The remaining third of the prescription: per-base normalisation

Production uses the tight window and a real abstention threshold, but an
**undivided** sum — the form with the known longer-template bias.
`sqrt_count_mean_llr` is the per-candidate-base normalised variant, and it had
never been benchmarked on the production path. One lever vs `m2base`, full
SIRV4 grid:

dSn +0.00, **dPr +0.32**, dF1 +0.08, 95% CI [−0.04, +0.21]; 60 cells up, 27
down, 69 unchanged.

Not significant — the interval includes zero. But the direction is consistent
in a way that a null effect usually is not: every corruption family has a
non-negative mean, recall does not move at all, and the gain is entirely
precision. That is the profile of a real but small effect, and it is the
opposite sign from this project's earlier raw-DTW finding that dividing by
length *worsens* things (that was a different scorer on a different task —
whole-template DTW, not a tight in-tie contrast — so the two are not in
conflict). A HEYA8 replication is running; nothing should change on the
strength of one grid and a CI that includes zero.

## Decision, fixed in advance

Report Sn and Pr for every cell, never F1 alone, and pair by cell with a
library-block bootstrap (seed 20260831).

* **`m2off` clearly worse** → the guided tie-break is earning its place as
  configured; the follow-up is calibration (an abstention threshold above
  1e-09), not removal.
* **`m2off` indistinguishable** → 3937 signal decisions per grid are buying
  nothing measurable. The honest conclusion is that the guided path's signal
  contribution is currently decorative, and the tight-window normalised metric
  becomes worth testing on its own merits rather than as a refinement.
* **`m2off` clearly better** → the wide-window no-abstention configuration is
  actively harmful, which would make the abstention threshold the first thing
  to fix, consistent with this project's own earlier finding that a calibrated
  abstention-only tie-breaker was the component with positive evidence while a
  hard signal override was not.

No default changes on the strength of one grid. SIRV4 is a spike-in with four
libraries; anything that survives here has to replicate on HEYA8, Sequin and
human before it is a candidate for a profile default, and a threshold chosen by
looking at these results would be tuning on the test set.
