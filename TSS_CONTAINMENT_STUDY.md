# Does the short transcript really exist? A data-driven TSS study

**Question.** Many reads align to a long transcript. A shorter transcript is
fully contained inside it. Is that short transcript real, or is its apparent
support just 5'-degraded reads of the long one? Can this be decided **from the
data alone**?

**Short answer.** Partly, and the honest system must say which part. There is a
sharp identifiability boundary, it can be located, and on the tractable side a
degradation-null hazard test detects real short transcripts with a measured
**0.000 false-positive rate** while **abstaining** rather than guessing when the
data cannot decide.

---

## 1. Why this is hard (and sometimes impossible)

Nanopore direct RNA (dRNA) sequences **3' → 5'**: adapter, poly(A) tail, then
the body. Consequences:

* the **3' end is reliable** — essentially every read has it, and poly(A) can
  confirm it;
* the **5' end is unreliable** — RNA degradation in vivo/in vitro plus early
  read termination truncate it.

The literature states this explicitly: with conventional dRNA "the 5′ and 3′
ends of poly(A) RNA cannot be identified unambiguously … due in part to RNA
degradation in vivo and in vitro that can obscure transcription start sites"
([Yan et al. 2022, *Genome Res* 32:162](https://genome.cshlp.org/content/32/1/162.long)).

So for a contained short transcript S inside long transcript L:

> If S shares L's splice chain **and** L's TES, and differs only by an internal
> TSS, then a read of S and a read of L that happened to stop at S's start are
> **the same observation**. If the degradation process is allowed an arbitrary
> profile, the two explanations are **mathematically unidentifiable**.

A system that always answers yes/no here is overclaiming. The implementation
therefore returns one of `supported` / `unsupported` / **`unidentifiable`**.

### The identifiability ladder

| Rung | S differs from L by | Discriminating evidence | Difficulty |
| --- | --- | --- | --- |
| 1 `own_junction` | a junction/exon L lacks | junction-specific reads settle it directly | easy |
| 2 `own_tes` | its own 3' end | dRNA's **reliable** end: a distinct 3' cluster, poly(A)-confirmable | tractable |
| 3 `tss_only` | internal TSS only | degradation-null probability test only | hard; abstention is a legitimate outcome |

**The rung decides which test runs.** This is not bookkeeping — applying the
5' hazard test to a rung-2 candidate asks the wrong question. Measured case:
SIRV303 is "contained" in SIRV301, but its start maps to offset **19** of the
parent, i.e. the two share a TSS and differ only at the 3' end. There is no
internal TSS to detect. Such candidates are routed to `evaluate_tes_support`,
which asks whether reads terminate at the candidate's own 3' end — evidence
that degradation cannot fake, because reads do not lose their 3' end.
Routing by rung moved SIRV303 and SIRV506 from abstention to correct
`supported` calls.

---

## 1a. How existing tools handle this, and what is missing

No mainstream long-read assembler decides *probabilistically* whether a
contained short transcript is real. They apply structural rules, require
orthogonal evidence, or defer the question:

| Tool | Handling of contained / incomplete models | Limitation for this question |
| --- | --- | --- |
| **IsoQuant** ([Prjibelski 2023, *Nat Biotechnol* 41:915](https://www.nature.com/articles/s41587-022-01565-y)) | Treats both ends as reliable only under `--fl_data`; otherwise 5'/3' truncation is modelled as a nuisance and truncated models are conservatively suppressed | No per-candidate test; a real short isoform and a truncation are handled by the same blanket rule |
| **Bambu** ([Chen 2023, *Nat Methods* 20:1187](https://doi.org/10.1038/s41592-023-01908-w)) | Machine-learning **NDR** (novel discovery rate) trained on annotation-derived labels, using full-length and equivalent-class read support | NDR is a global ranking calibrated to annotation, not a hypothesis test against a *degradation* null |
| **StringTie3** ([Shumate 2025, *Nat Methods*](https://doi.org/10.1038/s41592-026-03080-3)) | New **nascent mode** explicitly separates incomplete nascent RNA from mature isoforms | Targets nascent/pre-mRNA, not 5'-degradation of mature mRNA in dRNA |
| **ESPRESSO** ([Gao 2023, *Sci Adv* 9:eabq5072](https://pmc.ncbi.nlm.nih.gov/articles/PMC9858503/)) | Jointly considers all reads at a locus with per-read error profiles to fix splice sites | Corrects *junctions*; endpoint containment is not the target |
| **FLAIR / FLAIR2** ([docs](https://flair.readthedocs.io/en/latest/modules.html)) | `--check_splice` and short-read junction input raise splice confidence; users report UTR-localised "novel" models being likely false positives ([issue #104](https://github.com/BrooksLabUCSC/flair/issues/104)) | Explicitly recommends *orthogonal short reads*; no endpoint model |

The orthogonal-truth methods answer the question directly but need a different
library: **CAGE** and **RAMPAGE** ([Batut & Gingeras 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4372803/);
[ENCODE standards](https://www.encodeproject.org/data-standards/rampage/))
sequence 5'-complete cDNAs, and **ReCappable-seq**
([Yan 2022, *Genome Res* 32:162](https://genome.cshlp.org/content/32/1/162.long))
plus its nanopore form **ONT-cappable-seq** / Nanopore ReCappable
([Ugolini 2022, *NAR* 50:e59](https://doi.org/10.1093/nar/gkac144)) select
capped 5' ends. Standard dRNA carries **no cap signal**, which is precisely why
the question must be answered statistically here — and why §6 lists cap data
as the way to convert the hard rung from probabilistic to verifiable.

**The gap this study fills:** a per-candidate, three-way, FDR-controlled test
against an explicit *degradation* null, which reports `unidentifiable` instead
of guessing when dRNA physically cannot distinguish the two hypotheses.

## 2. The statistic: conditional termination hazard

Raw 5'-end density is the **wrong** object — it decays with depth simply
because fewer reads reach further 5'. The right object walks in the sequencing
direction (3' → 5') in spliced transcript coordinates:

```
offset d    = spliced distance from L's 5' end
at_risk(d)  = reads that actually reached d      (5' offset <= d)
ends_at(d)  = reads whose 5' end lies in bin d
hazard(d)   = ends_at(d) / at_risk(d)
```

Under degradation only, this hazard is smooth and slowly varying. A genuine
internal TSS at `d0` injects an **extra spike**, because reads transcribed from
S can only *begin* at `d0` — they are not the tail of a smooth decay.

Test: `H0` = the candidate bin terminates at the background hazard;
`H1` = an extra population starts exactly at `d0`. The p-value is an
**empirical bootstrap** under `H0`, not an asymptotic chi-square, because the
mixture weight is a boundary parameter (π ≥ 0) where the asymptotic null is
wrong.

Two design decisions were forced by the data, not chosen a priori:

1. **A global background is not safe.** Using one pooled hazard,
   degradation hotspots in single-isoform transcripts were called TSS at a
   **52%** rate. The null is now the **local neighbourhood median** hazard
   (candidate bin excluded), which asks whether the bin is anomalous *for the
   degradation regime it sits in*.
2. **No upstream reads is evidence, not a tie.** If nothing extends 5' of
   `d0`, degradation cannot be the explanation — degradation requires the
   parent to be present, and a smooth process would not deposit every read in
   one bin. Measured on SIRV, this pattern marks real nested transcripts. The
   first implementation called it `unidentifiable` and was wrong. It is
   recorded as `parent_unobserved`, because it supports **S's existence** and
   says nothing about whether **L** exists.
3. **The peak was CHOSEN, so a pointwise p-value is invalid.** EndpointRefine
   proposes a start because it already looks like a peak, then the test scores
   it on the very reads that made it look that way. The pointwise bootstrap
   p-value is therefore multiplied by the size of the search space that was
   scanned (**Bonferroni**): `p = min(1, p_pointwise * n_eligible_bins)`.
   A simulated max-statistic was tried first and rejected — it would need the
   risk set at *every* candidate bin, and reusing the candidate's `n_at_risk`
   for all of them understates the null's spread and is anti-conservative.
   Bonferroni needs only the search-space *size* and is conservative under any
   dependence between bins.

   The search space is **declared by the caller, not inferred from the data**.
   The default is 1, meaning "this offset was specified in advance". The
   pipeline passes the first exon's length in bins, because that is exactly
   the region EndpointRefine can propose a start in; deriving it from the
   observed read span would scale the penalty with sequencing depth, which is
   not a search space.

4. **Multiplicity is corrected at TWO levels, not one.** Within a parent,
   each alternative gets `p_within = min(1, p * k)` for the `k` alternatives
   tested on that parent, and the parent contributes the minimum of those
   upward, so one parent is one test. Across parents,
   **Benjamini-Yekutieli** rather than BH: alternative starts share reads,
   coverage and the degradation background, so positive dependence is not
   guaranteed and BY is valid under arbitrary dependence.

   A candidate stays `supported` only when **both** its own `p_within` and its
   parent's `q_value` clear alpha. Propagating only the parent q was a real
   bug: one strong alternative would certify its weak siblings. Both values
   are exposed per row.

   **Naming honesty:** in the pipeline the grouping key is the parent
   *candidate model*, not a gene. De novo assembly has no gene annotation, so
   this is **parent-model FDR**; only the SIRV study, which has a GTF, groups
   by real `gene_id`. Structural rules and abstentions carry `p_value = None`
   (JSON `null`, never `NaN`) and never join the family. **FDR applies to
   every calibrated HAZARD test** — which includes `own_tes` candidates that
   have a distinct TSS and are therefore hazard-tested — **and never to
   structural rules**: the 3'-end (TES-only) decision and the
   parent-unobserved dominant-start decision are threshold rules, not
   probabilities.

5. **Sharp peak, not broad hump.** `peak_mad` (median absolute deviation of
   in-bin 5' ends) is reported and a too-broad peak is rejected, so a wide
   degradation hotspot cannot pass as a tight TSS.

Implementation: `fin/analysis/tss_evidence.py`.

---

## 3. Ground-truth validation on SIRV

`experiments/prod_validation/tss_containment_study.py`, SIRV spike-ins with
absolute truth (84 SIRV transcripts, sample `SGNex_H9_directRNA_replicate3_run2`).

**Inventory:** 29 nested (short ⊂ long) pairs over 10 distinct contained
transcripts; 26 pairs have both members expressed; **28/29 differ at both
ends**, so SIRV populates rungs 1-2 and contains **no expressed `tss_only`
positive at all**.

**Negatives** — two kinds, and getting them right mattered more than the model:

* *pseudo-TSS*: random internal offsets with no real start nearby;
* *hard negatives*: the **maximum 5'-end pileup** position inside a
  single-isoform transcript — precisely the position most likely to fool a
  counting heuristic.

> **A label bug found by the data.** The first hard-negative run showed a 52%
> false-positive rate. Inspection revealed **13 of those 15 "false positives"
> were the method correctly finding the TSS of a *different real* SIRV
> transcript** overlapping the same locus (e.g. SIRV101 offset 350 = SIRV102's
> real start at 342). The labels were wrong, not the method. Negatives are now
> excluded from the neighbourhood of **any** real transcript start.

**Positives are counted once per short transcript.** A short model contained
in several long parents is one biological observation; the first run scored
each (short, parent) pair separately and inflated 4 positives into 7. The
deduplicated set keeps, per short transcript, the best-powered parent.

**Results are reported PER RUNG and PER LOCUS, never pooled.** The two rungs
are scored by different rules whose scores are not calibrated to a common
scale, so a pooled AUC would be meaningless; it is not computed at all.
Candidate rows from one locus are also correlated, so locus-grouped counts are
reported alongside row counts. Degradation background 0.0085 per 25 bp.

| Rung | unit | verdicts |
| --- | --- | --- |
| `own_tes` (positives) | 4 rows | 3 `supported`, 1 abstain, **0 wrong** |
| `own_tes` (positives) | **3 loci** | **3 `supported`** |
| `tss_only` (negatives) | 85 rows | 81 `unsupported`, 3 `supported`, 1 abstain |
| `tss_only` (negatives) | **8 loci** | 5 `unsupported`, 2 `supported`, 1 abstain |

After the full correction (Bonferroni selection -> within-parent -> BY across
loci), both views are emitted:

| view | positives | negatives |
| --- | --- | --- |
| **rows** | 3 `supported`, 1 abstain | **0 `supported`**, 84 `unsupported`, 1 abstain |
| **loci** | **3 `supported`** | **0 `supported`**, 7 `unsupported`, 1 abstain |

Aggregation rule, applied consistently and recorded in the JSON: a locus is
`supported` if ANY of its rows is supported, else `unsupported` if any row is
unsupported, else `unidentifiable`.

**Every SIRV false positive is eliminated, and all 3 positive loci remain
`supported`.** The fourth positive transcript **abstains** (`unidentifiable`,
zero reads at risk after 3'-end partitioning): it is retained, not dropped,
but it is unverified — abstention is not a verification. Row counts are
reported alongside locus counts because rows from one gene are correlated; the
locus is the GTF `gene_id`, attached explicitly to each row rather than
inferred from the candidate name.

**No model-vs-heuristic AUC is quoted, and the earlier "PR-AUC 1.000" is
withdrawn.** It was an artifact of mixing rungs plus non-tie-aware ranking.
The pooled AUC fields have been removed from the output rather than
annotated, because a caveat does not make an invalid metric valid.

## 4. The hard rung: ground-truth simulation

Because SIRV contains no expressed `tss_only` positive, the only way to
characterise the hardest case with known truth is to generate it. Reads are
simulated with the per-base degradation hazard derived from the measured
sample background; the short transcript shares the chain and the TES and
differs only by an internal TSS.

| depth | short = 10% | 25% | 50% | decided negatives |
| ---: | ---: | ---: | ---: | --- |
| 10 | *all abstain* | *all abstain* | *all abstain* | none decided |
| 20 | 0.097 | 0.741 | **1.000** | 31 TN, 0 FP |
| 50 | 0.625 | **1.000** | **1.000** | 40 TN, 0 FP |
| 100 | 0.650 | **1.000** | **1.000** | 40 TN, 0 FP |
| 200 | 0.550 | **1.000** | **1.000** | 40 TN, 0 FP |

**Totals: TP 368, FN 82, FP 0, TN 151**, plus 150 positive and 49 negative
abstentions. Recall among decided cases 0.818.

The FPR is **undefined at depth 10**, where every negative replicate abstains
and no negative is decided. The correct statement is "**zero false positives
among all 151 decided negatives**", not "FPR 0 at every depth".

Operating envelope, stated plainly:

* **< 20 reads at the locus → the method abstains.** It does not invent a
  transcript.
* **≥ 20 reads and the short isoform ≥ 25% of them → reliable detection.**
* **Minor isoforms at ~10% need ≥ 50 reads and are still only ~60% detected** —
  a genuine sensitivity limit, not a bug.
* Degradation alone was **never** called a TSS at any depth tested.

---

## 5. What is wired into the pipeline

`--tss-evidence-mode off|audit|require` (requires `--endpoint-refine`;
`off` by default, enabled by no profile):

* `off` — current behaviour, byte-identical output.
* `audit` — writes verdicts into `endpoint_refine.json` (`tss_verdicts`,
  `degradation_background_hazard`) and changes nothing else.
* `require` — only endpoint states whose alternative TSS is not `unsupported`
  are kept; if every alternative start of a model is refuted the model reverts
  to unsplit rather than emitting a lone primary state. **`unidentifiable`
  never drops a model**: insufficient evidence is not evidence of absence.

Inside the pipeline the same rung logic applies: a state sharing the primary's
3' end is tested as `tss_only`, while a state with its own 3' end first
restricts the read pool to reads ending there. Benjamini-Hochberg is applied
over every test in the run before `require` acts, so the pruning decision is
FDR-controlled rather than per-test.

The degradation background is estimated at run time from models that have no
contained sibling, so a locus already holding a shorter model cannot
contaminate its own null.

### Measured on a real human sample (H9 replicate2_run2, p00, precision profile)

Run-time background hazard **0.011874** per 25 bp, estimated from
nesting-free models. EndpointRefine proposed 150 splits / 300 endpoint states.

**Most real proposals are not TSS questions at all.** With production rung
routing, 103 of the 150 alternative states differ at the **3' end** (`own_tes`,
dRNA's reliable end), only **43** are the hard `tss_only` rung, and 4 are
unmappable. That distribution is itself a result: EndpointRefine is mostly
proposing 3'-end splits, where the evidence is strong.

| Verdict | n |
| --- | ---: |
| `supported` | **117** |
| `unidentifiable` | **33** (insufficient depth after 3'-end partitioning) |
| `unsupported` | **0** |

Only **14** of the 150 tests are selection-corrected hazard tests carrying a
p/q value; the other 136 are structural 3'-end or abstention calls with
`p_value: null`. The FDR statement therefore applies to those 14, and the
diagnostics say so per row rather than implying run-wide FDR control over
rules that have no calibrated probability.

Effect on output, verified byte-for-byte against the plain `--endpoint-refine`
run and against a repeat run:

| Property | Result |
| --- | --- |
| `audit` vs `--endpoint-refine` alone | **byte-identical** (`assembly.gtf`, `scores.tsv`) |
| `require` vs `audit` | **byte-identical** here (nothing was refuted) |
| transcripts | 10,362 in both |
| mass balance error | **1.46e-10** |
| repeat run | **byte-identical**, `tss_verdicts` identical |

Readings:

1. **`audit` is genuinely inert** and can be run on any sample at zero risk.
2. **Nothing was refuted on this sample**, so `require` changed nothing. An
   earlier build reported "2 unsupported"; that was an artifact of sending
   `own_tes` candidates through the 5' hazard test, and it disappeared once
   production used the same rung routing as the study. The corrected result is
   that EndpointRefine's existing guards already survive the statistical test.
3. **22% of proposals remain undecidable at this depth.** Those 33 states are
   labelled `unidentifiable` rather than silently asserted — the label marks
   exactly where deeper data or cap evidence would change the answer.

> **Provenance note.** The live `audit`/`require` artifacts, byte-identity,
> mass-conservation and determinism checks referenced above were produced
> before a final naming pass (`apply_locus_fdr` -> `apply_grouped_fdr`,
> `verdict_bh` -> `verdict_adjusted`) and the accompanying wording fixes.
> That pass changed **no runtime decision logic** — the SIRV study was re-run
> after it and reproduces the same verdicts — so the artifacts remain valid.

**Known limitation (verified by test):** `plan_endpoint_splits` re-cuts only
the **first exon**, so an alternative TSS in a downstream exon is out of scope —
that is a chain change, not an endpoint change.

---

## 6. What still needs your new data

* **A real `tss_only` positive set.** SIRV has none; the hard rung is currently
  characterised only by simulation.
* **Orthogonal TSS truth.** CAGE / RAMPAGE / ReCappable-seq (NRCeq) cap
  signal is the gold standard and would convert rung 3 from "probabilistic"
  to "verifiable". Standard dRNA has no cap signal.
* **Cross-sample recurrence.** A real alternative TSS recurs at the same
  coordinate across samples; degradation is stochastic. This is a cohort-only
  feature and is deliberately absent in single-sample runs.
* **Promotion decision.** `require` must be shown to reduce false contained
  models on real holdouts without losing genuine short isoforms before any
  profile enables it.

---

## 7. Reproduce

```bash
# SIRV ground-truth study + hardest-case simulation
python3 experiments/prod_validation/tss_containment_study.py \
  --gtf experiments/prod_validation/sirv4/_ref/full/annotation.gtf \
  --bam <sirv>/input.bam --nanocount <sirv>/nanocount.tsv \
  --out <outdir>

# audit verdicts on a real sample (output unchanged)
python3 -m fin.cli --profile real-drna-precision ... \
  --endpoint-refine --tss-evidence-mode audit
```
