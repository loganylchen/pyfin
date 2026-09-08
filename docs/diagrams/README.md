# The PyFIN pipeline, as it runs today

Five multi-panel figures. Everything below describes the **production default**
(`--profile real-drna`) at commit `aa21414`. Options that ship default-off are
named as options, never as part of the flow.

Two figure sets share this specification. The **illustrated set** (BioRender
style, icon-driven, for talks and the README) and the **schematic set** (Nature
Methods style, denser labels, for the methods section). Each prompt in
`prompts/` was reviewed by Codex against the code before generation — the first
pass returned "revise" on all four (inputs shown incomplete, "one pass"
overclaimed, the junction channel drawn as an additive cost, conditional outputs
drawn as unconditional, a result badge that read as a guarantee), the second
approved them. Known renderer liberties: `figB` added two accurate sub-labels
("transcript (annotated end)", "reads (with polyA)") not in the prompt, and draws
the read-supported extension as an appended block rather than a stretched last
exon.

| illustrated | covers |
| --- | --- |
| `figA_overview.png` | direct RNA → inputs → PyFIN (per interval / whole sample) → outputs; the polyA evidence ribbon |
| `figB_discover.png` | reads → chain grouping → family → collapse (shadow) → candidates; the extent sibling |
| `figC_assign.png` | one read: best-alignment tie, three evidence channels (junction = mask/prior, endpoints = cost), soft assignment, abundance; the endpoint kernel |
| `figD_finalize.png` | aggregate → select → snap → merge → refit; the snap-collision merge; the tested-grid badge |

| schematic | panels | what it covers |
| --- | --- | --- |
| `fig1_workflow.png` | a, b, c | the whole pipeline: sample prepass, per interval, finalize |
| `fig2_discovery.png` | a, b | chain clustering; the extent sibling |
| `fig3_evidence.png` | a, b, c | the endpoint kernel; the tie set; what each channel does |
| `fig4_merge_rule.png` | a, b | the snap-collision merge and its admission rule |
| `fig5_result.png` | a, b | why the merge cannot cost recall; the four-grid result |

Codex fact-checked this file against the code twice. The first pass found twelve
inaccuracies, the second found that nine were fixed, one over-corrected, two
under-corrected, plus five further imprecisions; all of those are folded in
below, and the two figures whose STRUCTURE was wrong were regenerated. Where a
statement is an empirical measurement rather than something the code enforces, it
is labelled as such.

**The default operating point** (verified against `fin/cli.py:176-248` and
`fin/pipeline/config.py:38-81`): `chain_cluster_discovery=True`,
`clustering="families"`, `quant_mode="m2_em"`, `canonical_gate=True`,
`junction_snap=True`, `snap_merge_chain=True`, `assignment_channels="ends"`,
`folded_tss_recovery_mode="require"`, `post_selection_refit=True` (and
*effective*, not merely requested), `abundance_feedback=False`,
`require_gtf_positive_mass=False`, `min_polya5p_reads=0`, `min_abundance=1.0`,
`min_gtf_abundance=1.0`, `threads=1`.

---

## Figure 1 — `fig1_workflow.png`: the three layers

*Panels: **a** sample prepass, **b** per interval, **c** finalize. Each box carries
its stage name and the two or three facts that matter for that stage.*

**Sample prepass (once).** `runner.run()` calls `generate_isolated_intervals()`.
Intervals are isolated by read coverage **and by annotation spans**, so an
annotated transcript can join two coverage islands into one interval. Isolation
is per chromosome and strand, not a claim of globally disjoint genomic
coordinates.

Then, when `candidate_existence_mode != "off"` **or**
`assignment_channels != "none"` — the second is what triggers it in this profile
— `_prepare_sample_evidence()` collects, per read, the krill polyA call and the
strand-aware genomic 3'/5' ends of its **primary** alignment, deduplicated by
read id. It is at most one krill pass: a usable cached artifact skips it. (The runner
tolerates a run with no signal — the 3' channel then abstains — but the
production CLI does not offer that path: a missing `--signal` is a fatal
missing-option error. Nor does the 5' term automatically survive it; the start
model needs enough annotation-derived observations to fit.) The result is written
to an identity-addressed, create-only artifact so every worker judges intervals
against the same evidence.

From it, `fin/analysis/sample_nuisance.py` fits the endpoint nuisance model: a
two-Laplace mixture (narrow jitter scale, broad error scale, broad weight `eta`)
over 3' distances, plus a one-sided 5' start model. The fit uses the input
annotation's claimed TES positions when at least two such distances survive
exclusion, otherwise split-half distances to the sample's own read modes. For
**assignment** the fit is leave-one-interval-out. The **discovery**
contradiction hook deliberately uses a pooled, non-excluded fit.

**Per interval.** Serial at the default `threads=1`; process parallelism is
opt-in. Discovery → gates → align → signal + ends → soft assign → interval
select. Figures 2 and 3.

**Finalize (once).** The band's AGGREGATE box is `aggregate_across_intervals()`,
and it — together with an `unify_extent_pairs()` step the figure does not draw —
runs in `run()`, *before* `_finalize_and_write()`. Everything from the SELECT box
on (`select_global` → `apply_junction_snap` → `refit_survivor_abundance` →
write) is inside `_finalize_and_write()`.

---

## Figure 2 — `fig2_discovery.png`: candidate discovery in one interval

*Panels: **a** chain clustering, **b** the extent sibling.*

`discover_candidates()` → `_chain_cluster_candidates()`.

1. **Group by intron chain**, ignoring endpoints. A read's chain comes from its
   CIGAR; annotated chains enter as hypotheses too, including annotation-only
   keys with no read.
2. **Cluster into families** by junction wobble (`chain_cluster_wobble_bp=6`),
   small cassette exons (`chain_cluster_cassette_max_exon_bp=70`) and
   containment.
3. **Then collapse within each family** — families are built *first*, folding
   happens *inside* them (this is the reverse of what the earlier draft said).
   An exact sub-chain folds into its longest *eligible retained* container and is
   kept as a dormant **shadow**, not discarded. Two exceptions matter: an **annotated chain never
   becomes a foldee**, and under the default span guard a read whose CIGAR
   exonically traverses an extra container intron stays on the shorter chain.
4. **Emit** per surviving member: the annotation model(s) owning that exact chain
   if the GTF has one (annotated geometry, pooled member reads as support), else
   a novel candidate spanning `min(read starts) … max(read ends)`. Single-exon
   members take a **separate route** — per-locus spatial grouping and overlap
   matching against single-exon annotations, not empty-chain equality.
5. **Extent siblings.** When an annotated model's reads CONTRADICT its claimed
   ends under the fitted kernels, a `<gtf_id>.ext` sibling is emitted so
   annotation and reads compete as two hypotheses. **Each contradicted endpoint** moves, and only
   **outward**; each **uncontradicted** endpoint keeps its annotated coordinate.
   Both ends can be contradicted in the same sibling. Contradiction is a narrow-component posterior test, reads that an
   *alternative* annotation explains are excluded, and the 5' side additionally
   requires a start CLUSTER. The sibling carries `source="gtf"`.

In panel **a**, the greyed shadow chain belongs to its surviving container —
it is provenance, not an independently emitted candidate. The three output arrows
are **conditional alternatives**, not three things emitted for every member: a
member becomes the annotated model *or* a novel candidate, and the extent sibling
appears only when a contradiction is found.

**Then the interval gates:** the junction registry is exported (no gate reads
that file, which is not to say junction evidence is unused elsewhere); the
canonical-motif gate drops novel candidates with a non-canonical junction
(GTF/fusion/mono exempt, abstains without genome sequence); the
junction-dominance gate exists but is default-off.

---

## Figures 3 — `fig3_evidence.png`: which candidate does a read belong to

*Panels: **a** the fitted endpoint kernel (narrow jitter spike on a broad error
hump, weight `eta`), **b** the per-read tie set with the three channels feeding
only the tied cells, **c** what each channel actually does.*

Per interval, inside `_quant_m2_em()` → `Assigner.assign()`.

1. **Align.** Every kept read is aligned with `mappy` against every candidate's
   spliced cDNA. A read's best-scoring candidates form **its own** tie set — a
   per-read set, not a globally tied column block.
2. **Add evidence the alignment does not have — to the tied cells only**, and
   **before** the assignment engine runs, not after it:
   - **Junction signal (M2).** Where tied candidates differ at a junction, the
     raw signal in the differing window is scored under krill by summed
     log-likelihood ratio. It does **not** simply add a cost: a sufficiently
     separated contrast becomes a **mask** over the losing structures, while an
     ambiguous read instead gets a covered-read vote prior or a flat fallback.
     It abstains when the tied candidates have unequal intron counts or no
     differing window.
   - **PolyA 3'.** For a read whose polyA call admits it, the strand-aware
     distance from its genomic 3' end to each candidate's claimed TES is scored
     under the fitted endpoint mixture and added as a cost. This is the
     direct-RNA-specific channel.
   - **Upstream 5'.** A read starting upstream of a candidate's claimed TSS is
     scored under the fitted one-sided start mixture — narrow jitter *and* broad
     error components, as a penalty relative to the likelihood at zero offset —
     and added as a cost. *(Correction: this
     assignment term applies per read; the read-CLUSTER requirement belongs to
     the discovery contradiction hook of Figure 2, not here.)*
   Both endpoint terms require at least two *eligible* candidate cells in the
   tie — cells not masked out, tested as `d_tx < 1e5` rather than by mathematical
   finiteness.
3. **Soft assign.** The engine is run at `beta=0` with
   `abundance_feedback=False`, so it is a **fixed-cost softmax**, not an
   iterative θ update — the θ branch exists but is disabled by default.
   Reported abundance is then the column sum of responsibilities — accumulated
   across all reads, not derived from the single read the figure follows. The
   boxes are labelled "SOFTMAX → RESPONSIBILITY → ABUNDANCE" for that reason.
4. **Still inside this stage, after the first assignment:** folded-TSS shadow
   recovery tests each dormant shadow's proposed TSS against the parent's
   read-5'-end hazard profile (with parent-mass, pure-5'-truncation and
   TES-agreement preconditions, so not every shadow is tested) and re-emits the
   supported ones, then reruns assignment **once** — and that second pass carries
   a recovery-specific positional likelihood, distinct from the ordinary upstream
   channel, which penalises the parent at the supported start bin or masks the
   recovered candidate for upstream starts. A recovered identity that already
   exists only gains support, without necessarily triggering a rerun. Post-EM
   mono resolution then runs through `select_m2_interval()`. **Interval selection also drops
   candidates here**, after assignment — global selection is not the only filter.

---

## Figures 4 and 5 — `fig4_merge_rule.png`, `fig5_result.png`: finalize, and the duplicate-model fix

*Figure 4 panels: **a** the discovery/snap/merge sequence, **b** the six admission
conditions and the three refusals. Figure 5 panels: **a** why removing a
duplicate cannot cost recall, **b** the four-grid result.*

Cross-interval:

1. `unify_extent_pairs` then `aggregate_across_intervals` (both in `run()`).
2. `select_global` — abundance floors, isoform-fraction and mono-exon rules.
   Both floors are numerically 1.0 but they are not the same rule: the profile
   sets `strict_novel_abundance_floor`, so a **novel** candidate needs abundance
   **> 1.0** while a **GTF** candidate passes at **≥ 1.0**. *(Correction:
   the strict-positive GTF gate `require_gtf_positive_mass` and the polyA+5'
   retention filter `min_polya5p_reads` are **not** active in this profile, and
   `enforce_gtf_positive_after_refit` is called with `enabled=False` here rather
   than skipped.)*
3. `apply_junction_snap` — novel multi-exon junctions within 6 bp of a more
   strongly supported read consensus are corrected (subject to a minimum support
   and a support-ratio threshold, with canonical target checks); then models with
   identical exon tuples are merged; then, new, models the snap made into one
   chain are merged.
4. `refit_survivor_abundance` — responsibilities are re-routed through the
   redirect map and renormalised over surviving weights, and abundance,
   confidence, `max_R` and the hard assignments are recomputed. Not every read
   takes that path: mono-resolved reads use forced target assignments, and reads
   with no viable survivor are explicitly accounted as selection-orphaned rather
   than reassigned somewhere arbitrary.
5. Write GTF + TSV.

**What Figure 4a shows.** *Discovery* emits two **different** chains that sit in
one family, differing at a junction (the circled gap) — not two variants of one
chain, which is precisely the case this merge refuses. *Snap* corrects the minority junction onto the
consensus, so the two now carry the same chain and differ only in where the first
exon begins. Note that the differing start is already present in the DISCOVERY
stage and the red circles appear there only: snapping moves **internal donor/acceptor coordinates only**, never the
outer endpoints. (The figure reads as a plus-strand example; the code is
strand-aware and takes the 3' end as `start` on the minus strand.) *Merge*
reconciles them into one model.

The old merge key was the FULL EXON TUPLE, so that last step did not happen: any
terminal difference kept the two apart. One human cell logged "snapped 249
junctions and merged 1 duplicate models" while carrying 141 redundant reference
matches *(measured, not enforced by code)*.

`snap_merge_chain` (now default) adds the second pass. It is deliberately narrow,
because one corrected intron chain is **not** one transcript — APA and
alternative-TSS isoforms share a chain. A merge requires all of: every member
novel and multi-exon, one chromosome, a known strand; the same post-snap chain;
the collision CREATED by the snap (per-survivor pre-snap chain sets non-empty and
pairwise disjoint, at least one corrected); ONE non-null discovery family across
every original member; identical strand-aware 3' ends; and a 5' spread within
25 bp as a group **diameter**, not transitive proximity. The survivor keeps its
own geometry.

**Measured** against the previous default (expression contract, paired bootstrap,
seed 20260831): SIRV4 dPr +1.45 / dF1 +0.70 over **156 cells**, HEYA8 +0.82 /
+0.32 over **130**, human GENCODE +0.30 / +0.043 over **7** (one sample, seven
conditions — not a full human grid), Sequin inert over **104**. Reported recall
is unchanged to two decimals on every condition of every grid and no evaluated
cell degraded. These are empirical results on those grids, **not guarantees the
merge algorithm enforces** — the merge criteria never consult reference matches.

---

## What the figures deliberately leave out

Three omissions a reader should know about, in order of importance:

1. **Post-assignment interval selection.** Structural, support and containment
   gates inside the interval remove candidates *after* assignment; the figures
   show only one "INTERVAL SELECT" box for all of it.
2. **The junction term's real mechanism.** Figure 3 draws three arrows into the
   tie as if all three were additive costs; the junction term can instead mask
   cells or supply a vote prior.
3. **Post-assignment identity changes.** Supported shadows add candidates and
   trigger a second assignment; annotation/extent alternatives resolve to one
   model; the refit reconciles final read mass.

## Wording that is operational, not biological

"A polyA call states this end is a polyadenylation site", "a read cannot begin
upstream of the transcript it came from", "supported shadow", "corrected
junction", "duplicate model" — each of these is a rule the code implements
(admitting a caller verdict, applying a penalty, passing a hazard test, meeting a
support ratio, satisfying the merge criteria). None is a demonstrated biological
fact about any individual molecule.

---

*Figures generated with `gpt-image-2` through the same OpenAI-compatible endpoint
the Codex CLI is configured against. They are schematic: this text is the
specification, the images are the map.*
