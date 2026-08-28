# PyFIN Algorithm Review

Scope: original review at `dev` `a2c4be4e`, followed by the verified profile
integration worktree with source SHA-256
`42e1210dea81860686e76aaf9fa058fc8a65808e0a31367570d8526f4da17954`.
Evidence includes formal tests, fresh SIRV tuning/holdout, full real p00 sweeps,
production manifests, and Advisor guidance.

## Resolution status

| Finding | Status | Verified remediation |
| --- | --- | --- |
| High 1: hidden SIRV defaults | Resolved with independent real samples | Balanced/precision real profiles plus SIRV; three real samples validate the new existence levers |
| High 2: secondary structural support | Resolved | Secondary records skipped before read ownership/generation; BAM regression test |
| High 3: one-score M2 certainty | Resolved | Two-score abstention plus per-read invalid-sequence abstention; one bad sequence no longer poisons a locus batch |
| High 4: summed LLR not wired | Resolved | `off|mean|summed_llr|auto`; real uses tight sum; SIRV uses SUM unguided and mean guided |
| Medium 1: family ownership lost before selection | Resolved | Stable discovery family IDs propagate through assignment/aggregation; family-aware fraction default passed two real paired tests |
| Medium 2: random candidate IDs | Resolved | Stable structural BLAKE2b IDs remove run-random output and equal-abundance tie breaks |
| Medium 4: M3 missingness semantics | Resolved by removal | Production M3 wiring removed; prototype preserved under `experiments/m3_coherence/` |
| Medium 7: no M2 route control | Partially resolved | Explicit metric/margins/auto route added; one-softmax semantics and inert options remain |
| Medium 8: config validation dead | Resolved | CLI calls `cfg.validate()` and reports Click usage errors |
| Medium 10/13: reproducibility/integration drift | Substantially resolved | Atomic source-hashed manifest and `PROFILE_OPTIMIZATION.md` ledger |
| High 5: endpoint isoforms | Open, measured lower impact | Only 62/4,477 (1.38%) missed T3 multi-exon truths were same-chain endpoint collapses with an emitted chain |
| High 6: post-drop abundance refit | Resolved for production beta=0 M2 | Named profiles renormalize sparse read responsibilities on final survivors and report explicit orphaned mass |

The findings below preserve the original evidence and rationale. Resolution
paragraphs state where the live source now supersedes the reviewed snapshot.

## Findings

### High 1 [Resolved]: Raw CLI defaults were unsuitable for the stated real-dRNA domain

Resolution: balanced `--profile real-drna` is now the CLI default;
`real-drna-precision`, `sirv`, and `custom` are explicit alternatives. Balanced
real disables full/polyA and soft-mass gates, uses a >1 read-equivalent floor,
requires five hard reads for novel mono models, applies read-supported junction
consensus, and uses tight summed LLR. The mono/consensus gains reproduced on two
additional H9 samples; the precision profile separately exposes the generation
support=2 trade-off instead of hiding its roughly 1 pp recall cost.

Evidence:

- `fin/cli.py` defaults `min_fulllen_fraction=0.1`, `min_polya5p_reads=1`,
  `min_abundance=3`, and `min_gtf_abundance=1`.
- `fin/pipeline/selection.py:359-361` labels abundance floors SIRV-tuned and
  states that they gut genuine low-abundance recall on real data;
  `selection.py:479-490` gives the corresponding warning for full-length gating.
- `fin/analysis/quantification.py:342-371` labels the 0.1 full-length threshold
  SIRV-tuned.
- `fin/analysis/quantification.py:426-462` gives the same warning for polyA/5'.
- `experiments/prod_validation/PRODUCTION_STATE.md` explicitly says the two
  SIRV finalize gates can halve real-dRNA recall and documents a different real
  profile with both disabled.

Why it matters:

The package describes itself as a direct-RNA transcriptome assembler, but a raw
invocation selects the synthetic benchmark profile. A user can obtain a
biologically different operating point without realizing it. Aggregate SIRV
F1 can improve while real low-expression, 5-prime-truncated, or polyA-unscored
isoforms disappear.

Recommendation:

- Make a first-class `real-drna` profile the default for the product domain.
- Require explicit `--profile sirv` for the synthetic gate set.
- Print the resolved profile and every active selection threshold in the run
  manifest/log.
- Keep tuning and final evaluation samples disjoint.

### High 2 [Resolved]: Candidate generation admitted secondary alignments

Resolution: discovery now skips secondary records before `all_read_ids`,
sequence capture, and structural generation. A primary+secondary same-query BAM
test proves only the primary intron chain contributes.

Evidence:

- Interval construction skips secondary records (`interval_manager.py:233-237`).
- Discovery refetches the interval and captures sequence only for primary records
  (`discovery.py:474-483`) but does not skip a secondary record before appending
  it to `non_fusion_reads` (`discovery.py:484-486`).
- Every appended alignment contributes a CIGAR intron chain
  (`discovery.py:250-260`).
- `cluster_families` and `collapse` document a precondition that per-variant
  read sets partition the family pool (`chain_cluster.py:577-581`). A read with
  two secondary/primary chains can appear in multiple variant sets.
- No discovery test exercises a secondary alignment.
- The production mapping script uses `--secondary=no`
  (`benchmarks/map_drna.sh:190`), which reduces exposure for that manifest but
  does not enforce the input contract for arbitrary BAMs.

Impact:

A multi-mapping read can support multiple structural variants, inflate discovery
support, create false family bridges, or satisfy ratios that assume independent
read ownership. The assignment matrix later has one row per read ID, so
structure generation and quantification operate on inconsistent populations.

Recommendation:

Skip secondary alignments before `all_read_ids`, `non_fusion_reads`, chain
extraction, and span collection. Add an explicit test where one query has primary
and secondary chains and prove it contributes to exactly one candidate family.

### High 3 [Resolved]: M2 converted partial technical scoring into infinite biological evidence

Resolution: both `tie_nll` and `m2_resolve_tie` require at least two successful
hypothesis scores. One/zero successful scores return no winner and preserve the
full M1 fallback tie. Candidate eventalign slices containing ambiguous bases now
abstain per read before batching; on H9 replicate3 this isolated 587 tied reads
with zero batch failures instead of degrading the final 8,486-read locus.

Evidence:

- `tie_nll` stores only candidate cells with a successful eventalign and at
  least one event (`assignment.py:197-203`, `227-233`).
- The coverage-gate path sorts the available scores. If only one candidate
  scored, `margin=float("inf")` (`assignment.py:403-410`).
- Any `margin >= 0.5` hard-selects that candidate and leaves all alternatives at
  `MISSING=1e6`, even when coverage was not established (`assignment.py:415-421`).
- A test in the older tie resolver explicitly calls one scored hypothesis an
  infinite-margin winner (`test_m2_resolve_tie.py:189-200`). The current
  gate-on tests do not cover partial scoring. This confirms that single-score
  infinity is a deliberate cross-generation design choice, not an accidental
  branch, while leaving the statistical objection unchanged.

Why it is unsound:

"Hypothesis B could not be scored" is not the same event as "the observed signal
strongly rejects B." Failure can come from slicing, event projection, backend
status, missing events, or model coverage. Treating missingness as certainty
creates non-random selection bias.

Recommendation:

A contrast must have at least two successfully scored hypotheses over an
explicitly comparable observation set. Otherwise abstain and retain the M1 tie
(or use a separately calibrated missing-hypothesis penalty backed by data).

### High 7 [Resolved]: One invalid M2 sequence poisoned a whole locus batch

Evidence from the independent H9 replicate3 run exposed one candidate slice
containing `N`. `krill.align_reads_variants` raised once for the full batch,
forcing every tied read in the final high-coverage locus into the slow singular
fallback. That fallback also lacks event positions, so `cover_by_read` became
empty and the locus silently changed from covered-vote priors to flat ambiguity.
Two pre-fix jobs exceeded 30 minutes without completing; the affected interval
contained 8,486 reads, 82 candidates, and 1,908 M1 ties.

Resolution: `tie_nll` validates every candidate payload as A/C/G/T/U before the
batch call. An invalid tie retains its full M1 fallback and contributes no M2
negative evidence; all valid reads remain in the batch. The fixed paired run
reported 587 abstaining reads, zero batch failures, completed in about 1,394
seconds, and has a regression test asserting that invalid payloads never enter
krill while valid tied reads still refine.

### High 4 [Resolved for production scope]: Default M2 omitted the implemented summed LLR

Resolution: summed LLR is wired into `m2_em`, limited to same-intron-count
contrasts, and independently parameterized. Real selected margin 1/flank 8. A
six-sample paired SIRV test found SUM +0.722 F1 when truly unguided but -0.284
with a full guide, so the SIRV profile resolves `auto` to SUM (1,4) for absent or
one-transcript stub GTFs and mean for usable guides. Both metrics score a common
genomic disagreement region while allowing each hypothesis its own eventalign.
SUM remains a tie-refinement metric, not an existence gate.

Evidence:

- Every tied candidate is eventaligned against its own sequence slice.
- The default genomic `gset` is shared, but `_mean_nll_in_gset` selects events
  after each candidate-specific event projection
  (`m2_junction_nll.py:1141-1201`). Event identities and counts can differ.
- `assignment.py:227-231` calls `read_cand_mean_nll` with no `reduce` argument,
  so default `m2_em` uses the legacy per-event mean over the wide class `gset`.
- The source itself explains why this is weak: mean averages a few
  discriminating events with many shared events (`m2_junction_nll.py:1357-1359`).
- Tight-window `summed_llr` is already implemented and unit-tested. It builds
  `diff_junction_windows`, passes `reduce="sum"`, and produces an undivided
  margin (`m2_junction_nll.py:1361-1365`, `1409-1452`). Its implementation
  comment states exactly that division by `n` dilutes the difference
  (`m2_junction_nll.py:1180-1182`).
- The only HEAD caller is `runner.py:699` inside `quant_mode="cluster"`. The
  winning/default `m2_em` path never calls it. A scan of all 33 historical
  `_tie_nll` versions found zero `reduce="sum"` or `summed_llr` uses.
- `PRODUCTION_STATE.md:57-59` says sum is active in the runner, while its own
  line 91 calls the sum/recovery work uncommitted. This is an integration claim,
  not evidence that default wiring exists.
- `--no-m2-tiebreak` also does not disable default `m2_em`: the flag and
  `m2_tiebreak_margin` control `_quant_argmax_keep`, while `m2_em` invokes
  `tie_nll` unconditionally.
- Corrected HEYA8 existence experiments over 6,737 candidates show structural
  `diff_decisive` AUROC 0.806 and read-support AUROC 0.774, while fixed broader
  private-fit signal metrics are roughly 0.49-0.51. That task is not the narrow
  summed wobble tie task.

Interpretation:

The user's recalled design is confirmed: use sum over tight differing events
because a wide mean dilutes the signal. It was implemented in `7420fe6a`, but
only for cluster mode; later `bc39ce8b` extracted the old mean path unchanged
into `assignment.py`. This was never a reversion because default `_tie_nll`
never used sum.

Event identity is intentionally not paired one-to-one. Alternative transcript
sequences can map the same raw read segment to different numbers and positions
of events; forcing event correspondence would discard part of the hypothesis
difference. The comparable object is the shared genomic disagreement window.
The remaining statistical risk is count-driven SUM scale when candidate-specific
event populations differ sharply. Current production guards this by restricting
SUM to same-intron-count contrasts, requiring two successful scores, using tight
windows and a separately validated margin; containment/cassette contrasts
abstain. Further work should calibrate candidate-specific event-set likelihoods
within shared genomic windows, not force equal event identities.

### High 5: The default generation model cannot represent same-chain endpoint isoforms

Evidence:

- `cluster_families` groups reads by exact intron chain and explicitly ignores
  3-prime position (`chain_cluster.py:438-456`).
- One candidate span is the union of all member-read extents
  (`discovery.py:343-352`).
- `three_prime_threshold` and canonical alternative expansion are bypassed by
  production family generation.
- Exact subchains are folded into longer chains (`chain_cluster.py:533-631`).
- Folded TSS recovery is implemented only inside `quant_mode="cluster"`
  (`runner.py:731-850`); default `m2_em` does not consume `CandidateSet.shadows`.
- Source comments record an earlier ~20% distinct-truth recall cost, 92% attributed
  to exact-subchain collapse (`config.py:29-35`). The span guard protects
  retained-intron evidence but not a genuine shorter terminal isoform with the
  same retained introns.

Impact:

Novel alternative TSS, TES, APA, and terminal-only isoforms sharing a splice
chain collapse into one union-span model. A genuine shorter 5-prime isoform can
also become dormant provenance with no default recovery path. This is a model
boundary, not a threshold bug. However, a direct real-p00 audit bounds its
current aggregate impact: among 4,477 missed T3 multi-exon truths, 95 belonged
to a multi-transcript same-chain group, and only 62 (1.38%) had that chain
already emitted by PyFIN. Endpoint representation therefore stays important for
biological completeness but is not the dominant present F1 gap; low-support
existence selection and junction-coordinate noise are larger.

Recommendation:

Use the existing `EndpointRefine` hook to define endpoint states within each
splice-chain family. Cluster strand-aware assigned-read TSS and TES modes, use
poly(A) evidence to corroborate TES/APA, and require a peak-over-degradation
model for TSS because dRNA 5-prime ends are frequently censored. Emit only
read-supported TSS/TES pairs rather than a Cartesian product, cap splitting,
and rerun abundance estimation after separation.

### High 6 [Resolved for production]: Candidate drops destroyed read mass

Resolution: named profiles now carry each interval's final beta=0 M2
responsibilities as a sparse read ledger while leaving the existence cascade on
its original evidence. Refit-enabled interval outputs are canonically ordered
before aggregation so worker completion order cannot move floating threshold
boundaries. After global selection and junction-snap merges, candidate
folds redirect complete soft columns, mono-resolved reads retain an explicit
mass-1 parent, and every other read is renormalized over final survivors. A read
with no surviving compatible column contributes one unit of selection-orphaned
mass instead of disappearing. The pass recomputes abundance, hard IDs/counts,
confidence, and `max_R`; `abundance_refit.json` separates alignment-unassigned
from selection-orphaned reads and asserts per-read and total conservation.

Paired live acceptance kept all structural record multisets and metrics
unchanged; equal-coordinate GTF groups may differ only in line order because
refit-on canonicalizes interval insertion while legacy off preserves worker
completion order. Independent release-source refit-on repeats produced
byte-identical GTF, TSV, and diagnostics artifacts. SIRV p00
recovered 54.56 read-mass with 78/1,310 assignable reads orphaned. Real tuning
recovered 25,622.35 read-mass; 49,795/294,111 assignable reads (16.93%) were
honestly unassigned because every compatible model was filtered. Independent
replicate3 recovered 52,651.97 mass with 109,027/594,431 reads (18.34%)
selection-orphaned. Refit-off
GTF/TSV remained byte-identical to v8/v6. Runtime and peak RSS were unchanged
within scheduler noise. Abundance-feedback and non-M2 modes warn-disable the
profile-implied pass; an explicit incompatible feedback/refit request errors.

Historical evidence:

- M2 cluster-recheck, junction-support, M2-support, containment-cluster, and
  wobble drops add columns to `drop_cols`; `selection.py:312-321` filters those
  columns without rerunning assignment.
- Global abundance, isoform, soft/hard, full-length, and polyA filters similarly
  remove results after assignment (`selection.py:370-523`).
- Only default post-EM mono resolution transfers hard reads/approximately one
  unit per read. Optional containment-collapse transfers soft abundance and hard
  IDs (`selection.py:241-281`).
- TPM is recomputed on surviving abundance, so it renormalizes the remaining
  mass rather than refitting reads to the survivor set.

Impact in the reviewed snapshot:

For structure-only assembly, filtering a model and discarding its assignments is
standard and can be acceptable. The High severity is scoped to the simultaneous
claim that the output provides abundance and TPM. For those fields, the result
is not a coherent estimator: reads assigned to a dropped shadow should often
support its nearest compatible survivor. Current output can undercount
survivors, lose total abundance, and then hide the loss by TPM normalization.
Sequential filters also use masses estimated in the presence of candidates that
later disappear. An alternative to requantification is to label these values
explicitly as model-filtered counts rather than transcript abundance.

Implemented strategy:

Separate existence selection from final quantification. After fixing the
survivor set, rebuild the compatibility matrix and rerun abundance estimation,
or transfer mass only through an explicit, tested fold mapping. Add invariants
for read-mass conservation and compare pre/post-selection assignments.

### Medium 1 [Resolved]: Discovery families were not propagated to selection

`cluster_families` already was the intended splice-family primitive, but family
identity ended at discovery and `isoform_fraction_drops` re-derived relatedness
as arbitrary same-strand genomic overlap. Stable BLAKE2b family IDs now cover
chromosome, strand, and sorted read/GTF chains; they propagate through
`TranscriptCandidate`, every quantification mode, cross-interval aggregation,
and junction merges. Attached GTF hypotheses share the family ID. Mono, fusion,
and legacy candidates remain family-less and use an overlap fallback.

`isoform_fraction_drops` now defaults to the family denominator;
`--isoform-fraction-locus overlap` reproduces the historical behavior. Paired
live validation passed the pre-registered 0.10-F1 guard: tuning T1 F1 +0.0076 and
T3 unchanged; independent replicate3 T1 +0.0218 and T3 +0.0059, with T3 recall
+0.0248 and precision -0.0056. SIRV p00/full were byte-identical between modes.
The earlier survivor-reconstruction audit (87/95 cross-family dominants) differed
from persistent discovery families because dropped bridge variants can preserve
a single-linkage family; this distinction is now explicit rather than hidden.

### Medium 2 [Resolved]: Random UUID candidate IDs broke reproducibility and affected ties

`_generate_novel_id()` used `uuid4()`, while `mono_resolve_drops` used
`candidate_id` to break equal-abundance host ties. Two identical precision runs
then produced one different transcript structure and a 0.028 T3-F1 difference,
confirming that the risk was observable rather than theoretical.

Novel IDs are now stable BLAKE2b hashes of chromosome, strand, endpoints, and
intron chain. GTF IDs remain untouched, and cluster-mode recovered candidates
use the same structural identity rule. A unit test fixes equality/change
invariants; duplicate final precision runs are required to be structurally and
metrically identical.

### Medium 3: Fusion-excluded records remain eligible for ordinary abundance assignment

Discovery adds read IDs/sequences before `is_fusion_read(rd)` skips the record
from generation (`discovery.py:474-486`). `CandidateSet.read_ids` then feeds every
quantification mode. When fusion is disabled, a long-soft-clipped/chimeric read
can still be assigned to an ordinary transcript; when enabled, it can compete
with both ordinary and fusion candidates.

This may be intentional for abundance salvage, but it is not documented as a
population decision. Define separate `ordinary_read_ids` and `fusion_read_ids`,
then explicitly choose whether a fusion read may contribute to an ordinary
isoform and under what alignment-coverage condition.

### Medium 4 [Resolved]: M3 missing comparisons encoded perfect similarity

The retired prototype filled cross-class and unscorable read pairs with distance
0, indistinguishable from perfect signal similarity in the coherence average.
Current named profiles never enabled it. Production configuration, CLI,
assignment, and ablation wiring have now been removed; the implementation and
focused tests are preserved under `experiments/m3_coherence/` for any future
masked/leave-one-out redesign. The generic EM engine remains, but production
M1/M2 calls it with `beta=0`. Real balanced/precision and SIRV p00/full outputs
were byte-identical before and after removal.

### Medium 5: Default structural gates are mostly tuned on synthetic truth

The 20 bp wobble/cassette/GTF guard improved synthetic HEYA8/SIRV honest F1 from
79.9 to 81.3, with recall restored by exact-junction support. That is good
mechanistic evidence. The report itself warns that dense engineered SIRV sites
are not representative and requires real-transcriptome confirmation.

The full novel-junction gate test shows near-zero mean F1 change: SIRV recall
-0.7 and precision +0.6, HEYA8 recall -0.4 and precision +0.6, with losses in
some corruptions. Yet the gate is default on. Containment-collapse ablations are
similarly mixed and remain off, correctly.

Keep synthetic-validated mechanisms available, but default-on promotion should
require an independent real-data acceptance gate with per-class recall, not only
aggregate F1.

### Medium 6: Canonical-only novel discovery excludes real noncanonical biology

The default canonical gate rejects every novel multi-exon transcript outside
GT-AG, GC-AG, AT-AC. GTF transcripts are trusted, so the bias applies
asymmetrically to novel discovery. Rare real noncanonical splice junctions,
reference errors, and sample-specific variants are unrepresentable.

Report canonical status as evidence and make strict removal a profile choice.
A high-support noncanonical chain should be allowed under a stronger read/signal
bar rather than categorically impossible.

### Medium 7 [Partially resolved]: `m2_em` exposes an iterative-EM interface but defaults to one softmax

A direct execution check confirmed:

```text
m3=False, abundance_feedback=False:
  iterations = 1
  max error versus one row-softmax = 0
  prior effect after update at beta=0 = 0
```

This is mathematically fine for the sparse rule matrix, but the names and exposed
knobs are misleading:

- `em_max_iter` and `em_tol` do no work on the default path.
- `em_sigma` is ignored; `m2_em` hard-codes sigma 1.
- `use_prior`, `prior_weight_cap`, and `score_alpha` are explicitly inert.
- `max_reads_per_interval_for_dtw` is unused.
- The M2 exact tie tolerance is hard-coded 1e-9; changing the exposed
  `m2_tiebreak_margin` does not widen the default tie set.
- `m2_tiebreak` itself is also ignored by default `m2_em`, so the intuitive
  CLI off-switch cannot perform a clean M2 ablation while retaining the mode.

Rename the default operation as calibrated sparse responsibility assignment, or
make the configuration match a real iterative model. Remove or reject inert
options rather than silently accepting them.

### Medium 8 [Resolved]: `PipelineConfig.validate()` was not on the production path

The CLI now validates the resolved config before runner setup and converts
validation failures to clean Click usage errors.

The validation function correctly checks paths and parallel invariants, and its
direct tests pass. `cli.main()` constructs the config and immediately calls
runner setup/run; no source/history call site invokes `validate()`.

Call it immediately after construction and translate errors to `click.UsageError`.
This is an implementation defect rather than a scientific issue, but it permits
invalid execution states.

### Medium 9: Observed-junction evidence can be incomplete without disabling gates

`compute_observed_junctions` fails open when BAM reading raises. Its own caveat
states that `BamReader.get_reads_in_region` may swallow a mid-iteration error and
return a partial non-empty list (`evidence.py:42-45`). Default support gates then
interpret undercounts as biological zero and may drop candidates.

The reader API should return completeness status. Hard negative evidence must
require a complete fetch; otherwise the gate should abstain.

### Medium 10: Signal model/profile metadata is insufficiently explicit

Krill and polyA default to RNA002. Pore model is a config field but is not a
normal CLI option in the 92-option interface, and input metadata is not validated
against it. Hidden environment knobs (`M1_MAX_INDEL_BP`, `MAPPY_R1_MIN_AS`,
`MAPPY_PRESET`, `KRILL_THREADS`) also change scientific results outside the
serialized config.

Promote all result-changing parameters to the run manifest, validate pore/signal
compatibility, and include them in reproducibility hashes.

### Medium 11: Fusion mode is an experimental hypothesis generator

Fusion clustering uses single-linkage within 500 bp on both breakpoints, so a
bridge read can join endpoints farther apart than 500 bp transitively.
Annotation arm variants are exempt from minimum support, and any annotation
participation can exempt a stitched combination from the normal read/read
intersection bar. This can expand candidate combinations beyond directly
observed fusion structures.

Fusion is off by default, which is appropriate. Before production promotion,
require direct breakpoint-spanning support per stitched structure, cap
single-linkage diameter, and validate precision on decoy/adaptor-rich reads.

### Medium 12: Family semantics still disagree with novel GTF gene IDs

`isoform_fraction_drops` now treats discovery-family siblings as isoforms, but
finalization still falls back to `gene_id=candidate_id` for every novel and
`write_gtf` groups records by that value (`finalize.py:73-80`,
`io_gtf.py:563-587`). Two family siblings can therefore compete as isoforms
during selection and then be emitted as two unrelated genes. This is not only a
presentation issue: downstream gene-level counting and interpretation see a
different ontology from the selection algorithm.

Assign a stable family/locus gene ID to novel siblings, or explicitly declare
that PyFIN does not perform novel gene grouping and stop using gene-like locus
semantics for its filters.

### Medium 13 [Resolved]: Recorded benchmark state had no reliable integration ledger

`PROFILE_OPTIMIZATION.md`, per-run source-hashed manifests, and preserved TSV
sweeps now connect evidence, resolved settings, production callers, and rejected
alternatives.

`PRODUCTION_STATE.md:83-91` mixes one landed change with four then-uncommitted
changes and says summed LLR is active in the runner. At current HEAD:

- span guard and mono locus splitting have landed;
- generation-time mono folding has been superseded by default post-EM mono
  resolution;
- shadow provenance has landed, but recovery only consumes it in cluster mode;
- summed LLR has landed as tested code, but not in default `m2_em`;
- `candidate_align`, `cluster_diff_regions`, and `explore_family_paths` are tested
  foundations with no production caller;
- the untracked HEYA8 research tree contains the latest signal/existence harness
  and 157 result/consultation artifacts under `diffsig/`.

There is no hidden branch or stash containing a complete integration. The extra
worktree is clean at an older detached refactor. `gpu_honesty.diff` is already
applied, while historical `wobble.diff` and `guided_gate.diff` have diverged and
cannot be treated as pending patches.

This process gap explains how a benchmark document came to describe a metric
that default execution never used. Maintain a versioned integration ledger with
one row per change: experimental evidence, implementation commit, production
caller, resolved defaults, threshold scale, superseding change, and holdout
validation status.

### Medium 14 [Resolved for real profiles]: Weak raw junction coordinates survived final selection

After mono/support selection, real p00 still contained 1,499 gffcompare `j`
models. Of these, 665 had the same intron count as a nearby reference chain with
all boundaries within 10 bp; for 568, every shifted truth junction was also the
strictly stronger local primary-read CIGAR mode. The correction evidence is
annotation-free even though GENCODE was used to label the audit.

Resolution: after global selection, real profiles scan primary BAM junctions
once. A novel multi-exon intron may move at most 6 bp only to a mode with at
least two reads and more than twice current-coordinate support; when canonical
gating is active, the target must pass the same motif set. GTF/fusion/mono
models are exempt. Fully identical corrected exon models merge with soft
abundance and hard read IDs conserved. Balanced T3 F1 improved by +0.229,
+0.567, and +0.146 on three H9 samples, with recall increasing in all three.
SIRV leaves this correction off.

### Low 1: Historical documents and source still disagree in other scientifically relevant ways

Examples:

- Older comments call argmax-first a production default; dispatcher uses
  argmax-keep and production defaults to m2_em.
- `record.md` calls primary M2/M3 assignment a dead end while HEAD still names
  m2_em as default, although the current path is a narrower exact-tie algorithm.

Generate an immutable run manifest and versioned algorithm description from code.
Do not use mutable narrative files as the only benchmark provenance.

## What is algorithmically sound

The review does not conclude that PyFIN is generally unreasonable. Several parts
are strong:

1. Coordinates and strand handling are explicit and consistently tested.
2. Family relations are simple, inspectable, and deterministic at the structural
   level.
3. Exact-subchain span guard encodes a biologically meaningful distinction:
   degradation does not span the missing intron; retained-intron evidence does.
4. GTF candidates cannot bridge read-derived families in the family primitive.
5. Assignment, selection, and output ownership have been separated into distinct
   modules, reducing accidental policy coupling.
6. PolyA scoring failure disables the filter instead of turning missing signal
   into negative evidence.
7. Exact-junction support protects real GTF siblings from synthetic wobble gates;
   the SIRV606 audit demonstrates why direct evidence is better than source labels.
8. Mono resolution is survivor-conditioned and empirically improved the real
   p00 operating point over generation-time folding.
9. The experiment suite often states its own SIRV/real-data caveats and includes
   honest expressed-truth denominators rather than scoring only easy transcripts.
10. Structural/basecall evidence is currently much stronger than private-region
    raw-signal evidence, and the repository is collecting the right diagnostics
    to discover that rather than hiding it.

## Overall judgment

| Capability | Judgment | Reason |
| --- | --- | --- |
| Structural candidate generation | Reasonable but incomplete | Strong de-fragmentation; endpoint states remain outside the model but explain only 1.38% of measured missed T3 multi-exon truth |
| SIRV/Sequin structure calling | Empirically strong | Six-sample mean/SUM test plus holdout; guide-conditional routing improves the aggregate |
| General real-dRNA defaults | Cross-sample validated and competitive | Mono support plus junction consensus improved three H9 samples; balanced mean T3 F1 38.883 vs IsoQuant 41.652 |
| M1 sequence evidence | Reasonable | Simple, inspectable, and empirically informative; exact-tie boundary is overly rigid |
| M2 raw-signal tie refinement | Integrated in a bounded niche | Real uses tight sum; SIRV routes by guide; same-count/two-score/invalid-read abstention constrain claims |
| Read-by-read coherence | Removed from production | M3 prototype retained under experiments; M4 has no production wiring |
| Abundance/TPM | Coherent on final survivors for production beta=0 M2 | Named profiles conserve each assignable read or report it as selection-orphaned; feedback/other modes remain unsupported |
| Fusion calling | Experimental | Off by default; support/cluster rules need precision validation |

The best concise description is:

> PyFIN is a competitive structural heuristic assembler with a narrow
> signal-assisted tie resolver and a mass-accounted final-survivor quantifier.
> It is not yet a unified, calibrated signal-generative transcript assembler.

## Prioritized remediation

### P0: correctness before further tuning

1. **Done:** remove secondary records from structural generation.
2. **Done:** M2 abstains unless at least two hypotheses score.
3. **Done:** mean/sum/off/auto routing, separate sum scale, and guide-aware SIRV A/B validation.
4. **Done:** balanced/precision real-dRNA and SIRV profiles with biological-domain default.
5. **Done:** isolate invalid eventalign payloads per read rather than poisoning a locus batch.
6. **Done:** replace UUID candidate IDs with stable structural BLAKE2b IDs.
7. **Done:** wire `PipelineConfig.validate()` into CLI execution.

### P1: make existence and abundance coherent

1. **Done for production beta=0 M2:** real mono/junction existence is validated and final survivors are requantified without re-filtering.
2. **Done:** mass-conservation, alignment-unassigned, selection-orphaned, forced-read, and TPM-shift diagnostics are written per run.
3. **Done:** propagate stable splice-family IDs and use them for isoform fraction.
4. Decide post-refit selection consistency. v9 requantifies the frozen survivor
   set, so the output is correct quantification under that set but is not a
   fixed point of the assignment-dependent gates. Offline audit of the absolute
   abundance floor found 0 violations on SIRV p00 and precision tuning and 2
   each on balanced tuning (0.018%) and replicate3 (0.013%), all novel mono
   models that lost a forced read. The isoform-fraction, soft-mass-ratio, mono
   hard-read, full-length, and poly(A) gates also consume refit-changed
   `num_reads`/assigned reads and remain unaudited, so the absolute-floor rate
   is not evidence that the full fixed-point gap is equally small. Thresholds
   were tuned against v8 diluted abundance and 7,252 of 11,281 balanced rows
   increased, so reusing them unchanged would silently alter filter strength.
   Compare freeze-once, a single post-refit filter, and a shrink-only
   fixed-point loop before promoting any of them.
5. Represent endpoint states inside a splice-chain family after higher-impact existence work; measured current impact is 1.38% of missed T3 multi-exon truth.
6. Calibrate candidate-specific M2 event-set likelihoods inside shared genomic
   disagreement windows; do not require one-to-one event correspondence.
7. Reconsider every default-on structural gate using independent real samples and
   per-isoform-class recall.

### P2: optional algorithms and reproducibility

1. Revisit the retired M3 prototype only with a missing-pair mask,
   leave-one-out coherence, and an enforced read cap.
2. Promote pore and hidden environment knobs into serialized configuration.
3. Consolidate or remove inert config and obsolete DTW/eventalign paths.
4. Harden fusion clustering/support and validate with decoys.
5. **Partially done:** source-hashed run manifest now records commit, resolved
   defaults, overrides, and environment; container/input hashes/backend identity
   remain to be added.

## Required validation program

### Unit/integration invariants

- One query with primary+secondary alignments contributes to exactly one family.
- Fusion-excluded read populations follow an explicit policy.
- One-scored M2 tie abstains; scored hypotheses use the same genomic
  disagreement windows while retaining hypothesis-specific event sets.
- A route-level test proves that mean, summed LLR, and off select the intended
  implementation; summed tests use tight diff windows and a calibrated threshold.
- Candidate selection followed by requantification conserves one unit per kept
  read (or reports deliberately unassigned mass).
- Stable IDs and byte-identical output across serial/parallel reruns.
- Endpoint-equivalent splice chains can retain distinct TSS/TES states.
- Isoform fraction never compares unrelated overlap-only loci.

### Statistical validation

- Freeze thresholds before final evaluation.
- Split by sample or chromosome so loci used to tune gates never score them.
- Keep SIRV/Sequin as corruption stress tests, not the sole production arbiter.
- For real data, combine cross-replicate reproducibility, short-read junction
  support, orthogonal TSS/TES evidence, and conservative annotation truth.
- Report metrics separately for GTF/novel, mono/multi, expression bins,
  canonical/noncanonical, containment, wobble, cassette, alternative ends, and
  fusion.
- For signal, report coverage/abstention, within-family pair accuracy, AUROC with
  cluster-bootstrap confidence intervals, and probability calibration.
- For abundance, report read-mass conservation, simulated mixture error, and
  rank/absolute error against spike-in concentrations.

## Advisor adjudication

The original adversarial pass highlighted M2 technical missingness, endpoint
loss, and abundance semantics. Follow-up execution resolved per-read missingness,
metric routing, family propagation, and production beta=0 post-selection
abundance refitting. Hypothesis-specific event sets inside a shared genomic
region remain intentional; SUM count-scale calibration is still open. The
endpoint audit narrowed current practical impact to 1.38% of missed T3
multi-exon truth, and M3 was removed from production.
