# v9 post-selection abundance refit — independent code review notes

Reviewed source hash `42e1210dea81860686e76aaf9fa058fc8a65808e0a31367570d8526f4da17954`.
No source change was made as a result of this review.

Shipped as commit `05e303ab99e04d1f6f783140c9eddb5a0df68973` on `dev`.

## Retired history backup

Before the `experiments/` data paths were stripped from 20 unpublished commits,
the full local history was preserved on branch `backup/dev-local-full` at
`ba4fb7a0c4d87cdc66af15ac50627c76e7f4e579`. After the rewritten history was
accepted and pushed, that branch was deleted and its unreachable objects were
garbage-collected to reclaim Git storage. The large files themselves were never
deleted; they remain on disk, untracked. Recording the SHA here keeps the
provenance chain readable even though the objects are gone.

## Verified correct

1. **Structural freeze is architectural, not incidental.** `min_gtf_abundance` /
   `floor_gtf_abundance` are referenced only in `fin/pipeline/selection.py`
   (pre-refit). `finalize_outputs()` explicitly drops nothing, and `write_gtf()`
   contains no filtering branch. Post-refit abundance therefore cannot change
   the emitted transcript set under any value.
2. **No global fold is missed.** The only `SelectionOutcome(action="fold")` in
   the codebase is `fin/pipeline/selection.py:273` (containment collapse inside
   `select_m2_interval`), which the ledger captures. `select_global()`
   (lines 324-539) contains zero `fold` occurrences, so discarding
   `_global_outcomes` in the runner is safe.
3. **Mono forced mass-1 is faithful to legacy semantics.**
   `mono_resolve_drops()` sets the parent `pq.abundance += 1.0` and the retained
   mono `mq.abundance = float(len(kept))`. Both folded and retained mono reads
   are mass-1 hard assignments in v8, so forcing mass 1 in the refit reproduces
   the legacy accounting rather than inflating it.
4. **Reads orphaned by a dropped mono parent are not force-pinned.** They have
   no post-selection owner, so they fall through to soft renormalization.
5. **`junction_snap` change is strictly additive.** `return_redirects=False`
   preserves the original return contract exactly.
6. **Config gating order is correct.** `profile_overrides` is populated at
   construction (`fin/cli.py:442`) before `validate()` (`fin/cli.py:551`).
7. **Programmatic misuse fails loudly.** Forcing
   `post_selection_refit_effective` without `m2_em` raises in
   `_order_interval_outputs()` rather than silently producing wrong numbers.

## Findings requiring narrowed wording (no numerical impact)

### A. Cross-interval reads are negligible but nonzero (resolved by measurement)

`extract_reads_for_interval()` fetches by region with **no strand filter** and
an overlap (not containment) test. Measured on the balanced tuning BAM using
production `generate_intervals_from_reads()` + `cluster_intervals()`
(15,190 intervals, matching the live run exactly):

| Intervals overlapped | Reads |
| ---: | ---: |
| 1 | 258,991 |
| 2 | 29,272 |
| 3 | 4,187 |
| 4 | 924 |
| >=5 | 37 |

**34,425 reads (11.73%) overlap two or more intervals.** That is an *upper
bound only*: geometric overlap means the read is fetched into more than one
interval's discovery, not that it reaches more than one `ResponsibilityLedger`.

The exact count was measured with `ledger_overlap_probe.py`, which monkeypatches
`fin.pipeline.runner.refit_survivor_abundance`, records ledger statistics, and
then delegates to the unmodified function. Production `fin/` source was not
touched. The probe run reproduced the release run exactly: 294,111 assignable
reads, 11,281 transcripts, `forced_target_conflicts = 3`, and `assembly.gtf` /
`scores.tsv` byte-identical to the release artifacts.

| Ledgers carrying the read | Reads |
| ---: | ---: |
| 1 | 294,100 |
| 2 | 9 |
| 3 | 1 |
| 4 | 1 |

**Only 11 reads (0.0037%) reach more than one ledger.** All 11 have differing
raw and effective candidate sets, none has identical normalized rows, and 8
differ in their surviving-candidate sets. The affected mass is bounded by 11 of
244,316 assigned units, i.e. 0.0045%.

Consequence: for the 294,100 single-interval reads, survivor renormalization is
exactly a beta=0 softmax restricted to surviving columns. For the 11
multi-interval reads the refit **sums the per-interval responsibility rows and
then renormalizes**, weighting each interval's opinion equally. Mass
conservation is unaffected (each read still contributes exactly 1.0; measured
per-read error <= 4.4e-16).

The correct claim is therefore: *survivor renormalization over merged local
responsibilities*, exact for single-interval reads and a documented
equal-weight merge for the 0.0037% carried by several intervals. It is not a
global EM or softmax rerun. `ALGORITHM.md` was updated to state this.

### B. Diagnostics counters are not a partition [FIXED]

`forced_reads + renormalized_reads + selection_orphaned_reads` exceeds
`assignable_reads` because a forced read whose target did not survive is counted
in both `forced_reads` and `selection_orphaned_reads`.

- real balanced: 26,653 + 224,045 + 49,795 = 300,493 vs 294,111 assignable
  (overlap 6,382)
- SIRV: 289 + 953 + 78 = 1,320 vs 1,310 assignable (overlap 10)

Mass balance is still exact. **Fixed:** `forced_orphaned_reads` now records the
intersection and `counter_semantics` documents it, so
`forced + renormalized + orphaned - forced_orphaned = assignable` holds exactly.
Measured live: SIRV 10, balanced 6,382, precision 4,657, replicate3 8,330.

### C. `alignment_unassigned_reads` is a scoped metric [FIXED]

It counts reads that entered interval quantification but received no
responsibility row. Intervals returning `None` early (no reads, no candidates)
contribute nothing, so it is not a BAM-global unaligned count. Real balanced
reports 851, whereas the BAM holds 295,252 distinct read names. **Fixed:**
`interval_quantification_unassigned_reads` is the scoped name;
`alignment_unassigned_reads` is retained as a compatibility alias.

### D. Audit mutates the artifact it audits [FIXED]

`abundance_refit_audit.py::audit()` writes `structural_identity: "passed"` back
into `abundance_refit.json` in place. The three-artifact determinism evidence
remains valid because both compared runs were stamped identically, but a
stamped run compared against an unstamped one would fail spuriously. **Fixed:**
the verdict is written only to `refit_identity.json`; `abundance_refit.json` is
never rewritten, and a regression test asserts byte equality across an audit.

### E. Stale comment in `fin/pipeline/parallel.py` [FIXED]

`run_parallel()` now states that aggregation is order-independent. Float
summation is order-sensitive, which is precisely why
`_order_interval_outputs()` exists. The comment could mislead a maintainer into
deleting the canonicalization. **Fixed:** the docstring now states that
aggregation is commutative but not bitwise associative and points at the
canonicalization.

### F. Dead statement [FIXED]

`fin/analysis/abundance_refit.py` updated `max_conservation_error` with a
constant `0.0` in the orphan branch; the statement had no effect. **Fixed:**
removed and replaced by the `forced_orphaned_reads` accounting.

## Empirical caveat worth surfacing

The largest single abundance shift is `novel_ab8d27c7bb3f7f2f`
(chrM:5590-10036, novel): `489.0 -> 8057.5`, i.e. `+7568.5`, which is 24% of the
entire 31,573.7 L1 abundance shift, leaving it at TPM 10,237. This is the
intended behaviour — it is the only surviving model for those reads — but it
shows the refit can concentrate very large mass onto one novel model at
high-coverage loci. No survivor ended at abundance 0 in either run.

## Open design question: should the gates rerun after the refit?

v9 requantifies a frozen survivor set. It is therefore *correct quantification
under that set*, not a *fixed point of the abundance gates*.

`post_refit_gate_audit.py` reapplies every assignment-dependent gate to the
refitted results inside a live run and reports what would be dropped, without
changing the run. All four audit runs produced `assembly.gtf` and `scores.tsv`
byte-identical to the release artifacts, proving zero perturbation. Full-length
support is recomputed from post-refit `assigned_read_ids`; cached `fulllen_frac`
is stale evidence and is not used.

| Profile | Survivors | Floor | Isoform fraction | Soft-mass | Mono | Full-length | Union | Rate |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| SIRV p00 | 81 | 0 | 0 | 0 | n/a | 0 | 0 | 0.000% |
| balanced tuning | 11,281 | 2 | 0 | n/a | 0 | n/a | 2 | 0.018% |
| balanced replicate3 | 15,656 | 2 | 3 | n/a | 0 | n/a | 5 | 0.032% |
| precision tuning | 10,217 | 0 | 0 | n/a | 0 | n/a | 0 | 0.000% |

`n/a` means the gate is inactive in that profile. `min_polya5p_reads` is 0
everywhere, so it is reported as active=false/unevaluated rather than as zero
violations; it needs the signal-derived polyA map to audit.

Every violator is a **novel multi-exon** model, not a mono model. The floor
violators all hold a single read at abundance `1.00` (balanced `chr22+`, 3
exons; replicate3 `chr7-` 3 exons and `chr10-` 9 exons) and fell from `2.00`
because a forced mono read moved to another parent. The replicate3
isoform-fraction violators are larger models (`chr15+` 4 exons at `20.0`,
`chr1-` 12 exons at `2.0`, `chr15+` 3 exons at `23.0`) whose family sibling grew
enough to push them under the 1% relative floor — a relative, not absolute,
effect. No candidate violates two gates at once.

Direction of change is otherwise strongly one-way: 7,252 balanced rows
increased, 5 decreased, 4,024 unchanged.

Simulated drop (remove the violators, rescore with the same gffcompare/truth
contract):

| Profile | Transcripts | T3 F1 | Delta | Released mass |
| --- | --- | --- | ---: | --- |
| balanced tuning | 11,281 -> 11,279 | 39.422 unchanged | 0.0000 | 2.0 of 244,316 |
| balanced replicate3 | 15,656 -> 15,651 | 38.6018 -> 38.5991 | -0.0027 | 47.0 of 485,404 |

So enforcing the fixed point is **not an accuracy win**: it is a consistency
change that costs a small amount of recall on replicate3 (2 of the 5 dropped
models matched truth) and changes nothing on balanced.

Why a naive rerun is wrong: every threshold was tuned against v8's diluted
abundance. Since refit almost always raises abundance, reusing the same numbers
changes effective filter strength. Structural/evidence gates (canonical,
junction evidence, containment, M2) do not depend on final abundance and should
not rerun.

### Measured with `post_refit_fixed_point_probe.py`

The union audit above is not by itself proof, because production
`select_global()` applies gates *sequentially* (each sees the previous gate's
survivors) and because "released mass" is not the same as an orphan-mass change
until the set is actually re-refitted. The fixed-point probe closes both gaps:
it applies the gates in production order and performs a real refit after every
shrink.

| Profile | freeze_once | once | fixed_point | Shrink rounds | Orphan mass delta |
| --- | ---: | ---: | ---: | ---: | ---: |
| SIRV p00 | 81 | 81 | 81 | 0 | 0.0 |
| balanced tuning | 11,281 | 11,279 | 11,279 | 1 | +2.0 |
| precision tuning | 10,217 | 10,217 | 10,217 | 0 | 0.0 |
| balanced replicate3 | 15,656 | 15,651 | 15,651 | 1 | +47.0 |

Results:

* **`once` and `fixed_point` produce identical survivor sets on all four
  profiles.** There is no cascade: the round after any drop removes nothing, so
  the loop converges immediately.
* Sequential gating selected exactly the same candidates as the independent
  union audit, so the earlier simulated-drop F1 numbers are valid for the true
  fixed point.
* **Released mass is not absorbed by surviving siblings; it becomes orphan
  mass.** Balanced assigned mass falls 244,316 -> 244,314 while orphan rises
  49,795 -> 49,797; replicate3 falls 485,404 -> 485,357 while orphan rises
  109,027 -> 109,074. Every released unit is accounted for, and the mass simply
  moves from "assigned" to "orphaned". This is the behaviour predicted for
  forced mono reads whose parent is deleted: they orphan immediately instead of
  flowing to another model.
* Structural F1 does not improve: balanced is unchanged and replicate3 loses
  `0.0027` (2 of its 5 dropped models matched truth).

### Decision

Against the pre-registered rule - drops confined to 2/5 transcripts, at most one
shrink round, no F1 gain, negligible orphan change - **production keeps
`freeze_once`.** No `post_refit_filter_mode` is added. Enforcing the fixed point
would delete a handful of real transcripts, convert their reads into orphan
mass, and slightly reduce recall, in exchange for internal tidiness alone.

The residual inconsistency is documented rather than fixed: the emitted set can
contain a small number of models that no longer satisfy the abundance or
isoform-fraction gate under refit accounting (0.000-0.032%, all novel
multi-exon). `post_refit_fixed_point_probe.py` is retained so the claim can be
re-measured whenever gates or thresholds change.

## Recommendation

Items A-F are wording, diagnostics, and provenance issues. None changes a
number, a structure, or an acceptance verdict. Item A was resolved by
measurement and corrected in `ALGORITHM.md` (documentation only).

Items B-F were fixed in a follow-up cleanup batch after `05e303ab` shipped,
rotating the source hash once to
`90c8117a958a2b0653ecef1235d6521ebc0a614315f5ea6e8cc55f1909e63eeb`. The cleanup
is provably numerically inert: SIRV p00, balanced tuning, precision tuning, and
replicate3 all reproduce `assembly.gtf` and `scores.tsv` byte-identically to the
`05e303ab` release. Only `abundance_refit.json` changes, by gaining
`forced_orphaned_reads`, `interval_quantification_unassigned_reads`, and
`counter_semantics`.

No correctness blocker was found. The two semantic checks requested during
review both cleared:

1. `select_global()` cannot emit a fold, so discarding `_global_outcomes` loses
   no redirect.
2. Multi-ledger reads exist but are 0.0037% of the population, bounding any
   deviation from restricted-softmax equivalence at 0.0045% of assigned mass.
