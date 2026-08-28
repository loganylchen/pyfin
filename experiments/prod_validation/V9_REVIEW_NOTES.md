# v9 post-selection abundance refit — independent code review notes

Reviewed source hash `42e1210dea81860686e76aaf9fa058fc8a65808e0a31367570d8526f4da17954`.
No source change was made as a result of this review.

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

### B. Diagnostics counters are not a partition

`forced_reads + renormalized_reads + selection_orphaned_reads` exceeds
`assignable_reads` because a forced read whose target did not survive is counted
in both `forced_reads` and `selection_orphaned_reads`.

- real balanced: 26,653 + 224,045 + 49,795 = 300,493 vs 294,111 assignable
  (overlap 6,382)
- SIRV: 289 + 953 + 78 = 1,320 vs 1,310 assignable (overlap 10)

Mass balance is still exact. A `forced_orphaned_reads` key would remove the
ambiguity.

### C. `alignment_unassigned_reads` is a scoped metric

It counts reads that entered interval quantification but received no
responsibility row. Intervals returning `None` early (no reads, no candidates)
contribute nothing, so it is not a BAM-global unaligned count. Real balanced
reports 851, whereas the BAM holds 295,252 distinct read names.

### D. Audit mutates the artifact it audits

`abundance_refit_audit.py::audit()` writes `structural_identity: "passed"` back
into `abundance_refit.json` in place. The three-artifact determinism evidence
remains valid because both compared runs were stamped identically, but a
stamped run compared against an unstamped one would fail spuriously. Writing
only `refit_identity.json` would remove the coupling.

### E. Stale comment in `fin/pipeline/parallel.py`

`run_parallel()` now states that aggregation is order-independent. Float
summation is order-sensitive, which is precisely why
`_order_interval_outputs()` exists. The comment could mislead a maintainer into
deleting the canonicalization.

### F. Dead statement

`fin/analysis/abundance_refit.py` updates `max_conservation_error` with a
constant `0.0` in the orphan branch; the statement has no effect.

## Empirical caveat worth surfacing

The largest single abundance shift is `novel_ab8d27c7bb3f7f2f`
(chrM:5590-10036, novel): `489.0 -> 8057.5`, i.e. `+7568.5`, which is 24% of the
entire 31,573.7 L1 abundance shift, leaving it at TPM 10,237. This is the
intended behaviour — it is the only surviving model for those reads — but it
shows the refit can concentrate very large mass onto one novel model at
high-coverage loci. No survivor ended at abundance 0 in either run.

## Open design question: should the gates rerun after the refit?

v9 requantifies a frozen survivor set. It is therefore *correct quantification
under that set*, not a *fixed point of the abundance gates*. Offline audit of
the absolute abundance floor only:

| Profile | Rows | Floor violations | Rate |
| --- | ---: | ---: | ---: |
| SIRV p00 (`>=3.0`, GTF floor raised to 3.0) | 81 | 0 | 0.000% |
| balanced tuning (`>1.0`) | 11,281 | 2 | 0.018% |
| balanced replicate3 (`>1.0`) | 15,656 | 2 | 0.013% |
| precision tuning (`>1.0`) | 10,217 | 0 | 0.000% |

The four balanced/replicate3 violations are novel mono models that fell from
`2.00` to `1.00` because a forced mono read moved to another parent. Direction
of change is otherwise strongly one-way: 7,252 rows increased, 5 decreased,
4,024 unchanged.

This rate covers the absolute floor only. `min_isoform_fraction`,
`max_soft_mass_ratio`, the mono hard-read floor, `min_fulllen_fraction`, and
`min_polya5p_reads` also consume refit-changed `num_reads`/assigned reads
(3,707 balanced rows changed `num_reads`) and have not been audited, so the
full fixed-point gap is not established.

Why a naive rerun is wrong: every threshold was tuned against v8's diluted
abundance. Since refit almost always raises abundance, reusing the same numbers
changes effective filter strength. Structural/evidence gates (canonical,
junction evidence, containment, M2) do not depend on final abundance and should
not rerun.

Registered as the next independent batch: audit every assignment-dependent gate
offline across all profiles, then compare freeze-once, a single post-refit
filter, and a shrink-only fixed-point loop (candidates only leave, so it
terminates, but it is path-dependent and cannot recover a first-round
mis-drop). Do not widen v9 to include it.

## Recommendation

Items A-F are wording, diagnostics, and provenance issues. None changes a
number, a structure, or an acceptance verdict. Item A was resolved by
measurement and corrected in `ALGORITHM.md` (documentation only). Fixing B-F
touches `fin/` source and would rotate the release hash, invalidating every
live acceptance artifact, so they belong in a follow-up batch rather than this
release.

No correctness blocker was found. The two semantic checks requested during
review both cleared:

1. `select_global()` cannot emit a fold, so discarding `_global_outcomes` loses
   no redirect.
2. Multi-ledger reads exist but are 0.0037% of the population, bounding any
   deviation from restricted-softmax equivalence at 0.0045% of assigned mass.
