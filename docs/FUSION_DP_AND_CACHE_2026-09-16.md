# Making the fusion adjudication 4x faster without changing a byte

Date: 2026-09-16. The question was the owner's: the per-molecule adjudication
was 6,277 s of a 6,631 s whole-library run — can an existing aligner replace
the hand-written DP, specifically mappy, which the codebase already depends on?

Result: **6,631 s → 1,643 s (4.0x) end to end, all 13 published artifacts
byte-identical**, and mappy is the wrong tool for reasons that have nothing to
do with speed.

## Two wrong guesses first

Worth recording, because both were structurally plausible and both were wrong.

**Guess 1: chain construction.** `build_fusion_chains` has an
`O(clusters × eligible)` subpath sweep. Measured: **1 second**, on 1,874
eligible paths and 595 clusters. Not the problem.

**Guess 2: a per-molecule cost model.** I added a timer per molecule and a
batch counter, and they contradicted each other — the counter said 22,000 of
23,403 molecules in 8 s, the timer emitted more records than the loop had
iterations left, with read ids in sorted order. I reverted both rather than
publish numbers I could not explain.

The explanation was not deep. `analyze_fusion_input` has **two** loops. Loop A
(`for index, molecule in enumerate(evidence.molecules)`) does reconcile,
coverage domain, path search and annotation for all 23,403 molecules in 8
seconds. Loop B (`for read_id in sorted(searches)`) does correction, proposal
construction and `validate_stitched_path` for the ~2,400 molecules that
produced candidate paths. I had put the timer's start in loop A and its end in
loop B, so every "per-molecule" duration was wall-clock since loop A finished,
and the sorted read ids were loop B's iteration order. Both probes were correct
about their own loop.

Lesson kept: a probe that disagrees with another probe is data about the probes.

## Finding 1: the genome cache was thrashing

`LazyGenomeFasta.__getitem__` returns a **whole chromosome** as a Python
string, and its LRU held **two**. That default is right for the assembly
pipeline, which sweeps one locus at a time. The fusion adjudication walks
candidate molecules in read-id order, so the chromosomes it needs arrive in
effectively random order.

Instrumented on a whole-genome run:

```
hits 137,042   misses 499   evictions 497
fetched_bases 66,637,687,660      fetch_seconds 246
names_seconds 0.003
```

**497 evictions on 499 misses** — nearly every miss threw out what the next one
wanted. 66.6 Gb fetched from a 3.1 Gb genome: the genome read 21 times over.
(`_names_set()`, rebuilt per miss, was innocent: 3 ms total.)

This also explains a puzzle from the acceptance record: why the two-locus
fixture and the chr9+chr22 run were fast while the whole genome was not. With a
two-chromosome reference **the cache cannot miss**.

The fix bounds the cache by **size** rather than count. A count bound cannot be
chosen safely without knowing the reference — 2 is right for one access pattern
and catastrophic for another, and 25 would be 3.1 GB on GRCh38 and far worse
elsewhere. A byte budget means the same thing on a spike-in FASTA (everything
fits) and on a genome larger than the budget. `cache_chroms=None` bounds by size
alone; both bounds apply when both are given; every existing caller's default is
untouched. An entry larger than the whole budget is **kept** — evicting it on
insert would make every access a miss, which is worse than no cache.

The tests state the limit too: a budget smaller than the working set thrashes
exactly as the count bound did. A size bound is not a better eviction policy; it
is a bound a caller can actually choose.

## Finding 2: mappy cannot replace the DP, and score equality is not a gate

Consulted Codex (gpt-6-astra), which read the module and ran cases through it.

**What `align_junction_window` actually is:** a three-state affine-gap
semi-global DP where the *query* is consumed end to end, the *target* is free at
both ends, and the traceback is **side-aware** — it partitions aligned and
matched base counts into left and right of a given junction offset, and counts
inserted bases, deleted bases and the longest single indel run.

**Why mappy fails, first thing first:** the function must return an optimal
alignment for *every* valid window, including one with no matching seeds. Codex
ran `query = "A"×75` against `target = "C"×75`, junction 35: score −154.0, 75
inserted bases, zero anchors, reason `large_indel`. A seed-chain-extend mapper
produces no hit there, and lowering seed length does not turn mapping into
exhaustive DP. Requiring `q_st == 0` and `q_en == len(query)` *filters* mapper
hits; it does not optimise over query-global, target-free alignments. Rescoring
a returned CIGAR with our formula does not establish it is optimal under that
formula. The existing mappy use in `validate_stitched_path` has a different
contract — it evaluates hits, tolerates partial coverage, and may return
`no_forward_template_alignment`.

**The sharper point:** score equality is not a sufficient gate for *any*
replacement. Two equally optimal intervals, `[3,5)` and `[2,4)`, carry
left/right anchor counts `0/2` and `1/1` — one rejects with
`missing_junction_anchor`, the other passes. A score-preserving replacement can
change the verdict.

Two contract details are easy to miss and both are load-bearing: a gap of length
L costs `gap_open + L*gap_extend` (4+2L, **not** 4+2(L−1)), and the recurrence
permits direct transitions between opposite gap states, charging another
opening. The contract is not "affine semi-global alignment"; it is **this
recurrence with this canonical traceback**.

Codex's ranking: compile the recurrence; then provably-exact shortcuts and
reuse; then dependency-correct numpy; then trial parasail with traceback
compatibility as a *blocking* requirement; reject mappy and edlib.

**Constraint this repo adds:** the production SIF has numpy and mappy (both
declared dependencies) but **not** parasail, Biopython or edlib, and numba is
present without being declared. So options 1 and 4 both mean a new dependency
and a rebuild of the frozen validation image.

## What was measured before being built

Codex's instruction was to measure the real call distribution, not guess.
`experiments/fusion/capture_junction_calls.py` records every call — and doubles
as the differential corpus.

| | BCR-ABL1 fixture | chr9+chr22 vs whole genome |
| --- | ---: | ---: |
| calls | 271 | 4,660 |
| exact repeats | 29.2% | **38.4%** |
| share of alignment time | 29.4% | — |
| query is an exact substring of target | 4.1% | — |
| DP cells, median | 5,320 (~73×73) | — |

That **killed Codex's own exact-substring shortcut** on value (4%) and justified
memoisation. Nothing else in the ranking is reachable without a new dependency.

## The change: memoisation, which needs no equivalence argument

The core is a pure function of the normalised window and the scoring parameters,
returning a frozen dataclass. The key is the **normalised window plus every
result-affecting parameter**, never the caller's junction coordinates — Codex's
trap: joint validation re-evaluates a junction with `tolerance=0` where the
first pass used `3`, and since padding is `flank + tolerance` the template
flanks go from 38 bases to 35. Same coordinates, different DP input, different
answer. Keying on the actual window covers this by construction.

`max_window` is deliberately *not* in the key, which is only safe because its
rejection happens before the cache is reached — asserted by a test.

## Numbers

Controlled profile, same job each time (chr9+chr22 molecules against the whole
genome, `--realign-workers 1` so cProfile sees everything):

| | baseline | + cache | + memo |
| --- | ---: | ---: | ---: |
| `analyze_fusion_input` | 667.2 s | 369.6 s | **213.4 s** |
| `align_junction_window` | 291.7 s | 259.4 s | **119.3 s** |
| `lazy_genome.__getitem__` | 275.8 s self | out of top 20 | — |
| `builtins.max` | 79.5M calls | 79.5M | **51.4M** |

Whole library, end to end: **6,631 s → 1,643 s (4.0x)**, wall clock 28 min 16 s.

Verification:
- **All 271 captured real calls replay field-for-field identical**, 0 mismatches.
- **All 13 published artifacts byte-identical** on chr9+chr22 *and* on the whole
  library, against the pre-optimisation run.
- 2,460 tests pass — and the suite itself went from 273 s to 66 s, because tests
  repeat windows too.

Codex then caught the capture harness measuring itself: run after the memo
landed, it would time a dict lookup for exactly the repeated windows whose
avoided work it exists to measure, understating the saving where it is largest.
It now clears the cache before each timed call. Cold-timed it reports 29.4%,
confirming the pre-memo 30.4% the decision was made on.

## What is left, and the decision it needs

The DP is now **55% of the remaining adjudication**. Getting more from it means
compiling the recurrence, which in this repo means adding numba or parasail and
rebuilding the SIF that HANDOFF section 5 freezes as the validation environment.
That is an owner decision, not a code decision.

Amdahl bounds the prize: at the pre-memo measurement, even an *infinitely* fast
alignment took the phase only from 370 s to 110 s. The memo has since taken it
to 213 s, so the remaining headroom is smaller than that.

Also unfixed by this work: `analyze_fusion_input`'s loop B is still serial.
Codex's structure for parallelising it — fork after the three junction
registries are complete, each worker reopening its own genome handle, results
reassembled in sorted read-id order — is recorded but not built. The genome
cache fix makes it cheaper than it would have been: forked children share large
chromosome strings copy-on-write.

Separately fixed while measuring: the elapsed clock started *after*
`GTFReader.parse()`, so the pipeline reported `annotation indexed [0s]` for tens
of seconds of real work and every later phase quoted a time the wall clock
disagreed with (1,643 s reported against 1,696 s measured outside the process).
