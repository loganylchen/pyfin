# Making the fusion analysis work on real direct RNA

Date: 2026-09-15. Five Codex (gpt-6-astra) rounds. The positive control is
BCR-ABL1 in K562 direct RNA (SG-NEx replicate4 run1), the textbook fusion.

## The result that started this

The independent fusion analysis, validated on synthetic error-free reads
(`FUSION_VALIDATION_2026-09-12.md`) and extended with the fused-transcript
layer (`FUSION_TRANSCRIPT_MODELS_2026-09-15.md`), found **nothing** on the
real positive control. A locus-restricted BAM (98 records, 11 nominated
molecules, 7 SA-tagged chr22↔chr9) against the full GRCh38 with GENCODE v44
produced **0 chains at the defaults**. At a relaxed `min_read_coverage=0.8`
it produced ONE chain, `insufficient_full_chain_support`, gene path
`[[BCR],[]]` — the ABL1 arm unlabelled and the selected structure wrong.

Lowering the coverage threshold "worked" in the sense that it rescued a
wrong explanation. The four real defects were these.

## 1. An SA tag is a locator, not a geometry

minimap2 writes an SA CIGAR as `min(query span, reference span)M` plus one
`I`/`D` for the difference: every intron and indel position is discarded.
The BCR arm of read `97ea8286` is a spliced `...2M717N26M...` in its
supplementary record and `175M723D` in the SA tag — one block spanning 898
reference bases. `segment_from_cigar` splits blocks on `N` only, so the SA
summary entered the analysis as an independent observed structure, competing
with (and beating) the real one.

No deletion-length rule can repair this: the intron's *position* is gone, so
turning `723D` into `723N` would place a splice after all 175 aligned bases.

`reconcile_sa_locators` (`fin/fusion/records.py`) removes an SA-only geometry
when a detailed alignment already represents the same placement: same read,
chromosome, strand and read length; reciprocal query overlap ≥ 0.9; query
endpoints within 10 bp; outer reference endpoints within 10 bp. Query overlap
alone is not enough — it would collapse a same-chromosome paralogous
placement onto the wrong locus. Every distinct *detailed* geometry survives;
nothing is chosen by MAPQ, source or input order. An SA entry with no
corroboration is kept as it is: it still nominates the molecule and the
search windows. Supersessions are recorded per read with their correspondence
metrics (`sa_supersessions` in `fusion_reads.jsonl`).

## 2. A deletion-free anchor rejects nanopore reads

The old rule required every one of the 12 template positions on each side of
a boundary to be aligned. A direct-RNA read with a 2-bp deletion five bases
from a real junction failed it. The correct spliced BCR-ABL1 path was
rejected as `missing_junction_anchor` at query coverage 0.939.

The rule is now, per side, over a window of `2 × min_anchor`: at least
`min_anchor` template positions aligned, at least `ceil(min_identity × aligned)`
of them matching, and no more than 4 consecutive unaligned template positions
*immediately adjoining* the boundary. Deletions deeper in the flank
contribute nothing but do not disqualify; they remain subject to the global
indel and template-run limits. Note `min_identity × aligned`, not
`min_identity × min_anchor` — with 24 aligned positions the latter would
accept a 37.5%-identity flank.

`validate_stitched_path` now records `(aligned, matched, required, edge gap)`
per flank in `anchor_detail`, so a rejection says what it failed on.

## 3. A template can carry sequence the read never aligned

The wrong unspliced path **validated** as `whole_transcript_anchored`. Its
template contained the 717-bp intron; the read simply failed to align exon
13, and nothing checked that. Boundaries come from block transitions, so an
invented span *inside* one block was invisible.

Every template block must now be represented by the same forward alignment:
at least `min(min_anchor, L)` aligned positions and at least 70% of its
length; no unaligned template run longer than `max_indel` inside the hit; no
unrepresented template prefix or suffix longer than `max_indel`. Per-block
counts are published in `block_support`. This does not ban intron retention —
a retained intron supported throughout by aligned query sequence passes; an
unsupported one does not.

## 4. Direct-RNA molecules are not genomic end to end

These reads carry ~80 bp of non-genomic 3′ sequence (adapter and poly(A)
carry-through; in this library the tail is not even A-rich) plus a few 5′
bases. `min_read_coverage=0.9` over the full read rejects genuine short
fusions for sequence that is not supposed to map.

`coverage_domain` computes, once per molecule, the terminal bases that **no
detailed alignment could place** — judged on the union of every detailed
mapping including low-MAPQ, normal and paralogous alternatives, so an
independently placed terminal arm is never excluded. Bounds: ≤100 bp per end,
≤120 bp total, ≤20% of the molecule; a whole eligible overhang or nothing,
never a slice of a longer one, never anything internal. The same domain is
used by the path search and the validator, so every competing hypothesis is
judged on one denominator. Raw and effective coverage are both reported.

Separately, `find_read_paths` no longer discards paths that cover a few bases
fewer than the best one. Covered bases are not a substitute for the sequence
adjudication that follows, and a slightly longer wrong geometry must not
erase the right hypothesis.

Nomination: `min_clip_bases` 100 → 250. At 100 the adapter/poly(A) tail
nominates **8.6% of a real library** (85,528 of 993,703 primaries on K562
rep4), which trips the 100k cap; 250 leaves 17,287. SA and supplementary
records still nominate regardless of clip length.

## The consequence for the transcript layer

The chain layer keys on exact topology, so five molecules of one fusion that
lost different numbers of 5′ exons are five chains of support 1 — nothing is
"supported" at `min_support=2`. At the transcript level they are one
transcript observed to different extents.

Models are therefore **aggregate envelopes**: structures whose middle arms and
breakpoints agree and whose outer arms stand in a 5′-truncation / 3′-truncation
relation are unioned, taking the longer extent on each side. A model is
reported when at least `min_support` distinct molecules fit inside its
envelope. Because no single molecule need have observed the whole envelope,
every model carries `envelope_status` (`single_molecule_spanned` when one
observation reaches both ends, else `aggregate`) and the per-end
`five/three_prime_observed_read_ids`.

Membership is **exact**. A molecule whose corrector decisively settled on a
placement 2 bp away is recorded as `nearby_read_ids` — visible, not counted.
Coordinate proximity is not evidence. Two models merge only when (a) some
molecule is an exact member of both, i.e. it holds both placements itself,
and (b) the whole group stays within tolerance pairwise across all members —
a bounded diameter, not single linkage, so a bridge molecule cannot chain two
distinguishable populations. The merged record is the placement annotation
agrees with on the most arms, then the longest; every member's complete
placement, and every alternative any member already carried, is listed in
`coordinate_alternatives`. A regenerated envelope gets a new id (aliases are
resolved transitively), re-resolved annotation, recomputed endpoint reach and
end statuses, and keeps `unresolved_microhomology_same_sequence` only when
the record *and* every member spell one sequence.

## The positive control now

Frozen as `tests/data/fusion_k562_bcr_abl1/` — 95 real records (unmodified
sequences, SA tags including `175M723D`, hard clips, indels, MAPQ) on two
300-kb reference slices with every BCR/ABL1 GENCODE v44 transcript whole,
1.3 MB, seconds to run. `experiments/fusion/build_bcr_abl1_fixture.py`
rebuilds it; reads whose SA partner lies outside the two loci are left out
entirely rather than stripped of their partner.

| | before | after |
| --- | --- | --- |
| chains at defaults | 0 | 7 |
| fused-transcript models | — | **1** |
| junction | — | **BCR:e14>ABL1:e2** (K562's canonical b3a2) |
| arm classes | — | `arm_fsm` / `arm_fsm` |
| transcripts | — | ENST00000305877 / ENST00000359540 → ENST00000318560 / ENST00000372348 |
| length | — | 7172 bp |
| unique support | — | 4 molecules |
| nearby (reported, not counted) | — | 1 |
| envelope | — | `single_molecule_spanned` |

The 5th molecule holds a decisive 2-bp-shifted donor and is `nearby`. The
6th has an unannotated chr9 exon (support 1, not reported). The 7th — the
`175M723D` read — still fails, and the reason is now legible: its left flank
upstream of the first BCR splice boundary has 23 aligned template positions
but only 15 matching, against `ceil(0.7 × 23) = 17`. That is local identity,
not coverage, not template representation, not the fusion edge.

## Throughput

The per-molecule genome search was single-threaded: mappy holds the GIL in
`map()`, and `n_threads` only speeds the index build. A full-library run sat
at 17k molecules for 7 hours. `FusionGenomeAligner.map_many` builds the index
once, forks a pool, and each worker maps chunks with the inherited
copy-on-write index; results return keyed by read id and are consumed in the
caller's deterministic order. A test asserts serial and 3-worker results are
equal per read.

Operational limits, stated rather than hidden: not thread-safe or reentrant
(`_POOL_ALIGNER` is process-global); fork-based, so Linux only; an ordinary
worker exception propagates and the pool is cleaned up, but abrupt native
worker death would leave the blocking map awaiting a task — crash supervision
is separate work; workers touch only the mappy index, never the genome
mapping, and one that needs the genome must reopen it. Codex verified that a
forked child closing its inherited flock descriptor cannot release the
parent's lock (same open file description; only an explicit `LOCK_UN` would),
so the publication lock is safe under this ownership.

## What this does not claim

The two-locus fixture proves the regression, not genome-wide uniqueness of
the mapping. Windowed anchors deliberately tolerate more local error than the
previous exact rule: there is no internal-domain bypass for an adapter
chimera (an internal adapter cannot be clipped out of one contiguous hit
while both arms are retained, and a >10 bp insertion fails `large_indel`),
but these rules are not an adapter detector, and adapter-bearing negative
controls are still owed. `supported_rna_fusion` remains RNA-level chain
evidence, not a demonstrated genomic rearrangement.
