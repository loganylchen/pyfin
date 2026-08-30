# PyFIN Algorithm

Snapshot: `dev` at `a2c4be4e` plus the verified profile-integration worktree,
source SHA-256
`42e1210dea81860686e76aaf9fa058fc8a65808e0a31367570d8526f4da17954`,
2026-08-26. This document describes what the current source actually executes.
Older README, production-state, and experiment notes are evidence but do not
override the live source and run manifests.

## 1. Problem formulation

PyFIN takes nanopore direct-RNA reads and constructs a set of transcript
hypotheses. It then assigns reads to those hypotheses, removes candidates that
look like alignment/degradation artifacts, and reports transcript structure and
abundance.

Inputs:

- A spliced BAM/SAM alignment to the genome.
- An optional GTF/GFF annotation.
- A genome FASTA.
- Basecalled read sequences in FASTQ/BAM.
- Optional SLOW5/BLOW5 or POD5 raw signal.

Outputs:

- Assembled known, novel, and optional fusion transcripts in GTF.
- Soft abundance, hard read count, confidence, TPM, and diagnostic scores in TSV.
- Optional fusion breakpoints in BEDPE.

The implementation has three conceptually separate decisions:

1. `Generation`: which transcript structures are allowed to exist as candidates.
2. `Existence/selection`: which generated structures survive evidence gates.
3. `Abundance`: how read mass is divided among surviving or competing structures.

Current production code does not solve all three as one joint probabilistic
model. It uses structural generation, a narrow signal-assisted assignment, and
sequential heuristic selection.

## 2. Coordinates and core state

Internal genomic spans are 0-based, half-open `[start, end)`. GTF is converted
at the I/O boundary to/from 1-based inclusive coordinates.

Core objects (`fin/candidates/dataclasses.py`):

- `IntronChain`: immutable tuple of genomic `(donor, acceptor)` introns.
- `TranscriptCandidate`: ID, chain, 3-prime coordinate, spliced sequence, source,
  discovery reads, strand/span, persistent splice-family ID, and optional fusion breakpoints.
- `CandidateSet`: one genomic interval's candidates, read IDs/sequences,
  generation clusters, folded-shadow provenance, and read spans.
- `QuantResult`: soft abundance, hard assignment count/IDs, confidence, source,
  exon structure, persistent family ID, scores, full-length fraction, and fusion metadata.

Candidate source is exactly `gtf`, `novel`, or `fusion`.

## 3. Configuration profiles

The CLI resolves `--profile real-drna|real-drna-precision|sirv|custom` before
constructing `PipelineConfig`. Balanced `real-drna` is the product default. Explicit command-line
values override profile values, including boolean `--no-*` flags. Python callers
use `PipelineConfig.from_profile`; direct construction retains programmatic
backward-compatible defaults.

Shared named-profile settings:

```text
chain_cluster_discovery = True
clustering              = families
quant_mode               = m2_em
canonical_gate           = True
m2_diff_cover_gate       = True
m2_cluster_recheck       = True
novel_junction_min_reads = 2
m2_support_gate          = True
containment_cluster      = True
mono_resolve_post_em     = True
abundance_feedback       = False
min_gtf_abundance        = 1.0
min_isoform_fraction     = 0.01
isoform_fraction_locus   = family
post_selection_refit     = True (named profiles)
use_gpu                  = True
threads                  = 1
gpu_workers              = 0
```

Domain-specific values selected by the recorded optimization:

| Setting | `sirv` | `real-drna` |
| --- | ---: | ---: |
| `m2_metric` | auto: SUM unguided, mean guided | `summed_llr` |
| summed margin / flank | 1 / 4 | 1 / 8 |
| novel abundance boundary | inclusive >=3.0 | strict >1.0 |
| `floor_gtf_abundance` | true | false |
| `max_soft_mass_ratio` | 2.0 | 0.0 |
| `min_fulllen_fraction` | 0.1 | 0.0 |
| `min_polya5p_reads` | 0 | 0 |
| novel mono hard-read floor | off | 5 |
| read-supported junction snap | off | tolerance 6, support 2, ratio 2 |

`real-drna-precision` inherits balanced real settings and raises generation-time
`min_novel_reads` from 1 to 2. Across three real samples this gains about 1.1 T3
F1 beyond balanced selection but costs about 1 pp T1/T3 recall. `custom` applies
no profile overlay; pre-family historical reproduction additionally requires
`--isoform-fraction-locus overlap`. Every CLI run validates paths/settings and atomically writes
`run_manifest.json` with the resolved config, explicit overrides, source hash,
Git commit when available, and result-changing environment variables.

## 4. Interval construction

Entry: `fin/io/interval_manager.py:generate_isolated_intervals`.

1. Scan BAM alignments.
2. Skip unmapped and secondary alignments while constructing interval seeds.
3. Mark supplementary or terminal-soft-clipped reads (clip >=250 bp) as fusion
   candidates and exclude them from ordinary interval seeds.
4. Convert each primary non-fusion alignment to `(chrom, start, end, strand)`.
5. Cluster overlapping/nearby spans independently by chromosome and strand.
   `max_gap=0` means overlap or adjacency is required.
6. Build transcript spans from optional GTF and merge them with read intervals.
7. Return deterministic, non-overlapping strand-separated intervals plus the
   fusion read IDs/records.

Intervals carry no gene ID. Overlapping same-strand genes may share an interval.

## 5. Read collection inside an interval

Entry: `fin/candidates/discovery.py:discover_candidates`.

The BAM is fetched again for the interval:

1. Keep mapped, non-secondary alignments on the interval strand.
2. Add each retained query ID to `all_read_ids`.
3. Capture sequence only from non-supplementary alignments.
4. Classify supplementary or >=250 bp terminal-soft-clipped records as fusion
   records and omit them from structural generation.
5. Pass the remaining records to candidate generation.

Secondary CIGAR chains cannot contribute structural support, so one query cannot
populate several variants through primary/secondary mappings.

## 6. Intron-chain extraction

`fin/candidates/intron_chains.py:extract_intron_chain` walks CIGAR operations
from `reference_start`. Each `N` operation creates one genomic intron. The
production family path uses those CIGAR coordinates directly:

- No 3-prime clustering is used to define a multi-exon family.
- No canonical-site search/snapping is used in the family branch.
- Canonical motif validation occurs later as a drop gate.

The legacy generation path can expand each junction over nearby canonical sites
and cap the Cartesian product with `max_chains_per_read=16`; production family
generation bypasses that path.

## 7. Family construction

Entry: `fin/candidates/chain_cluster.py:cluster_families`.

### 7.1 Exact variants

Reads are grouped by their exact intron tuple. Chain-less reads enter one mono
holding bucket. For each exact multi-exon chain, `variant_reads[chain]` records
its discovery read IDs.

### 7.2 Structural relation graph

Distinct chains become vertices. An undirected edge exists if any relation is
true (`chain_cluster.py:123-197`):

- `wobble(a,b)`: same intron count and, for every pair,
  `abs(donor_a-donor_b) <= 6` and `abs(acceptor_a-acceptor_b) <= 6`.
- `containment(short,long)`: `short` is a strict contiguous subchain of `long`,
  with both coordinates of every aligned intron within 6 bp.
- `cassette(a,b)`: one chain has exactly one extra intron, the outer splice
  boundaries match within 6 bp, all other introns match, and the skipped exon is
  `0 < length < 70` bp.

Union-find connected components define families. This is single-linkage: A can
join C through B even when A and C are not directly related.

The function itself only groups. It does not fold, choose final read ownership,
or synthesize a new path. The discovery adapter attaches each multi-exon GTF to
at most one related read family without bridging families; unmatched GTF chains
form GTF-only families. Every family receives a stable BLAKE2b ID over chromosome,
strand, and its sorted read/GTF intron chains.

## 8. Exact-subchain collapse

Entry: `fin/candidates/chain_cluster.py:collapse`.

Within each family, variants are sorted longest-first, then by support and
coordinates. A strict exact contiguous subchain folds into the first retained
longest container.

With the default span guard:

1. Identify introns present in the container but outside the subchain match.
2. For each subchain read with span `[s,e)`, test whether it covers an extra
   container intron completely: `s <= donor` and `e >= acceptor`.
3. A spanning read is positive evidence for a retained-intron/alternative
   structure, so it remains on a surviving subchain candidate.
4. A non-spanning or span-missing read is treated as degradation and folds into
   the longer container.
5. Folded chains and read IDs remain dormant provenance under the container.

Wobble, cassette, and non-exact containment siblings are not folded. No new
chain is synthesized. The read sets of folded variants are pooled into the
container before assignment.

This operation ignores transcript endpoints as independent isoform states.
Folded 5-prime alternatives can only be recovered by the cluster quantification
mode; default `m2_em` does not invoke the TSS recovery code.

## 9. Candidate emission

Entry: `fin/candidates/discovery.py:_chain_cluster_candidates`.

### 9.1 GTF candidates

All same-strand GTF transcripts overlapping the interval are emitted with their
annotation span and spliced sequence. Multi-exon GTF hypotheses inherit their
attached or GTF-only family ID. `gtf_by_chain.setdefault` chooses the first GTF
transcript for merging discovery reads when several share an exact chain; the
other GTF candidates remain available.

### 9.2 Multi-exon novel candidates

For every surviving chain member:

- Candidate start is the minimum start of its member reads.
- Candidate end is the maximum end of its member reads.
- Exons are the outer span minus the introns.
- Sequence is stitched from the reference and reverse-complemented on `-`.
- An exact-chain GTF match receives the reads instead of creating a novel.
- Otherwise a stable structure-derived `novel_<blake2b>` ID is assigned.
- The candidate inherits its discovery family's stable `family_id`.

All reads with the same exact intron chain create one novel model; alternative
TSS/TES/APA states with the same chain are not separately represented.

### 9.3 Mono-exon candidates

The production runner disables generation-time mono folding when
`mono_resolve_post_em=True`. Chain-less reads are spatially split into
non-overlapping loci. A locus is matched to the best-overlap mono GTF transcript
or emitted as one novel mono candidate spanning the union of its reads.

### 9.4 Minimum discovery support

If `min_novel_reads > 1`, low-support novel candidates are removed after
collapse. The current default is 1, but later junction and abundance gates are
stricter.

## 10. Optional de-novo generation paths

Off by default:

- `denovo_graph`: clusters observed junctions, builds an intron graph, and
  enumerates read-supported paths with TSS braking.
- Legacy canonical alternative expansion: searches nearby canonical motifs and
  emits bounded per-read chain alternatives.
- `explore_family_paths`: family-local path completion. It has unit tests but no
  production call site.
- `cluster_diff_regions`: emits fine comparable allele regions with shared
  genomic flanks. It has unit tests and experiments but no production call site.
- `candidate_align`: builds one multi-candidate mappy index and calculates
  coverage/identity/3-prime goodness. It is not production-wired.

## 11. Fusion generation

Fusion is off by default. The enabled path is:

### F1: chimeric read resolution

`fin/fusion/chimeric.py` takes primary alignments with long terminal clips,
re-aligns the longest clipped segment against a genome-wide mappy index, and
constructs two genomic arms. Reads with an internal unmapped gap >=30 bp are
rejected as likely adapter chimeras. Missing mappy/backend state fails closed to
no fusion call rather than aborting the pipeline.

### F2: arm assembly

`fin/fusion/arm_assembly.py` single-linkage clusters reads with equal
`(chromA,strandA,chromB,strandB)` and both breakpoints within 500 bp. Each arm
receives read-derived chains plus optional overlapping annotation variants.
Read-derived variants require at least 2 reads; annotation variants are exempt.

### F3: stitching

`fin/fusion/stitch.py` forms arm-A x arm-B combinations. Read/read pairs require
an intersection of at least `fusion_min_support` read IDs. Combinations involving
annotation are support-exempt and inherit support from the read arm or cluster.
The result is a fusion `TranscriptCandidate` with joined sequence and two
breakpoints.

Fusion candidates are appended to the ordinary candidate set and compete in the
selected quantification mode. They are exempt from most later drop filters.

## 12. Pre-assignment gates

### 12.1 Canonical gate

`fin/pipeline/selection.py:canonical_gate_select`, default on:

- Applies only to novel multi-exon candidates.
- Keeps a candidate only when every intron has one of GT-AG, GC-AG, AT-AC in
  transcript orientation.
- GTF, fusion, and mono candidates are exempt.
- Missing chromosome sequence makes the gate a no-op.

### 12.2 Junction dominance

`junction_dominance_select`, default off:

For every novel intron, count same-strand primary-read CIGAR junctions within
`junction_dominance_tol_bp=2`. Drop the candidate if its junction support is
below 2 or a distinct junction within a 20 bp window has strictly greater
support. GTF/fusion/mono candidates are exempt.

## 13. M1 sequence score

`fin/scoring/mappy_score.py:score_hit` reconstructs a map-ont alignment score:

```text
AS = 2 * matches - 4 * mismatches
     - sum_over_indels min(4 + 2*k, 24 + k)

mismatches = max(0, hit.NM - total_indel_bases)
```

Any single insertion/deletion longer than `M1_MAX_INDEL_BP` (environment
default 50) invalidates the hit. Per candidate, the best valid hit is retained.

`fin/scoring/mappy_distance.py` transforms scores for `m1_em`:

```text
d_M1[i,j] = row_max_AS[i] - AS[i,j]
```

A missing candidate hit receives `row_max - (row_min - 1)`, the worst distance.
A read with no candidate hit has an all-zero uninformative row.

The production `m2_em` path does not use this dense distance. It uses M1 only to:

- Keep reads with at least one valid `AS > 0` and above hidden
  `MAPPY_R1_MIN_AS` if set.
- Define an exact best-AS tie set:
  `T_i = {j: AS_ij >= max_AS_i - 1e-9}`.
- Mask every candidate outside `T_i`.

Because reconstructed AS is effectively discrete, this is an exact tie, not a
near-tie margin.

## 14. M2 junction signal score

Entry: `fin/pipeline/assignment.py:tie_nll` and
`fin/scoring/m2_junction_nll.py`.

Only reads with at least two exact M1 winners are eligible. `m2_metric` selects
`off`, legacy `mean`, or tight `summed_llr`; it is independent of the argmax-only
`m2_tiebreak` switch. A valid signal contrast requires at least two successfully
scored hypotheses. Otherwise M2 abstains and preserves the full M1 tie.

### 14.1 Legacy/guided-SIRV mean class window

For each tie class:

1. Locate internal structural disagreement regions.
2. For each candidate and disagreement anchor, find its nearest transcript-frame
   junction.
3. Take `K=10` transcript bases on each side.
4. Project those bases back to genomic coordinates.
5. Union all candidate positions into one genomic `gset`.

### 14.2 Mean event score

For each read-candidate tie cell:

1. Reuse the best mappy hit to slice the candidate sequence.
2. Run krill eventalign against that candidate hypothesis.
3. Project each candidate-frame event position back to the genome.
4. Keep events whose genomic position belongs to the shared class `gset`.
5. Compute candidate-specific per-event mean NLL:

```text
nll_event      = 0.5 * z^2 + log(model_stdv)
M2_mean(i,j)   = mean(nll_event for events projected into class gset)
```

When `m2_metric="mean"`, `tie_nll` uses the legacy `reduce="mean"`. The same
genomic set is requested for all hypotheses, but each hypothesis has its own
event alignment and may contribute a different event set/count. Before batched
eventalign, every candidate slice must contain only A/C/G/T/U. A read with an
ambiguous candidate slice abstains individually and retains its M1 tie; it can no
longer poison the entire locus and force a per-read fallback.

### 14.3 Coverage gate and sparse distance

With `m2_diff_cover_gate=True` (both named profiles):

- Build structural wobble spans for the tie.
- A span is called covered when at least one candidate projection straddles it.
- Require at least two successful candidate scores; one/zero means abstain.
- Sort successful scores and compute `margin = second_best - best`.
- Use the metric-specific margin: mean 0.5, or real-profile summed LLR 1.0.
- If the margin passes, set winner distance to 0 and every other tie cell to
  `MISSING=1e6` (hard assignment after softmax).
- Only covered, distinguishing reads add one vote to the winner's
  `covered_vote`.
- Reads with no NLL or margin <0.5 are fuzzy.
- Fuzzy reads split according to covered votes among their tie candidates:
  `d_ij = -log(vote_j / sum_votes)`.
- With no votes, fuzzy reads use distance 0 for every tied candidate, yielding
  a uniform `1/K` split.

With the coverage gate off:

- Singleton ties get distance 0.
- Unscorable ties are flat.
- Scored cells use `NLL_j - min_NLL`.
- Unscored cells use `(max_NLL-min_NLL)+1`.
- Positive distances are divided by their median before EM.

### 14.4 Tight-window summed LLR

`m2_metric="summed_llr"` is wired into default `m2_em` assignment and is the
real-dRNA profile default:

1. Require tied candidates to have the same intron count. Containment/cassette
   contrasts abstain because their candidate-private event populations are
   asymmetric.
2. `diff_junction_windows` keeps tight windows around differing junction
   boundaries (`flank=8` in the real profile).
3. Candidate event NLL uses `reduce="sum"`, not a per-event mean.
4. The margin is the undivided summed-NLL gap between the two best successfully
   scored candidates (`>=1` in the optimized real profile).

The design reason from the original experiment remains: mean division dilutes a
few discriminating events among many shared events. Historical code had this
metric only in cluster mode; the profile integration added equivalent batch and
fallback reducers to `fin.pipeline.assignment.tie_nll`.

Real p00 validation selected margin 1/flank 8 over eight other combinations.
Against wide mean, tight sum improved T1/T3 F1 slightly and reduced assignment
runtime from 1391 to about 599 seconds before other profile filters. A later
six-sample SIRV paired test found a guide interaction: true p00 SUM gained 0.722
mean F1 (4 wins, 1 tie, 1 loss), while full-guide SUM lost 0.284 (2 wins, 4
losses). The SIRV profile therefore resolves `auto` to SUM (margin 1/flank 4)
when no usable GTF exists and mean when a GTF has at least two transcripts.

Summation is not promoted beyond its validated niche. Same-intron-count scope,
two-score abstention, and distinct sum thresholds prevent technical missingness
or containment event-count asymmetry from becoming false biological certainty.

## 15. Retired coherence prototypes and PolyA

### 15.1 M3 junction coherence [removed]

M3 read-by-read junction-window DTW was removed from production configuration,
CLI, assignment, and ablation rows. Current named profiles never enabled it, and
its missing-pair semantics were not production-ready. The implementation and
focused tests remain under `experiments/m3_coherence/` for possible redesign;
the generic EM engine remains available, while production M1/M2 calls it with
`beta=0`. Post-removal real balanced/precision and SIRV p00/full GTF/TSV outputs
are byte-identical to the pre-removal runs.

### 15.2 M4 difference-region DTW [experimental only]

`fin/scoring/diff_region_dtw.py` can construct class-partitioned read-read
coherence from structural difference regions. It has no production config or
assignment wiring after M3 removal.

### 15.3 PolyA

`fin/scoring/polya.py` performs whole-read krill eventalign with `polya=True`
and returns `(polya_length, polya_qc)`. Batch failure disables the entire polyA
filter to avoid interpreting a technical failure as negative evidence. The
configured pore model defaults to RNA002.

## 16. Quantification modes

### 16.1 `argmax`

The dispatcher uses `_quant_argmax_keep`, not the older `_quant_argmax_first`.
Every read is aligned to every candidate. Exact best-AS ties split mass `1/K`;
non-tied cells receive zero. There is no iterative EM.

The unused `_quant_argmax_first` hard-selects the lowest candidate index on an
AS tie. Since discovery orders GTF before novel candidates, that method contains
an implicit GTF prior, but it is not the current dispatcher path.

### 16.2 `m1_em`

Build the dense M1 distance matrix and call `em_with_coherence` with `beta=0`,
distance temperature `em_sigma`, and optional abundance feedback/length
normalization. M2 and read-by-read coherence are absent. Optional krill tiebreak can rewrite
ambiguous responsibilities afterward.

### 16.3 `m2_em` (default)

Build the sparse exact-tie M2 distance described above. Call
`em_with_coherence` with `sigma=1.0` and `beta=0`, not `config.em_sigma`. Then quantify via
`fin/analysis/quantification.py:quantify_transcripts` and run interval-level
selection. `m2_metric` controls this path directly: `off` keeps the exact M1
tie, `mean` uses the wide class window, and `summed_llr` uses the tight same-
intron-count window. The separate `m2_tiebreak` flag controls argmax mode only.

With `abundance_feedback=False`, this mode is mathematically one row-softmax
rather than an iterative abundance estimator.

### 16.4 `cluster`

Each generation cluster is processed independently:

1. Align cluster reads only to cluster members.
2. Treat all members within `cluster_m1_tie_margin=20` of the best AS as tied.
3. A unique best becomes a certain compatibility anchor.
4. A tied read that spans a differing junction can use summed-LLR M2; a winner
   requires at least two scored hypotheses and a strict score difference.
5. Otherwise the read remains compatible with every tied member.
6. An internal compatibility EM repeatedly splits ambiguous reads in proportion
   to current member abundance while certain reads stay fixed.
7. Keep members with abundance >=1.
8. Optionally re-test folded subchain TSS proposals against weighted 5-prime read
   end peaks and recover accepted short isoforms.

Repository experiments report this mode lost to `m2_em`; it is off by default.

## 17. EM/fixed-point equations

Entry: `fin/analysis/assignments.py:em_with_coherence`.

Inputs are read-candidate distance `D`, read-read distance `Q`, temperature
`sigma`, and coherence weight `beta`.

Initialization:

```text
R_ij <- exp(-(D_ij - min_j D_ij) / sigma)
R_i  <- R_i / sum_j R_ij
```

An optional `prior_weights` vector multiplies initialization only. Production
assembly does not pass it.

At each iteration:

```text
cluster_mass_j = max(sum_k R_kj, 1e-10)
C_ij           = sum_k Q_ik * R_kj / cluster_mass_j
E_ij           = D_ij + beta * C_ij
```

If abundance feedback is enabled:

```text
count_j = cluster_mass_j
if length normalization: count_j <- count_j / effective_length_j
theta_j = count_j / sum_l count_l
E_ij   <- E_ij - sigma * log(theta_j + 1e-12)
```

Responsibility update:

```text
R_ij <- exp(-(E_ij - min_j E_ij) / sigma)
R_i  <- R_i / sum_j R_ij
```

Stop when `max(abs(R_old-R_new)) < 1e-4` or after 1000 iterations. The reported
trace is the surrogate

```text
L = -sum_ij R_ij * E_ij + sigma * H(R)
```

not a likelihood derived from an explicitly specified generative model.

Default `m2_em` has `Q=0`, `beta=0`, no feedback, and `sigma=1`; initialization
already equals every subsequent update. A direct execution check produced one
iteration and exact equality to a single softmax.

## 18. Interval-level selection

Entry: `fin/pipeline/selection.py:select_m2_interval`. It runs only after default
`m2_em` assignment. The cascade uses the original responsibility matrix and
candidate list.

1. `m2_cluster_recheck` (default on): cluster structurally similar candidates
   with 20 bp wobble and <70 bp cassette differences. Pick the highest soft-mass
   anchor. Drop a novel sibling when `abundance/anchor < 0.15`. A low GTF sibling
   is droppable only when enabled and its distinguishing junction has fewer than
   1 exact supporting read.
2. `novel_junction_min_reads=2`: drop a novel multi-exon candidate when any
   junction has fewer than 2 same-strand primary CIGAR observations within 2 bp.
3. Optional guided junction gate: analogous rule for GTF candidates; default off.
4. `m2_support_gate` (default on): a multi-exon GTF/novel candidate survives only
   if it is at least one read's sole M1 winner or an M2 minimum. Fusion and mono
   are exempt. With tie acceptance, equal M2 minima can support all tied winners.
5. `containment_cluster` (default on): drop a novel chain that is a contiguous
   <=6 bp subchain of a longer candidate when both
   `soft_ab_shadow <= 0.3 * soft_ab_parent` and
   `discovery_reads_shadow <= 0.3 * discovery_reads_parent`, provided the shadow
   has at most 10 discovery reads.
6. Optional `containment_collapse` (default off): fold a strand-consistent exact
   3-prime suffix/prefix novel truncation into a higher-abundance parent. This
   operation transfers soft abundance and hard read IDs.
7. `mono_resolve_post_em` (default on): for each mono hard-assigned read, find
   surviving multi-exon candidates containing its span inside one exon. Fold to
   the sole host or highest-EM-abundance host. The mono retains uncovered reads
   and is dropped if fewer than 2 remain.
8. Remove every column marked by a non-fold drop. Interval selection does not
   re-estimate those columns; named profiles carry the original responsibility
   ledger to the final survivor-set refit in Section 20.

Observed junction evidence skips secondary/supplementary/unmapped records and is
memoized per interval. Missing/unreadable evidence disables support gates.

## 19. Global selection

After commutative interval aggregation, `select_global` applies filters in this
order:

1. Source-aware abundance floor:
   - SIRV novel uses inclusive >=3.0;
   - balanced/precision real novel uses strict >1.0;
   - GTF >=1.0 unless `floor_gtf_abundance` raises it;
   - fusion always exempt.
2. Isoform fraction: for every novel multi-exon candidate, find the maximum
   abundance among candidates in its persisted splice family and drop when
   `abundance/family_max < 0.01`. Attached GTF can set the maximum but is exempt;
   mono/fusion are family-less. A family-less legacy candidate falls back to
   same-strand genomic overlap, and `--isoform-fraction-locus overlap` restores
   the historical denominator globally.
3. Soft/hard ratio: SIRV drops a novel multi-exon candidate with zero hard
   reads or `soft_abundance / hard_read_count >= 2.0`; real profiles disable it.
4. Novel mono gate: balanced and precision real profiles require >=5 hard reads;
   SIRV/custom leave it off unless requested.
5. Full-length fraction: for novel candidates with at least 3 exons and at least
   4 scored hard reads, count reads whose strand-aware 5-prime and 3-prime ends
   are both within 25 bp of candidate ends. Drop when the fraction is <0.1.
6. PolyA/5-prime support: when signal is available, require at least one hard
   read with krill polyA QC PASS, tail length >10, and genomic 5-prime end within
   25 bp of the candidate. Fusion and, by default, GTF are exempt. Empty/failed
   polyA scoring disables the filter.

These are sequential: later locus-relative filters see only earlier survivors.
The survivor set is intentionally frozen on these pre-refit values; refitted
abundance is not fed back through the existence filters.

### 19.1 Finalized-model junction consensus

Both real profiles then scan primary BAM CIGARs once and count exact junctions by
chromosome/strand. For each novel multi-exon survivor, an intron may move <=6 bp
to the strongest nearby mode only when target support is >=2 and strictly more
than twice current-coordinate support. With the canonical gate on, the target
must pass the same configured motif set. GTF, fusion, and mono models are exempt.
Models that become identical in complete exon coordinates merge; soft abundance
is summed, hard read IDs are unioned, confidence is hard-read weighted, and
`max_R` takes the maximum. The merge also records absorbed-to-representative IDs
so the later responsibility ledger follows the same structural identity. This is
annotation-free and runs after existence selection, so it cannot change the EM
competition set. It gained T3 F1 on three real samples while retaining or
increasing recall. SIRV leaves it off.

## 19b. Candidate evidence layer and the frozen evidence ranker (optional)

After global selection and before junction snapping/refit, the pipeline can
compute one observable-feature row per survivor
(`fin/analysis/candidate_evidence.py`): weakest/median BAM junction support,
family share/rank from persisted discovery families, sibling sub/superchain
geometry, strand-aware 5'/3'/full-length read-end agreement, per-junction
canonical fraction, mono/exon-count/length, and EM quantities. One whole-BAM
pass (`collect_ranking_bam_evidence`) supplies junction support and primary
read ends together and carries an explicit completeness flag. Sentinel -1
means "evidence source unavailable", never "low". `--candidate-evidence`
writes `candidate_evidence.tsv` (read-only diagnostic; outputs unchanged).

`--ranking-mode filter` (EXPERIMENTAL, off by default and enabled by no
profile) scores each row with a frozen L2-logistic model
(`fin/analysis/candidate_ranking.py`, exact constants mirrored bit-for-bit in
`experiments/prod_validation/models/candidate_ranker_v1.json`) and removes
NOVEL candidates whose raw logit falls below the frozen operating point;
GTF/fusion are always exempt, an incomplete BAM scan refuses to filter, and
an explicitly requested filter that cannot compute evidence raises instead of
silently keeping everything. Filtering runs BEFORE snapping and the
survivor-abundance refit, so released read mass follows the ordinary refit
accounting (renormalized or explicitly orphaned). The score is a ranking
logit, not a calibrated probability; the `confidence` column is untouched.

Provenance: trained on the H9 r2r2 tuning sample only (labels = gffcompare
exact-match vs GENCODE, used only as labels); chromosome-grouped CV AUC
0.809; threshold frozen on the tuning frontier under a hard T1-not-lower
constraint. Frozen evaluation: r3r1 T1 38.044->38.600 / T3 39.732->42.077
(both above StringTie3/IsoQuant); one-shot r4r2 audit T1 38.580->38.256 /
T3 39.788->40.758 (T3 gap to IsoQuant halved but not closed, own T1
slightly lower) - which is why the mode ships experimental rather than
default, pending the reserved unopened holdouts.

## 19c. EndpointRefine (EXPERIMENTAL, off by default)

`--endpoint-refine` splits a novel multi-exon survivor into at most
`--endpoint-max-splits` endpoint states when strand-aware read-end modes
support distinct (TSS, TES) pairs: 25 bp mode clustering of assigned reads'
genomic ends, pair support >= 3 reads and >= 15% of end-mapped reads,
interior-TSS (5'-degradation-direction) modes requiring 2x support, stable
BLAKE2b endpoint IDs, and mono/GTF/fusion exempt. Post-split
requantification is mandatory and structural: the split plan emits per-read
routes that `refit_survivor_abundance` consumes (`split_routes` /
`split_primary`), so every read's mass is re-dealt over the endpoint states
under the same conservation invariants, and `validate()` rejects
`endpoint_refine` without the effective refit. Poly(A)-supported TES
strengthening is a declared v2 hook (`polya_read_ids`); v1 is signal-free.
Measured upper bound of the addressed error class: 62/4,477 (1.38%) of
missed T3 multi-exon truths. Off in every profile pending holdout data.

## 19d. TSS evidence for contained models (EXPERIMENTAL, off by default)

`--tss-evidence-mode audit|require` (needs `--endpoint-refine`) decides whether
a CONTAINED shorter model is a real transcript or a 5'-degradation artifact of
the longer one. dRNA reads 3'->5', so the 3' end is reliable and the 5' end is
degraded; when a short model shares the parent's chain AND TES and differs only
by an internal TSS, the two explanations are mathematically unidentifiable, so
the verdict is three-way: `supported` / `unsupported` / `unidentifiable`.

Evidence is routed by identifiability rung. A state with its OWN 3' end is
decided by `evaluate_tes_support` on the reliable end (reads terminating at the
candidate's 3' cluster, which degradation cannot fake) after restricting the
read pool to reads ending there. A state sharing the parent's 3' end has no
internal TSS signal available except the conditional termination hazard
`ends_at(d)/at_risk(d)` in spliced coordinates, tested against a LOCAL
neighbourhood-median null (a single global hazard mislabelled degradation
hotspots at a 52% rate).

Three statistical guards, all of which change the answer. (i) SELECTION: the
candidate bin was chosen because it looked extreme and is then tested on the
same reads, so the pointwise bootstrap p is Bonferroni-adjusted by the size of
the search space actually scanned -- in the pipeline the FIRST EXON's length
in bins, which is deterministic and is exactly the region EndpointRefine may
propose a start in. A simulated max-statistic was implemented first and
rejected: it needs the risk set at every candidate bin, and reusing the
candidate's own `n_at_risk` understates the null spread. (ii) The bootstrap is
empirical, not asymptotic chi-square, because the mixture weight is a boundary
parameter. (iii) MULTIPLICITY, two levels: each alternative gets
`p_within = min(1, p * k)` for the k alternatives on its parent model, the
parent contributes the minimum of those upward, and Benjamini-Yekutieli runs
across parent models -- BY, not BH, because alternatives share reads, coverage
and the degradation background so positive dependence is not guaranteed. A
candidate survives only if BOTH its own `p_within` and its parent q clear
alpha, so one strong start cannot certify its weak siblings.

This is parent-MODEL FDR, not gene-locus FDR: de novo assembly has no gene
annotation, and the grouping key is the parent candidate id. Structural decisions and abstentions carry
`p_value = None` (JSON null, never NaN) and never enter the family, so FDR is
claimed for every calibrated HAZARD test (including `own_tes` candidates with
a distinct TSS) and never for the structural TES-only or parent-unobserved
rules. Peak sharpness (`peak_mad`) rejects broad humps,
including on the no-upstream-reads path. `unidentifiable` never drops a model.

Validation, reported per rung and per locus; no pooled AUC is computed because
the two rungs are scored by rules that are not calibrated to a common scale.
SIRV, deduplicated per short transcript, after the full correction: all 3
positive loci `supported`, the fourth positive transcript abstaining (retained
but unverified), and ZERO false positives among the 85 negative rows
(84 `unsupported`, 1 abstention). Ground-truth simulation of the hardest rung:
zero false positives among all 151 decided negatives (the FPR is undefined at
depth 10, where every negative abstains); recall 0.097 / 0.741 / 1.000 at
depth 20 for 10% / 25% / 50% short fractions, 0.625 / 1.000 / 1.000 at depth
50, and abstention below 20 reads. Live human sample: `audit` byte-identical
to `--endpoint-refine` alone, mass error 1.46e-10, repeat runs byte-identical.
Full study: `TSS_CONTAINMENT_STUDY.md`.

## 20. Abundance, confidence, aggregation, and TPM

For responsibility matrix `R` and hard labels `h`:

```text
abundance_j        = sum_i R_ij
hard_count_j       = count_i(h_i == j)
confidence_j       = mean_i(R_ij where h_i == j), or 0
assigned_read_ids  = {read_i: h_i == j}
```

Across intervals, `fin/analysis/quantification.py:aggregate_across_intervals`
unions IDs. Abundance is deduplicated heuristically:

```text
abundance_out = sum_interval_abundance
                * unique_hard_read_count / abundance_weight
abundance_weight = sum_interval max(hard_count, 1)
```

When no hard IDs exist, aggregation falls back to the naive abundance sum.
Confidence and diagnostic scores are hard-read-weighted averages. `max_R` and
full-length fraction use maxima across intervals.

After all filtering, `fin/analysis/quantification.py:compute_tpm` calculates
TSV TPM as:

```text
RPK_j = abundance_j / (spliced_length_j / 1000)
TPM_j = RPK_j / sum_k RPK_k * 1e6
```

Named profiles first order interval `(QuantResult, ledger)` pairs by a
strand-aware genomic key before legacy aggregation, then run
`fin/analysis/abundance_refit.py` after global selection and junction snapping.
Per interval, the final beta=0 responsibility matrix is
stored sparsely as `read -> {candidate: weight}`. For each unique read, interval
weights are summed, containment and snap merges redirect candidate IDs, and the
surviving weights are renormalized:

```text
S_i = sum_{j in final survivors} w_ij
R'_ij = w_ij / S_i                         when S_i > 0
unassigned_i = 1                           when S_i == 0
```

For a read carried by a single interval this is exactly the beta=0 softmax
result after deleting non-survivor columns. A read whose responsibilities come
from several intervals is not a single softmax row: its per-interval rows are
summed first, which weights each interval's opinion equally, and only then
renormalized. Instrumented measurement on real balanced tuning found 11 of
294,111 assignable reads (0.0037%) in more than one ledger, bounding the
affected mass at 0.0045% of assigned mass.
Mono-resolved fragment reads instead retain their explicit mass-1 parent mapping;
if different intervals force different parents, the highest pre-refit final
abundance wins with candidate ID as deterministic tie-break. The pass recomputes
abundance, hard IDs/counts, confidence, and `max_R`. Forced mono reads contribute
to abundance/hard count but not the EM-confidence mean. `abundance_refit.json`
separates reads that never received a responsibility row
(`interval_quantification_unassigned_reads`, aliased by the legacy
`alignment_unassigned_reads`; both are scoped to intervals that reached
quantification and are not BAM-wide unaligned counts) from selection-orphaned
reads, and asserts `assigned mass + unassigned mass = assignable reads`. The
read counters deliberately overlap: a forced read whose parent did not survive
is both forced and orphaned, so
`forced + renormalized + orphaned - forced_orphaned = assignable`. Paired on/off structural
checks compare normalized GTF record multisets because equal `(chrom,start)`
gene groups retain their upstream insertion order; independent refit-on repeats
must still be byte-identical. `custom` and explicit
`--no-post-selection-refit` retain the historical model-filtered values.

GTF gene IDs are resolved after selection. Novel candidates default to
`gene_id=candidate_id`, so novel isoforms are not grouped into a shared output
gene unless another layer assigns one.

## 21. Parallel execution

Genome access is lazy by default (`fin/io/lazy_genome.py`): each process
opens the indexed FASTA and holds at most `genome_cache_chroms` (2)
chromosomes, replacing the eager whole-genome dict that dominated worker
memory (measured 67.1 -> 40.9 GB whole-job peak on the tuning sample with
byte-identical outputs; `lazy_genome=False` restores the historical eager
load).

`fin/pipeline/parallel.py:run_parallel` uses multiprocessing `spawn`:

- The parent generates a deterministic interval list.
- Each worker rebuilds that list and receives interval indices.
- Exactly `gpu_workers` processes request GPU; the rest are CPU.
- One `PipelineRunner` is created and reused per worker.
- Results are collected by interval index and aggregated commutatively.
- Spawn avoids inherited BLAS pools and non-fork-safe CUDA/pysam/POD5 handles.

Structural/abundance determinism is tested when large fixtures are available.
Novel transcript IDs are stable BLAKE2b hashes of chromosome, strand, genomic
endpoints, and intron chain; equal-abundance tie-breaks and output identity no
longer depend on `uuid4()`.

## 22. Default paths in one page

Shared path:

```text
BAM primary (non-secondary) interval seeds
  -> strand-separated overlap intervals
  -> exact CIGAR intron chains (3' endpoints ignored)
  -> family union-find: wobble 6 / cassette <70 / containment
  -> exact-subchain collapse with span guard
  -> union-span novel candidates + same-strand GTF candidates
  -> canonical novel gate
  -> M1 align every read to every candidate
  -> keep AS>0; retain exact best-AS tie set only
  -> profile-specific M2 exact-tie refinement
  -> one-softmax assignment (read coherence absent, abundance feedback off)
  -> interval recheck/junction/support/containment gates
  -> post-EM mono read folding
  -> cross-interval aggregation and profile-specific finalization
  -> optional finalized-model junction consensus + mass-preserving merge
  -> survivor-set responsibility refit + orphaned-mass diagnostics
  -> GTF / scoring TSV / refit audit / manifest / optional BEDPE
```

`real-drna` (raw CLI default):

```text
tight summed LLR (margin 1, flank 8, same intron count only)
  -> strict abundance >1 read-equivalent
  -> novel mono hard-read floor 5
  -> isoform fraction 1%
  -> read-supported junction consensus (6 bp / support 2 / ratio 2)
  -> survivor-set abundance refit
  -> soft/hard, full-length, and polyA gates OFF

real-drna-precision adds generation support >=2; balanced real keeps 1.
```

`sirv`:

```text
SUM (margin 1/flank 4) when unguided; wide mean (margin 0.5) when guided
  -> novel/GTF abundance 3/1 with GTF floor raised to 3
  -> isoform fraction 1%, soft/hard ratio 2x, full-length 10%
  -> survivor-set abundance refit
  -> polyA gate OFF
```

The algorithm remains a structural heuristic assembler with conservative,
exact-tie signal refinement rather than a dense raw-signal generative model.
Profiles make the empirically different synthetic and biological operating
points explicit instead of hiding them in one default.

## 23. Integration status of recorded algorithm changes

`experiments/prod_validation/PRODUCTION_STATE.md:83-91` was written around an
uncommitted working state. Its rows do not all describe current HEAD equally:

| Recorded change | Current HEAD status | Classification |
| --- | --- | --- |
| SIRV vs real operating point | Named, validated profiles plus manifests | Landed in CLI and Python API |
| fusion soft-clip 50 -> 250 | Present in interval/discovery fusion classification | Landed |
| generation-time mono fold | Available, but disabled when default post-EM mono resolution is on | Superseded by a later policy |
| exact-subchain span guard | Present and default on in family discovery | Landed |
| mono split by locus/GTF overlap | Present in production family discovery | Landed |
| folded shadows | Preserved in `CandidateSet.shadows` | Landed as provenance |
| stable novel IDs | Content hash of chromosome/strand/endpoints/chain | Landed for reproducibility and deterministic ties |
| summed-LLR M2 | Wired into `m2_em`; real default; SIRV auto-selects it only when unguided | Landed and guide-selected |
| 5-prime recovery scaffold | Implemented inside cluster quantification | Off-default; winning mode abandoned it |
| `candidate_align` | Implemented and tested | Foundation, no production caller |
| `cluster_diff_regions` | Implemented and tested | Foundation, no production caller |
| `explore_family_paths` | Implemented and tested | Foundation, no production caller |
| Read-by-read signal coherence | M3 retired; M4 experimental | Prototype preserved, no production wiring |
| three de-novo precision levers plan | Most concepts later landed with changed semantics/defaults | Historical plan, not an applicable patch |

Current untracked `experiments/wobble_heya8/` scripts preserve another research
layer. `explore_diff_signal.py` explicitly uses summed NLL because the user
rejected mean dilution. Its `diffsig/codex_consult.md` separates well-posed
wobble contrasts from containment, where candidate-private event sets make a
symmetric sum invalid. Later 6,737-candidate existence experiments test broader
private-region signal features; their near-random signal AUROCs do not directly
invalidate the narrower summed wobble tie resolver, but they do argue against
promoting signal to a general existence gate.

There is no hidden integration branch or stash at this snapshot. Local branches
are `dev` and `main`; the only additional worktree is clean and detached at the
older `b01bcf61` refactor. Historical patch files are not a clean source of
missing code: `gpu_honesty.diff` is already applied, while `wobble.diff` and
`guided_gate.diff` have diverged from HEAD. The profile/summed-M2 gap was
integrated directly from the tested implementation and a new reproducible sweep.
Remaining unintegrated state is limited to optional research foundations without
winning end-to-end evidence, not a hidden Git branch. See
`PROFILE_OPTIMIZATION.md` for decisions and metrics.
