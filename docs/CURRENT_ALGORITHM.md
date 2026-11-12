# What PyFIN actually does

The shipped pipeline, stage by stage, as `--profile real-drna` executes it.

This document is about **mechanism**. It contains no benchmark results and no
performance claims; those belong to the manuscript. Where a rule exists
because a measurement said so, the rule is stated and the measurement is not.

Scope: the assembly path in full, plus the acceptance and abstention rules of
the independent fusion analysis. To *use* the tool see [`USAGE.md`](USAGE.md).

Source of truth: `fin/`. Every non-obvious claim below carries a `file:line`.
When this document and `pyfin_out/run_manifest.json` disagree about a setting,
the manifest is right -- it records what a given run resolved.

---

## 0. The production profile

`real-drna` sets 27 keys explicitly (`fin/pipeline/config.py`); everything
else is the raw default. The ones that shape the result:

| Setting | Value | What it turns on |
|---|---|---|
| `m2_metric` | `summed_llr` | signal junction scoring, summed log-likelihood-ratio |
| `m2_summed_llr_flank` / `_margin` | 8 / 1.0 | scoring window, structure-selection margin |
| `assignment_channels` | `ends` | endpoint evidence participates in assignment |
| `m2_support_gate_attest_reads` | **1** | the attestation exemption in §5.2 |
| `canonical_gate` | on | novel multi-exon candidates need canonical junctions |
| `novel_junction_min_reads` | 2 | each novel junction needs two observing reads |
| `m2_cluster_recheck` | on | post-EM wobble-shadow removal |
| `containment_cluster` | on | low-support sub-chain removal |
| `folded_tss_recovery_mode` | `require` | folded short chains may be recovered |
| `junction_snap` (+ `snap_merge_chain`) | on | novel junctions snap to stronger observations |
| `post_selection_refit` | on | abundance is re-fit over survivors |
| `min_abundance` / `min_gtf_abundance` | 1.0 / 1.0 | soft-mass floors |
| `strict_novel_abundance_floor` | on | novel floor is strict `>`, not `>=` |
| `min_isoform_fraction` | 0.01 | minor-isoform suppression for novel |
| `drop_mono_exon_novel` + `min_mono_exon_reads` | on, 5 | weak novel single-exon suppression |
| `max_soft_mass_ratio` | **0.0** | *disabled* |
| `min_fulllen_fraction` | **0.0** | *disabled* |
| `min_polya5p_reads` | **0** | *disabled* |

The last three matter: they are SIRV-derived levers that the direct-RNA
profile turns **off**. Reading their help text as shipped behaviour is a
common mistake.

---

## 1. Input and sample-level preparation

**Intervals.** Reads are grouped into working intervals by overlap, separated
by chromosome and strand (`fin/io/interval_manager.py:449`). Reads carrying a
fusion signature are held out of the interval read pool entirely
(`fin/io/interval_manager.py:641`), so they never enter candidate generation
or the quantification axis.

**Poly(A) is signal evidence; endpoints are not.** Whether a read has a
poly(A) tail is decided from the raw signal. Where its ends *are* comes from
the BAM. The two are kept apart deliberately.

**Endpoint nuisance.** The pipeline fits, per sample, how much transcript ends
move on their own -- so that a later endpoint disagreement can be judged
against that sample's own noise rather than a fixed tolerance. The fit uses
annotation distances where an annotation is available and a split-half
fallback otherwise (`fin/analysis/sample_nuisance.py:608`).

Discovery uses the pooled fit. **Assignment refits it leave-one-interval-out**
-- excluding the interval's own reads -- so an interval cannot set the
yardstick it is then measured by (`fin/pipeline/runner.py:2127`,
`fin/pipeline/runner.py:1994`).

Sample-level evidence is cached and shared across runs of the same sample.
What is cached is the common evidence, never a verdict about a candidate.

---

## 2. Discovery: generating candidates inside one interval

1. **CIGAR to intron chain.** Each read's spliced alignment becomes an ordered
   intron chain.
2. **Exact grouping first.** Reads are grouped on the chain alone. Endpoints
   do not participate in identity at this stage -- two reads of the same
   isoform truncated to different extents are the same structure.
3. **Families.** Chains that need to be compared with each other are grouped
   into a family. A family is a comparison scope, **not** a consensus: member
   coordinates are never rewritten to a family consensus.
4. **Collapse.** A chain that is an exact contiguous sub-chain of a longer one
   is folded into it as a *shadow*, longest-first, then by read support and
   coordinates (`fin/candidates/chain_cluster.py:748`). Two qualifications,
   both active by default:
   - `chain_cluster_fold_span_guard` is on. A read that aligns **exonically
     across an extra intron of the proposed container** contradicts the fold;
     the shorter chain survives, carrying those reads.
   - certain annotation-only chains cannot act as containers
     (`fin/candidates/chain_cluster.py:745`), and an annotated chain is exempt
     from being folded.
5. **Emit.** Candidates get sequences and endpoint hypotheses, including
   `.ext` siblings. An `.ext` sibling is a **read-supported extent that
   contradicts the annotated one** -- its coordinates come from `_read_extent`
   (`fin/candidates/discovery.py:760`, `:787`), not from the annotation. Its
   `source="gtf"` does not make its new endpoints annotated.

**Compatibility is computed here**, because the profile asks for it
(`fin/pipeline/runner.py:754`, `fin/candidates/discovery.py:591`). §5.2 uses it.

---

## 3. Pre-assignment gates

These run before any read is assigned, and they remove *candidates* without
removing those candidates' reads from the competition pool. The novel-junction
read-support requirement is **not** one of them: it runs after assignment, in
interval selection (§5.2).

- **Canonical gate** (on): a novel multi-exon candidate with any non-canonical
  junction is dropped. GTF, fusion and mono candidates are exempt
  (`fin/pipeline/selection.py:1055`). A candidate removed here cannot be
  repaired by the later junction snap -- it no longer exists.
- **Junction dominance** is an additional conditional branch
  (`fin/pipeline/runner.py:815`) and is **off** by default
  (`fin/pipeline/config.py:990`).

**The junction registry records evidence; it is not a third gate.** No gate
re-reads the serialised `junction_support.tsv`. But the in-memory
`CandidateSet.junction_support` is not passive either -- it feeds the
attestation veto in §5.2 (`fin/pipeline/selection.py:369`).

---

## 4. Assignment: first the competition set, then the evidence

### 4.1 Structural mask (M1)

Every read is re-mapped to every candidate cDNA **in the interval**, producing
a read x candidate alignment-score matrix. Families govern *discovery*; they
are not the boundary of assignment competition -- the runner hands the whole
interval candidate set to `Assigner.assign` and the matrix is built over all
of it (`fin/pipeline/runner.py:2531`, `fin/pipeline/assignment.py:959`,
`:986`). Exact AS ties are what the signal
layer exists to break; a read with a single best candidate is already decided.

### 4.2 Junction scoring (M2)

Two windows, easily confused:

- the **summed-score window flank is 8** (`m2_summed_llr_flank`);
- the **coverage-span flank is 2** (`_M2_EM_FLANK`,
  `fin/pipeline/assignment.py:40`).

*Coverage* means every disagreement span is straddled on at least one
candidate alignment. A sufficient structure-level score margin can select a
structure **without** verified coverage; verified coverage is additionally
required for a candidate to contribute a `covered_vote`
(`fin/pipeline/assignment.py:458`, `:577`, `:1122`).

### 4.3 Endpoint channels

With `assignment_channels = ends`, poly(A) 3' evidence adds an NLL term only
on cells that are still eligible ties, and the 5' channel is computed per read
without requiring a start cluster.

### 4.4 Turning costs into responsibilities

A softmax over the combined costs gives each read a distribution over
candidates. **This is not an abundance EM** -- abundance comes later.

The no-vote fallback has a specific shape worth stating, because it is easy to
get backwards (`fin/pipeline/assignment.py:1150`):

- **With covered votes:** normalise over *structures*, then divide each
  structure's share among its tied candidates.
- **Without any covered vote:** set **equal costs on every tied candidate**.
  A structure represented by two endpoint siblings therefore receives twice
  the initial mass of a singleton structure, before endpoint terms apply.

### 4.5 Signal abstention vs. dependency failure

These are different and the distinction matters when reading a run:

| Situation | Behaviour |
|---|---|
| krill not installed | **raises** at the unguarded import (`fin/pipeline/assignment.py:365`) |
| krill imports, backend yields no aligner | AS ties preserved, no junction scores (`:396`) |
| invalid sequence payload, missing result, too few scored hypotheses | abstain (`:487`, `:529`, `:597`) |
| batch failure | per-read fallback |
| poly(A) unavailable | its own path (`fin/scoring/polya.py:44`) |

The CLI requiring `--signal` (`fin/cli.py:456`) does not mean every read has
usable signal evidence, and a signal failure does not generally abort the run.

---

## 5. After assignment, still inside the interval

### 5.1 Folded-TSS shadow recovery

A chain folded away in §2.4 may be recovered -- but only after asking whether
the short chain has a start of its own, rather than being a truncation of the
container.

### 5.2 Interval selection, and the support gate

Extent is resolved first, then the wobble-shadow recheck, then novel junction
support, then the support gate, then containment, then mono
(`fin/pipeline/selection.py:257`, `:276`, `:416`, `:448`, `:529`).

**Novel junction support** (2 reads): a novel multi-exon candidate is dropped
if any of its junctions is observed by fewer than two same-strand primary
CIGAR introns within 2 bp on both boundaries (`fin/pipeline/selection.py:266`,
`fin/scoring/m2_junction_nll.py:1096`). GTF, fusion and mono candidates are
exempt. It is evaluated here, after assignment, not with the §3 gates.

**The ordinary support rule.** A multi-exon candidate (GTF or novel) survives
the support gate if it is some read's **sole** M1 best, or some read's **M2
minimum** within a scored tie. Fusion and mono candidates are exempt. With tie
acceptance, equal M2 minima support all tied winners.

The gate also **fails open**: when both support-evidence dictionaries are
empty it drops nothing (`fin/scoring/m2_junction_nll.py:578`), so "survives
only if" would be too absolute.

**The shipped default adds a conditional exemption on top of that rule, and
the rule alone no longer describes production.** `real-drna` sets
`m2_support_gate_attest_reads = 1` (`fin/pipeline/config.py:137`), and
selection subtracts qualifying annotated candidates from the drop set
(`fin/pipeline/selection.py:362`):

1. **Eligibility.** Only **annotated** multi-exon candidates qualify, and at
   least K=1 distinct compatible read is required. The eligible read
   population is same-strand, mapped, non-secondary and non-fusion, **with
   MAPQ >= 1** -- compatibility depth is counted only from reads that entered
   `indels_by_read` (`fin/candidates/discovery.py:522`, `:592`), so a
   geometrically compatible MAPQ-0 read cannot supply attestation.
2. **Compatible is not uniquely-supporting.** A read's non-empty chain must
   exactly equal the candidate chain, or be a strand-aware contiguous **3'
   suffix** of it -- ascending-coordinate suffix on `+`, prefix on `-`.
   Junction wobble is not accepted in this match. Read ids are deduplicated
   per model, but **one read can attest several nested models**
   (`fin/candidates/discovery.py:162`).
3. **Span conditions.** A shortened read must not extend across the omitted
   upstream intron beyond a 10 bp slop; its 3' end may overshoot the annotated
   end by at most 50 bp. That is an overshoot limit, not a requirement to
   finish near the TES. No per-read indel-cleanliness requirement applies by
   default (`fin/candidates/discovery.py:157`, `:230`).
4. **The veto is `novel`, not merely `jitter`,** with a 10 bp window
   (`fin/pipeline/config.py:478`, `fin/scoring/m2_junction_nll.py:521`). An
   exemption is refused when:
   - an unobserved candidate junction has an observed neighbour within the
     window on **both** boundaries; or
   - an observed junction absent from the interval's annotated-junction
     universe shares **either** boundary, within the window, with an
     unobserved candidate junction; or
   - such a novel observed intron lies wholly inside a candidate exon and is
     longer than the veto window.
5. **The annotation universe includes already-masked candidates** -- novelty
   is judged against all annotated candidates in the interval, not only
   survivors (`fin/pipeline/selection.py:374`). Missing junction observations
   or invalid geometry do not establish an exemption
   (`fin/scoring/m2_junction_nll.py:494`).
6. **An `.ext` sibling inherits its annotation's compatibility depth**, which
   was measured against the original extent (`fin/candidates/discovery.py:719`).
   Extending the TES does not retroactively convert overshooting reads into
   attestation support.
7. **It changes this gate only.** It does not resurrect candidates dropped by
   earlier gates (they stay masked, `fin/pipeline/selection.py:367`), does not
   add candidates, does not move read mass, and does not bypass the global
   abundance cascade (`fin/pipeline/runner.py:352`).

> **The distinction to carry away:** under the shipped profile an annotated
> model may survive the support gate **without winning any individual M1 or
> M2 competition**, via K=1 compatibility attestation subject to the full
> novel-junction veto. That is neither blanket annotation retention, nor
> evidence that the compatible read singles out that isoform.

---

## 6. Sample-level finalization

1. **Unify extent pairs, then aggregate** across intervals.
2. **Global selection** operates on the evidence as it stood **before** the
   refit -- selection does not read its own output.
3. **Junction snap** (on): a novel junction may be replaced by a more strongly
   observed nearby one; `snap_merge_chain` then merges models that became
   structurally identical, preserving read and abundance mass. GTF and fusion
   are exempt. The snap is per junction and has no minority-support veto.
4. **Post-selection refit** (on): abundance is re-fit over the survivor set,
   redistributing the mass of dropped candidates and reporting what was
   orphaned.

### 6.1 Quantification and TPM

Abundance is the fractional EM mass. `num_reads` is the hard argmax count.
TPM by default **divides abundance by transcript length**
(`tpm_length_normalize=True`, `fin/pipeline/config.py:509`), which is the
short-read convention. `--no-tpm-length-normalize` **still normalises to one
million**, it simply stops dividing by length:
`TPM = 1e6 * abundance / sum(abundances)` (`fin/analysis/quantification.py:554`).
That is the direct-RNA convention. Disabling length normalisation does not
disable million-sum normalisation, so the TPM column is not generally the raw
abundance.

**No structural decision reads TPM.** The flag changes the reported column and
the refit diagnostics only.

### 6.2 Default outputs

`assembly.gtf` (abundance as an attribute, no TPM), `scores.tsv` (per
candidate, including TPM), `junction_support.tsv`, `abundance_refit.json`,
`run_manifest.json`, and a folded-TSS recovery report where applicable
(`fin/pipeline/runner.py:327`, `:338`, `:517`, `fin/cli.py:756`). Artifacts
are conditional -- an empty or failed run does not produce all of them. Field
definitions are in [`USAGE.md`](USAGE.md).

---

## 7. Independent fusion analysis (`--fusion`, opt-in)

Fusion runs **after an unchanged assembly** and writes to its own directory.
`fin/cli.py:593` passes `fusion_enabled=False` to the assembly config, so the
old coupled path is unreachable from the CLI and the assembly GTF never
contains a fusion transcript.

Acceptance, in order: complete multi-arm read paths are enumerated (every
coverage-qualified prefix is a hypothesis, not only leaves); junctions are
nominated within +/-3 bp and adjudicated by a fixed-window affine DP with a
margin; paths are validated; surviving molecules are chained; chains with
enough whole-chain support become fused-transcript models.

Rules worth knowing because they change what is reported:

- **A circular reference's seam is not a rearrangement.** A read crossing the
  origin of a linearised contig (chrM by default, 50 bp tolerance) is
  classified `circular_reference_origin_spanning`, and its chain never becomes
  a supported fusion however many molecules cross it.
- **`chain_support_status`** reports the support gate alone, for every event
  class -- `status` is a legacy composite whose meaning depends on which
  branch fired.
- **`event_family_id`** links reported models that describe the same junction
  geometry within tolerance. It is a geometric grouping, not a claim that the
  members are isoforms of one biological event.
- **`junction_dinucleotides`** are the reference bases at the reported
  placement. They are an annotation of coordinates, not evidence of splicing,
  and their distribution is subject to motif-selection bias because the
  nominator searches for canonical motifs.
- **`fusion_junction_evidence.jsonl`** answers, per reported junction, which
  molecule on which retained path put it exactly there. Model membership is
  **not** precise-coordinate support and the export never infers one from the
  other.

Fusion has **not** been benchmarked against any fusion caller. Signal
validation of fusion junctions is available but off by default.

---

## 8. What is off by default

Reading an option's help text as shipped behaviour is the most common way to
misdescribe this pipeline. Under `real-drna` the following are **not** active:

| Option | Default |
|---|---|
| `max_soft_mass_ratio`, `min_fulllen_fraction`, `min_polya5p_reads` | disabled by the profile |
| `junction_dominance_filter` | off (`fin/pipeline/config.py:990`) |
| `m2_support_gate_gtf_exempt` | `False` -- a *different* mechanism from §5.2 attestation |
| `read_extent_mode_tol_bp` | 0 -- ends use the outward extreme, not a mode (`:564`) |
| `tss_contradiction_min_share`, `tss_contradiction_auto_gate`, both distance floors | 0 / off (`:565`) |
| `m2_scramble_seed` | 0 -- no signal-score permutation |
| `selection_audit` | off -- candidate-fate tracing is opt-in (`:1065`) |
| `denovo_graph` | off, **and unreachable under production discovery** (`fin/candidates/discovery.py:1041` returns first) |
| `--fusion` | off |

---

## Provenance

Written against `dev` at `4932747`, 2026-09-21, with the stage structure and
worked examples of [`PIPELINE_WALKTHROUGH.md`](PIPELINE_WALKTHROUGH.md)
(2026-09-08) re-checked against the current source rather than carried over.
That document's early stages -- interval handling, poly(A) and nuisance,
discovery and collapse, the pre-assignment gates -- were re-verified as still
accurate; its account of M2 windows, the flat-tie fallback, abstention
branches, outputs and above all the support gate had gone stale and is
superseded here.
