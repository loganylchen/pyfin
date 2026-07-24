# pyfin — current production state (2026-07-22)

Snapshot of the DEFAULT `fin` CLI pipeline as it stands on `dev` (HEAD `828565a` + the
uncommitted generation/quant body about to be committed). Written to give a faithful
picture of "what production is" after this session's changes.

## 0. The two operating points (IMPORTANT)

The `fin` CLI defaults are **SIRV-tuned**. Several finalize gates that lift synthetic-SIRV
Tx-F1 **halve real-dRNA recall** (dRNA is 3'→5' so genuine isoforms are 5'-truncated and
often lack a clean polyA 5' end). So there are effectively two profiles:

| profile | how invoked | min-fulllen-fraction | min-polya5p-reads | min-abundance |
|---|---|---|---|---|
| **SIRV/default** | raw `fin` CLI defaults | 0.1 | 1 | 3 |
| **real dRNA** (our p00 benchmark) | `--min-fulllen-fraction 0 --min-polya5p-reads 0` | 0 | 0 | 3 (or 0) |

All p00 numbers this session are the **real-dRNA** operating point (SIRV gates off). The
generation/structural changes below are default-ON and help **both** profiles.

## 1. Pipeline flow (default quant_mode = m2_em)

```
reads(BAM) ─▶ intervals ─▶ GENERATION ─▶ pre-EM gates ─▶ M1/M2/EM quant
                                                             ─▶ post-EM structural levers
                                                             ─▶ finalize filters ─▶ GTF
```

### 1a. Intervals & fusion exclusion  (io/interval_manager.py)
- Strand-separated per-locus intervals.
- `is_fusion_read` excludes a read from generation if supplementary OR a terminal
  soft-clip **≥ 250 bp** (was 50; raised this session — commit `828565a` — dRNA poly-A/
  adapter tails run to ~200bp and were mis-flagged). Same threshold gates the `--fusion` caller.

### 1b. Generation  (candidates/discovery.py `_chain_cluster_candidates` + chain_cluster.py) — `chain_cluster_discovery=True`
1. Group reads by **exact intron chain** (3' ends ignored).
2. **Exact-subchain fold** into the longest container, WITH **span-guard** (`fold_span_guard=True`,
   new): a read folds only if its span does NOT run exonically across a container's EXTRA
   intron; a read that spans it (retained-intron / alt isoform) is kept as its own candidate.
3. **Mono-fold** (`fold_monoexon_contained=True`, new): a single-exon read wholly inside a
   multi-exon candidate's exon folds into it (5'/3' degradation fragment); intronic/uncontained
   mono reads stay separate.
4. **Union-find cluster** survivors by wobble/cassette(≤70bp)/containment (bp=6); ALL members kept.
5. **Mono reads split by locus** (new, `_spatial_read_clusters`) and each locus overlap-matched
   to a single-exon GTF (else novel) — distinct single-exon genes no longer collapse.
6. Novel candidate span = union of its reads' extents; sequence stitched from genome.
   GTF transcripts added as `source="gtf"` candidates. Retained fold shadows carried for recovery.
- `canonical_search_bp` (=4) is IGNORED when chain-cluster is on.

### 1c. Pre-EM gate  (runner.py)
- **canonical_gate** ON: drop a NOVEL multi-exon candidate if any junction motif ∉
  {GT-AG, GC-AG, AT-AC}. GTF/fusion/mono exempt.
- junction_dominance_filter: OFF.

### 1d. Quant  (quant_mode="m2_em", runner.py + scoring/m2_junction_nll.py)
- **M1**: mappy AS argmax → per-read best-AS tie set + mappability mask.
- **M2**: junction-NLL tie resolution (`m2_tiebreak` ON, margin 1e-9; `m2_diff_cover_gate`
  ON margin 0.5; `m2_tie_scoregate_split` ON). The **summed-LLR** metric (sum per-event NLL
  over tight differing-junction windows) is used by the runner's M2 path.
- **EM**: seeded by the M2 junction-NLL distance; RSEM-style. `m3_coherence` OFF,
  `abundance_feedback` OFF. `em_max_iter` 1000, `em_tol` 1e-4.

### 1e. Post-EM structural levers  (runner.py, all default-ON, m2_em only)
- **m2_cluster_recheck** (bp=20, fraction=0.15, cassette=70): cluster multi-exon by structure;
  drop a NOVEL wobble shadow whose EM abundance < fraction×anchor. `novel_displaces_gtf` ON
  (a GTF sibling is droppable only if its distinguishing junction has < 1 exact read; jct_tol=0).
- **novel_junction_min_reads=2**: drop a NOVEL multi-exon candidate if any junction is spliced
  by < 2 directly-observed reads.
- **m2_support_gate** ON (tie-accept): keep a multi-exon candidate iff it is some read's M1
  sole-best-AS OR M2-best.
- **containment_cluster** ON (ab_ratio 0.3 AND read_ratio 0.3, cap 10 reads): drop a NOVEL
  sub-chain truncation/exon-skip shadow that is a low-support minor of a longer candidate.
- containment_collapse: OFF.

### 1f. Finalize filters  (runner.py `_finalize_and_write`, in order; `fin` CLI defaults)
1. **min_abundance=3** (NOVEL EM abundance) / min_gtf_abundance=1 (GTF soft-mass); floor_gtf_abundance OFF.
2. **min_isoform_fraction=0.01** (NOVEL multi-exon vs dominant locus isoform).
3. **max_soft_mass_ratio=2.0** (NOVEL multi-exon soft/hard read ratio).
4. **min_fulllen_fraction=0.1** (NOVEL multi-exon full-length read fraction; window 25, min_reads 4). ← SIRV gate, off for dRNA.
5. **min_polya5p_reads=1** (NOVEL polyA+5' evidence; needs signal; polya5p_exempt_gtf ON). ← SIRV gate, off for dRNA.
- drop_mono_exon_novel OFF; guided_junction_min_reads 0 (OFF).

## 2. This session's changes (all default-ON generation/quant)

| change | status | p00 effect (m2_em, dRNA gates-off) |
|---|---|---|
| fusion soft-clip 50→250 | **LANDED** `828565a` | =6278→6323, Pr 62.6→66.0 |
| mono-fold (fold single-exon fragments in) | uncommitted, default-ON | =6323→6414, Pr 66.0→69.9 |
| span-guard (keep retained-intron reads separate) | uncommitted, default-ON | =6414→6602, Pr 69.9→71.7 |
| mono-split (per-locus mono + overlap GTF match) | uncommitted, default-ON | benchmarking (job 211274) |
| shadows + summed-LLR M2 + recovery scaffold | uncommitted | (infra; summed-LLR active in runner) |

Cumulative p00: `=` 6278→6602 (+324), Pr 62.6→71.7 (+9.1pt), output shrunk (de-fragmentation).

## 3. Default-OFF experiments (present in code, not in production)

- `quant_mode="cluster"` (per-cluster M3/EM + 5'-TSS recovery; cluster_quant.py) — LOST to m2_em.
- `denovo_wobble_tol` (wobble merge), `denovo_graph` (intron-graph assembly), junction_graph.py — de-novo experiments.
- `junction_dominance_filter`, `containment_collapse`, `abundance_feedback`, `m3_coherence`, `guided_junction_min_reads` — off.

## 4. Known limitations

- **Exact-subchain fold of genuine short isoforms** (nested truncations with truncated reads):
  span-guard fixes the retained-intron subset; genuine short isoforms whose reads are ALSO
  truncated still fold (retained as shadows; downstream 5'-TSS recovery is the intended path).
- **Selection is the recall bottleneck**: p00 generation ceiling ~11310 `=` vs output ~6602 —
  the `min_abundance=3` floor kills low-abundance isolated TPs (63% have <3 reads). The proper
  fix is §14-style candidate-set model selection (not yet built).
- Finalize gates 4/5 are SIRV-tuned (see §0).
