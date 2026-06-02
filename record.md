# pyfin — experiment record (M3 / false-positive removal)

> Detailed experiment log kept in-repo to keep project memory lean.
> **All numbers are gffcompare Transcript-level F1** (`experiments/gffcompare_sweep.py`,
> SIRV p00, NO-GTF de novo — the mandated competitor metric), NOT the internal proxy.
> The proxy historically UNDERSTATES gffcompare wins, so every go/no-go below was decided
> on gffcompare. Every threshold is SIRV-tuned; nothing is hardcoded into production.
> **Production code is UNTOUCHED.** Everything here is read-only diagnostics + GTF outputs.

Competitor bar (SIRV no-GTF de novo): production baseline peak **44.7–44.8**;
ablation M1-first champion **45.2**; FLAIR 44.0 > ESPRESSO 38.5 > TALON 27.1 >
IsoQuant 25.5 > StringTie3 20.7 > Bambu 12.7.

---

## 1. Headline result — gffcompare Tx-F1 ladder

| T (abundance) | baseline | M2-margin (thr −0.05) | locus-relabund (<0.1) | stacked (both) |
|---|---|---|---|---|
| 3  | 40.8 | 42.5 | 45.2 | **46.5** |
| 5  | 44.5 | 45.9 | 45.0 | **46.2** |
| 7  | 43.8 | 45.1 | 43.7 | 44.8 |
| 10 | 44.8 | 45.1 | 44.4 | 44.6 |

- **M2-margin filter alone:** +1.4 Tx-F1 at peak (44.5→45.9 @T5), lifts EVERY T. Beats
  FLAIR (44.0) and the ablation champion (45.2). Precision-dominant: at T5 Pr +7.4,
  Sn −1.7 (~3 real TP lost on a 176-tx truth). Dropped 45 of 523 novel multi-exon cands.
- **locus-relabund alone:** strongest at low T (45.2 @T3) — minor-isoform suppression.
- **Stacked (M2-margin + relabund<0.1):** peak 46.5 @T3.
- **★ SUPERSEDED by §3c:** loosening locus-relative abundance ALONE (signal-free) reaches peak
  **47.8 @T3 at zero recall cost (relabund<0.2)**, climbing to **49.6 @T3 at relabund<0.4** — and
  M2-margin signal turns out to be REDUNDANT with it. See §3c for the full ladder + SIRV caveat.

---

## 2. What the precision-hurting FP actually are

The remaining FP that cost precision (T≥5 slice: 49 TP / 63 FP; 59/63 are M2-margin-reachable):

- **NOT assembly chimeras.** Full-chain phasing AUC = 0.56 (≈ coin flip). `full_frac`
  (fraction of assigned reads carrying the full chain) = 0.27 for FP — IDENTICAL to TP's 0.27.
  They are carried by real reads.
- They are **structurally-real-looking MINOR 2-intron alt-splice paths**: wrong *combinations*
  of individually-real junctions, or genuine-but-minor isoforms that should be suppressed.
- **All per-structure filters plateau at AUC ~0.6:** per-junction raw support 0.62,
  consecutive-pair support 0.63, full-chain phasing 0.56. A per-junction-support test cannot
  see them *by construction* (each junction IS supported).

Only **signal** and **abundance** separate them, and they separate DIFFERENT FP (orthogonal):
- **M2-margin (signal)** — AUC 0.80.
- **locus-relative abundance** — AUC 0.695.

---

## 3. The two working levers (detail)

### 3a. M2-margin M3 post-filter — AUC 0.80
Per NOVEL multi-exon candidate C:
- cluster wobble siblings (`cluster_candidates_by_chain`, wobble=20),
- build shared K=10 junction window (`class_junction_window_set`),
- per C-assigned read: `margin = NLL(best sibling) − NLL(C)` (>0 ⇒ read prefers C),
- `mean_margin` over ≥2 scored reads.
- **Drop** reachable C with `mean_margin < thr` (thr = −0.05). **Unreachable NEVER dropped**
  (singleton / no shared window / <2 scored reads).

Reachability ceiling: ~44% of high-count FP are disjoint/singleton → invisible to this method.
That is its hard limit; it still nets +1.4 because the reachable FP it does remove are clean.

### 3b. Locus-relative abundance — AUC 0.695 (pure GTF, no signal)
`relabund(C) = C.abundance / max(abundance over overlapping novel candidates)`.
- FP median 0.226 vs TP median 0.519 → minor-isoform suppression.
- Drop novel multi-exon C with `relabund < 0.1`.
- Orthogonal to M2-margin (signal vs count) → stacks cleanly.

---

## 3c. FOLLOW-UP RESEARCH — "how else to remove these FP" (signal-free feature mining)

Tool `benchmarks/diag_fp_gtf_features.py` (read-only, host, NO f5c): label every novel
multi-exon candidate in the baseline `assembly.gtf` TP/FP by intron-chain vs SIRV truth, mine
features derivable from the GTF alone (abundance + num_reads), measure AUC, emit filtered GTFs,
measure gffcompare.

**Tested axes (full set 64 TP / 976 FP, then HARD slice abundance≥T):**
| feature | full-set AUC | T≥3 AUC | verdict |
|---|---|---|---|
| `relabund` (locus-relative abundance) | 0.918 | 0.876 | **dominant lever** |
| `mean_weight` = abundance/num_reads (argmax_keep mean 1/K fractional weight) | 0.889 | 0.796 | NOT orthogonal — relabund dominates it |
| `abundance` (raw count floor) | 0.915 | — | the gffcompare T sweep itself |
| `n_introns` | 0.505 | — | coin flip (confirms user-rejected gate) |

- **`mean_weight` is NOT a new lever.** Looked strong on the full set, but at T≥3 `relabund<0.2`
  removes 52 FP / 0 TP whereas `mean_weight<0.3` removes 49 FP but COSTS 5 TP. It's the same
  minor/borrowed-candidate mechanism in count form; relabund strictly dominates. Joint AND is worse
  than relabund alone. Dead end as an independent axis.

- **★ `relabund` is FAR stronger than previously credited and it is SIGNAL-FREE.** gffcompare ladder
  (drop novel multi-exon with relabund < thr; pure GTF, no f5c):

  | relabund thr | T3 | T5 | T7 | T10 | Sn@T3 | dropped |
  |---|---|---|---|---|---|---|
  | baseline | 40.8 | 44.5 | 43.8 | 44.8 | 40.9 | 0 |
  | **<0.2** | 47.8 | 47.0 | 45.6 | 44.8 | **40.9 (ZERO recall cost)** | 945 |
  | <0.3 | 49.1 | 48.2 | 46.7 | 45.7 | 39.8 | 961 |
  | **<0.4 (SIRV peak)** | **49.6** | 48.7 | 47.2 | 46.1 | 38.6 | 972 |
  | <0.5 | 48.1 | 47.1 | 45.5 | 44.4 | 36.4 (recall breaks) | 980 |
  | M2-margin ∩ relabund<0.2 | 48.6 | 47.7 | 46.2 | 45.1 | 39.2 | — |

  At `<0.2`, **Sn is byte-identical to baseline at every T → zero recall cost on SIRV**, pure
  precision (T3 Pr 40.7→57.6), peak **47.8**. The F1 curve climbs to **49.6 @T3 at <0.4**, then
  falls at <0.5 as recall finally breaks. relabund's value is that it is LOCUS-relative, so it
  decouples precision from the absolute count floor: it lets you run at low T (high recall) while
  suppressing minor-isoform FP — biggest gains are at low T (+7→+9 @T3), converging to baseline at
  T10.

- **★ M2-margin (f5c signal) is REDUNDANT with relabund.** Stacking M2-margin on relabund<0.2
  (48.6) is BELOW relabund<0.3/<0.4 alone (49.1/49.6). Both filters target the same minor/borrowed
  candidates; the expensive signal filter is dominated by the free count filter. The previously
  recorded "stacked 46.5" was M2-margin + relabund**<0.1** (too-tight relabund); loosening relabund
  alone surpasses it without any signal.

- **★★ SIRV-ARTIFACT WARNING (strongest of any lever found).** The F1 curve peaks at relabund<0.4
  (drops 972 of 1129 candidates) precisely BECAUSE SIRV has almost no genuine minor isoforms —
  suppressing minor isoforms is nearly free here. The recall break only appears at <0.5. On a REAL
  transcriptome with real minor/low-abundance isoforms this same threshold would be catastrophic for
  recall. relabund is the most powerful AND the most SIRV-overfit lever found → if ever wired:
  configurable, default conservative (≤0.2) or OFF, re-tuned on real data. Do NOT hardcode.

- **★ PRIOR ART — this is a standard, biologically-motivated filter, NOT a pyfin invention.**
  "Minimum isoform fraction relative to the dominant isoform at a locus" is a 20-year-old textbook
  heuristic. `relabund` is literally the same quantity:

  | tool | parameter | default | meaning |
  |---|---|---|---|
  | Cufflinks | `--min-isoform-fraction` / `-F` | **0.10** | drop isoforms below 10% of the locus's major isoform |
  | StringTie (current) | `-f` | **0.01** | same, 1% |
  | StringTie (legacy v1.3.x) | `-f` | 0.10 | |
  | StringTie `--conservative` | `-f` | 0.05 | conservative mode raises to 5% |

  **Biological rationale (from the Cufflinks/StringTie docs themselves):** the low-fraction tail is
  enriched for NON-real transcripts — primarily **incompletely-spliced precursors (pre-mRNA)** that
  look like an isoform missing/adding an intron, plus RT/template-switching artifacts, mapping
  errors, RNA-degradation fragments, and stochastic splicing noise (Pickrell 2010; Melamud & Moult).
  So suppressing the minor-isoform tail is principled, not a metric hack.

- **★★ RECOMMENDED THRESHOLD vs the SIRV-optimal threshold — the key distinction.** The *strategy*
  is sound; the SIRV-*optimal operating point* is overfit. Real tools use **1% (StringTie) to 10%
  (Cufflinks)**. The SIRV optimum here is `relabund<0.4` (40%) — 4–40× more aggressive — and even
  `<0.2` is "free". That is only possible BECAUSE SIRV has no genuine low-abundance isoform tail
  (synthetic, few isoforms/gene). On a real transcriptome the 20–40% band is full of REAL
  low-abundance functional isoforms (tissue/condition-specific, NMD/PTC targets that are natively
  low-abundance, regulatory minor isoforms) → 0.4 would gut recall. The "F1 keeps rising to 0.4"
  curve is itself EVIDENCE of the SIRV artifact (it measures SIRV's missing tail, not filter
  quality). **If wired to production: adopt the CONCEPT + the literature threshold (default ~0.01–0.1,
  aligned with StringTie `-f`/Cufflinks `-F`), configurable, NEVER the SIRV-tuned 0.4.** On real data
  the F1 peak is expected to land in the 1–10% band with a far more modest gain.

- **Irreducible prize unchanged:** at T≥10 even relabund<0.4 only reaches 46.1; the high-count
  structural FP (alt-assemblies of individually-real junctions) survive every count/structure/signal
  lever. Only a splice-site PWM/MaxEnt (real-data, §4) is likely to touch them.

**Net research answer:** beyond M2-margin, the one lever that materially removes these FP is
**locus-relative abundance** (signal-free, dominant, but heavily SIRV-tuned). mean_weight and
n_introns are dead ends; M2-margin is redundant with relabund. New artifacts:
`experiments/diag_m3_margin_filter/assembly_relabund0{2,3,4,5}.gtf`, `assembly_m2_relab02.gtf`.

## 4. Real-data lever flag (untestable on SIRV)

The structurally-real minor alt-paths are exactly what a **splice-site strength model
(PWM / MaxEnt)** would score down — a spurious donor/acceptor that is canonical-dinucleotide
(passes the GT-AG gate) but weak in the full splice-site PWM context. This CANNOT be tested on
synthetic SIRV (its splice sites are designed, not drawn from a genomic PWM). Flag for when a
real transcriptome arrives: splice-site PWM/MaxEnt scoring at discovery is the most likely
generalizable precision lever beyond the canonical-dinucleotide gate.

---

## 5. Artifacts

- Script (read-only diagnostic, drives production candidate path):
  `benchmarks/diag_m3_margin_filter.py`
- GTFs in `experiments/diag_m3_margin_filter/`:
  - `assembly.gtf` — baseline (1129 cands, 523 novel multi-exon)
  - `assembly_filtered.gtf` — M2-margin filter applied (−45 cands)
  - `assembly_relabund.gtf` — locus-relabund < 0.1 applied
  - `assembly_both.gtf` — stacked (both filters)
- Earlier read-only diagnostics: `benchmarks/diag_m2_margin_cluster.py`,
  `benchmarks/diag_m3_event_region.py`, `benchmarks/diag_m3_coherence.py`

---

## 6. Status / open decisions

- M2-margin gffcompare go/no-go = **POSITIVE (+1.4 peak)**. Stacked with relabund = +1.7.
- Optional next step before wiring: sweep thresholds (M2-margin thr ∈ {−0.02,−0.05,−0.08,0};
  relabund ∈ {0.05,0.1,0.15}) to locate the gffcompare-optimal operating point.

### ★ relabund WIRED TO PRODUCTION as `min_isoform_fraction` (default 0.01, configurable)
User go-ahead given. Wired as the standard Cufflinks `--min-isoform-fraction`/StringTie `-f`
heuristic, NOT the SIRV-tuned 0.4.
- **Config:** `PipelineConfig.min_isoform_fraction: float = 0.01` (0.0 disables). Default =
  StringTie's shipped default, recall-safe for long reads.
- **Code:** `fin/analysis/quantification.py::isoform_fraction_drops()` (pure, bucketed by
  chrom/strand, locus-dominant never dropped); applied in `runner._finalize_and_write` as a
  post-quant candidate FILTER gated by `enable_score_filter AND min_isoform_fraction>0`, AFTER
  combined_score filter, BEFORE gene_id resolution. Exempts gtf/fusion/mono. NEVER touches
  EM/assignment. Tests: `tests/unit/test_isoform_fraction.py` (11). Codex review: **APPROVE**.
- **★ HONEST gffcompare result at the literature thresholds (SIRV no-GTF, baseline = current
  prod assembly):**

  | min_isoform_fraction | T3 | T5 | T7 | T10 | vs baseline |
  |---|---|---|---|---|---|
  | baseline (off) | 40.8 | 44.5 | 43.8 | **44.8** | — |
  | **0.01 (default)** | 40.8 | 44.5 | 43.8 | 44.8 | **identical — no-op on SIRV** |
  | 0.05 | 40.8 | 44.5 | 43.8 | 44.8 | identical — no-op |
  | 0.10 | 43.9 | 44.7 | 43.8 | 44.8 | +3.1 @T3, +0.2 @T5; peak unchanged |

  **At 0.01–0.05 the SIRV gffcompare gain is ZERO.** The 604 candidates dropped at 0.01 are all
  abundance<3, so the gffcompare T≥3 abundance floor already culls them → the filter is REDUNDANT
  with the count floor at these operating points. The large gains (47.8 @<0.2, 49.6 @<0.4 in §3c)
  require the SIRV-overfit regime. So `min_isoform_fraction=0.01` is shipped as a **real-data
  correctness/insurance feature** (catches a high-abundance minor that's still <1% of a very
  dominant neighbor — a case the absolute count floor cannot see), NOT a SIRV benchmark win. It is
  gffcompare-neutral on SIRV by design (does no harm), and the user may raise it (e.g. 0.10 for the
  visible low-T SIRV gain, or re-tune on real data). NEVER set the SIRV-optimal 0.4 in production.

### USER DECISION — do NOT re-propose
- `n_introns ≥ N` structural gate is REJECTED (+0.8 on SIRV but a SIRV artifact; real genes
  legitimately have many introns).
- M2/M3-as-primary-assignment/quant signal is a settled DEAD END in every form (full-matrix,
  tie-restricted, EM, coherence, M1-gated cascade; GTF-in and GTF-free) — all lose to plain M1
  argmax. M3 is only ever a post-hoc FILTER.
