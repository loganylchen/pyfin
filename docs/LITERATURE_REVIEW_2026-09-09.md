# Literature review: raw signal, guided emission, degradation-aware assignment, evaluation contracts (2026-09-09)

Commissioned for the "best on every axis" campaign (HANDOFF 4o). Produced by a
web-search agent from primary sources; every claim carries its URL. Items are
tagged [A] actionable / [B] background. Read alongside
`PRODUCTION_JUNCTION_ACCURACY_BENCHMARK.md` (our own raw-DTW result) and
`PRODUCTION_CLUSTER_DTW_REPLAY_2026-09-09.md`.

# Signal, guided-emission, and evaluation literature for PyFIN

Legend: **[A]** = directly actionable · **[B]** = background.

---

## 1. Raw signal for splice/structure decisions in dRNA (2022–2026)

**NanoSplicer** is still the only published tool that lets raw signal *decide* a splice junction. It extracts the junction-within-read (JWR) squiggle, enumerates canonical GT-AG candidates within ±10 bp, predicts candidate squiggles with Tombo, and DTW-aligns observed→candidate. ([Bioinformatics 2022](https://academic.oup.com/bioinformatics/article/38/15/3741/6594111), [PMC9344838](https://ncbi.nlm.nih.gov/pmc/articles/PMC9344838), [github.com/shimlab/NanoSplicer](https://github.com/shimlab/NanoSplicer))

**[A] How it avoids the "longer template wins" bias — three mechanisms PyFIN does not currently use:**
1. Alignment quality is the **average log-likelihood over the nucleotides of the candidate squiggle**, not a sum (Suppl. S1.7). Per-base normalisation makes candidates of different length comparable. PyFIN's *summed* junction NLL is exactly the un-normalised form.
2. Scores are converted to a **posterior over the candidate set**, output as `prob_uniform_prior` and `prob_seq_pattern_prior` (a motif prior) plus `best_prob` — i.e. a normalised competition among a closed hypothesis set, with an explicit prior slot.
3. An abstention gate: **SIQ** = max alignment quality across candidates; JWRs below threshold (docs suggest SIQ < −0.8) are refused, and the accuracy claim is stated only for SIQ > 0.8 **and** best posterior > 0.8. This mirrors your calibrated-abstention tie-breaker result.

**[B] Caveats:** validated on **cDNA / R9 only**; the authors defer dRNA and R10 to future work. An AMSI student project probing dRNA failure modes found minimap2 ≥ NanoSplicer on 498 classified dRNA reads, attributing failures to RNA's larger dwell-time variance and Tombo resquiggling ([srs.amsi.org.au](https://srs.amsi.org.au/student-blog/nanosplicer/)). **This is independent corroboration of your internal DTW-vs-BAM result.** Tombo is deprecated and RNA002/R9-only.

**Signal-space mappers do not cover this problem. [B]** RawHash2 ([Bioinformatics 2024](https://academic.oup.com/bioinformatics/article/40/8/btae478/7723993)) and Sigmoni ([Bioinformatics 2024](https://academic.oup.com/bioinformatics/article/40/Supplement_1/i287/7700860)) are DNA-only; RawHash2 lists RNA explicitly as *future work* and has no splice-aware chaining ([arXiv v4](https://arxiv.org/html/2309.05771v4)). No basecalling-free isoform caller exists — the nearest analogue is tRNA **isodecoder** discrimination from ionic current ([bioRxiv 2025](https://www.biorxiv.org/content/10.64898/2025.12.27.696718v2.full)).

**Length-bias handling in the wider signal literature [A]:**
- Normalised DTW `nMetric = Metric / gold_length`, because absolute DTW distance grows proportionally to template length ([PLoS Comput Biol 2021](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1009350)).
- Thresholding a raw best-match score has "very low specificity at higher sensitivity… especially when matching against a longer reference" — argues for an **LLR against a null/background model on the same signal**, so per-sample terms cancel ([arXiv 1705.01620](https://arxiv.org/pdf/1705.01620)).
- seq2squiggle normalised its DTW metric for read-length differences ([Bioinformatics 2025](https://academic.oup.com/bioinformatics/article/41/1/btae744/7930676)).

**Duration/dwell models [A/B]:** **SegPore** (eLife 2026, [104618](https://elifesciences.org/articles/104618)) uses a **hierarchical HMM** with a "jiggling translocation" hypothesis, segmenting into base *and transition* blocks — the only current model that treats dwell explicitly, and it beats Nanopolish/Tombo; RNA004-capable. **Uncalled4** (Nat Methods 22:681, 2025, [s41592-025-02631-4](https://www.nature.com/articles/s41592-025-02631-4)) provides basecaller-guided DTW and a reproducible **de novo pore-model training** method — the practical replacement for Tombo. **f5c ≥1.4** supports `--pore rna004` ([hasindu2008.github.io/f5c](https://hasindu2008.github.io/f5c/)). Nanopolish/f5c HMMs carry explicit `P(stay)`/`P(skip)` transition priors ([CRANE, arXiv](https://arxiv.org/pdf/2603.20420)).

**Simulators as hypothesis generators [B]:** Squigulator ([Genome Res 34:778, 2024](https://genome.cshlp.org/content/34/5/778)) is the only simulator supporting **RNA004** and has a noise-free "ideal" mode that emits an expected squiggle for any candidate sequence — but no scoring module; NanoSimFormer reports Squigulator at 81.3% median accuracy on HG002-RNA004 vs its own 99.3% ([PMC13312124](https://pmc.ncbi.nlm.nih.gov/articles/PMC13312124/)). **[A]** ONT's RNA004 model is a **9-mer, in standard (normalised) units, stored 5'→3' while signal is 3'→5'** ([nanoporetech/kmer_models](https://github.com/nanoporetech/kmer_models)); Poregen offers an open de novo alternative ([Bioinformatics 2025](https://academic.oup.com/bioinformatics/article/41/4/btaf111/8078598)).

**No source found** for any tool repurposing Nanocompore/xPore/m6Anet-style comparative signal tests for transcript *structure*.

---

## 2. 5' end / cap / TSS from dRNA signal

**[A] The 5' end of a dRNA read is structurally untrustworthy.** ~11 nt are lost when the motor releases the strand — confirmed on ERCC and IVT controls, so it is not degradation ([eLife 49658](https://elifesciences.org/articles/49658); [PMC8906549](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8906549/)). Structure-probing work sees a matching 10–15 nt shift. Standard DRS has **no marker for the native 5' end**, so "full-length" cannot be called from a plain read.

**[B] Cap detection is chemical, not signal-based.** NRCeq swaps m7G for a 45-nt adapter; NanoCount-style screening of capped+polyadenylated reads gave 294,107 full-length scaffolds from ~4M GM12878 reads (~7%) ([NAR 50:3475](https://academic.oup.com/nar/article/50/6/3475/6542485)). CapTrap-seq is the cDNA analogue ([Nat Commun 2024](https://www.nature.com/articles/s41467-024-49523-3)). Raw-signal m7G work targets *internal* m7G, not the cap. **No signal-only cap/TSS classifier found.**

**[A] The adapter/motor region *is* a distinguishable signal segment** (open-channel → discrete adapter segments → translocation) and is what read-start segmenters key on ([review PMC11388133](https://pmc.ncbi.nlm.nih.gov/articles/PMC11388133/)) — but nobody has published using it to grade 5' completeness. **[A]** The mature signal 3'-end analogue is Dorado `--estimate-poly-a` (`pt:i` tag; −1 = no anchor, 0 = anchor but no length), benchmarked as the recommended tool on synthetic standards, with a known spurious ~0–10 nt peak from fragmentation/mispriming ([Dorado docs](https://software-docs.nanoporetech.com/dorado/latest/basecaller/polya_estimation/); [GigaScience 2025](https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giaf098/8246622)). Note **version-to-version instability** in Dorado polyA lengths ([dorado#1562](https://github.com/nanoporetech/dorado/issues/1562)).

---

## 3. How guided tools emit low-support annotated transcripts

Two paradigms.

**Hard annotation bypass:**
- **TALON**: `talon_filter_transcripts` whitelists **(a) Known, or (b)** novel with ≥5 reads in all datasets; internal-priming `--maxFracA 0.5` applies to novel only. **Known transcripts skip the filter — 1 read suffices** ([github](https://github.com/mortazavilab/TALON)).
- **StringTie `-e`**: quantify-only over `-G`. v2.2.1 and **v3.0.0 (2025-01-07)** fixed low-coverage transcripts not being reported; **v3.0.2 (2025-10-16)** fixed `-G` transcripts with *zero* coverage not being emitted ([releases](https://github.com/gpertea/stringtie/releases)).
- **Bambu**: annotated transcripts carry `NDR = NA` and are never NDR-filtered ([Nat Methods 2023](https://www.nature.com/articles/s41592-023-01908-w)).
- **FLAIR2 `--annotation_reliant`**: aligns reads to the annotated transcriptome *first*, discovering novel isoforms only from leftovers ([Genome Biol 25:173](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03301-y)).

**[A] Lower-but-nonzero bar — the exact rule you asked for:**
- **IsoQuant 3.x** `default_ont`: **`min_known_count=1`, `min_nonfl_count=3`, `min_novel_count=3`** (`isoquant.py::set_model_construction_options`). Known non-FL spliced isoforms additionally need ≥1 left- and ≥1 right-support read and presence in the intron graph. `sensitive_ont` = 1/2/3; `pacbio_all` is documented as recovering "lowly-covered or **5'-truncated** genes". Truncation is a first-class event class (`ism_5/3`, `fake_terminal_exon_5`, `major_exon_elongation_5`); `--transcript_quantification` defaults to **`unique_only`**, matching delta **6** for ONT ([docs](https://ablab.github.io/IsoQuant/cmd.html), [formats](https://ablab.github.io/IsoQuant/formats.html)).
- **ESPRESSO**: `-N 2`, `-R 0` gate **novel** isoforms only; FSM/annotated detection has no equivalent count cutoff, ends need not match (±`SJ_dist 35`).
- **Mandalorion** `-R 5`; **isONform** `--iso_abundance 3` (recommended; shipped default not verified).

**[A] StringTie's ordinary `-G` guided assembly is the outlier — guides get *no* coverage relaxation.** Defaults `-c 1` (multi-exon), `-s 4.75` single-exon (`-L` forces `-s 1.5 -g 0`), `-f 0.01`, `-j 1` ([manual](https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual)). StringTie3 adds nascent-aware `-N` and a reworked long-read module removing poly(A)-priming artifacts, and v3.0.2 credits "a feature that uses unaligned poly(A) tails to better anchor transcript ends" — **this, not an annotation prior, is the likeliest source of StringTie3's recall edge on 5'-degraded reads.**

**Benchmarks [A]:** LRGASP found tools with high reference-TSS/TTS match but **low CAGE/QuantSeq support**, concluding pipelines are "heavily guided by reference annotations when completing transcript sequences" — the published statement of "free recall". Mandalorion had the most consistent CAGE support. LRGASP also concedes its own SIRV benchmark "might be biased" because the annotation was available ([Nat Methods 21:1349](https://www.nature.com/articles/s41592-024-02298-3)). FLAIR/FLAMES are the most annotation-sensitive: true-positive reads get "corrected" away when the guide is incomplete ([Nat Commun 2024](https://www.nature.com/articles/s41467-024-48117-3)).

---

## 4. Degradation-aware assignment

**[A] LIQA is the only published truncation-aware likelihood.** A read is "complete" if its end is within 20 bp of an annotated 5' end, else **right-censored**; read-length CDF estimated by **Kaplan–Meier**, giving emission `f(L_r | iso=i) = f(L_r) / Σ_{l_j > L_r} θ_j` — a **risk set over only isoforms longer than the read**, weighted by current abundance. This *is* the short-vs-truncated-long identifiability model ([Genome Biol 2021](https://link.springer.com/article/10.1186/s13059-021-02399-8)).

**[A] oarfish** frames your exact figure — reads compatible with T1 and T2 where "T2 is a proper prefix of T1" — and adds a coverage term `p(δ_i,a)=1/(1+e^{−aδ_i})`, `a=2`, over binned segments (`--model-coverage`). **But the sign rewards *under*-covered segments: it is a coverage-uniformity prior, not a 5'-decay prior**, and it explicitly carries no effective length ([Bioinformatics 41:i304, 2025](https://europepmc.org/article/MED/40662837)).

**[A] NanoCount has no effective-length term at all** — compatibility is initialised `1/n_alignments` and updated to `θ_ref / Σθ`; `transcript_length` is metadata only. Options are `max_dist_3_prime=50`, `max_dist_5_prime=-1` (off), `min_query_fraction_aligned=0.5`. **There is no `--three_prime_shift`.** ([NAR 50:e19](https://academic.oup.com/nar/article/50/4/e19/6439677)). `salmon --ont` **suppresses** length correction; **lr-kallisto** is the only tool estimating effective length empirically (`mean aligned length − k`) ([PLoS Comput Biol](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1013692)). IsoQuant detects truncation ambiguity but splits ambiguous reads with equal weights. **TranSigner** weights compatibility by 3'/5' read-end position on the transcript (`--psw`) ([Genome Biol 26:257, 2025](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03723-2)).

**[A] The gap is real and narrow:** no published EM conditions on 5'-start position given a fitted per-sample degradation profile. **INDEGRA** supplies a ready-made estimator — a fragmentation rate κ (cuts/nt) → per-transcript DTI index separating biological decay from technical bias ([bioRxiv 2024](https://www.biorxiv.org/content/10.1101/2024.12.12.627949v1)). Degradation magnitude and RIN-based correction: [NARGAB lqad060](https://academic.oup.com/nargab/article/5/2/lqad060/7192649).

**[A] Do not validate with stock NanoSim** — it places read starts **uniformly**, containing none of the 5' bias in question; use the LRGASP fork with empirical 5'/3' truncation distributions ([lrgasp-simulation](https://github.com/LRGASP/lrgasp-simulation/)).

---

## 5. Evaluation contracts

**LRGASP** disclaims real-sample truth outright: "Without a ground truth in cell-line data, we gauged method accuracy using spike-ins, simulated data and GENCODE manual curation." Four truth sets: SIRV-Set 4 (E0 only, + 15 long SIRVs, + 92 ERCCs; **no Sequins**), a *newer* GENCODE than participants received, 50 undisclosed manually annotated genes, and simulations. Orthogonal evidence = CAGE (5'), QuantSeq (3'), Illumina junctions, ±50 nt to reference TSS/TTS. **Evaluation used SQANTI3; gffcompare is never mentioned.** Metrics include **redundancy**, **non-redundant precision**, PDR ([Nat Methods 21:1349](https://www.nature.com/articles/s41592-024-02298-3)).

**SG-NEx** (Nat Methods 22:801, 2025, [s41592-025-02623-4](https://www.nature.com/articles/s41592-025-02623-4)) — includes your HEYA8/H9 lines; truth = six spike-in sets (Sequin V1/V2 Mix A, ERCC, **SIRV E0 and E2**, long SIRVs; direct RNA excluded from spike-in comparisons). For biological samples it states plainly there is a "lack of a ground truth"; orthogonal anchors are matched Illumina, in-silico fragmentation, and **RT-qPCR/dPCR on 13 discordant genes in MCF7** (dPCR r = 0.97 long-read vs −0.6 short-read).

**No consensus read-support cutoff for "expressed" exists** — no source found. The 2025 Genome Research review states SQANTI3-style curation "cannot provide formal benchmarking capabilities" ([PMC12047247](https://pmc.ncbi.nlm.nih.gov/articles/PMC12047247/)).

**[A] gffcompare `=` critiques — documented, from primary sources:**
- `=` is intron-chain match and **"matching transcripts are allowed to differ on the length of the first and last exons"** — TSS/TES-blind *by design* ([F1000Res 2020, PMC7222033](https://pmc.ncbi.nlm.nih.gov/articles/PMC7222033/); [manual](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)).
- An opt-in fix exists and is **off by default**: `--strict-match` requires ends within `-e` (default 100 bp), else emits `~` ([gffcompare.cpp](https://github.com/gpertea/gffcompare/blob/master/gffcompare.cpp)).
- **Sibling collapse is explicit**: "GffCompare **removes isoforms with the same SJC but different ends within each sample before comparison**" ([Volden & Vollmers, Genome Biol 2023, PMC10351160](https://pmc.ncbi.nlm.nih.gov/articles/PMC10351160/)); identical-chain query siblings are deduplicated (precision unpenalised) while `c`-contained siblings are retained. `.refmap` is one row per reference with a *list* of queries.
- Intron-chain Sn/Pr "completely ignores single-exon transcripts."
- A 2026 barley dRNA benchmark found 5' truncation was the dominant failure mode, caught by **SQANTI3 ISM proportions and not by gffcompare** ([PMC12967217](https://pmc.ncbi.nlm.nih.gov/articles/PMC12967217/)).

**[A] Alternatives:** SQANTI3 "Reference Match" = FSM **plus both ends within 50 bp**; ISM is split into 5'-fragment / 3'-fragment / internal ([Nat Methods 2024, PMC11093726](https://pmc.ncbi.nlm.nih.gov/articles/PMC11093726/)). **TUSCO** (Nat Commun 2026, [s41467-026-72089-1](https://www.nature.com/articles/s41467-026-72089-1)) is a drop-in contract: curated single-isoform endogenous genes, TP = Reference Match only, PTP = FSM_non-RM or ISM, reporting Sn, redundant precision, PDR, 1−FDR, inverse redundancy, non-redundant precision. **No head-to-head gffcompare-vs-SQANTI3 benchmark exists.**

---

## Top actionable items

1. **Normalise the junction NLL per candidate base and convert to a posterior with an explicit prior** (NanoSplicer S1.7) instead of a summed LLR — this is the published fix for the exact "longer template wins" failure you measured, and it comes with a built-in abstention gate (SIQ) matching your calibrated tie-breaker.
2. **IsoQuant's ONT rule is `min_known_count=1` vs `min_novel_count=3`**, plus ≥1 left/right support for non-FL known models. StringTie3's guided assembly has *no* annotation relaxation — its edge on degraded reads likely comes from the v3.0.2 poly(A)-tail end-anchoring, which is reproducible for you from Dorado `pt:i`.
3. **LIQA's censored risk-set denominator `f(L_r)/Σ_{l_j>L_r}θ_j`** is the only published short-vs-truncated-long likelihood; a 5'-start-conditioned version (fitting the profile with INDEGRA's κ) is genuinely unpublished territory.
4. **Report `--strict-match` and LRGASP/TUSCO redundancy + non-redundant precision alongside `=`** — the sibling collapse you already documented internally is a stated gffcompare behaviour, not a harness bug.
5. **Signal is a dead end for 5' completeness**: no cap/TSS signal classifier exists, and the ~11 nt motor loss is chemistry-level. Your 3' polyA channel being saturated is consistent with the literature.
