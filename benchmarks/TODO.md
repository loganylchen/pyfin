# Benchmark TODO — pyfin signal-based assembly

Goal (per `.omc/plans/ralplan-pyfin-sig-assembly.md` AC-26/AC-27, R8 opt B3):
**pyfin must strictly outperform at least one established baseline on at least one metric** for both the assembly task and the fusion-detection task, on a public direct-RNA-seq dataset.

---

## 1. Datasets

### 1.1 Primary — assembly + quantification
- **SG-NEx** (Singapore Nanopore Expression project) — recommended primary dataset.
  - URL: https://github.com/GoekeLab/sg-nex-data
  - Pick 1 cell line (e.g. **HCT116** or **K562**), one DRS run, ~10–20M reads.
  - Available: POD5/FAST5 + FASTQ + reference BAM + truth GTF (GENCODE) + spike-ins (SIRVs).
  - Why: standard reference for long-read isoform tools, has matched short-read RNA-seq for cross-validation.

### 1.2 Spike-in truth — calibration set
- **SIRV-Set 4** (Lexogen): 69 synthetic isoforms with known abundance.
- **Sequin Mix A/B** (Garvan): tiered concentration, ground-truth GTF.
- Why: SG-NEx contains SIRVs already; gives an absolute ground truth for sensitivity/precision and abundance accuracy.

### 1.3 Fusion-positive — fusion benchmark
- **Universal Human Reference RNA + fusion cell lines** — published DRS data with validated fusions:
  - **K562**: BCR–ABL1 (validated)
  - **MCF7**: BCAS4–BCAS3, ARFGEF2–SULF2
  - **NCI-H660**: TMPRSS2–ERG
- Sources to consider: SG-NEx (K562, MCF7), ENA studies cited by JAFFAL paper (PRJEB48143), LongGF supplement.
- Why: real fusion truth set with orthogonal RT-PCR / short-read confirmation.

### 1.4 Negative control — no-fusion
- **GM12878** (lymphoblastoid, no known recurrent fusions) DRS from SG-NEx.
- Why: bound false-positive rate of the fusion caller.

---

## 2. Tools to compare against

### 2.1 Assembly + quantification
| Tool | Version target | Binary / install | Notes |
|---|---|---|---|
| **Bambu** | ≥3.4 | R + Bioconductor; needs aligned BAM | Reference-guided, widely used baseline |
| **IsoQuant** | ≥3.6 | `pip install isoquant` | De-novo + reference-guided |
| **FLAIR** (optional) | latest | `pip install flair-brookslab` | Add if time permits |
| **StringTie2** (optional) | ≥2.2 | `conda install stringtie` | Long-read mode `-L` |

### 2.2 Fusion detection
| Tool | Version target | Notes |
|---|---|---|
| **JAFFAL** | latest | Long-read fusion caller, requires bpipe |
| **LongGF** | ≥0.1.2 | Compiled C++ binary |
| **Genion** (optional) | latest | If runtime budget allows |

All four are already enumerated in `benchmarks/run_benchmark.sh`; missing binaries skip gracefully.

---

## 3. Metrics

### 3.1 Assembly
- **Transcript-level sensitivity** = TP / (TP + FN) vs GENCODE truth set, matched by exon-chain.
- **Transcript-level precision** = TP / (TP + FP).
- **F1**.
- **Novel-isoform detection rate** (intron chains absent from GTF but supported by ≥3 reads).
- **SIRV recovery rate** at thresholds {0.1, 1, 10} TPM.

### 3.2 Quantification
- **Spearman / Pearson correlation** of estimated TPM vs SIRV/Sequin truth.
- **Median |log2FC|** between predicted and known abundance.
- **Resolution of overlapping isoforms** (per-transcript NRMSE on isoform groups).

### 3.3 Fusion
- **Recall** on validated fusion set.
- **Precision** = validated / total reported.
- **F1**.
- **Breakpoint accuracy**: fraction of calls within ±100 bp of validated coordinates.
- **False-positive rate** on GM12878 negative control (fusions per million reads).

### 3.4 Engineering
- **Wall-clock** end-to-end (single GPU vs CPU-only).
- **Peak RSS** memory.
- **GPU memory** peak (for `pyfin --use-gpu`).
- **Reproducibility**: 3 reruns, max relative deviation in TPM and fusion call set.

---

## 4. Benchmark Plan (workflow)

### Phase A — Environment
- [ ] Provision GPU node (≥24 GB VRAM) + CPU-only node for fallback runs.
- [ ] Install all baseline tools: Bambu (R), IsoQuant, JAFFAL, LongGF.
  - Capture exact versions in `benchmarks/versions.txt`.
- [ ] Lock pyfin commit SHA used for the run.

### Phase B — Data prep
- [ ] Download SG-NEx HCT116/K562/MCF7 DRS (POD5 or BLOW5 + FASTQ).
- [ ] Align reads with minimap2 `-ax splice -uf -k14` (matched protocol across tools).
- [ ] Index BAM, validate read counts, drop secondary/duplicate.
- [ ] Stage truth GTF (GENCODE v45) + SIRV/Sequin annotation.
- [ ] Curate validated fusion call set (BCR-ABL1, BCAS4-BCAS3, TMPRSS2-ERG, …).

### Phase C — Run pyfin
- [ ] `pyfin --bam … --gtf … --signal … --fusion --output-dir out/pyfin/` on each dataset.
- [ ] Record GPU + CPU runtimes; verify CPU/GPU result equivalence (already covered by `tests/integration/test_cpu_fallback.py` for synthetic data — repeat at scale).

### Phase D — Run baselines
- [ ] Bambu — `prepareAnnotations` + `bambu()` on the same BAM.
- [ ] IsoQuant — `isoquant.py --reference … --genedb … --bam …`.
- [ ] JAFFAL — long-read pipeline on FASTQ.
- [ ] LongGF — on the BAM.
- [ ] Each writes into `out/<tool>/`; harness already supports this layout.

### Phase E — Score
- [ ] Extend `benchmarks/compare_results.py`:
  - Parse each tool's GTF/TSV.
  - Match transcripts to truth via exon-chain hash.
  - Compute Section 3 metrics; emit `comparison.tsv` + per-metric plots (matplotlib, optional).
- [ ] Add fusion comparator: parse BEDPE / tool-specific fusion output, intersect with truth set within ±100 bp.

### Phase F — Decide pass/fail
- Per AC-27, **at least one** of the following must be true on the primary dataset:
  - pyfin F1 > Bambu F1 OR pyfin F1 > IsoQuant F1 (assembly), AND/OR
  - pyfin fusion recall > JAFFAL recall OR pyfin fusion F1 > LongGF F1.
- If neither holds: tune `score_alpha`, `prior_weight_cap`, `em_beta` on a held-out subset, then rerun.

### Phase G — Publish
- [ ] Write `benchmarks/RESULTS.md` with table of metrics, hardware, versions, commit SHA, command lines.
- [ ] Add 1–2 figures (precision-recall + fusion confusion matrix).
- [ ] Move into README "Benchmarks" section.

---

## 5. Open Decisions (need user input before Phase A)

- [ ] **Primary dataset**: SG-NEx HCT116 default — confirm or pick alternative.
- [ ] **Hardware budget**: 1 GPU (24 GB) + 32-core CPU acceptable? Otherwise tweak `max_reads_per_interval_for_dtw`.
- [ ] **Time budget**: full 4-tool sweep on 1 dataset is ~1–2 days wall-clock; OK?
- [ ] **Fusion negative control**: GM12878 OK or prefer iPSC line?

---

## 6. Out-of-scope for v1 benchmark (defer)

- Multi-sample DE analysis (Bambu/IsoQuant strength; pyfin's `quantify` mode handles known-tx case but DE is downstream).
- RNA modification benchmarks (m6A, etc.) — separate effort.
- Short-read concordance — only if a reviewer explicitly asks.

---

## 7. Acceptance for "benchmark complete"

- [ ] `benchmarks/RESULTS.md` exists with concrete numbers on ≥1 public dataset.
- [ ] At least 1 strict win recorded (assembly **or** fusion).
- [ ] `comparison.tsv` reproducible from raw outputs via `compare_results.py`.
- [ ] Hardware + tool versions + commit SHA captured.
- [ ] CI smoke (`tests/integration/test_benchmark_smoke.py`) still green.
