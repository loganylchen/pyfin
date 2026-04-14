# PyFIN Tutorial: Nanopore Transcript Quantification

This tutorial walks through a complete analysis using PyFIN — from raw data preparation to final quantification results.

## Table of Contents

1. [Prerequisites](#1-prerequisites)
2. [Installation](#2-installation)
3. [Input Data Preparation](#3-input-data-preparation)
4. [Running Quantification](#4-running-quantification)
5. [Understanding the Output](#5-understanding-the-output)
6. [Advanced Usage](#6-advanced-usage)
7. [Docker Usage](#7-docker-usage)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Prerequisites

### Hardware

- **CPU mode**: Any modern x86_64 machine, ≥16 GB RAM recommended
- **GPU mode** (recommended): NVIDIA GPU with compute capability ≥7.0 (Volta or newer), CUDA Toolkit 12.x installed

### Software

| Tool | Version | Purpose |
|------|---------|---------|
| Python | ≥3.8 | Runtime |
| f5c | ≥1.3 | Signal-level alignment (eventalign) |
| samtools | ≥1.15 | BAM indexing and sorting |
| minimap2 | ≥2.24 | Read-to-transcript alignment (used internally via mappy) |

Install external tools with conda:

```bash
conda install -c bioconda f5c samtools minimap2
```

---

## 2. Installation

### Option A: From source (recommended for development)

```bash
git clone https://github.com/loganylchen/pyfin.git
cd pyfin

# CPU only
pip install -e .

# With GPU acceleration (requires nvcc on PATH)
pip install -e ".[gpu]"
```

The build system auto-detects your GPU's compute capability. If no GPU is found, the CUDA DTW extension is skipped and PyFIN falls back to CPU.

### Option B: Docker

```bash
# GPU image (requires nvidia-docker)
docker pull <dockerhub_user>/pyfin:gpu-latest

# CPU image
docker pull <dockerhub_user>/pyfin:cpu-latest
```

### Verify installation

```bash
# Check package
python -c "import fin; print(f'PyFIN v{fin.__version__}')"

# Check CUDA DTW
python -c "from fin._dtw import is_available; print(f'CUDA DTW: {is_available()}')"

# Check CuPy (GPU EM)
python -c "from fin.analysis.assignments import CUPY_AVAILABLE; print(f'CuPy: {CUPY_AVAILABLE}')"

# Check external tools
which f5c samtools
```

---

## 3. Input Data Preparation

PyFIN requires **four input files per sample** plus **two shared reference files**.

### Per-sample files

| File | Format | Description |
|------|--------|-------------|
| BAM | `.bam` (indexed) | Reads aligned to genome with minimap2 |
| FASTQ | `.fq` / `.fastq` | Raw basecalled reads |
| Signal | `.blow5` / `.slow5` / `.pod5` | Raw nanopore signal |

### Shared reference files

| File | Format | Description |
|------|--------|-------------|
| GTF | `.gtf` | Transcript annotations (e.g., GENCODE, Ensembl) |
| Genome FASTA | `.fa` (indexed) | Reference genome sequence |

### Step-by-step data preparation

#### 3.1 Basecall with Dorado

```bash
dorado basecaller rna004_130bps_sup@v5.1.0 pod5_dir/ \
    --emit-fastq > sample1.fastq

# Or if you already have basecalled data, skip this step
```

#### 3.2 Convert signal to BLOW5 (if starting from POD5)

```bash
# Option 1: Use BLOW5 format
blue-crab p2s pod5_dir/ -o sample1.blow5

# Option 2: Use POD5 directly (set --signal-format pod5 later)
```

#### 3.3 Align reads to genome

```bash
minimap2 -ax splice -uf -k14 --secondary=no \
    genome.fa sample1.fastq \
    | samtools sort -o sample1.bam
samtools index sample1.bam
```

> **Note**: Use `-ax splice -uf -k14` for direct RNA-seq. The `-uf` flag handles the RNA strand orientation.

#### 3.4 Prepare genome FASTA index

```bash
samtools faidx genome.fa
```

#### 3.5 Verify all files are ready

```bash
# Check all files exist
ls -lh genome.fa genome.fa.fai annotations.gtf
ls -lh sample1.bam sample1.bam.bai sample1.fastq sample1.blow5
```

Your directory should look like:

```
project/
├── genome.fa            # Reference genome
├── genome.fa.fai        # FASTA index
├── annotations.gtf      # Gene annotations (GENCODE/Ensembl)
├── sample1.bam          # Aligned reads
├── sample1.bam.bai      # BAM index
├── sample1.fastq        # Raw reads
├── sample1.blow5        # Raw signal
├── sample2.bam          # (if multiple samples)
├── sample2.bam.bai
├── sample2.fastq
└── sample2.blow5
```

---

## 4. Running Quantification

### 4.1 Single sample

```bash
fin quantify \
    --gtf annotations.gtf \
    --genome genome.fa \
    --sample "sample1:sample1.bam:sample1.fastq:sample1.blow5" \
    --output-dir ./results \
    --signal-format slow5
```

### 4.2 Multiple samples

Specify `--sample` multiple times:

```bash
fin quantify \
    --gtf annotations.gtf \
    --genome genome.fa \
    --sample "ctrl_rep1:ctrl1.bam:ctrl1.fastq:ctrl1.blow5" \
    --sample "ctrl_rep2:ctrl2.bam:ctrl2.fastq:ctrl2.blow5" \
    --sample "treat_rep1:treat1.bam:treat1.fastq:treat1.blow5" \
    --sample "treat_rep2:treat2.bam:treat2.fastq:treat2.blow5" \
    --output-dir ./results \
    --signal-format slow5
```

### 4.3 Using POD5 signal files

```bash
fin quantify \
    --gtf annotations.gtf \
    --genome genome.fa \
    --sample "sample1:sample1.bam:sample1.fastq:sample1.pod5" \
    --output-dir ./results \
    --signal-format pod5
```

### 4.4 CPU-only mode

```bash
fin quantify \
    --gtf annotations.gtf \
    --genome genome.fa \
    --sample "sample1:sample1.bam:sample1.fastq:sample1.blow5" \
    --output-dir ./results \
    --no-gpu
```

### 4.5 Custom tool paths

If `f5c` or `samtools` are not on your `$PATH`:

```bash
fin quantify \
    --gtf annotations.gtf \
    --genome genome.fa \
    --sample "sample1:sample1.bam:sample1.fastq:sample1.blow5" \
    --output-dir ./results \
    --f5c-path /path/to/f5c \
    --samtools-path /path/to/samtools
```

### 4.6 Verbose output

Add `-v` for debug-level logging:

```bash
fin quantify \
    --gtf annotations.gtf \
    --genome genome.fa \
    --sample "sample1:sample1.bam:sample1.fastq:sample1.blow5" \
    --output-dir ./results \
    -v
```

### Sample name format

The `--sample` argument follows this format:

```
--sample "NAME:BAM_PATH:FASTQ_PATH:SIGNAL_PATH"
```

| Field | Description |
|-------|-------------|
| `NAME` | Sample identifier (used in output column headers and BAM filenames) |
| `BAM_PATH` | Path to the indexed BAM file |
| `FASTQ_PATH` | Path to the FASTQ file |
| `SIGNAL_PATH` | Path to BLOW5/SLOW5/POD5 signal file |

---

## 5. Understanding the Output

After a successful run, the output directory contains:

```
results/
├── counts.tsv           # Raw abundance matrix
├── tpm.tsv              # TPM-normalized expression matrix
├── sample1.bam          # Per-sample read-to-transcript assignment BAM
├── sample1.bam.bai      # BAM index
├── sample2.bam          # (one per sample)
├── sample2.bam.bai
└── work/                # Intermediate files (can be deleted)
    ├── sample1/
    │   ├── chr1_100_5000/
    │   │   ├── ENST001/
    │   │   │   ├── candidate.fa
    │   │   │   └── aligned.bam
    │   │   └── ENST002/
    │   │       ├── candidate.fa
    │   │       └── aligned.bam
    │   └── ...
    └── sample2/
        └── ...
```

### 5.1 counts.tsv

Raw probability-weighted read counts per transcript per sample:

```
transcript_id    ctrl_rep1    ctrl_rep2    treat_rep1    treat_rep2
ENST00000456328  12.3400      11.8900      5.2100        4.9800
ENST00000450305  0.0000       0.0000       0.0000        0.0000
ENST00000488147  45.6700      48.2300      42.1000       43.5600
```

- Each row is a transcript (identified by GTF `transcript_id`)
- Values are probability-weighted counts from the EM algorithm, not simple integer read counts

### 5.2 tpm.tsv

TPM (Transcripts Per Million) normalized expression values:

```
transcript_id    ctrl_rep1    ctrl_rep2    treat_rep1    treat_rep2
ENST00000456328  5234.1200    5102.3400    2104.5600     2011.2300
ENST00000488147  18923.4500   19234.5600   17845.2300    18102.3400
```

- TPM normalizes for transcript length and sequencing depth
- Columns within each sample sum to 1,000,000
- Suitable for cross-sample comparison

### 5.3 Per-sample assignment BAMs

Each `<sample_name>.bam` contains the final read-to-transcript mapping after EM assignment:

```bash
# View assignment BAM
samtools view results/sample1.bam | head

# Check which transcripts have assigned reads
samtools idxstats results/sample1.bam

# Extract reads assigned to a specific transcript
samtools view results/sample1.bam ENST00000456328
```

The BAM references are transcript sequences (not the genome), with each read placed on its hard-assigned transcript. This is useful for:

- Visualizing read coverage per transcript
- Downstream analysis of read-level assignments
- Validating quantification results

### 5.4 Work directory

The `work/` directory contains per-candidate intermediate files from the scoring pipeline. These can be deleted after a successful run to save disk space:

```bash
rm -rf results/work/
```

---

## 6. Advanced Usage

### 6.1 Python API

For programmatic access or integration into custom pipelines:

```python
from fin.pipeline.quantify_runner import QuantifyRunner, SampleInput

samples = [
    SampleInput(
        name="ctrl_rep1",
        bam_path="ctrl1.bam",
        fastq_path="ctrl1.fastq",
        signal_path="ctrl1.blow5",
    ),
    SampleInput(
        name="treat_rep1",
        bam_path="treat1.bam",
        fastq_path="treat1.fastq",
        signal_path="treat1.blow5",
    ),
]

runner = QuantifyRunner(
    gtf_path="annotations.gtf",
    genome_fasta_path="genome.fa",
    samples=samples,
    output_dir="./results",
    signal_format="slow5",
    use_gpu=True,
    em_sigma=1.0,       # EM temperature (lower = harder assignments)
    em_beta=0.5,         # Coherence weight (0 = ignore read-to-read similarity)
    em_max_iter=1000,
    em_tol=1e-4,
)

try:
    runner.setup()
    runner.run()
finally:
    runner.cleanup()
```

### 6.2 EM parameters

The EM algorithm has two key parameters that control read assignment behavior:

| Parameter | Default | Effect |
|-----------|---------|--------|
| `em_sigma` | 1.0 | Temperature for converting distances to probabilities. **Lower** values make assignments sharper (more confident). **Higher** values make assignments softer (more uncertain). |
| `em_beta` | 0.5 | Weight for read-to-read coherence regularization. `0.0` = pure distance-based assignment. `1.0` = heavily weight read similarity (reads with similar signals should be assigned together). |

Tuning guidance:

- **High-confidence annotations** (well-annotated organism): Use default `sigma=1.0, beta=0.5`
- **Many similar isoforms** (complex alternative splicing): Try `sigma=0.5, beta=0.7` for sharper distinction
- **Noisy data** (low quality reads): Try `sigma=2.0, beta=0.3` to be more permissive

### 6.3 Using individual modules

```python
# Candidate discovery only
from fin.candidates.discovery import discover_gtf_only
from fin.io.interval_manager import GenomicInterval
from fin.io.io_gtf import GTFReader

interval = GenomicInterval(chrom="chr1", start=1000, end=50000, strand="+")
with GTFReader("annotations.gtf") as gtf:
    gtf.parse()
    candidates = discover_gtf_only(interval, "reads.bam", gtf, genome_fasta_seq)

# EM assignment only
from fin.analysis.assignments import em_with_coherence
R, hard_assignments, log_likelihoods = em_with_coherence(
    dist_read_to_tx=dist_matrix,       # (n_reads, n_transcripts)
    dist_read_to_read=dtw_matrix,      # (n_reads, n_reads)
    sigma=1.0,
    beta=0.5,
    use_gpu=True,
)
# R: soft assignment matrix (n_reads, n_transcripts)
# hard_assignments: array of length n_reads, each value is the assigned transcript index
```

---

## 7. Docker Usage

### 7.1 GPU image

```bash
docker run --gpus all --rm \
    -v /path/to/data:/data \
    <dockerhub_user>/pyfin:gpu-latest \
    fin quantify \
        --gtf /data/annotations.gtf \
        --genome /data/genome.fa \
        --sample "sample1:/data/sample1.bam:/data/sample1.fastq:/data/sample1.blow5" \
        --output-dir /data/results
```

> **Note**: `--gpus all` requires [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html).

### 7.2 CPU image

```bash
docker run --rm \
    -v /path/to/data:/data \
    <dockerhub_user>/pyfin:cpu-latest \
    fin quantify \
        --gtf /data/annotations.gtf \
        --genome /data/genome.fa \
        --sample "sample1:/data/sample1.bam:/data/sample1.fastq:/data/sample1.blow5" \
        --output-dir /data/results \
        --no-gpu
```

### 7.3 Branch-specific images

Development branch images are tagged with the branch name:

```bash
# Latest stable
docker pull <dockerhub_user>/pyfin:gpu-latest
docker pull <dockerhub_user>/pyfin:cpu-latest

# Development branch
docker pull <dockerhub_user>/pyfin:gpu-dev
docker pull <dockerhub_user>/pyfin:cpu-dev
```

---

## 8. Troubleshooting

### "Missing external tools: f5c, samtools"

Install the required tools:

```bash
conda install -c bioconda f5c samtools
```

Or point to custom paths:

```bash
fin quantify --f5c-path /custom/path/f5c --samtools-path /custom/path/samtools ...
```

### "Failed to build mappy index"

A candidate transcript has an empty sequence. This usually means the GTF annotations don't match the genome FASTA. Verify:

```bash
# Check that chromosome names match between GTF and genome
grep "^>" genome.fa | head
grep -v "^#" annotations.gtf | cut -f1 | sort -u | head
```

Common mismatch: GTF uses `chr1` while genome uses `1` (or vice versa).

### CUDA errors

```bash
# Verify CUDA is accessible
nvidia-smi
nvcc --version

# Test PyFIN CUDA
python -c "from fin._dtw import is_available; print(is_available())"

# Fall back to CPU if GPU issues persist
fin quantify --no-gpu ...
```

### Out of memory (GPU)

For large datasets, the DTW pairwise computation may exceed GPU memory. PyFIN will raise a descriptive `MemoryError`. Solutions:

1. Use `--no-gpu` to fall back to CPU
2. Process a genome region with fewer reads (split your BAM)
3. Use a GPU with more VRAM

### Out of memory (RAM)

Signal files (BLOW5/POD5) can be large. Ensure sufficient RAM:

- ~2x the signal file size is a reasonable estimate
- For 10 million reads, expect 32-64 GB RAM usage

### Empty output / no transcripts quantified

Check that:

1. The BAM has mapped reads in regions covered by GTF annotations
2. The GTF contains transcript entries (not just gene-level)
3. Chromosome naming is consistent across all files

```bash
# Check BAM has mapped reads
samtools flagstat sample1.bam

# Check GTF has transcripts
grep -c "transcript" annotations.gtf

# Check overlap: reads in regions with annotations
samtools view -c sample1.bam chr1:1000-50000
```

---

## Complete Example

End-to-end from raw data to results:

```bash
# 1. Prepare reference
samtools faidx genome.fa

# 2. Basecall (if needed)
dorado basecaller rna004_130bps_sup@v5.1.0 pod5_dir/ --emit-fastq > reads.fastq

# 3. Convert signal
blue-crab p2s pod5_dir/ -o reads.blow5

# 4. Align to genome
minimap2 -ax splice -uf -k14 --secondary=no genome.fa reads.fastq \
    | samtools sort -o reads.bam
samtools index reads.bam

# 5. Run PyFIN quantification
fin quantify \
    --gtf gencode.v44.annotation.gtf \
    --genome genome.fa \
    --sample "my_sample:reads.bam:reads.fastq:reads.blow5" \
    --output-dir ./pyfin_results \
    -v

# 6. Inspect results
head pyfin_results/counts.tsv
head pyfin_results/tpm.tsv
samtools idxstats pyfin_results/my_sample.bam | sort -k3 -nr | head

# 7. Clean up intermediate files (optional)
rm -rf pyfin_results/work/
```
