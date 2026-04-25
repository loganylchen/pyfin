# PyFIN - Signal-Based Reference-Based Transcriptome Assembly

A Python package for nanopore Direct RNA-seq **transcriptome assembly with fusion
detection and per-transcript signal-based scoring**. Reference-based discovery
of known and novel isoforms; SA-tag fusion calling; f5c eventalign + DTW
read-to-read clustering produce three quality scores per candidate, which feed
back into EM as a prior to bias quantification toward signal-coherent isoforms.

## Pipeline Overview

```
BAM + GTF + POD5/SLOW5/BLOW5 + FASTQ + genome FASTA
                |
   +------------+--------------+
   v                           v
Phase 1                     Phase 1.5
Candidate Discovery         Fusion Detection (--fusion)
(GTF + novel isoforms)      (SA-tags -> breakpoints
   |                         -> fusion candidates)
   +------------+------------+
                v
        Phase 2: minimap2 + f5c eventalign scoring
                |
                v
        Phase 3: Read-to-read DTW (CUDA / CPU)
                |
                v
        Phase 4: Composite scoring per candidate
                  - coherence (within-cluster signal homogeneity)
                  - discrimination (LL margin vs next-best)
                  - combined = coherence^alpha * discrimination^(1-alpha)
                |
                v
        Phase 5: EM with combined-score prior (CuPy / numpy)
                |
                v
        Phase 6: Quantification + outputs
                  - assembly.gtf  (with score attributes)
                  - scores.tsv    (per-candidate metrics + TPM)
                  - fusions.bedpe (when --fusion)
```

### Phase Details

1. **Candidate Discovery** - Reference-based: GTF transcripts + novel isoforms via intron-chain extraction + 3' end clustering
1.5. **Fusion Detection** *(optional, `--fusion`)* - Parses SA tags, clusters breakpoints, builds spliced fusion candidates merged into the same `CandidateSet`
2. **External Scoring** - minimap2 + f5c eventalign produce per-read signal-to-reference distances
3. **Signal DTW** - Pairwise read-to-read signal similarity (CUDA C++ extension, CPU fallback)
4. **Composite Scoring** - Coherence + discrimination + combined geometric mean per candidate
5. **EM Assignment** - Probabilistic read-to-transcript assignment with coherence regularization and combined-score prior (CuPy GPU)
6. **Quantification** - Probability-weighted abundance + TPM, written to GTF / TSV / BEDPE

## Installation

### Basic (CPU only)

```bash
git clone https://github.com/loganylchen/pyfin.git
cd pyfin
pip install -e .
```

### With GPU acceleration

Requires CUDA Toolkit (nvcc on PATH):

```bash
# Install with CUDA DTW extension
pip install -e .

# Install CuPy for GPU-accelerated EM algorithm
pip install -e ".[gpu]"
```

The build system auto-detects your GPU's compute capability via `nvidia-smi`. If no GPU is detected, it falls back to PTX compute_70 for forward compatibility.

### External tool dependencies

```bash
# Required for the scoring pipeline
conda install -c bioconda minimap2 f5c samtools

# Or install individually
# minimap2: https://github.com/lh3/minimap2
# f5c: https://github.com/hasindu2008/f5c
# samtools: https://github.com/samtools/samtools
```

### Verify installation

```bash
# Check Python package
python -c "import fin; print(f'PyFIN v{fin.__version__}')"

# Check CUDA DTW availability
python -c "from fin._dtw import is_available; print(f'CUDA DTW: {is_available()}')"

# Check CuPy availability
python -c "from fin.analysis.assignments import CUPY_AVAILABLE; print(f'CuPy: {CUPY_AVAILABLE}')"

# Check external tools
which minimap2 f5c samtools
```

## Quick Start

### CLI — assembly only

```bash
pyfin \
    --bam reads.bam \
    --gtf annotations.gtf \
    --genome genome.fa \
    --fastq reads.fq \
    --signal signal.blow5 \
    --signal-format slow5 \
    --output-dir pyfin_out/
```

Outputs:
- `pyfin_out/assembly.gtf` — assembled transcripts (gtf + novel) with `coherence_score`, `discrimination_score`, `combined_score`, `tpm` attributes
- `pyfin_out/scores.tsv` — per-candidate metrics table

### CLI — assembly + fusion detection

```bash
pyfin \
    --bam reads.bam --gtf annotations.gtf --genome genome.fa \
    --fastq reads.fq --signal signal.blow5 \
    --output-dir pyfin_out/ \
    --fusion --min-support 2 --max-dist 500 --flank-bp 500
```

Adds `pyfin_out/fusions.bedpe` listing called fusion breakpoints with read support.

### CLI — quantify known transcripts across samples

```bash
pyfin quantify \
    --gtf annotations.gtf --genome genome.fa \
    --sample s1:s1.bam:s1.fq:s1.blow5 \
    --sample s2:s2.bam:s2.fq:s2.blow5 \
    --output-dir pyfin_quant/
```

### Python API

```python
from fin.pipeline import PipelineConfig, PipelineRunner

config = PipelineConfig(
    bam_path="reads.bam",
    gtf_path="annotations.gtf",
    genome_fasta_path="genome.fa",
    fastq_path="reads.fq",
    signal_path="signal.blow5",
    signal_format="slow5",
    work_dir="./pyfin_output",
    output_gtf="./pyfin_output/assembly.gtf",
    output_tsv="./pyfin_output/scores.tsv",
    output_bedpe="./pyfin_output/fusions.bedpe",
    fusion_enabled=True,        # SA-tag fusion calling
    use_gpu=True,
    use_prior=True,             # Apply combined_score prior to EM
    score_alpha=0.5,            # coherence vs discrimination weight
    em_sigma=1.0,
    em_beta=0.5,
)

runner = PipelineRunner(config)
runner.setup()
results = runner.run()
runner.cleanup()

for cid, qr in sorted(results.items(), key=lambda x: -x[1].abundance):
    print(f"{qr.candidate_id}\t{qr.abundance:.2f}\t{qr.combined_score:.3f}\t{qr.source}")
```

## Per-Transcript Scoring

Every candidate (gtf, novel, fusion) gets three scores in `[0, 1]` derived from
the f5c eventalign distances and read-to-read DTW clustering:

| Score | Meaning | Computed from |
|---|---|---|
| **coherence** | How tightly the reads assigned to this candidate cluster in signal space | Read-to-read DTW distances weighted by EM responsibilities |
| **discrimination** | How much better this candidate is than the next-best for its reads | Sigmoid-normalized log-likelihood gap on the eventalign distance matrix |
| **combined** | `coherence^alpha * discrimination^(1-alpha)` (`alpha = score_alpha`, default 0.5) | Geometric mean of the two |

`combined_score` is also fed back into the EM as a prior weight (clip-to-cap
normalization at `prior_weight_cap`, default 10×), biasing assignments toward
signal-coherent isoforms. Disable with `--no-prior` or `use_prior=False`.

## Fusion Detection

Enabled with `--fusion`. The pipeline:

1. Parses BAM `SA` tags to extract supplementary-alignment breakpoints.
2. Clusters breakpoints within `--max-dist` bp on both sides (default 500 bp).
3. Drops clusters with fewer than `--min-support` supporting reads (default 2).
4. Builds spliced fusion candidates by joining `±flank-bp` of genome sequence
   around each breakpoint (default 500 bp).
5. Fusion candidates flow through the same f5c → DTW → composite scoring → EM
   pipeline as regular transcripts; they appear in `assembly.gtf` with
   `transcript_source "fusion"` and in `fusions.bedpe`.

### Using individual modules

```python
# Candidate discovery only
from fin.candidates import discover_candidates
from fin.io.interval_manager import GenomicInterval
from fin.io.io_gtf import GTFReader

interval = GenomicInterval(chrom="chr1", start=1000, end=50000, strand="+")
with GTFReader("annotations.gtf") as gtf:
    gtf.parse()
    genome_seq = open("genome.fa").read()  # simplified
    candidates = discover_candidates(interval, "reads.bam", gtf, genome_seq)

# Eventalign parsing only
from fin.scoring import parse_eventalign_tsv, build_distance_matrix
scores = parse_eventalign_tsv("eventalign.tsv", candidate_lengths={"tx1": 1000})
dist = build_distance_matrix(scores, read_ids=["r1", "r2"], candidate_ids=["tx1", "tx2"])

# EM assignment only
from fin.analysis.assignments import em_with_coherence
R, assignments, lls = em_with_coherence(
    dist_read_to_tx, dist_read_to_read,
    sigma=1.0, beta=0.5, use_gpu=True
)
```

## GPU Acceleration

PyFIN uses GPU acceleration in three areas:

| Component | Technology | Fallback |
|-----------|-----------|----------|
| DTW pairwise distances | CUDA C++ extension (`fin._dtw`) | scipy/numpy CPU |
| EM algorithm (matrix ops) | CuPy (`cupy.matmul`, `cupy.exp`) | numpy CPU |
| DTW variable-length pairs | CUDA `dtw_distance()` per pair | numpy row-vectorized |

### GPU memory safety

The DTW pairwise batch function checks available GPU memory before allocation. For large datasets (>10,000 reads), it will raise a `MemoryError` with a descriptive message rather than crashing with a CUDA OOM error.

### Disabling GPU

```python
# Pipeline level
config = PipelineConfig(..., use_gpu=False)

# EM algorithm level
R, assignments, lls = em_with_coherence(
    dist_read_to_tx, dist_read_to_read,
    use_gpu=False  # Force numpy CPU
)
```

## Module Structure

```
fin/
  candidates/               # Transcript candidate discovery
    dataclasses.py          # IntronChain, TranscriptCandidate, CandidateSet
    intron_chains.py        # CIGAR -> intron chain, 3' clustering
    discovery.py            # Full candidate discovery per interval
  scoring/                   # External tool scoring
    external_tools.py       # minimap2/f5c subprocess wrappers
    eventalign_parser.py    # f5c TSV -> distance matrix
    signal_dtw.py           # Signal extraction + pairwise DTW
  pipeline/                  # Orchestration
    config.py               # PipelineConfig dataclass
    runner.py               # PipelineRunner (Phase 1-5)
  analysis/                  # Statistical methods
    assignments.py          # EM with coherence (CuPy GPU)
    quantification.py       # Probability-weighted abundance
    clustering.py           # 3' end position clustering
  io/                        # I/O for bioinformatics formats
    io_bam.py               # BAM reader (pysam)
    io_gtf.py               # GTF reader
    io_pod5.py              # POD5 signal reader
    io_slow5.py             # SLOW5/BLOW5 signal reader
    io_fasta.py             # FASTA reader
    io_fastq.py             # FASTQ reader
    interval_manager.py     # Genomic interval management
  _dtw/                      # CUDA DTW extension
    __init__.py             # Python API (dtw_distance, dtw_pairwise)
    dtw_api.cpp             # C++/CUDA bindings
    dtw.hpp                 # CUDA DTW kernel
  utils/                     # Utilities
    log_config.py           # Logging configuration
    sequences.py            # Sequence utilities
```

## Configuration Reference

### PipelineConfig parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `bam_path` | str | required | Path to BAM file |
| `gtf_path` | str | None | Path to GTF annotation file |
| `genome_fasta_path` | str | required | Path to genome FASTA |
| `fastq_path` | str | required | Path to reads FASTQ |
| `signal_path` | str | required | Path to signal file (SLOW5/BLOW5/POD5) |
| `signal_format` | str | "slow5" | Signal format: "slow5" or "pod5" |
| `work_dir` | str | "./pyfin_work" | Working directory for intermediate files |
| `output_gtf` | str | None | Output assembly GTF path |
| `output_tsv` | str | None | Output per-candidate scoring TSV path |
| `output_bedpe` | str | None | Output fusion BEDPE path (with `fusion_enabled`) |
| `three_prime_threshold` | int | 24 | 3' end clustering distance (bp) |
| `em_sigma` | float | 1.0 | EM temperature (lower = harder assignments) |
| `em_beta` | float | 0.5 | Coherence weight (0 = ignore read similarity) |
| `em_max_iter` | int | 1000 | Maximum EM iterations |
| `em_tol` | float | 1e-4 | EM convergence tolerance |
| `score_alpha` | float | 0.5 | Coherence vs discrimination weight in `combined_score` |
| `prior_weight_cap` | float | 10.0 | Max multiplicative boost from `combined_score` prior |
| `use_prior` | bool | True | Feed `combined_score` into EM as prior weight |
| `use_gpu` | bool | True | Enable GPU acceleration (DTW + EM) |
| `max_reads_per_interval_for_dtw` | int | 2000 | DTW subsampling cap per interval |
| `fusion_enabled` | bool | False | Enable SA-tag fusion detection |
| `fusion_min_support` | int | 2 | Min supporting reads per fusion breakpoint |
| `fusion_max_dist` | int | 500 | Breakpoint clustering distance (bp) |
| `fusion_flank_bp` | int | 500 | Flank on each side of breakpoint for fusion sequence |
| `f5c_path` | str | "f5c" | Path to f5c binary |
| `samtools_path` | str | "samtools" | Path to samtools binary |
| `max_reads` | int | None | Limit total reads processed |

## Testing

### Run all unit tests

```bash
# With the correct Python environment (needs pysam, numpy, scipy)
pytest tests/unit/ -v

# Run specific test modules
pytest tests/unit/test_intron_chains.py -v      # Intron chain extraction
pytest tests/unit/test_eventalign_parser.py -v   # Eventalign TSV parsing
pytest tests/unit/test_quantification.py -v      # Quantification
pytest tests/unit/test_assignments.py -v         # EM algorithm
pytest tests/unit/test_dtw.py -v                 # DTW (CUDA tests skipped without GPU)
```

### Test with GPU

```bash
# DTW tests require CUDA GPU
pytest tests/unit/test_dtw.py -v  # CUDA tests auto-skip if no GPU

# Verify GPU is being used
python -c "
from fin._dtw import is_available
from fin.analysis.assignments import CUPY_AVAILABLE
print(f'CUDA DTW: {is_available()}')
print(f'CuPy EM:  {CUPY_AVAILABLE}')
"
```

### Test data

Test data is located in `tests/testdata/`:
- `RNA004.test.pod5` - POD5 signal file
- `test.fa` - Reference FASTA
- `RNA004.test.bam` - Aligned reads

## Development

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run tests with coverage
pytest tests/unit/ --cov=fin --cov-report=term-missing

# Lint
ruff check fin/
```

## Benchmarking

A benchmark harness is provided in `benchmarks/` to compare pyfin against
established long-read isoform tools (Bambu, IsoQuant) and fusion callers
(JAFFAL, LongGF). Missing tools skip gracefully.

```bash
bash benchmarks/run_benchmark.sh \
    --dataset /path/to/data.bam \
    --tools pyfin,bambu,isoquant,jaffal,longgf \
    --output-dir out/
```

See `benchmarks/README.md` for the harness and `benchmarks/TODO.md` for the
recommended datasets (SG-NEx, SIRV/Sequin, fusion-positive cell lines), the
metric set (sensitivity/precision/F1, abundance correlation, fusion recall +
breakpoint accuracy), and acceptance criteria.

## License

MIT License
