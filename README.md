# PyFIN - Nanopore Direct RNA-seq Analysis Pipeline

A Python package for nanopore Direct RNA-seq transcript discovery, assignment, and quantification. Uses external tools (minimap2, f5c) for proven alignment accuracy and focuses on multi-evidence probabilistic transcript assignment.

## Pipeline Overview

```
BAM + GTF + POD5/SLOW5 + FASTQ + genome FASTA
         |
   +-----------+
   v           v
Phase 1     Signal
Candidate   Loading
Discovery     |
   |          |
   v          |
Phase 2       |
minimap2 +    |
f5c scoring   |
   |          |
   +----------+
   v          v
Phase 3: Signal DTW (read-to-read)
   |
   v
Phase 4: EM Assignment (eventalign + DTW coherence)
   |
   v
Phase 5: Quantification (probability-weighted)
```

### Phase Details

1. **Candidate Discovery** - Identifies transcript candidates from GTF annotations and novel isoforms via intron chain extraction + 3' end clustering
2. **External Scoring** - Aligns reads to candidates with minimap2, scores with f5c eventalign
3. **Signal DTW** - Computes pairwise read-to-read signal similarity using GPU-accelerated DTW
4. **EM Assignment** - Probabilistic read-to-transcript assignment with coherence regularization (CuPy GPU-accelerated)
5. **Quantification** - Probability-weighted transcript abundance estimation

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
    use_gpu=True,        # Enable GPU for DTW + EM
    em_sigma=1.0,        # EM temperature
    em_beta=0.5,         # Coherence weight
)

runner = PipelineRunner(config)
runner.setup()
results = runner.run()
runner.cleanup()

# Print results
for cid, qr in sorted(results.items(), key=lambda x: -x[1].abundance):
    print(f"{qr.candidate_id}\t{qr.abundance:.2f}\t{qr.confidence:.3f}\t{qr.source}")
```

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
| `three_prime_threshold` | int | 24 | 3' end clustering distance (bp) |
| `em_sigma` | float | 1.0 | EM temperature (lower = harder assignments) |
| `em_beta` | float | 0.5 | Coherence weight (0 = ignore read similarity) |
| `em_max_iter` | int | 1000 | Maximum EM iterations |
| `use_gpu` | bool | True | Enable GPU acceleration |
| `max_reads` | int | None | Limit number of reads processed |

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

## License

MIT License
