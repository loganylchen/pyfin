# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**fin** is a Python bioinformatics tool for detecting RNA modifications using nanopore Direct RNA-seq data. It compares signal differences between native RNA and whole-transcriptomic in vitro transcribed products to identify RNA modifications.

This is a hybrid Python-C/CUDA project with GPU acceleration support. The Python package is named `fin` and the repository is at `https://github.com/loganylchen/pyfin`.

## Development Commands

### Installation (Development)
```bash
# Install dependencies
python -m pip install --upgrade pip setuptools wheel
pip install -r requirements.txt
pip install -e .
```

### Testing
```bash
# Run tests with coverage
pytest tests/ -v --cov=fin --cov-report=xml

# Run a single test
pytest tests/path/to/test_file.py::test_function_name -v
```

### Linting and Code Formatting
```bash
# Run flake8 linter
flake8 fin/ --count --select=E9,F63,F7,F82 --show-source --statistics
flake8 fin/ --count --exit-zero --max-complexity=10 --max-line-length=88 --statistics

# Check code formatting with black
black --check fin/

# Auto-format code with black
black fin/
```

## Architecture

### Hybrid Python-C/CUDA Structure
This project uses Python as the primary interface with computationally intensive operations implemented in C and CUDA:

1. **Python Package (`fin/`)**: High-level API and I/O handling
2. **C Extension (`fin/_f5c/`)**: f5c library for signal processing and event alignment
3. **CUDA Extension (`fin/_opendba/`)**: GPU-accelerated Dynamic Time Warping (DTW) and clustering

### Key Components

**Core Modules (`fin/core/`):**
- `signal_processor.py`: Signal processing algorithms for nanopore data
- `modification_detector.py`: RNA modification detection logic
- `eventalign.py`: Event alignment functionality (interfaces with f5c)
- `event_detection.py`: Event detection algorithms
- `dtw_gpu.py`: GPU-accelerated Dynamic Time Warping
- `f5c_wrapper.py`: Python wrapper for the f5c C library

**I/O Modules (`fin/io/`):**
- `io_signal.py`: Signal file parsing (FAST5, POD5, SLOW5 formats)
- `io_bam.py`: BAM alignment file parsing
- `io_manager.py`: Centralized I/O management
- `sequence_extractor.py`: Sequence extraction utilities

**C/CUDA Components:**
- `fin/_f5c/`: f5c C library (17 C source files) for core signal processing, event alignment, and HMM models
- `fin/_opendba/`: OpenDBA CUDA library for DTW Barycenter Averaging with GPU acceleration

### Build System

The project uses setuptools with custom extensions:

- **setup.py**: Main build configuration with two C/C++ extensions
- **f5c extension**: Pure C extension for signal processing
- **OpenDBA extension**: CUDA extension (falls back to CPU if PyTorch not available)

GPU support is optional and autodetected via PyTorch's build system. The extensions are built with `-O3` optimization and require numpy headers.

### Data Flow

1. **Input**: Nanopore sequencing data in FAST5/POD5/SLOW5 formats + BAM alignments
2. **Signal Processing**: Raw signals processed through f5c library for event detection and alignment
3. **GPU Acceleration**: DTW and clustering operations accelerated via OpenDBA CUDA kernels
4. **Output**: RNA modification detection results and optionally consensus signals

### Key Dependencies

**Scientific Computing:**
- numpy, pandas, scipy (core numerical operations)
- h5py (HDF5 file handling)

**Bioinformatics:**
- pysam (BAM file parsing)
- pod5, slow5, pyf5 (nanopore signal formats)
- ont-pyguppy-client-lib (ONT basecalling integration)

**Build Requirements:**
- numpy (headers for C extensions)
- PyTorch (optional, for CUDA build tools)
- CUDA Toolkit (optional, for GPU acceleration)

### CI/CD

GitHub Actions workflow (`.github/workflows/python-package.yml`):
- Tests on Ubuntu and macOS
- Python versions 3.8-3.11
- Runs pytest, flake8, black
- Builds and pushes Docker images on main branch

## Algorithm Pipeline

The core algorithm processes gene regions by integrating isoform identification with nanopore signal analysis. For each targeted gene (e.g., `TP53`):

### 1. Region Isolation and Input Preparation

**Input Files:**
- BAM file: Contains read alignments
- GTF file: Contains reference isoform annotations

**Processing:**
- Each gene region is processed independently
- Extract all isoform sequences from GTF annotations
- Extract all overlapping reads from BAM file

### 2. Sequence Extraction

**Annotation-Driven Sequences:**
- Parse GTF file to extract all annotated isoforms for the target gene
- Each isoform is represented as a splice variant with exon coordinates
- Transform genomic coordinates to transcript sequences

**Read-Derived Sequences:**
- Parse BAM alignments to infer read-level isoforms
- Each aligned read represents a potential isoform based on its splicing pattern
- Extract the read sequence and its alignment structure

### 3. Signal Alignment via Eventalign

**Event Alignment Process:**
- For each read, extract its nanopore raw signal (from FAST5/POD5/SLOW5)
- For each isoform sequence (both annotated and read-derived), perform event alignment using f5c
- Eventalign generates a sequence of events mapping signal to genomic positions

**Completeness Calculation:**
- For each read-isoform pair, calculate alignment completeness score
- Completeness measures how well the read's signal aligns to the isoform sequence
- High completeness indicates the read likely originated from that isoform

### 4. Pairwise DTW Comparison

**Distance Matrix Construction:**
- For all reads mapping to the region, compute pairwise DTW (Dynamic Time Warping) distances
- Use GPU-accelerated DTW via OpenDBA for performance
- Result: `N x N` distance matrix where N is the number of reads

### 5. Integration Matrix

**Combining Metrics:**
- **Matrix 1**: Pairwise DTW distance matrix (N x N)
- **Matrix 2**: Eventalign completeness scores (reads x isoforms)
- Integrate both matrices to identify:
  - Which isoforms are well-supported by read signals
  - Which reads cluster together based on signal similarity
  - Which annotated isoforms are real vs. false

**Validation Logic:**
- High completeness reads for an isoform provide strong evidence
- DTW clusters confirm signal-level similarity among isoform-supporting reads
- Annotated isoforms with sufficient high-completeness read support are validated as "real"

### 6. Read Assignment and Output

**Isoform Assignment:**
- For each read, assign to the most likely isoform based on:
  - Eventalign completeness score
  - Signal similarity (DTW distance to other reads)
- Reads without clear isoform support are marked as ambiguous or novel

**Output Formats:**

**GTF Format:**
- Validated isoforms with coordinates and metadata
- Isoform features: exon boundaries, splice junctions, transcript structure
- Support metrics: number of supporting reads, average completeness scores

**BED12 Format:**
- Read-level mapping information
- For each assigned read:
  - Chromosome, start, end coordinates
  - Block sizes and offsets (exon structure)
  - Assigned isoform ID
  - Mapping quality and completeness score

### 7. Key Features

**Novel Isoform Discovery:**
- Read-derived isoforms with strong signal support but no annotation are flagged as novel
- These represent potential alternatively spliced variants not in reference databases

**Quality Scoring:**
- Each isoform receives a confidence score based on:
  - Number and quality of supporting reads
  - Signal coherence within read cluster
  - Completeness of eventalign for all supporting reads

**Performance Optimization:**
- Regional processing allows parallelization across genes or chromosomes
- GPU acceleration for DTW reduces computation time from hours to minutes
- Memory-efficient streaming of BAM and signal files

## Important Notes

1. **GPU Support**: CUDA acceleration for DTW is optional but recommended. The pipeline will run on CPU if no GPU is available, though significantly slower for large datasets.

2. **Input Requirements**:
   - Nanopore direct RNA sequencing data (FAST5, POD5, or SLOW5 files)
   - BAM file from minimap2 alignment
   - Reference genome in FASTA format
   - Gene annotation GTF file (optional, for targeted analysis)

3. **Configuration**:
   - `min_completeness_threshold`: Minimum eventalign completeness to consider read-isoform match
   - `dtw_similarity_threshold`: Maximum DTW distance for reads to be in same isoform cluster
   - `min_read_support`: Minimum number of reads required to validate an isoform

4. **Missing Root README**: There's no README.md at the project root, only in the OpenDBA subdirectory.

5. **Empty Tests Directory**: The `fin/tests/` directory exists but contains no actual tests.

6. **No Docker Build**: Despite Docker configuration in CI, there's no Dockerfile in the repository.

## External Resources

- **OpenDBA Documentation**: Comprehensive README at `fin/_opendba/README.md` with details on GPU acceleration, nanopore data processing, and clustering algorithms
- **f5c**: C library adapted for signal processing (originally from nanopolish)
