# FIN - RNA Isoform Detection

A Python tool for detecting RNA modifications using nanopore Direct RNA-seq data by comparing signal differences between native RNA and whole-transcriptomic in vitro transcribed products.

## Installation

### Development Installation

```bash
# Clone the repository
git clone https://github.com/loganylchen/pyfin.git
cd pyfin

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

### Docker Images

FIN provides pre-built Docker images for easy deployment:

**CPU-only version (recommended for most users):**
```bash
docker pull loganylchen/pyfin:latest
docker run -it --rm -v $(pwd):/data loganylchen/pyfin:latest python3 -m fin.cli --help
```

**GPU-accelerated version (for high-performance computing):**
```bash
docker pull loganylchen/pyfin:latest-gpu
docker run -it --rm --gpus all -v $(pwd):/data loganylchen/pyfin:latest-gpu python3 -m fin.cli --help
```

## Usage

```bash
# Basic usage
FIN.py --input <input.bam> --reference <reference.fa> --output <output_dir>

# See all options
FIN.py --help
```

## Features

- **RNA modification detection** using nanopore Direct RNA-seq
- **GPU-accelerated DTW** (Dynamic Time Warping) via OpenDBA
- **Event alignment** using f5c for signal processing
- **Multiple input formats**: FAST5, POD5, SLOW5
- **Flexible deployment**: Native Python or Docker

## Architecture

FIN is a hybrid Python-C/CUDA project:

- **Python layer**: High-level API and I/O handling (`fin/`)
- **f5c C extension**: Event detection and alignment (`fin/_f5c/`)
- **OpenDBA CUDA extension**: GPU-accelerated DTW (`fin/_opendba/`)

## Requirements

- Python 3.7+
- NumPy, Pandas, SciPy
- Bioinformatics libraries: pysam, h5py
- Nanopore signal libraries: pod5, slow5, pyf5
- CUDA Toolkit (optional, for GPU acceleration)

## License

MIT License

## Contributing

Contributions are welcome! Please see the project repository for guidelines.
