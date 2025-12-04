# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**fin** is a Python bioinformatics tool for assemble transcriptome using nanopore Direct RNA-seq data. It correctes the alignments based on raw current signals.

This is a hybrid Python-C/CUDA project with GPU acceleration support. The Python package is named `fin` and the repository is at `https://github.com/loganylchen/pyfin`.



## Architecture

### Hybrid Python-C/CUDA Structure
This project uses Python as the primary interface with computationally intensive operations implemented in C and CUDA:

1. **Python Package (`fin/`)**: High-level API and I/O handling
2. **C Extension (`fin/_f5c/`)**: f5c library for signal processing and event alignment
3. **CUDA Extension (`fin/_opendba/`)**: GPU-accelerated Dynamic Time Warping (DTW) and clustering

