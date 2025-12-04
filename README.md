# py-fin

A Python package for nanopore Direct RNA-seq data analysis and transcriptome assembly.

## Overview

`py-fin` is a bioinformatics tool designed for analyzing nanopore Direct RNA sequencing data. It provides efficient I/O handling for various genomic file formats and tools for organizing and processing RNA-seq reads.

## Installation

### From Source (Development)

```bash
git clone https://github.com/loganylchen/pyfin.git
cd pyfin
pip install -e .
```

### Install with Optional Dependencies

```bash
# Install with GPU support
pip install -e ".[gpu]"

# Install with development tools
pip install -e ".[dev]"

# Install everything
pip install -e ".[all]"
```

## Features

### File Format Readers

- **FASTA Reader** (`FASTAReader`, `FASTARecord`): Read and manipulate FASTA files
- **BAM/SAM Reader** (`BamReader`): Read alignment files with pysam
- **GTF/GFF Reader** (`GTFReader`): Parse genome annotations
- **BED Reader** (`BEDReader`): Read genomic intervals
- **Signal Files**: Support for nanopore signal formats:
  - `Fast5Reader` (ONT FAST5)
  - `Slow5Reader` (SLOW5/BLOW5)
  - `Pod5Reader` (POD5)

### Read Subset Manager

The `ReadSubsetManager` organizes reads from BAM files into logical subsets based on GTF annotations:

- Group reads by transcriptomic intervals
- Identify fusion candidate reads (supplementary alignments, chimeric reads)
- Separate unannotated reads
- Lazy generation of subset BAM files
- Generate data bundles for processing

Example:
```python
from fin.io import create_subset_manager

# Create manager and analyze reads
manager = create_subset_manager(
    bam_path='reads.bam',
    gtf_path='annotation.gtf',
    fasta_path='genome.fa'
)

# Iterate over subsets
for subset in manager.iterate_subsets():
    # Get data bundle with reads and sequences
    bundle = manager.get_data_bundle(subset)

    # Process reads
    for read in bundle.reads:
        print(read['query_name'])

    # Optionally write subset to BAM
    if subset.num_reads > 0:
        manager.write_subset_bam(subset, f'{subset.subset_id}.bam')
```

## Quick Start

### Reading FASTA Files

```python
from fin.io import FASTAReader

with FASTAReader('genome.fa') as reader:
    for record in reader.iterate_records():
        print(f"{record.id}: {record.length} bp, {record.gc_content:.1f}% GC")
```

### Reading BAM Files

```python
from fin.io import BamReader

with BamReader('alignments.bam') as reader:
    # Get stats
    stats = reader.get_file_stats()
    print(f"Total reads: {stats['total_reads']}")
    print(f"Mapped: {stats['mapped_pct']:.1f}%")

    # Fetch reads in region
    for read in reader.fetch(region='chr1:1000-2000'):
        print(read['query_name'])
```

### Reading GTF Annotations

```python
from fin.io import GTFReader

with GTFReader('annotation.gtf') as reader:
    reader.parse()

    # Get gene
    gene = reader.get_gene('GENE001')
    print(f"{gene.gene_id}: {gene.num_transcripts} transcripts")

    # Get transcripts in region
    transcripts = reader.get_transcripts_in_region('chr1', 1000, 5000)
    for tx in transcripts:
        print(f"{tx.transcript_id}: {tx.num_exons} exons")
```

## Requirements

- Python >= 3.8
- See `pyproject.toml` for dependencies

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.
