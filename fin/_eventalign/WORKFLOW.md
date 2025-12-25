# Event Alignment Workflow Guide

## Overview

Event alignment requires **two** critical data sources for each read:

1. **Raw signal data** - From POD5/FAST5/SLOW5 files (electrical current measurements)
2. **Basecalled read sequence** - From FASTQ/BAM files (nucleotide sequence)

**Important**: POD5 files do NOT contain basecalled sequences. You must provide them separately.

## Quick Start

### Option 1: Using FASTQ (Recommended for new users)

If you have basecalled FASTQ files from Guppy, Dorado, or another basecaller:

```python
from fin.io.io_fastq import FastqReader
from fin.io.io_pod5 import Pod5Reader
from fin._eventalign import run_eventalign, MODEL_RNA002

# Pair reads from FASTQ with signals from POD5
# See: examples/load_from_fastq_pod5.py
```

### Option 2: Using BAM (Recommended for aligned reads)

If you have aligned BAM files:

```python
from fin.io.io_bam import BamReader
from fin.io.io_pod5 import Pod5Reader
from fin._eventalign import run_eventalign, MODEL_RNA002

# Pair reads from BAM with signals from POD5
# See: examples/load_from_bam_pod5.py
```

## Data Flow Diagram

```
┌─────────────────┐     ┌─────────────────┐
│   POD5 File     │     │  FASTQ/BAM File │
│ (Raw Signal)    │     │ (Sequences)     │
└────────┬────────┘     └────────┬────────┘
         │                       │
         │ 1. Load signal        │ 2. Load sequence
         │                       │
         ▼                       ▼
┌─────────────────────────────────────┐
│     Pair by read_id                 │
│  (signal + sequence for each read)  │
└────────┬────────────────────────────┘
         │
         │ 3. Run eventalign
         ▼
┌─────────────────────────────────────┐
│     Aligned Events                  │
│  (event-to-reference mapping)       │
└─────────────────────────────────────┘
```

## Why Both Files Are Needed

### POD5 Files
- **Contain**: Raw electrical signal (picoamperes vs time)
- **Do NOT contain**: Basecalled nucleotide sequences
- **Purpose**: High-fidelity storage of nanopore measurements

### FASTQ/BAM Files
- **Contain**: Basecalled nucleotide sequences (A, C, G, T)
- **Purpose**: Sequence interpretation of the signal

**The alignment algorithm matches patterns in the signal to k-mers in the sequence.** Without the correct sequence, alignment is impossible because the algorithm doesn't know what pattern to look for.

## Common Workflows

### Workflow 1: Basecalling with Guppy

```bash
# 1. Basecall POD5 to FASTQ
guppy_basecaller \
  -i /path/to/pod5_dir \
  -s /path/to/output \
  --flowcell FLO-MIN106 \
  --kit SQK-RNA002 \
  --cpu_threads_per_caller 8

# 2. Run event alignment using Python
python examples/load_from_fastq_pod5.py
```

### Workflow 2: Alignment with Minimap2

```bash
# 1. Basecall (if not already done)
guppy_basecaller -i pod5/ -s output/ --flowcell FLO-MIN106 --kit SQK-RNA002

# 2. Align to reference
minimap2 -ax map-ont reference.fa output/basecalled.fastq | \
  samtools sort -o aligned.bam

# 3. Index BAM
samtools index aligned.bam

# 4. Run event alignment using Python
python examples/load_from_bam_pod5.py
```

## Troubleshooting

### Alignment Returns "no_alignment"

**Symptom**: Events are detected but alignment returns 0 pairs

**Root Cause**: The read sequence doesn't match the signal

**Check**:
1. Are read IDs matching between FASTQ/BAM and POD5?
2. Was the FASTQ generated from the same POD5 file?
3. Is the read sequence the basecalled result of THIS specific signal?

**Solution**: Verify your data sources
```python
# Check if read IDs match
from fin.io.io_pod5 import Pod5Reader
from fin.io.io_fastq import FastqReader

with Pod5Reader("signal.pod5") as pod5:
    pod5_ids = set(pod5.read_ids)

with FastqReader("sequences.fastq") as fastq:
    fastq_ids = set(read_id for read_id, _, _ in fastq.iter_reads())

print(f"POD5 read IDs: {len(pod5_ids)}")
print(f"FASTQ read IDs: {len(fastq_ids)}")
print(f"Matching: {len(pod5_ids & fastq_ids)}")
```

### Poor Quality Alignments

**Symptom**: Alignment succeeds but avg_log_emission is low (< -5.0)

**Possible causes**:
1. Poor signal quality (low variance, extreme drift)
2. Wrong pore model (RNA002 vs RNA004)
3. Degraded RNA or DNA sample
4. Incorrect basecalling parameters

**Debug**: Use the debug script
```bash
python examples/debug_alignment.py
```

## File Format Summary

| Format | Contains | Use For |
|--------|----------|---------|
| POD5 | Raw signal only | Signal extraction |
| FAST5 | Raw signal (+ optional sequence) | Legacy format |
| FASTQ | Basecalled sequences | Simple sequence storage |
| BAM | Aligned sequences | Mapping + sequences |

## API Reference

### run_eventalign()

```python
result = run_eventalign(
    read_ids,      # List[str]: Read identifiers
    read_seqs,     # List[str]: Basecalled sequences (A,C,G,T)
    ref_seqs,      # List[str]: Reference sequences
    ref_names,     # List[str]: Reference names
    ref_lens,      # List[int]: Reference lengths
    signals,       # List[np.ndarray]: Raw signals (float32, pA)
    sample_rates,  # List[float]: Sample rates (Hz)
    model_id=MODEL_RNA002,
)
```

**Returns**: Dictionary with keys
- `summary`: Processing summary
- `scalings`: Normalization parameters for each read
- `events`: Detected events for each read
- `full`: Full alignment for each (read, reference) pair
- `mapping`: Base-to-event mapping for each (read, reference) pair

## Example Scripts

- `examples/load_from_fastq_pod5.py` - Load from FASTQ + POD5
- `examples/load_from_bam_pod5.py` - Load from BAM + POD5
- `examples/debug_alignment.py` - Diagnose alignment failures
- `examples/run_eventalign_example.py` - Basic usage with synthetic data
