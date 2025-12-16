# PyFin Examples

This directory contains example scripts demonstrating various features of the PyFin package.

## Prerequisites

Build the extensions first:
```bash
cd <project-root>
pip install -e .
```

Or build extensions in-place:
```bash
python setup.py build_ext --inplace
```

## Quick Test

Test all examples at once:
```bash
python examples/test_all_examples.py
```

## Examples by Category

### Event Detection

**test_event.py** - Basic event detection with visualization
- Generates synthetic signal
- Detects events using f5c algorithm
- Visualizes raw signal and events

**test_trimming_fix.py** - Adapter trimming demonstration
- Shows MAD-based adapter detection
- Visualizes trimming regions

### Event Alignment (Eventalign)

**test_profile_hmm.py** ⭐ **RECOMMENDED**
- Full f5c Profile HMM eventalign
- Compares ABEA vs Profile HMM
- Shows HMM state distribution

**raw_signal_alignment_example.py**
- Complete pipeline: signal → events → alignment
- Uses real pore models
- Demonstrates 3-state HMM

**test_eventalign.py**
- Comprehensive eventalign visualization
- Event-to-kmer mapping
- Requires matplotlib

**debug_eventalign.py**
- Quick diagnostic tool for eventalign
- Shows detailed output for debugging

### DTW (Dynamic Time Warping)

**dtw_nanopore_example.py**
- CUDA-accelerated DTW for signal comparison
- Requires CUDA toolkit

**test_cuda_dtw.py**
- CUDA DTW benchmarking
- Requires CUDA and test data

**dtw_pairwise_example.py**
- Pairwise DTW distance calculation
- Requires CUDA

**benchmark_dtw_*.py**
- Performance benchmarking scripts
- Require CUDA and test data

### I/O and Utilities

**test_package.py**
- Test FASTA reader and utilities
- No special requirements

**interval_workflow_example.py**
- Interval/region management
- Requires BAM files

### Complete Workflows

**region_transcript_analysis_workflow.py** ⭐ **COMPREHENSIVE WORKFLOW**
- Complete pipeline for region-based transcript analysis
- Inputs: BAM, genome FASTA, transcriptome FASTA, GTF, POD5
- Features:
  1. Separate reads into isolated genomic regions
  2. Get candidate transcripts per region
  3. Eventalign each read to all candidate transcripts
  4. Compute DTW pairwise distances among reads
- Run with `--demo` for synthetic data demonstration
- Example:
  ```bash
  python region_transcript_analysis_workflow.py \
      --bam reads.bam \
      --genome genome.fa \
      --transcriptome transcripts.fa \
      --gtf annotation.gtf \
      --pod5 signals.pod5 \
      --output results/
  ```

### Comparison Tools

**compare_with_f5c.py**
- Compare PyFin eventalign with f5c
- Requires f5c installed

**compare_profile_hmm_with_f5c.py**
- Detailed Profile HMM vs f5c comparison
- Requires f5c output TSV files

### Deprecated/Old Examples

**eventalign_example.py**
- Old EventAligner class (deprecated)
- Use test_profile_hmm.py instead

**test_adapter_trimming.py**
- Standalone adapter trimming demo
- Deprecated (now integrated in event detection)

## Common Issues

### "Error: f5c extensions not available"

Build the extensions:
```bash
cd <project-root>
pip install -e .
```

### "ModuleNotFoundError: No module named 'fin'"

Install the package:
```bash
pip install -e .
```

Or add to PYTHONPATH:
```bash
export PYTHONPATH=<project-root>:$PYTHONPATH
```

### "CUDA not available"

DTW examples require CUDA toolkit. Either:
1. Install CUDA toolkit
2. Skip CUDA examples (package still works without CUDA)

### Import errors

Make sure you're using the correct import paths:
```python
# Correct
from fin._f5c._event import detect_events
from fin._f5c._eventalign import eventalign, profile_hmm_eventalign

# Wrong (old)
from fin import detect_events, eventalign
```

## Example Workflow

1. **Start simple**:
   ```bash
   python examples/test_event.py
   ```

2. **Try Profile HMM eventalign** (most important):
   ```bash
   python examples/test_profile_hmm.py
   ```

3. **Explore visualizations**:
   ```bash
   python examples/raw_signal_alignment_example.py
   ```

4. **If you have f5c installed**, compare results:
   ```bash
   python examples/compare_with_f5c.py --synthetic
   ```

## Getting Help

- Check the docstrings in each example file
- See main README.md for package documentation
- See PROFILE_HMM_README.md for eventalign details
- Run examples with `--help` flag (if supported)

## Contributing

When adding new examples:
1. Add clear docstrings
2. Include error handling for missing dependencies
3. Add to the appropriate category above
4. Test with `test_all_examples.py`
