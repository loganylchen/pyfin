# Eventalign: Event-to-Sequence Alignment

## Overview

The `eventalign` function aligns detected nanopore events to a reference DNA/RNA sequence, creating a mapping between k-mers in the sequence and events in the raw signal. This is a key step in nanopore signal analysis, enabling downstream applications like:

- **Methylation calling**: Identifying modified bases from signal deviations
- **Basecalling QC**: Validating basecaller accuracy
- **Signal analysis**: Understanding nanopore sequencing characteristics
- **Adapter detection**: Locating adapter sequences in reads

## Installation

The eventalign functionality is included in the pyfin package. Build from source:

```bash
cd /home/logan/Projects/pyfin
pip install -e . --force-reinstall --no-deps
```

This compiles the C extensions including `_eventalign`.

## Basic Usage

```python
import numpy as np
from fin import detect_events, eventalign

# Load or generate raw signal (float32 numpy array)
raw_signal = np.random.randn(10000).astype(np.float32)

# Define reference sequence
sequence = "ACGTACGTACGTACGT"

# Align events to sequence
result = eventalign(
    raw_signal=raw_signal,
    sequence=sequence,
    is_rna=False,      # Use RNA detection parameters if True
    kmer_size=5        # K-mer size (default: 5)
)

# Access results
print(f"Detected {result['n_events']} events")
print(f"Scaling: scale={result['scaling']['scale']:.4f}, shift={result['scaling']['shift']:.2f}")

# Examine base-to-event mapping
for i, mapping in enumerate(result['base_to_event_map']):
    kmer = mapping['kmer']
    start = mapping['start']
    stop = mapping['stop']
    print(f"K-mer {kmer} (position {i}): events {start}-{stop}")
```

## Function Signature

```python
eventalign(
    raw_signal: np.ndarray,    # 1D float32 array of raw signal
    sequence: str,             # Reference DNA/RNA sequence
    is_rna: bool = False,      # Use RNA detection parameters
    kmer_size: int = 5,        # K-mer size for alignment
    model: dict = None         # Optional k-mer model (uses default if None)
) -> dict
```

## Return Value

The function returns a dictionary with:

- **`base_to_event_map`**: List of dictionaries, one per k-mer in the sequence:
  - `kmer`: The k-mer sequence (string)
  - `start`: First event index mapped to this k-mer (int, -1 if unmapped)
  - `stop`: Last event index mapped to this k-mer (int, -1 if unmapped)

- **`scaling`**: Dictionary with calibration parameters:
  - `scale`: Scale factor for signal normalization
  - `shift`: Shift offset for signal normalization

- **`n_events`**: Total number of detected events (int)

- **`n_aligned_pairs`**: Number of event-to-kmer pairs in the alignment (int)

## Algorithm Overview

The eventalign process consists of:

1. **Event Detection**: Raw signal → events (using Scrappie algorithm)
   - Adapter trimming (first 200 samples, last 10 samples)
   - T-test based segmentation
   - Returns: event mean, stdv, start index, length

2. **Scaling Estimation**: Calculate signal normalization parameters
   - Method of moments estimation
   - Matches event levels to model k-mer levels
   - Returns: scale and shift parameters

3. **Event-to-K-mer Alignment**: Map events to k-mers
   - Dynamic programming alignment (simplified in Python wrapper)
   - Creates base-to-event mapping
   - Each k-mer gets a range of events

4. **Post-processing**: Create output data structures

## Examples

### Example 1: Synthetic Signal

```python
from fin import eventalign
import numpy as np

# Generate synthetic signal (50 samples per base)
def generate_signal(sequence, samples_per_base=50):
    base_levels = {'A': 100, 'C': 95, 'G': 105, 'T': 90}
    signal = []
    for base in sequence:
        level = base_levels[base]
        signal.extend(np.random.normal(level, 2.5, samples_per_base))
    return np.array(signal, dtype=np.float32)

sequence = "ACGTACGT" * 5
raw_signal = generate_signal(sequence)

result = eventalign(raw_signal, sequence, is_rna=False, kmer_size=5)
print(f"Aligned {result['n_aligned_pairs']} pairs")
```

### Example 2: Methylation Detection

```python
# After eventalign, analyze signal deviations to detect modifications
base_to_event = result['base_to_event_map']
scaling = result['scaling']

for i, mapping in enumerate(base_to_event):
    if mapping['start'] == -1:
        continue
    
    kmer = mapping['kmer']
    # Check if k-mer contains potential modification site (e.g., CpG)
    if 'CG' in kmer:
        # Analyze events for this k-mer
        start, stop = mapping['start'], mapping['stop']
        # ... compare event levels to model expectations
```

### Example 3: Complete Pipeline

See `examples/test_eventalign.py` for a comprehensive example with:
- Synthetic signal generation
- Event detection
- Eventalign
- Visualization of results

Run it:
```bash
cd examples
python test_eventalign.py
```

## Parameters

### `kmer_size`

The k-mer size affects alignment granularity:
- **k=3**: Very coarse alignment, less accurate
- **k=5**: Standard for DNA/RNA (default)
- **k=6**: More specific but requires more events
- **k>6**: May suffer from sparse k-mer coverage

Trade-offs:
- Larger k → more specific but fewer events per k-mer
- Smaller k → less specific but more events per k-mer

### `is_rna`

Use `is_rna=True` for direct RNA sequencing:
- Different event detection parameters
- Adjusted for RNA signal characteristics
- May affect event detection sensitivity

### `model`

Custom k-mer models (optional):
```python
model = {
    0: {'level_mean': 100.0, 'level_stdv': 2.5},  # kmer rank 0
    1: {'level_mean': 95.0, 'level_stdv': 2.3},   # kmer rank 1
    # ... for all 4^k kmers
}
result = eventalign(signal, sequence, model=model)
```

Default: Uses built-in model based on `kmer_size`.

## Performance

- **Event detection**: ~1-2 ms per 10k samples (C implementation)
- **Alignment**: ~5-10 ms for 40-base sequence with 2000 samples
- **Memory**: O(n_events + n_kmers)

## Limitations

1. **Current implementation**: Simplified alignment algorithm
   - Uses uniform event distribution heuristic
   - Full adaptive banded alignment coming in future release

2. **Model support**: Default models only
   - Custom model support is basic
   - Future: Load models from file (ONT format)

3. **Accuracy**: Best for sequences > 20 bases
   - Shorter sequences may have poor alignment
   - Adapter trimming removes first ~200 samples

## Troubleshooting

### "eventalign extension not available"

The C extension wasn't compiled. Rebuild:
```bash
pip install -e . --force-reinstall
```

### "sequence too short for kmer_size"

Sequence must be ≥ kmer_size. Either:
- Increase sequence length
- Decrease `kmer_size`

### Poor alignment quality

- Check signal quality (adequate samples per base)
- Verify sequence matches signal origin
- Try different `kmer_size`
- Ensure proper `is_rna` setting

## Related Functions

- **`detect_events()`**: Event detection only
- **f5c eventalign**: Full f5c implementation (BAM input required)

## References

- **f5c**: https://github.com/hasindu2008/f5c
- **Nanopolish**: https://github.com/jts/nanopolish
- **Scrappie**: https://github.com/nanoporetech/scrappie (event detection)

## Future Enhancements

- [ ] Full adaptive banded alignment (from f5c/align.c)
- [ ] Model file loading (ONT format)
- [ ] HMM-based realignment
- [ ] GPU acceleration
- [ ] Methylation calling integration
- [ ] Multi-threading for large datasets
