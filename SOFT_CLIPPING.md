# Soft-Clipping in Eventalign

## Overview

The eventalign module now includes **soft-clipping** functionality inspired by f5c/nanopolish. This allows the aligner to automatically handle untrimmed adapter sequences and low-quality regions at the start and end of reads.

## What is Soft-Clipping?

Soft-clipping is a probabilistic mechanism that allows the HMM aligner to **skip events** at the beginning and end of the signal that don't align well to the reference sequence. These are typically:

- **Adapter sequences** (if not removed during event detection)
- **Poly-A tails** in RNA
- **Open pore regions** (low variation signal)
- **Stall regions** at the start

## How It Works

### Pre-Flanking States (Start of Read)

The aligner maintains a probability distribution `pre_flank[i]` representing the log probability of skipping the first `i` events as background/adapter signal:

```
pre_flank[0] = log(1 - TRANS_START_TO_CLIP)  # No skipping
pre_flank[1] = log(TRANS_START_TO_CLIP) + background_prob + log(1 - TRANS_CLIP_SELF)
pre_flank[i] = log(TRANS_CLIP_SELF) + background_prob + pre_flank[i-1]
```

### Post-Flanking States (End of Read)

Similarly, `post_flank[i]` represents the probability that event `i` was the last aligned event, with all remaining events treated as background:

```
post_flank[n-1] = log(1 - TRANS_START_TO_CLIP)  # All aligned
post_flank[i] = log(TRANS_CLIP_SELF) + background_prob + post_flank[i+1]
```

### Transition Probabilities

Two key parameters control soft-clipping behavior:

- **`TRANS_START_TO_CLIP = 0.5`**: Probability of entering the clipping state (50% chance)
- **`TRANS_CLIP_SELF = 0.9`**: Probability of staying in the clipping state (90% chance)

These values mean:
- Once in the clipping state, there's a 90% chance of staying (continuing to skip events)
- Each event contributes a background probability penalty of `-3.0` (log scale)

## Benefits

### 1. Robust to Trimming Quality

Even if adapter trimming fails or is incomplete during event detection, the aligner won't be affected:

```python
from fin._f5c import eventalign
import numpy as np

# Signal with untrimmed adapter at start (first 1000 samples are adapter)
raw_signal = np.concatenate([
    np.random.randn(1000) * 2 + 100,  # Adapter: low variation
    actual_signal                      # Real signal
])

# Aligner will automatically skip the adapter events
result = eventalign(raw_signal, sequence, model, is_rna=1)
# Only real signal events are aligned to sequence
```

### 2. Better Alignment Quality

- Low-quality events at boundaries don't corrupt the alignment
- The HMM naturally finds the best alignment boundaries
- Reduces false positive alignments in noisy regions

### 3. Computational Efficiency

- Only high-confidence alignments are output
- Adapter events are implicitly handled without explicit detection
- Cleaner output with fewer spurious alignments

## Implementation Details

### Viterbi Algorithm with Flanking

The alignment uses a modified Viterbi algorithm:

1. **Initialization**: First event can align to any k-mer with pre-flanking probability
2. **Recursion**: Each cell considers:
   - Staying at same k-mer (multiple events per k-mer)
   - Moving to next k-mer (progression along sequence)
   - Pre-flanking (for first k-mer only)
3. **Termination**: Best ending considers post-flanking probabilities
4. **Traceback**: Reconstruct alignment path, excluding flanked regions

### DP Table Structure

```
           k-mer 0    k-mer 1    k-mer 2    ...    k-mer n
event 0   [score]    [score]    [score]    ...    [score]
event 1   [score]    [score]    [score]    ...    [score]
   ...       ...        ...        ...      ...       ...
event m   [score]    [score]    [score]    ...    [score]
```

Each cell stores:
- Maximum log probability of alignment up to that point
- Traceback pointer to previous cell

### Emission Probabilities

Events are scored against k-mer model using log-normal PDF:

```c
log_prob = log_normal_pdf(scaled_event_level, 
                          model_mean, 
                          model_stdv, 
                          model_log_stdv)
```

Where `scaled_event_level` accounts for drift and scaling.

## Comparison with Simple Alignment

### Before (Simple Heuristic):

```python
# Uniform distribution of events across k-mers
events_per_kmer = n_events / n_kmers
# All events assigned, including adapters
```

**Problems:**
- Adapter events assigned to real k-mers (incorrect)
- No probabilistic model
- Poor alignment quality

### After (Soft-Clipping HMM):

```python
# HMM with flanking states
# pre_flank[i] and post_flank[i] allow skipping
# Viterbi finds optimal alignment
# Only high-probability alignments output
```

**Benefits:**
- Adapter events automatically skipped
- Probabilistic scoring
- High-quality alignments only

## Usage Example

```python
import numpy as np
from fin._f5c import eventalign

# Load your data
raw_signal = np.load('signal_with_adapters.npy')  # May contain adapters
sequence = "ACGTACGTACGT..."  # Reference sequence
model = load_kmer_model()  # 5-mer model

# Align with automatic soft-clipping
result = eventalign(
    raw_signal=raw_signal,
    sequence=sequence,
    model=model,
    is_rna=1,
    kmer_size=5
)

# Only clean alignments are returned
print(f"Events detected: {result['n_events']}")
print(f"Events aligned: {result['n_aligned_pairs']}")  # May be less if adapters present
print(f"Soft-clipped: {result['n_events'] - result['n_aligned_pairs']}")

# Base-to-event map only contains aligned events
for kmer_data in result['base_to_event_map']:
    if kmer_data['start'] != -1:  # Has aligned events
        print(f"K-mer {kmer_data['kmer']}: events {kmer_data['start']}-{kmer_data['stop']}")
```

## Performance Notes

- **Time Complexity**: O(n_events × n_kmers) for the DP table
- **Space Complexity**: O(n_events × n_kmers) for DP and traceback
- **Optimization**: Banded alignment can reduce to O(n_events × band_width) if needed

## References

1. **f5c**: https://github.com/hasindu2008/f5c
2. **Nanopolish**: https://github.com/jts/nanopolish
3. **Scrappie**: https://github.com/nanoporetech/scrappie (original event detection)

The soft-clipping algorithm is adapted from f5c's eventalign module, which itself is based on nanopolish's HMM alignment approach.

## Configuration

The soft-clipping behavior can be tuned by modifying the constants in `fin/_f5c/eventalign.c`:

```c
#define TRANS_START_TO_CLIP 0.5f  // Increase to clip more aggressively
#define TRANS_CLIP_SELF 0.9f       // Increase to stay in clipping longer
```

Higher values = more aggressive clipping (more events skipped)
Lower values = more conservative clipping (fewer events skipped)

## Troubleshooting

### Too Many Events Clipped

If legitimate signal is being clipped:
- Check adapter trimming in event detection
- Verify k-mer model is appropriate for your data
- Consider adjusting `TRANS_START_TO_CLIP` to lower value (e.g., 0.3)

### Too Few Events Clipped

If adapters are still being aligned:
- Event detection may be failing to detect adapter boundaries
- Consider increasing `TRANS_START_TO_CLIP` (e.g., 0.7)
- Check that MAD-based trimming is working correctly

### Alignment Quality Issues

- Verify scaling parameters are reasonable
- Check k-mer model matches your chemistry (DNA vs RNA)
- Ensure reference sequence matches your sample
