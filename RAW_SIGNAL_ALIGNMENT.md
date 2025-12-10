# Raw Signal to Sequence Alignment

Complete guide to aligning nanopore raw signal directly to reference sequences using real pore models and f5c's 3-state HMM algorithm.

## Overview

The `eventalign` function provides an end-to-end pipeline for aligning nanopore raw signal to a reference sequence:

```
Raw Signal → Event Detection → 3-State HMM Alignment → Base-to-Event Mapping
```

## Features

### 1. Real Pore Models

Uses authentic nanopore pore models from f5c/nanopolish:

**RNA Models:**
- **RNA R9.4 5-mer** (default): `MODEL_ID_RNA_R9_NUCLEOTIDE`
  - Kmer size: 5
  - Model: `rna002_model_builtin_data`
  - 1,024 kmers (4^5)
  
- **RNA004 9-mer**: `MODEL_ID_RNA_RNA004_NUCLEOTIDE`
  - Kmer size: 9
  - Model: `rna004_model_builtin_data`
  - 262,144 kmers (4^9)

**DNA Models:**
- **DNA R9.4 5-mer** (default): `MODEL_ID_DNA_R9_NUCLEOTIDE`
  - Kmer size: 5
  - Model: `r9_4_450bps_cpg_6mer_template_model_builtin_data`
  - 1,024 kmers (4^5)

Each model contains:
- `level_mean`: Expected signal level for each kmer
- `level_stdv`: Standard deviation of signal
- `level_log_stdv`: Pre-computed log(stdv) for performance

### 2. Event Detection

MAD-based event detection with adapter trimming (from `event_detection_simple.c`):

1. **MAD Trimming**: Removes adapter regions using Median Absolute Deviation
   - Calculates median of signal
   - Computes MAD: `median(|signal - median|)`
   - Trims regions where signal deviates > threshold from median

2. **Peak Detection**: Identifies event boundaries using t-statistic
   - Computes running t-statistic between adjacent windows
   - Detects peaks where t-stat exceeds threshold
   - Uses sophisticated peak detection from f5c

3. **Event Segmentation**: Creates events with mean, stdv, start, length
   - Each event represents a detected signal level
   - Filtered for quality (removes very short/noisy events)

### 3. Three-State HMM Alignment

Full f5c HMM implementation (from `align.c`):

**States:**
```
PSR9_MATCH (0)      - Event matches kmer model
PSR9_BAD_EVENT (1)  - Noisy event to skip
PSR9_KMER_SKIP (2)  - Kmer with no event
```

**Transitions (10 types):**
```
MATCH → MATCH (same kmer)    : lp_mm_self  [stay probability]
MATCH → MATCH (next kmer)    : lp_mm_next  [move probability]
MATCH → BAD_EVENT            : lp_mb       [0.001]
MATCH → KMER_SKIP            : lp_mk       [0.0025]
BAD_EVENT → BAD_EVENT        : lp_bb       [stay probability]
BAD_EVENT → KMER_SKIP        : lp_bk       [0.0025]
BAD_EVENT → MATCH (next)     : lp_bm_next  [move from bad]
BAD_EVENT → MATCH (same)     : lp_bm_self  [0.01]
KMER_SKIP → KMER_SKIP        : lp_kk       [stay probability]
KMER_SKIP → MATCH            : lp_km       [move from skip]
```

**Dynamic Probabilities:**
```c
float events_per_base = (float)n_events / (float)n_kmers;
float p_stay = 1.0f - (1.0f / events_per_base);
float p_skip = 0.0025f;
float p_bad = 0.001f;
```

**Emission Model:**
```c
// F5C scaling model
float gp_mean = scaling.scale * model[kmer].level_mean + scaling.shift;
float gp_stdv = model[kmer].level_stdv * scaling.var;
float log_prob = log_normal_pdf(event.mean, gp_mean, gp_stdv, gp_log_stdv);
```

### 4. Soft-Clipping

Handles untrimmed adapters and low-quality regions:

**Pre-flanking** (start of read):
```c
TRANS_START_TO_CLIP = 0.5   // 50% prob to enter clipping
TRANS_CLIP_SELF = 0.9        // 90% prob to stay clipping
```

**Post-flanking** (end of read):
- Computed backwards from read end
- Allows natural termination in clipped region

**Benefit:**
- Automatically skips adapter events
- No need for perfect adapter trimming
- Events in adapter regions don't affect alignment

## Usage

### Python API

```python
from fin._f5c import eventalign
import numpy as np

# Raw signal from FAST5/POD5 file
raw_signal = np.array([...], dtype=np.float32)

# Reference sequence
sequence = "ATGCGATACGTAGCTAGCTAGCTA"

# Align for DNA
result = eventalign(raw_signal, sequence, is_rna=0, kmer_size=5)

# Align for RNA
result = eventalign(raw_signal, sequence, is_rna=1, kmer_size=5)

# Use RNA004 9-mer model
result = eventalign(raw_signal, sequence, is_rna=1, kmer_size=9)
```

### Parameters

- **raw_signal**: 1D numpy float32 array of raw signal samples
- **sequence**: Reference DNA/RNA sequence string (ACGT or ACGU)
- **model**: Optional (currently unused, uses built-in models)
- **is_rna**: `1` for RNA, `0` for DNA (default: 0)
- **kmer_size**: `5` or `9` (default: auto-selected based on model)

### Return Value

Dictionary with:
```python
{
    'base_to_event_map': [
        [event_idx1, event_idx2, ...],  # Events for kmer 0
        [event_idx3, event_idx4, ...],  # Events for kmer 1
        ...
    ],
    'scaling': {
        'scale': 1.234,
        'shift': 56.789
    },
    'n_events': 1500,
    'n_aligned_pairs': 1200
}
```

**base_to_event_map**: List of lists, one per kmer position
- Each entry contains event indices aligned to that kmer
- Empty list `[]` means kmer was skipped (no events)
- Only includes MATCH states (BAD_EVENT filtered out)

**scaling**: Piecewise linear scaling parameters
- `scale`: Multiplicative factor
- `shift`: Additive offset
- Transform: `signal_level = scale * model_level + shift`

**n_events**: Total events detected (after MAD trimming)

**n_aligned_pairs**: Number of event-kmer pairs (excludes soft-clipped)

## Examples

### Example 1: Basic RNA Alignment

```python
import numpy as np
from fin._f5c import eventalign

# RNA sequence
sequence = "AUGCGAUACGUAGCUAGCUA"

# Load raw signal (from FAST5/POD5)
raw_signal = load_raw_signal(...)  # float32 array

# Align using RNA R9.4 5-mer model
result = eventalign(raw_signal, sequence, is_rna=1, kmer_size=5)

print(f"Detected {result['n_events']} events")
print(f"Aligned {result['n_aligned_pairs']} pairs")

# Access base-to-event mapping
for kmer_idx in range(len(result['base_to_event_map'])):
    event_indices = result['base_to_event_map'][kmer_idx]
    kmer_seq = sequence[kmer_idx:kmer_idx+5]
    print(f"Kmer {kmer_idx} ({kmer_seq}): {len(event_indices)} events")
```

### Example 2: DNA with Untrimmed Adapters

```python
# DNA sequence
sequence = "ATGCGATACGTAGCTAGCTA"

# Raw signal with adapters (not trimmed)
# Soft-clipping will handle this automatically
raw_signal = load_untrimmed_signal(...)

# Align - adapters handled by soft-clipping
result = eventalign(raw_signal, sequence, is_rna=0, kmer_size=5)

# Adapter events are automatically skipped
# Only proper sequence events are in base_to_event_map
```

### Example 3: RNA004 9-mer Model

```python
# Longer RNA sequence for 9-mer
sequence = "AUGCGAUACGUAGCUAGCUAGCUAGCUGCUAGCUA"

# Use RNA004 9-mer model for better accuracy
result = eventalign(raw_signal, sequence, is_rna=1, kmer_size=9)

# 9-mer models provide finer resolution
for kmer_idx in range(len(result['base_to_event_map'])):
    event_indices = result['base_to_event_map'][kmer_idx]
    kmer_seq = sequence[kmer_idx:kmer_idx+9]
    print(f"9-mer {kmer_idx} ({kmer_seq}): {len(event_indices)} events")
```

## Algorithm Details

### Scaling Estimation

Piecewise linear scaling between events and model:

```c
simple_scalings_t estimate_scalings(
    const char *sequence,
    int32_t seq_len,
    simple_model_t *model,
    uint32_t kmer_size,
    event_table events)
{
    // 1. Match events to kmers (rough alignment)
    // 2. Compute linear regression: event_mean = scale * model_mean + shift
    // 3. Return scaling parameters
}
```

### DP Table Structure

3-state HMM requires expanded table:

```
Old (simplified):  dp[n_events][n_kmers]
New (3-state):     dp[n_events+1][n_kmers * 3]

For each kmer block:
  [block * 3 + 0] = MATCH state
  [block * 3 + 1] = BAD_EVENT state
  [block * 3 + 2] = KMER_SKIP state
```

### Traceback

```
1. Find best ending state (check all 3 states * all kmers + post-flank)
2. Walk backwards recording states:
   - MATCH: Output event-kmer pair, decrement row
   - BAD_EVENT: Record but don't output, decrement row
   - KMER_SKIP: Don't decrement row (no event consumed)
3. Reverse path to get forward alignment
4. Filter: Only output MATCH states
```

## Performance

### Memory Usage

```
Event detection: O(n_samples)
Model: O(4^k) where k is kmer_size
  - k=5: 1,024 kmers
  - k=9: 262,144 kmers
Alignment: O(n_events * n_kmers * 3)
```

### Time Complexity

```
Event detection: O(n_samples)
Scaling: O(n_events * n_kmers)
Alignment: O(n_events * n_kmers)
Total: O(n_samples + n_events * n_kmers)
```

### GPU Acceleration

GPU version available (align.cu):
- Currently has simplified 2-state algorithm
- TODO: Update to match 3-state CPU implementation
- Enable with: `#define CUDA_ENABLED`

## Comparison with F5C

| Feature | PyFin | F5C |
|---------|-------|-----|
| Event detection | ✓ MAD + t-stat | ✓ Same |
| Pore models | ✓ RNA/DNA R9.4, RNA004 | ✓ Same |
| HMM states | ✓ 3 states | ✓ 3 states |
| Transitions | ✓ 10 types, dynamic | ✓ Same |
| Emission model | ✓ F5C scaling | ✓ Same |
| Soft-clipping | ✓ Pre/post flanking | ✓ Same |
| GPU support | ⚠ Needs update | ✓ Full |

## Troubleshooting

### "Event detection failed"
- Check raw_signal is float32
- Ensure signal has reasonable values (not all zeros)
- Try different is_rna setting

### "Alignment failed"
- Sequence too short for kmer_size
- No events detected (check signal quality)
- Model mismatch (wrong is_rna setting)

### Few aligned pairs
- Normal if sequence has adapters (soft-clipping works)
- Check sequence matches expected chemistry
- Verify is_rna matches data type

### Model kmer size mismatch
- Specify correct kmer_size for model
- RNA R9.4: kmer_size=5
- RNA004: kmer_size=9
- DNA R9.4: kmer_size=5

## References

1. **F5C**: https://github.com/hasindu2008/f5c
2. **Nanopolish**: https://github.com/jts/nanopolish
3. **Paper**: Simpson et al., "Detecting DNA cytosine methylation using nanopore sequencing" (2017)
