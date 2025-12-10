# Event Alignment Module Structure

This directory contains the modular implementation of event-to-sequence alignment with soft-clipping support.

## File Organization

### Core Alignment Files

1. **`align_common.h`** - Common definitions and structures
   - Shared data structures (`simple_model_t`, `simple_scalings_t`, `simple_aligned_pair_t`)
   - Constants (`TRANS_START_TO_CLIP`, `TRANS_CLIP_SELF`)
   - HMM state definitions
   - Utility functions (`get_rank`, `get_kmer_rank`)
   - Function prototypes for alignment implementations

2. **`align.c`** - CPU implementation
   - `align_with_flanking_cpu()` - Main alignment function
   - `make_pre_flanking()` - Calculate pre-flanking probabilities
   - `make_post_flanking()` - Calculate post-flanking probabilities
   - `log_normal_pdf()` - Emission probability calculation
   - Viterbi algorithm with dynamic programming

3. **`align.cu`** - GPU (CUDA) implementation
   - `align_with_flanking_gpu()` - GPU-accelerated alignment
   - CUDA kernels for parallel processing:
     - `kernel_make_pre_flanking` - Parallel flanking probability calculation
     - `kernel_make_post_flanking` - Parallel post-flanking calculation
     - `kernel_init_dp_row` - Initialize first row of DP table
     - `kernel_fill_dp` - Fill DP table with event-kmer scores
   - Memory management for GPU (cudaMalloc/cudaFree)
   - Host-device data transfer

4. **`eventalign.c`** - Python wrapper
   - Python module interface (`PyInit__eventalign`)
   - `py_eventalign()` - Main Python-callable function
   - Event detection integration
   - Model parameter estimation
   - Automatic selection between CPU/GPU implementations

### Event Detection Files

5. **`event_detection_simple.h`** - Event detection header
   - Event and event table structures
   - Detector parameter structures
   - Function prototypes

6. **`event_detection_simple.c`** - Event detection implementation
   - `detect_events_simple()` - Main event detection
   - `trim_raw_by_mad()` - MAD-based adapter trimming
   - `compute_tstat()` - T-statistic calculation for event boundaries
   - Peak detection algorithms

## Compilation

### CPU-only Build

```bash
# Build with CPU alignment only
python setup.py build_ext --inplace
```

This compiles:
- `align.c` → CPU alignment
- `eventalign.c` → Python wrapper (uses CPU)
- `event_detection_simple.c` → Event detection

### GPU-enabled Build

```bash
# Build with CUDA support
CUDA_ENABLED=1 python setup.py build_ext --inplace
```

This compiles:
- `align.c` → CPU alignment (fallback)
- `align.cu` → GPU alignment (CUDA)
- `eventalign.c` → Python wrapper (auto-selects GPU if available)
- `event_detection_simple.c` → Event detection

The build system automatically:
1. Detects if CUDA is available
2. Compiles `.cu` files with `nvcc`
3. Links CUDA runtime libraries
4. Defines `CUDA_ENABLED` preprocessor flag

## Usage

### From Python

```python
from fin._f5c import eventalign
import numpy as np

# Load your data
raw_signal = np.load('signal.npy')  # float32 array
sequence = "ACGTACGTACGT..."  # Reference sequence
model = load_kmer_model()  # 5-mer model dictionary

# Align with automatic soft-clipping
result = eventalign(
    raw_signal=raw_signal,
    sequence=sequence,
    model=model,
    is_rna=1,
    kmer_size=5
)

# Results
print(f"Events detected: {result['n_events']}")
print(f"Events aligned: {result['n_aligned_pairs']}")
print(f"Soft-clipped: {result['n_events'] - result['n_aligned_pairs']}")
```

### CPU vs GPU Selection

The module automatically uses GPU if compiled with CUDA support:

```c
#ifdef CUDA_ENABLED
    // Use GPU implementation
    int n_pairs = align_with_flanking_gpu(...);
#else
    // Use CPU implementation
    int n_pairs = align_with_flanking_cpu(...);
#endif
```

You can force CPU-only by:
1. Not defining `CUDA_ENABLED` at compile time
2. Building without CUDA toolkit installed

## Algorithm Details

### Soft-Clipping with Flanking States

Both CPU and GPU implementations use the same algorithm:

1. **Pre-Flanking** (`pre_flank[i]`)
   - Probability of skipping first `i` events
   - Allows adapter events at start to be ignored
   - Based on transition probabilities: `TRANS_START_TO_CLIP`, `TRANS_CLIP_SELF`

2. **Post-Flanking** (`post_flank[i]`)
   - Probability that event `i` was last aligned
   - Allows trailing events (poly-A tails, adapters) to be ignored
   - Same transition model as pre-flanking

3. **Viterbi DP Table**
   - `dp[event][kmer]` = max log probability of alignment
   - Traceback to reconstruct optimal path
   - Only outputs high-confidence alignments

4. **Emission Probabilities**
   - Log-normal PDF of event level vs k-mer model
   - Scaled event levels account for drift
   - Model parameters: mean, stdv

### CPU Implementation

- Sequential processing of events
- Nested loops over events and k-mers
- Memory efficient (uses 2D arrays)
- Suitable for:
  - Single reads
  - Small batches
  - Systems without GPU

### GPU Implementation

- Parallel processing of k-mers
- Sequential processing of events (dependencies)
- Memory transferred to/from device
- Kernels launched per event
- Suitable for:
  - Large batches of reads
  - High-throughput analysis
  - Systems with CUDA GPU

## Performance Comparison

| Metric | CPU | GPU |
|--------|-----|-----|
| **Single read** | Fast (ms) | Overhead from transfers |
| **Batch (100s)** | Linear scaling | Massive parallelism |
| **Memory** | Host RAM only | Device RAM + transfers |
| **Compatibility** | All systems | CUDA-capable only |

**Recommendation:**
- Use CPU for: Single reads, interactive analysis, debugging
- Use GPU for: High-throughput pipelines, large datasets

## Extending the Module

### Adding New Alignment Algorithm

1. Create `align_new.c` or `align_new.cu`
2. Implement with signature:
   ```c
   int32_t align_with_flanking_new(
       simple_aligned_pair_t *out,
       const char *sequence,
       int32_t seq_len,
       event_table events,
       simple_model_t *model,
       uint32_t kmer_size,
       simple_scalings_t scaling);
   ```
3. Add declaration to `align_common.h`
4. Add selection logic in `eventalign.c`:
   ```c
   #ifdef USE_NEW_ALGORITHM
       int n_pairs = align_with_flanking_new(...);
   #else
       int n_pairs = align_with_flanking_cpu(...);
   #endif
   ```

### Optimizations

**CPU:**
- Banded DP (limit k-mer range per event)
- SIMD vectorization for emission calculations
- Multi-threading for batch processing

**GPU:**
- Shared memory for model parameters
- Coalesced memory access patterns
- Larger thread blocks for better occupancy
- Streams for concurrent kernel execution

## Testing

```bash
# Run alignment tests
python examples/eventalign_example.py

# Run event detection tests  
python examples/test_event.py

# Benchmark CPU vs GPU
python examples/benchmark_alignment.py
```

## Troubleshooting

### Compilation Errors

**CPU:**
- `undefined reference to align_with_flanking_cpu`
  - Ensure `align.c` is included in build
  - Check `align_common.h` is in include path

**GPU:**
- `nvcc not found`
  - Install CUDA toolkit
  - Set `PATH` to include CUDA bin directory
- `CUDA_CHECK failed`
  - Check GPU is available: `nvidia-smi`
  - Verify CUDA runtime version matches compile version

### Runtime Errors

- **Segmentation fault**
  - Check array bounds in alignment code
  - Verify event table is not empty
  - Ensure sequence length >= kmer_size

- **Poor alignment quality**
  - Check scaling parameters (`shift`, `scale`)
  - Verify k-mer model matches your data (DNA vs RNA)
  - Inspect pre/post flanking probabilities

## References

- f5c: https://github.com/hasindu2008/f5c
- Nanopolish: https://github.com/jts/nanopolish
- CUDA Programming Guide: https://docs.nvidia.com/cuda/

## License

This implementation is adapted from f5c and nanopolish, which are licensed under MIT.
