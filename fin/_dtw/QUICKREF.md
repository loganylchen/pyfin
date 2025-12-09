# CUDA DTW Quick Reference

## Import

```python
from fin._dtw import dtw_distance, is_available, cleanup
```

## Basic Usage

```python
import numpy as np

# Create sequences (must be 1D, float32)
seq1 = np.array([1.0, 2.0, 3.0], dtype=np.float32)
seq2 = np.array([1.5, 2.5, 3.5], dtype=np.float32)

# Compute distance
distance = dtw_distance(seq1, seq2)
```

## With Options

```python
# Open boundaries (for subsequence matching)
distance = dtw_distance(
    seq1, seq2,
    use_open_start=True,  # Allow alignment to start anywhere
    use_open_end=True     # Allow alignment to end anywhere
)
```

## Check Availability

```python
if is_available():
    print("CUDA DTW available")
else:
    print("CUDA DTW not available")
```

## Cleanup

```python
# After processing (frees GPU memory)
cleanup()
```

## Input Requirements

- **Type**: `np.float32` (will auto-convert)
- **Dimensions**: 1D arrays only
- **Length**: Any positive length (can be different)
- **Memory**: Must fit in GPU memory

## Common Patterns

### Compare Multiple Sequences

```python
reference = np.random.randn(100).astype(np.float32)
sequences = [np.random.randn(100).astype(np.float32) for _ in range(10)]

distances = [dtw_distance(reference, seq) for seq in sequences]
best_match = np.argmin(distances)
```

### Batch Processing with Cleanup

```python
results = []
for i, seq in enumerate(sequences):
    dist = dtw_distance(reference, seq)
    results.append(dist)
    
    if i % 100 == 0:  # Periodic cleanup
        cleanup()

cleanup()  # Final cleanup
```

### Handle Missing CUDA

```python
try:
    from fin._dtw import dtw_distance
    distance = dtw_distance(seq1, seq2)
except ImportError:
    # Fallback to CPU implementation
    distance = cpu_dtw(seq1, seq2)
```

## Performance Tips

1. **Pre-convert to float32**: Avoid repeated conversions
   ```python
   seq = seq.astype(np.float32)  # Do once
   ```

2. **Ensure contiguous arrays**: Better memory access
   ```python
   seq = np.ascontiguousarray(seq)
   ```

3. **Batch similar sizes**: GPU efficiency
   ```python
   # Group sequences by length before processing
   ```

4. **Call cleanup() periodically**: Free GPU memory
   ```python
   if batch_count % 1000 == 0:
       cleanup()
   ```

## Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| "RuntimeError: CUDA DTW not available" | CUDA not installed | Install CUDA Toolkit |
| "ValueError: must be 1-dimensional" | 2D/3D array | Use 1D arrays only |
| "TypeError: must be float32" | Wrong dtype | Convert with `.astype(np.float32)` |
| "ValueError: cannot be empty" | Zero length | Check input length |
| "RuntimeError: CUDA DTW failed" | GPU error | Check `nvidia-smi`, reduce size |

## Examples

See:
- `examples/test_cuda_dtw.py` - Full test suite
- `examples/dtw_nanopore_example.py` - Nanopore signal analysis
- `fin/_dtw/README.md` - Complete documentation
