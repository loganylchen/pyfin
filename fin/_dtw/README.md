# CUDA-Accelerated Dynamic Time Warping (DTW)

This module provides GPU-accelerated Dynamic Time Warping distance calculation using CUDA.

## Features

- **GPU Acceleration**: Uses NVIDIA CUDA for fast DTW computation
- **NumPy Integration**: Works seamlessly with NumPy arrays
- **Flexible Boundaries**: Supports open start and open end boundary conditions
- **Easy to Use**: Simple Python API with automatic type conversion

## Requirements

- NVIDIA GPU with CUDA support (Compute Capability 8.0 or higher for best performance)
- CUDA Toolkit (11.0 or higher recommended)
- NumPy

## Installation

The CUDA extension is built automatically during package installation if CUDA toolkit is detected:

```bash
pip install -e .
```

If CUDA is not available, the extension will be skipped and the package will install without GPU support.

## Usage

### Basic Example

```python
import numpy as np
from fin._dtw import dtw_distance

# Create two sequences
seq1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
seq2 = np.array([1.5, 2.5, 3.5, 4.5, 5.5], dtype=np.float32)

# Compute DTW distance
distance = dtw_distance(seq1, seq2)
print(f"DTW distance: {distance}")
```

### With Open Boundaries

```python
# Enable open start and end boundaries
distance = dtw_distance(
    seq1, seq2,
    use_open_start=True,
    use_open_end=True
)
```

### Check Availability

```python
from fin._dtw import is_available

if is_available():
    print("CUDA DTW is available!")
else:
    print("CUDA DTW is not available")
```

### Resource Cleanup

```python
from fin._dtw import cleanup

# After you're done computing DTW distances
cleanup()  # Free GPU resources
```

## API Reference

### `dtw_distance(seq1, seq2, use_open_start=False, use_open_end=False)`

Compute DTW distance between two sequences.

**Parameters:**
- `seq1` (array-like): First sequence (will be converted to float32)
- `seq2` (array-like): Second sequence (will be converted to float32)
- `use_open_start` (bool, optional): Enable open start boundary (default: False)
- `use_open_end` (bool, optional): Enable open end boundary (default: False)

**Returns:**
- `float`: DTW distance

**Raises:**
- `RuntimeError`: If CUDA extension is not available
- `ValueError`: If input sequences are invalid

### `cleanup()`

Reset CUDA device and free all GPU resources.

### `is_available()`

Check if CUDA DTW extension is available.

**Returns:**
- `bool`: True if CUDA extension is available

## Examples

### Time Series Comparison

```python
import numpy as np
from fin._dtw import dtw_distance

# Two time series with slight temporal shift
t = np.linspace(0, 2*np.pi, 100)
signal1 = np.sin(t).astype(np.float32)
signal2 = np.sin(t + 0.5).astype(np.float32)  # Phase shift

distance = dtw_distance(signal1, signal2)
print(f"Distance between shifted signals: {distance}")
```

### Batch Processing

```python
from fin._dtw import dtw_distance, cleanup
import numpy as np

# Reference sequence
reference = np.random.randn(100).astype(np.float32)

# Compare against multiple sequences
sequences = [np.random.randn(100).astype(np.float32) for _ in range(10)]

distances = []
for seq in sequences:
    dist = dtw_distance(reference, seq)
    distances.append(dist)

print(f"Mean distance: {np.mean(distances):.4f}")
print(f"Std distance: {np.std(distances):.4f}")

# Clean up GPU resources
cleanup()
```

### Different Length Sequences

```python
# DTW naturally handles different length sequences
seq1 = np.random.randn(50).astype(np.float32)
seq2 = np.random.randn(100).astype(np.float32)

distance = dtw_distance(seq1, seq2)
print(f"Distance: {distance}")
```

## Performance Notes

- Input arrays are automatically converted to float32 and made contiguous
- For best performance, pre-allocate arrays as float32 with C-contiguous layout
- The GPU achieves best performance with sequences of length 100-10000
- For very small sequences (<50 elements), CPU implementations may be faster due to overhead

## Testing

Run the test suite:

```bash
python examples/test_cuda_dtw.py
```

This will run basic functionality tests and performance benchmarks.

## Troubleshooting

### Import Error

If you get an import error:
1. Check that CUDA toolkit is installed: `nvcc --version`
2. Verify the extension was built: Look for `_cuda_dtw` in the build output
3. Check GPU availability: `nvidia-smi`

### Runtime Errors

- **"CUDA error"**: Check that your GPU has sufficient memory
- **"Invalid input"**: Ensure sequences are 1D and not empty
- **"nvcc not found"**: Install CUDA toolkit and ensure nvcc is in PATH

## Technical Details

### Architecture

The module uses the OpenDBA DTW kernel adapted for pairwise distance calculation:
- Wavefront algorithm for efficient GPU parallelization
- Shared memory optimization for cache efficiency
- Support for open boundary conditions

### Compute Capability

The default build targets Ampere architecture (compute_80). To target different architectures, modify the `--generate-code` flag in `setup.py`.

## License

This module is part of the py-fin package and follows the same license.
