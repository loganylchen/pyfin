# Adapter Trimming Implementation for Nanopore Signals

## Overview

I've successfully implemented MAD (Median Absolute Deviation) based adapter trimming in your `event_detection_simple.c` file, following the approach used in f5c and Scrappie.

## What Was Added

### 1. Statistical Helper Functions

- **`medianf()`**: Calculates the median of a float array
- **`quantilef()`**: Calculates quantiles (percentiles) of data
- **`madf()`**: Calculates Median Absolute Deviation with scaling factor (1.4826)

### 2. Core Trimming Functions

#### `trim_raw_by_mad()`
- Divides raw signal into non-overlapping chunks (default: 100 samples)
- Calculates MAD for each chunk to measure signal variability
- Identifies low-variation regions (adapters) at start and end
- Trims chunks where MAD ≤ threshold

#### `trim_and_segment_raw()`
- Applies MAD-based trimming first
- Adds fixed trimming (200 samples from start, 10 from end)
- Returns the trimmed signal range

### 3. Parameters Used (from f5c/Scrappie defaults)

```c
int trim_start = 200;      // Additional fixed trim from start
int trim_end = 10;          // Additional fixed trim from end
int varseg_chunk = 100;     // Chunk size for MAD calculation
float varseg_thresh = 0.0;  // Percentile threshold (0.0 = median)
```

### 4. Updated Functions

- **`compute_sum_sumsq()`**: Now accepts a start offset to work with trimmed ranges
- **`detect_events_simple()`**: Uses trimmed signal range (`work_start`, `work_end`)
- **`getevents_simple()`**: Applies adapter trimming before event detection

## How It Works

### Step-by-Step Process:

1. **Chunk Division**: Raw signal is divided into 100-sample chunks
2. **MAD Calculation**: For each chunk, calculate MAD = median(|x - median(x)|) × 1.4826
3. **Threshold Determination**: Use percentile of MADs as threshold
4. **Start Trimming**: Remove chunks from beginning where MAD ≤ threshold
5. **End Trimming**: Remove chunks from end where MAD ≤ threshold
6. **Fixed Trimming**: Additional 200 samples from start, 10 from end
7. **Event Detection**: Detect events only in the retained signal range
8. **Position Adjustment**: Adjust event positions to account for trimming

## Why This Works

- **Adapter regions** have **low variation** (stable, predictable signal)
- **DNA/RNA sequence** has **high variation** (different k-mers produce different levels)
- **Open pore** also has **low variation** (no molecule present)

MAD-based trimming effectively identifies and removes these low-variation regions.

## Usage Example

```c
// In your Python wrapper
event_table getevents_simple(size_t nsample, float *rawptr, int is_rna);

// The function now automatically:
// 1. Trims adapters using MAD-based method
// 2. Applies fixed trimming
// 3. Detects events on clean signal
// 4. Returns events with corrected positions
```

## Benefits

1. **Automatic adapter detection**: No manual specification needed
2. **Robust to signal variations**: Uses MAD (robust statistic)
3. **Consistent with f5c/nanopolish**: Same methodology
4. **Improved accuracy**: Events detected only on actual sequenced molecule

## Testing

Run the demo script to visualize the trimming:

```bash
python examples/test_adapter_trimming.py
```

This will generate a plot showing:
- Original signal with trim regions highlighted
- MAD per chunk (showing low variation in adapters)
- Trimmed signal (clean DNA/RNA sequence)

## Key Changes to Your Code

1. Added 5 new helper functions for statistics and trimming
2. Modified `compute_sum_sumsq()` to accept start offset
3. Updated `detect_events_simple()` to use trimmed ranges
4. Enhanced `getevents_simple()` with automatic adapter trimming
5. Added DNA parameter defaults in addition to RNA

## References

- f5c implementation: https://github.com/hasindu2008/f5c/blob/main/src/events.c
- Scrappie event detection: https://github.com/nanoporetech/scrappie
- MAD statistic: Robust measure of variability, less sensitive to outliers than standard deviation
