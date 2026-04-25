#!/usr/bin/env python3
"""
Benchmark comparison: CUDA DTW vs dtaidistance vs fastdtw vs pure Python vs Numba

This script compares the performance and accuracy of multiple DTW implementations:
1. fin._dtw (CUDA GPU-accelerated)
2. dtaidistance (optimized CPU with C extensions)
3. fastdtw (approximate DTW algorithm)
4. Pure Python (naive implementation)
5. Numba JIT (JIT-compiled Python)

Requirements:
    pip install dtaidistance fastdtw numba
"""

import numpy as np
import time
from typing import List, Tuple
import sys


# =============================================================================
# Pure Python Implementation
# =============================================================================

def dtw_python(seq1: np.ndarray, seq2: np.ndarray) -> float:
    """
    Pure Python DTW implementation (naive, unoptimized).
    
    This is the reference implementation - slow but easy to understand.
    Uses squared Euclidean distance to match CUDA implementation.
    """
    n, m = len(seq1), len(seq2)
    
    # Initialize cost matrix with infinity
    dtw_matrix = np.full((n + 1, m + 1), np.inf)
    dtw_matrix[0, 0] = 0.0
    
    # Fill the cost matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # Use squared difference to match CUDA implementation
            diff = seq1[i-1] - seq2[j-1]
            cost = diff * diff
            dtw_matrix[i, j] = cost + min(
                dtw_matrix[i-1, j],      # insertion
                dtw_matrix[i, j-1],      # deletion
                dtw_matrix[i-1, j-1]     # match
            )
    
    return dtw_matrix[n, m]


# =============================================================================
# Numba-accelerated Implementation
# =============================================================================

try:
    from numba import jit
    
    @jit(nopython=True)
    def dtw_numba(seq1: np.ndarray, seq2: np.ndarray) -> float:
        """
        Numba JIT-compiled DTW implementation.
        
        Uses @jit decorator to compile to machine code at runtime.
        Uses squared Euclidean distance to match CUDA implementation.
        """
        n, m = len(seq1), len(seq2)
        
        # Initialize cost matrix
        dtw_matrix = np.full((n + 1, m + 1), np.inf)
        dtw_matrix[0, 0] = 0.0
        
        # Fill the cost matrix
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                # Use squared difference to match CUDA implementation
                diff = seq1[i-1] - seq2[j-1]
                cost = diff * diff
                
                # Find minimum of three predecessors
                min_prev = dtw_matrix[i-1, j-1]
                if dtw_matrix[i-1, j] < min_prev:
                    min_prev = dtw_matrix[i-1, j]
                if dtw_matrix[i, j-1] < min_prev:
                    min_prev = dtw_matrix[i, j-1]
                
                dtw_matrix[i, j] = cost + min_prev
        
        return dtw_matrix[n, m]
    
    NUMBA_AVAILABLE = True
    print("✓ Numba available")
except ImportError:
    NUMBA_AVAILABLE = False
    print("✗ Numba not available (install: pip install numba)")


# =============================================================================
# External Library Implementations
# =============================================================================

# Try to import all DTW implementations
implementations = {}

# 1. Pure Python (always available)
implementations['Pure Python'] = dtw_python
print("✓ Pure Python DTW available")

# 2. Numba JIT
if NUMBA_AVAILABLE:
    implementations['Numba JIT'] = dtw_numba

# 3. CUDA DTW (our implementation)
try:
    from fin._dtw import dtw_distance as cuda_dtw, is_available as cuda_available, cleanup
    if cuda_available():
        implementations['CUDA (GPU)'] = cuda_dtw
        print("✓ CUDA DTW available")
    else:
        print("✗ CUDA DTW not available")
except ImportError as e:
    print(f"✗ CUDA DTW import failed: {e}")

# 4. dtaidistance
try:
    from dtaidistance import dtw as dtaidistance_dtw
    implementations['dtaidistance (C)'] = lambda x, y, **kwargs: dtaidistance_dtw.distance(x, y)
    print("✓ dtaidistance available")
except ImportError:
    print("✗ dtaidistance not available (install: pip install dtaidistance)")

# 5. fastdtw
try:
    from fastdtw import fastdtw
    from scipy.spatial.distance import euclidean
    
    def fastdtw_wrapper(x, y, **kwargs):
        # Ensure arrays are 1D and reshape if needed
        x_1d = np.asarray(x).ravel()
        y_1d = np.asarray(y).ravel()
        # fastdtw expects each element to be a sequence, so reshape to (n, 1)
        x_2d = x_1d.reshape(-1, 1)
        y_2d = y_1d.reshape(-1, 1)
        return fastdtw(x_2d, y_2d, dist=euclidean)[0]
    
    implementations['fastdtw (approx)'] = fastdtw_wrapper
    print("✓ fastdtw available")
except ImportError:
    print("✗ fastdtw not available (install: pip install fastdtw)")

print()


def benchmark_dtw(
    impl_name: str,
    dtw_func,
    seq1: np.ndarray,
    seq2: np.ndarray,
    n_runs: int = 10,
    warmup: int = 2
) -> Tuple[float, float]:
    """
    Benchmark a DTW implementation.
    
    Returns:
        Tuple of (distance, avg_time_ms)
    """
    # Warm-up runs
    for _ in range(warmup):
        try:
            _ = dtw_func(seq1, seq2)
        except:
            pass
    
    # Actual benchmark
    times = []
    distance = None
    
    for _ in range(n_runs):
        start = time.perf_counter()
        try:
            distance = dtw_func(seq1, seq2)
        except Exception as e:
            print(f"  Error in {impl_name}: {e}")
            return None, None
        end = time.perf_counter()
        times.append((end - start) * 1000)  # Convert to ms
    
    avg_time = np.mean(times)
    std_time = np.std(times)
    
    return distance, avg_time, std_time


def run_benchmark_suite(sequence_lengths: List[int], n_runs: int = 10):
    """
    Run benchmarks for different sequence lengths.
    
    For very slow implementations (Pure Python), reduce number of runs
    and skip very long sequences.
    """
    print("="*80)
    print("DTW Implementation Benchmark Suite")
    print("="*80)
    print(f"Number of runs per test: {n_runs}")
    print(f"Implementations available: {len(implementations)}")
    print()
    
    # Warn about Pure Python performance
    if 'Pure Python' in implementations and max(sequence_lengths) > 500:
        print("⚠ WARNING: Pure Python implementation is very slow for long sequences")
        print("  Consider reducing sequence lengths or skipping Pure Python for large tests\n")
    
    if len(implementations) == 0:
        print("ERROR: No DTW implementations available!")
        print("Please install at least one: pip install dtaidistance fastdtw")
        return
    
    results = {name: [] for name in implementations.keys()}
    
    for length in sequence_lengths:
        print(f"\n{'='*80}")
        print(f"Sequence Length: {length}")
        print(f"{'='*80}")
        
        # Generate random sequences
        np.random.seed(42)  # For reproducibility
        seq1 = np.random.randn(length).astype(np.float32)
        seq2 = np.random.randn(length).astype(np.float32)
        
        length_results = {}
        
        for impl_name, dtw_func in implementations.items():
            # Skip Pure Python for very long sequences (too slow)
            if impl_name == 'Pure Python' and length > 500:
                print(f"\nSkipping {impl_name} (too slow for length {length})...")
                results[impl_name].append((length, None, None))
                continue
            
            # Reduce runs for Pure Python to save time
            test_runs = max(3, n_runs // 3) if impl_name == 'Pure Python' and length > 200 else n_runs
            
            print(f"\nTesting {impl_name}...")
            distance, avg_time, std_time = benchmark_dtw(
                impl_name, dtw_func, seq1, seq2, n_runs=test_runs
            )
            
            if distance is not None:
                print(f"  Distance: {distance:.4f}")
                print(f"  Avg time: {avg_time:.2f} ± {std_time:.2f} ms")
                length_results[impl_name] = {
                    'distance': distance,
                    'time': avg_time,
                    'std': std_time
                }
                results[impl_name].append((length, avg_time, distance))
            else:
                print("  Benchmark failed")
                results[impl_name].append((length, None, None))
        
        # Print comparison
        if len(length_results) > 1:
            print(f"\n{'-'*80}")
            print("Comparison:")
            
            # Find fastest
            fastest = min(length_results.items(), key=lambda x: x[1]['time'])
            print(f"  Fastest: {fastest[0]} ({fastest[1]['time']:.2f} ms)")
            
            # Compare distances
            distances = [r['distance'] for r in length_results.values()]
            if len(set([f"{d:.4f}" for d in distances])) == 1:
                print("  Distance agreement: ✓ All implementations agree")
            else:
                print("  Distance agreement: ✗ Implementations differ")
                for name, res in length_results.items():
                    print(f"    {name}: {res['distance']:.4f}")
            
            # Speedup comparison
            if 'CUDA (GPU)' in length_results and len(length_results) > 1:
                cuda_time = length_results['CUDA (GPU)']['time']
                print("\n  Speedup vs CUDA:")
                for name, res in length_results.items():
                    if name != 'CUDA (GPU)':
                        speedup = res['time'] / cuda_time
                        print(f"    {name}: {speedup:.2f}x slower")
    
    # Summary
    print(f"\n\n{'='*80}")
    print("BENCHMARK SUMMARY")
    print(f"{'='*80}\n")
    
    for impl_name in implementations.keys():
        print(f"{impl_name}:")
        impl_results = results[impl_name]
        
        if all(r[1] is not None for r in impl_results):
            print(f"  {'Length':<10} {'Time (ms)':<15} {'Distance':<15}")
            print(f"  {'-'*40}")
            for length, time_ms, distance in impl_results:
                print(f"  {length:<10} {time_ms:<15.2f} {distance:<15.4f}")
        else:
            print("  Some benchmarks failed")
        print()
    
    # Overall winner
    print(f"{'='*80}")
    print("OVERALL PERFORMANCE")
    print(f"{'='*80}")
    
    avg_times = {}
    for impl_name, impl_results in results.items():
        valid_times = [t for _, t, _ in impl_results if t is not None]
        if valid_times:
            avg_times[impl_name] = np.mean(valid_times)
    
    if avg_times:
        sorted_impls = sorted(avg_times.items(), key=lambda x: x[1])
        print("\nRanking (by average time across all sequence lengths):")
        for i, (name, avg_time) in enumerate(sorted_impls, 1):
            print(f"  {i}. {name:<25} {avg_time:>10.2f} ms")
        
        if len(sorted_impls) > 1:
            fastest_time = sorted_impls[0][1]
            print(f"\nSpeedup factors (relative to {sorted_impls[0][0]}):")
            for name, avg_time in sorted_impls[1:]:
                speedup = avg_time / fastest_time
                print(f"  {name:<25} {speedup:>10.2f}x slower")
        
        # Explain CUDA performance characteristics
        if 'CUDA (GPU)' in avg_times and 'Numba JIT' in avg_times:
            print(f"\n{'='*80}")
            print("WHY NUMBA MAY BE FASTER THAN CUDA:")
            print(f"{'='*80}")
            print("CUDA has fixed overhead (~100-500μs) per call:")
            print("  • 6× cudaMalloc allocations")
            print("  • 2× Host→Device memory transfers")
            print("  • 1× Device→Host memory transfer")
            print("  • 6× cudaFree deallocations")
            print("  • Kernel launch overhead")
            print("\nNumba advantages for short sequences:")
            print("  • Zero memory transfer overhead")
            print("  • Direct CPU execution (cache-friendly)")
            print("  • LLVM-optimized native code")
            print("\nCUDA becomes faster when:")
            print("  • Sequence length > 5000 (parallelism >> overhead)")
            print("  • Batch processing (amortize overhead over many DTW calls)")
            print("  • Very long sequences > 10k (GPU dominates)")
            print("\nTo optimize CUDA for short sequences:")
            print("  • Implement memory pooling (reuse cudaMalloc)")
            print("  • Add batch API (compute multiple DTW at once)")
            print("  • Use pinned memory for faster transfers")


def test_accuracy():
    """
    Test accuracy: Compare results on identical sequences and simple cases.
    """
    print(f"\n{'='*80}")
    print("ACCURACY TEST")
    print(f"{'='*80}\n")
    
    if len(implementations) < 2:
        print("Need at least 2 implementations to compare accuracy")
        return
    
    # Test 1: Identical sequences (should give distance ~0)
    print("Test 1: Identical sequences")
    seq = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
    print(f"  Sequence: {seq}")
    
    for name, func in implementations.items():
        try:
            dist = func(seq, seq)
            print(f"  {name:<25} distance: {dist:.6f}")
        except Exception as e:
            print(f"  {name:<25} error: {e}")
    
    # Test 2: Simple different sequences
    print("\nTest 2: Different sequences")
    seq1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
    seq2 = np.array([1.5, 2.5, 3.5, 4.5, 5.5], dtype=np.float32)
    print(f"  Seq1: {seq1}")
    print(f"  Seq2: {seq2}")
    
    distances = {}
    for name, func in implementations.items():
        try:
            dist = func(seq1, seq2)
            distances[name] = dist
            print(f"  {name:<25} distance: {dist:.6f}")
        except Exception as e:
            print(f"  {name:<25} error: {e}")
    
    # Check agreement
    if len(distances) > 1:
        dist_values = list(distances.values())
        max_diff = max(dist_values) - min(dist_values)
        print(f"\n  Maximum difference: {max_diff:.6f}")
        if max_diff < 0.01:
            print("  ✓ All implementations agree (within 0.01)")
        else:
            print("  ⚠ Implementations show some variation")


def print_implementation_info():
    """
    Print information about each DTW implementation.
    """
    print(f"\n{'='*80}")
    print("IMPLEMENTATION DETAILS")
    print(f"{'='*80}\n")
    
    if 'Pure Python' in implementations:
        print("1. Pure Python")
        print("   - Naive implementation using nested loops")
        print("   - O(n*m) time complexity, O(n*m) space complexity")
        print("   - Slowest but easiest to understand")
        print("   - Reference implementation for correctness\n")
    
    if 'Numba JIT' in implementations:
        print("2. Numba JIT")
        print("   - JIT-compiled Python code using LLVM")
        print("   - Same algorithm as Pure Python but compiled")
        print("   - 10-100x faster than Pure Python")
        print("   - No GPU required, works on any CPU")
        print("   - Best for short-medium sequences (< 5000 points)")
        print("   - Zero memory transfer overhead\n")
    
    if 'CUDA (GPU)' in implementations:
        print("3. CUDA (GPU)")
        print("   - GPU-accelerated using NVIDIA CUDA")
        print("   - Parallel wavefront algorithm")
        print("   - 100-1000x faster than Pure Python for long sequences")
        print("   - Requires NVIDIA GPU and CUDA toolkit")
        print("   - Overhead: ~100-500μs per call (cudaMalloc/cudaMemcpy)")
        print("   - Best for: sequences > 5000 points or batch processing")
        print("   - Note: May be slower than Numba for short sequences\n")
    
    if 'dtaidistance (C)' in implementations:
        print("4. dtaidistance (C)")
        print("   - Optimized C implementation with Python bindings")
        print("   - Highly optimized CPU code")
        print("   - 50-200x faster than Pure Python")
        print("   - No GPU required\n")
    
    if 'fastdtw (approx)' in implementations:
        print("5. fastdtw (approx)")
        print("   - Approximate DTW using FastDTW algorithm")
        print("   - Reduces complexity using multilevel approach")
        print("   - May not give exact DTW distance")
        print("   - Good for very long sequences\n")
    
    print("PERFORMANCE EXPECTATIONS:")
    print("  • Sequences < 2000:    Numba ≈ dtaidistance > CUDA > Pure Python")
    print("  • Sequences 2000-5000: CUDA ≈ Numba ≈ dtaidistance > Pure Python")
    print("  • Sequences > 5000:    CUDA > dtaidistance ≈ Numba > Pure Python")
    print("  • Batch processing:    CUDA >> all others (amortizes overhead)\n")


def main():
    """
    Main benchmark script.
    """
    if len(implementations) == 0:
        print("\nERROR: No DTW implementations available!")
        print("\nTo install missing implementations:")
        print("  pip install numba         # JIT-compiled Python")
        print("  pip install dtaidistance  # Optimized CPU DTW")
        print("  pip install fastdtw       # Approximate DTW")
        print("\nFor CUDA DTW, rebuild the package:")
        print("  cd /path/to/pyfin && pip install -e .")
        sys.exit(1)
    
    print(f"Found {len(implementations)} DTW implementation(s)\n")
    
    # Print implementation details
    print_implementation_info()
    
    # Run accuracy tests first
    test_accuracy()
    
    # Run performance benchmarks
    # Choose sequence lengths to show crossover point where CUDA becomes faster
    if 'CUDA (GPU)' in implementations and 'Pure Python' not in implementations:
        # Show CUDA advantage for long sequences
        sequence_lengths = [100, 500, 1000, 2000, 5000, 10000]
        print("NOTE: Including long sequences to demonstrate GPU advantage")
    elif 'Pure Python' in implementations:
        # Limit max length because Pure Python is very slow
        sequence_lengths = [100, 500, 1000, 2000]
        print("NOTE: Limited sequence lengths due to Pure Python being slow")
    else:
        # Default range
        sequence_lengths = [100, 500, 1000, 2000, 5000]
    
    n_runs = 10
    
    print(f"\n{'='*80}")
    print("Starting performance benchmarks...")
    print(f"Sequence lengths: {sequence_lengths}")
    print(f"Runs per length: {n_runs}")
    print(f"{'='*80}\n")
    
    run_benchmark_suite(sequence_lengths, n_runs=n_runs)
    
    # Cleanup CUDA if used
    if 'CUDA (GPU)' in implementations:
        try:
            cleanup()
            print("\n✓ CUDA resources cleaned up")
        except:
            pass
    
    print("\nBenchmark complete!")


if __name__ == "__main__":
    main()
