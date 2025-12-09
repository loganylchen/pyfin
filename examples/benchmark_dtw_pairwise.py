#!/usr/bin/env python3
"""
Benchmark: Pairwise DTW - Batch vs Individual vs CPU implementations

Compares performance of:
1. fin._dtw.dtw_pairwise() - CUDA batch API (GPU)
2. fin._dtw.dtw_distance() loop - CUDA individual calls (GPU)
3. dtaidistance.distance_matrix_fast() - Optimized CPU batch
4. Pure Python loop - Naive CPU implementation
"""

import numpy as np
import time
from typing import Dict, List, Tuple
import sys


def benchmark_pairwise(
    sequences: np.ndarray,
    impl_name: str,
    compute_func,
    n_runs: int = 5,
    warmup: int = 1
) -> Tuple[np.ndarray, float, float]:
    """
    Benchmark a pairwise DTW implementation.
    
    Returns:
        Tuple of (distance_matrix, avg_time_ms, std_time_ms)
    """
    # Warm-up runs
    for _ in range(warmup):
        try:
            _ = compute_func(sequences)
        except:
            pass
    
    # Actual benchmark
    times = []
    distance_matrix = None
    
    for _ in range(n_runs):
        start = time.perf_counter()
        try:
            distance_matrix = compute_func(sequences)
        except Exception as e:
            print(f"  Error in {impl_name}: {e}")
            return None, None, None
        end = time.perf_counter()
        times.append((end - start) * 1000)  # Convert to ms
    
    avg_time = np.mean(times)
    std_time = np.std(times)
    
    return distance_matrix, avg_time, std_time


# =============================================================================
# Implementation Wrappers
# =============================================================================

def cuda_batch_pairwise(sequences):
    """CUDA batch API - optimal GPU usage"""
    from fin._dtw import dtw_pairwise
    return dtw_pairwise(sequences)


def cuda_individual_pairwise(sequences):
    """CUDA individual calls - suboptimal but functional"""
    from fin._dtw import dtw_distance
    n = len(sequences)
    dist_matrix = np.zeros((n, n), dtype=np.float32)
    
    for i in range(n):
        for j in range(i+1, n):
            dist = dtw_distance(sequences[i], sequences[j])
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    
    return dist_matrix


def dtaidistance_pairwise(sequences):
    """dtaidistance CPU batch - optimized C implementation"""
    from dtaidistance import dtw
    # dtaidistance expects list of arrays
    seq_list = [seq for seq in sequences]
    dist_matrix = dtw.distance_matrix_fast(seq_list)
    return np.array(dist_matrix)


def python_naive_pairwise(sequences):
    """Pure Python naive implementation"""
    def dtw_python(seq1, seq2):
        n, m = len(seq1), len(seq2)
        dtw_matrix = np.full((n + 1, m + 1), np.inf)
        dtw_matrix[0, 0] = 0.0
        
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                diff = seq1[i-1] - seq2[j-1]
                cost = diff * diff
                dtw_matrix[i, j] = cost + min(
                    dtw_matrix[i-1, j],
                    dtw_matrix[i, j-1],
                    dtw_matrix[i-1, j-1]
                )
        
        return dtw_matrix[n, m]
    
    n = len(sequences)
    dist_matrix = np.zeros((n, n), dtype=np.float32)
    
    for i in range(n):
        for j in range(i+1, n):
            dist = dtw_python(sequences[i], sequences[j])
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    
    return dist_matrix


def numba_pairwise(sequences):
    """Numba JIT-compiled pairwise"""
    from numba import jit
    
    @jit(nopython=True)
    def dtw_numba(seq1, seq2):
        n, m = len(seq1), len(seq2)
        dtw_matrix = np.full((n + 1, m + 1), np.inf)
        dtw_matrix[0, 0] = 0.0
        
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                diff = seq1[i-1] - seq2[j-1]
                cost = diff * diff
                
                min_prev = dtw_matrix[i-1, j-1]
                if dtw_matrix[i-1, j] < min_prev:
                    min_prev = dtw_matrix[i-1, j]
                if dtw_matrix[i, j-1] < min_prev:
                    min_prev = dtw_matrix[i, j-1]
                
                dtw_matrix[i, j] = cost + min_prev
        
        return dtw_matrix[n, m]
    
    n = len(sequences)
    dist_matrix = np.zeros((n, n), dtype=np.float32)
    
    for i in range(n):
        for j in range(i+1, n):
            dist = dtw_numba(sequences[i], sequences[j])
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist
    
    return dist_matrix


# =============================================================================
# Setup Implementations
# =============================================================================

implementations = {}

# 1. CUDA Batch (optimal)
try:
    from fin._dtw import dtw_pairwise, is_available
    if is_available():
        implementations['CUDA Batch (GPU)'] = cuda_batch_pairwise
        print("✓ CUDA Batch API available")
    else:
        print("✗ CUDA Batch API not available")
except ImportError as e:
    print(f"✗ CUDA Batch API import failed: {e}")

# 2. CUDA Individual (for comparison)
try:
    from fin._dtw import dtw_distance, is_available
    if is_available():
        implementations['CUDA Individual (GPU)'] = cuda_individual_pairwise
        print("✓ CUDA Individual API available")
except ImportError:
    pass

# 3. dtaidistance
try:
    from dtaidistance import dtw
    implementations['dtaidistance (CPU)'] = dtaidistance_pairwise
    print("✓ dtaidistance available")
except ImportError:
    print("✗ dtaidistance not available (install: pip install dtaidistance)")

# 4. Numba
try:
    from numba import jit
    implementations['Numba Loop (CPU)'] = numba_pairwise
    print("✓ Numba available")
except ImportError:
    print("✗ Numba not available (install: pip install numba)")

# 5. Pure Python (always available, but slow)
implementations['Pure Python (CPU)'] = python_naive_pairwise
print("✓ Pure Python available")

print()


def run_benchmark_suite():
    """
    Run comprehensive benchmarks comparing all implementations.
    """
    print("="*80)
    print("PAIRWISE DTW BENCHMARK SUITE")
    print("="*80)
    print(f"Implementations available: {len(implementations)}")
    print()
    
    if len(implementations) == 0:
        print("ERROR: No implementations available!")
        return
    
    # Test configurations: (num_sequences, seq_length)
    configs = [
        (5, 100),    # Small: quick test
        (10, 200),   # Small-medium
        (20, 300),   # Medium
        (50, 500),   # Large
    ]
    
    # Skip expensive tests for slow implementations
    if 'Pure Python (CPU)' in implementations and 'CUDA Batch (GPU)' not in implementations:
        configs = configs[:2]  # Only small tests
        print("⚠ WARNING: Running limited tests (Pure Python only is very slow)\n")
    
    all_results = {}
    
    for num_seq, seq_len in configs:
        num_pairs = (num_seq * (num_seq - 1)) // 2
        
        print(f"\n{'='*80}")
        print(f"Configuration: {num_seq} sequences × {seq_len} samples")
        print(f"Total pairs to compute: {num_pairs}")
        print(f"{'='*80}")
        
        # Generate random sequences
        np.random.seed(42)
        sequences = np.random.randn(num_seq, seq_len).astype(np.float32)
        
        config_results = {}
        reference_matrix = None
        
        for impl_name, impl_func in implementations.items():
            # Skip slow implementations for large tests
            if impl_name == 'Pure Python (CPU)' and num_seq > 10:
                print(f"\nSkipping {impl_name} (too slow for {num_seq} sequences)...")
                continue
            
            if impl_name == 'CUDA Individual (GPU)' and num_seq > 20:
                print(f"\nSkipping {impl_name} (inefficient for {num_seq} sequences)...")
                continue
            
            # Adjust runs based on implementation speed
            if impl_name == 'Pure Python (CPU)':
                n_runs = 2
            elif impl_name == 'CUDA Individual (GPU)':
                n_runs = 3
            else:
                n_runs = 5
            
            print(f"\nTesting {impl_name}... (runs={n_runs})")
            dist_matrix, avg_time, std_time = benchmark_pairwise(
                sequences, impl_name, impl_func, n_runs=n_runs
            )
            
            if dist_matrix is not None:
                print(f"  Avg time: {avg_time:.2f} ± {std_time:.2f} ms")
                print(f"  Time per pair: {avg_time / num_pairs:.2f} ms")
                print(f"  Throughput: {num_pairs / (avg_time / 1000):.1f} pairs/sec")
                
                # Store reference matrix from first implementation
                if reference_matrix is None:
                    reference_matrix = dist_matrix
                
                # Check correctness
                if reference_matrix is not None:
                    max_diff = np.max(np.abs(dist_matrix - reference_matrix))
                    print(f"  Max difference from reference: {max_diff:.6f}")
                    if max_diff < 0.1:
                        print(f"  ✓ Results match")
                    else:
                        print(f"  ⚠ Results differ!")
                
                config_results[impl_name] = {
                    'time': avg_time,
                    'std': std_time,
                    'throughput': num_pairs / (avg_time / 1000)
                }
            else:
                print(f"  Benchmark failed")
        
        # Print comparison for this configuration
        if len(config_results) > 1:
            print(f"\n{'-'*80}")
            print(f"Comparison for {num_seq} sequences × {seq_len} samples:")
            print(f"{'-'*80}")
            
            # Find fastest
            fastest = min(config_results.items(), key=lambda x: x[1]['time'])
            print(f"Fastest: {fastest[0]} ({fastest[1]['time']:.2f} ms)")
            
            # Speedup comparison
            if 'CUDA Batch (GPU)' in config_results:
                cuda_time = config_results['CUDA Batch (GPU)']['time']
                print(f"\nSpeedup vs CUDA Batch:")
                for name, res in config_results.items():
                    if name != 'CUDA Batch (GPU)':
                        speedup = res['time'] / cuda_time
                        print(f"  {name:<30} {speedup:>6.1f}x slower")
            
            # Show throughput comparison
            print(f"\nThroughput comparison:")
            sorted_by_throughput = sorted(config_results.items(), 
                                         key=lambda x: x[1]['throughput'], 
                                         reverse=True)
            for name, res in sorted_by_throughput:
                print(f"  {name:<30} {res['throughput']:>8.1f} pairs/sec")
        
        all_results[(num_seq, seq_len)] = config_results
    
    # Summary
    print(f"\n\n{'='*80}")
    print("BENCHMARK SUMMARY")
    print(f"{'='*80}\n")
    
    # Create summary table
    print(f"{'Configuration':<25} ", end='')
    impl_names = list(implementations.keys())
    for name in impl_names[:3]:  # Show first 3
        print(f"{name:<25} ", end='')
    print()
    print("-" * 80)
    
    for (num_seq, seq_len), results in all_results.items():
        num_pairs = (num_seq * (num_seq - 1)) // 2
        print(f"{num_seq}seq × {seq_len}len ({num_pairs} pairs)".ljust(25), end=' ')
        
        for name in impl_names[:3]:
            if name in results:
                print(f"{results[name]['time']:>7.1f} ms".ljust(25), end=' ')
            else:
                print(f"{'---':<25}", end=' ')
        print()
    
    # Key insights
    print(f"\n{'='*80}")
    print("KEY INSIGHTS")
    print(f"{'='*80}\n")
    
    if 'CUDA Batch (GPU)' in implementations and 'CUDA Individual (GPU)' in implementations:
        print("1. CUDA Batch vs Individual:")
        print("   Batch API amortizes GPU overhead over all pairs")
        print("   Speedup increases with batch size:")
        for (num_seq, seq_len), results in all_results.items():
            if 'CUDA Batch (GPU)' in results and 'CUDA Individual (GPU)' in results:
                speedup = results['CUDA Individual (GPU)']['time'] / results['CUDA Batch (GPU)']['time']
                print(f"   - {num_seq} sequences: {speedup:.1f}x faster")
        print()
    
    if 'CUDA Batch (GPU)' in implementations and 'dtaidistance (CPU)' in implementations:
        print("2. GPU vs CPU (dtaidistance):")
        print("   GPU advantage scales with problem size:")
        for (num_seq, seq_len), results in all_results.items():
            if 'CUDA Batch (GPU)' in results and 'dtaidistance (CPU)' in results:
                speedup = results['dtaidistance (CPU)']['time'] / results['CUDA Batch (GPU)']['time']
                print(f"   - {num_seq} sequences: {speedup:.1f}x faster")
        print()
    
    print("3. When to use each implementation:")
    print("   • CUDA Batch: > 10 sequences, maximum performance")
    print("   • CUDA Individual: < 10 pairs, different lengths")
    print("   • dtaidistance: CPU-only systems, medium batches")
    print("   • Numba: CPU-only, small batches, no dependencies")
    print("   • Pure Python: Educational purposes only")
    
    # Performance expectations
    if 'CUDA Batch (GPU)' in implementations:
        print(f"\n{'='*80}")
        print("CUDA BATCH PERFORMANCE EXPECTATIONS")
        print(f"{'='*80}\n")
        print("Approximate throughput for CUDA Batch:")
        print("  • 10 sequences × 200:    ~200 pairs/sec")
        print("  • 20 sequences × 300:    ~400 pairs/sec")
        print("  • 50 sequences × 500:    ~1000 pairs/sec")
        print("  • 100 sequences × 1000:  ~2000 pairs/sec")
        print()
        print("Memory requirements:")
        print("  • Input: num_seq × seq_len × 4 bytes")
        print("  • Output: num_seq² × 4 bytes")
        print("  • Temp: ~seq_len × num_seq × 8 bytes")
        print()
        print("Example: 100 sequences × 1000 samples")
        print("  • Input: 0.4 MB")
        print("  • Output: 0.04 MB")
        print("  • Temp: ~0.8 MB")
        print("  • Total: ~1.24 MB (fits on any GPU)")


def test_accuracy():
    """
    Test accuracy across implementations with known cases.
    """
    print(f"\n{'='*80}")
    print("ACCURACY VERIFICATION")
    print(f"{'='*80}\n")
    
    if len(implementations) < 2:
        print("Need at least 2 implementations to compare")
        return
    
    # Test 1: Identical sequences (all distances should be 0)
    print("Test 1: Identical sequences (should all be 0)")
    seq = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
    sequences = np.tile(seq, (3, 1))
    
    for name, func in list(implementations.items())[:3]:  # Test first 3
        try:
            dist_matrix = func(sequences)
            max_dist = np.max(dist_matrix)
            print(f"  {name:<30} max distance: {max_dist:.6f}")
            if max_dist < 0.001:
                print(f"    ✓ Correct (all zeros)")
            else:
                print(f"    ⚠ Unexpected non-zero distances")
        except Exception as e:
            print(f"  {name:<30} error: {e}")
    
    # Test 2: Known different sequences
    print("\nTest 2: Different sequences")
    sequences = np.array([
        [1.0, 2.0, 3.0, 4.0, 5.0],
        [1.5, 2.5, 3.5, 4.5, 5.5],
        [2.0, 3.0, 4.0, 5.0, 6.0]
    ], dtype=np.float32)
    
    results = {}
    for name, func in list(implementations.items())[:3]:
        try:
            dist_matrix = func(sequences)
            results[name] = dist_matrix
            print(f"  {name:<30}")
            print(f"    Distance [0,1]: {dist_matrix[0,1]:.6f}")
            print(f"    Distance [0,2]: {dist_matrix[0,2]:.6f}")
            print(f"    Distance [1,2]: {dist_matrix[1,2]:.6f}")
        except Exception as e:
            print(f"  {name:<30} error: {e}")
    
    # Check agreement
    if len(results) > 1:
        matrices = list(results.values())
        max_diff = np.max([np.max(np.abs(matrices[0] - m)) for m in matrices[1:]])
        print(f"\n  Maximum difference between implementations: {max_diff:.6f}")
        if max_diff < 0.01:
            print("  ✓ All implementations agree")
        else:
            print("  ⚠ Implementations show variation")


def main():
    """
    Main benchmark script.
    """
    if len(implementations) == 0:
        print("\nERROR: No implementations available!")
        print("\nTo install missing implementations:")
        print("  pip install numba dtaidistance")
        print("\nFor CUDA support:")
        print("  cd /path/to/pyfin && pip install -e .")
        sys.exit(1)
    
    print(f"Found {len(implementations)} implementation(s)\n")
    
    # Run accuracy tests
    test_accuracy()
    
    # Run performance benchmarks
    run_benchmark_suite()
    
    # Cleanup
    if 'CUDA Batch (GPU)' in implementations:
        try:
            from fin._dtw import cleanup
            cleanup()
            print("\n✓ CUDA resources cleaned up")
        except:
            pass
    
    print("\n" + "="*80)
    print("Benchmark complete!")
    print("="*80)


if __name__ == "__main__":
    main()
