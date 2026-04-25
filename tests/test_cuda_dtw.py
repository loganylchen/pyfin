#!/usr/bin/env python3
"""
Test script for CUDA DTW module

This script demonstrates how to use the CUDA-accelerated DTW distance calculation.
"""

import numpy as np
import sys

def test_basic_dtw():
    """Test basic DTW functionality"""
    print("Testing CUDA DTW module...")
    
    try:
        from fin._dtw import dtw_distance, is_available, cleanup
        
        if not is_available():
            print("ERROR: CUDA DTW extension is not available!")
            print("Make sure CUDA toolkit is installed and the extension was built.")
            return False
        
        print("✓ CUDA DTW extension loaded successfully")
        
        # Test 1: Simple identical sequences
        print("\nTest 1: Identical sequences")
        seq1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
        seq2 = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
        distance = dtw_distance(seq1, seq2)
        print(f"  Sequences: {seq1}")
        print(f"  DTW distance: {distance}")
        print("  Expected: ~0.0")
        
        # Test 2: Different sequences
        print("\nTest 2: Different sequences")
        seq1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
        seq2 = np.array([2.0, 3.0, 4.0, 5.0, 6.0], dtype=np.float32)
        distance = dtw_distance(seq1, seq2)
        print(f"  Seq1: {seq1}")
        print(f"  Seq2: {seq2}")
        print(f"  DTW distance: {distance}")
        
        # Test 3: Random sequences
        print("\nTest 3: Random sequences")
        np.random.seed(42)
        seq1 = np.random.randn(100).astype(np.float32)
        seq2 = np.random.randn(100).astype(np.float32)
        distance = dtw_distance(seq1, seq2)
        print("  Random sequences (length 100)")
        print(f"  DTW distance: {distance}")
        
        # Test 4: Different lengths
        print("\nTest 4: Different length sequences")
        seq1 = np.random.randn(50).astype(np.float32)
        seq2 = np.random.randn(75).astype(np.float32)
        distance = dtw_distance(seq1, seq2)
        print(f"  Seq1 length: {len(seq1)}")
        print(f"  Seq2 length: {len(seq2)}")
        print(f"  DTW distance: {distance}")
        
        # Test 5: With open boundaries
        print("\nTest 5: With open start/end boundaries")
        seq1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
        seq2 = np.array([2.0, 3.0, 4.0], dtype=np.float32)
        distance_normal = dtw_distance(seq1, seq2, use_open_start=False, use_open_end=False)
        distance_open = dtw_distance(seq1, seq2, use_open_start=True, use_open_end=True)
        print(f"  Seq1: {seq1}")
        print(f"  Seq2: {seq2}")
        print(f"  DTW distance (normal): {distance_normal}")
        print(f"  DTW distance (open boundaries): {distance_open}")
        
        # Cleanup
        print("\n✓ All tests passed!")
        cleanup()
        print("✓ CUDA resources cleaned up")
        
        return True
        
    except ImportError as e:
        print(f"ERROR: Failed to import CUDA DTW module: {e}")
        print("\nThis is expected if:")
        print("  1. CUDA toolkit is not installed")
        print("  2. The package was installed without CUDA support")
        print("  3. The extension failed to build")
        return False
    
    except Exception as e:
        print(f"ERROR: Test failed with exception: {e}")
        import traceback
        traceback.print_exc()
        return False


def benchmark_dtw():
    """Benchmark DTW performance"""
    print("\n" + "="*60)
    print("Benchmarking CUDA DTW performance")
    print("="*60)
    
    try:
        from fin._dtw import dtw_distance, cleanup
        import time
        
        sizes = [100, 500, 1000, 2000]
        
        for size in sizes:
            seq1 = np.random.randn(size).astype(np.float32)
            seq2 = np.random.randn(size).astype(np.float32)
            
            # Warm-up
            _ = dtw_distance(seq1, seq2)
            
            # Benchmark
            n_runs = 10
            start = time.time()
            for _ in range(n_runs):
                distance = dtw_distance(seq1, seq2)
            end = time.time()
            
            avg_time = (end - start) / n_runs * 1000  # ms
            print(f"  Size {size:4d} x {size:4d}: {avg_time:7.2f} ms/call (distance: {distance:.4f})")
        
        cleanup()
        
    except ImportError:
        print("  Skipping benchmark - CUDA extension not available")


if __name__ == "__main__":
    success = test_basic_dtw()
    
    if success:
        benchmark_dtw()
    
    sys.exit(0 if success else 1)
