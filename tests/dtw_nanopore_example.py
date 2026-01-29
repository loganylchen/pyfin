#!/usr/bin/env python3
"""
Simple example demonstrating CUDA DTW usage for nanopore signal analysis
"""

import numpy as np
from fin._dtw import dtw_distance, is_available, cleanup

def main():
    print("CUDA DTW Example for Nanopore Signals")
    print("=" * 50)
    
    # Check if CUDA is available
    if not is_available():
        print("\nWARNING: CUDA DTW extension is not available.")
        print("This example requires CUDA support.")
        return
    
    print("\n✓ CUDA DTW extension is available\n")
    
    # Simulate nanopore signal data (normalized current values)
    # In real usage, these would come from fast5/pod5 files
    
    # Example 1: Compare similar signals
    print("Example 1: Comparing similar signals")
    print("-" * 50)
    
    # Reference signal (e.g., from a known base sequence)
    reference_signal = np.array([
        1.2, 1.5, 1.8, 2.1, 2.3, 2.5, 2.4, 2.2, 2.0, 1.8,
        1.6, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.9, 1.0
    ], dtype=np.float32)
    
    # Test signal with slight noise
    test_signal = reference_signal + np.random.normal(0, 0.1, len(reference_signal)).astype(np.float32)
    
    distance = dtw_distance(reference_signal, test_signal)
    print(f"Reference signal length: {len(reference_signal)}")
    print(f"Test signal length: {len(test_signal)}")
    print(f"DTW distance: {distance:.4f}")
    print("(Lower distance indicates more similar signals)\n")
    
    # Example 2: Compare different signals
    print("Example 2: Comparing different signals")
    print("-" * 50)
    
    signal1 = np.array([1.0, 1.5, 2.0, 2.5, 3.0], dtype=np.float32)
    signal2 = np.array([3.0, 2.5, 2.0, 1.5, 1.0], dtype=np.float32)
    
    distance = dtw_distance(signal1, signal2)
    print(f"Signal 1: {signal1}")
    print(f"Signal 2: {signal2}")
    print(f"DTW distance: {distance:.4f}\n")
    
    # Example 3: Handle different length signals
    print("Example 3: Different length signals")
    print("-" * 50)
    
    # Longer reference signal
    long_signal = np.random.randn(100).astype(np.float32)
    # Shorter test signal
    short_signal = np.random.randn(50).astype(np.float32)
    
    distance = dtw_distance(long_signal, short_signal)
    print(f"Long signal length: {len(long_signal)}")
    print(f"Short signal length: {len(short_signal)}")
    print(f"DTW distance: {distance:.4f}")
    print("(DTW naturally handles different lengths)\n")
    
    # Example 4: Using open boundaries
    print("Example 4: Open boundary alignment")
    print("-" * 50)
    
    # Subsequence matching with open boundaries
    main_signal = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], dtype=np.float32)
    sub_signal = np.array([3, 4, 5, 6, 7], dtype=np.float32)
    
    dist_normal = dtw_distance(main_signal, sub_signal, 
                               use_open_start=False, use_open_end=False)
    dist_open = dtw_distance(main_signal, sub_signal,
                            use_open_start=True, use_open_end=True)
    
    print(f"Main signal: {main_signal}")
    print(f"Sub signal: {sub_signal}")
    print(f"DTW distance (normal): {dist_normal:.4f}")
    print(f"DTW distance (open boundaries): {dist_open:.4f}")
    print("(Open boundaries help with subsequence matching)\n")
    
    # Example 5: Batch comparison
    print("Example 5: Batch signal comparison")
    print("-" * 50)
    
    reference = np.random.randn(100).astype(np.float32)
    test_signals = [np.random.randn(100).astype(np.float32) for _ in range(5)]
    
    print(f"Comparing {len(test_signals)} signals against reference...")
    distances = [dtw_distance(reference, sig) for sig in test_signals]
    
    for i, dist in enumerate(distances, 1):
        print(f"  Signal {i}: distance = {dist:.4f}")
    
    best_match = np.argmin(distances)
    print(f"\nBest match: Signal {best_match + 1} (distance: {distances[best_match]:.4f})")
    
    # Clean up CUDA resources
    print("\n" + "=" * 50)
    cleanup()
    print("✓ CUDA resources cleaned up")
    print("\nExample completed successfully!")

if __name__ == "__main__":
    main()
