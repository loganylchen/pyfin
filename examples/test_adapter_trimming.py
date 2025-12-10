#!/usr/bin/env python3
"""
Test script to demonstrate adapter trimming functionality
"""

import numpy as np
import matplotlib.pyplot as plt

# Mock implementation for demonstration
# In reality, this would call your compiled C library

def simulate_nanopore_signal(n_samples=10000):
    """
    Simulate a nanopore signal with:
    - Adapter region at start (low variation)
    - Actual DNA sequence (high variation)
    - Adapter region at end (low variation)
    """
    np.random.seed(42)
    
    # Start adapter: low variation, stable signal
    adapter_start = np.random.normal(100, 2, 500)
    
    # DNA sequence: high variation, different k-mer levels
    dna_signal = []
    kmer_levels = [90, 100, 110, 95, 105, 115, 92, 108]
    for _ in range(8000 // len(kmer_levels)):
        for level in kmer_levels:
            # Each k-mer has some noise but distinct mean
            kmer_signal = np.random.normal(level, 8, len(kmer_levels))
            dna_signal.extend(kmer_signal)
    
    dna_signal = np.array(dna_signal[:8000])
    
    # End adapter: low variation, stable signal
    adapter_end = np.random.normal(105, 2, 1500)
    
    # Combine all parts
    signal = np.concatenate([adapter_start, dna_signal, adapter_end])
    
    return signal.astype(np.float32)


def calculate_mad_per_chunk(signal, chunk_size=100):
    """Calculate MAD for each chunk"""
    n_chunks = len(signal) // chunk_size
    mads = []
    
    for i in range(n_chunks):
        chunk = signal[i * chunk_size:(i + 1) * chunk_size]
        median = np.median(chunk)
        mad = np.median(np.abs(chunk - median)) * 1.4826
        mads.append(mad)
    
    return np.array(mads)


def visualize_adapter_trimming():
    """Visualize the adapter trimming process"""
    signal = simulate_nanopore_signal()
    chunk_size = 100
    
    # Calculate MAD per chunk
    mads = calculate_mad_per_chunk(signal, chunk_size)
    threshold = np.percentile(mads, 0)  # 0th percentile (minimum)
    
    # Find trim points
    trim_start = 0
    for i, mad in enumerate(mads):
        if mad > threshold:
            trim_start = i * chunk_size
            break
    
    trim_end = len(signal)
    for i in range(len(mads) - 1, -1, -1):
        if mads[i] > threshold:
            trim_end = (i + 1) * chunk_size
            break
    
    # Add fixed trimming (200 from start, 10 from end)
    trim_start += 200
    trim_end -= 10
    
    # Create visualization
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    
    # Plot 1: Raw signal with trim regions highlighted
    ax = axes[0]
    time = np.arange(len(signal))
    ax.plot(time, signal, 'b-', linewidth=0.5, alpha=0.7)
    ax.axvspan(0, trim_start, alpha=0.3, color='red', label='Trimmed (Adapter)')
    ax.axvspan(trim_end, len(signal), alpha=0.3, color='red')
    ax.axvspan(trim_start, trim_end, alpha=0.2, color='green', label='Retained (DNA)')
    ax.set_xlabel('Sample Index')
    ax.set_ylabel('Signal (pA)')
    ax.set_title('Nanopore Signal with Adapter Trimming')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: MAD per chunk
    ax = axes[1]
    chunk_positions = np.arange(len(mads)) * chunk_size + chunk_size / 2
    ax.bar(chunk_positions, mads, width=chunk_size * 0.8, alpha=0.7)
    ax.axhline(threshold, color='r', linestyle='--', label=f'Threshold: {threshold:.2f}')
    ax.axvline(trim_start, color='g', linestyle='-', linewidth=2, label='Trim Start')
    ax.axvline(trim_end, color='g', linestyle='-', linewidth=2, label='Trim End')
    ax.set_xlabel('Sample Index')
    ax.set_ylabel('MAD (Median Absolute Deviation)')
    ax.set_title('Signal Variation per Chunk (MAD-based Adapter Detection)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Trimmed signal only
    ax = axes[2]
    trimmed_signal = signal[trim_start:trim_end]
    trimmed_time = np.arange(len(trimmed_signal))
    ax.plot(trimmed_time, trimmed_signal, 'g-', linewidth=0.5)
    ax.set_xlabel('Sample Index (after trimming)')
    ax.set_ylabel('Signal (pA)')
    ax.set_title('Trimmed Signal (Adapters Removed)')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('/Users/logan/Projects/pyfin/examples/adapter_trimming_demo.png', dpi=150)
    print(f"Visualization saved to: /Users/logan/Projects/pyfin/examples/adapter_trimming_demo.png")
    
    # Print statistics
    print(f"\nAdapter Trimming Statistics:")
    print(f"  Original signal length: {len(signal)} samples")
    print(f"  Trimmed from start: {trim_start} samples")
    print(f"  Trimmed from end: {len(signal) - trim_end} samples")
    print(f"  Retained signal length: {trim_end - trim_start} samples")
    print(f"  Percentage retained: {100 * (trim_end - trim_start) / len(signal):.1f}%")
    print(f"\n  MAD threshold: {threshold:.3f}")
    print(f"  Mean MAD in adapters: {np.mean(mads[:5]):.3f}")
    print(f"  Mean MAD in DNA region: {np.mean(mads[5:-15]):.3f}")


if __name__ == "__main__":
    print("Demonstrating Nanopore Adapter Trimming (MAD-based method)")
    print("=" * 60)
    visualize_adapter_trimming()
    print("\nAdapter trimming removes regions with low signal variation")
    print("(adapters, open pore) from the beginning and end of reads.")
