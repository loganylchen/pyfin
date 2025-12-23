#!/usr/bin/env python
"""
Example usage of the fin._eventalign API.

This example demonstrates:
1. Event detection from raw nanopore signal
2. Loading pore models (RNA002 and RNA004)
"""

import numpy as np
from fin._eventalign import getevents, set_model, MODEL_RNA002, MODEL_RNA004


def example_getevents():
    """Demonstrate event detection from raw signal."""
    print("=" * 60)
    print("Example 1: Event Detection")
    print("=" * 60)

    # Generate synthetic signal (simulating nanopore current)
    # In practice, this would be actual raw signal data from a FAST5 file
    np.random.seed(42)
    n_samples = 10000
    signal = np.random.randn(n_samples).astype(np.float32) * 10 + 120

    print(f"Input signal: {n_samples} samples")
    print(f"Signal range: [{signal.min():.2f}, {signal.max():.2f}]")

    # Detect events
    events = getevents(signal)

    print(f"\nDetected {events['n_events']} events")
    print(f"\nFirst 5 events:")
    print(f"  {'Start':>10} {'Length':>10} {'Mean':>10} {'Stdv':>10}")
    print(f"  {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
    for i in range(min(5, events['n_events'])):
        print(f"  {events['starts'][i]:10d} {events['lengths'][i]:10.1f} "
              f"{events['means'][i]:10.2f} {events['stdvs'][i]:10.2f}")


def example_set_model():
    """Demonstrate loading pore models."""
    print("\n" + "=" * 60)
    print("Example 2: Loading Pore Models")
    print("=" * 60)

    # Load RNA002 model (k=5)
    print("\nLoading RNA002 model (k=5)...")
    model_002 = set_model(MODEL_RNA002)
    print(f"  K-mer size: {model_002['kmer_size']}")
    print(f"  Number of k-mers: {model_002['num_kmer']}")
    print(f"  Mean level range: [{model_002['level_means'].min():.2f}, "
          f"{model_002['level_means'].max():.2f}]")
    print(f"  Stdv range: [{model_002['level_stdvs'].min():.4f}, "
          f"{model_002['level_stdvs'].max():.4f}]")

    # Load RNA004 model (k=9)
    print("\nLoading RNA004 model (k=9)...")
    model_004 = set_model(MODEL_RNA004)
    print(f"  K-mer size: {model_004['kmer_size']}")
    print(f"  Number of k-mers: {model_004['num_kmer']}")
    print(f"  Mean level range: [{model_004['level_means'].min():.2f}, "
          f"{model_004['level_means'].max():.2f}]")
    print(f"  Stdv range: [{model_004['level_stdvs'].min():.4f}, "
          f"{model_004['level_stdvs'].max():.4f}]")


def example_kmer_lookup():
    """Demonstrate looking up k-mer model values."""
    print("\n" + "=" * 60)
    print("Example 3: K-mer Model Lookup")
    print("=" * 60)

    # Load RNA002 model
    model = set_model(MODEL_RNA002)

    # Convert kmer string to index
    # RNA002 uses k=5, bases are A=0, C=1, G=2, T=3
    # Index = sum(base * 4^position) for position in 0..k-1
    def kmer_to_index(kmer: str) -> int:
        base_to_val = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'U': 0}
        index = 0
        for i, base in enumerate(kmer):
            index = index * 4 + base_to_val[base.upper()]
        return index

    # Example kmers
    kmers = ["AAAAA", "CCCCC", "GGGGG", "UUUUU", "ACGUU"]

    print(f"\nRNA002 model values for example k-mers:")
    print(f"  {'K-mer':>10} {'Index':>10} {'Mean':>10} {'Stdv':>10}")
    print(f"  {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
    for kmer in kmers:
        idx = kmer_to_index(kmer)
        mean = model['level_means'][idx]
        stdv = model['level_stdvs'][idx]
        print(f"  {kmer:>10} {idx:10d} {mean:10.2f} {stdv:10.4f}")


if __name__ == "__main__":
    example_getevents()
    example_set_model()
    example_kmer_lookup()
    print("\n" + "=" * 60)
    print("All examples completed successfully!")
    print("=" * 60)
