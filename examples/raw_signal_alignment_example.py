#!/usr/bin/env python3
"""
Example: Raw Signal to Sequence Alignment

Demonstrates the complete pipeline for aligning nanopore raw signal
directly to a reference sequence using real pore models.

Pipeline:
1. Raw signal → Event detection (MAD-based adapter trimming)
2. Events → 3-state HMM alignment with soft-clipping
3. Output: Base-to-event mapping with proper handling of adapters

This uses f5c's full 3-state HMM algorithm with:
- MATCH state: Event matches kmer
- BAD_EVENT state: Noisy events to skip
- KMER_SKIP state: Kmers with no events
- Dynamic transition probabilities
- Real pore models (RNA R9.4, RNA004, DNA R9.4)
"""

import numpy as np
from fin._f5c import eventalign

def generate_synthetic_signal(sequence, kmer_size=5, is_rna=False):
    """
    Generate synthetic nanopore signal for a given sequence.
    
    Args:
        sequence: DNA/RNA sequence string
        kmer_size: Size of kmers (5 or 9)
        is_rna: Whether this is RNA (True) or DNA (False)
    
    Returns:
        1D numpy array of float32 raw signal values
    """
    # Simplified model levels (these would come from real pore models)
    base_to_level = {'A': 100.0, 'C': 95.0, 'G': 105.0, 'T': 90.0, 'U': 90.0}
    
    signal = []
    
    # Add adapter region (noisy signal at start)
    adapter_samples = 200
    signal.extend(np.random.normal(80, 10, adapter_samples))
    
    # Generate signal for each kmer (with some noise)
    n_kmers = len(sequence) - kmer_size + 1
    for i in range(n_kmers):
        kmer = sequence[i:i+kmer_size]
        
        # Average level based on bases in kmer
        level = np.mean([base_to_level.get(b, 100.0) for b in kmer])
        
        # Each kmer generates ~10 samples (with variation)
        n_samples = np.random.randint(8, 15)
        kmer_signal = np.random.normal(level, 2.5, n_samples)
        signal.extend(kmer_signal)
    
    # Add adapter region (noisy signal at end)
    signal.extend(np.random.normal(85, 10, adapter_samples))
    
    return np.array(signal, dtype=np.float32)


def example_rna_alignment():
    """Example: RNA sequence alignment with R9.4 5-mer model"""
    print("=" * 70)
    print("Example 1: RNA Alignment (R9.4 5-mer model)")
    print("=" * 70)
    
    # RNA sequence (note: U instead of T)
    sequence = "AUGCGAUACGUAGCUAGCUAGCUAGCUGCUAGCUAGCUA"
    print(f"\nSequence: {sequence}")
    print(f"Length: {len(sequence)} bases")
    
    # Generate synthetic signal
    raw_signal = generate_synthetic_signal(sequence, kmer_size=5, is_rna=True)
    print(f"Raw signal: {len(raw_signal)} samples")
    
    # Align raw signal to sequence
    result = eventalign(raw_signal, sequence, is_rna=1, kmer_size=5)
    
    # Display results
    print(f"\nAlignment Results:")
    print(f"  Events detected: {result['n_events']}")
    print(f"  Aligned pairs: {result['n_aligned_pairs']}")
    print(f"  Scaling: scale={result['scaling']['scale']:.3f}, shift={result['scaling']['shift']:.3f}")
    
    # Show base-to-event mapping for first 10 kmers
    print(f"\nBase-to-Event Mapping (first 10 kmers):")
    base_to_event = result['base_to_event_map']
    for i in range(min(10, len(base_to_event))):
        kmer_pos = i
        kmer_seq = sequence[i:i+5]
        event_indices = base_to_event[i]
        n_events = len(event_indices) if event_indices else 0
        print(f"  Kmer {kmer_pos:3d} ({kmer_seq}): {n_events} events")
    
    return result


def example_dna_alignment():
    """Example: DNA sequence alignment with R9.4 5-mer model"""
    print("\n" + "=" * 70)
    print("Example 2: DNA Alignment (R9.4 5-mer model)")
    print("=" * 70)
    
    # DNA sequence
    sequence = "ATGCGATACGTAGCTAGCTAGCTAGCTGCTAGCTAGCTA"
    print(f"\nSequence: {sequence}")
    print(f"Length: {len(sequence)} bases")
    
    # Generate synthetic signal
    raw_signal = generate_synthetic_signal(sequence, kmer_size=5, is_rna=False)
    print(f"Raw signal: {len(raw_signal)} samples")
    
    # Align raw signal to sequence
    result = eventalign(raw_signal, sequence, is_rna=0, kmer_size=5)
    
    # Display results
    print(f"\nAlignment Results:")
    print(f"  Events detected: {result['n_events']}")
    print(f"  Aligned pairs: {result['n_aligned_pairs']}")
    print(f"  Scaling: scale={result['scaling']['scale']:.3f}, shift={result['scaling']['shift']:.3f}")
    
    # Show base-to-event mapping
    print(f"\nBase-to-Event Mapping (first 10 kmers):")
    base_to_event = result['base_to_event_map']
    for i in range(min(10, len(base_to_event))):
        kmer_pos = i
        kmer_seq = sequence[i:i+5]
        event_indices = base_to_event[i]
        n_events = len(event_indices) if event_indices else 0
        print(f"  Kmer {kmer_pos:3d} ({kmer_seq}): {n_events} events")
    
    return result


def example_rna004_alignment():
    """Example: RNA sequence alignment with RNA004 9-mer model"""
    print("\n" + "=" * 70)
    print("Example 3: RNA Alignment (RNA004 9-mer model)")
    print("=" * 70)
    
    # Longer RNA sequence for 9-mer model
    sequence = "AUGCGAUACGUAGCUAGCUAGCUAGCUGCUAGCUAGCUAGGCUAGCUAGCUA"
    print(f"\nSequence: {sequence}")
    print(f"Length: {len(sequence)} bases")
    
    # Generate synthetic signal
    raw_signal = generate_synthetic_signal(sequence, kmer_size=9, is_rna=True)
    print(f"Raw signal: {len(raw_signal)} samples")
    
    # Align raw signal to sequence (will auto-select RNA004 9-mer model)
    result = eventalign(raw_signal, sequence, is_rna=1, kmer_size=9)
    
    # Display results
    print(f"\nAlignment Results:")
    print(f"  Events detected: {result['n_events']}")
    print(f"  Aligned pairs: {result['n_aligned_pairs']}")
    print(f"  Scaling: scale={result['scaling']['scale']:.3f}, shift={result['scaling']['shift']:.3f}")
    
    # Show base-to-event mapping
    print(f"\nBase-to-Event Mapping (first 10 kmers):")
    base_to_event = result['base_to_event_map']
    for i in range(min(10, len(base_to_event))):
        kmer_pos = i
        kmer_seq = sequence[i:i+9]
        event_indices = base_to_event[i]
        n_events = len(event_indices) if event_indices else 0
        print(f"  Kmer {kmer_pos:3d} ({kmer_seq}): {n_events} events")
    
    return result


def example_with_adapters():
    """Example showing how soft-clipping handles untrimmed adapters"""
    print("\n" + "=" * 70)
    print("Example 4: Soft-Clipping with Untrimmed Adapters")
    print("=" * 70)
    
    sequence = "ATGCGATACGTAGCTAGCTAGCTAGCTGCTAGCTAGCTA"
    print(f"\nSequence: {sequence}")
    
    # Generate signal with long adapter regions
    raw_signal = generate_synthetic_signal(sequence, kmer_size=5, is_rna=False)
    
    # Add extra noisy adapter signal
    pre_adapter = np.random.normal(75, 15, 500)  # Long pre-adapter
    post_adapter = np.random.normal(80, 15, 500)  # Long post-adapter
    raw_signal = np.concatenate([pre_adapter, raw_signal, post_adapter]).astype(np.float32)
    
    print(f"Raw signal: {len(raw_signal)} samples (includes long adapters)")
    
    # Align - soft-clipping should handle adapters automatically
    result = eventalign(raw_signal, sequence, is_rna=0, kmer_size=5)
    
    print(f"\nAlignment Results:")
    print(f"  Events detected: {result['n_events']}")
    print(f"  Aligned pairs: {result['n_aligned_pairs']}")
    print(f"  Scaling: scale={result['scaling']['scale']:.3f}, shift={result['scaling']['shift']:.3f}")
    print(f"\nNote: Soft-clipping automatically skipped adapter events!")
    print(f"      (TRANS_START_TO_CLIP=0.5, TRANS_CLIP_SELF=0.9)")
    
    return result


if __name__ == "__main__":
    print("\nRaw Signal to Sequence Alignment Examples")
    print("Using f5c's 3-state HMM with real pore models\n")
    
    # Run examples
    try:
        example_rna_alignment()
        example_dna_alignment()
        example_rna004_alignment()
        example_with_adapters()
        
        print("\n" + "=" * 70)
        print("All examples completed successfully!")
        print("=" * 70)
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
