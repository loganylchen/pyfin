#!/usr/bin/env python3
"""
Example usage of the new fin._eventalign Python wrapper

This script demonstrates how to use the eventalign module to align
nanopore signals to reference sequences, with both CPU and GPU support.
"""

import numpy as np
import sys

# Try to import the eventalign module
try:
    from fin._eventalign import EventAligner, align_single, align_batch, MODEL_RNA002, MODEL_RNA004
    print("✓ Successfully imported fin._eventalign module")
except ImportError as e:
    print(f"✗ Failed to import fin._eventalign: {e}")
    print("\nTo build the module, run:")
    print("  pip install -e .")
    sys.exit(1)


def example_single_alignment():
    """Example: Align a single read to a single reference (auto-detect GPU/CPU)"""
    print("=" * 70)
    print("Example 1: Single Read Alignment (RNA002 model, k=5)")
    print("=" * 70)

    # Create a synthetic RNA signal (simulating nanopore current)
    # In real usage, this would be raw signal from fast5/slow5 files
    signal_length = 10000
    np.random.seed(42)
    signal = np.random.randn(signal_length).astype(np.float32) * 10 + 120

    # Define read and reference
    read_name = "test_read_001"
    ref_sequence = "ACGUACGUACGU"  # 12 bases -> 8 kmers for k=5
    ref_name = "chr1"
    sample_rate = 4000.0  # Hz, typical for ONT data

    print("\nInput:")
    print(f"  Signal length: {len(signal)} samples")
    print(f"  Signal stats: mean={signal.mean():.2f}, std={signal.std():.2f}")
    print(f"  Read name: {read_name}")
    print(f"  Reference sequence: {ref_sequence} (length={len(ref_sequence)})")
    print(f"  Reference name: {ref_name}")
    print(f"  Sample rate: {sample_rate} Hz")

    # Method 1: Using the convenience function (auto-detects GPU)
    print("\n--- Method 1: Using convenience function (auto GPU detection) ---")
    result = align_single(signal, read_name, ref_sequence, ref_name,
                         model=MODEL_RNA002, sample_rate=sample_rate)

    print_result(result)

    # Method 2: Using the EventAligner class directly
    print("\n--- Method 2: Using EventAligner class (auto GPU detection) ---")
    aligner = EventAligner(model=MODEL_RNA002)  # Auto-detect GPU
    print(f"  Using GPU: {aligner.using_gpu}")

    result2 = aligner.align_read_single_ref(signal, read_name, ref_sequence,
                                           ref_name, sample_rate)

    print_result(result2)


def example_force_cpu_gpu():
    """Example: Force CPU or GPU usage"""
    print("\n" + "=" * 70)
    print("Example 2: Force CPU or GPU Backend")
    print("=" * 70)

    # Create synthetic signal
    signal_length = 10000
    np.random.seed(43)
    signal = np.random.randn(signal_length).astype(np.float32) * 10 + 120

    read_name = "test_read_002"
    ref_sequence = "ACGUACGUACGU"
    ref_name = "chr1"

    print("\nTrying to force GPU usage...")
    try:
        aligner_gpu = EventAligner(model=MODEL_RNA002, use_gpu=True)
        print(f"  ✓ GPU backend available: {aligner_gpu.using_gpu}")
        result_gpu = aligner_gpu.align_read_single_ref(signal, read_name,
                                                      ref_sequence, ref_name)
        print(f"  Alignments: {result_gpu['n_alignments']}")
    except ImportError as e:
        print(f"  ✗ GPU backend not available: {e}")

    print("\nForcing CPU usage...")
    try:
        aligner_cpu = EventAligner(model=MODEL_RNA002, use_gpu=False)
        print(f"  ✓ CPU backend available: {not aligner_cpu.using_gpu}")
        result_cpu = aligner_cpu.align_read_single_ref(signal, read_name,
                                                      ref_sequence, ref_name)
        print(f"  Alignments: {result_cpu['n_alignments']}")
    except ImportError as e:
        print(f"  ✗ CPU backend not available: {e}")


def example_single_alignment_rna004():
    """Example: Align with RNA004 model (k=9)"""
    print("\n" + "=" * 70)
    print("Example 3: Single Read Alignment (RNA004 model, k=9)")
    print("=" * 70)

    # RNA004 uses k=9, so we need a longer sequence
    signal_length = 15000
    np.random.seed(44)
    signal = np.random.randn(signal_length).astype(np.float32) * 10 + 120

    read_name = "test_read_003"
    ref_sequence = "ACGUACGUACGUACGUACGU"  # 20 bases -> 12 kmers for k=9
    ref_name = "chr2"

    print("\nInput:")
    print(f"  Signal length: {len(signal)} samples")
    print(f"  Read name: {read_name}")
    print(f"  Reference sequence: {ref_sequence} (length={len(ref_sequence)})")
    print("  Model: RNA004 (k=9)")

    aligner = EventAligner(model=MODEL_RNA004)
    print(f"  Using GPU: {aligner.using_gpu}")
    print(f"  K-mer size: {aligner.kmer_size}")

    result = aligner.align_read_single_ref(signal, read_name, ref_sequence,
                                           ref_name)

    print_result(result)


def example_batch_alignment():
    """Example: Batch align multiple reads to multiple references"""
    print("\n" + "=" * 70)
    print("Example 4: Batch Alignment")
    print("=" * 70)

    # Create multiple synthetic signals
    np.random.seed(45)
    n_reads = 3
    n_refs = 2

    signals = [np.random.randn(10000 + i * 1000).astype(np.float32) * 10 + 120
               for i in range(n_reads)]
    read_names = [f"read_{i:03d}" for i in range(n_reads)]

    ref_sequences = ["ACGUACGUACGU", "UGCAUGCAUGCA"]
    ref_names = ["chr1", "chr2"]
    ref_lengths = [len(seq) for seq in ref_sequences]

    print("\nInput:")
    print(f"  Number of reads: {n_reads}")
    print(f"  Number of references: {n_refs}")
    print(f"  Total alignments: {n_reads * n_refs}")

    # Run batch alignment with auto GPU detection
    aligner = EventAligner(model=MODEL_RNA002)
    print(f"  Using GPU: {aligner.using_gpu}")

    results = aligner.align_batch(signals, read_names, ref_sequences,
                                 ref_names, ref_lengths)

    print("\nResults:")
    for key, result in results.items():
        read_name, ref_name = key
        print(f"  ({read_name}, {ref_name}): "
              f"success={result['success']}, "
              f"n_alignments={result['n_alignments']}, "
              f"events_per_base={result['events_per_base']:.2f}")


def print_result(result):
    """Helper function to print alignment result"""
    print("\nResult:")
    print(f"  Success: {result['success']}")
    print(f"  Events detected: {result['n_events']}")
    print(f"  Aligned pairs: {result['n_alignments']}")
    print(f"  Events per base: {result['events_per_base']:.2f}")

    if result['n_alignments'] > 0:
        ref_pos = result['ref_positions']
        read_pos = result['read_positions']
        print("  First 5 alignments:")
        for i in range(min(5, len(ref_pos))):
            print(f"    Ref pos {ref_pos[i]:3d} <-> Event {read_pos[i]:3d}")


def main():
    """Main function"""
    print("\n" + "=" * 70)
    print("FIN EventAlign Wrapper Examples")
    print("=" * 70)

    # Run examples
    try:
        example_single_alignment()
        example_force_cpu_gpu()
        example_single_alignment_rna004()
        example_batch_alignment()

        print("\n" + "=" * 70)
        print("All examples completed successfully!")
        print("=" * 70)

        print("\nUsage summary:")
        print("  1. Import: from fin._eventalign import align_single, align_batch")
        print("  2. Single alignment: align_single(signal, read_name, ref_seq, ref_name)")
        print("  3. Batch alignment: align_batch(signals, read_names, ref_seqs, ref_names, lengths)")
        print("  4. Create aligner: aligner = EventAligner(model=1)")
        print("  5. Force GPU: aligner = EventAligner(model=1, use_gpu=True)")
        print("  6. Force CPU: aligner = EventAligner(model=1, use_gpu=False)")

    except Exception as e:
        print(f"\n✗ Error running examples: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
