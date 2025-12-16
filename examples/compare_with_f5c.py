#!/usr/bin/env python3
"""
Compare PyFin eventalign results with f5c reference implementation

This script helps debug alignment issues by:
1. Running both PyFin and f5c on the same data
2. Comparing the results in detail
3. Identifying specific differences

Usage:
    python compare_with_f5c.py --fast5 <file> --fastq <file> --reference <file>

Or for synthetic data:
    python compare_with_f5c.py --synthetic
"""

import numpy as np
import argparse
import subprocess
import tempfile
import os
from pathlib import Path

try:
    from fin._f5c._event import detect_events
    from fin._f5c._eventalign import eventalign
except ImportError:
    print("Error: f5c extensions not available")
    print("Build with: python setup.py build_ext --inplace")
    import sys

    sys.exit(1)


def generate_simple_test_sequence():
    """Generate a simple test sequence for debugging"""
    # Use a sequence with clear patterns
    return "AAAACCCGGGTTT" * 10  # Repetitive pattern for easy debugging


def generate_test_signal(sequence, events_per_base=10):
    """
    Generate very simple synthetic signal for testing
    Each base has a distinct level, minimal noise
    """
    base_levels = {"A": 100.0, "C": 95.0, "G": 105.0, "T": 90.0}

    signal = []
    # Small adapter
    signal.extend(np.random.normal(80, 5, 100))

    # Clear signal for each base
    for base in sequence:
        level = base_levels.get(base, 100.0)
        # Fixed number of samples per base for predictability
        base_signal = np.random.normal(level, 1.0, events_per_base)
        signal.extend(base_signal)

    # Small tail
    signal.extend(np.random.normal(85, 5, 100))

    return np.array(signal, dtype=np.float32)


def analyze_pyfin_alignment(raw_signal, sequence, kmer_size=5):
    """Run PyFin eventalign and analyze results (RNA-only mode)"""
    print("\n" + "=" * 80)
    print("PYFIN ALIGNMENT ANALYSIS (RNA-only mode)")
    print("=" * 80)

    # Detect events first to see what we're working with
    print("\n[1] Event Detection:")
    events = detect_events(raw_signal)  # RNA-only mode
    print(f"    Total events detected: {len(events)}")
    print(f"    Signal length: {len(raw_signal)} samples")
    print(f"    Events per sample: {len(events) / len(raw_signal):.4f}")

    # Show first few events
    print(f"\n    First 5 events:")
    for i, e in enumerate(events[:5]):
        print(
            f"      Event {i}: mean={e['mean']:.2f}, stdv={e['stdv']:.2f}, "
            f"start={e['start']}, length={e['length']}"
        )

    # Run alignment
    print(f"\n[2] Event Alignment:")
    print(f"    Sequence length: {len(sequence)} bases")
    print(f"    K-mer size: {kmer_size}")
    print(f"    Expected k-mers: {len(sequence) - kmer_size + 1}")

    result = eventalign(raw_signal, sequence, kmer_size=kmer_size)  # RNA-only

    print(f"\n    Results:")
    print(f"      Total events: {result['n_events']}")
    print(f"      Aligned pairs: {result['n_aligned_pairs']}")
    print(f"      Alignment rate: {result['n_aligned_pairs']/result['n_events']*100:.1f}%")
    print(f"      Scale: {result['scaling']['scale']:.4f}")
    print(f"      Shift: {result['scaling']['shift']:.2f}")

    # Analyze base-to-event mapping
    print(f"\n[3] Base-to-Event Mapping Analysis:")
    base_to_event_map = result["base_to_event_map"]

    # Check how many k-mers have events
    kmers_with_events = sum(1 for mapping in base_to_event_map if mapping["start"] != -1)
    print(f"    K-mers with events: {kmers_with_events}/{len(base_to_event_map)}")
    print(f"    K-mers without events: {len(base_to_event_map) - kmers_with_events}")

    # Count events per k-mer
    events_per_kmer = []
    for mapping in base_to_event_map:
        if mapping["start"] != -1 and mapping["stop"] != -1:
            n = mapping["stop"] - mapping["start"] + 1
            events_per_kmer.append(n)

    if events_per_kmer:
        print(
            f"    Events per k-mer (aligned): {np.mean(events_per_kmer):.2f} ± {np.std(events_per_kmer):.2f}"
        )
        print(f"    Min: {np.min(events_per_kmer)}, Max: {np.max(events_per_kmer)}")

    # Show detailed mapping for first 10 k-mers
    print(f"\n[4] Detailed Mapping (first 10 k-mers):")
    print(
        f"    {'Pos':<5} {'K-mer':<10} {'Start':<8} {'Stop':<8} {'N Events':<10} {'Expected':<10}"
    )
    print(f"    {'-'*60}")

    expected_events_per_kmer = len(events) / len(base_to_event_map)

    for i in range(min(10, len(base_to_event_map))):
        mapping = base_to_event_map[i]
        kmer = mapping["kmer"]
        start = mapping["start"]
        stop = mapping["stop"]

        if start != -1 and stop != -1:
            n_events = stop - start + 1
            status = (
                "✓" if abs(n_events - expected_events_per_kmer) < expected_events_per_kmer else "⚠"
            )
        else:
            n_events = 0
            status = "✗"

        print(
            f"    {i:<5} {kmer:<10} {start:<8} {stop:<8} {n_events:<10} "
            f"{expected_events_per_kmer:.1f} {status}"
        )

    # Check for issues
    print(f"\n[5] Issue Detection:")
    issues = []

    if result["n_aligned_pairs"] < result["n_events"] * 0.5:
        issues.append("⚠ Less than 50% of events aligned - possible alignment failure")

    if kmers_with_events < len(base_to_event_map) * 0.8:
        issues.append(
            f"⚠ Only {kmers_with_events}/{len(base_to_event_map)} k-mers have events - possible gaps"
        )

    if events_per_kmer and np.std(events_per_kmer) > np.mean(events_per_kmer):
        issues.append("⚠ High variance in events per k-mer - uneven alignment")

    # Check for consecutive k-mers without events
    no_event_runs = []
    current_run = 0
    for mapping in base_to_event_map:
        if mapping["start"] == -1:
            current_run += 1
        else:
            if current_run > 0:
                no_event_runs.append(current_run)
            current_run = 0

    if no_event_runs and max(no_event_runs) > 5:
        issues.append(
            f"⚠ Long gaps detected: max {max(no_event_runs)} consecutive k-mers without events"
        )

    if not issues:
        print("    ✓ No major issues detected")
    else:
        for issue in issues:
            print(f"    {issue}")

    return result


def compare_alignments_detailed(pyfin_result, sequence, kmer_size):
    """
    Detailed analysis of alignment quality
    """
    print("\n" + "=" * 80)
    print("ALIGNMENT QUALITY ANALYSIS")
    print("=" * 80)

    base_to_event_map = pyfin_result["base_to_event_map"]
    n_kmers = len(sequence) - kmer_size + 1

    print(f"\n[1] Coverage Statistics:")
    print(f"    Expected k-mers: {n_kmers}")
    print(f"    K-mers in result: {len(base_to_event_map)}")

    if len(base_to_event_map) != n_kmers:
        print(f"    ⚠ WARNING: K-mer count mismatch!")
        print(f"      This suggests the alignment may be truncated or incorrect")

    # Check if k-mers match sequence
    print(f"\n[2] K-mer Sequence Validation:")
    mismatches = 0
    for i in range(min(len(base_to_event_map), n_kmers)):
        expected_kmer = sequence[i : i + kmer_size]
        actual_kmer = base_to_event_map[i]["kmer"]
        if expected_kmer != actual_kmer:
            if mismatches < 5:  # Show first 5 mismatches
                print(f"    ✗ Position {i}: expected '{expected_kmer}', got '{actual_kmer}'")
            mismatches += 1

    if mismatches == 0:
        print(f"    ✓ All k-mers match sequence")
    else:
        print(f"    ✗ {mismatches} k-mer mismatches found!")

    # Check event index ordering
    print(f"\n[3] Event Index Ordering:")
    prev_stop = -1
    ordering_issues = 0
    for i, mapping in enumerate(base_to_event_map):
        if mapping["start"] != -1:
            if prev_stop != -1 and mapping["start"] <= prev_stop:
                if ordering_issues < 3:
                    print(f"    ✗ Position {i}: start={mapping['start']} <= prev_stop={prev_stop}")
                ordering_issues += 1
            prev_stop = mapping["stop"]

    if ordering_issues == 0:
        print(f"    ✓ Event indices are monotonically increasing")
    else:
        print(f"    ✗ {ordering_issues} event index ordering violations!")

    # Distribution analysis
    print(f"\n[4] Distribution Analysis:")
    event_counts = []
    for mapping in base_to_event_map:
        if mapping["start"] != -1 and mapping["stop"] != -1:
            event_counts.append(mapping["stop"] - mapping["start"] + 1)

    if event_counts:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(12, 4))

        plt.subplot(131)
        plt.hist(event_counts, bins=20, edgecolor="black", alpha=0.7)
        plt.xlabel("Events per K-mer")
        plt.ylabel("Count")
        plt.title("Distribution of Events per K-mer")
        plt.axvline(
            np.mean(event_counts),
            color="red",
            linestyle="--",
            label=f"Mean: {np.mean(event_counts):.1f}",
        )
        plt.legend()

        plt.subplot(132)
        plt.plot(event_counts, marker="o", markersize=2, alpha=0.6)
        plt.xlabel("K-mer Position")
        plt.ylabel("Number of Events")
        plt.title("Events per K-mer Across Sequence")
        plt.axhline(np.mean(event_counts), color="red", linestyle="--", alpha=0.5)

        plt.subplot(133)
        coverage = [1 if m["start"] != -1 else 0 for m in base_to_event_map]
        plt.plot(np.cumsum(coverage) / np.arange(1, len(coverage) + 1) * 100, linewidth=2)
        plt.xlabel("K-mer Position")
        plt.ylabel("Coverage (%)")
        plt.title("Cumulative K-mer Coverage")
        plt.ylim([0, 105])
        plt.axhline(100, color="green", linestyle="--", alpha=0.3)

        plt.tight_layout()
        plt.savefig("alignment_analysis.png", dpi=150, bbox_inches="tight")
        print(f"    ✓ Distribution plots saved to 'alignment_analysis.png'")
        plt.show()


def test_synthetic_alignment():
    """Test alignment with simple synthetic data"""
    print("\n" + "=" * 80)
    print("SYNTHETIC DATA TEST")
    print("=" * 80)

    # Create simple test case
    sequence = generate_simple_test_sequence()
    print(f"\nSequence: {sequence[:50]}... (length: {len(sequence)})")

    # Generate signal
    raw_signal = generate_test_signal(sequence, events_per_base=10)
    print(f"Signal: {len(raw_signal)} samples")

    # Test with RNA-only mode
    print(f"\n--- Testing RNA alignment ---")
    result = analyze_pyfin_alignment(raw_signal, sequence, kmer_size=5)

    # Detailed comparison
    compare_alignments_detailed(result, sequence, kmer_size=5)


def main():
    parser = argparse.ArgumentParser(description="Compare PyFin with f5c alignment")
    parser.add_argument("--synthetic", action="store_true", help="Test with synthetic data")
    parser.add_argument("--fast5", type=str, help="FAST5 file for real data test")
    parser.add_argument("--reference", type=str, help="Reference sequence file (FASTA)")

    args = parser.parse_args()

    if args.synthetic:
        test_synthetic_alignment()
    else:
        print("Real data comparison not yet implemented.")
        print("Use --synthetic flag for synthetic data test")


if __name__ == "__main__":
    main()
