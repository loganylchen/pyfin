#!/usr/bin/env python3
"""
Quick debug script to check eventalign results

This will show you exactly what's coming out of eventalign
and help identify issues.
"""

import numpy as np

try:
    from fin._f5c._event import detect_events
    from fin._f5c._eventalign import eventalign
except ImportError:
    print("Error: f5c extensions not available")
    print("Build with: python setup.py build_ext --inplace")
    import sys

    sys.exit(1)


def debug_eventalign():
    """Run eventalign with simple test data and show detailed output"""

    # Very simple test: short sequence with clear signal
    sequence = "AAACCCGGGTTT"
    print(f"Sequence: {sequence} (length: {len(sequence)})")
    print(
        f"Expected 5-mers: {len(sequence) - 5 + 1} = {', '.join([sequence[i:i+5] for i in range(len(sequence)-4)])}"
    )

    # Generate very clear signal
    base_levels = {"A": 100.0, "C": 95.0, "G": 105.0, "T": 90.0}
    signal = []

    # Small adapter
    signal.extend([80] * 50)

    # 10 samples per base, clear levels
    for base in sequence:
        level = base_levels[base]
        signal.extend([level] * 10)  # Exactly 10 samples, no noise

    # Small tail
    signal.extend([85] * 50)

    signal = np.array(signal, dtype=np.float32)
    print(f"\nSignal: {len(signal)} samples")
    print(f"Expected: {50 + len(sequence)*10 + 50} = {50 + len(sequence)*10 + 50}")

    # Detect events
    print("\n" + "=" * 70)
    print("EVENT DETECTION")
    print("=" * 70)
    events = detect_events(signal)  # RNA-only mode
    print(f"Detected {len(events)} events")
    print("\nFirst 10 events:")
    for i in range(min(10, len(events))):
        e = events[i]
        print(
            f"  Event {i:2d}: mean={e['mean']:6.2f}, stdv={e['stdv']:5.2f}, "
            f"start={e['start']:4d}, length={e['length']:3d}"
        )

    # Run alignment
    print("\n" + "=" * 70)
    print("ALIGNMENT")
    print("=" * 70)

    print("\nRunning eventalign...")
    result = eventalign(signal, sequence, kmer_size=5)  # RNA-only mode

    print(f"\nResults:")
    print(f"  n_events: {result['n_events']}")
    print(f"  n_aligned_pairs: {result['n_aligned_pairs']}")
    print(
        f"  scaling: scale={result['scaling']['scale']:.4f}, shift={result['scaling']['shift']:.2f}"
    )

    # Check base_to_event_map
    print("\n" + "=" * 70)
    print("BASE-TO-EVENT MAP")
    print("=" * 70)

    base_to_event_map = result["base_to_event_map"]
    print(f"Length: {len(base_to_event_map)}")
    print(f"Expected: {len(sequence) - 5 + 1}")

    if len(base_to_event_map) == 0:
        print("\n❌ ERROR: base_to_event_map is EMPTY!")
        print("   This means alignment completely failed.")
        return

    print(f"\nDetailed mapping:")
    print(f"{'Idx':<4} {'K-mer':<8} {'Start':<8} {'Stop':<8} {'N Events':<10}")
    print("-" * 50)

    for i, mapping in enumerate(base_to_event_map):
        kmer = mapping["kmer"]
        start = mapping["start"]
        stop = mapping["stop"]
        n_events = stop - start + 1 if start != -1 else 0

        # Expected k-mer
        expected_kmer = sequence[i : i + 5] if i + 5 <= len(sequence) else "?"
        match = "✓" if kmer == expected_kmer else "✗"

        print(
            f"{i:<4} {kmer:<8} {start:<8} {stop:<8} {n_events:<10} {match} (expected: {expected_kmer})"
        )

    # Analysis
    print("\n" + "=" * 70)
    print("ANALYSIS")
    print("=" * 70)

    # Check k-mer correctness
    correct_kmers = sum(
        1 for i, m in enumerate(base_to_event_map) if m["kmer"] == sequence[i : i + 5]
    )
    print(f"K-mer correctness: {correct_kmers}/{len(base_to_event_map)}")

    # Check coverage
    covered = sum(1 for m in base_to_event_map if m["start"] != -1)
    print(
        f"Coverage: {covered}/{len(base_to_event_map)} k-mers have events ({100*covered/len(base_to_event_map):.1f}%)"
    )

    # Check event assignment
    total_assigned_events = sum(
        m["stop"] - m["start"] + 1 for m in base_to_event_map if m["start"] != -1
    )
    print(
        f"Events assigned: {total_assigned_events}/{result['n_events']} ({100*total_assigned_events/result['n_events']:.1f}%)"
    )

    # Check for issues
    print("\n" + "=" * 70)
    print("ISSUES DETECTED")
    print("=" * 70)

    issues = []

    if len(base_to_event_map) != len(sequence) - 5 + 1:
        issues.append(
            f"❌ Wrong number of k-mers: got {len(base_to_event_map)}, expected {len(sequence) - 5 + 1}"
        )

    if correct_kmers != len(base_to_event_map):
        issues.append(
            f"❌ K-mer sequence mismatch: {len(base_to_event_map) - correct_kmers} incorrect"
        )

    if covered < len(base_to_event_map) * 0.8:
        issues.append(
            f"⚠️  Low coverage: only {covered}/{len(base_to_event_map)} k-mers have events"
        )

    if total_assigned_events < result["n_events"] * 0.8:
        issues.append(
            f"⚠️  Many unassigned events: {result['n_events'] - total_assigned_events} events not mapped"
        )

    # Check ordering
    prev_stop = -1
    for i, m in enumerate(base_to_event_map):
        if m["start"] != -1:
            if m["start"] < prev_stop:
                issues.append(
                    f"❌ Event ordering violation at k-mer {i}: start={m['start']} < prev_stop={prev_stop}"
                )
                break
            prev_stop = m["stop"]

    if not issues:
        print("✅ No major issues detected!")
    else:
        for issue in issues:
            print(f"  {issue}")

    return result


if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("EVENTALIGN DEBUG")
    print("=" * 70)
    debug_eventalign()
