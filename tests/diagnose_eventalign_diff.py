#!/usr/bin/env python3
"""
Diagnostic script to compare pyfin eventalign vs f5c eventalign results.

This script helps identify why eventalign results differ from f5c by:
1. Checking event ordering (reversed vs raw)
2. Comparing signal position mappings
3. Analyzing k-mer alignment differences
"""

import numpy as np
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    from fin._f5c import detect_events, profile_hmm_eventalign
    from fin.io import Pod5Reader

    print("✓ Loaded pyfin eventalign modules")
except ImportError as e:
    print(f"✗ Failed to load pyfin modules: {e}")
    sys.exit(1)


def diagnose_event_order(signal: np.ndarray):
    """
    Diagnose event ordering by checking event detection.

    Issues to check:
    1. Are events in raw signal order or reversed?
    2. Does event[0] correspond to start or end of signal?
    """
    print("\n" + "=" * 70)
    print("DIAGNOSTIC 1: Event Ordering")
    print("=" * 70)

    # Detect events
    events = detect_events(signal)
    n_events = len(events)

    print(f"\nDetected {n_events} events from signal of length {len(signal)}")
    print("\nFirst 3 events (should be at START of raw signal):")
    for i in range(min(3, n_events)):
        evt = events[i]
        print(
            f"  Event {i}: start={evt['start']}, length={evt['length']:.1f}, mean={evt['mean']:.2f}"
        )

    print("\nLast 3 events (should be at END of raw signal):")
    for i in range(max(0, n_events - 3), n_events):
        evt = events[i]
        print(
            f"  Event {i}: start={evt['start']}, length={evt['length']:.1f}, mean={evt['mean']:.2f}"
        )

    # Check: first event should have low start index, last event should have high start
    first_start = events[0]["start"]
    last_start = events[-1]["start"]

    print("\n📊 Event Order Analysis:")
    print(f"  First event start: {first_start}")
    print(f"  Last event start:  {last_start}")
    print(f"  Signal length:     {len(signal)}")

    if first_start < last_start:
        print("  ✓ Events are in RAW SIGNAL order (first event at start of signal)")
        print("    This means: event[0] = 3' end of RNA, event[n-1] = 5' end of RNA")
        return "raw_order"
    else:
        print("  ✗ Events appear REVERSED (first event at end of signal)")
        print("    This means: event[0] = 5' end of RNA, event[n-1] = 3' end of RNA")
        return "reversed"


def diagnose_eventalign_mapping(signal: np.ndarray, sequence: str):
    """
    Diagnose eventalign event-to-kmer mapping.

    Issues to check:
    1. Do alignment positions match raw signal coordinates?
    2. Are event indices consistent with event detection?
    3. Do k-mers map to correct signal regions?
    """
    print("\n" + "=" * 70)
    print("DIAGNOSTIC 2: Eventalign Event Mapping")
    print("=" * 70)

    # Run eventalign
    try:
        result = profile_hmm_eventalign(raw_signal=signal, sequence=sequence, kmer_size=5)
    except Exception as e:
        print(f"✗ Eventalign failed: {e}")
        return None

    alignment = result.get("alignment", [])
    n_aligned = len([a for a in alignment if a.get("event_idx", -1) >= 0])

    print("\nAlignment summary:")
    print(f"  Total events detected: {result.get('n_events', 0)}")
    print(f"  Events aligned: {n_aligned}")
    print(f"  Sequence length: {len(sequence)}")
    print("  K-mer size: 5")
    print(f"  Expected k-mers: {len(sequence) - 4}")

    # Get valid alignments
    valid_alns = [a for a in alignment if a.get("event_idx", -1) >= 0]

    if not valid_alns:
        print("✗ No valid alignments found")
        return None

    # Sort by sequence position
    valid_alns.sort(key=lambda x: x.get("ref_position", 0))

    print("\n📊 First 5 aligned k-mers:")
    for i, aln in enumerate(valid_alns[:5]):
        pos = aln.get("ref_position", -1)
        kmer = aln.get("ref_kmer", "?")
        event_idx = aln.get("event_idx", -1)
        sig_start = aln.get("signal_start", -1)
        sig_len = aln.get("signal_length", 0)

        # Expected kmer from sequence
        if 0 <= pos < len(sequence) - 4:
            expected_kmer = sequence[pos : pos + 5]
        else:
            expected_kmer = "???"

        match = "✓" if kmer == expected_kmer else "✗"
        print(f"  Pos {pos:3d}: {kmer} {match} (expected: {expected_kmer})")
        print(f"           event_idx={event_idx}, signal[{sig_start}:{sig_start+int(sig_len)}]")

    print("\n📊 Last 5 aligned k-mers:")
    for i, aln in enumerate(valid_alns[-5:]):
        pos = aln.get("ref_position", -1)
        kmer = aln.get("ref_kmer", "?")
        event_idx = aln.get("event_idx", -1)
        sig_start = aln.get("signal_start", -1)
        sig_len = aln.get("signal_length", 0)

        if 0 <= pos < len(sequence) - 4:
            expected_kmer = sequence[pos : pos + 5]
        else:
            expected_kmer = "???"

        match = "✓" if kmer == expected_kmer else "✗"
        print(f"  Pos {pos:3d}: {kmer} {match} (expected: {expected_kmer})")
        print(f"           event_idx={event_idx}, signal[{sig_start}:{sig_start+int(sig_len)}]")

    # Check signal position ordering
    print("\n📊 Signal Position Ordering:")
    first_sig_start = valid_alns[0].get("signal_start", 0)
    last_sig_start = valid_alns[-1].get("signal_start", 0)

    print(f"  First k-mer (pos 0, 5' end): signal_start={first_sig_start}")
    print(f"  Last k-mer (3' end):         signal_start={last_sig_start}")

    if first_sig_start > last_sig_start:
        print("  ✓ Signal positions DECREASE from 5'→3' (expected for RNA)")
        print("    This is correct: 5' end aligns to END of signal, 3' end to START")
    else:
        print("  ✗ Signal positions INCREASE from 5'→3' (unexpected)")
        print("    This suggests event reversal may be incorrect")

    return result


def compare_with_f5c(pyfin_result, f5c_records):
    """
    Compare pyfin eventalign with f5c eventalign results.

    Args:
        pyfin_result: Result from profile_hmm_eventalign
        f5c_records: List of f5c eventalign records (from read_f5c_eventalign_table)
    """
    print("\n" + "=" * 70)
    print("DIAGNOSTIC 3: Pyfin vs F5C Comparison")
    print("=" * 70)

    if not pyfin_result or not f5c_records:
        print("✗ Missing data for comparison")
        return

    pyfin_alns = [a for a in pyfin_result.get("alignment", []) if a.get("event_idx", -1) >= 0]

    print("\nAlignment counts:")
    print(f"  Pyfin: {len(pyfin_alns)} events aligned")
    print(f"  F5C:   {len(f5c_records)} events aligned")

    # Sort both by sequence position
    pyfin_alns.sort(key=lambda x: x.get("ref_position", 0))
    f5c_sorted = sorted(f5c_records, key=lambda x: x["position"])

    # Compare first few k-mers
    print("\n📊 Comparing first 5 k-mers:")
    for i in range(min(5, len(pyfin_alns), len(f5c_sorted))):
        pyfin = pyfin_alns[i]
        f5c = f5c_sorted[i]

        pyfin_kmer = pyfin.get("ref_kmer", "?")
        f5c_kmer = f5c["reference_kmer"]
        pyfin_pos = pyfin.get("ref_position", -1)
        f5c_pos = f5c["position"]
        pyfin_sig = pyfin.get("signal_start", -1)
        f5c_sig = f5c["start_idx"]

        kmer_match = "✓" if pyfin_kmer == f5c_kmer else "✗"
        pos_match = "✓" if pyfin_pos == f5c_pos else "✗"
        sig_diff = abs(pyfin_sig - f5c_sig) if pyfin_sig >= 0 else float("inf")
        sig_match = "✓" if sig_diff < 10 else "✗"

        print(f"  Position {i}:")
        print(f"    K-mer: pyfin={pyfin_kmer:5s} f5c={f5c_kmer:5s} {kmer_match}")
        print(f"    Pos:   pyfin={pyfin_pos:4d} f5c={f5c_pos:4d} {pos_match}")
        print(f"    Signal: pyfin={pyfin_sig:5d} f5c={f5c_sig:5d} diff={sig_diff:.0f} {sig_match}")


def main():
    """Run diagnostics on eventalign differences."""
    import argparse

    parser = argparse.ArgumentParser(description="Diagnose eventalign vs f5c differences")
    parser.add_argument("--pod5", required=True, help="POD5 file")
    parser.add_argument("--read-id", required=True, help="Read ID to analyze")
    parser.add_argument("--sequence", required=True, help="Reference sequence (5'->3')")
    parser.add_argument("--f5c-table", help="F5C eventalign table (optional)")

    args = parser.parse_args()

    # Load signal
    print(f"Loading read: {args.read_id}")
    with Pod5Reader(args.pod5) as reader:
        result = reader.get_calibrated_signal(args.read_id)
        if result is None:
            print(f"✗ Read not found: {args.read_id}")
            return 1
        signal, metadata = result
        signal = np.array(signal, dtype=np.float32)

    print(f"✓ Loaded signal: {len(signal)} samples")
    print(f"✓ Sequence length: {len(args.sequence)} bases")

    # Run diagnostics
    event_order = diagnose_event_order(signal)
    alignment_result = diagnose_eventalign_mapping(signal, args.sequence)

    # Compare with f5c if provided
    if args.f5c_table and alignment_result:
        try:
            import gzip
            import csv

            # Load f5c data
            is_gzipped = args.f5c_table.endswith(".gz")
            open_func = gzip.open if is_gzipped else open
            mode = "rt" if is_gzipped else "r"

            f5c_records = []
            with open_func(args.f5c_table, mode) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    if row["read_name"] == args.read_id:
                        f5c_records.append(
                            {
                                "position": int(row["position"]),
                                "reference_kmer": row["reference_kmer"],
                                "start_idx": int(row["start_idx"]),
                                "end_idx": int(row["end_idx"]),
                                "event_index": int(row["event_index"]),
                                "event_level_mean": float(row["event_level_mean"]),
                            }
                        )

            print(f"\n✓ Loaded {len(f5c_records)} f5c records")
            compare_with_f5c(alignment_result, f5c_records)

        except Exception as e:
            print(f"✗ Failed to load f5c data: {e}")

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"\nEvent ordering: {event_order}")
    print("Check the analysis above for specific discrepancies.")
    print("\nSee EVENTALIGN_F5C_DIFF_ANALYSIS.md for detailed explanation.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
