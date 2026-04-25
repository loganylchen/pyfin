#!/usr/bin/env python3
"""
Analyze soft-clipping behavior at 3' vs 5' ends.

This script checks which events are being soft-clipped (excluded from alignment)
at each end of the sequence to understand the asymmetry between pyfin and f5c.
"""



def analyze_alignment_coverage(pyfin_alignment, f5c_records, sequence_length):
    """
    Compare which sequence positions are covered by alignments.

    Args:
        pyfin_alignment: List of alignment dicts from pyfin
        f5c_records: List of f5c alignment records
        sequence_length: Length of reference sequence
    """
    print("\n" + "=" * 70)
    print("ALIGNMENT COVERAGE ANALYSIS")
    print("=" * 70)

    # Get covered positions for pyfin
    pyfin_positions = set()
    for aln in pyfin_alignment:
        if aln.get("event_idx", -1) >= 0:
            pos = aln.get("ref_position", -1)
            if pos >= 0:
                pyfin_positions.add(pos)

    # Get covered positions for f5c
    f5c_positions = set()
    for rec in f5c_records:
        pos = rec.get("position", -1)
        if pos >= 0:
            f5c_positions.add(pos)

    print(f"\nSequence length: {sequence_length} bases")
    print(f"Expected k-mers (5-mer): {sequence_length - 4}")

    print("\n📊 Pyfin Coverage:")
    print(f"  Aligned positions: {len(pyfin_positions)}")
    if pyfin_positions:
        print(f"  First aligned: position {min(pyfin_positions)} (5' end)")
        print(f"  Last aligned:  position {max(pyfin_positions)} (3' end)")
        pyfin_5p_missing = min(pyfin_positions)
        pyfin_3p_missing = (sequence_length - 5) - max(pyfin_positions)
        print(f"  Unaligned at 5' end: {pyfin_5p_missing} positions")
        print(f"  Unaligned at 3' end: {pyfin_3p_missing} positions")

    print("\n📊 F5C Coverage:")
    print(f"  Aligned positions: {len(f5c_positions)}")
    if f5c_positions:
        print(f"  First aligned: position {min(f5c_positions)} (5' end)")
        print(f"  Last aligned:  position {max(f5c_positions)} (3' end)")
        f5c_5p_missing = min(f5c_positions)
        f5c_3p_missing = (sequence_length - 5) - max(f5c_positions)
        print(f"  Unaligned at 5' end: {f5c_5p_missing} positions")
        print(f"  Unaligned at 3' end: {f5c_3p_missing} positions")

    if pyfin_positions and f5c_positions:
        print("\n🔍 Difference:")
        print(f"  5' end: pyfin missing {pyfin_5p_missing}, f5c missing {f5c_5p_missing}")
        print(f"          Δ = {pyfin_5p_missing - f5c_5p_missing} (positive = pyfin skips more)")
        print(f"  3' end: pyfin missing {pyfin_3p_missing}, f5c missing {f5c_3p_missing}")
        print(f"          Δ = {pyfin_3p_missing - f5c_3p_missing} (positive = pyfin skips more)")

        if pyfin_5p_missing > f5c_5p_missing:
            print(
                f"\n  ⚠️  Pyfin skips {pyfin_5p_missing - f5c_5p_missing} more positions at 5' end"
            )
            print("      This is the END of raw signal (high event indices)")
            print("      Affected by PRE-flanking soft-clipping")
        elif pyfin_5p_missing < f5c_5p_missing:
            print(f"\n  ⚠️  F5C skips {f5c_5p_missing - pyfin_5p_missing} more positions at 5' end")

        if pyfin_3p_missing > f5c_3p_missing:
            print(
                f"\n  ⚠️  Pyfin skips {pyfin_3p_missing - f5c_3p_missing} more positions at 3' end"
            )
            print("      This is the START of raw signal (low event indices)")
            print("      Affected by POST-flanking soft-clipping")
        elif pyfin_3p_missing < f5c_3p_missing:
            print(f"\n  ⚠️  F5C skips {f5c_3p_missing - pyfin_3p_missing} more positions at 3' end")


def analyze_event_coverage(pyfin_result, f5c_records, n_events_total):
    """
    Compare which events are included in alignments.
    """
    print("\n" + "=" * 70)
    print("EVENT COVERAGE ANALYSIS")
    print("=" * 70)

    # Get pyfin event indices
    pyfin_events = set()
    for aln in pyfin_result.get("alignment", []):
        evt_idx = aln.get("event_idx", -1)
        if evt_idx >= 0:
            pyfin_events.add(evt_idx)

    # Get f5c event indices
    f5c_events = set()
    for rec in f5c_records:
        evt_idx = rec.get("event_index", -1)
        if evt_idx >= 0:
            f5c_events.add(evt_idx)

    print(f"\nTotal events detected: {n_events_total}")

    print("\n📊 Pyfin Event Usage:")
    print(f"  Events aligned: {len(pyfin_events)}")
    if pyfin_events:
        print(f"  Event index range: {min(pyfin_events)} to {max(pyfin_events)}")
        print(f"  Unused at low indices: {min(pyfin_events)}")
        print(f"  Unused at high indices: {n_events_total - 1 - max(pyfin_events)}")

    print("\n📊 F5C Event Usage:")
    print(f"  Events aligned: {len(f5c_events)}")
    if f5c_events:
        print(f"  Event index range: {min(f5c_events)} to {max(f5c_events)}")
        print(f"  Unused at low indices: {min(f5c_events)}")
        print(f"  Unused at high indices: {n_events_total - 1 - max(f5c_events)}")

    print("\n💡 Note: Event indices may be in different coordinate systems!")
    print("   - If events are reversed for alignment, event[0] = 5' end")
    print("   - If events are in raw order, event[0] = 3' end")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Analyze soft-clipping asymmetry")
    parser.add_argument("--pyfin-result", help="Path to saved pyfin alignment result (pickle)")
    parser.add_argument("--f5c-table", help="F5C eventalign TSV file")
    parser.add_argument("--read-id", help="Read ID to analyze")
    parser.add_argument("--sequence", required=True, help="Reference sequence")

    args = parser.parse_args()

    # This is a template - actual usage would load results from files
    print("=" * 70)
    print("SOFT-CLIPPING ASYMMETRY ANALYZER")
    print("=" * 70)

    print("\n📋 To use this script:")
    print("  1. Run pyfin eventalign and save the result")
    print("  2. Run f5c eventalign and save the TSV output")
    print("  3. Pass both results to this script for comparison")

    print("\n📊 What this analyzes:")
    print("  - Which sequence positions are missing from alignment at each end")
    print("  - Which event indices are unused at each end")
    print("  - Whether the asymmetry is in pre-flanking (5' end) or post-flanking (3' end)")

    print("\n🔧 Potential fixes:")
    print("  1. Adjust TRANS_START_TO_CLIP (lower = less clipping)")
    print("  2. Make pre/post-flanking asymmetric (different clip rates)")
    print("  3. Adjust adapter trimming (trim_start vs trim_end)")
    print("  4. Check event reversal is applied consistently")

    print("\n💡 Key insight:")
    print("  For RNA, after event reversal:")
    print("    - PRE-flanking affects 5' end (high signal indices)")
    print("    - POST-flanking affects 3' end (low signal indices)")
    print("  So if pyfin skips more 5' events → increase TRANS_START_TO_CLIP for pre-flanking")
    print("  Or if pyfin skips more 3' events → increase TRANS_START_TO_CLIP for post-flanking")


if __name__ == "__main__":
    main()
