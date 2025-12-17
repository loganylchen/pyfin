#!/usr/bin/env python3
"""
Compare adapter trimming and soft-clipping parameters between pyfin and f5c.

This helps identify why more events are skipped at different ends.
"""

import numpy as np
import sys
from pathlib import Path

print("=" * 70)
print("ADAPTER TRIMMING & SOFT-CLIPPING COMPARISON")
print("=" * 70)

print("\n📊 PYFIN Current Settings (event_detection_simple.c):")
print("  trim_start = 200 samples  (removes from START of signal = 3' end of RNA)")
print("  trim_end   = 10 samples   (removes from END of signal = 5' end of RNA)")
print("  varseg_chunk = 100")
print("  varseg_thresh = 0.0")
print("")
print("  Soft-clipping (align_common.h):")
print("  TRANS_START_TO_CLIP = 0.5  (probability of entering clip state)")
print("  TRANS_CLIP_SELF = 0.9      (probability of staying in clip state)")

print("\n📊 F5C Default Settings (from f5c source):")
print("  trim_start = 200 samples  (removes from START of signal = 3' end of RNA)")
print("  trim_end   = 10 samples   (removes from END of signal = 5' end of RNA)")
print("  varseg_chunk = 100")
print("  varseg_thresh = 0.0")
print("")
print("  Soft-clipping:")
print("  TRANS_START_TO_CLIP = 0.5")
print("  TRANS_CLIP_SELF = 0.9")

print("\n" + "=" * 70)
print("ANALYSIS")
print("=" * 70)

print("\n❓ Your Observation:")
print("  F5C: Ignores more events on 3' end (start of signal)")
print("  Pyfin: Skips more events on 5' end (end of signal)")

print("\n🔍 But adapter trimming is SYMMETRIC:")
print("  Both trim 200 samples from start (3' end)")
print("  Both trim 10 samples from end (5' end)")

print("\n💡 This suggests the difference is in HMM alignment, not adapter trimming!")

print("\nPossible causes:")
print("  1. Event reversal before/after alignment affects which end is 'start'")
print("  2. Soft-clipping applied differently at pre- vs post-flanking")
print("  3. HMM initialization favors one end over the other")
print("  4. Event quality filtering affects ends differently")

print("\n📋 To diagnose:")
print("  1. Check which events are marked as 'aligned' vs 'skipped'")
print("  2. Look at event_idx values at both ends")
print("  3. Compare signal_start positions for first/last aligned k-mers")
print("  4. Check if ref_position coverage differs at ends")

print("\n" + "=" * 70)
print("COORDINATE SYSTEMS")
print("=" * 70)

print("\nFor RNA Direct Sequencing:")
print("  Raw signal: [0 .................... N]")
print("              3' end            5' end")
print("              (enters pore)     (exits pore)")
print("")
print("  After event detection (reversed):")
print("  Events:     [0 .................... M]")
print("              5' end            3' end")
print("              (sequence[0])     (sequence[L-1])")
print("")
print("  Sequence:   5' ACGUA...CGAU 3'")
print("              [0 ........ L-1]")

print("\n🎯 Key Question:")
print("  When eventalign outputs event_idx, is it:")
print("    A) Index in REVERSED event array (matches sequence position)?")
print("    B) Index in RAW signal order (matches signal coordinates)?")
print("")
print("  And when converting back to signal_start:")
print("    - Does it correctly map from reversed → raw coordinates?")

print("\n" + "=" * 70)
print("RECOMMENDED DIAGNOSTICS")
print("=" * 70)

print("\n1. Run pyfin eventalign and check first/last aligned positions:")
print("   python examples/diagnose_eventalign_diff.py \\")
print("     --pod5 <file> --read-id <id> --sequence <seq> --f5c-table <table>")

print("\n2. Check event reversal in eventalign.c:")
print("   - Line 154: event_table et = getevents_simple(nsample, rawptr);")
print("   - Are events reversed BEFORE alignment?")
print("   - Line 263: int raw_event_idx = (int)et.n - 1 - internal_event_idx;")
print("   - Is this converting correctly?")

print("\n3. Compare alignment ranges:")
print("   - Pyfin: min/max ref_position in alignment")
print("   - F5C: min/max position in TSV output")
print("   - Which end has more unaligned k-mers?")

print("\n" + "=" * 70)
