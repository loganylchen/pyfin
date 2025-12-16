#!/usr/bin/env python3
"""
Test Profile HMM eventalign (full f5c version)
Compare with simple ABEA alignment
"""

import numpy as np
import sys
from pathlib import Path

# Add parent directory to path if needed
parent_dir = Path(__file__).parent.parent
if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

try:
    from fin._f5c import detect_events, eventalign, profile_hmm_eventalign
except ImportError as e:
    print(f"Error importing eventalign extensions: {e}")
    print("Make sure to build the extensions first:")
    print(f"  cd {parent_dir}")
    print("  pip install -e .")
    print("Or try: python setup.py build_ext --inplace")
    sys.exit(1)


def test_profile_hmm():
    """Test Profile HMM eventalign with synthetic data"""

    # Simple test sequence (RNA)
    sequence = "ACGUACGUACGUACGUACGU"

    # Generate synthetic signal
    # For RNA, typical current levels ~100-200 pA
    np.random.seed(42)
    n_events = 60  # ~3 events per base
    signal = np.random.randn(n_events * 50).astype(np.float32) * 10 + 120

    print("=" * 70)
    print("Profile HMM Eventalign Test (Full f5c Version)")
    print("=" * 70)
    print(f"Sequence: {sequence}")
    print(f"Length: {len(sequence)} bases")
    print(f"Signal samples: {len(signal)}")
    print(f"Expected events: ~{(len(sequence) - 4) * 3}")
    print()

    # Test 1: Simple ABEA alignment (original)
    print("Test 1: Simple ABEA Alignment")
    print("-" * 70)
    try:
        result_abea = eventalign(raw_signal=signal, sequence=sequence, kmer_size=5)
        print(f"✓ ABEA alignment successful")
        print(f"  Events detected: {result_abea['n_events']}")
        print(f"  Aligned pairs: {result_abea['n_aligned_pairs']}")
        print(
            f"  Scaling: scale={result_abea['scaling']['scale']:.3f}, shift={result_abea['scaling']['shift']:.1f}"
        )

        # Show base_to_event_map sample
        print(f"  Base-to-event map (first 5):")
        for i, mapping in enumerate(result_abea["base_to_event_map"][:5]):
            print(f"    {i}: kmer={mapping['kmer']}, events={mapping['start']}-{mapping['stop']}")
    except Exception as e:
        print(f"✗ ABEA alignment failed: {e}")
        import traceback

        traceback.print_exc()

    print()

    # Test 2: Profile HMM alignment (full f5c)
    print("Test 2: Profile HMM Alignment (Full f5c)")
    print("-" * 70)
    try:
        result_hmm = profile_hmm_eventalign(
            raw_signal=signal, sequence=sequence, kmer_size=5, events_per_base=3.0
        )
        print(f"✓ Profile HMM alignment successful")
        print(f"  Events detected: {result_hmm['n_events']}")
        print(f"  Aligned records: {result_hmm['n_aligned']}")
        print(f"  Events per base: {result_hmm['events_per_base']:.2f}")
        print(
            f"  Scaling: scale={result_hmm['scaling']['scale']:.3f}, shift={result_hmm['scaling']['shift']:.1f}"
        )

        # Show alignment sample
        print(f"  Alignment records (first 10):")
        print(
            f"    {'Pos':<4} {'Kmer':<6} {'Evt':<4} {'State':<5} {'EvtMean':<8} {'ModelMean':<10} {'ScaledMean':<11}"
        )
        for i, aln in enumerate(result_hmm["alignment"][:10]):
            evt_str = str(aln["event_idx"]) if aln["event_idx"] >= 0 else "-"
            print(
                f"    {aln['ref_position']:<4} {aln['ref_kmer']:<6} {evt_str:<4} "
                f"{aln['hmm_state']:<5} {aln['event_mean']:<8.1f} "
                f"{aln['model_mean']:<10.3f} {aln['scaled_model_mean']:<11.1f}"
            )

        # Count HMM states
        states = {}
        for aln in result_hmm["alignment"]:
            state = aln["hmm_state"]
            states[state] = states.get(state, 0) + 1

        print(f"  HMM state distribution:")
        for state, count in sorted(states.items()):
            state_name = {"M": "MATCH", "K": "KMER_SKIP", "B": "BAD_EVENT"}.get(state, state)
            print(f"    {state_name}: {count}")

    except Exception as e:
        print(f"✗ Profile HMM alignment failed: {e}")
        import traceback

        traceback.print_exc()

    print()
    print("=" * 70)
    print("Test Complete")
    print("=" * 70)
    print()
    print("Key Differences:")
    print("  ABEA (eventalign):")
    print("    - Simple adaptive banded alignment")
    print("    - Returns (event_idx, kmer_idx) pairs")
    print("    - Fast but less detailed")
    print()
    print("  Profile HMM (profile_hmm_eventalign):")
    print("    - Full Viterbi HMM alignment")
    print("    - Returns event_alignment_t with HMM states")
    print("    - Matches f5c eventalign output")
    print("    - Includes MATCH (M), BAD_EVENT (B), KMER_SKIP (K) states")


if __name__ == "__main__":
    test_profile_hmm()
