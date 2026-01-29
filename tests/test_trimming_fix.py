"""
Test to demonstrate the double-trimming fix

This script shows that MAD-based trimming now works correctly
without losing additional signal due to fixed trimming.
"""

import numpy as np
import matplotlib.pyplot as plt

try:
    from fin._f5c._event import detect_events
except ImportError:
    print("Error: f5c event detection not available")
    print("Build with: python setup.py build_ext --inplace")
    import sys

    sys.exit(1)


def test_adapter_trimming():
    """Test that adapter trimming works correctly for different adapter lengths"""

    print("=" * 80)
    print("Testing Adapter Trimming Fix".center(80))
    print("=" * 80)

    # Test cases with different adapter lengths
    test_cases = [
        ("Short adapter (150 samples)", 150),
        ("Medium adapter (250 samples)", 250),
        ("Long adapter (400 samples)", 400),
    ]

    fig, axes = plt.subplots(len(test_cases), 2, figsize=(14, 10))
    fig.suptitle(
        "Adapter Trimming Test: MAD-Based Method (No Double Trimming)",
        fontsize=14,
        fontweight="bold",
    )

    for idx, (label, adapter_length) in enumerate(test_cases):
        print(f"\n[Test {idx+1}] {label}")
        print("-" * 80)

        # Generate synthetic signal
        # Adapter: low variation (stable signal)
        adapter = np.random.normal(loc=150.0, scale=2.0, size=adapter_length).astype(np.float32)

        # Good DNA signal: high variation (different k-mers)
        dna_signal_length = 2000
        dna_signal = np.random.normal(loc=100.0, scale=10.0, size=dna_signal_length).astype(
            np.float32
        )

        # Tail: low variation
        tail = np.random.normal(loc=140.0, scale=3.0, size=100).astype(np.float32)

        # Combine
        full_signal = np.concatenate([adapter, dna_signal, tail])

        print(f"  Total signal length: {len(full_signal)} samples")
        print(f"  Adapter: 0 to {adapter_length}")
        print(f"  DNA signal: {adapter_length} to {adapter_length + dna_signal_length}")
        print(f"  Tail: {adapter_length + dna_signal_length} to {len(full_signal)}")

        # Detect events
        try:
            events = detect_events(full_signal)  # RNA-only mode

            if len(events) > 0:
                first_event_start = events[0]["start"]
                last_event_end = events[-1]["start"] + events[-1]["length"]

                print(f"\n  Events detected: {len(events)}")
                print(f"  First event starts at sample: {first_event_start}")
                print(f"  Last event ends at sample: {last_event_end}")
                print(f"  Trimmed from start: ~{first_event_start} samples")
                print(f"  Trimmed from end: ~{len(full_signal) - last_event_end} samples")

                # Check if trimming was appropriate
                expected_trim = adapter_length
                actual_trim = first_event_start

                if abs(actual_trim - expected_trim) < 200:  # Within 200 samples tolerance
                    status = "✓ GOOD"
                    color = "green"
                else:
                    status = "✗ BAD (likely double-trimming)"
                    color = "red"

                print(f"\n  Status: {status}")
                print(f"  Expected trimming: ~{expected_trim} samples")
                print(f"  Actual trimming: {actual_trim} samples")
                print(f"  Difference: {abs(actual_trim - expected_trim)} samples")

                # Plot 1: Full signal with regions marked
                ax1 = axes[idx, 0]
                ax1.plot(full_signal, linewidth=0.5, alpha=0.7, color="darkblue")
                ax1.axvspan(0, adapter_length, alpha=0.3, color="red", label="Adapter")
                ax1.axvspan(
                    adapter_length,
                    adapter_length + dna_signal_length,
                    alpha=0.2,
                    color="green",
                    label="DNA signal",
                )
                ax1.axvspan(
                    adapter_length + dna_signal_length,
                    len(full_signal),
                    alpha=0.3,
                    color="orange",
                    label="Tail",
                )

                # Mark first event
                ax1.axvline(
                    first_event_start,
                    color=color,
                    linestyle="--",
                    linewidth=2,
                    label=f"First event (sample {first_event_start})",
                )

                ax1.set_xlabel("Sample Index")
                ax1.set_ylabel("Current (pA)")
                ax1.set_title(f"{label} - Full Signal")
                ax1.legend(loc="upper right", fontsize=8)
                ax1.grid(True, alpha=0.3)

                # Plot 2: Zoomed to trim boundary
                ax2 = axes[idx, 1]
                zoom_start = max(0, adapter_length - 100)
                zoom_end = min(len(full_signal), adapter_length + 300)
                ax2.plot(
                    range(zoom_start, zoom_end),
                    full_signal[zoom_start:zoom_end],
                    linewidth=0.8,
                    color="darkblue",
                )
                ax2.axvline(
                    adapter_length,
                    color="purple",
                    linestyle=":",
                    linewidth=2,
                    label=f"True adapter end ({adapter_length})",
                )
                ax2.axvline(
                    first_event_start,
                    color=color,
                    linestyle="--",
                    linewidth=2,
                    label=f"First event ({first_event_start})",
                )

                # Shade the adapter region
                ax2.axvspan(zoom_start, adapter_length, alpha=0.3, color="red")

                ax2.set_xlabel("Sample Index")
                ax2.set_ylabel("Current (pA)")
                ax2.set_title(f"{label} - Trim Boundary (Zoomed)")
                ax2.legend(fontsize=8)
                ax2.grid(True, alpha=0.3)

            else:
                print("  ✗ ERROR: No events detected")
                axes[idx, 0].text(
                    0.5, 0.5, "No events detected", ha="center", va="center", fontsize=12
                )
                axes[idx, 1].text(
                    0.5, 0.5, "No events detected", ha="center", va="center", fontsize=12
                )

        except Exception as e:
            print(f"  ✗ ERROR: {e}")
            import traceback

            traceback.print_exc()
            axes[idx, 0].text(
                0.5, 0.5, f"Error: {e}", ha="center", va="center", fontsize=10, color="red"
            )
            axes[idx, 1].text(
                0.5, 0.5, f"Error: {e}", ha="center", va="center", fontsize=10, color="red"
            )

    plt.tight_layout()
    plt.savefig("adapter_trimming_test.png", dpi=150, bbox_inches="tight")
    print("\n" + "=" * 80)
    print(f"Figure saved to: adapter_trimming_test.png")
    print("=" * 80)
    print("\nExpected results:")
    print("  ✓ First event should start close to adapter end (within ~100 samples)")
    print("  ✓ No 200-sample fixed offset should be visible")
    print("  ✓ MAD-based trimming adapts to actual adapter length")

    plt.show()


if __name__ == "__main__":
    test_adapter_trimming()
