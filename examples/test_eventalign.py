"""
Example: Nanopore Event-to-Sequence Alignment (Eventalign)

This example demonstrates:
1. Generating synthetic nanopore signal for a known sequence
2. Detecting events from the raw signal
3. Aligning events to the reference sequence (eventalign)
4. Visualizing the alignment results
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

try:
    from fin._f5c._event import detect_events
    from fin._f5c._eventalign import eventalign
except ImportError:
    print("Error: f5c extensions not available")
    print("Build with: python setup.py build_ext --inplace")
    import sys

    sys.exit(1)
import pickle


def visualize_eventalign(raw_signal, sequence, result, kmer_size=5):
    """
    Create comprehensive visualization of eventalign results.

    Args:
        raw_signal: Raw nanopore signal
        sequence: Reference sequence
        result: Eventalign result dictionary
        kmer_size: K-mer size used for alignment
    """
    base_to_event_map = result["base_to_event_map"]
    scaling = result["scaling"]

    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(4, 2, hspace=0.3, wspace=0.3)

    # --- Panel 1: Raw signal overview ---
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(raw_signal, color="darkblue", linewidth=0.5, alpha=0.7, label="Raw Signal")
    ax1.axvspan(0, 200, alpha=0.3, color="red", label="Trimmed (adapter)")
    ax1.axvspan(
        len(raw_signal) - 10, len(raw_signal), alpha=0.3, color="orange", label="Trimmed (tail)"
    )
    ax1.set_xlabel("Sample Index", fontsize=11)
    ax1.set_ylabel("Current (pA)", fontsize=11)
    ax1.set_title("Raw Nanopore Signal", fontsize=13, fontweight="bold")
    ax1.legend(loc="upper right")
    ax1.grid(True, alpha=0.3)

    # --- Panel 2: Event-to-kmer alignment visualization ---
    ax2 = fig.add_subplot(gs[1, :])

    # Color map for different bases
    base_colors = {"A": "red", "C": "blue", "G": "green", "T": "orange", "U": "orange"}

    y_offset = 0
    for i, mapping in enumerate(base_to_event_map):
        kmer = mapping["kmer"]
        start_event = mapping["start"]
        stop_event = mapping["stop"]

        if start_event == -1 or stop_event == -1:
            continue

        # Get the middle base of the k-mer for coloring
        middle_base = kmer[kmer_size // 2]
        color = base_colors.get(middle_base, "gray")

        # Draw rectangle for this k-mer's events
        width = stop_event - start_event + 1
        rect = Rectangle(
            (start_event, y_offset),
            width,
            0.8,
            facecolor=color,
            edgecolor="black",
            alpha=0.6,
            linewidth=0.5,
        )
        ax2.add_patch(rect)

        # Add k-mer label for first few
        if i < 15 or i % 5 == 0:
            ax2.text(
                start_event + width / 2,
                y_offset + 0.4,
                kmer,
                ha="center",
                va="center",
                fontsize=6,
                fontweight="bold",
            )

        y_offset += 1

    ax2.set_xlim(0, result["n_events"])
    ax2.set_ylim(0, len(base_to_event_map))
    ax2.set_xlabel("Event Index", fontsize=11)
    ax2.set_ylabel("K-mer Position in Sequence", fontsize=11)
    ax2.set_title(f"Event-to-K-mer Alignment Map (k={kmer_size})", fontsize=13, fontweight="bold")
    ax2.grid(True, alpha=0.3, axis="x")

    # --- Panel 3: Events per k-mer distribution ---
    ax3 = fig.add_subplot(gs[2, 0])
    events_per_kmer = []
    for mapping in base_to_event_map:
        if mapping["start"] != -1 and mapping["stop"] != -1:
            events_per_kmer.append(mapping["stop"] - mapping["start"] + 1)
        else:
            events_per_kmer.append(0)

    ax3.hist(events_per_kmer, bins=30, color="purple", alpha=0.7, edgecolor="black")
    ax3.set_xlabel("Events per K-mer", fontsize=11)
    ax3.set_ylabel("Count", fontsize=11)
    ax3.set_title("Events per K-mer Distribution", fontsize=12, fontweight="bold")
    ax3.axvline(
        np.mean(events_per_kmer),
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"Mean: {np.mean(events_per_kmer):.1f}",
    )
    ax3.legend()
    ax3.grid(True, alpha=0.3)

    # --- Panel 4: Scaling parameters ---
    ax4 = fig.add_subplot(gs[2, 1])
    ax4.axis("off")

    info_text = f"""
    Eventalign Summary
    {'='*40}
    
    Sequence length:     {len(sequence)} bases
    K-mer size:          {kmer_size}
    Number of k-mers:    {len(base_to_event_map)}
    
    Events detected:     {result['n_events']}
    Aligned pairs:       {result['n_aligned_pairs']}
    Events per k-mer:    {np.mean(events_per_kmer):.2f} ± {np.std(events_per_kmer):.2f}
    
    Scaling Parameters:
    -------------------
    Scale:               {scaling['scale']:.4f}
    Shift:               {scaling['shift']:.2f} pA
    
    Coverage:
    ---------
    K-mers with events:  {sum(1 for e in events_per_kmer if e > 0)} / {len(base_to_event_map)}
    Coverage:            {100*sum(1 for e in events_per_kmer if e > 0)/len(base_to_event_map):.1f}%
    """

    ax4.text(
        0.1,
        0.9,
        info_text,
        transform=ax4.transAxes,
        fontsize=10,
        verticalalignment="top",
        family="monospace",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
    )

    # --- Panel 5: Sequence view with k-mer boundaries ---
    ax5 = fig.add_subplot(gs[3, :])

    # Display sequence with k-mer boundaries
    seq_display = list(sequence[:50])  # Show first 50 bases
    x_positions = np.arange(len(seq_display))

    for i, base in enumerate(seq_display):
        color = base_colors.get(base, "gray")
        ax5.text(i, 0, base, ha="center", va="center", fontsize=12, color=color, fontweight="bold")

        # Mark k-mer boundaries
        if i % kmer_size == 0:
            ax5.axvline(i - 0.5, color="black", linestyle="--", alpha=0.3)

    ax5.set_xlim(-1, len(seq_display))
    ax5.set_ylim(-0.5, 0.5)
    ax5.set_xlabel("Base Position", fontsize=11)
    ax5.set_title(
        f"Reference Sequence (first 50 bases, k-mer boundaries marked)",
        fontsize=12,
        fontweight="bold",
    )
    ax5.set_yticks([])
    ax5.grid(True, alpha=0.2, axis="x")

    return fig


def main():
    """Main execution function."""

    print("=" * 80)
    print("Nanopore Eventalign (Event-to-Sequence Alignment) Example".center(80))
    print("=" * 80)

    # 1. Define a reference sequence
    print("\n[1] Setting up reference sequence...")
    sequence = "AAACTGCAAATAAAATAATAAAAAGGAAAGGAAACTTTGAACCTTATGTACATGAAATAATGCAGGTCTAGCTAAACATAATGCTAGTCCTAGATTACTTATTGATTTAAAAACAAAAAACAAAAAATAGTAAAATAAAACAATTAATGTTTTCATAGACCTAAGGAAAAAGAATTTCTTGACAAGTACAAAAAAAAATTTAAAGCATTCCTTTTTAATTTGTAATTTTTTTTACTGTGGAATAAATTTGCAGAATGTCAGTTTTGTTTAGTAACAGAATTTGATATATTGAGCAGGAAACGCAATATTTGGATTATAAAATTCTTGCTTTAATAAAAATTCCTTAAACAGTGAAAAAAAT"  # 40 bases
    print(f"    Sequence: {sequence}")
    print(f"    Length: {len(sequence)} bases")

    # 2. Generate synthetic signal from sequence
    print("\n[2] Generating synthetic nanopore signal from sequence...")
    # raw_signal = generate_signal_from_sequence(sequence, samples_per_base=50, noise_level=2.5)
    with open("examples/signal.pkl", "rb") as f:
        raw_signal = pickle.load(f)
    print(f"    ✓ Generated {len(raw_signal)} samples")
    print(f"    ✓ Signal range: {raw_signal.min():.2f} - {raw_signal.max():.2f} pA")
    print(f"    ✓ Average samples per base: ~50")

    # 3. Detect events
    print("\n[3] Detecting events from raw signal...")
    try:
        events = detect_events(raw_signal, is_rna=False)
        print(f"    ✓ Detected {len(events)} events")
        event_means = [e["mean"] for e in events]
        print(
            f"    ✓ Event mean current: {np.mean(event_means):.2f} ± {np.std(event_means):.2f} pA"
        )
    except Exception as e:
        print(f"    ✗ Error: {e}")
        return

    # 4. Align events to sequence (eventalign)
    print("\n[4] Aligning events to reference sequence...")
    kmer_size = 5
    try:
        result = eventalign(raw_signal, sequence[::-1], is_rna=True, kmer_size=kmer_size)
        print(f"    ✓ Alignment complete!")
        print(f"    ✓ K-mer size: {kmer_size}")
        print(f"    ✓ Number of k-mers: {len(result['base_to_event_map'])}")
        print(f"    ✓ Aligned pairs: {result['n_aligned_pairs']}")
        print(
            f"    ✓ Scaling - scale: {result['scaling']['scale']:.4f}, "
            f"shift: {result['scaling']['shift']:.2f} pA"
        )
    except Exception as e:
        print(f"    ✗ Error: {e}")
        import traceback

        traceback.print_exc()
        return

    # 5. Display alignment details
    print("\n[5] Base-to-Event Mapping Details:")
    print("-" * 80)
    print(f"{'K-mer':<8} {'Position':<10} {'Event Start':<12} {'Event Stop':<12} {'N Events':<10}")
    print("-" * 80)

    base_to_event_map = result["base_to_event_map"]
    for i in range(min(10, len(base_to_event_map))):  # Show first 10
        mapping = base_to_event_map[i]
        kmer = mapping["kmer"]
        start = mapping["start"]
        stop = mapping["stop"]
        n_events = stop - start + 1 if start != -1 and stop != -1 else 0
        print(f"{kmer:<8} {i:<10} {start:<12} {stop:<12} {n_events:<10}")

    if len(base_to_event_map) > 10:
        print(f"... ({len(base_to_event_map) - 10} more k-mers)")

    # 6. Create visualization
    print("\n[6] Creating visualization...")
    fig = visualize_eventalign(raw_signal, sequence, result, kmer_size=kmer_size)

    # Save figure
    output_path = "eventalign_example.png"
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    print(f"    ✓ Figure saved to: {output_path}")

    plt.show()

    print("\n" + "=" * 80)
    print("Eventalign example completed successfully!".center(80))
    print("=" * 80)
    print("\nNext steps:")
    print("  - Try with different k-mer sizes (3, 5, 6)")
    print("  - Test with real nanopore data (FAST5/POD5 files)")
    print("  - Use alignment for downstream analysis (methylation calling, etc.)")


if __name__ == "__main__":
    main()
