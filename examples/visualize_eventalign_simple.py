#!/usr/bin/env python
"""
Simple visualization of eventalign results from f5c.

This script uses only the eventalign TSV file (which contains embedded raw samples)
to visualize:
- Raw signal with detected events
- Event mean levels colored by k-mer
- Model vs observed comparison
"""

import argparse
import gzip
import sys

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# Set style
sns.set_style("whitegrid")


def read_eventalign_tsv(tsv_path: str, max_events: int = None) -> list[dict]:
    """
    Read eventalign TSV file.

    Returns list of dicts with parsed data including raw samples.
    """
    events = []
    open_func = gzip.open if tsv_path.endswith(".gz") else open

    with open_func(tsv_path, "rt") as f:
        for i, line in enumerate(f):
            if max_events and i >= max_events:
                break
            parts = line.strip().split("\t")
            if len(parts) < 15:
                continue

            # Parse comma-separated raw samples (column 16, 0-indexed)
            raw_samples = (
                np.array([float(x) for x in parts[15].split(",")])
                if len(parts) > 15 and parts[15]
                else np.array([])
            )

            event = {
                "ref_name": parts[0],
                "ref_pos": int(parts[1]),
                "ref_kmer": parts[2],
                "read_id": parts[3],
                "strand": parts[4],
                "event_idx": int(parts[5]),
                "event_mean": float(parts[6]),
                "event_stdv": float(parts[7]),
                "event_length": float(parts[8]),
                "model_kmer": parts[9],
                "model_mean": float(parts[10]),
                "model_stdv": float(parts[11]),
                "standardized_mean": float(parts[12]) if parts[12] != "inf" else 0.0,
                "start_idx": int(parts[13]),
                "end_idx": int(parts[14]),
                "raw_samples": raw_samples,
            }
            events.append(event)

    return events


def kmer_to_color(kmer: str) -> tuple[float, float, float]:
    """Map k-mer to a unique color based on nucleotide composition."""
    base_colors = {
        "A": (0.8, 0.2, 0.2),   # Red
        "C": (0.2, 0.2, 0.8),   # Blue
        "G": (0.2, 0.8, 0.2),   # Green
        "T": (0.8, 0.8, 0.2),   # Yellow
        "U": (0.8, 0.2, 0.2),   # Red (same as A)
        "N": (0.5, 0.5, 0.5),   # Gray
    }

    # Average color for k-mer
    colors = [base_colors.get(base.upper(), (0.5, 0.5, 0.5)) for base in kmer]
    return tuple(sum(c) / len(c) for c in zip(*colors))


def plot_eventalign_visualization(
    events: list[dict],
    output_path: str = None,
    show_n_events: int = 100,
):
    """
    Create visualization of eventalign results.
    """
    # Limit events for display
    display_events = events[:show_n_events]

    # Reconstruct signal from embedded raw samples
    all_samples = []
    for ev in display_events:
        if len(ev["raw_samples"]) > 0:
            all_samples.extend(ev["raw_samples"].tolist())
    signal = np.array(all_samples)

    # Build sample index for each event
    sample_positions = []
    pos = 0
    for ev in display_events:
        sample_positions.append((pos, pos + len(ev["raw_samples"])))
        pos += len(ev["raw_samples"])

    # Create figure
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, 1.2], hspace=0.35)

    # =========================================================================
    # Panel 1: Raw signal with event boundaries
    # =========================================================================
    ax1 = fig.add_subplot(gs[0])

    # Plot raw signal
    x_signal = np.arange(len(signal))
    ax1.plot(x_signal, signal, "k-", alpha=0.4, linewidth=0.8, label="Raw signal")

    # Overlay events with their mean levels and colors
    for i, ev in enumerate(display_events):
        color = kmer_to_color(ev["ref_kmer"])
        start_pos, end_pos = sample_positions[i]

        # Draw event mean as horizontal line
        if end_pos > start_pos:
            ax1.hlines(
                ev["event_mean"],
                start_pos,
                end_pos,
                colors=color,
                linewidths=2.5,
                alpha=0.9,
            )

        # Draw event boundary lines
        ax1.axvline(start_pos, color=color, linewidth=0.5, alpha=0.3)
        ax1.axvline(end_pos, color=color, linewidth=0.5, alpha=0.3)

    ax1.set_ylabel("Current (pA)", fontsize=11)
    ax1.set_title("Raw Signal with Detected Events (colored by k-mer)", fontsize=12, fontweight="bold")
    ax1.set_xlim(0, len(signal))
    ax1.legend(loc="upper right")

    # =========================================================================
    # Panel 2: Event mean levels colored by k-mer
    # =========================================================================
    ax2 = fig.add_subplot(gs[1])

    # Plot event means as colored bars
    x_centers = []
    bar_widths = []
    for i, (start_pos, end_pos) in enumerate(sample_positions):
        x_centers.append((start_pos + end_pos) / 2)
        bar_widths.append(end_pos - start_pos)

    for i, ev in enumerate(display_events):
        color = kmer_to_color(ev["ref_kmer"])
        ax2.bar(
            x_centers[i],
            ev["event_mean"],
            width=bar_widths[i],
            color=color,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.5,
        )

    ax2.set_ylabel("Event Mean (pA)", fontsize=11)
    ax2.set_title("Event Mean Levels by K-mer", fontsize=12, fontweight="bold")
    ax2.set_xlim(0, len(signal))

    # =========================================================================
    # Panel 3: Event alignment with model comparison
    # =========================================================================
    ax3 = fig.add_subplot(gs[2])

    # Get data
    ref_positions = [ev["ref_pos"] for ev in display_events]
    event_means = [ev["event_mean"] for ev in display_events]
    model_means = [ev["model_mean"] for ev in display_events]
    z_scores = [ev["standardized_mean"] for ev in display_events]

    x = np.arange(len(display_events))
    width = 0.35

    # Color by k-mer
    colors = [kmer_to_color(ev["ref_kmer"]) for ev in display_events]

    # Create twin axis for z-scores
    ax3_twin = ax3.twinx()

    # Plot event means (bars)
    bars1 = ax3.bar(x - width/2, event_means, width, label="Observed Event Mean",
                   color=colors, alpha=0.8, edgecolor="black", linewidth=0.5)

    # Plot model means (gray bars for comparison)
    bars2 = ax3.bar(x + width/2, model_means, width, label="Model Expected Mean",
                   color="lightgray", alpha=0.7, edgecolor="black", linewidth=0.5)

    # Plot z-scores as line on twin axis
    line = ax3_twin.plot(x, z_scores, "o-", color="darkred", markersize=4,
                        linewidth=1.5, label="Z-Score", alpha=0.8)
    ax3_twin.axhline(0, color="black", linestyle="--", alpha=0.3)
    ax3_twin.axhline(2, color="red", linestyle="--", alpha=0.3)
    ax3_twin.axhline(-2, color="red", linestyle="--", alpha=0.3)

    # Add reference k-mers as text labels (sparse)
    for i, ev in enumerate(display_events):
        if i % 10 == 0:  # Label every 10th event
            ax3.text(i, ax3.get_ylim()[0] + 5, ev["ref_kmer"],
                    ha="center", va="bottom", fontsize=7, rotation=45,
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))

    ax3.set_xlabel("Event Index", fontsize=11)
    ax3.set_ylabel("Current (pA)", fontsize=11)
    ax3_twin.set_ylabel("Z-Score", fontsize=11, color="darkred")
    ax3_twin.tick_params(axis="y", labelcolor="darkred")
    ax3.set_title("Event Alignment: Observed vs Model (with Z-Scores)", fontsize=12, fontweight="bold")

    # Combine legends
    bars = [bars1, bars2] + line
    labels = [b.get_label() for b in bars]
    ax3.legend(bars, labels, loc="upper right")

    # =========================================================================
    # Add k-mer color legend
    # =========================================================================
    unique_kmers = sorted(set(ev["ref_kmer"] for ev in display_events))
    # Limit to top 20 most common for readability
    from collections import Counter
    kmer_counts = Counter(ev["ref_kmer"] for ev in display_events)
    top_kmers = [k for k, _ in kmer_counts.most_common(20)]

    legend_elements = [mpatches.Patch(color=kmer_to_color(k), label=f"{k} ({kmer_counts[k]})")
                      for k in top_kmers]
    fig.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(0.92, 0.5),
              title="K-mer (count)", fontsize=9)

    plt.tight_layout(rect=[0, 0, 0.88, 1])

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to {output_path}")

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Visualize f5c eventalign results"
    )
    parser.add_argument(
        "eventalign",
        nargs="?",
        type=str,
        default="examples/test_data/one_read.eventalign.tsv.gz",
        help="Path to eventalign TSV file (can be gzipped)",
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None,
        help="Path to save output figure",
    )
    parser.add_argument(
        "--n-events", "-n",
        type=int,
        default=100,
        help="Number of events to display (default: 100)",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Maximum events to read from file (default: all)",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("F5C EventAlign Visualization")
    print("=" * 60)

    # Read eventalign results
    print(f"\nReading: {args.eventalign}")
    events = read_eventalign_tsv(args.eventalign, max_events=args.max_events)
    print(f"  Total events: {len(events)}")

    if len(events) == 0:
        print("Error: No events found in eventalign file")
        sys.exit(1)

    # Get read info
    read_id = events[0]["read_id"]
    ref_name = events[0]["ref_name"]
    print(f"  Read ID: {read_id}")
    print(f"  Reference: {ref_name}")
    print(f"  Reference positions: {events[0]['ref_pos']} - {events[-1]['ref_pos']}")

    # Count events with raw samples
    n_with_samples = sum(1 for ev in events if len(ev["raw_samples"]) > 0)
    print(f"  Events with raw samples: {n_with_samples}/{len(events)}")

    # Create visualization
    print(f"\nCreating visualization (showing first {args.n_events} events)...")
    plot_eventalign_visualization(
        events=events,
        output_path=args.output,
        show_n_events=args.n_events,
    )

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
