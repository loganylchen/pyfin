#!/usr/bin/env python
"""
Visualization of eventalign results from f5c.

This script uses only the eventalign TSV file (which contains embedded raw samples)
to visualize:
- Raw signal with detected events
- Event mean levels vs model mean levels (different colors)
- Model comparison with z-scores
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
    show_n_events: int = None,
):
    """
    Create visualization of eventalign results showing all events.
    """
    # Show all events if not specified
    display_events = events if show_n_events is None else events[:show_n_events]

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
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1.2, 1], hspace=0.4)

    # =========================================================================
    # Panel 1: Raw signal with ALL events overlaid
    # =========================================================================
    ax1 = fig.add_subplot(gs[0])

    # Plot raw signal (downsample if too long for performance)
    x_signal = np.arange(len(signal))
    step = max(1, len(signal) // 50000)  # Downsample to at most 50000 points
    ax1.plot(x_signal[::step], signal[::step], "k-", alpha=0.5, linewidth=0.5, label="Raw signal")

    # Overlay ALL events with their mean levels (one color for observed, one for model)
    for i, ev in enumerate(display_events):
        start_pos, end_pos = sample_positions[i]

        if end_pos > start_pos:
            # Event mean (observed) - BLUE
            ax1.hlines(
                ev["event_mean"],
                start_pos,
                end_pos,
                colors="#1f77b4",  # Blue for event mean
                linewidths=2,
                alpha=0.8,
            )

            # Model mean (expected) - ORANGE
            ax1.hlines(
                ev["model_mean"],
                start_pos,
                end_pos,
                colors="#ff7f0e",  # Orange for model mean
                linewidths=2,
                alpha=0.6,
                linestyles="dashed",
            )

    ax1.set_ylabel("Current (pA)", fontsize=11)
    ax1.set_title(f"Raw Signal with All {len(display_events)} Events (Blue=Observed, Orange=Model)",
                 fontsize=12, fontweight="bold")
    ax1.set_xlim(0, len(signal))

    # Custom legend
    legend_elements = [
        mpatches.Patch(color="#1f77b4", label="Event Mean (Observed)"),
        mpatches.Patch(color="#ff7f0e", label="Model Mean (Expected)"),
    ]
    ax1.legend(handles=legend_elements, loc="upper right", fontsize=9)

    # =========================================================================
    # Panel 2: Side-by-side comparison of Event Mean vs Model Mean
    # =========================================================================
    ax2 = fig.add_subplot(gs[1])

    # Create bar positions
    x_centers = []
    bar_widths = []
    for i, (start_pos, end_pos) in enumerate(sample_positions):
        x_centers.append((start_pos + end_pos) / 2)
        bar_widths.append(max(1, end_pos - start_pos))

    # Plot event means in BLUE
    for i, ev in enumerate(display_events):
        ax2.bar(
            x_centers[i],
            ev["event_mean"],
            width=bar_widths[i],
            color="#1f77b4",  # Blue
            alpha=0.7,
            edgecolor="none",
            linewidth=0,
        )

    # Overlay model means in ORANGE (narrower bars on top)
    for i, ev in enumerate(display_events):
        ax2.bar(
            x_centers[i],
            ev["model_mean"],
            width=bar_widths[i] * 0.5,  # Narrower for visibility
            color="#ff7f0e",  # Orange
            alpha=0.9,
            edgecolor="white",
            linewidth=0.5,
        )

    ax2.set_ylabel("Current (pA)", fontsize=11)
    ax2.set_title("Event Mean Levels: Blue=Observed, Orange=Model (narrower bars)",
                 fontsize=12, fontweight="bold")
    ax2.set_xlim(0, len(signal))

    # =========================================================================
    # Panel 3: Event alignment by reference position (all events)
    # =========================================================================
    ax3 = fig.add_subplot(gs[2])

    # Get data for ALL events
    ref_positions = [ev["ref_pos"] for ev in display_events]
    event_means = [ev["event_mean"] for ev in display_events]
    model_means = [ev["model_mean"] for ev in display_events]

    x = np.arange(len(display_events))

    # Scatter plot for better visualization of all events
    # Event means (BLUE scatter)
    ax3.scatter(x, event_means, c="#1f77b4", alpha=0.6, s=10,
               label="Event Mean (Observed)", zorder=2)

    # Model means (ORANGE scatter)
    ax3.scatter(x, model_means, c="#ff7f0e", alpha=0.6, s=10,
               label="Model Mean (Expected)", zorder=2)

    # Connect with lines to show the relationship
    for i in range(0, len(display_events), max(1, len(display_events) // 500)):
        ax3.plot([i, i], [event_means[i], model_means[i]], "gray",
                alpha=0.2, linewidth=0.5, zorder=1)

    ax3.set_xlabel("Event Index", fontsize=11)
    ax3.set_ylabel("Current (pA)", fontsize=11)
    ax3.set_title(f"All {len(display_events)} Events: Observed vs Model (scatter plot)",
                 fontsize=12, fontweight="bold")
    ax3.legend(loc="upper right", fontsize=9)

    # =========================================================================
    # Panel 4: Z-scores showing deviation from model
    # =========================================================================
    ax4 = fig.add_subplot(gs[3])

    z_scores = [ev["standardized_mean"] for ev in display_events]

    ax4.plot(x, z_scores, "o-", color="#d62728", markersize=2, linewidth=0.8, alpha=0.7)
    ax4.axhline(0, color="black", linestyle="-", alpha=0.5, linewidth=1)
    ax4.axhline(2, color="red", linestyle="--", alpha=0.5, linewidth=1, label="+2σ")
    ax4.axhline(-2, color="red", linestyle="--", alpha=0.5, linewidth=1, label="-2σ")

    ax4.set_xlabel("Event Index", fontsize=11)
    ax4.set_ylabel("Z-Score", fontsize=11)
    ax4.set_title("Standardized Event Levels (deviation from model)",
                 fontsize=12, fontweight="bold")
    ax4.legend(loc="upper right", fontsize=9)

    # Color outliers
    outlier_mask = np.abs(z_scores) > 2
    if np.any(outlier_mask):
        outlier_x = np.array(x)[outlier_mask]
        outlier_y = np.array(z_scores)[outlier_mask]
        ax4.scatter(outlier_x, outlier_y, c="red", s=20, alpha=0.8, zorder=3)

    # =========================================================================
    # Add statistics text
    # =========================================================================
    stats_text = (
        f"Total events: {len(display_events)}\n"
        f"Signal samples: {len(signal)}\n"
        f"Mean event mean: {np.mean(event_means):.2f} pA\n"
        f"Mean model mean: {np.mean(model_means):.2f} pA\n"
        f"RMSD: {np.sqrt(np.mean((np.array(event_means) - np.array(model_means))**2)):.2f} pA\n"
        f"Outliers (|z|>2): {np.sum(outlier_mask)}/{len(z_scores)}"
    )
    fig.text(0.02, 0.5, stats_text, fontsize=9, verticalalignment="center",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.3))

    plt.tight_layout(rect=[0, 0, 0.98, 1])

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to {output_path}")

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Visualize f5c eventalign results - showing all events"
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
        default=None,
        help="Number of events to display (default: all)",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=None,
        help="Maximum events to read from file (default: all)",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("F5C EventAlign Visualization (All Events)")
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

    # Determine number of events to show
    n_show = args.n_events if args.n_events else len(events)

    # Create visualization
    print(f"\nCreating visualization (showing {n_show} events)...")
    plot_eventalign_visualization(
        events=events,
        output_path=args.output,
        show_n_events=n_show,
    )

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
