#!/usr/bin/env python
"""
Visualize eventalign results from f5c with raw signal and events.

This script loads:
1. Raw signal from POD5 file
2. Reference sequence from FASTA
3. Event alignment results from eventalign TSV

Then creates an interactive visualization showing:
- Raw signal with detected events overlaid
- Event-to-reference alignment
- K-mer level expectations vs observed event means
"""

import argparse
import gzip
import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import seaborn as sns

# Set style
sns.set_style("whitegrid")


def read_fasta(fasta_path: str) -> tuple[str, str]:
    """Read reference sequence from FASTA file."""
    with open(fasta_path, "r") as f:
        header = f.readline().strip().lstrip(">")
        seq = "".join(line.strip() for line in f)
    return header, seq


def read_eventalign_tsv(tsv_path: str) -> list[dict]:
    """
    Read eventalign TSV file.

    Columns (nanopolish/f5c format):
    1. ref_name
    2. ref_pos (0-based)
    3. ref_kmer
    4. read_id
    5. strand
    6. event_idx
    7. event_mean
    8. event_stdv
    9. event_length
    10. model_kmer
    11. model_mean
    12. model_stdv
    13. standardized_mean
    14. start_idx (sample position)
    15. end_idx (sample position)
    16. raw_samples (comma-separated)

    Returns list of dicts with parsed data.
    """
    events = []
    open_func = gzip.open if tsv_path.endswith(".gz") else open

    with open_func(tsv_path, "rt") as f:
        for line in f:
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


def read_pod5_signal(pod5_path: str, read_id: str) -> tuple[np.ndarray, float]:
    """Read raw signal from POD5 file."""
    try:
        import pod5
    except ImportError:
        raise ImportError(
            "pod5 package required. Install with: pip install pod5"
        )

    with pod5.Reader(pod5_path) as reader:
        for read in reader:
            if str(read.read_id) == read_id:
                signal = read.signal
                sample_rate = read.run_info.sample_rate
                return signal, sample_rate

    raise ValueError(f"Read ID {read_id} not found in {pod5_path}")


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
    signal: np.ndarray,
    events: list[dict],
    ref_seq: str,
    sample_rate: float = 4000.0,
    output_path: str = None,
    show_first_n_events: int = 100,
):
    """
    Create visualization of eventalign results.

    Parameters
    ----------
    signal : np.ndarray
        Raw signal data
    events : list of dict
        Parsed eventalign events
    ref_seq : str
        Reference sequence
    sample_rate : float
        Signal sample rate in Hz
    output_path : str, optional
        Path to save figure
    show_first_n_events : int
        Number of events to show in detail view
    """
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1.5, 1], hspace=0.3)

    # Limit events for display
    display_events = events[:show_first_n_events]

    # Calculate time axis
    time_samples = np.arange(len(signal)) / sample_rate

    # =========================================================================
    # Panel 1: Raw signal with event boundaries
    # =========================================================================
    ax1 = fig.add_subplot(gs[0])

    # Plot raw signal (downsample if too long)
    step = max(1, len(signal) // 20000)
    ax1.plot(time_samples[::step], signal[::step], "k-", alpha=0.3, linewidth=0.5)

    # Overlay events with their mean levels
    for i, ev in enumerate(display_events):
        color = kmer_to_color(ev["ref_kmer"])
        start_time = ev["start_idx"] / sample_rate
        end_time = ev["end_idx"] / sample_rate

        # Draw event mean as horizontal line
        ax1.hlines(
            ev["event_mean"],
            start_time,
            end_time,
            colors=color,
            linewidths=2,
            alpha=0.8,
        )

        # Draw event boundaries
        ax1.axvline(start_time, color=color, linewidth=0.5, alpha=0.3)
        ax1.axvline(end_time, color=color, linewidth=0.5, alpha=0.3)

    ax1.set_ylabel("Current (pA)")
    ax1.set_title("Raw Signal with Detected Events")
    ax1.set_xlim(0, max(len(signal) / sample_rate, display_events[-1]["end_idx"] / sample_rate))

    # =========================================================================
    # Panel 2: Event mean levels colored by k-mer
    # =========================================================================
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    # Plot event means as colored bars
    for ev in display_events:
        color = kmer_to_color(ev["ref_kmer"])
        start_time = ev["start_idx"] / sample_rate
        end_time = ev["end_idx"] / sample_rate
        mid_time = (start_time + end_time) / 2

        ax2.bar(
            mid_time,
            ev["event_mean"],
            width=end_time - start_time,
            color=color,
            alpha=0.8,
            edgecolor="none",
        )

    ax2.set_ylabel("Event Mean (pA)")
    ax2.set_title("Event Mean Levels by K-mer")

    # =========================================================================
    # Panel 3: Event alignment to reference with model comparison
    # =========================================================================
    ax3 = fig.add_subplot(gs[2])

    # Create alignment visualization
    ref_positions = [ev["ref_pos"] for ev in display_events]
    event_means = [ev["event_mean"] for ev in display_events]
    model_means = [ev["model_mean"] for ev in display_events]

    # Plot event means and model means
    x = np.arange(len(display_events))
    width = 0.35

    # Color by k-mer
    colors = [kmer_to_color(ev["ref_kmer"]) for ev in display_events]

    bars1 = ax3.bar(x - width/2, event_means, width, label="Event Mean", color=colors, alpha=0.8)
    bars2 = ax3.bar(x + width/2, model_means, width, label="Model Mean", color="gray", alpha=0.5)

    # Add reference sequence as text labels
    for i, ev in enumerate(display_events):
        if i % 5 == 0:  # Don't crowd the x-axis
            ax3.text(i, ax3.get_ylim()[0] + 2, ev["ref_kmer"],
                    ha="center", va="bottom", fontsize=8, rotation=45)

    ax3.set_xlabel("Event Index")
    ax3.set_ylabel("Current (pA)")
    ax3.set_title("Event Means vs Model Expectations")
    ax3.legend(loc="upper right")

    # =========================================================================
    # Panel 4: Standardized event levels (z-scores)
    # =========================================================================
    ax4 = fig.add_subplot(gs[3], sharex=ax3)

    z_scores = [ev["standardized_mean"] for ev in display_events]
    ax4.plot(x, z_scores, "o-", markersize=3, linewidth=1)
    ax4.axhline(0, color="black", linestyle="--", alpha=0.5)
    ax4.axhline(2, color="red", linestyle="--", alpha=0.3, label="+2σ")
    ax4.axhline(-2, color="red", linestyle="--", alpha=0.3, label="-2σ")

    ax4.set_xlabel("Event Index")
    ax4.set_ylabel("Z-Score")
    ax4.set_title("Standardized Event Levels (Model Deviations)")
    ax4.legend(loc="upper right")

    # =========================================================================
    # Add k-mer color legend
    # =========================================================================
    # Create unique k-mers for legend
    unique_kmers = sorted(set(ev["ref_kmer"] for ev in display_events))[:20]
    legend_elements = [mpatches.Patch(color=kmer_to_color(k), label=k) for k in unique_kmers]
    fig.legend(handles=legend_elements, loc="center left", bbox_to_anchor=(0.92, 0.5), title="K-mer")

    plt.tight_layout(rect=[0, 0, 0.92, 1])

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to {output_path}")

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Visualize f5c eventalign results with raw signal"
    )
    parser.add_argument(
        "--pod5",
        type=str,
        default="examples/test_data/one_read.pod5",
        help="Path to POD5 file containing raw signal",
    )
    parser.add_argument(
        "--fasta",
        type=str,
        default="examples/test_data/one_read.fa",
        help="Path to reference FASTA file",
    )
    parser.add_argument(
        "--eventalign",
        type=str,
        default="examples/test_data/one_read.eventalign.tsv.gz",
        help="Path to eventalign TSV file (can be gzipped)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to save output figure",
    )
    parser.add_argument(
        "--n-events",
        type=int,
        default=100,
        help="Number of events to display in detail (default: 100)",
    )
    parser.add_argument(
        "--sample-rate",
        type=float,
        default=4000.0,
        help="Signal sample rate in Hz (default: 4000)",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("F5C EventAlign Visualization")
    print("=" * 60)

    # Read reference sequence
    print(f"\nReading reference from: {args.fasta}")
    ref_header, ref_seq = read_fasta(args.fasta)
    print(f"  Reference: {ref_header}")
    print(f"  Length: {len(ref_seq)} bp")

    # Read eventalign results
    print(f"\nReading eventalign results from: {args.eventalign}")
    events = read_eventalign_tsv(args.eventalign)
    print(f"  Total events: {len(events)}")

    if len(events) == 0:
        print("Error: No events found in eventalign file")
        return

    # Get read ID from events
    read_id = events[0]["read_id"]
    print(f"  Read ID: {read_id}")
    print(f"  Reference positions: {events[0]['ref_pos']} - {events[-1]['ref_pos']}")

    # Read raw signal
    print(f"\nReading raw signal from: {args.pod5}")
    try:
        signal, sample_rate = read_pod5_signal(args.pod5, read_id)
        print(f"  Signal length: {len(signal)} samples")
        print(f"  Sample rate: {sample_rate} Hz")
        print(f"  Duration: {len(signal) / sample_rate:.2f} seconds")
    except Exception as e:
        print(f"  Warning: Could not read POD5 file: {e}")
        print("  Using signal samples from eventalign data...")
        # Concatenate all raw samples from events
        all_samples = []
        for ev in events:
            if len(ev["raw_samples"]) > 0:
                all_samples.extend(ev["raw_samples"].tolist())
        signal = np.array(all_samples)
        sample_rate = args.sample_rate
        print(f"  Reconstructed signal length: {len(signal)} samples")

    # Create visualization
    print(f"\nCreating visualization...")
    plot_eventalign_visualization(
        signal=signal,
        events=events,
        ref_seq=ref_seq,
        sample_rate=sample_rate,
        output_path=args.output,
        show_first_n_events=args.n_events,
    )

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
