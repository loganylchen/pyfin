#!/usr/bin/env python
"""
Visualization comparing f5c eventalign results vs getevents output.

This script:
1. Reads raw signal from POD5 file
2. Runs getevents to detect events
3. Reads f5c eventalign results from TSV
4. Compares the two event detection methods
"""

import argparse
import gzip
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# Set style
sns.set_style("whitegrid")


def read_pod5_signal(pod5_path: str, read_id: str = None) -> tuple[np.ndarray, float, str]:
    """Read raw signal from POD5 file."""
    try:
        import pod5
    except ImportError:
        raise ImportError(
            "pod5 package required. Install with: pip install pod5"
        )

    with pod5.Reader(pod5_path) as reader:
        reads = list(reader.reads())

        if read_id:
            # Find specific read
            for read in reads:
                if str(read.read_id) == read_id:
                    signal = read.signal
                    sample_rate = read.run_info.sample_rate
                    return signal, sample_rate, str(read.read_id)
            raise ValueError(f"Read ID {read_id} not found in {pod5_path}")
        else:
            # Return first read
            if len(reads) == 0:
                raise ValueError(f"No reads found in {pod5_path}")
            read = reads[0]
            signal = read.signal
            sample_rate = read.run_info.sample_rate
            return signal, sample_rate, str(read.read_id)


def read_eventalign_tsv(tsv_path: str, max_events: int = None) -> list[dict]:
    """Read eventalign TSV file from f5c."""
    events = []
    open_func = gzip.open if tsv_path.endswith(".gz") else open

    with open_func(tsv_path, "rt") as f:
        for i, line in enumerate(f):
            if max_events and i >= max_events:
                break
            parts = line.strip().split("\t")
            if len(parts) < 15:
                continue

            event = {
                "ref_name": parts[0],
                "ref_pos": int(parts[1]),
                "ref_kmer": parts[2],
                "read_id": parts[3],
                "event_idx": int(parts[5]),
                "event_mean": float(parts[6]),
                "event_stdv": float(parts[7]),
                "model_mean": float(parts[10]),
                "model_stdv": float(parts[11]),
                "standardized_mean": float(parts[12]) if parts[12] != "inf" else 0.0,
            }
            events.append(event)

    return events


def run_getevents(signal: np.ndarray):
    """Run getevents from fin._eventalign module."""
    try:
        from fin._eventalign import getevents
    except ImportError:
        raise ImportError(
            "fin._eventalign module not available. Please build the package first."
        )

    # Convert signal to float32 if needed
    if signal.dtype != np.float32:
        signal = signal.astype(np.float32)

    return getevents(signal)


def plot_comparison(
    signal: np.ndarray,
    f5c_events: list[dict],
    getevents_result: dict,
    output_path: str = None,
):
    """
    Create visualization comparing f5c events vs getevents output.
    """
    # Parse getevents output
    ge_starts = getevents_result["starts"]
    ge_ends = ge_starts + getevents_result["lengths"]
    ge_means = getevents_result["means"]

    # Create figure
    fig = plt.figure(figsize=(18, 10))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.2, 1, 1], hspace=0.35)

    # Downsample signal for plotting
    x_signal = np.arange(len(signal))
    step = max(1, len(signal) // 50000)

    # =========================================================================
    # Panel 1: Raw signal with both event detections
    # =========================================================================
    ax1 = fig.add_subplot(gs[0])

    # Plot raw signal
    ax1.plot(x_signal[::step], signal[::step], "k-", alpha=0.4, linewidth=0.5, label="Raw signal")

    # Overlay f5c events (BLUE)
    for ev in f5c_events:
        # Use ref_pos as a proxy for position (not exactly sample position but ordered)
        idx = ev["event_idx"]
        if idx < len(ge_starts):
            start_pos = int(ge_starts[idx])
            end_pos = int(ge_ends[idx])
        else:
            continue

        if end_pos > start_pos:
            # F5C event mean (BLUE solid)
            ax1.hlines(
                ev["event_mean"],
                start_pos,
                end_pos,
                colors="#1f77b4",
                linewidths=2,
                alpha=0.8,
            )

    # Overlay getevents results (GREEN)
    for i in range(len(ge_means)):
        start_pos = int(ge_starts[i])
        end_pos = int(ge_ends[i])

        if end_pos > start_pos:
            # getevents mean (GREEN dashed)
            ax1.hlines(
                ge_means[i],
                start_pos,
                end_pos,
                colors="#2ca02c",
                linewidths=2,
                alpha=0.5,
                linestyles="dashed",
            )

    ax1.set_ylabel("Current (pA)", fontsize=11)
    ax1.set_title("Raw Signal with Event Detection: Blue=f5c, Green=getevents (ours)",
                 fontsize=12, fontweight="bold")
    ax1.set_xlim(0, len(signal))

    legend_elements = [
        mpatches.Patch(color="#1f77b4", label="f5c Event Mean"),
        mpatches.Patch(color="#2ca02c", label="getevents (ours)"),
    ]
    ax1.legend(handles=legend_elements, loc="upper right", fontsize=9)

    # =========================================================================
    # Panel 2: Side-by-side comparison of event means
    # =========================================================================
    ax2 = fig.add_subplot(gs[1])

    # Align by event index
    n_compare = min(len(f5c_events), len(ge_means))
    x = np.arange(n_compare)

    f5c_means = [ev["event_mean"] for ev in f5c_events[:n_compare]]
    ge_means_arr = ge_means[:n_compare]

    # Scatter plot comparison
    ax2.scatter(x, f5c_means, c="#1f77b4", alpha=0.6, s=15,
               label="f5c Event Mean", zorder=2)
    ax2.scatter(x, ge_means_arr, c="#2ca02c", alpha=0.6, s=15,
               label="getevents Mean", zorder=2)

    # Connect with lines to show differences
    for i in range(0, n_compare, max(1, n_compare // 200)):
        ax2.plot([i, i], [f5c_means[i], ge_means_arr[i]], "gray",
                alpha=0.2, linewidth=0.5, zorder=1)

    ax2.set_ylabel("Current (pA)", fontsize=11)
    ax2.set_title(f"Event Means Comparison ({n_compare} events)", fontsize=12, fontweight="bold")
    ax2.legend(loc="upper right", fontsize=9)

    # =========================================================================
    # Panel 3: Difference plot and statistics
    # =========================================================================
    ax3 = fig.add_subplot(gs[2])

    # Calculate differences
    f5c_arr = np.array(f5c_means)
    ge_arr = np.array(ge_means_arr)
    differences = f5c_arr - ge_arr
    relative_diff = differences / f5c_arr * 100

    # Plot differences
    ax3.plot(x, differences, "o-", color="#d62728", markersize=3, linewidth=0.8, alpha=0.7)
    ax3.axhline(0, color="black", linestyle="-", alpha=0.5, linewidth=1)

    # Color outliers
    outlier_mask = np.abs(differences) > 10
    if np.any(outlier_mask):
        outlier_x = np.array(x)[outlier_mask]
        outlier_y = np.array(differences)[outlier_mask]
        ax3.scatter(outlier_x, outlier_y, c="red", s=30, alpha=0.8, zorder=3,
                   label=f"Large diff (>10 pA): {np.sum(outlier_mask)}")

    ax3.set_xlabel("Event Index", fontsize=11)
    ax3.set_ylabel("Difference (f5c - getevents) [pA]", fontsize=11)
    ax3.set_title("Event Mean Differences", fontsize=12, fontweight="bold")
    ax3.legend(loc="upper right", fontsize=9)

    # =========================================================================
    # Add statistics box
    # =========================================================================
    mean_diff = np.mean(differences)
    std_diff = np.std(differences)
    rmsd = np.sqrt(np.mean(differences**2))
    corr = np.corrcoef(f5c_arr, ge_arr)[0, 1]

    stats_text = (
        f"Statistics:\n"
        f"f5c events: {len(f5c_events)}\n"
        f"getevents: {getevents_result['n_events']}\n"
        f"Compared: {n_compare}\n\n"
        f"Mean diff: {mean_diff:.2f} pA\n"
        f"Std diff: {std_diff:.2f} pA\n"
        f"RMSD: {rmsd:.2f} pA\n"
        f"Correlation: {corr:.4f}\n"
        f"Large diffs (>10pA): {np.sum(outlier_mask)}"
    )
    fig.text(0.02, 0.5, stats_text, fontsize=9, verticalalignment="center",
            bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.5))

    plt.tight_layout(rect=[0, 0, 0.96, 1])

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to {output_path}")

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Compare f5c eventalign events vs getevents output"
    )
    parser.add_argument(
        "--pod5",
        type=str,
        default="examples/test_data/one_read.pod5",
        help="Path to POD5 file",
    )
    parser.add_argument(
        "--eventalign",
        type=str,
        default="examples/test_data/one_read.eventalign.tsv.gz",
        help="Path to f5c eventalign TSV file",
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None,
        help="Path to save output figure",
    )

    args = parser.parse_args()

    print("=" * 60)
    print("F5C vs getevents Comparison")
    print("=" * 60)

    # Read POD5 signal
    print(f"\nReading POD5: {args.pod5}")
    signal, sample_rate, read_id = read_pod5_signal(args.pod5)
    print(f"  Read ID: {read_id}")
    print(f"  Signal length: {len(signal)} samples")
    print(f"  Sample rate: {sample_rate} Hz")

    # Run getevents
    print(f"\nRunning getevents...")
    getevents_result = run_getevents(signal)
    print(f"  Detected {getevents_result['n_events']} events")

    # Read f5c eventalign results
    print(f"\nReading f5c eventalign: {args.eventalign}")
    f5c_events = read_eventalign_tsv(args.eventalign)
    print(f"  Total events: {len(f5c_events)}")

    # Create visualization
    print(f"\nCreating comparison visualization...")
    plot_comparison(
        signal=signal,
        f5c_events=f5c_events,
        getevents_result=getevents_result,
        output_path=args.output,
    )

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
