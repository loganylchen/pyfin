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
        raise ImportError("pod5 package required. Install with: pip install pod5")

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
            signal = read.signal_pa
            sample_rate = read.run_info.sample_rate
            return signal, sample_rate, str(read.read_id)


def read_eventalign_tsv(tsv_path: str, max_events: int = None) -> list[dict]:
    """Read eventalign TSV file from f5c.

    TSV format columns:
    0: ref_name, 1: ref_pos, 2: ref_kmer, 3: read_id, 4: strand,
    5: event_idx, 6: event_mean, 7: event_stdv, 8: event_length,
    9: model_kmer, 10: model_mean, 11: model_stdv,
    12: standardized_mean, 13: start_idx (signal position),
    14: end_idx (signal position)
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
                "start_idx": int(parts[13]),  # Signal start position
                "end_idx": int(parts[14]),  # Signal end position
            }
            events.append(event)

    return events


def run_getevents(signal: np.ndarray):
    """Run getevents from fin._eventalign module."""
    try:
        from fin._eventalign import getevents
    except ImportError:
        raise ImportError("fin._eventalign module not available. Please build the package first.")

    # Convert signal to float32 if needed
    if signal.dtype != np.float32:
        signal = signal.astype(np.float32)

    return getevents(signal)


def load_rna002_model():
    """Load the RNA002 pore model."""
    try:
        from fin._eventalign import set_model, MODEL_RNA002
    except ImportError:
        raise ImportError("fin._eventalign module not available. Please build the package first.")

    return set_model(MODEL_RNA002)


def kmer_to_index(kmer: str, kmer_size: int = 5) -> int:
    """
    Convert k-mer string to model index.
    RNA002 uses k=5, bases are A=0, C=1, G=2, T=3
    Index = sum(base * 4^(k-1-position)) for position in 0..k-1
    """
    base_to_val = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 0}
    index = 0
    for base in kmer:
        index = index * 4 + base_to_val[base.upper()]
    return index


def find_matching_events(f5c_events: list[dict], ge_starts: np.ndarray, ge_ends: np.ndarray):
    """
    Find matching events between f5c and getevents based on signal position overlap.

    An f5c event matches a getevents event if they overlap significantly.
    Returns a list of (f5c_idx, ge_idx) tuples for matched pairs.
    """
    matches = []
    ge_matched = set()

    for f5c_idx, f5c_ev in enumerate(f5c_events):
        f5c_start = f5c_ev["start_idx"]
        f5c_end = f5c_ev["end_idx"]

        # Find the best matching getevents event based on overlap
        best_ge_idx = None
        best_overlap = 0

        for ge_idx in range(len(ge_starts)):
            if ge_idx in ge_matched:
                continue  # Skip already matched events

            ge_start = ge_starts[ge_idx]
            ge_end = ge_ends[ge_idx]

            # Calculate overlap
            overlap_start = max(f5c_start, ge_start)
            overlap_end = min(f5c_end, ge_end)
            overlap = max(0, overlap_end - overlap_start)

            # Calculate overlap fraction (overlap relative to both events)
            f5c_len = f5c_end - f5c_start
            ge_len = ge_end - ge_start
            min_len = min(f5c_len, ge_len)
            overlap_ratio = overlap / min_len if min_len > 0 else 0

            # Match if overlap ratio is > 0.3
            if overlap_ratio > 0.3 and overlap_ratio > best_overlap:
                best_overlap = overlap_ratio
                best_ge_idx = ge_idx

        if best_ge_idx is not None:
            matches.append((f5c_idx, best_ge_idx))
            ge_matched.add(best_ge_idx)

    return matches


def plot_comparison(
    signal: np.ndarray,
    f5c_events: list[dict],
    getevents_result: dict,
    rna002_model: dict = None,
    output_path: str = None,
):
    """
    Create visualization comparing f5c events vs getevents output.
    Events are aligned based on their signal position overlap, not by event index.
    """
    # Parse getevents output
    ge_starts = getevents_result["starts"]
    ge_ends = ge_starts + getevents_result["lengths"]
    ge_means = getevents_result["means"]

    # Find matching events based on signal position overlap
    matches = find_matching_events(f5c_events, ge_starts, ge_ends)

    print(
        f"  Matched {len(matches)} pairs out of {len(f5c_events)} f5c events and {len(ge_means)} getevents"
    )

    # Create figure with 4 panels if model is provided, 3 otherwise
    n_panels = 4 if rna002_model is not None else 3
    height_ratios = [1.2, 1, 1] if rna002_model is None else [1.2, 1, 1, 1]
    fig = plt.figure(figsize=(18, 12 if rna002_model else 10))
    gs = fig.add_gridspec(n_panels, 1, height_ratios=height_ratios, hspace=0.35)

    # Downsample signal for plotting
    x_signal = np.arange(len(signal))
    step = max(1, len(signal) // 50000)

    # =========================================================================
    # Panel 1: Raw signal with both event detections
    # =========================================================================
    ax1 = fig.add_subplot(gs[0])

    # Plot raw signal
    ax1.plot(x_signal[::step], signal[::step], "k-", alpha=0.4, linewidth=0.5, label="Raw signal")

    # Overlay f5c events (BLUE) and f5c model means (ORANGE) - use actual signal positions from TSV
    for ev in f5c_events:
        start_pos = ev["start_idx"]
        end_pos = ev["end_idx"]

        if end_pos > start_pos and end_pos < len(signal):
            # F5C event mean (BLUE solid)
            ax1.hlines(
                ev["event_mean"],
                start_pos,
                end_pos,
                colors="#1f77b4",
                linewidths=2,
                alpha=0.8,
            )
            # F5C model mean (ORANGE dotted)
            ax1.hlines(
                ev["model_mean"],
                start_pos,
                end_pos,
                colors="#ff7f0e",
                linewidths=1.5,
                alpha=0.6,
                linestyles="dotted",
            )

    # Overlay getevents results (GREEN)
    for i in range(len(ge_means)):
        start_pos = int(ge_starts[i])
        end_pos = int(ge_ends[i])

        if end_pos > start_pos and end_pos < len(signal):
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
    ax1.set_title(
        "Raw Signal with Event Detection: Blue=f5c Event, Orange=f5c Model, Green=getevents (ours)",
        fontsize=12,
        fontweight="bold",
    )
    ax1.set_xlim(0, len(signal))

    legend_elements = [
        mpatches.Patch(color="#1f77b4", label="f5c Event Mean"),
        mpatches.Patch(color="#ff7f0e", label="f5c Model Mean"),
        mpatches.Patch(color="#2ca02c", label="getevents (ours)"),
    ]
    ax1.legend(handles=legend_elements, loc="upper right", fontsize=9)

    # =========================================================================
    # Panel 2: Side-by-side comparison of matched event means
    # =========================================================================
    ax2 = fig.add_subplot(gs[1])

    # Use matched pairs for comparison
    n_compare = len(matches)
    x = np.arange(n_compare)

    f5c_means = [f5c_events[f5c_idx]["event_mean"] for f5c_idx, _ in matches]
    f5c_model_means = [f5c_events[f5c_idx]["model_mean"] for f5c_idx, _ in matches]
    ge_means_matched = [ge_means[ge_idx] for _, ge_idx in matches]

    # Scatter plot comparison
    ax2.scatter(x, f5c_means, c="#1f77b4", alpha=0.6, s=15, label="f5c Event Mean", zorder=3)
    ax2.scatter(x, f5c_model_means, c="#ff7f0e", alpha=0.6, s=15, label="f5c Model Mean", zorder=2)
    ax2.scatter(x, ge_means_matched, c="#2ca02c", alpha=0.6, s=15, label="getevents Mean", zorder=3)

    # Connect with lines to show differences
    for i in range(0, n_compare, max(1, n_compare // 200)):
        ax2.plot(
            [i, i], [f5c_means[i], ge_means_matched[i]], "gray", alpha=0.2, linewidth=0.5, zorder=1
        )

    ax2.set_ylabel("Current (pA)", fontsize=11)
    ax2.set_title(
        f"Matched Event Means Comparison ({n_compare} matched pairs)",
        fontsize=12,
        fontweight="bold",
    )
    ax2.legend(loc="upper right", fontsize=9)

    # =========================================================================
    # Panel 3: Difference plot and statistics
    # =========================================================================
    ax3 = fig.add_subplot(gs[2])

    # Calculate differences for matched pairs
    f5c_arr = np.array(f5c_means)
    f5c_model_arr = np.array(f5c_model_means)
    ge_arr = np.array(ge_means_matched)
    differences = f5c_arr - ge_arr

    # Plot differences
    ax3.plot(x, differences, "o-", color="#d62728", markersize=3, linewidth=0.8, alpha=0.7, label="f5c Event - getevents")
    ax3.plot(x, f5c_model_arr - ge_arr, "o-", color="#ff7f0e", markersize=3, linewidth=0.8, alpha=0.5, label="f5c Model - getevents")
    ax3.axhline(0, color="black", linestyle="-", alpha=0.5, linewidth=1)

    # Color outliers
    outlier_mask = np.abs(differences) > 10
    if np.any(outlier_mask):
        outlier_x = np.array(x)[outlier_mask]
        outlier_y = np.array(differences)[outlier_mask]
        ax3.scatter(
            outlier_x,
            outlier_y,
            c="red",
            s=30,
            alpha=0.8,
            zorder=4,
            label=f"Large diff (>10 pA): {np.sum(outlier_mask)}",
        )

    ax3.set_xlabel("Matched Event Pair Index", fontsize=11)
    ax3.set_ylabel("Difference [pA]", fontsize=11)
    ax3.set_title(
        "Event Mean Differences (Matched by Signal Position)", fontsize=12, fontweight="bold"
    )
    ax3.legend(loc="upper right", fontsize=8)

    # =========================================================================
    # Panel 4: K-mer model comparison (f5c model_mean vs raw pore model)
    # =========================================================================
    if rna002_model is not None:
        ax4 = fig.add_subplot(gs[3])

        # Get unique kmers and their model means
        unique_kmers = {}
        for ev in f5c_events:
            kmer = ev["ref_kmer"]
            if kmer not in unique_kmers:
                unique_kmers[kmer] = ev["model_mean"]

        # Get raw model values for each kmer
        kmer_list = list(unique_kmers.keys())
        f5c_model_values = [unique_kmers[k] for k in kmer_list]
        raw_model_values = []
        raw_model_indices = []
        for kmer in kmer_list:
            idx = kmer_to_index(kmer, rna002_model["kmer_size"])
            raw_model_values.append(rna002_model["level_means"][idx])
            raw_model_indices.append(idx)

        f5c_model_arr = np.array(f5c_model_values)
        raw_model_arr = np.array(raw_model_values)

        # Scatter plot comparison
        ax4.scatter(raw_model_arr, f5c_model_arr, c="#9467bd", alpha=0.6, s=20, edgecolors="black", linewidth=0.5)

        # Add diagonal line for perfect match
        min_val = min(f5c_model_arr.min(), raw_model_arr.min())
        max_val = max(f5c_model_arr.max(), raw_model_arr.max())
        ax4.plot([min_val, max_val], [min_val, max_val], "k--", alpha=0.3, linewidth=1, label="Perfect match")

        # Calculate statistics
        model_diff = f5c_model_arr - raw_model_arr
        mean_model_diff = np.mean(model_diff)
        std_model_diff = np.std(model_diff)
        rmsd_model = np.sqrt(np.mean(model_diff**2))
        corr_model = np.corrcoef(f5c_model_arr, raw_model_arr)[0, 1] if len(f5c_model_arr) > 1 else 0.0

        ax4.set_xlabel("Raw Pore Model Mean (pA)", fontsize=11)
        ax4.set_ylabel("f5c Model Mean (pA)", fontsize=11)
        ax4.set_title(
            f"K-mer Model Comparison ({len(kmer_list)} unique kmers)\n"
            f"Mean diff: {mean_model_diff:.2f} pA, RMSD: {rmsd_model:.2f} pA, Corr: {corr_model:.4f}",
            fontsize=12,
            fontweight="bold",
        )
        ax4.legend(loc="upper left", fontsize=8)
        ax4.grid(True, alpha=0.3)

        # Update statistics text box to include k-mer model comparison
        stats_text = (
            f"Statistics:\n"
            f"f5c events: {len(f5c_events)}\n"
            f"getevents: {getevents_result['n_events']}\n"
            f"Matched pairs: {n_compare}\n\n"
            f"f5c Event vs getevents:\n"
            f"  Mean diff: {mean_diff:.2f} pA\n"
            f"  Std diff: {std_diff:.2f} pA\n"
            f"  RMSD: {rmsd:.2f} pA\n"
            f"  Correlation: {corr:.4f}\n\n"
            f"f5c Model vs getevents:\n"
            f"  Mean diff: {mean_diff_model:.2f} pA\n"
            f"  Std diff: {std_diff_model:.2f} pA\n"
            f"  RMSD: {rmsd_model:.2f} pA\n"
            f"  Correlation: {corr_model:.4f}\n\n"
            f"K-mer Model (f5c vs raw):\n"
            f"  Mean diff: {mean_model_diff:.2f} pA\n"
            f"  RMSD: {rmsd_model:.2f} pA\n"
            f"  Correlation: {corr_model:.4f}\n\n"
            f"Large diffs (>10pA): {np.sum(outlier_mask)}"
        )
    else:
        # =========================================================================
        # Add statistics box (no k-mer model comparison)
        # =========================================================================
        stats_text = (
            f"Statistics:\n"
            f"f5c events: {len(f5c_events)}\n"
            f"getevents: {getevents_result['n_events']}\n"
            f"Matched pairs: {n_compare}\n\n"
            f"f5c Event vs getevents:\n"
            f"  Mean diff: {mean_diff:.2f} pA\n"
            f"  Std diff: {std_diff:.2f} pA\n"
            f"  RMSD: {rmsd:.2f} pA\n"
            f"  Correlation: {corr:.4f}\n\n"
            f"f5c Model vs getevents:\n"
            f"  Mean diff: {mean_diff_model:.2f} pA\n"
            f"  Std diff: {std_diff_model:.2f} pA\n"
            f"  RMSD: {rmsd_model:.2f} pA\n"
            f"  Correlation: {corr_model:.4f}\n\n"
            f"Large diffs (>10pA): {np.sum(outlier_mask)}"
        )

    fig.text(
        0.02,
        0.5,
        stats_text,
        fontsize=7 if rna002_model else 8,
        verticalalignment="center",
        bbox=dict(boxstyle="round", facecolor="lightblue", alpha=0.5),
    )

    plt.tight_layout(rect=[0, 0, 0.96, 1])

    if output_path:
        plt.savefig(output_path, dpi=600, bbox_inches="tight")
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
        "--output",
        "-o",
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

    # Load RNA002 model for k-mer comparison
    print(f"\nLoading RNA002 pore model...")
    rna002_model = load_rna002_model()
    print(f"  K-mer size: {rna002_model['kmer_size']}")
    print(f"  Number of k-mers: {rna002_model['num_kmer']}")

    # Create visualization
    print(f"\nCreating comparison visualization...")
    plot_comparison(
        signal=signal,
        f5c_events=f5c_events,
        getevents_result=getevents_result,
        rna002_model=rna002_model,
        output_path=args.output,
    )

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
