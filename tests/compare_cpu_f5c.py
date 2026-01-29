#!/usr/bin/env python3
"""
Compare CPU eventalign result with f5c eventalign result.

This script:
1. Loads f5c eventalign results from a file
2. Runs CPU eventalign on the same data
3. Plots signal with events from both tools, showing matching alignments
4. Compares alignment statistics and identifies differences
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from fin._eventalign import run_eventalign, MODEL_RNA002

# For loading signal from POD5 files
try:
    from fin.io.io_pod5 import Pod5Reader
    POD5_AVAILABLE = True
except ImportError:
    POD5_AVAILABLE = False
    print("Warning: pod5 not available, will use synthetic signal")


def load_reference_from_fasta(fasta_path: str) -> tuple:
    """Load reference sequence from FASTA file."""
    with open(fasta_path, "r") as f:
        lines = f.readlines()

    ref_name = None
    ref_sequence = []
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            ref_name = line[1:].split()[0]
        else:
            ref_sequence.append(line)

    ref_sequence = "".join(ref_sequence).upper()
    ref_length = len(ref_sequence)

    return ref_name, ref_sequence, ref_length


def load_signal_from_pod5(pod5_path: str, read_id: str = None) -> tuple:
    """Load signal data from POD5 file."""
    if not POD5_AVAILABLE:
        # Generate synthetic signal
        print("Generating synthetic signal...")
        read_id = "synthetic_read"
        np.random.seed(42)
        n_samples = 100000
        signal = np.random.randn(n_samples).astype(np.float32) * 10 + 120
        sample_rate = 4000.0
        return read_id, signal, sample_rate

    with Pod5Reader(pod5_path) as reader:
        if read_id is None:
            read_id = reader.read_ids[0]

        read = reader.get_read(read_id)
        if read is None:
            raise ValueError(f"Read {read_id} not found in POD5 file")

        signal = read.signal_pa.astype(np.float32)
        sample_rate = float(read.run_info.sample_rate)

        return str(read.read_id), signal, sample_rate


def load_f5c_eventalign_results(eventalign_path: str) -> dict:
    """
    Load f5c eventalign results from file.

    The f5c eventalign output format (tsv):
    chromosome    position    reference_kmer    read_index    strand    event_idx    event_mean    event_stdv    model_mean    model_stdv    scaled_mean    scaled_stdv    start_idx    end_idx
    """
    alignments = []

    with open(eventalign_path, "r") as f:
        for line in f:
            if line.startswith("chromosome"):
                continue  # Skip header

            parts = line.strip().split("\t")
            if len(parts) < 14:
                continue

            alignment = {
                "ref_position": int(parts[1]),
                "ref_kmer": parts[2],
                "event_idx": int(parts[4]),
                "event_mean": float(parts[5]),
                "event_stdv": float(parts[6]),
                "model_mean": float(parts[7]),
                "model_stdv": float(parts[8]),
                "start": int(parts[12]),
                "end": int(parts[13]),
            }
            alignments.append(alignment)

    return alignments


def f5c_alignments_to_dict(f5c_alignments: list) -> list:
    """Convert f5c alignments to a format matching CPU eventalign output."""
    result = []
    for align in f5c_alignments:
        result.append({
            "ref_position": align["ref_position"],
            "ref_kmer": align["ref_kmer"],
            "event_idx": align["event_idx"],
            "rc": False,  # f5c may have strand info
            "model_kmer": align["ref_kmer"],  # f5c uses ref_kmer as model_kmer
            "hmm_state": "C",  # f5c doesn't output HMM state
            "event_mean": align["event_mean"],
            "event_stdv": align["event_stdv"],
        })
    return result


def plot_comparison(
    signal: np.ndarray,
    cpu_result: dict,
    f5c_alignments: list,
    output_path: str = None
):
    """
    Plot signal with events from CPU and f5c eventalign.

    Shows:
    - Raw signal trace
    - CPU eventalign results as points
    - f5c eventalign results as points
    - Matching events highlighted
    """
    events = cpu_result['events'][0]
    cpu_alignments = cpu_result['full'][0][0]

    # Convert f5c alignments to dict format
    f5c_aligns = f5c_alignments_to_dict(f5c_alignments)

    # Create position mappings for comparison
    cpu_positions = {(a['ref_position'], a['event_idx']): a for a in cpu_alignments}
    f5c_positions = {(a['ref_position'], a['event_idx']): a for a in f5c_aligns}

    # Find matching and different alignments
    matching = []
    cpu_only = []
    f5c_only = []

    for key in cpu_positions:
        if key in f5c_positions:
            matching.append((cpu_positions[key], f5c_positions[key]))
        else:
            cpu_only.append(cpu_positions[key])

    for key in f5c_positions:
        if key not in cpu_positions:
            f5c_only.append(f5c_positions[key])

    # Create figure with subplots
    fig, axes = plt.subplots(3, 1, figsize=(16, 12), sharex=True)

    # Determine signal range to plot
    max_samples = min(50000, len(signal))

    # Plot 1: CPU eventalign
    ax = axes[0]
    sample_indices = np.arange(max_samples)
    ax.plot(sample_indices, signal[:max_samples], color='lightgray', alpha=0.5, linewidth=0.5, label='Raw Signal')

    # Plot CPU alignments as scatter points
    for align in cpu_alignments:
        event_idx = align['event_idx']
        if event_idx < len(events['starts']):
            start = events['starts'][event_idx]
            if start < max_samples:
                ax.scatter(start, events['means'][event_idx], c='blue', s=20, alpha=0.7, edgecolors='none')

    ax.set_ylabel('Current (pA)', fontsize=11, fontweight='bold')
    ax.set_title(f'CPU Eventalign ({len(cpu_alignments)} alignments)', fontsize=12, fontweight='bold')
    ax.legend(['Signal', 'Events'], loc='upper right')
    ax.grid(True, alpha=0.3)

    # Plot 2: f5c eventalign
    ax = axes[1]
    ax.plot(sample_indices, signal[:max_samples], color='lightgray', alpha=0.5, linewidth=0.5, label='Raw Signal')

    # Plot f5c alignments as scatter points
    for align in f5c_aligns:
        event_idx = align['event_idx']
        if event_idx < len(events['starts']):
            start = events['starts'][event_idx]
            if start < max_samples:
                ax.scatter(start, align['event_mean'], c='green', s=20, alpha=0.7, edgecolors='none')

    ax.set_ylabel('Current (pA)', fontsize=11, fontweight='bold')
    ax.set_title(f'F5C Eventalign ({len(f5c_aligns)} alignments)', fontsize=12, fontweight='bold')
    ax.legend(['Signal', 'Events'], loc='upper right')
    ax.grid(True, alpha=0.3)

    # Plot 3: Comparison overlay
    ax = axes[2]
    ax.plot(sample_indices, signal[:max_samples], color='lightgray', alpha=0.5, linewidth=0.5, label='Raw Signal')

    # Plot matching events
    for cpu_a, f5c_a in matching:
        event_idx = cpu_a['event_idx']
        if event_idx < len(events['starts']):
            start = events['starts'][event_idx]
            if start < max_samples:
                ax.scatter(start, events['means'][event_idx], c='lime', s=30, alpha=0.8,
                          marker='o', edgecolors='black', linewidth=0.5, label='Matching' if matching.index((cpu_a, f5c_a)) == 0 else "")

    # Plot CPU-only events
    for align in cpu_only[:50]:  # Limit to avoid overcrowding
        event_idx = align['event_idx']
        if event_idx < len(events['starts']):
            start = events['starts'][event_idx]
            if start < max_samples:
                ax.scatter(start, events['means'][event_idx], c='blue', s=20, alpha=0.5,
                          marker='x', label='CPU only' if cpu_only.index(align) == 0 else "")

    # Plot f5c-only events
    for align in f5c_only[:50]:  # Limit to avoid overcrowding
        event_idx = align['event_idx']
        if event_idx < len(events['starts']):
            start = events['starts'][event_idx]
            if start < max_samples:
                ax.scatter(start, align['event_mean'], c='red', s=20, alpha=0.5,
                          marker='+', label='F5C only' if f5c_only.index(align) == 0 else "")

    ax.set_ylabel('Current (pA)', fontsize=11, fontweight='bold')
    ax.set_xlabel('Sample Index', fontsize=11, fontweight='bold')
    ax.set_title(f'Comparison: {len(matching)} matching, {len(cpu_only)} CPU only, {len(f5c_only)} F5C only',
                 fontsize=12, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.suptitle('CPU vs F5C Eventalign Comparison', fontsize=14, fontweight='bold')
    plt.tight_layout()

    if output_path is None:
        output_path = Path(__file__).parent / "cpu_f5c_comparison.png"

    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Comparison plot saved to: {output_path}")

    return len(matching), len(cpu_only), len(f5c_only)


def plot_alignment_trace(
    signal: np.ndarray,
    cpu_result: dict,
    f5c_alignments: list,
    output_path: str = None
):
    """
    Plot alignment trace showing reference positions on the signal.

    Creates a visualization similar to f5c's eventalign output,
    showing which reference k-mers align to which events.
    """
    events = cpu_result['events'][0]
    cpu_alignments = cpu_result['full'][0][0]
    f5c_aligns = f5c_alignments_to_dict(f5c_alignments)

    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)

    max_samples = min(50000, len(signal))
    sample_indices = np.arange(max_samples)

    # Plot 1: CPU alignments
    ax = axes[0]
    ax.plot(sample_indices, signal[:max_samples], color='lightgray', alpha=0.3, linewidth=0.5)

    # Group alignments by reference position
    for align in cpu_alignments:
        event_idx = align['event_idx']
        if event_idx < len(events['starts']):
            start = events['starts'][event_idx]
            end = start + events['lengths'][event_idx]
            pos = align['ref_position']

            if start < max_samples:
                # Draw line from reference position to event
                ax.plot([start, end], [pos, pos], 'b-', alpha=0.3, linewidth=1)
                ax.scatter(start, pos, c='blue', s=10, alpha=0.5)

    ax.set_ylabel('Reference Position', fontsize=11, fontweight='bold')
    ax.set_title('CPU Eventalign - Reference to Event Mapping', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Plot 2: f5c alignments
    ax = axes[1]
    ax.plot(sample_indices, signal[:max_samples], color='lightgray', alpha=0.3, linewidth=0.5)

    for align in f5c_aligns:
        event_idx = align['event_idx']
        if event_idx < len(events['starts']):
            start = events['starts'][event_idx]
            end = start + events['lengths'][event_idx]
            pos = align['ref_position']

            if start < max_samples:
                ax.plot([start, end], [pos, pos], 'g-', alpha=0.3, linewidth=1)
                ax.scatter(start, pos, c='green', s=10, alpha=0.5)

    ax.set_ylabel('Reference Position', fontsize=11, fontweight='bold')
    ax.set_xlabel('Sample Index', fontsize=11, fontweight='bold')
    ax.set_title('F5C Eventalign - Reference to Event Mapping', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)

    plt.suptitle('Alignment Trace: Reference Positions', fontsize=14, fontweight='bold')
    plt.tight_layout()

    if output_path is None:
        output_path = Path(__file__).parent / "cpu_f5c_trace.png"

    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Trace plot saved to: {output_path}")


def print_statistics(cpu_result: dict, f5c_alignments: list, n_matching: int, n_cpu_only: int, n_f5c_only: int):
    """Print detailed comparison statistics."""
    cpu_alignments = cpu_result['full'][0][0]
    f5c_aligns = f5c_alignments_to_dict(f5c_alignments)

    print("\n" + "=" * 70)
    print("ALIGNMENT COMPARISON STATISTICS")
    print("=" * 70)

    print(f"\nCPU eventalign: {len(cpu_alignments)} alignments")
    print(f"F5C eventalign: {len(f5c_aligns)} alignments")

    print(f"\nMatching alignments: {n_matching}")
    print(f"CPU only: {n_cpu_only}")
    print(f"F5C only: {n_f5c_only}")

    if len(cpu_alignments) > 0:
        match_rate = n_matching / len(cpu_alignments) * 100
        print(f"\nMatch rate (vs CPU): {match_rate:.2f}%")

    if len(f5c_aligns) > 0:
        match_rate = n_matching / len(f5c_aligns) * 100
        print(f"Match rate (vs F5C): {match_rate:.2f}%")

    # Compare scalings
    cpu_scalings = cpu_result['scalings'][0]
    print(f"\nCPU scalings: scale={cpu_scalings['scale']:.4f}, shift={cpu_scalings['shift']:.4f}, var={cpu_scalings['var']:.4f}")

    # Compare event statistics
    events = cpu_result['events'][0]
    print(f"\nEvents detected: {len(events['starts'])}")
    print(f"  Mean length: {np.mean(events['lengths']):.2f} samples")
    print(f"  Mean signal: {np.mean(events['means']):.2f} pA")
    print(f"  Mean stdv: {np.mean(events['stdvs']):.4f}")

    # Show some detailed comparisons for matching alignments
    if n_matching > 0:
        print("\n" + "-" * 70)
        print("SAMPLE MATCHING ALIGNMENTS (first 10)")
        print("-" * 70)
        print(f"{'RefPos':>8} {'CPU_Event':>10} {'CPU_Kmer':>10} {'F5C_Event':>10} {'F5C_Kmer':>10}")
        print("-" * 60)

        count = 0
        for cpu_a in cpu_alignments:
            key = (cpu_a['ref_position'], cpu_a['event_idx'])
            if key in [(c['ref_position'], c['event_idx']) for c in f5c_aligns]:
                f5c_a = next(f for f in f5c_aligns if f['ref_position'] == cpu_a['ref_position'] and f['event_idx'] == cpu_a['event_idx'])
                print(f"{cpu_a['ref_position']:8d} {cpu_a['event_idx']:10d} {cpu_a['ref_kmer']:>10} {f5c_a['event_idx']:10d} {f5c_a['ref_kmer']:>10}")
                count += 1
                if count >= 10:
                    break


def main():
    """Main comparison function."""
    print("=" * 70)
    print("CPU vs F5C Eventalign Comparison")
    print("=" * 70)

    # Paths to test data
    test_dir = Path(__file__).parent / "test_data"
    fasta_path = test_dir / "one_read.fa"
    pod5_path = test_dir / "one_read.pod5"
    eventalign_path = test_dir / "one_read.eventalign"

    # Check if files exist
    if not fasta_path.exists():
        print(f"FASTA file not found: {fasta_path}")
        print("Please ensure test data is available")
        return

    if not eventalign_path.exists():
        print(f"F5C eventalign file not found: {eventalign_path}")
        print("Will skip F5C comparison")
        return

    # Load reference
    print(f"\nLoading reference from: {fasta_path}")
    ref_name, ref_seq, ref_len = load_reference_from_fasta(fasta_path)
    print(f"  Reference: {ref_name}")
    print(f"  Length: {ref_len} bp")

    # Load signal
    print(f"\nLoading signal from: {pod5_path}")
    read_id, signal, sample_rate = load_signal_from_pod5(pod5_path)
    print(f"  Read ID: {read_id}")
    print(f"  Signal length: {len(signal)} samples")
    print(f"  Sample rate: {sample_rate} Hz")

    # Use reference as read sequence (placeholder)
    read_seq = ref_seq

    # Run CPU eventalign
    print("\n" + "-" * 70)
    print("RUNNING CPU EVENTALIGN")
    print("-" * 70)
    cpu_result = run_eventalign(
        read_ids=[read_id],
        read_seqs=[read_seq],
        ref_seqs=[ref_seq],
        ref_names=[ref_name],
        ref_lens=[ref_len],
        signals=[signal],
        sample_rates=[sample_rate],
        model_id=MODEL_RNA002,
    )
    print(f"CPU complete: {len(cpu_result['full'][0][0])} alignments")

    # Load f5c eventalign results
    print("\n" + "-" * 70)
    print("LOADING F5C EVENTALIGN RESULTS")
    print("-" * 70)
    f5c_alignments = load_f5c_eventalign_results(eventalign_path)
    print(f"F5C loaded: {len(f5c_alignments)} alignments")

    # Generate comparison plots
    print("\n" + "-" * 70)
    print("GENERATING PLOTS")
    print("-" * 70)

    n_matching, n_cpu_only, n_f5c_only = plot_comparison(signal, cpu_result, f5c_alignments)
    plot_alignment_trace(signal, cpu_result, f5c_alignments)

    # Print statistics
    print_statistics(cpu_result, f5c_alignments, n_matching, n_cpu_only, n_f5c_only)

    print("\n" + "=" * 70)
    print("COMPARISON COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
