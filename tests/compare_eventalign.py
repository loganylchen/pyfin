#!/usr/bin/env python3
"""
Compare CPU eventalign, GPU eventalign, and original f5c eventalign results.

This script:
1. Runs eventalign using CPU implementation
2. Runs eventalign using GPU (CUDA) implementation
3. Loads original f5c eventalign results
4. Plots events with raw signals and connects alignment locations
5. Compares performance (timing) between methods
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import time
import subprocess

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from fin._eventalign import run_eventalign, run_eventalign_cuda, MODEL_RNA002, CUDA_AVAILABLE

# For loading signal from POD5 files
try:
    from fin.io.io_pod5 import Pod5Reader
    POD5_AVAILABLE = True
except ImportError:
    POD5_AVAILABLE = False
    print("Warning: pod5 not available")


def run_f5c_eventalign(bam_path, reference_path, output_path):
    """
    Run original f5c eventalign and measure time.

    Returns:
        tuple: (elapsed_time, return_code)
    """
    cmd = [
        "f5c",
        "eventalign",
        "--bam", bam_path,
        "--reference", reference_path,
        "--signal-index",
        "--scale-events",
        "--summary",
        "--out", output_path
    ]

    print(f"Running f5c eventalign...")
    print(f"  Command: {' '.join(cmd)}")

    start_time = time.time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        elapsed_time = time.time() - start_time
        return elapsed_time, result.returncode
    except FileNotFoundError:
        print("  f5c command not found, skipping f5c comparison")
        return None, -1
    except subprocess.TimeoutExpired:
        print("  f5c timed out after 300 seconds")
        return None, -1
    except Exception as e:
        print(f"  Error running f5c: {e}")
        return None, -1


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
        # Generate synthetic signal for testing
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
    Load original f5c eventalign results from .tsv file.

    Returns:
        dict with event alignments organized by reference position
    """
    alignments = {}
    with open(eventalign_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            # Parse f5c eventalign output format
            contig = parts[0]
            ref_position = int(parts[1])
            ref_kmer = parts[2]
            strand = parts[3]
            event_idx = int(parts[4])
            event_mean = float(parts[5])
            event_stdv = float(parts[6])
            model_kmer = parts[7]
            hmm_state = parts[8]

            key = (contig, ref_position)
            alignments[key] = {
                'ref_position': ref_position,
                'ref_kmer': ref_kmer,
                'event_idx': event_idx,
                'event_mean': event_mean,
                'event_stdv': event_stdv,
                'model_kmer': model_kmer,
                'hmm_state': hmm_state,
            }

    return alignments


def plot_comparison(signal, cpu_result, gpu_result, f5c_alignments, events):
    """
    Plot comparison of CPU, GPU, and f5c eventalign results.

    Args:
        signal: raw signal array
        cpu_result: results from run_eventalign (CPU)
        gpu_result: results from run_eventalign_cuda (GPU)
        f5c_alignments: results from original f5c
        events: detected events
    """
    n_events = len(events['means'])

    # Create figure with subplots
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(4, 1, hspace=0.3, height_ratios=[1, 1, 1, 0.5])

    # Sample indices for plotting
    sample_starts = events['starts']
    sample_ends = sample_starts + events['lengths']

    # Plot 1: Raw signal with events (CPU)
    ax1 = fig.add_subplot(gs[0])
    ax1.plot(signal, color='lightgray', alpha=0.5, linewidth=0.5, label='Raw Signal')

    # Mark events
    cpu_full = cpu_result['full'][0][0] if cpu_result['full'][0][0] else []
    cpu_event_positions = {}  # event_idx -> ref_position

    for align in cpu_full:
        event_idx = align['event_idx']
        ref_pos = align['ref_position']
        cpu_event_positions[event_idx] = ref_pos

        # Highlight aligned events
        if event_idx < n_events:
            start = sample_starts[event_idx]
            end = min(sample_ends[event_idx], len(signal))
            ax1.axvspan(start, end, alpha=0.3, color='blue', linewidth=0)

    # Plot event means
    event_means = events['means']
    event_positions = (sample_starts + sample_ends) / 2
    ax1.scatter(event_positions, event_means, c='blue', s=10, alpha=0.6, label='Event Means', zorder=3)

    ax1.set_ylabel('Signal (pA)', fontsize=10)
    ax1.set_title('CPU Eventalign - Events on Raw Signal', fontsize=12, fontweight='bold')
    ax1.legend(loc='upper right', fontsize=8)
    ax1.set_xlim(0, min(len(signal), 50000))

    # Plot 2: Raw signal with events (GPU)
    ax2 = fig.add_subplot(gs[1])
    ax2.plot(signal, color='lightgray', alpha=0.5, linewidth=0.5, label='Raw Signal')

    gpu_full = gpu_result['full'][0][0] if gpu_result['full'][0][0] else []
    gpu_event_positions = {}

    for align in gpu_full:
        event_idx = align['event_idx']
        ref_pos = align['ref_position']
        gpu_event_positions[event_idx] = ref_pos

        if event_idx < n_events:
            start = sample_starts[event_idx]
            end = min(sample_ends[event_idx], len(signal))
            ax2.axvspan(start, end, alpha=0.3, color='green', linewidth=0)

    ax2.scatter(event_positions, event_means, c='green', s=10, alpha=0.6, label='Event Means', zorder=3)

    ax2.set_ylabel('Signal (pA)', fontsize=10)
    ax2.set_title('GPU Eventalign - Events on Raw Signal', fontsize=12, fontweight='bold')
    ax2.legend(loc='upper right', fontsize=8)
    ax2.set_xlim(0, min(len(signal), 50000))

    # Plot 3: Raw signal with events (f5c)
    ax3 = fig.add_subplot(gs[2])
    ax3.plot(signal, color='lightgray', alpha=0.5, linewidth=0.5, label='Raw Signal')

    f5c_event_positions = {}
    for align in f5c_alignments.values():
        event_idx = align['event_idx']
        ref_pos = align['ref_position']
        f5c_event_positions[event_idx] = ref_pos

        if event_idx < n_events:
            start = sample_starts[event_idx]
            end = min(sample_ends[event_idx], len(signal))
            ax3.axvspan(start, end, alpha=0.3, color='red', linewidth=0)

    ax3.scatter(event_positions, event_means, c='red', s=10, alpha=0.6, label='Event Means', zorder=3)

    ax3.set_ylabel('Signal (pA)', fontsize=10)
    ax3.set_xlabel('Sample Index', fontsize=10)
    ax3.set_title('F5C Eventalign - Events on Raw Signal', fontsize=12, fontweight='bold')
    ax3.legend(loc='upper right', fontsize=8)
    ax3.set_xlim(0, min(len(signal), 50000))

    # Plot 4: Alignment comparison table
    ax4 = fig.add_subplot(gs[3])
    ax4.axis('off')

    # Create comparison table
    comparison_data = []
    max_rows = 10  # Show first 10 alignments

    # Get CPU alignments
    cpu_aligns = cpu_result['full'][0][0] if cpu_result['full'][0][0] else []
    gpu_aligns = gpu_result['full'][0][0] if gpu_result['full'][0][0] else []

    for i in range(min(max_rows, len(cpu_aligns))):
        cpu_align = cpu_aligns[i] if i < len(cpu_aligns) else None
        gpu_align = gpu_aligns[i] if i < len(gpu_aligns) else None

        if cpu_align:
            comparison_data.append([
                i,
                cpu_align['event_idx'],
                cpu_align['ref_position'],
                cpu_align['ref_kmer'],
                gpu_align['ref_position'] if gpu_align else 'N/A',
                'MATCH' if gpu_align and gpu_align['ref_position'] == cpu_align['ref_position'] else 'DIFF',
            ])

    # Create table
    table = ax4.table(cellText=comparison_data,
                      colLabels=['Index', 'Event', 'CPU RefPos', 'RefKmer', 'GPU RefPos', 'Match'],
                      cellLoc='center',
                      loc='center',
                      bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)

    # Color code matches
    for i in range(len(comparison_data)):
        if comparison_data[i][5] == 'MATCH':
            table[(i+1, 5)].set_facecolor('#90EE90')  # Light green
        else:
            table[(i+1, 5)].set_facecolor('#FFB6C1')  # Light red

    plt.suptitle('Event Alignment Comparison: CPU vs GPU vs F5C', fontsize=14, fontweight='bold')
    plt.tight_layout()

    # Save figure
    output_path = Path(__file__).parent / "eventalign_comparison.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {output_path}")

    return fig


def plot_alignment_traces(signal, cpu_result, gpu_result, f5c_alignments, events):
    """
    Plot alignment traces showing event-to-reference mappings.
    """
    fig, axes = plt.subplots(3, 1, figsize=(16, 12), sharex=True)

    sample_starts = events['starts']
    sample_ends = sample_starts + events['lengths']
    event_means = events['means']
    event_positions = (sample_starts + sample_ends) / 2
    n_events = len(sample_starts)

    # CPU alignment trace
    ax = axes[0]
    ax.plot(signal, color='lightgray', alpha=0.3, linewidth=0.5)

    cpu_full = cpu_result['full'][0][0] if cpu_result['full'][0][0] else []
    for align in cpu_full:
        event_idx = align['event_idx']
        ref_pos = align['ref_position']

        if event_idx < n_events:
            start = sample_starts[event_idx]
            end = min(sample_ends[event_idx], len(signal))

            # Draw vertical line to reference position
            ax.plot([start, end], [event_means[event_idx], event_means[event_idx]],
                   color='blue', linewidth=2, alpha=0.7)
            ax.text(start, event_means[event_idx] + 5, str(ref_pos),
                   fontsize=6, color='blue', rotation=90)

    ax.scatter(event_positions, event_means, c='blue', s=15, zorder=3)
    ax.set_ylabel('Signal (pA)', fontsize=10)
    ax.set_title('CPU: Event Reference Positions', fontsize=12, fontweight='bold')
    ax.set_xlim(0, min(len(signal), 50000))

    # GPU alignment trace
    ax = axes[1]
    ax.plot(signal, color='lightgray', alpha=0.3, linewidth=0.5)

    gpu_full = gpu_result['full'][0][0] if gpu_result['full'][0][0] else []
    for align in gpu_full:
        event_idx = align['event_idx']
        ref_pos = align['ref_position']

        if event_idx < n_events:
            start = sample_starts[event_idx]
            end = min(sample_ends[event_idx], len(signal))

            ax.plot([start, end], [event_means[event_idx], event_means[event_idx]],
                   color='green', linewidth=2, alpha=0.7)
            ax.text(start, event_means[event_idx] + 5, str(ref_pos),
                   fontsize=6, color='green', rotation=90)

    ax.scatter(event_positions, event_means, c='green', s=15, zorder=3)
    ax.set_ylabel('Signal (pA)', fontsize=10)
    ax.set_title('GPU: Event Reference Positions', fontsize=12, fontweight='bold')
    ax.set_xlim(0, min(len(signal), 50000))

    # F5C alignment trace
    ax = axes[2]
    ax.plot(signal, color='lightgray', alpha=0.3, linewidth=0.5)

    for align in f5c_alignments.values():
        event_idx = align['event_idx']
        ref_pos = align['ref_position']

        if event_idx < n_events:
            start = sample_starts[event_idx]
            end = min(sample_ends[event_idx], len(signal))

            ax.plot([start, end], [event_means[event_idx], event_means[event_idx]],
                   color='red', linewidth=2, alpha=0.7)
            ax.text(start, event_means[event_idx] + 5, str(ref_pos),
                   fontsize=6, color='red', rotation=90)

    ax.scatter(event_positions, event_means, c='red', s=15, zorder=3)
    ax.set_ylabel('Signal (pA)', fontsize=10)
    ax.set_xlabel('Sample Index', fontsize=10)
    ax.set_title('F5C: Event Reference Positions', fontsize=12, fontweight='bold')
    ax.set_xlim(0, min(len(signal), 50000))

    plt.suptitle('Event Alignment Traces: Reference Positions on Signal', fontsize=14, fontweight='bold')
    plt.tight_layout()

    output_path = Path(__file__).parent / "eventalign_traces.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Traces plot saved to: {output_path}")


def compare_results(cpu_result, gpu_result, f5c_alignments, cpu_time, gpu_time, f5c_time):
    """Print detailed comparison statistics."""
    print("\n" + "="*70)
    print("ALIGNMENT COMPARISON STATISTICS")
    print("="*70)

    cpu_full = cpu_result['full'][0][0] if cpu_result['full'][0][0] else []
    gpu_full = gpu_result['full'][0][0] if gpu_result['full'][0][0] else []

    print(f"\nCPU alignments: {len(cpu_full)}")
    print(f"GPU alignments: {len(gpu_full)}")
    print(f"F5C alignments: {len(f5c_alignments)}")

    # Compare CPU vs GPU
    cpu_gpu_matches = 0
    for cpu_align in cpu_full:
        for gpu_align in gpu_full:
            if (cpu_align['event_idx'] == gpu_align['event_idx'] and
                cpu_align['ref_position'] == gpu_align['ref_position']):
                cpu_gpu_matches += 1
                break

    print(f"\nCPU vs GPU matches: {cpu_gpu_matches}/{len(cpu_full)} ({100*cpu_gpu_matches/max(len(cpu_full),1):.1f}%)")

    # Compare CPU vs F5C
    cpu_f5c_matches = 0
    for cpu_align in cpu_full:
        for f5c_align in f5c_alignments.values():
            if (cpu_align['event_idx'] == f5c_align['event_idx'] and
                cpu_align['ref_position'] == f5c_align['ref_position']):
                cpu_f5c_matches += 1
                break

    print(f"CPU vs F5C matches: {cpu_f5c_matches}/{len(cpu_full)} ({100*cpu_f5c_matches/max(len(cpu_full),1):.1f}%)")

    # Scalings comparison
    print("\n" + "-"*70)
    print("SCALINGS COMPARISON")
    print("-"*70)

    cpu_scalings = cpu_result['scalings'][0]
    gpu_scalings = gpu_result['scalings'][0]

    print(f"\nCPU: scale={cpu_scalings['scale']:.4f}, shift={cpu_scalings['shift']:.4f}, var={cpu_scalings['var']:.4f}")
    print(f"GPU: scale={gpu_scalings['scale']:.4f}, shift={gpu_scalings['shift']:.4f}, var={gpu_scalings['var']:.4f}")

    # Timing comparison
    print("\n" + "-"*70)
    print("PERFORMANCE COMPARISON")
    print("-"*70)

    print(f"\nCPU time:  {cpu_time:.4f} seconds")
    print(f"GPU time:  {gpu_time:.4f} seconds")
    if f5c_time is not None:
        print(f"F5C time:  {f5c_time:.4f} seconds")

        # Calculate speedups
        if gpu_time > 0:
            gpu_speedup = cpu_time / gpu_time
            print(f"\nGPU vs CPU speedup: {gpu_speedup:.2f}x")
        if f5c_time > 0:
            cpu_vs_f5c = f5c_time / cpu_time
            gpu_vs_f5c = f5c_time / gpu_time
            print(f"CPU vs F5C speedup: {cpu_vs_f5c:.2f}x")
            print(f"GPU vs F5C speedup: {gpu_vs_f5c:.2f}x")
    else:
        if gpu_time > 0:
            gpu_speedup = cpu_time / gpu_time
            print(f"\nGPU vs CPU speedup: {gpu_speedup:.2f}x")


def plot_timing_comparison(cpu_time, gpu_time, f5c_time):
    """Plot timing comparison bar chart."""
    fig, ax = plt.subplots(figsize=(10, 6))

    methods = []
    times = []
    colors = []

    if cpu_time is not None:
        methods.append('CPU')
        times.append(cpu_time)
        colors.append('#3498db')  # Blue

    if gpu_time is not None:
        methods.append('GPU (CUDA)')
        times.append(gpu_time)
        colors.append('#2ecc71')  # Green

    if f5c_time is not None:
        methods.append('F5C')
        times.append(f5c_time)
        colors.append('#e74c3c')  # Red

    bars = ax.bar(methods, times, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)

    # Add value labels on bars
    for bar, time in zip(bars, times):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
               f'{time:.4f}s',
               ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax.set_ylabel('Time (seconds)', fontsize=12, fontweight='bold')
    ax.set_title('Event Alignment Performance Comparison', fontsize=14, fontweight='bold')
    ax.grid(axis='y', alpha=0.3)

    # Add speedup annotations
    if cpu_time and gpu_time:
        speedup = cpu_time / gpu_time
        ax.annotate(f'GPU is\n{speedup:.1f}x faster',
                   xy=(1, gpu_time),
                   xytext=(0.5, (cpu_time + gpu_time) / 2),
                   arrowprops=dict(arrowstyle='->', color='black', lw=1.5),
                   fontsize=10, ha='center', fontweight='bold')

    plt.tight_layout()

    output_path = Path(__file__).parent / "eventalign_timing.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Timing plot saved to: {output_path}")


def main():
    """Main comparison function."""
    print("="*70)
    print("EVENTALIGN COMPARISON: CPU vs GPU vs F5C")
    print("="*70)

    # Paths to test data
    test_dir = Path(__file__).parent / "test_data"
    fasta_path = test_dir / "one_read.fa"
    pod5_path = test_dir / "one_read.pod5"
    eventalign_path = test_dir / "one_read.eventalign.tsv"

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

    # Use reference as read sequence (in practice, use actual basecalled sequence)
    read_seq = ref_seq

    print("\n" + "-"*70)
    print("RUNNING CPU EVENTALIGN")
    print("-"*70)
    cpu_start = time.time()
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
    cpu_time = time.time() - cpu_start
    print(f"CPU complete: {len(cpu_result['full'][0][0])} alignments in {cpu_time:.4f}s")

    print("\n" + "-"*70)
    print("RUNNING GPU EVENTALIGN")
    print("-"*70)
    gpu_time = None
    if CUDA_AVAILABLE:
        gpu_start = time.time()
        gpu_result = run_eventalign_cuda(
            read_ids=[read_id],
            read_seqs=[read_seq],
            ref_seqs=[ref_seq],
            ref_names=[ref_name],
            ref_lens=[ref_len],
            signals=[signal],
            sample_rates=[sample_rate],
            model_id=MODEL_RNA002,
        )
        gpu_time = time.time() - gpu_start
        print(f"GPU complete: {len(gpu_result['full'][0][0])} alignments in {gpu_time:.4f}s")
    else:
        print("CUDA not available, using CPU result as placeholder")
        gpu_result = cpu_result
        gpu_time = cpu_time

    # Load f5c results
    print("\n" + "-"*70)
    print("LOADING F5C EVENTALIGN RESULTS")
    print("-"*70)
    f5c_time = None
    if eventalign_path.exists():
        f5c_alignments = load_f5c_eventalign_results(eventalign_path)
        print(f"F5C complete: {len(f5c_alignments)} alignments")
    else:
        print(f"F5C eventalign file not found: {eventalign_path}")
        print("Will skip F5C comparison")
        f5c_alignments = {}

    # Optionally run f5c and measure time
    test_dir = Path(__file__).parent / "test_data"
    bam_path = test_dir / "one_read.bam"
    if bam_path.exists() and fasta_path.exists():
        print("\n" + "-"*70)
        print("RUNNING F5C EVENTALIGN (for timing)")
        print("-"*70)
        f5c_time, _ = run_f5c_eventalign(str(bam_path), str(fasta_path), str(test_dir / "f5c_output.eventalign"))

    # Get events
    events = cpu_result['events'][0]
    n_events = len(events['starts'])
    print(f"\nTotal events detected: {n_events}")

    # Print comparison statistics
    compare_results(cpu_result, gpu_result, f5c_alignments, cpu_time, gpu_time, f5c_time)

    # Generate plots
    print("\n" + "-"*70)
    print("GENERATING PLOTS")
    print("-"*70)

    # Limit signal for plotting
    plot_signal = signal[:50000]

    # Filter events for plotting region
    plot_mask = events['starts'] < 50000
    plot_events = {
        'starts': events['starts'][plot_mask],
        'lengths': events['lengths'][plot_mask],
        'means': events['means'][plot_mask],
        'stdvs': events['stdvs'][plot_mask],
    }

    plot_comparison(plot_signal, cpu_result, gpu_result, f5c_alignments, plot_events)
    plot_alignment_traces(plot_signal, cpu_result, gpu_result, f5c_alignments, plot_events)
    plot_timing_comparison(cpu_time, gpu_time, f5c_time)

    print("\n" + "="*70)
    print("COMPARISON COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()
