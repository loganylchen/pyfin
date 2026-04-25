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
import gzip

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
    Load f5c eventalign results from gzipped TSV file.

    The f5c eventalign output format (tsv):
    chromosome    position    reference_kmer    read_index    strand    event_idx
    event_mean    event_stdv    log_prob    model_kmer    model_mean    model_stdv
    scale    start_idx    end_idx    signal_samples

    Format (0-indexed, 16 columns total):
    chromosome(0), position(1), reference_kmer(2), read_index(3), strand(4),
    event_idx(5), event_mean(6), event_stdv(7), log_prob(8), model_kmer(9),
    model_mean(10), model_stdv(11), scale(12), start_idx(13), end_idx(14),
    signal_samples(15)
    """
    alignments = []

    # Open file (handle both .gz and plain text)
    open_func = gzip.open if str(eventalign_path).endswith('.gz') else open
    mode = 'rt' if str(eventalign_path).endswith('.gz') else 'r'

    with open_func(eventalign_path, mode) as f:
        for line in f:
            if line.startswith("chromosome"):
                continue  # Skip header

            parts = line.strip().split("\t")
            if len(parts) < 15:
                continue

            alignment = {
                "ref_position": int(parts[1]),
                "ref_kmer": parts[2],
                "event_idx": int(parts[5]),
                "event_mean": float(parts[6]),
                "event_stdv": float(parts[7]),
                "log_prob": float(parts[8]),  # Log probability
                "model_kmer": parts[9],  # String, not float
                "model_mean": float(parts[10]),
                "model_stdv": float(parts[11]),
                "scale": float(parts[12]),  # Scale/drift value
                "start": int(parts[13]),
                "end": int(parts[14]) if len(parts) > 14 else None,
            }
            alignments.append(alignment)

    return alignments


def f5c_alignments_to_dict(f5c_alignments: list, events: dict) -> list:
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
    f5c_aligns = f5c_alignments_to_dict(f5c_alignments, events)

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

    # Plot CPU-only events (limited to avoid overcrowding)
    for i, align in enumerate(cpu_only[:50]):
        event_idx = align['event_idx']
        if event_idx < len(events['starts']):
            start = events['starts'][event_idx]
            if start < max_samples:
                ax.scatter(start, events['means'][event_idx], c='blue', s=20, alpha=0.5,
                          marker='x', label='CPU only' if i == 0 else "")

    # Plot f5c-only events (limited to avoid overcrowding)
    for i, align in enumerate(f5c_only[:50]):
        event_idx = align['event_idx']
        if event_idx < len(events['starts']):
            start = events['starts'][event_idx]
            if start < max_samples:
                ax.scatter(start, align['event_mean'], c='red', s=20, alpha=0.5,
                          marker='+', label='F5C only' if i == 0 else "")

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
    Plot alignment trace showing reference positions vs event index.

    Creates a visualization showing which reference k-mers align to which events,
    with both CPU and f5c results overlaid in different colors.
    """
    events = cpu_result['events'][0]
    cpu_alignments = cpu_result['full'][0][0]
    f5c_aligns = f5c_alignments_to_dict(f5c_alignments, events)

    fig, ax = plt.subplots(1, 1, figsize=(16, 10))

    # Plot CPU alignments as blue points
    for align in cpu_alignments:
        event_idx = align['event_idx']
        pos = align['ref_position']
        ax.scatter(event_idx, pos, c='blue', s=20, alpha=0.6, marker='o', edgecolors='none', label='CPU' if align == cpu_alignments[0] else "")

    # Plot f5c alignments as green points
    for align in f5c_aligns:
        event_idx = align['event_idx']
        pos = align['ref_position']
        ax.scatter(event_idx, pos, c='green', s=20, alpha=0.6, marker='x', label='F5C' if align == f5c_aligns[0] else "")

    # Find and highlight matching alignments
    cpu_positions = {(a['ref_position'], a['event_idx']): a for a in cpu_alignments}
    f5c_positions = {(a['ref_position'], a['event_idx']): a for a in f5c_aligns}

    for key in cpu_positions:
        if key in f5c_positions:
            pos, event_idx = key
            ax.scatter(event_idx, pos, c='red', s=30, alpha=0.8, marker='o', edgecolors='black', linewidth=1, label='Matching' if key == list(cpu_positions.keys())[0] else "")

    ax.set_xlabel('Event Index', fontsize=11, fontweight='bold')
    ax.set_ylabel('Reference Position', fontsize=11, fontweight='bold')
    ax.set_title('Alignment Trace: Reference Position vs Event Index', fontsize=12, fontweight='bold')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if output_path is None:
        output_path = Path(__file__).parent / "cpu_f5c_trace.png"

    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Trace plot saved to: {output_path}")


def print_statistics(cpu_result: dict, f5c_alignments: list, n_matching: int, n_cpu_only: int, n_f5c_only: int):
    """Print detailed comparison statistics."""
    cpu_alignments = cpu_result['full'][0][0]
    f5c_aligns = f5c_alignments_to_dict(f5c_alignments, cpu_result['events'][0])

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

    # Get f5c scale info
    if f5c_aligns:
        f5c_scales = [a.get('scale', 0) for a in f5c_alignments]
        print(f"F5C scales: mean={np.mean(f5c_scales):.4f}, std={np.std(f5c_scales):.4f}")

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
            f5c_match = next((f for f in f5c_aligns if f['ref_position'] == cpu_a['ref_position'] and f['event_idx'] == cpu_a['event_idx']), None)
            if f5c_match:
                print(f"{cpu_a['ref_position']:8d} {cpu_a['event_idx']:10d} {cpu_a['ref_kmer']:>10} {f5c_match['event_idx']:10d} {f5c_match['ref_kmer']:>10}")
                count += 1
                if count >= 10:
                    break


def evaluate_alignment_quality(
    cpu_result: dict,
    f5c_alignments: list,
    signal: np.ndarray,
) -> dict:
    """
    Evaluate which eventalign result is better based on multiple metrics.

    Metrics:
    1. Coverage: How many reference positions are covered
    2. Continuity: How continuous the alignment is (fewer gaps)
    3. Consistency: Log probability scores (higher is better)
    4. Event utilization: How well events are used
    """
    cpu_alignments = cpu_result['full'][0][0]
    f5c_aligns = f5c_alignments_to_dict(f5c_alignments, cpu_result['events'][0])
    events = cpu_result['events'][0]

    metrics = {
        'cpu': {},
        'f5c': {},
        'winner': {}
    }

    # 1. Coverage - number of unique reference positions covered
    cpu_ref_positions = set(a['ref_position'] for a in cpu_alignments)
    f5c_ref_positions = set(a['ref_position'] for a in f5c_aligns)

    metrics['cpu']['coverage'] = len(cpu_ref_positions)
    metrics['f5c']['coverage'] = len(f5c_ref_positions)

    if metrics['cpu']['coverage'] > metrics['f5c']['coverage']:
        metrics['winner']['coverage'] = 'CPU'
        metrics['winner']['coverage_margin'] = metrics['cpu']['coverage'] - metrics['f5c']['coverage']
    elif metrics['f5c']['coverage'] > metrics['cpu']['coverage']:
        metrics['winner']['coverage'] = 'F5C'
        metrics['winner']['coverage_margin'] = metrics['f5c']['coverage'] - metrics['cpu']['coverage']
    else:
        metrics['winner']['coverage'] = 'TIE'

    # 2. Continuity - measure gaps in alignment
    if len(cpu_alignments) > 1:
        cpu_sorted = sorted(cpu_alignments, key=lambda x: x['ref_position'])
        cpu_gaps = [cpu_sorted[i+1]['ref_position'] - cpu_sorted[i]['ref_position']
                    for i in range(len(cpu_sorted)-1)]
        metrics['cpu']['mean_gap'] = np.mean(cpu_gaps)
        metrics['cpu']['max_gap'] = max(cpu_gaps) if cpu_gaps else 0
        metrics['cpu']['gap_std'] = np.std(cpu_gaps)

    if len(f5c_aligns) > 1:
        f5c_sorted = sorted(f5c_aligns, key=lambda x: x['ref_position'])
        f5c_gaps = [f5c_sorted[i+1]['ref_position'] - f5c_sorted[i]['ref_position']
                    for i in range(len(f5c_sorted)-1)]
        metrics['f5c']['mean_gap'] = np.mean(f5c_gaps)
        metrics['f5c']['max_gap'] = max(f5c_gaps) if f5c_gaps else 0
        metrics['f5c']['gap_std'] = np.std(f5c_gaps)

    # Lower mean gap is better (more continuous)
    if 'mean_gap' in metrics['cpu'] and 'mean_gap' in metrics['f5c']:
        if metrics['cpu']['mean_gap'] < metrics['f5c']['mean_gap']:
            metrics['winner']['continuity'] = 'CPU'
        elif metrics['f5c']['mean_gap'] < metrics['cpu']['mean_gap']:
            metrics['winner']['continuity'] = 'F5C'
        else:
            metrics['winner']['continuity'] = 'TIE'

    # 3. Consistency - log probability scores (higher is better)
    if f5c_aligns and 'log_prob' in f5c_alignments[0]:
        f5c_log_probs = [a.get('log_prob', -float('inf')) for a in f5c_alignments]
        metrics['f5c']['mean_log_prob'] = np.mean(f5c_log_probs)
        metrics['f5c']['min_log_prob'] = min(f5c_log_probs)
        metrics['f5c']['max_log_prob'] = max(f5c_log_probs)

    # CPU log probabilities from alignment results
    if cpu_alignments and 'log_prob' in cpu_alignments[0]:
        cpu_log_probs = [a.get('log_prob', -float('inf')) for a in cpu_alignments]
        metrics['cpu']['mean_log_prob'] = np.mean(cpu_log_probs)
        metrics['cpu']['min_log_prob'] = min(cpu_log_probs)
        metrics['cpu']['max_log_prob'] = max(cpu_log_probs)

        # Compare log probabilities between CPU and F5C
        if 'mean_log_prob' in metrics['f5c']:
            if metrics['cpu']['mean_log_prob'] > metrics['f5c']['mean_log_prob']:
                metrics['winner']['log_prob'] = 'CPU'
                metrics['winner']['log_prob_margin'] = metrics['cpu']['mean_log_prob'] - metrics['f5c']['mean_log_prob']
            elif metrics['f5c']['mean_log_prob'] > metrics['cpu']['mean_log_prob']:
                metrics['winner']['log_prob'] = 'F5C'
                metrics['winner']['log_prob_margin'] = metrics['f5c']['mean_log_prob'] - metrics['cpu']['mean_log_prob']
            else:
                metrics['winner']['log_prob'] = 'TIE'

    # CPU doesn't output log_prob directly, so we compare event quality
    # Lower event stdv typically means cleaner events
    if len(events['stdvs']) > 0:
        metrics['cpu']['mean_event_stdv'] = np.mean(events['stdvs'])
        metrics['cpu']['median_event_stdv'] = np.median(events['stdvs'])

    # 4. Event utilization - unique events used / total events
    cpu_events_used = set(a['event_idx'] for a in cpu_alignments)
    f5c_events_used = set(a['event_idx'] for a in f5c_aligns)
    total_events = len(events['starts'])

    metrics['cpu']['event_utilization'] = len(cpu_events_used) / total_events if total_events > 0 else 0
    metrics['f5c']['event_utilization'] = len(f5c_events_used) / total_events if total_events > 0 else 0

    if metrics['cpu']['event_utilization'] > metrics['f5c']['event_utilization']:
        metrics['winner']['event_utilization'] = 'CPU'
    elif metrics['f5c']['event_utilization'] > metrics['cpu']['event_utilization']:
        metrics['winner']['event_utilization'] = 'F5C'
    else:
        metrics['winner']['event_utilization'] = 'TIE'

    # 5. Alignment density - alignments per reference position
    if cpu_ref_positions:
        metrics['cpu']['alignments_per_position'] = len(cpu_alignments) / len(cpu_ref_positions)
    if f5c_ref_positions:
        metrics['f5c']['alignments_per_position'] = len(f5c_aligns) / len(f5c_ref_positions)

    # 6. K-mer consistency - check if aligned k-mers match reference
    cpu_kmer_match = sum(1 for a in cpu_alignments if a['ref_kmer'] == a.get('model_kmer', a['ref_kmer']))
    f5c_kmer_match = sum(1 for a in f5c_aligns if a['ref_kmer'] == a.get('model_kmer', a['ref_kmer']))

    if len(cpu_alignments) > 0:
        metrics['cpu']['kmer_consistency'] = cpu_kmer_match / len(cpu_alignments) * 100
    if len(f5c_aligns) > 0:
        metrics['f5c']['kmer_consistency'] = f5c_kmer_match / len(f5c_aligns) * 100

    return metrics


def print_evaluation_report(metrics: dict):
    """Print a detailed evaluation report."""
    print("\n" + "=" * 70)
    print("ALIGNMENT QUALITY EVALUATION")
    print("=" * 70)

    print("\n1. COVERAGE (Reference positions covered)")
    print("-" * 70)
    print(f"  CPU:     {metrics['cpu']['coverage']} positions")
    print(f"  F5C:     {metrics['f5c']['coverage']} positions")
    winner = metrics['winner'].get('coverage', 'N/A')
    if winner != 'TIE':
        margin = metrics['winner'].get('coverage_margin', 0)
        print(f"  Winner:  {winner} (by {margin} positions)")
    else:
        print("  Winner:  TIE")

    print("\n2. CONTINUITY (Gap analysis)")
    print("-" * 70)
    if 'mean_gap' in metrics['cpu']:
        print(f"  CPU:     mean gap = {metrics['cpu']['mean_gap']:.2f}, "
              f"max gap = {metrics['cpu']['max_gap']}, std = {metrics['cpu']['gap_std']:.2f}")
    if 'mean_gap' in metrics['f5c']:
        print(f"  F5C:     mean gap = {metrics['f5c']['mean_gap']:.2f}, "
              f"max gap = {metrics['f5c']['max_gap']}, std = {metrics['f5c']['gap_std']:.2f}")
    if 'continuity' in metrics['winner']:
        print(f"  Winner:  {metrics['winner']['continuity']} (lower gap is better)")

    print("\n3. EVENT UTILIZATION")
    print("-" * 70)
    print(f"  CPU:     {metrics['cpu']['event_utilization']*100:.1f}% of events used")
    print(f"  F5C:     {metrics['f5c']['event_utilization']*100:.1f}% of events used")
    print(f"  Winner:  {metrics['winner']['event_utilization']}")

    print("\n4. ALIGNMENT DENSITY")
    print("-" * 70)
    if 'alignments_per_position' in metrics['cpu']:
        print(f"  CPU:     {metrics['cpu']['alignments_per_position']:.2f} alignments/position")
    if 'alignments_per_position' in metrics['f5c']:
        print(f"  F5C:     {metrics['f5c']['alignments_per_position']:.2f} alignments/position")

    print("\n5. K-MER CONSISTENCY")
    print("-" * 70)
    if 'kmer_consistency' in metrics['cpu']:
        print(f"  CPU:     {metrics['cpu']['kmer_consistency']:.1f}% consistent")
    if 'kmer_consistency' in metrics['f5c']:
        print(f"  F5C:     {metrics['f5c']['kmer_consistency']:.1f}% consistent")

    print("\n6. LOG PROBABILITIES")
    print("-" * 70)
    if 'mean_log_prob' in metrics['cpu']:
        print(f"  CPU:     mean = {metrics['cpu']['mean_log_prob']:.4f}, "
              f"min = {metrics['cpu']['min_log_prob']:.4f}, "
              f"max = {metrics['cpu']['max_log_prob']:.4f}")
    if 'mean_log_prob' in metrics['f5c']:
        print(f"  F5C:     mean = {metrics['f5c']['mean_log_prob']:.4f}, "
              f"min = {metrics['f5c']['min_log_prob']:.4f}, "
              f"max = {metrics['f5c']['max_log_prob']:.4f}")
    if 'log_prob' in metrics['winner']:
        winner = metrics['winner'].get('log_prob', 'N/A')
        if winner != 'TIE':
            margin = metrics['winner'].get('log_prob_margin', 0)
            print(f"  Winner:  {winner} (by {margin:.4f})")
        else:
            print("  Winner:  TIE")
        print("  (Higher log probability indicates better alignment quality)")

    # Overall verdict
    print("\n" + "=" * 70)
    print("OVERALL VERDICT")
    print("=" * 70)

    wins = {'CPU': 0, 'F5C': 0, 'TIE': 0}
    for key in ['coverage', 'continuity', 'event_utilization']:
        winner = metrics['winner'].get(key, 'TIE')
        wins[winner] = wins.get(winner, 0) + 1

    print(f"\nCPU wins:   {wins['CPU']} categories")
    print(f"F5C wins:   {wins['F5C']} categories")
    print(f"Ties:       {wins['TIE']} categories")

    if wins['CPU'] > wins['F5C']:
        print("\n  CPU eventalign appears better overall.")
    elif wins['F5C'] > wins['CPU']:
        print("\n  F5C eventalign appears better overall.")
    else:
        print("\n  Both methods perform similarly.")

    print("=" * 70)


def main():
    """Main comparison function."""
    print("=" * 70)
    print("CPU vs F5C Eventalign Comparison")
    print("=" * 70)

    # Paths to test data
    test_dir = Path(__file__).parent / "test_data"
    fasta_path = test_dir / "one_read.fa"
    pod5_path = test_dir / "one_read.pod5"
    eventalign_path = test_dir / "one_read.eventalign.tsv.gz"

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

    # Evaluate alignment quality
    print("\n" + "-" * 70)
    print("EVALUATING ALIGNMENT QUALITY")
    print("-" * 70)
    metrics = evaluate_alignment_quality(cpu_result, f5c_alignments, signal)
    print_evaluation_report(metrics)

    print("\n" + "=" * 70)
    print("COMPARISON COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
