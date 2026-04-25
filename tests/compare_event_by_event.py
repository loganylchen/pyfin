#!/usr/bin/env python3
"""
Detailed event-by-event comparison between f5c and pyfin eventalign results.

This script:
1. Loads f5c reference TSV output
2. Runs pyfin eventalign
3. Compares event-by-event with detailed statistics
4. Creates visualizations
"""

import numpy as np
import gzip
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def load_f5c_eventalign_tsv(tsv_path: str) -> List[Dict]:
    """Load f5c eventalign TSV output."""
    alignments = []
    opener = gzip.open if tsv_path.endswith('.gz') else open

    with opener(tsv_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 16:
                continue

            aln = {
                'reference_name': parts[0],
                'reference_position': int(parts[1]),
                'reference_kmer': parts[2],
                'read_id': parts[3],
                'strand': parts[4],
                'event_idx': int(parts[5]),
                'event_mean': float(parts[6]),
                'event_stdv': float(parts[7]),
                'duration': float(parts[8]),
                'model_kmer': parts[9],
                'model_mean': float(parts[10]),
                'model_stdv': float(parts[11]),
                'scaled_mean': float(parts[12]),
                'start_idx': int(parts[13]),
                'end_idx': int(parts[14]),
            }

            if len(parts) > 15:
                aln['raw_samples'] = np.array([float(x) for x in parts[15].split(',')])

            alignments.append(aln)

    return alignments


def load_reference_from_fasta(fasta_path: str) -> Tuple[str, str, int]:
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


def load_read_sequence_from_fastq(fastq_path: str) -> Tuple[str, str]:
    """Load read sequence from FASTQ file."""
    with open(fastq_path, "r") as f:
        lines = f.readlines()

    read_id = None
    read_sequence = None
    for i, line in enumerate(lines):
        line = line.strip()
        if i == 0 and line.startswith("@"):
            # Parse header: @read_id optional_fields...
            read_id = line[1:].split()[0]
        elif i == 1:
            read_sequence = line.upper()
            break  # Got what we need

    if read_sequence is None:
        raise ValueError(f"Could not parse FASTQ file: {fastq_path}")

    return read_id, read_sequence


def load_signal_from_pod5(pod5_path: str) -> Tuple[str, np.ndarray, float]:
    """Load signal data from POD5 file."""
    try:
        from fin.io.io_pod5 import Pod5Reader
        with Pod5Reader(pod5_path) as reader:
            read_id = reader.read_ids[0]
            read = reader.get_read(read_id)
            signal = read.signal_pa.astype(np.float32)
            sample_rate = float(read.run_info.sample_rate)
            return str(read.read_id), signal, sample_rate
    except Exception as e:
        print(f"POD5 load failed: {e}")
        np.random.seed(42)
        n_samples = 100000
        signal = np.random.randn(n_samples).astype(np.float32) * 10 + 120
        sample_rate = 4000.0
        return "synthetic_read", signal, sample_rate


def get_kmer_rank(kmer: str) -> int:
    """
    Calculate lexicographic rank of a k-mer.

    Rank is computed from last base to first: A=0, C=1, G=2, T=3
    This matches the f5c/nanopolish k-mer ranking scheme.

    Example:
        AAAAA -> 0
        AAAAC -> 1
        AAAAT -> 3
        TTTTT -> 1023 (for k=5)
    """
    base_rank = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    r = 0
    # From last base to first (same as get_kmer_rank in align.cpp)
    for i, base in enumerate(reversed(kmer)):
        r += base_rank.get(base, 0) << (i << 1)
    return r


def export_pyfin_to_f5c_tsv(pyfin_result: Dict, ref_name: str, ref_seq: str,
                            output_path: str, signal: Optional[np.ndarray] = None,
                            include_raw_samples: bool = False,
                            pore_model: Optional[Dict] = None):
    """
    Export pyfin eventalign results in f5c TSV format.

    Args:
        pyfin_result: Result dictionary from run_eventalign
        ref_name: Reference sequence name
        ref_seq: Reference sequence (for k-mer lookup)
        output_path: Output TSV file path (.gz for compressed)
        signal: Raw signal array (optional, for raw_samples column)
        include_raw_samples: Whether to include raw samples column (large files!)
        pore_model: Pore model dict from set_model() with level_means/stdvs arrays
    """
    from fin._eventalign import set_model, MODEL_RNA002

    # Extract data from pyfin result
    pyfin_align = pyfin_result["full"][0][0]
    pyfin_events = pyfin_result["events"][0]
    scalings = pyfin_result["scalings"][0]
    read_id = pyfin_result.get("read_ids", ["unknown"])[0]

    # Load pore model if not provided
    if pore_model is None:
        pore_model = set_model(MODEL_RNA002)

    level_means = pore_model['level_means']
    level_stdvs = pore_model['level_stdvs']

    with (gzip.open if output_path.endswith('.gz') else open)(output_path, 'wt') as f:
        for aln in pyfin_align:
            ref_pos = aln['ref_position']
            event_idx = aln['event_idx']
            ref_kmer = aln['ref_kmer']
            model_kmer = aln['model_kmer']
            hmm_state = aln['hmm_state']

            # Get event statistics
            event_mean = float(pyfin_events['means'][event_idx])
            event_stdv = float(pyfin_events['stdvs'][event_idx])
            event_length = float(pyfin_events['lengths'][event_idx])
            event_start = int(pyfin_events['starts'][event_idx])
            event_end = event_start + int(event_length)

            # Get model mean/stdv from pore model using k-mer rank
            kmer_rank = get_kmer_rank(model_kmer)
            model_mean = float(level_means[kmer_rank])
            model_stdv = float(level_stdvs[kmer_rank])

            # Apply scaling to get scaled mean
            scale = scalings['scale']
            shift = scalings['shift']
            scaled_mean = (event_mean - shift) / scale

            # Strand (for RNA, this is typically '-' for complement)
            strand = '+'

            # Build TSV line (f5c format has 16 columns)
            # Format: ref_name pos ref_kmer read_id strand event_idx event_mean
            #         event_stdv duration model_kmer model_mean model_stdv
            #         scaled_mean start_idx end_idx [raw_samples]
            parts = [
                ref_name,
                str(ref_pos),
                ref_kmer,
                read_id,
                strand,
                str(event_idx),
                f"{event_mean:.4f}",
                f"{event_stdv:.4f}",
                f"{event_length:.4f}",
                model_kmer,
                f"{model_mean:.4f}",
                f"{model_stdv:.4f}",
                f"{scaled_mean:.4f}",
                str(event_start),
                str(event_end),
            ]

            # Add raw samples if requested and signal is available
            if include_raw_samples and signal is not None:
                raw_samples = signal[event_start:event_end]
                raw_samples_str = ','.join(f"{x:.2f}" for x in raw_samples)
                parts.append(raw_samples_str)

            f.write('\t'.join(parts) + '\n')

    print(f"  Exported {len(pyfin_align)} event alignments to: {output_path}")


def index_f5c_by_position(alignments: List[Dict]) -> Dict[int, List[Dict]]:
    """Index f5c alignments by reference position."""
    pos_map = {}
    for aln in alignments:
        pos = aln['reference_position']
        if pos not in pos_map:
            pos_map[pos] = []
        pos_map[pos].append(aln)
    return pos_map


def compare_event_by_event(f5c_alignments: List[Dict], pyfin_result: Dict,
                           ref_seq: str) -> Dict:
    """
    Perform event-by-event comparison between f5c and pyfin.

    Returns detailed comparison results.
    """
    # Index f5c by position
    f5c_by_pos = index_f5c_by_position(f5c_alignments)

    # Extract pyfin alignments
    pyfin_align = pyfin_result["full"][0][0]
    pyfin_events = pyfin_result["events"][0]

    # Index pyfin by position
    pyfin_by_pos = {}
    for aln in pyfin_align:
        pos = aln['ref_position']
        if pos not in pyfin_by_pos:
            pyfin_by_pos[pos] = []
        pyfin_by_pos[pos].append(aln)

    # Get all positions
    all_positions = sorted(set(f5c_by_pos.keys()) | set(pyfin_by_pos.keys()))

    comparison = {
        'position_details': [],
        'f5c_only_positions': [],
        'pyfin_only_positions': [],
        'matched_positions': [],
        'event_mean_diff': [],
        'event_count_diff': [],
        'statistics': {}
    }

    print(f"\n{'='*80}")
    print("Event-by-Event Comparison")
    print(f"{'='*80}")

    # Compare position by position
    for pos in all_positions:
        f5c_events = f5c_by_pos.get(pos, [])
        pyfin_events_at_pos = pyfin_by_pos.get(pos, [])

        pos_detail = {
            'position': pos,
            'ref_kmer': ref_seq[pos:pos+5] if pos < len(ref_seq) - 4 else 'N/A',
            'f5c_count': len(f5c_events),
            'pyfin_count': len(pyfin_events_at_pos),
            'f5c_means': [e['event_mean'] for e in f5c_events],
            'pyfin_event_indices': [e['event_idx'] for e in pyfin_events_at_pos],
        }

        # Get event means from pyfin
        if pyfin_events_at_pos:
            event_indices = [e['event_idx'] for e in pyfin_events_at_pos]
            event_means = pyfin_events['means'][event_indices]
            pos_detail['pyfin_means'] = event_means
        else:
            pos_detail['pyfin_means'] = []

        comparison['position_details'].append(pos_detail)

        # Categorize
        if len(f5c_events) > 0 and len(pyfin_events_at_pos) > 0:
            comparison['matched_positions'].append(pos)

            # Calculate mean difference
            f5c_mean = np.mean(pos_detail['f5c_means'])
            pyfin_mean = np.mean(pos_detail['pyfin_means'])
            comparison['event_mean_diff'].append(f5c_mean - pyfin_mean)
            comparison['event_count_diff'].append(len(f5c_events) - len(pyfin_events_at_pos))

        elif len(f5c_events) > 0:
            comparison['f5c_only_positions'].append(pos)
        else:
            comparison['pyfin_only_positions'].append(pos)

    # Calculate statistics
    n_f5c = len(f5c_alignments)
    n_pyfin = len(pyfin_align)
    n_matched_pos = len(comparison['matched_positions'])

    comparison['statistics'] = {
        'f5c_total_alignments': n_f5c,
        'pyfin_total_alignments': n_pyfin,
        'f5c_unique_positions': len(f5c_by_pos),
        'pyfin_unique_positions': len(pyfin_by_pos),
        'matched_positions': n_matched_pos,
        'f5c_only_positions': len(comparison['f5c_only_positions']),
        'pyfin_only_positions': len(comparison['pyfin_only_positions']),
    }

    if comparison['event_mean_diff']:
        comparison['statistics']['mean_diff_mean'] = np.mean(comparison['event_mean_diff'])
        comparison['statistics']['mean_diff_std'] = np.std(comparison['event_mean_diff'])
        comparison['statistics']['mean_diff_abs'] = np.mean(np.abs(comparison['event_mean_diff']))

    if comparison['event_count_diff']:
        comparison['statistics']['count_diff_mean'] = np.mean(comparison['event_count_diff'])
        comparison['statistics']['count_diff_std'] = np.std(comparison['event_count_diff'])

    return comparison


def print_detailed_comparison(f5c_alignments: List[Dict], comparison: Dict):
    """Print detailed comparison results."""
    stats = comparison['statistics']

    print("\nOverall Statistics:")
    print(f"  F5C total alignments: {stats['f5c_total_alignments']}")
    print(f"  PyFin total alignments: {stats['pyfin_total_alignments']}")
    print(f"  F5C unique positions: {stats['f5c_unique_positions']}")
    print(f"  PyFin unique positions: {stats['pyfin_unique_positions']}")
    print(f"  Matched positions: {stats['matched_positions']}")
    print(f"  F5C only positions: {stats['f5c_only_positions']}")
    print(f"  PyFin only positions: {stats['pyfin_only_positions']}")

    if 'mean_diff_mean' in stats:
        print("\n  Event mean difference (F5C - PyFin):")
        print(f"    Mean: {stats['mean_diff_mean']:.4f}")
        print(f"    Std: {stats['mean_diff_std']:.4f}")
        print(f"    Abs mean: {stats['mean_diff_abs']:.4f}")

    if 'count_diff_mean' in stats:
        print("\n  Event count difference per position (F5C - PyFin):")
        print(f"    Mean: {stats['count_diff_mean']:.2f}")
        print(f"    Std: {stats['count_diff_std']:.2f}")

    # Print first 50 matched positions
    print(f"\n{'='*80}")
    print("First 50 Matched Positions (Event-by-Event)")
    print(f"{'='*80}")
    print(f"{'Pos':>6} {'RefKmer':>10} {'F5C_N':>7} {'PyFin_N':>8} {'F5C_Mean':>10} {'PyFin_Mean':>11} {'Diff':>8}")
    print(f"{'-'*6} {'-'*10} {'-'*7} {'-'*8} {'-'*10} {'-'*11} {'-'*8}")

    matched_details = [d for d in comparison['position_details']
                      if d['position'] in comparison['matched_positions']]

    for i, detail in enumerate(matched_details[:50]):
        f5c_mean = np.mean(detail['f5c_means']) if len(detail['f5c_means']) > 0 else 0
        pyfin_mean = np.mean(detail['pyfin_means']) if len(detail['pyfin_means']) > 0 else 0
        diff = f5c_mean - pyfin_mean

        print(f"{detail['position']:6d} {detail['ref_kmer']:>10} "
              f"{detail['f5c_count']:7d} {detail['pyfin_count']:8d} "
              f"{f5c_mean:10.2f} {pyfin_mean:11.2f} {diff:8.4f}")

    # Print F5C-only positions
    if comparison['f5c_only_positions']:
        print(f"\n{'='*80}")
        print(f"F5C-Only Positions (first 20 of {len(comparison['f5c_only_positions'])})")
        print(f"{'='*80}")
        print(f"{'Pos':>6} {'RefKmer':>10} {'N_Events':>9} {'Event_Means':>30}")
        print(f"{'-'*6} {'-'*10} {'-'*9} {'-'*30}")

        for i, pos in enumerate(comparison['f5c_only_positions'][:20]):
            detail = comparison['position_details'][pos]  # Get detail
            means_str = ','.join(f"{m:.1f}" for m in detail['f5c_means'][:3])
            if len(detail['f5c_means']) > 3:
                means_str += '...'
            print(f"{pos:6d} {detail['ref_kmer']:>10} {detail['f5c_count']:9d} {means_str:>30}")

    # Print PyFin-only positions
    if comparison['pyfin_only_positions']:
        print(f"\n{'='*80}")
        print(f"PyFin-Only Positions (first 20 of {len(comparison['pyfin_only_positions'])})")
        print(f"{'='*80}")
        print(f"{'Pos':>6} {'RefKmer':>10} {'N_Events':>9} {'Event_Indices':>30}")
        print(f"{'-'*6} {'-'*10} {'-'*9} {'-'*30}")

        for i, pos in enumerate(comparison['pyfin_only_positions'][:20]):
            detail = comparison['position_details'][pos]
            indices_str = ','.join(str(idx) for idx in detail['pyfin_event_indices'][:3])
            if len(detail['pyfin_event_indices']) > 3:
                indices_str += '...'
            print(f"{pos:6d} {detail['ref_kmer']:>10} {detail['pyfin_count']:9d} {indices_str:>30}")


def create_text_visualization(f5c_alignments: List[Dict], comparison: Dict,
                              output_path: str):
    """Create a text-based visualization of the comparison."""
    with open(output_path, 'w') as f:
        f.write("# Event-by-Event Comparison Report\n")
        f.write("="*80 + "\n\n")

        stats = comparison['statistics']
        f.write("## Overall Statistics\n\n")
        f.write(f"- F5C total alignments: {stats['f5c_total_alignments']}\n")
        f.write(f"- PyFin total alignments: {stats['pyfin_total_alignments']}\n")
        f.write(f"- F5C unique positions: {stats['f5c_unique_positions']}\n")
        f.write(f"- PyFin unique positions: {stats['pyfin_unique_positions']}\n")
        f.write(f"- Matched positions: {stats['matched_positions']}\n")
        f.write(f"- F5C only positions: {stats['f5c_only_positions']}\n")
        f.write(f"- PyFin only positions: {stats['pyfin_only_positions']}\n\n")

        if 'mean_diff_mean' in stats:
            f.write("## Event Mean Difference (F5C - PyFin)\n\n")
            f.write(f"- Mean: {stats['mean_diff_mean']:.4f}\n")
            f.write(f"- Std: {stats['mean_diff_std']:.4f}\n")
            f.write(f"- Absolute mean: {stats['mean_diff_abs']:.4f}\n\n")

        # Position coverage map
        f.write("## Position Coverage Map\n\n")
        f.write("Legend:\n")
        f.write("  = Both F5C and PyFin have events at this position\n")
        f.write("  F Only F5C has events\n")
        f.write("  P Only PyFin has events\n")
        f.write("  . Neither has events\n\n")

        all_positions = sorted(comparison['statistics']['matched_positions'] +
                               comparison['f5c_only_positions'] +
                               comparison['pyfin_only_positions'])

        if all_positions:
            pos_range = (min(all_positions), max(all_positions))
            f.write(f"Position range: {pos_range[0]} - {pos_range[1]}\n\n")

            # Create 100-character wide map
            map_width = 100
            pos_per_char = max(1, (pos_range[1] - pos_range[0]) // map_width)

            f.write("Coverage Map:\n")
            for start in range(pos_range[0], pos_range[1], map_width * pos_per_char):
                end = min(start + map_width * pos_per_char, pos_range[1])
                line = []
                for p in range(start, end, pos_per_char):
                    if p in comparison['matched_positions']:
                        line.append('=')
                    elif p in comparison['f5c_only_positions']:
                        line.append('F')
                    elif p in comparison['pyfin_only_positions']:
                        line.append('P')
                    else:
                        line.append('.')
                f.write(f"{start:6d}: {''.join(line)}\n")

        # Detailed event table
        f.write("\n## Detailed Event Table (First 100 Matched Positions)\n\n")
        f.write(f"{'Pos':>6} {'RefKmer':>10} {'F5C_N':>7} {'PyFin_N':>8} "
               f"{'F5C_Mean':>10} {'PyFin_Mean':>11} {'Diff':>8}\n")
        f.write(f"{'-'*6} {'-'*10} {'-'*7} {'-'*8} {'-'*10} {'-'*11} {'-'*8}\n")

        matched_details = [d for d in comparison['position_details']
                          if d['position'] in comparison['matched_positions']]

        for detail in matched_details[:100]:
            f5c_mean = np.mean(detail['f5c_means']) if detail['f5c_means'] else 0
            pyfin_mean = np.mean(detail['pyfin_means']) if detail['pyfin_means'] else 0
            diff = f5c_mean - pyfin_mean

            f.write(f"{detail['position']:6d} {detail['ref_kmer']:>10} "
                   f"{detail['f5c_count']:7d} {detail['pyfin_count']:8d} "
                   f"{f5c_mean:10.2f} {pyfin_mean:11.2f} {diff:8.4f}\n")

    print(f"\nText visualization saved to: {output_path}")


def create_matplotlib_visualization(f5c_alignments: List[Dict], comparison: Dict,
                                    output_path: str):
    """Create matplotlib visualizations."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
    except ImportError:
        print("\nmatplotlib not available, skipping visualizations")
        print("Install with: pip install matplotlib")
        return

    stats = comparison['statistics']
    matched_details = [d for d in comparison['position_details']
                      if d['position'] in comparison['matched_positions']]

    if not matched_details:
        print("\nNo matched positions to visualize")
        return

    # Create figure with subplots
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.3, wspace=0.3)

    # 1. Event count comparison
    ax1 = fig.add_subplot(gs[0, 0])
    positions = [d['position'] for d in matched_details]
    f5c_counts = [d['f5c_count'] for d in matched_details]
    pyfin_counts = [d['pyfin_count'] for d in matched_details]

    ax1.scatter(f5c_counts, pyfin_counts, alpha=0.5, s=10)
    ax1.plot([0, max(max(f5c_counts), max(pyfin_counts))],
             [0, max(max(f5c_counts), max(pyfin_counts))],
             'r--', label='y=x')
    ax1.set_xlabel('F5C Event Count')
    ax1.set_ylabel('PyFin Event Count')
    ax1.set_title('Event Count per Position')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. Event mean comparison
    ax2 = fig.add_subplot(gs[0, 1])
    f5c_means = [np.mean(d['f5c_means']) for d in matched_details]
    pyfin_means = [np.mean(d['pyfin_means']) for d in matched_details]

    ax2.scatter(f5c_means, pyfin_means, alpha=0.5, s=10)
    min_mean = min(min(f5c_means), min(pyfin_means))
    max_mean = max(max(f5c_means), max(pyfin_means))
    ax2.plot([min_mean, max_mean], [min_mean, max_mean], 'r--', label='y=x')
    ax2.set_xlabel('F5C Event Mean')
    ax2.set_ylabel('PyFin Event Mean')
    ax2.set_title('Event Mean per Position')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. Position coverage along reference
    ax3 = fig.add_subplot(gs[1, :])

    all_positions = sorted(set(positions))
    f5c_coverage = [1 if p in [d['position'] for d in matched_details] or
                    p in comparison['f5c_only_positions'] else 0
                    for p in all_positions]
    pyfin_coverage = [1 if p in [d['position'] for d in matched_details] or
                      p in comparison['pyfin_only_positions'] else 0
                      for p in all_positions]

    ax3.plot(all_positions, f5c_coverage, label='F5C', linewidth=2)
    ax3.plot(all_positions, pyfin_coverage, label='PyFin', linewidth=2)
    ax3.set_xlabel('Reference Position')
    ax3.set_ylabel('Coverage (1=has events)')
    ax3.set_title('Position Coverage Along Reference')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_yticks([0, 1])
    ax3.set_yticklabels(['No', 'Yes'])

    # 4. Event count distribution
    ax4 = fig.add_subplot(gs[2, 0])
    ax4.hist(f5c_counts, bins=20, alpha=0.5, label='F5C', color='blue')
    ax4.hist(pyfin_counts, bins=20, alpha=0.5, label='PyFin', color='orange')
    ax4.set_xlabel('Event Count per Position')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Event Count Distribution')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # 5. Event mean difference
    ax5 = fig.add_subplot(gs[2, 1])
    mean_diffs = [f5c_means[i] - pyfin_means[i] for i in range(len(f5c_means))]
    ax5.hist(mean_diffs, bins=50, edgecolor='black', alpha=0.7)
    ax5.axvline(0, color='red', linestyle='--', linewidth=2, label='Zero difference')
    ax5.axvline(np.mean(mean_diffs), color='green', linestyle='--',
               linewidth=2, label=f'Mean: {np.mean(mean_diffs):.4f}')
    ax5.set_xlabel('Event Mean Difference (F5C - PyFin)')
    ax5.set_ylabel('Frequency')
    ax5.set_title('Event Mean Difference Distribution')
    ax5.legend()
    ax5.grid(True, alpha=0.3)

    plt.suptitle('F5C vs PyFin Event Alignment Comparison', fontsize=16, y=0.995)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nVisualization saved to: {output_path}")
    plt.close()


def main():
    """Main comparison function."""
    test_dir = Path(__file__).parent / "test_data"
    tsv_path = test_dir / "one_read.eventalign.tsv.gz"
    fasta_path = test_dir / "one_read.fa"
    pod5_path = test_dir / "one_read.pod5"
    fastq_path = test_dir / "one_read.fq"

    print("="*80)
    print("Detailed Event-by-Event Comparison: F5C vs PyFin")
    print("="*80)

    # Load F5C output
    print(f"\nLoading F5C output from: {tsv_path}")
    f5c_alignments = load_f5c_eventalign_tsv(str(tsv_path))
    print(f"  Loaded {len(f5c_alignments)} alignments")

    # Load reference
    print(f"\nLoading reference from: {fasta_path}")
    ref_name, ref_seq, ref_len = load_reference_from_fasta(str(fasta_path))
    print(f"  Reference: {ref_name}")
    print(f"  Length: {ref_len} bp")

    # Load signal and run pyfin
    print(f"\nLoading signal from: {pod5_path}")
    read_id, signal, sample_rate = load_signal_from_pod5(pod5_path)
    print(f"  Read ID: {read_id}")
    print(f"  Signal length: {len(signal)} samples")

    # Load read sequence from FASTQ if available
    if fastq_path.exists():
        print(f"\nLoading basecalled read sequence from: {fastq_path}")
        fastq_read_id, read_seq = load_read_sequence_from_fastq(str(fastq_path))
        print(f"  FASTQ read ID: {fastq_read_id}")
        print(f"  Read sequence length: {len(read_seq)} bp")

        # Verify read IDs match
        if fastq_read_id != read_id:
            print("  WARNING: Read ID mismatch!")
            print(f"    POD5 read_id: {read_id}")
            print(f"    FASTQ read_id: {fastq_read_id}")
    else:
        print("\n  WARNING: FASTQ not found, using reference sequence as placeholder")
        print("  For proper comparison, use basecalled read sequence from FASTQ/BAM")
        read_seq = ref_seq

    try:
        from fin._eventalign import run_eventalign, MODEL_RNA002

        print("\nRunning PyFin eventalign...")
        result = run_eventalign(
            read_ids=[read_id],
            read_seqs=[read_seq],
            ref_seqs=[ref_seq],
            ref_names=[ref_name],
            ref_lens=[ref_len],
            signals=[signal],
            sample_rates=[sample_rate],
            model_id=MODEL_RNA002,
        )

        pyfin_align = result["full"][0][0]

        if len(pyfin_align) == 0:
            print("\nPyFin alignment FAILED")
            print(f"  Status: {result['mapping'][0][0].get('status', 'unknown')}")
            print(f"  Events detected: {result['events'][0]['starts'].shape[0]}")
            print("\n  Cannot perform event-by-event comparison due to alignment failure.")
            print("  This is expected when using reference sequence as read_seq.")
            print("  For proper comparison, provide basecalled read sequence from FASTQ/BAM.")
            return

        print(f"  PyFin alignment SUCCESS: {len(pyfin_align)} alignments")

        # Export pyfin results in f5c TSV format
        output_dir = Path(__file__).parent
        tsv_output = output_dir / "pyfin_eventalign.tsv.gz"
        print("\nExporting PyFin results in f5c TSV format...")
        # Add read_id to result for export
        result['read_ids'] = [read_id]
        export_pyfin_to_f5c_tsv(result, ref_name, ref_seq, str(tsv_output),
                                 signal=signal, include_raw_samples=False)

        # Perform detailed comparison
        comparison = compare_event_by_event(f5c_alignments, result, ref_seq)

        # Print detailed results
        print_detailed_comparison(f5c_alignments, comparison)

        # Create visualizations
        output_dir = Path(__file__).parent
        text_output = output_dir / "detailed_comparison_report.txt"
        plot_output = output_dir / "comparison_visualization.png"

        create_text_visualization(f5c_alignments, comparison, str(text_output))
        create_matplotlib_visualization(f5c_alignments, comparison, str(plot_output))

        print("\n" + "="*80)
        print("Summary")
        print("="*80)

        stats = comparison['statistics']
        print(f"\nF5C: {stats['f5c_total_alignments']} alignments at {stats['f5c_unique_positions']} positions")
        print(f"PyFin: {stats['pyfin_total_alignments']} alignments at {stats['pyfin_unique_positions']} positions")
        print(f"Matched positions: {stats['matched_positions']}")

        if 'mean_diff_abs' in stats:
            print(f"Mean event difference: {stats['mean_diff_abs']:.4f} pA")

        if stats['matched_positions'] > 0:
            match_pct = 100 * stats['matched_positions'] / stats['f5c_unique_positions']
            print(f"Position overlap: {match_pct:.1f}%")

    except ImportError as e:
        print(f"\nERROR: {e}")
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
