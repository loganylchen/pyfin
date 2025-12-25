#!/usr/bin/env python3
"""
Compare pyfin eventalign results with f5c reference output.

This script:
1. Loads f5c eventalign TSV output (reference)
2. Runs pyfin eventalign on the same data
3. Compares and visualizes the results

Usage:
    python compare_with_f5c.py
"""

import numpy as np
import gzip
from pathlib import Path
from typing import Dict, List, Tuple


def load_f5c_eventalign_tsv(tsv_path: str) -> List[Dict]:
    """
    Load f5c eventalign TSV output.

    Format (tab-separated):
    1. reference_name
    2. reference_position (0-based)
    3. reference_kmer
    4. read_id
    5. strand (t/f)
    6. event_idx
    7. event_mean
    8. event_stdv
    9. duration (sum of weights)
    10. model_kmer
    11. model_mean
    12. model_stdv
    13. scaled_mean
    14. start_idx (in raw signal)
    15. end_idx (in raw signal)
    16. raw_samples (comma-separated)
    """
    alignments = []

    # Handle .gz files
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

            # Parse raw samples (comma-separated)
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


def analyze_f5c_output(alignments: List[Dict]) -> Dict:
    """Analyze f5c eventalign output and extract statistics."""
    if not alignments:
        return {}

    # Group by reference position
    pos_to_events: Dict[int, List[Dict]] = {}
    for aln in alignments:
        pos = aln['reference_position']
        if pos not in pos_to_events:
            pos_to_events[pos] = []
        pos_to_events[pos].append(aln)

    # Calculate statistics
    means = [a['event_mean'] for a in alignments]
    stdvs = [a['event_stdv'] for a in alignments]
    events_per_pos = [len(events) for events in pos_to_events.values()]

    stats = {
        'total_alignments': len(alignments),
        'unique_positions': len(pos_to_events),
        'reference_name': alignments[0]['reference_name'],
        'read_id': alignments[0]['read_id'],
        'events_per_position_mean': np.mean(events_per_pos),
        'events_per_position_std': np.std(events_per_pos),
        'event_mean_mean': np.mean(means),
        'event_mean_std': np.std(means),
        'event_mean_min': np.min(means),
        'event_mean_max': np.max(means),
        'event_stdv_mean': np.mean(stdvs),
        'event_stdv_std': np.std(stdvs),
        'position_range': (
            min(a['reference_position'] for a in alignments),
            max(a['reference_position'] for a in alignments)
        ),
        'positions': set(pos_to_events.keys()),
    }

    return stats


def print_f5c_summary(stats: Dict):
    """Print summary of f5c eventalign output."""
    print("=" * 70)
    print("F5C Eventalign Output Summary")
    print("=" * 70)

    print(f"\nReference: {stats['reference_name']}")
    print(f"Read ID: {stats['read_id']}")
    print(f"\nTotal alignments: {stats['total_alignments']}")
    print(f"Unique reference positions: {stats['unique_positions']}")
    print(f"Position range: {stats['position_range'][0]} - {stats['position_range'][1]}")

    print(f"\nEvent Statistics:")
    print(f"  Mean level: {stats['event_mean_mean']:.2f} +/- {stats['event_mean_std']:.2f}")
    print(f"    Range: [{stats['event_mean_min']:.2f}, {stats['event_mean_max']:.2f}]")
    print(f"  Stdv level: {stats['event_stdv_mean']:.4f} +/- {stats['event_stdv_std']:.4f}")

    print(f"\nEvents per position:")
    print(f"  Mean: {stats['events_per_position_mean']:.2f}")
    print(f"  Std: {stats['events_per_position_std']:.2f}")


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
        # Generate synthetic signal for testing
        print("Generating synthetic signal...")
        np.random.seed(42)
        n_samples = 100000
        signal = np.random.randn(n_samples).astype(np.float32) * 10 + 120
        sample_rate = 4000.0
        return "synthetic_read", signal, sample_rate


def run_pyfin_and_compare(f5c_stats: Dict, ref_seq: str, pod5_path: str) -> Dict:
    """Run pyfin eventalign and compare with f5c results."""
    print("\n" + "=" * 70)
    print("Running PyFin Eventalign")
    print("=" * 70)

    try:
        from fin._eventalign import run_eventalign, MODEL_RNA002

        # Load signal from POD5
        read_id, signal, sample_rate = load_signal_from_pod5(pod5_path)
        print(f"  Read ID: {read_id}")
        print(f"  Signal length: {len(signal)} samples")
        print(f"  Sample rate: {sample_rate} Hz")

        # IMPORTANT: For proper comparison, we need the basecalled read sequence
        # that matches the signal. The reference sequence is NOT the same as the read.
        # For this comparison, we'll use the reference as read_seq, but alignment may fail.

        print(f"\n  WARNING: Using reference sequence as placeholder")
        print(f"  For proper comparison, use actual basecalled read sequence from FASTQ/BAM")
        print(f"  See: load_from_fastq_pod5.py or load_from_bam_pod5.py")

        # Run eventalign
        result = run_eventalign(
            read_ids=[read_id],
            read_seqs=[ref_seq],  # Using reference as placeholder
            ref_seqs=[ref_seq],
            ref_names=["SIRV101"],
            ref_lens=[len(ref_seq)],
            signals=[signal],
            sample_rates=[sample_rate],
            model_id=MODEL_RNA002,
        )

        # Extract pyfin results
        pyfin_align = result["full"][0][0]
        pyfin_mapping = result["mapping"][0][0]
        pyfin_events = result["events"][0]

        comparison = {
            'pyfin_success': len(pyfin_align) > 0,
            'pyfin_n_alignments': len(pyfin_align),
            'pyfin_events_per_base': pyfin_mapping.get('events_per_base', 0),
            'pyfin_n_events': pyfin_events['starts'].shape[0],
            'f5c_n_alignments': f5c_stats['total_alignments'],
            'f5c_unique_positions': f5c_stats['unique_positions'],
            'status': pyfin_mapping.get('status', 'unknown'),
        }

        if len(pyfin_align) > 0:
            print(f"\n  PyFin alignment SUCCESS: {len(pyfin_align)} alignments")
            print(f"  Events per base: {comparison['pyfin_events_per_base']:.2f}")

            # Show first few pyfin alignments
            print(f"\n  First 10 PyFin alignments:")
            print(f"    {'RefPos':>8} {'EventIdx':>10} {'RefKmer':>10} {'State':>6}")
            print(f"    {'-'*8} {'-'*10} {'-'*10} {'-'*6}")
            for i in range(min(10, len(pyfin_align))):
                ea = pyfin_align[i]
                print(f"    {ea['ref_position']:8d} {ea['event_idx']:10d} "
                      f"{ea['ref_kmer']:>10} {ea['hmm_state']:>6}")

            # Compare position coverage
            pyfin_positions = set(a['ref_position'] for a in pyfin_align)

            comparison['pyfin_unique_positions'] = len(pyfin_positions)
            comparison['position_overlap'] = len(f5c_stats['positions'] & pyfin_positions)

            print(f"\n  Position Coverage Comparison:")
            print(f"    F5C unique positions: {comparison['f5c_unique_positions']}")
            print(f"    PyFin unique positions: {comparison['pyfin_unique_positions']}")
            print(f"    Overlap: {comparison['position_overlap']} "
                  f"({100*comparison['position_overlap']/max(comparison['f5c_unique_positions'], 1):.1f}%)")

            # Event statistics comparison
            pyfin_means = pyfin_events['means']
            print(f"\n  Event Statistics Comparison:")
            print(f"    F5C event mean: {f5c_stats['event_mean_mean']:.2f} +/- {f5c_stats['event_mean_std']:.2f}")
            print(f"    PyFin event mean: {np.mean(pyfin_means):.2f} +/- {np.std(pyfin_means):.2f}")

        else:
            print(f"\n  PyFin alignment FAILED")
            print(f"  Status: {comparison['status']}")
            print(f"  Events detected: {comparison['pyfin_n_events']}")

            # Show diagnostic info
            print(f"\n  Diagnostics:")
            print(f"    F5C alignments: {comparison['f5c_n_alignments']}")
            print(f"    F5C unique positions: {comparison['f5c_unique_positions']}")
            print(f"    PyFin events detected: {comparison['pyfin_n_events']}")

            # Calculate expected events per kmer
            n_kmers = len(ref_seq) - 5 + 1
            expected_epk = comparison['pyfin_n_events'] / n_kmers
            print(f"    Events/k-mer ratio: {expected_epk:.2f} (typical: 2-4 for RNA)")

        return comparison

    except ImportError as e:
        print(f"\n  ERROR: fin._eventalign module not available: {e}")
        return {'error': str(e)}
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        return {'error': str(e)}


def save_comparison_report(f5c_stats: Dict, comparison: Dict, output_path: str):
    """Save comparison report to file."""
    with open(output_path, 'w') as f:
        f.write("# Event Alignment Comparison Report\n\n")
        f.write("## F5C Reference Results\n\n")
        f.write(f"- Total alignments: {f5c_stats['total_alignments']}\n")
        f.write(f"- Unique positions: {f5c_stats['unique_positions']}\n")
        f.write(f"- Position range: {f5c_stats['position_range'][0]} - {f5c_stats['position_range'][1]}\n\n")

        f.write("### Event Statistics\n\n")
        f.write(f"- Event mean: {f5c_stats['event_mean_mean']:.2f} +/- {f5c_stats['event_mean_std']:.2f}\n")
        f.write(f"- Event stdv: {f5c_stats['event_stdv_mean']:.4f} +/- {f5c_stats['event_stdv_std']:.4f}\n")
        f.write(f"- Events per position: {f5c_stats['events_per_position_mean']:.2f} +/- {f5c_stats['events_per_position_std']:.2f}\n\n")

        f.write("## PyFin Comparison\n\n")
        if 'error' in comparison:
            f.write(f"ERROR: {comparison['error']}\n")
        elif comparison.get('pyfin_success', False):
            f.write(f"- PyFin alignments: {comparison['pyfin_n_alignments']}\n")
            f.write(f"- Events per base: {comparison['pyfin_events_per_base']:.2f}\n")
            f.write(f"- Position overlap: {comparison['position_overlap']}/{comparison['f5c_unique_positions']} "
                    f"({100*comparison['position_overlap']/comparison['f5c_unique_positions']:.1f}%)\n")
            f.write(f"- Status: SUCCESS\n")
        else:
            f.write(f"- Status: FAILED ({comparison.get('status', 'unknown')})\n")
            f.write(f"- Events detected: {comparison.get('pyfin_n_events', 'N/A')}\n")
            f.write(f"\n### Notes\n\n")
            f.write(f"This is expected when using reference sequence as read_seq.\n")
            f.write(f"For proper comparison:\n")
            f.write(f"1. Basecall POD5 with Guppy/Dorado to get FASTQ\n")
            f.write(f"2. Or use BAM file with aligned reads\n")
            f.write(f"3. See load_from_fastq_pod5.py or load_from_bam_pod5.py\n")

    print(f"\nComparison report saved to: {output_path}")


def main():
    """Main comparison function."""
    test_dir = Path(__file__).parent / "test_data"
    tsv_path = test_dir / "one_read.eventalign.tsv.gz"
    fasta_path = test_dir / "one_read.fa"
    pod5_path = test_dir / "one_read.pod5"
    output_path = Path(__file__).parent / "comparison_report.txt"

    print("=" * 70)
    print("F5C vs PyFin Eventalign Comparison")
    print("=" * 70)

    # Load F5C output
    print(f"\nLoading F5C output from: {tsv_path}")
    f5c_alignments = load_f5c_eventalign_tsv(str(tsv_path))
    print(f"  Loaded {len(f5c_alignments)} alignments")

    # Show first few f5c alignments
    print(f"\n  First 10 F5C alignments:")
    print(f"    {'Pos':>6} {'EventIdx':>8} {'RefKmer':>10} {'EventMean':>10} {'ModelMean':>10}")
    print(f"    {'-'*6} {'-'*8} {'-'*10} {'-'*10} {'-'*10}")
    for i, aln in enumerate(f5c_alignments[:10]):
        print(f"    {aln['reference_position']:6d} {aln['event_idx']:8d} {aln['reference_kmer']:>10} "
              f"{aln['event_mean']:10.2f} {aln['model_mean']:10.2f}")

    # Load reference
    print(f"\nLoading reference from: {fasta_path}")
    ref_name, ref_seq, ref_len = load_reference_from_fasta(str(fasta_path))
    print(f"  Reference: {ref_name}")
    print(f"  Length: {ref_len} bp")

    # Analyze F5C output
    f5c_stats = analyze_f5c_output(f5c_alignments)
    print_f5c_summary(f5c_stats)

    # Compare with PyFin
    comparison = run_pyfin_and_compare(f5c_stats, ref_seq, str(pod5_path))

    # Save report
    save_comparison_report(f5c_stats, comparison, str(output_path))

    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"\nF5C reference output: {len(f5c_alignments)} event alignments")
    print(f"  Position range: {f5c_stats['position_range'][0]} - {f5c_stats['position_range'][1]}")
    print(f"  Events per position: {f5c_stats['events_per_position_mean']:.2f} +/- {f5c_stats['events_per_position_std']:.2f}")

    if 'error' in comparison:
        print(f"\nPyFin: ERROR - {comparison['error']}")
    elif comparison.get('pyfin_success', False):
        print(f"\nPyFin: {comparison['pyfin_n_alignments']} event alignments")
        print(f"  Position overlap: {comparison['position_overlap']}/{comparison['f5c_unique_positions']} "
              f"({100*comparison['position_overlap']/comparison['f5c_unique_positions']:.1f}%)")
        print(f"\n  Status: Results comparable!" if comparison['position_overlap'] > 0 else "  Status: No position overlap")
    else:
        print(f"\nPyFin: Alignment failed")
        print(f"  Status: {comparison['status']}")
        print(f"  Events detected: {comparison.get('pyfin_n_events', 'N/A')}")
        print(f"\n  Note: This is expected when using reference sequence as read_seq")
        print(f"        For proper comparison, use basecalled read sequence from FASTQ/BAM")


if __name__ == "__main__":
    main()
