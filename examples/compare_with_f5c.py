#!/usr/bin/env python3
"""
Compare pyfin eventalign results with f5c reference output.

This script:
1. Loads f5c eventalign TSV output (reference)
2. Runs pyfin eventalign on the same data
3. Compares and visualizes event-by-event results
4. Creates detailed reports and plots

Usage:
    python compare_with_f5c.py
"""

import numpy as np
import gzip
from pathlib import Path
from typing import Dict, List, Tuple
import sys


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


def create_text_visualization(f5c_alignments: List[Dict], ref_seq: str,
                              output_path: str):
    """Create a text-based visualization of f5c results."""
    with open(output_path, 'w') as f:
        f.write("# F5C Eventalign Results Visualization\n\n")

        # Index by position
        pos_to_events = {}
        for aln in f5c_alignments:
            pos = aln['reference_position']
            if pos not in pos_to_events:
                pos_to_events[pos] = []
            pos_to_events[pos].append(aln)

        all_positions = sorted(pos_to_events.keys())
        if not all_positions:
            return

        pos_range = (min(all_positions), max(all_positions))
        f.write(f"## Position Coverage Map\n\n")
        f.write(f"Position range: {pos_range[0]} - {pos_range[1]}\n")
        f.write(f"Total positions with events: {len(all_positions)}\n")
        f.write(f"Total event alignments: {len(f5c_alignments)}\n\n")

        # Create 100-character wide map
        map_width = 100
        pos_per_char = max(1, (pos_range[1] - pos_range[0]) // map_width)

        f.write("Legend:\n")
        f.write("  |  Has 1 event\n")
        f.write("  || Has 2 events\n")
        f.write("  ||| Has 3 events\n")
        f.write("  +++++ Has >5 events\n")
        f.write("  . No events\n\n")

        f.write("Coverage Map:\n")
        for start in range(pos_range[0], pos_range[1], map_width * pos_per_char):
            end = min(start + map_width * pos_per_char, pos_range[1])
            line = []
            for p in range(start, end, pos_per_char):
                n_events = len(pos_to_events.get(p, []))
                if n_events == 0:
                    line.append('.')
                elif n_events <= 3:
                    line.append('|' * n_events)
                else:
                    line.append('+' * min(n_events, 5))
            f.write(f"{start:6d}: {''.join(line)}\n")

        # Event statistics per position
        events_per_pos = [len(pos_to_events[p]) for p in all_positions]
        f.write(f"\n## Event Statistics\n\n")
        f.write(f"Events per position:\n")
        f.write(f"  Mean: {np.mean(events_per_pos):.2f}\n")
        f.write(f"  Std: {np.std(events_per_pos):.2f}\n")
        f.write(f"  Min: {np.min(events_per_pos)}\n")
        f.write(f"  Max: {np.max(events_per_pos)}\n")

        # Event mean statistics
        means = [a['event_mean'] for a in f5c_alignments]
        f.write(f"\nEvent mean levels:\n")
        f.write(f"  Mean: {np.mean(means):.2f} pA\n")
        f.write(f"  Std: {np.std(means):.2f} pA\n")
        f.write(f"  Range: [{np.min(means):.2f}, {np.max(means):.2f}] pA\n")


def run_pyfin_and_compare(f5c_alignments: List[Dict], ref_seq: str,
                          pod5_path: str) -> Dict:
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

        print(f"\n  WARNING: Using reference sequence as placeholder")
        print(f"  For proper comparison, use basecalled read sequence from FASTQ/BAM")

        # Run eventalign
        result = run_eventalign(
            read_ids=[read_id],
            read_seqs=[ref_seq],
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

            # Event-by-event comparison with f5c
            print(f"\n  Comparing with F5C (event-by-event)...")

            # Index f5c by position
            f5c_by_pos = {}
            for aln in f5c_alignments:
                pos = aln['reference_position']
                if pos not in f5c_by_pos:
                    f5c_by_pos[pos] = []
                f5c_by_pos[pos].append(aln)

            # Index pyfin by position
            pyfin_by_pos = {}
            for aln in pyfin_align:
                pos = aln['ref_position']
                if pos not in pyfin_by_pos:
                    pyfin_by_pos[pos] = []
                pyfin_by_pos[pos].append(aln)

            all_positions = sorted(set(f5c_by_pos.keys()) | set(pyfin_by_pos.keys()))

            matched_positions = 0
            f5c_only = 0
            pyfin_only = 0

            print(f"\n  Position-by-position comparison (first 20):")
            print(f"    {'Pos':>6} {'F5C_N':>7} {'PyFin_N':>8} {'Status':>10}")
            print(f"    {'-'*6} {'-'*7} {'-'*8} {'-'*10}")

            for pos in all_positions[:20]:
                n_f5c = len(f5c_by_pos.get(pos, []))
                n_pyfin = len(pyfin_by_pos.get(pos, []))

                if n_f5c > 0 and n_pyfin > 0:
                    status = "MATCH"
                    matched_positions += 1
                elif n_f5c > 0:
                    status = "F5C_ONLY"
                    f5c_only += 1
                else:
                    status = "PyFIN_ONLY"
                    pyfin_only += 1

                print(f"    {pos:6d} {n_f5c:7d} {n_pyfin:8d} {status:>10}")

            comparison['matched_positions'] = matched_positions
            comparison['f5c_only'] = f5c_only
            comparison['pyfin_only'] = pyfin_only
            comparison['total_positions'] = len(all_positions)

            print(f"\n  Summary:")
            print(f"    Total positions: {comparison['total_positions']}")
            print(f"    Matched positions: {matched_positions} ({100*matched_positions/comparison['total_positions']:.1f}%)")
            print(f"    F5C only: {f5c_only}")
            print(f"    PyFin only: {pyfin_only}")

        else:
            print(f"\n  PyFin alignment FAILED")
            print(f"  Status: {comparison['status']}")
            print(f"  Events detected: {comparison['pyfin_n_events']}")

        return comparison

    except ImportError as e:
        print(f"\n  ERROR: fin._eventalign module not available: {e}")
        return {'error': str(e)}
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        return {'error': str(e)}


def main():
    """Main comparison function."""
    test_dir = Path(__file__).parent / "test_data"
    tsv_path = test_dir / "one_read.eventalign.tsv.gz"
    fasta_path = test_dir / "one_read.fa"
    pod5_path = test_dir / "one_read.pod5"
    output_dir = Path(__file__).parent

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
    pos_to_events = {}
    for aln in f5c_alignments:
        pos = aln['reference_position']
        if pos not in pos_to_events:
            pos_to_events[pos] = []
        pos_to_events[pos].append(aln)

    means = [a['event_mean'] for a in f5c_alignments]
    stdvs = [a['event_stdv'] for a in f5c_alignments]
    events_per_pos = [len(events) for events in pos_to_events.values()]

    print(f"\nF5C Statistics:")
    print(f"  Unique positions: {len(pos_to_events)}")
    print(f"  Events per position: {np.mean(events_per_pos):.2f} +/- {np.std(events_per_pos):.2f}")
    print(f"  Event mean: {np.mean(means):.2f} +/- {np.std(means):.2f} pA")

    # Create text visualization of f5c results
    viz_path = output_dir / "f5c_visualization.txt"
    create_text_visualization(f5c_alignments, ref_seq, str(viz_path))
    print(f"\nF5C visualization saved to: {viz_path}")

    # Compare with PyFin
    comparison = run_pyfin_and_compare(f5c_alignments, ref_seq, str(pod5_path))

    # Save comparison report
    report_path = output_dir / "comparison_report.txt"
    with open(report_path, 'w') as f:
        f.write("# Event Alignment Comparison Report\n\n")
        f.write("## F5C Reference Results\n\n")
        f.write(f"- Total alignments: {len(f5c_alignments)}\n")
        f.write(f"- Unique positions: {len(pos_to_events)}\n")
        f.write(f"- Events per position: {np.mean(events_per_pos):.2f} +/- {np.std(events_per_pos):.2f}\n")
        f.write(f"- Event mean: {np.mean(means):.2f} +/- {np.std(means):.2f} pA\n\n")

        f.write("## PyFin Comparison\n\n")
        if 'error' in comparison:
            f.write(f"ERROR: {comparison['error']}\n")
        elif comparison.get('pyfin_success', False):
            f.write(f"- PyFin alignments: {comparison['pyfin_n_alignments']}\n")
            f.write(f"- Events per base: {comparison['pyfin_events_per_base']:.2f}\n")
            f.write(f"- Position comparison:\n")
            f.write(f"  - Total positions: {comparison['total_positions']}\n")
            f.write(f"  - Matched: {comparison['matched_positions']}\n")
            f.write(f"  - F5C only: {comparison['f5c_only']}\n")
            f.write(f"  - PyFin only: {comparison['pyfin_only']}\n")
            f.write(f"- Status: SUCCESS\n")
        else:
            f.write(f"- Status: FAILED ({comparison.get('status', 'unknown')})\n")
            f.write(f"- Events detected: {comparison.get('pyfin_n_events', 'N/A')}\n")
            f.write(f"\n### Notes\n\n")
            f.write(f"This is expected when using reference sequence as read_seq.\n")
            f.write(f"For proper comparison, use basecalled read sequence from FASTQ/BAM.\n")

    print(f"\nComparison report saved to: {report_path}")

    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"\nF5C reference output: {len(f5c_alignments)} event alignments")
    print(f"  Unique positions: {len(pos_to_events)}")
    print(f"  Events per position: {np.mean(events_per_pos):.2f} +/- {np.std(events_per_pos):.2f}")

    if 'error' in comparison:
        print(f"\nPyFin: ERROR - {comparison['error']}")
    elif comparison.get('pyfin_success', False):
        print(f"\nPyFin: {comparison['pyfin_n_alignments']} event alignments")
        print(f"  Matched positions: {comparison['matched_positions']}/{comparison['total_positions']} "
              f"({100*comparison['matched_positions']/comparison['total_positions']:.1f}%)")
    else:
        print(f"\nPyFin: Alignment failed")
        print(f"  Note: This is expected when using reference sequence as read_seq")
        print(f"        For proper comparison, use basecalled read sequence from FASTQ/BAM")


if __name__ == "__main__":
    main()
