#!/usr/bin/env python3
"""
Example: Load data from FASTQ (sequences) and POD5 (signals) for event alignment.

This script demonstrates the typical workflow where:
1. Basecalled read sequences are loaded from a FASTQ file
2. Raw signals are loaded from a POD5 file
3. They are paired by read ID for event alignment

FASTQ files are typically produced by basecallers like Guppy, Dorado, or Bonafider.
"""

import numpy as np
from pathlib import Path
from typing import List, Tuple

try:
    from fin.io.io_fastq import FastqReader
    from fin.io.io_pod5 import Pod5Reader
    from fin._eventalign import run_eventalign, MODEL_RNA002
except ImportError as e:
    print(f"Import error: {e}")
    print("Please ensure pyfin is installed correctly")
    exit(1)


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


def pair_reads_with_signals(
    fastq_path: str, pod5_path: str, max_reads: int = None
) -> Tuple[List[str], List[str], List[np.ndarray], List[float]]:
    """
    Pair basecalled sequences from FASTQ with signals from POD5.

    Args:
        fastq_path: Path to FASTQ file with basecalled reads
        pod5_path: Path to POD5 file with raw signals
        max_reads: Maximum number of reads to load (None = all)

    Returns:
        Tuple of (read_ids, read_seqs, signals, sample_rates)
    """
    read_ids = []
    read_seqs = []
    signals = []
    sample_rates = []

    # Build a mapping of read_id to signal from POD5
    print(f"Loading signals from POD5: {pod5_path}")
    with Pod5Reader(pod5_path) as pod5_reader:
        pod5_read_ids = pod5_reader.read_ids
        print(f"  Found {len(pod5_read_ids)} reads in POD5")

        # Pre-load signals into a dictionary for efficient lookup
        signal_map = {}
        for read_id in pod5_read_ids:
            read = pod5_reader.get_read(read_id)
            if read is not None:
                signal_map[str(read.read_id)] = {
                    "signal": read.signal_pa.astype(np.float32),
                    "sample_rate": float(read.run_info.sample_rate),
                }

    # Read basecalled sequences from FASTQ and pair with signals
    print(f"Loading sequences from FASTQ: {fastq_path}")
    with FastqReader(fastq_path) as fastq_reader:
        n_loaded = 0
        n_missing = 0
        n_skipped = 0

        for i, (read_id, sequence, quality) in enumerate(fastq_reader.iter_reads()):
            if max_reads is not None and i >= max_reads:
                break

            # Skip sequences with non-standard bases
            if any(base not in "ACGTU" for base in sequence.upper()):
                n_skipped += 1
                continue

            # Pair with signal from POD5
            if read_id in signal_map:
                read_ids.append(read_id)
                read_seqs.append(sequence.upper().replace("U", "T"))
                signals.append(signal_map[read_id]["signal"])
                sample_rates.append(signal_map[read_id]["sample_rate"])
                n_loaded += 1
            else:
                n_missing += 1

        print(f"  Loaded {n_loaded} paired reads")
        if n_missing > 0:
            print(f"  Warning: {n_missing} reads in FASTQ not found in POD5")
        if n_skipped > 0:
            print(f"  Warning: {n_skipped} reads skipped due to non-standard bases")

    return read_ids, read_seqs, signals, sample_rates


def main():
    """Main function demonstrating FASTQ+POD5 workflow."""
    print("=" * 70)
    print("Event Alignment from FASTQ (sequences) + POD5 (signals)")
    print("=" * 70)

    # Test data paths
    test_dir = Path(__file__).parent / "test_data"
    fastq_path = test_dir / "basecalled.fastq"
    pod5_path = test_dir / "one_read.pod5"
    fasta_path = test_dir / "one_read.fa"

    # Check if files exist
    if not fastq_path.exists():
        print(f"\nFASTQ file not found: {fastq_path}")
        print("\nTo create a FASTQ file with basecalled reads:")
        print("1. Basecall POD5 with Guppy:")
        print("   guppy_basecaller -i /path/to/pod5_dir -s /path/to/output")
        print("                      --flowcell FLO-MIN106 --kit SQK-RNA002")
        print("                      --cpu_threads_per_caller 8")
        print("2. Or use Dorado (newer):")
        print("   dorado basecaller DNA_MODEL /path/to/pod5_dir > basecalled.fastq")
        return

    if not fasta_path.exists():
        print(f"\nReference FASTA not found: {fasta_path}")
        return

    # Load reference
    print(f"\nLoading reference from: {fasta_path}")
    ref_name, ref_seq, ref_len = load_reference_from_fasta(fasta_path)
    print(f"  Reference: {ref_name}")
    print(f"  Length: {ref_len} bp")

    # Load reads (sequences from FASTQ, signals from POD5)
    print("\n" + "-" * 70)
    read_ids, read_seqs, signals, sample_rates = pair_reads_with_signals(
        str(fastq_path), str(pod5_path), max_reads=10
    )

    if len(read_ids) == 0:
        print("\nNo reads loaded. Check that:")
        print("1. FASTQ file contains basecalled sequences")
        print("2. POD5 file contains the same read IDs")
        print("3. Read IDs match between FASTQ and POD5")
        return

    print(f"\nLoaded {len(read_ids)} reads with paired sequences and signals")

    # Show summary of first few reads
    print("\nRead summary:")
    for i in range(min(3, len(read_ids))):
        print(f"  {i+1}. {read_ids[i]}")
        print(f"     Sequence length: {len(read_seqs[i])} bp")
        print(f"     Signal length: {len(signals[i])} samples")
        print(f"     Sample rate: {sample_rates[i]} Hz")

    # Run event alignment
    print("\n" + "=" * 70)
    print("Running Event Alignment")
    print("=" * 70)

    result = run_eventalign(
        read_ids=read_ids,
        read_seqs=read_seqs,
        ref_seqs=[ref_seq],
        ref_names=[ref_name],
        ref_lens=[ref_len],
        signals=signals,
        sample_rates=sample_rates,
        model_id=MODEL_RNA002,
    )

    # Display results
    print("\n" + "=" * 70)
    print("Results")
    print("=" * 70)

    print(f"\nSummary:")
    summary = result["summary"]
    print(f"  Reads processed: {summary['num_reads']}")
    print(f"  References: {summary['num_refs']}")

    # Show alignment results for each read
    for i, read_id in enumerate(read_ids):
        print(f"\n{read_id}:")
        full_align = result["full"][i][0]
        mapping = result["mapping"][i][0]

        if len(full_align) > 0:
            print(f"  Alignment SUCCESS: {len(full_align)} aligned events")
            print(f"  Events per base: {mapping['events_per_base']:.2f}")

            # Show first few aligned events
            print(f"\n  First 5 aligned events:")
            print(f"    {'RefPos':>8} {'EventIdx':>10} {'RefKmer':>10} {'State':>6}")
            print(f"    {'-'*8} {'-'*10} {'-'*10} {'-'*6}")
            for k in range(min(5, len(full_align))):
                ea = full_align[k]
                print(f"    {ea['ref_position']:8d} {ea['event_idx']:10d} {ea['ref_kmer']:>10} {ea['hmm_state']:>6}")
        else:
            print(f"  Alignment FAILED")
            print(f"  Status: {mapping.get('status', 'unknown')}")


if __name__ == "__main__":
    main()
