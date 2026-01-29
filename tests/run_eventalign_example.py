#!/usr/bin/env python3
"""
Example usage of the run_eventalign function.

This demonstrates how to use run_eventalign to perform pair-wise
event alignment between reads and references.
"""

import numpy as np
from pathlib import Path
from fin._eventalign import run_eventalign, MODEL_RNA002

# For loading signal from POD5 files
try:
    from fin.io.io_pod5 import Pod5Reader

    POD5_AVAILABLE = True
except ImportError:
    POD5_AVAILABLE = False
    print("Warning: pod5 not available, will use synthetic signal")


def load_reference_from_fasta(fasta_path: str) -> tuple:
    """
    Load reference sequence from FASTA file.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Tuple of (ref_name, ref_sequence, ref_length)
    """
    with open(fasta_path, "r") as f:
        lines = f.readlines()

    # Parse FASTA
    ref_name = None
    ref_sequence = []
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            ref_name = line[1:].split()[0]  # Get first part of header
        else:
            ref_sequence.append(line)

    ref_sequence = "".join(ref_sequence).upper()
    ref_length = len(ref_sequence)

    return ref_name, ref_sequence, ref_length


def load_signal_from_pod5(pod5_path: str, read_id: str = None) -> tuple:
    """
    Load signal data from POD5 file.

    Args:
        pod5_path: Path to POD5 file
        read_id: Specific read ID to load (None = first read)

    Returns:
        Tuple of (read_id, signal, sample_rate)
    """
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
            # Get first read
            read_id = reader.read_ids[0]

        read = reader.get_read(read_id)
        if read is None:
            raise ValueError(f"Read {read_id} not found in POD5 file")

        signal = read.signal_pa.astype(np.float32)
        sample_rate = float(read.run_info.sample_rate)

        return str(read.read_id), signal, sample_rate


def example_single_read_single_ref():
    """Example: Single read against single reference"""
    print("=" * 70)
    print("Example 1: Single Read vs Single Reference")
    print("=" * 70)

    # Paths to test data
    test_dir = Path(__file__).parent
    fasta_path = test_dir / "test_data" / "one_read.fa"
    pod5_path = test_dir / "test_data" / "one_read.pod5"

    # Check if files exist
    if not fasta_path.exists():
        print(f"FASTA file not found: {fasta_path}")
        print("Please ensure test data is available")
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

    # IMPORTANT: We need the basecalled read sequence that corresponds to this signal!
    # POD5 files only contain raw signal, not sequences.
    #
    # For this example, we're using the reference as a placeholder, but in practice:
    # 1. Basecall the POD5 with Guppy/Dorado to get FASTQ
    # 2. Or use a BAM file with aligned reads
    # See: load_from_fastq_pod5.py or load_from_bam_pod5.py for proper workflow
    print(f"\n  WARNING: Using reference sequence as placeholder")
    print(f"  For production use, provide actual basecalled read sequence!")
    read_seq = ref_seq  # Use first 500bp as read sequence

    print(f"\nRunning eventalign...")
    print(f"  Model: RNA002 (k=5)")
    print(f"  Read length: {len(read_seq)} bp")

    # Run eventalign
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

    # Display results
    print("\n" + "=" * 70)
    print("Results:")
    print("=" * 70)

    print(f"\nSummary:")
    summary = result["summary"]
    print(f"  Reads processed: {summary['num_reads']}")
    print(f"  References: {summary['num_refs']}")

    # Scalings
    print(f"\nScalings:")
    for i, sc in enumerate(result["scalings"]):
        print(f"  Read {i}: scale={sc['scale']:.4f}, shift={sc['shift']:.4f}, var={sc['var']:.4f}")

    # Events
    print(f"\nDetected Events:")
    events = result["events"][0]
    n_ev = len(events["starts"])
    print(f"  Total events: {n_ev}")
    print(f"\n  ALL events:")
    print(f"    {'Start':>10} {'Length':>10} {'Mean':>10} {'Stdv':>10}")
    print(f"    {'-'*10} {'-'*10} {'-'*10} {'-'*10}")
    for i in range(n_ev):
        print(
            f"    {events['starts'][i]:10d} "
            f"{events['lengths'][i]:10.1f} "
            f"{events['means'][i]:10.2f} "
            f"{events['stdvs'][i]:10.4f}"
        )

    # Full alignment results
    print(f"\nAlignment Results:")
    for i in range(len(result["full"])):  # For each read
        for j in range(len(result["full"][i])):  # For each reference
            full_align = result["full"][i][j]
            mapping = result["mapping"][i][j]

            if len(full_align) > 0:
                print(f"  Read {i} vs Ref {j}: {len(full_align)} aligned events")
                print(f"\n  First 5 aligned events:")
                print(
                    f"    {'RefPos':>8} {'EventIdx':>10} {'RefKmer':>10} {'ModelKmer':>10} "
                    f"{'State':>6}"
                )
                print(f"    {'-'*8} {'-'*10} {'-'*10} {'-'*10} {'-'*6}")
                for k in range(min(5, len(full_align))):
                    ea = full_align[k]
                    print(
                        f"    {ea['ref_position']:8d} "
                        f"{ea['event_idx']:10d} "
                        f"{ea['ref_kmer']:>10} "
                        f"{ea['model_kmer']:>10} "
                        f"{ea['hmm_state']:>6}"
                    )

                # Mapping info
                start = mapping["start"]
                stop = mapping["stop"]
                print(f"\n  Base-to-Event Mapping:")
                print(f"    Events per base: {mapping['events_per_base']:.2f}")
                print(f"    First 10 bases: start=[{', '.join(map(str, start[:10]))}]")
                print(f"                    stop=[{', '.join(map(str, stop[:10]))}]")

                # Diagnostic info
                print(f"\n  Alignment Statistics:")
                print(f"    Status: {mapping.get('status', 'unknown')}")
                print(f"    Aligned pairs: {mapping.get('n_aligned_pairs', 'N/A')}")
                print(f"    Event alignments: {mapping.get('n_event_alignment', 'N/A')}")
            else:
                print(f"  Read {i} vs Ref {j}: No alignment")
                # Show diagnostic info
                mapping = result["mapping"][i][j]
                print(f"    Status: {mapping.get('status', 'unknown')}")
                print(f"    Events detected: {mapping.get('n_events', 'N/A')}")
                print(f"    K-mers in reference: {mapping.get('n_kmers', 'N/A')}")
                print(f"    Reference length: {mapping.get('ref_len', 'N/A')} bp")
                print(f"    Read length: {mapping.get('read_len', 'N/A')} bp")


def example_single_read_multi_ref():
    """Example: Single read against multiple references"""
    print("\n" + "=" * 70)
    print("Example 2: Single Read vs Multiple References")
    print("=" * 70)

    test_dir = Path(__file__).parent
    fasta_path = test_dir / "one_read.fa"
    pod5_path = test_dir / "one_read.pod5"

    if not fasta_path.exists():
        print(f"FASTA file not found: {fasta_path}")
        return

    # Load reference
    ref_name, ref_seq, ref_len = load_reference_from_fasta(fasta_path)

    # Load signal
    read_id, signal, sample_rate = load_signal_from_pod5(pod5_path)

    # Create a shorter read sequence
    read_seq = ref_seq[:500]

    # Create multiple references (e.g., different regions)
    ref_seqs = [
        ref_seq[:1000],  # Reference 1: first 1000bp
        ref_seq[500:1500],  # Reference 2: 500-1500bp
        ref_seq[1000:2000],  # Reference 3: 1000-2000bp
    ]
    ref_names = [f"{ref_name}_region1", f"{ref_name}_region2", f"{ref_name}_region3"]
    ref_lens = [len(r) for r in ref_seqs]

    print(f"\nRunning pair-wise alignment:")
    print(f"  1 read vs {len(ref_seqs)} references")

    result = run_eventalign(
        read_ids=[read_id],
        read_seqs=[read_seq],
        ref_seqs=ref_seqs,
        ref_names=ref_names,
        ref_lens=ref_lens,
        signals=[signal],
        sample_rates=[sample_rate],
        model_id=MODEL_RNA002,
    )

    print(f"\nResults:")
    for j, ref_name in enumerate(ref_names):
        full_align = result["full"][0][j]
        mapping = result["mapping"][0][j]
        if len(full_align) > 0:
            print(
                f"  {ref_name}: {len(full_align)} aligned events, "
                f"events/base={mapping['events_per_base']:.2f}"
            )
        else:
            print(f"  {ref_name}: No alignment")


def example_multi_read_multi_ref():
    """Example: Multiple reads against multiple references"""
    print("\n" + "=" * 70)
    print("Example 3: Multiple Reads vs Multiple References")
    print("=" * 70)

    test_dir = Path(__file__).parent
    fasta_path = test_dir / "one_read.fa"
    pod5_path = test_dir / "one_read.pod5"

    if not fasta_path.exists():
        print(f"FASTA file not found: {fasta_path}")
        return

    # Load reference
    ref_name, ref_seq, ref_len = load_reference_from_fasta(fasta_path)

    # Load signal
    _, signal_template, sample_rate = load_signal_from_pod5(pod5_path)

    # Create multiple "reads" by slicing the signal
    # In practice, these would be different actual reads
    read_signals = [
        signal_template[:50000],  # First 50k samples
        signal_template[50000:100000],  # Next 50k samples
    ]
    read_ids = ["read_001", "read_002"]
    read_seqs = [
        ref_seq[:500],
        ref_seq[500:1000],
    ]

    # Create multiple references
    ref_seqs = [
        ref_seq[:1000],
        ref_seq[500:1500],
    ]
    ref_names = [f"{ref_name}_region1", f"{ref_name}_region2"]
    ref_lens = [len(r) for r in ref_seqs]

    print(f"\nRunning pair-wise alignment:")
    print(f"  {len(read_ids)} reads vs {len(ref_seqs)} references")

    result = run_eventalign(
        read_ids=read_ids,
        read_seqs=read_seqs,
        ref_seqs=ref_seqs,
        ref_names=ref_names,
        ref_lens=ref_lens,
        signals=read_signals,
        sample_rates=[sample_rate] * len(read_ids),
        model_id=MODEL_RNA002,
    )

    print(f"\nPair-wise results:")
    for i, read_id in enumerate(read_ids):
        print(f"\n  {read_id}:")
        for j, ref_name in enumerate(ref_names):
            full_align = result["full"][i][j]
            mapping = result["mapping"][0][j]
            if len(full_align) > 0:
                print(
                    f"    vs {ref_name}: {len(full_align)} aligned events, "
                    f"events/base={mapping['events_per_base']:.2f}"
                )
            else:
                print(f"    vs {ref_name}: No alignment")


def main():
    """Main function"""
    print("run_eventalign Example")
    print("=" * 70)

    try:
        # Example 1: Single read vs single reference
        example_single_read_single_ref()

        # Example 2: Single read vs multiple references
        # example_single_read_multi_ref()

        # Example 3: Multiple reads vs multiple references
        # example_multi_read_multi_ref()

        print("\n" + "=" * 70)
        print("All examples completed successfully!")
        print("=" * 70)

    except Exception as e:
        print(f"\nError: {e}")
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
