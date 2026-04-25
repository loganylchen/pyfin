#!/usr/bin/env python3
"""
Debug script to diagnose alignment failures.
"""

import numpy as np
from pathlib import Path
from fin._eventalign import run_eventalign, MODEL_RNA002

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


def load_signal_from_pod5(pod5_path: str) -> tuple:
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
        # Generate synthetic signal
        print("Generating synthetic signal...")
        np.random.seed(42)
        n_samples = 100000
        signal = np.random.randn(n_samples).astype(np.float32) * 10 + 120
        sample_rate = 4000.0
        return "synthetic_read", signal, sample_rate


def main():
    """Debug alignment failure."""
    print("=" * 70)
    print("Debugging Alignment Failure")
    print("=" * 70)

    # Paths to test data
    test_dir = Path(__file__).parent / "test_data"
    fasta_path = test_dir / "one_read.fa"
    pod5_path = test_dir / "one_read.pod5"

    # Load reference
    print(f"\nLoading reference from: {fasta_path}")
    ref_name, ref_seq, ref_len = load_reference_from_fasta(fasta_path)
    print(f"  Reference: {ref_name}")
    print(f"  Length: {ref_len} bp")
    print(f"  First 50bp: {ref_seq[:50]}")

    # Load signal
    print(f"\nLoading signal from: {pod5_path}")
    read_id, signal, sample_rate = load_signal_from_pod5(pod5_path)
    print(f"  Read ID: {read_id}")
    print(f"  Signal length: {len(signal)} samples")
    print(f"  Sample rate: {sample_rate} Hz")
    print(f"  Signal stats: mean={signal.mean():.2f}, std={signal.std():.2f}, min={signal.min():.2f}, max={signal.max():.2f}")

    # IMPORTANT: The read sequence MUST match the signal!
    # Using the reference sequence as a substitute will NOT work because
    # the signal was generated from a different molecule/sequence.
    #
    # For real data, you need to provide the basecalled read sequence
    # that corresponds to this specific signal. Load it from:
    # - BAM file (see load_from_bam_pod5.py)
    # - FASTQ file (see load_from_fastq_pod5.py)
    #
    # For this debug script, we'll use the reference as a placeholder
    # but alignment will likely fail since it doesn't match the signal.
    print("\n" + "!" * 70)
    print("WARNING: Using reference sequence as read sequence!")
    print("This will likely cause alignment failure because the read_seq")
    print("must match the actual basecalled sequence from the signal.")
    print("!")
    print("For proper workflow, see:")
    print("  - load_from_fastq_pod5.py (use FASTQ from basecaller)")
    print("  - load_from_bam_pod5.py (use BAM with aligned reads)")
    print("!" * 70)

    read_seq = ref_seq

    print("\nRunning eventalign...")
    print("  Model: RNA002 (k=5)")
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

    # Display detailed diagnostics
    print("\n" + "=" * 70)
    print("Detailed Diagnostics")
    print("=" * 70)

    # Summary
    summary = result["summary"]
    print("\nSummary:")
    print(f"  Reads processed: {summary['num_reads']}")
    print(f"  References: {summary['num_refs']}")

    # Scalings
    print("\nScalings:")
    for i, sc in enumerate(result["scalings"]):
        print(f"  Read {i}:")
        print(f"    scale:  {sc['scale']:.6f}")
        print(f"    shift:  {sc['shift']:.6f}")
        print(f"    var:    {sc['var']:.6f}")

    # Events
    print("\nEvents:")
    events = result["events"][0]
    starts = events['starts']
    n_ev = len(starts)
    print(f"  Total events: {n_ev}")
    print(f"  Events/kmer ratio: {n_ev / (ref_len - 5 + 1):.2f}")

    # Event statistics
    lengths = events['lengths']
    means = events['means']
    stdvs = events['stdvs']

    print("\n  Event Statistics:")
    print(f"    Start positions: min={starts.min()}, max={starts.max()}")
    print(f"    Lengths: mean={lengths.mean():.2f}, std={lengths.std():.2f}, min={lengths.min():.2f}, max={lengths.max():.2f}")
    print(f"    Means:   mean={means.mean():.2f}, std={means.std():.2f}, min={means.min():.2f}, max={means.max():.2f}")
    print(f"    Stdv:    mean={stdvs.mean():.4f}, std={stdvs.std():.4f}, min={stdvs.min():.4f}, max={stdvs.max():.4f}")

    # Show first and last 10 events
    print("\n  First 10 events:")
    print(f"    {'Idx':>4} {'Start':>10} {'Length':>8} {'Mean':>8} {'Stdv':>8}")
    print(f"    {'-'*4} {'-'*10} {'-'*8} {'-'*8} {'-'*8}")
    for i in range(min(10, n_ev)):
        print(f"    {i:4d} {starts[i]:10d} {lengths[i]:8.1f} {means[i]:8.2f} {stdvs[i]:8.4f}")

    print("\n  Last 10 events:")
    print(f"    {'Idx':>4} {'Start':>10} {'Length':>8} {'Mean':>8} {'Stdv':>8}")
    print(f"    {'-'*4} {'-'*10} {'-'*8} {'-'*8} {'-'*8}")
    for i in range(max(0, n_ev - 10), n_ev):
        print(f"    {i:4d} {starts[i]:10d} {lengths[i]:8.1f} {means[i]:8.2f} {stdvs[i]:8.4f}")

    # Alignment results
    print("\nAlignment Results:")
    full_align = result["full"][0][0]
    mapping = result["mapping"][0][0]

    if len(full_align) > 0:
        print(f"  Alignment SUCCESS: {len(full_align)} aligned events")
        print("\n  First 5 aligned events:")
        print(f"    {'RefPos':>8} {'EventIdx':>10} {'RefKmer':>10} {'State':>6}")
        print(f"    {'-'*8} {'-'*10} {'-'*10} {'-'*6}")
        for k in range(min(5, len(full_align))):
            ea = full_align[k]
            print(f"    {ea['ref_position']:8d} {ea['event_idx']:10d} {ea['ref_kmer']:>10} {ea['hmm_state']:>6}")

        start = mapping["start"]
        stop = mapping["stop"]
        print("\n  Base-to-Event Mapping:")
        print(f"    Events per base: {mapping['events_per_base']:.2f}")
        print(f"    First 10 bases: start={[int(s) for s in start[:10]]}")
        print(f"    Last 10 bases:  start={[int(s) for s in start[-10:]]}")
    else:
        print("  Alignment FAILED")
        print(f"    Status: {mapping.get('status', 'unknown')}")
        print(f"    Events detected: {mapping.get('n_events', 'N/A')}")
        print(f"    K-mers in reference: {mapping.get('n_kmers', 'N/A')}")
        print(f"    Reference length: {mapping.get('ref_len', 'N/A')} bp")
        print(f"    Read length: {mapping.get('read_len', 'N/A')} bp")

        # Additional debugging
        print("\n  Possible Failure Reasons:")
        n_kmers = mapping.get('n_kmers', ref_len - 5 + 1)

        # Check events per kmer ratio
        epk = n_ev / n_kmers
        print(f"    1. Events/k-mer ratio: {epk:.2f} (typical: 2-4 for RNA)")

        # Check scaling quality
        sc = result["scalings"][0]
        # Expected: scale around 1.0, shift around 0 (or adjusted for signal level)
        print(f"    2. Scaling: scale={sc['scale']:.4f} (expected ~1.0)")
        print(f"               shift={sc['shift']:.4f} (adjusts signal level)")
        print("               (After scaling, event means should match model levels)")

        # Check if model levels match scaled event means
        from fin._eventalign import set_model
        model = set_model(MODEL_RNA002)
        model_means = model['level_means']
        print(f"    3. Model level range: [{model_means.min():.2f}, {model_means.max():.2f}]")
        scaled_means = (means - sc['shift']) / sc['scale']
        print(f"       Scaled event mean range: [{scaled_means.min():.2f}, {scaled_means.max():.2f}]")

        # Check signal quality
        print(f"    4. Signal quality: std={signal.std():.2f} (should be > 5 for good data)")
        print(f"       Signal drift: first_mean={signal[:1000].mean():.2f}, last_mean={signal[-1000:].mean():.2f}")

        print("\n  Common Issues:")
        print("    - MISSING READ SEQUENCE: POD5 files don't contain basecalled sequences!")
        print("      Solution: Use FASTQ (from Guppy/Dorado) or BAM files")
        print("      See: load_from_fastq_pod5.py or load_from_bam_pod5.py")
        print("    - Wrong read sequence: read_seq must match the actual signal")
        print("      The reference sequence is NOT the same as the read sequence!")
        print("    - Synthetic/random signal: Will NOT produce good alignment")
        print("    - Poor signal quality: Low variance or extreme drift")
        print("    - Wrong pore model: RNA002 vs RNA004")


if __name__ == "__main__":
    main()
