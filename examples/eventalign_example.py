#!/usr/bin/env python3
"""
Example usage of f5c eventalign functionality

This demonstrates how to use the EventAligner class to align nanopore
events to reference k-mers.
"""

import sys
from pathlib import Path
import logging

# Set debug level (choose one of these options):

# Option 1: Set global debug level
# logging.basicConfig(level=logging.DEBUG)

# Option 2: Use the package logger with debug level
from fin.utils.log_config import setup_logger
logger = setup_logger(__name__, level='DEBUG', log_file='eventalign_debug.log')

try:
    from fin._f5c import EventAligner, _F5C_AVAILABLE
except ImportError as e:
    print("Error: f5c module not available")
    print(f"Import error: {e}")
    print("\nTo fix this:")
    print("  1. Make sure you have dependencies: sudo apt-get install zlib1g-dev")
    print("  2. Install with: pip install -e .")
    sys.exit(1)


def example_usage():
    """Example usage of EventAligner"""
    print("=" * 70)
    print("f5c EventAlign Example")
    print("=" * 70)

    # Note: These are placeholder paths - replace with your actual files
    bam_path = "path/to/your/aligned_reads.bam"
    fasta_path = "path/to/your/reference.fa"
    slow5_path = "path/to/your/signals.blow5"  # Optional but recommended

    print(f"\nInput files:")
    print(f"  BAM: {bam_path}")
    print(f"  FASTA: {fasta_path}")
    print(f"  SLOW5: {slow5_path}")

    # Check if files exist
    for path in [bam_path, fasta_path]:
        if not Path(path).exists():
            print(f"\n✗ File not found: {path}")
            print("  Please update the paths in this example script")
            return

    if slow5_path and not Path(slow5_path).exists():
        print(f"\n⚠ Warning: SLOW5 file not found: {slow5_path}")
        print("  Event alignment will work but may be slower")
        slow5_path = None

    try:
        # Create EventAligner
        print("\nInitializing EventAligner...")
        aligner = EventAligner(bam_path, fasta_path, slow5_path)

        print("✓ EventAligner initialized successfully")

        # Process reads in batches
        print("\nAligning reads (showing first 5 reads)...")
        print("-" * 70)

        read_count = 0
        total_events = 0

        for read_id, alignments in aligner.align_reads():
            read_count += 1
            total_events += len(alignments)

            if read_count <= 5:  # Show first 5 reads
                print(f"\nRead {read_count}: {read_id}")
                print(f"  Events aligned: {len(alignments)}")

                # Show first few alignments
                if alignments:
                    print("\n  First 3 event alignments:")
                    print("  " + "-" * 65)
                    print("  {:<8} {:<8} {:<10} {:>8} {:>8} {:>8}".format(
                        "Event#", "KmerIdx", "Kmer", "Eventμ", "Modelμ", "Prob"
                    ))
                    print("  " + "-" * 65)

                    for i, event in enumerate(alignments[:3]):
                        print("  {event_idx:<8} {kmer_idx:<8} {kmer:<10} "
                              "{event_mean:>8.2f} {model_mean:>8.2f} "
                              "{posterior_probability:>8.3f}".format(**event))

            # Stop after processing a reasonable number
            if read_count >= 100:
                print(f"\n... processed {read_count} reads (stopping demo)")
                break

        print("-" * 70)
        print(f"\nSummary:")
        print(f"  Reads processed: {read_count}")
        print(f"  Total events aligned: {total_events}")
        print(f"  Average events per read: {total_events / read_count:.1f}")

    except Exception as e:
        print(f"\n✗ Error during event alignment: {e}")
        print("\nTroubleshooting:")
        print("  1. Make sure BAM file is sorted and indexed")
        print("  2. Verify FASTA reference matches BAM alignments")
        print("  3. Check that SLOW5 file contains signals for the reads")
        raise


def example_batch_processing():
    """Example of processing alignments in batches"""
    print("\n" + "=" * 70)
    print("Batch Processing Example")
    print("=" * 70)

    bam_path = "path/to/your/aligned_reads.bam"
    fasta_path = "path/to/your/reference.fa"

    try:
        aligner = EventAligner(bam_path, fasta_path)

        print("Processing reads in batches...")
        batch_num = 0

        while True:
            batch_num += 1
            alignments = aligner.align_batch()

            if not alignments:
                break

            print(f"  Batch {batch_num}: {len(alignments)} alignments")

            # Process the alignments
            for alignment in alignments:
                # Your processing logic here
                pass

            if batch_num >= 10:  # Demo only
                print(f"\n... stopping demo after {batch_num} batches")
                break

    except Exception as e:
        print(f"Error: {e}")


def example_simple_interface():
    """Example using the simple eventalign() generator function"""
    print("\n" + "=" * 70)
    print("Simple Interface Example")
    print("=" * 70)

    bam_path = "test_reads.bam"
    fasta_path = "test_ref.fa"
    slow5_path = "test_signals.blow5"

    try:
        # Using the generator function (less control but simpler)
        for read_id, alignments in EventAligner.eventalign(
            bam_path, fasta_path, slow5_path, max_reads=10
        ):
            print(f"\nRead: {read_id}")
            for event in alignments[:3]:  # Show first 3
                print(f"  {event['kmer']}: μ={event['event_mean']:.2f}, "
                      f"model μ={event['model_mean']:.2f}, "
                      f"prob={event['posterior_probability']:.3f}")

    except Exception as e:
        print(f"Error: {e}")


def check_compilation():
    """Check if the module compiled correctly"""
    print("=" * 70)
    print("Compilation Check")
    print("=" * 70)

    if _F5C_AVAILABLE:
        print("✓ f5c module is available")
        print(f"  Module: {EventAligner.__module__}")
        print(f"  Class: {EventAligner.__name__}")

        # Show available methods
        print("\n✓ Available methods:")
        for method in ['__init__', 'align_batch', 'align_reads', 'close']:
            if hasattr(EventAligner, method):
                print(f"  - {method}")

        print("\n✓ Ready to use!")
        return True
    else:
        print("✗ f5c module not available")
        print("\nTo compile:")
        print("  1. Install dependencies:")
        print("     sudo apt-get install zlib1g-dev libhts-dev")
        print("  2. Build the module:")
        print("     pip install -e .")
        print("\nTroubleshooting:")
        print("  - Check logs/build.log for compilation errors")
        print("  - Make sure third_party/f5c/ exists (git submodule)")
        return False


def main():
    """Main function"""
    print("f5c EventAlign Example Script")
    print("=" * 70)
    print("\nThis script demonstrates how to use the f5c eventalign functionality")
    print("from Python.\n")

    # Check compilation first
    if not check_compilation():
        return

    print("\n" + "=" * 70)
    print("Examples (tests with placeholder paths):")
    print("=" * 70)
    print("\n1. Basic usage (EventAligner class)")
    print("2. Batch processing")
    print("3. Simple generator interface")
    print("4. Run with real files (modify paths first)\n")

    # Show examples without running (since we don't have real test data)
    print("Example 1: Basic usage")
    print("-" * 70)
    print("""
    from fin._f5c import EventAligner

    aligner = EventAligner(
        bam_path="aligned_reads.bam",
        fasta_path="reference.fa",
        slow5_path="signals.blow5"  # Optional
    )

    for read_id, alignments in aligner.align_reads(max_reads=100):
        print(f"Read: {read_id}, Events: {len(alignments)}")
        for aln in alignments:
            print(f"  {aln['kmer']}: {aln['event_mean']:.2f}")
    """)

    print("\nExample 2: Using context manager")
    print("-" * 70)
    print("""
    with EventAligner(bam_path, fasta_path, slow5_path) as aligner:
        for read_id, alignments in aligner.align_reads():
            process_alignments(read_id, alignments)
    # Resources automatically released
    """)

    print("\n" + "=" * 70)
    print("Next Steps:")
    print("-" * 70)
    print("1. Update the file paths in this script")
    print("2. Run: python example_eventalign.py")
    print("3. Integrate into your analysis pipeline")
    print("\nFor more information, see the docstrings in fin/_f5c/__init__.py")
    print("=" * 70)


if __name__ == "__main__":
    main()
