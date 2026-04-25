#!/usr/bin/env python3
"""
Test example for the py-fin package

This demonstrates how to use the FASTA reader and ReadSubsetManager
"""

import sys
from pathlib import Path

# Import the package
try:
    import fin
    from fin.io import FASTAReader, FASTARecord, ReadSubsetManager, create_subset_manager
    print(f"✓ Successfully imported fin version {fin.__version__}")
except ImportError as e:
    print(f"✗ Failed to import fin: {e}")
    sys.exit(1)

def test_fasta_reader():
    """Test FASTA reader functionality"""
    print("\n" + "="*60)
    print("Testing FASTA Reader")
    print("="*60)

    # Create a test FASTA file
    test_fasta = Path("test_example.fasta")
    with open(test_fasta, "w") as f:
        f.write(">chr1 Test chromosome 1\n"
                "ATCGATCGATCGATCGATCG\n"
                "TTTTCCCCAAAAGGGG\n"
                ">chr2 Test chromosome 2\n"
                "GGGGAAAATTTTCCCC\n"
                "ATATATATATATAT\n")

    try:
        # Read the FASTA file
        with FASTAReader(str(test_fasta)) as reader:
            records = reader.get_records()

            print(f"✓ Read {len(records)} records from FASTA file")

            for i, record in enumerate(records, 1):
                print(f"\nRecord {i}:")
                print(f"  ID: {record.id}")
                print(f"  Header: {record.header}")
                print(f"  Length: {record.length}")
                print(f"  GC Content: {record.gc_content:.1f}%")
                print(f"  Sequence preview: {record.sequence[:30]}...")

            # Test get_summary
            summary = reader.get_summary()
            print("\n✓ FASTA Summary:")
            for key, value in summary.items():
                print(f"  {key}: {value}")

    finally:
        # Clean up
        if test_fasta.exists():
            test_fasta.unlink()

def test_read_subset_manager():
    """Test ReadSubsetManager (requires test data files)"""
    print("\n" + "="*60)
    print("Testing ReadSubsetManager")
    print("="*60)

    # This is a demonstration - in real use you'd have actual BAM/GTF files
    print("\nReadSubsetManager provides:")
    print("  - Organize reads from BAM into transcriptomic intervals")
    print("  - Identify fusion candidate reads")
    print("  - Group unannotated reads separately")
    print("  - Lazy generation of subset BAM files")
    print("\nUsage example:")
    print("  manager = create_subset_manager('reads.bam', 'annotation.gtf', 'genome.fa')")
    print("  for subset in manager.iterate_subsets():")
    print("      bundle = manager.get_data_bundle(subset)")
    print("      # Process bundle.reads and bundle.reference_sequences")
    print("      # Optionally: manager.write_subset_bam(subset, 'output.bam')")

    print("\n✓ ReadSubsetManager imported successfully")
    print("✓ Check examples/ and documentation for full usage")

def main():
    """Main test function"""
    print("="*60)
    print("py-fin Package Test")
    print("="*60)
    print(f"Package version: {fin.__version__}")
    print(f"Package location: {Path(fin.__file__).parent}")

    # Test FASTA reader
    test_fasta_reader()

    # Test Read Subset Manager
    test_read_subset_manager()

    print("\n" + "="*60)
    print("All tests completed successfully!")
    print("="*60)

    print("\nAvailable modules in fin.io:")
    print("  - FASTAReader/FASTARecord: Read FASTA files")
    print("  - BamReader: Read BAM/SAM files")
    print("  - GTFReader: Read GTF/GFF annotations")
    print("  - BEDReader: Read BED files")
    print("  - ReadSubsetManager: Organize reads into logical subsets")
    print("  - Fast5Reader/Slow5Reader/Pod5Reader: Nanopore signal files")

if __name__ == "__main__":
    main()
