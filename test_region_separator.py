#!/usr/bin/env python3
"""
Region Separator Testing Script

This script demonstrates how to use the RegionSeparator class to split BAM files
by gene regions and filter out fusion-like reads.
"""

import sys
import os
from pathlib import Path
import logging

# Add the project root to Python path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

from fin.io.region_separator import RegionSeparator, FusionReadDetector, separate_regions


def setup():
    """Setup logging and imports."""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    print("Set up logging and imports completed.\n")


def example_basic_usage():
    """
    Example 1: Basic Usage

    Start with a basic example of separating reads by gene regions.
    """
    print("=" * 60)
    print("Example 1: Basic Usage")
    print("=" * 60)

    # Use test data files
    bam_file = "examples/test_data/test_reads.bam.gz"  # Gzipped BAM file
    gtf_file = "examples/test_data/test_annotation.gtf.gz"  # Gzipped GTF file
    output_dir = "output/regions"  # Output directory

    if os.path.exists(bam_file) and os.path.exists(gtf_file):
        separator = RegionSeparator(
            bam_file=bam_file,
            gtf_file=gtf_file,
            output_dir=output_dir,
            detect_fusions=True  # Enable fusion read detection
        )

        # Extract gene regions (returns iterator)
        print("Extracting gene regions...")
        for i, region in enumerate(separator.extract_gene_regions(min_reads=5)):
            print(f"Gene {i+1}: {region.gene_id} ({region.gene_name or 'N/A'})")
            print(f"  Location: {region.chrom}:{region.start}-{region.end}")
            print(f"  Strand: {region.strand}")
            print(f"  Reads: {region.read_count}")
            print(f"  Transcripts: {len(region.transcripts)}")
            print()

            # Only show first 3 genes to avoid too much output
            if i >= 2:
                print("... (truncated after 3 genes)")
                break
    else:
        print("Please update the file paths above with your actual BAM and GTF files!")

    print()


def example_fusion_detection():
    """
    Example 2: Fusion Read Detection

    Explore how the fusion read detection works.
    """
    print("=" * 60)
    print("Example 2: Fusion Read Detection")
    print("=" * 60)

    bam_file = "examples/test_data/test_reads.bam.gz"  # Test BAM file

    # Create a fusion detector with custom parameters
    fusion_detector = FusionReadDetector(
        min_alignment_blocks=2,      # Consider reads with ≥2 alignment blocks
        max_block_gap=100000,        # Flag if gap > 100kb between blocks
        min_block_size=50,           # Minimum block size
        max_trimmed_bases=20         # Max soft-clipped bases at both ends
    )

    if os.path.exists(bam_file):
        import pysam
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            fusion_count = 0
            total_count = 0

            # Check first 1000 reads
            for i, record in enumerate(bam):
                if i >= 1000:
                    break

                total_count += 1
                if fusion_detector.is_fusion_like(record):
                    fusion_count += 1

            print(f"Analyzed {total_count} reads")
            print(f"Detected {fusion_count} fusion-like reads ({fusion_count/total_count*100:.2f}%)")
    else:
        print("Please update the bam_file path above with your actual BAM file!")

    print()


def example_writing_output():
    """
    Example 3: Writing Output BAM Files

    Write separate BAM files for each gene region with sufficient read coverage.
    """
    print("=" * 60)
    print("Example 3: Writing Output BAM Files")
    print("=" * 60)

    bam_file = "examples/test_data/test_reads.bam.gz"
    gtf_file = "examples/test_data/test_annotation.gtf.gz"
    output_dir = "output/regions"

    if os.path.exists(bam_file) and os.path.exists(gtf_file):
        separator = RegionSeparator(
            bam_file=bam_file,
            gtf_file=gtf_file,
            output_dir=output_dir,
            detect_fusions=True
        )

        # Write BAM files for genes with at least 10 reads
        print("Writing BAM files to separate regions...")
        created_files = separator.write_region_bams(min_reads=10, prefix="gene")

        print(f"\nCreated {len(created_files)} BAM files")
        print("\nFile locations:")
        for f in created_files:
            print(f"  {f}")
    else:
        print("Please update the file paths above with your actual files!")

    print()


def example_convenient_function():
    """
    Example 4: Convenient Function

    Use the convenient `separate_regions` function for a one-step solution.
    """
    print("=" * 60)
    print("Example 4: Convenient Function")
    print("=" * 60)

    bam_file = "examples/test_data/test_reads.bam.gz"
    gtf_file = "examples/test_data/test_annotation.gtf.gz"
    output_dir = "output/regions"

    if os.path.exists(bam_file) and os.path.exists(gtf_file):
        print("Separating reads by gene regions...")
        print(f"Input BAM: {bam_file}")
        print(f"Input GTF: {gtf_file}")
        print(f"Output directory: {output_dir}")

        # Call the convenient function
        created_files = separate_regions(
            bam_file=bam_file,
            gtf_file=gtf_file,
            output_dir=output_dir,
            min_reads=5,               # Minimum reads per gene
            detect_fusions=True,       # Filter fusion reads
            write_bams=True,           # Write individual BAMs
            region_list="output/regions.tsv"
        )

        print(f"\nSuccessfully created {len(created_files)} BAM files")
        print(f"(Saved to: {output_dir}.)")
    else:
        print("Please update the file paths above with your actual files!")

    print()


def example_region_statistics():
    """
    Example 5: Examining Region Statistics

    Get detailed statistics about the separated regions.
    """
    print("=" * 60)
    print("Example 5: Examining Region Statistics")
    print("=" * 60)

    bam_file = "examples/test_data/test_reads.bam.gz"
    gtf_file = "examples/test_data/test_annotation.gtf.gz"
    output_dir = "output/regions"

    if os.path.exists(bam_file) and os.path.exists(gtf_file):
        separator = RegionSeparator(
            bam_file=bam_file,
            gtf_file=gtf_file,
            output_dir=output_dir,
            detect_fusions=True
        )

        # Get statistics
        separator.print_statistics()
        stats = separator.get_statistics()

        # Access detailed stats
        print("\nDetailed Statistics:")
        print(f"- Total reads processed: {stats['total_reads']}")
        print(f"- Fusion-like reads removed: {stats['fusion_reads']}")
        print(f"- Genes with reads: {stats['genes_with_reads']}/{stats['total_genes']}")

        # Find genes with most reads
        sorted_genes = sorted(stats['reads_per_gene'].items(),
                              key=lambda x: x[1], reverse=True)[:10]
        print("\nTop 10 genes by read count:")
        for rank, (gene_id, read_count) in enumerate(sorted_genes, 1):
            print(f"{rank}. {gene_id}: {read_count} reads")
    else:
        print("Please update the file paths above with your actual files!")

    print()


def example_custom_fusion_detector():
    """
    Example 6: Creating a Custom Fusion Detector

    Customize fusion detection parameters for your specific needs.
    """
    print("=" * 60)
    print("Example 6: Custom Fusion Detector")
    print("=" * 60)

    bam_file = "examples/test_data/test_reads.bam.gz"

    # Create a stricter fusion detector
    strict_fusion_detector = FusionReadDetector(
        min_alignment_blocks=3,       # Require at least 3 blocks
        max_block_gap=50000,          # More sensitive to gaps
        min_block_size=100,           # Larger minimum block size
        max_trimmed_bases=10          # Stricter on clipping
    )

    # Create a lenient fusion detector
    lenient_fusion_detector = FusionReadDetector(
        min_alignment_blocks=2,       # Only require 2 blocks
        max_block_gap=200000,         # Allow larger gaps
        min_block_size=30,            # Smaller minimum block
        max_trimmed_bases=50          # Allow more clipping
    )

    if os.path.exists(bam_file):
        import pysam
        with pysam.AlignmentFile(bam_file, 'rb') as bam:
            # Count with strict detector
            strict_count = sum(1 for r in bam if strict_fusion_detector.is_fusion_like(r))

            # Reset iterator
            bam.reset()

            # Count with lenient detector
            lenient_count = sum(1 for r in bam if lenient_fusion_detector.is_fusion_like(r))

            print(f"Strict fusion detector: {strict_count} reads flagged")
            print(f"Lenient fusion detector: {lenient_count} reads flagged")
    else:
        print("Please update the bam_file path above with your actual BAM file!")

    print()


def main():
    """Run all examples."""
    print("\n" + "=" * 60)
    print("Region Separator Testing Examples")
    print("=" * 60)
    print()

    setup()

    # Note: These examples use placeholder file paths.
    # Update the bam_file and gtf_file variables in each function
    # with your actual file paths to run the examples.

    example_basic_usage()
    example_fusion_detection()
    example_writing_output()
    example_convenient_function()
    example_region_statistics()
    example_custom_fusion_detector()

    print("\n" + "=" * 60)
    print("Testing script completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
