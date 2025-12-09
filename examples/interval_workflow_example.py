#!/usr/bin/env python3
"""
Example workflow for using the new interval manager approach

This demonstrates how to:
1. Generate isolated, non-overlapping intervals from BAM/GTF
2. Skip fusion reads during interval generation
3. Extract reads and annotation on-demand for specific intervals
4. Strand-separated intervals with read count tracking
"""

import sys
from pathlib import Path

# Import the package
from fin.io import (
    generate_isolated_intervals,
    extract_reads_for_interval,
    extract_annotation_for_interval,
    GenomicInterval
)
from fin.utils.log_config import get_package_logger


def example_workflow():
    """
    Example workflow showing how to use the new interval-based approach

    Note: This uses placeholder file paths - replace with actual files
    """
    logger = get_package_logger(__name__)

    print("\n" + "=" * 60)
    print("Interval Generation Workflow Example")
    print("=" * 60)

    # Step 1: Generate isolated intervals
    # ====================================
    print("\nStep 1: Generating isolated intervals...")
    print("-" * 60)

    bam_path = "path/to/your/reads.bam"
    gtf_path = "path/to/your/annotation.gtf"

    result = generate_isolated_intervals(
        bam_path=bam_path,
        gtf_path=gtf_path,
        max_gap=1000,  # Merge intervals within 1kb
        max_reads=None  # Process all reads
    )

    intervals = result['intervals']
    fusion_read_ids = result['fusion_read_ids']

    print(f"✓ Generated {len(intervals)} isolated intervals")
    print(f"✓ Identified {len(fusion_read_ids)} fusion reads")
    print(f"✓ Processed {result['num_reads_processed']} total reads")

    # Show strand separation
    strand_counts = {}
    for interval in intervals:
        strand_counts[interval.strand] = strand_counts.get(interval.strand, 0) + 1
    print("\nInterval distribution by strand:")
    for strand, count in sorted(strand_counts.items()):
        print(f"  {strand if strand else 'unstranded'}: {count} intervals")

    # Step 2: Process intervals in batches
    # ====================================
    print("\nStep 2: Processing intervals...")
    print("-" * 60)

    batch_size = 100  # Process 100 intervals at a time

    # Show interval stats
    print(f"\nProcessing {len(intervals)} intervals in batches of {batch_size}")
    print("=" * 60)

    total_extracted_reads = 0

    for batch_start in range(0, len(intervals), batch_size):
        batch_end = min(batch_start + batch_size, len(intervals))
        batch_intervals = intervals[batch_start:batch_end]

        print(f"\n[Batch {batch_start // batch_size + 1}]")
        print(f"Processing intervals {batch_start}-{batch_end - 1}")
        print(f"\n{'Index':<8} {'Interval':<25} {'Strand':<10} {'Interval Reads':<15} {'Extracted':<12}")
        print("-" * 72)

        # Step 3: Extract data for each interval on-demand
        # ================================================
        for i, interval in enumerate(batch_intervals):
            idx = batch_start + i

            # Extract reads for this interval (skipping fusions)
            reads = extract_reads_for_interval(
                bam_path=bam_path,
                interval=interval,
                fusion_read_ids=fusion_read_ids
            )

            # Extract annotation for this interval
            annotation = extract_annotation_for_interval(
                gtf_path=gtf_path,
                interval=interval
            )

            total_extracted_reads += len(reads)

            # Display interval info with strand separation
            strand_display = interval.strand if interval.strand else "N/A"
            print(f"{idx:<8} {interval.region_string:<25} {strand_display:<10} "
                  f"{interval.read_count:<15} {len(reads):<12}")

            # Your processing logic here...
            # For example:
            # - Align reads to transcripts
            # - Detect modifications
            # - Quantify expression

    print(f"\n✓ Processed all {len(intervals)} intervals")
    print(f"✓ Total extracted reads: {total_extracted_reads}")
    print(f"✓ Average reads per interval: {total_extracted_reads / len(intervals):.1f}")

    # Alternative: Process specific intervals by chromosome
    # =====================================================
    print("\nAlternative: Process intervals by chromosome...")
    print("-" * 60)

    # Group intervals by chromosome
    chrom_intervals = {}
    for interval in intervals:
        if interval.chrom not in chrom_intervals:
            chrom_intervals[interval.chrom] = []
        chrom_intervals[interval.chrom].append(interval)

    # Process each chromosome separately
    for chrom, chrom_interval_list in chrom_intervals.items():
        print(f"\nProcessing chromosome {chrom}: {len(chrom_interval_list)} intervals")

        # You could parallelize this per chromosome
        # or save chromosome-specific data

        # Example: Save chromosome-level summary
        total_reads = 0
        for interval in chrom_interval_list:
            reads = extract_reads_for_interval(bam_path, interval, fusion_read_ids)
            total_reads += len(reads)

        print(f"  Total reads on {chrom}: {total_reads}")

    print("\n✓ Completed chromosome-level processing")


def example_parallel_processing():
    """
    Example showing how to use parallel processing with intervals
    """
    print("\n" + "=" * 60)
    print("Parallel Processing Example")
    print("=" * 60)

    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp

    def process_interval(interval_info):
        """Process a single interval (for parallel execution)"""
        bam_path, interval, fusion_read_ids = interval_info

        # Extract reads for this interval
        reads = extract_reads_for_interval(
            bam_path=bam_path,
            interval=interval,
            fusion_read_ids=fusion_read_ids
        )

        # Return results (must be pickle-able)
        return {
            'interval': interval.region_string,
            'num_reads': len(reads)
        }

    # Generate intervals (same as before)
    bam_path = "path/to/your/reads.bam"
    result = generate_isolated_intervals(bam_path=bam_path)
    intervals = result['intervals']
    fusion_read_ids = result['fusion_read_ids']

    # Process intervals in parallel
    num_workers = min(mp.cpu_count(), len(intervals))

    # Prepare data for parallel processing
    process_data = [
        (bam_path, interval, fusion_read_ids)
        for interval in intervals
    ]

    print(f"Processing {len(intervals)} intervals in parallel using {num_workers} workers...")

    results = []
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Process intervals in chunks to avoid memory issues
        chunk_size = 100
        for i in range(0, len(process_data), chunk_size):
            chunk = process_data[i:i + chunk_size]
            chunk_results = executor.map(process_interval, chunk)
            results.extend(chunk_results)

            print(f"  Processed {len(results)} intervals so far...")

    print(f"\n✓ Completed parallel processing of {len(results)} intervals")

    # Summarize results
    total_reads = sum(r['num_reads'] for r in results)
    print(f"✓ Total reads across all intervals: {total_reads}")


def main():
    """Main example function"""
    logger = get_package_logger(__name__)

    print("\n" + "=" * 60)
    print("Interval Manager Workflow Examples")
    print("=" * 60)
    print("\nNote: These examples use placeholder file paths.")
    print("Replace with actual files to run the examples.")

    # Uncomment one of these examples to try:

    # Example 1: Basic workflow
    # example_workflow()

    # Example 2: Parallel processing
    # example_parallel_processing()

    print("\n" + "=" * 60)
    print("Key Benefits of This Approach:")
    print("-" * 60)
    print("1. ✓ Intervals are isolated and non-overlapping")
    print("2. ✓ Strand-separated (+ and - strands kept separate)")
    print("3. ✓ Fusion reads are identified in first pass")
    print("4. ✓ Read counts tracked per interval (no gene IDs)")
    print("5. ✓ Reads/annotation extracted on-demand")
    print("6. ✓ Memory-efficient (no caching required)")
    print("7. ✓ Easy to parallelize per interval")
    print("8. ✓ Flexible processing order")
    print("=" * 60)


if __name__ == "__main__":
    main()
