#!/usr/bin/env python3
"""
Generate example data for testing the region separator.

This script creates:
1. A small GTF file with a few test genes
2. A corresponding BAM file with artificial reads
3. Some fusion-like reads to test the filtering
"""

import os
import sys
from pathlib import Path

try:
    import pysam
    import numpy as np
except ImportError as e:
    print(f"Error: {e}")
    print("Please install required packages: pip install pysam numpy")
    sys.exit(1)


def generate_gtf(output_path: str):
    """Generate a simple GTF file with test genes."""
    gtf_content = """chr1\tunknown\tgene\t1000\t5000\t.\t+\t.\tgene_id "gene_A"; gene_name "TEST1";\nchr1\tunknown\texon\t1000\t1500\t.\t+\t.\tgene_id "gene_A"; transcript_id "transcript_A1"; exon_number "1";\nchr1\tunknown\texon\t2000\t2500\t.\t+\t.\tgene_id "gene_A"; transcript_id "transcript_A1"; exon_number "2";\nchr1\tunknown\texon\t4000\t5000\t.\t+\t.\tgene_id "gene_A"; transcript_id "transcript_A1"; exon_number "3";\nchr1\tunknown\texon\t1000\t1500\t.\t+\t.\tgene_id "gene_A"; transcript_id "transcript_A2"; exon_number "1";\nchr1\tunknown\texon\t3500\t5000\t.\t+\t.\tgene_id "gene_A"; transcript_id "transcript_A2"; exon_number "2";\nchr1\tunknown\tgene\t5500\t11000\t.\t-\t.\tgene_id "gene_B"; gene_name "TEST2";\nchr1\tunknown\texon\t5500\t6200\t.\t-\t.\tgene_id "gene_B"; transcript_id "transcript_B1"; exon_number "3";\nchr1\tunknown\texon\t8000\t8500\t.\t-\t.\tgene_id "gene_B"; transcript_id "transcript_B1"; exon_number "2";\nchr1\tunknown\texon\t10000\t11000\t.\t-\t.\tgene_id "gene_B"; transcript_id "transcript_B1"; exon_number "1";\nchr2\tunknown\tgene\t5000\t9000\t.\t+\t.\tgene_id "gene_C"; gene_name "TEST3";\nchr2\tunknown\texon\t5000\t6000\t.\t+\t.\tgene_id "gene_C"; transcript_id "transcript_C1"; exon_number "1";\nchr2\tunknown\texon\t7000\t7500\t.\t+\t.\tgene_id "gene_C"; transcript_id "transcript_C1"; exon_number "2";\nchr2\tunknown\texon\t8000\t9000\t.\t+\t.\tgene_id "gene_C"; transcript_id "transcript_C1"; exon_number "3";\n"""

    with open(output_path, 'w') as f:
        f.write(gtf_content)

    print(f"Generated GTF file: {output_path}")


def create_bam_header():
    """Create a BAM file header."""
    header = {
        'HD': {'VN': '1.0', 'SO': 'unsorted'},
        'SQ': [
            {'SN': 'chr1', 'LN': 249250621},
            {'SN': 'chr2', 'LN': 243199373},
        ]
    }
    return header


def generate_normal_reads(bam_file, read_prefix):
    """Generate normal alignment reads for testing."""
    # Normal read for gene A (chr1:1000-5000, + strand)
    r1 = pysam.AlignedSegment(bam_file.header)
    r1.query_name = f"{read_prefix}_001"
    r1.query_sequence = "A" * 200  # 200bp read
    r1.flag = 0  # positive strand
    r1.reference_id = 0  # chr1
    r1.reference_start = 1100
    r1.mapping_quality = 60
    r1.cigartuples = [(0, 200)]  # 200M
    r1.next_reference_id = -1
    r1.next_reference_start = -1
    r1.template_length = 0
    r1.query_qualities = pysam.qualitystring_to_array("I" * 200)

    # Another normal read for gene A (spans multiple exons)
    r2 = pysam.AlignedSegment(bam_file.header)
    r2.query_name = f"{read_prefix}_002"
    r2.query_sequence = "C" * 500
    r2.flag = 16  # negative strand
    r2.reference_id = 0
    r2.reference_start = 4800
    r2.mapping_quality = 60
    r2.cigartuples = [(0, 500)]
    r2.query_qualities = pysam.qualitystring_to_array("I" * 500)

    # Normal read for gene B (chr1:5500-11000, - strand)
    r3 = pysam.AlignedSegment(bam_file.header)
    r3.query_name = f"{read_prefix}_003"
    r3.query_sequence = "G" * 300
    r3.flag = 16
    r3.reference_id = 0
    r3.reference_start = 10500
    r3.mapping_quality = 55
    r3.cigartuples = [(0, 300)]
    r3.query_qualities = pysam.qualitystring_to_array("I" * 300)

    # Normal read for gene C (chr2:5000-9000, + strand)
    r4 = pysam.AlignedSegment(bam_file.header)
    r4.query_name = f"{read_prefix}_004"
    r4.query_sequence = "T" * 400
    r4.flag = 0
    r4.reference_id = 1  # chr2
    r4.reference_start = 5200
    r4.mapping_quality = 58
    r4.cigartuples = [(0, 400)]
    r4.query_qualities = pysam.qualitystring_to_array("I" * 400)

    return [r1, r2, r3, r4]


def generate_fusion_reads(bam_file, read_prefix):
    """Generate fusion-like reads for testing."""
    fusion_reads = []

    # Fusion-like: read with multiple alignment blocks (split alignment)
    # This simulates a chimeric read
    r1 = pysam.AlignedSegment(bam_file.header)
    r1.query_name = f"{read_prefix}_f001"
    r1.query_sequence = "A" * 400
    r1.flag = 2048  # supplementary alignment flag
    r1.reference_id = 0  # chr1
    r1.reference_start = 1100
    r1.mapping_quality = 20  # low mapq
    # 150M + 200N (large gap) + 250M - looks like a fusion!
    r1.cigartuples = [(0, 150), (3, 200), (0, 250)]  # 150 + 250 = 400
    r1.query_qualities = pysam.qualitystring_to_array("I" * 400)
    fusion_reads.append(r1)

    # Fusion-like: read with extreme soft clipping on both ends
    r2 = pysam.AlignedSegment(bam_file.header)
    r2.query_name = f"{read_prefix}_f002"
    r2.query_sequence = "C" * 1000
    r2.flag = 0
    r2.reference_id = 0
    r2.reference_start = 5000
    r2.mapping_quality = 5  # very low mapq
    # 400S + 200M + 400S - most bases are clipped! (1000 = 400+200+400)
    r2.cigartuples = [(4, 400), (0, 200), (4, 400)]
    r2.query_qualities = pysam.qualitystring_to_array("I" * 1000)
    fusion_reads.append(r2)

    # Normal spliced read (NOT fusion - this is typical for RNA-seq)
    r3 = pysam.AlignedSegment(bam_file.header)
    r3.query_name = f"{read_prefix}_003_spliced"
    r3.query_sequence = "G" * 500
    r3.flag = 0
    r3.reference_id = 1  # chr2
    r3.reference_start = 5200
    r3.mapping_quality = 60  # High quality
    # Normal RNA spliced alignment with reasonable exon sizes
    # 100M + 30N + 100M + 30N + 100M + 30N + 100M + 30N + 100M = 500 total
    r3.cigartuples = [(0, 100), (3, 30), (0, 100), (3, 30), (0, 100), (3, 30), (0, 100), (3, 30), (0, 100)]
    r3.query_qualities = pysam.qualitystring_to_array("I" * 500)
    fusion_reads.append(r3)

    return fusion_reads


def generate_bam_file(gtf_path: str, bam_path: str):
    """Generate a test BAM file with normal and fusion reads."""
    # Read GTF to get reference names
    print(f"Reading GTF: {gtf_path}")

    # Create header
    header = create_bam_header()

    # Write BAM file
    with pysam.AlignmentFile(bam_path, 'wb', header=header) as out_bam:
        # Generate normal reads
        print("Generating normal reads...")
        normal_reads = generate_normal_reads(out_bam, "read")
        for read in normal_reads:
            out_bam.write(read)

        # Generate fusion-like reads
        print("Generating fusion-like reads...")
        fusion_reads = generate_fusion_reads(out_bam, "read")
        for read in fusion_reads:
            out_bam.write(read)

        print(f"Generated {len(normal_reads)} normal reads and {len(fusion_reads)} fusion-like reads")

    # Index the BAM file
    print(f"Indexing BAM file: {bam_path}")
    pysam.index(bam_path)
    print(f"Created BAM index: {bam_path}.bai")


def main():
    """Main function to generate example data."""
    # Create output directory
    example_dir = Path("/Users/logan/Projects/pyfin/examples/test_data")
    example_dir.mkdir(parents=True, exist_ok=True)

    # File paths
    gtf_path = example_dir / "test_annotation.gtf"
    bam_path = example_dir / "test_reads.bam"

    # Generate files
    print("=== Generating Example Test Data ===\n")
    generate_gtf(str(gtf_path))
    generate_bam_file(str(gtf_path), str(bam_path))

    print("\n=== Summary ===")
    print(f"GTF file: {gtf_path}")
    print(f"  - 3 genes (gene_A, gene_B, gene_C)")
    print(f"  - Multiple transcripts per gene")
    print(f"BAM file: {bam_path}")
    print(f"  - 4 normal reads (1 per transcript region)")
    print(f"  - 3 fusion-like reads for testing fusion detection")
    print(f"  - 6 total reads")
    print("\nYou can test the region_separator with these files!")
    print(f"\nExample usage:")
    print(f"from fin.io.region_separator import separate_regions")
    print(f'separate_regions("{bam_path}", "{gtf_path}", "output/regions")')


if __name__ == "__main__":
    main()
