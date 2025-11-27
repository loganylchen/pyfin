#!/usr/bin/env python3
"""
Demonstration of sequence extraction from BED and BAM files.

This script shows how to extract sequences from:
1. Reference FASTA files using BED coordinates
2. BAM alignment files
"""

from fin.io.sequence_extractor import (
    ReferenceExtractor,
    BamSequenceExtractor,
    GenomicRegion,
    extract_sequence_from_bed_and_fasta,
    get_sequence_at_position
)

def demo_reference_extraction():
    """Demonstrate reference sequence extraction."""
    print("=" * 60)
    print("1. REFERENCE EXTRACTION FROM FASTA")
    print("=" * 60)

    # Example with reference FASTA file
    fasta_file = "reference.fasta"  # Replace with your file

    try:
        with ReferenceExtractor(fasta_file) as extractor:
            # Extract specific region
            region = GenomicRegion(
                chrom="chr1",
                start=10000,
                end=10100,
                name="test_region",
                strand="+"
            )

            seq = extractor.extract_sequence(region)
            print(f"Region: {region}")
            print(f"Sequence length: {len(seq)} bp")
            print(f"Sequence (first 50 bp): {seq[:50]}...")
            print()

            # Extract with flanking context
            seq_with_context = extractor.extract_region_with_context(
                chrom="chr1",
                pos=10500,
                flank=20
            )
            print(f"Region at chr1:10500 with ±20 bp flank")
            print(f"Sequence: {seq_with_context[:50]}...")
            print()

            # Get contig info
            contigs = extractor.get_contigs()
            print(f"Available contigs: {contigs[:5]}...")  # First 5
            print()

    except FileNotFoundError:
        print(f"File not found: {fasta_file}")
        print("Please provide a valid FASTA file\n")


def demo_bed_extraction():
    """Demonstrate BED file parsing with reference."""
    print("=" * 60)
    print("2. BED FILE SEQUENCE EXTRACTION")
    print("=" * 60)

    bed_file = "regions.bed"  # Replace with your file
    fasta_file = "reference.fasta"  # Replace with your file

    # Example BED format:
    # chr1\t1000\t1100\tregion1\t100\t+
    # chr1\t2000\t2100\tregion2\t100\t-

    try:
        # Extract all sequences from BED
        sequences = extract_sequence_from_bed_and_fasta(bed_file, fasta_file)

        print(f"Extracted {len(sequences)} sequences from BED file")
        for name, seq in list(sequences.items())[:3]:  # Show first 3
            print(f"  {name}: {len(seq)} bp")
        print()

        # Extract with context around a specific position
        seq, region = get_sequence_at_position(
            chrom="chr1",
            pos=5000,
            fasta_file=fasta_file,
            flank=30
        )
        print(f"Sequence at chr1:5000 ±30 bp:")
        print(f"  Region: {region}")
        print(f"  Sequence: {seq}")
        print()

    except FileNotFoundError as e:
        print(f"File not found: {e}")
        print("Please provide valid BED and FASTA files\n")


def demo_bam_extraction():
    """Demonstrate BAM file sequence extraction."""
    print("=" * 60)
    print("3. BAM FILE SEQUENCE EXTRACTION")
    print("=" * 60)

    bam_file = "alignments.bam"  # Replace with your file

    try:
        with BamSequenceExtractor(bam_file) as extractor:
            # Extract all reads in a region
            print("Extracting reads from chr1:10000-11000...")
            count = 0
            for read_id, seq, is_reverse in extractor.extract_sequences_in_region(
                chrom="chr1",
                start=10000,
                end=11000,
                min_mapq=20
            ):
                if count < 3:  # Show first 3 reads
                    strand = "-" if is_reverse else "+"
                    print(f"  Read: {read_id[:30]}...")
                    print(f"  Strand: {strand}, Length: {len(seq)} bp")
                    print(f"  Sequence: {seq[:50]}...")
                    print()
                count += 1

            print(f"Total reads extracted: {count}")
            print()

            # Extract soft-clipped sequences
            print("Extracting soft-clipped sequences...")
            clip_count = 0
            for read_id, clipped_seq, position in extractor.extract_soft_clipped_sequences(
                min_clip_len=10
            ):
                if clip_count < 3:  # Show first 3
                    print(f"  Read: {read_id[:20]}...")
                    print(f"  Position: {position}, Length: {len(clipped_seq)} bp")
                    print(f"  Clip: {clipped_seq}")
                    print()
                clip_count += 1

            print(f"Total soft clips extracted: {clip_count}")
            print()

    except FileNotFoundError:
        print(f"File not found: {bam_file}")
        print("Please provide a valid BAM file\n")


def demo_transcript_extraction():
    """Demonstrate spliced transcript extraction."""
    print("=" * 60)
    print("4. SPLICED TRANSCRIPT EXTRACTION")
    print("=" * 60)

    fasta_file = "reference.fasta"  # Replace with your file

    try:
        with ReferenceExtractor(fasta_file) as extractor:
            # Define exons for a transcript
            # These would typically come from a GTF/GFF file
            exons = [
                (1000, 1100),   # exon 1
                (1200, 1350),   # exon 2
                (1500, 1600),   # exon 3
            ]

            transcript_seq = extractor.extract_transcript_sequence(
                chrom="chr1",
                exons=exons,
                strand="+"
            )

            print(f"Transcript from {len(exons)} exons:")
            for i, (s, e) in enumerate(exons, 1):
                print(f"  Exon {i}: {s}-{e} ({e-s} bp)")
            print(f"Total transcript length: {len(transcript_seq)} bp")
            print(f"Sequence: {transcript_seq[:100]}...")
            print()

    except FileNotFoundError:
        print(f"File not found: {fasta_file}")
        print("Please provide a valid FASTA file\n")


def demo_rna_coord_transform():
    """Demonstrate RNA coordinate transformations."""
    print("=" * 60)
    print("5. RNA COORDINATE TRANSFORMATIONS")
    print("=" * 60)

    from fin.utils.sequences import reverse_complement

    # Short example sequence
    seq = "ATCGATCGATCG"

    print(f"Original sequence: {seq}")
    print(f"Reverse complement: {reverse_complement(seq)}")
    print()

    print("Note: For RNA mapping:")
    print("- RNA is sequenced 5'→3'")
    print("- Nanopore signals are generated 3'→5'")
    print("- When read is forward strand in BAM:")
    print("  - Signal is reverse complement of reference")
    print("  - Need to reverse complement to match signal")
    print("- When read is reverse strand in BAM:")
    print("  - Signal already matches reference orientation")
    print()


def main():
    """Run all demonstrations."""
    print("\n" + "=" * 60)
    print("PYBALEEN SEQUENCE EXTRACTION DEMONSTRATIONS")
    print("=" * 60 + "\n")

    demo_rna_coord_transform()
    demo_reference_extraction()
    demo_bed_extraction()
    demo_bam_extraction()
    demo_transcript_extraction()

    print("=" * 60)
    print("For more details, see the documentation:")
    print("- ReferenceExtractor: Extract from FASTA")
    print("- BamSequenceExtractor: Extract from BAM")
    print("- IOManager: Unified interface for both")
    print("=" * 60)


if __name__ == "__main__":
    main()
