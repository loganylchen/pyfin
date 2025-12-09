"""
Sequence extraction utilities for reference genomes and BAM files.

Provides functions to extract DNA/RNA sequences from:
- Reference FASTA files using BED coordinates (supports .gz)
- BAM alignment files
- Reference genomes with coordinate transformations
"""

import pysam
from typing import List, Dict, Optional, Iterator, Tuple, Union
import logging
import gzip
import tempfile
import os
from pathlib import Path
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class GenomicRegion:
    """Represents a genomic region."""
    chrom: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
    name: Optional[str] = None
    strand: str = "+"

    @property
    def length(self) -> int:
        """Region length in bases."""
        return self.end - self.start

    def __str__(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}({self.strand})"


class ReferenceExtractor:
    """
    Extract sequences from reference FASTA/FASTQ files.

    Handles coordinate transformations and strand-specific extraction.
    Supports both plain text and gzipped (.gz) FASTA files.
    """

    def __init__(self, reference_file: str):
        """
        Initialize reference extractor.

        Args:
            reference_file: Path to reference FASTA file (plain or .gz)
        """
        self.reference_file = reference_file
        self.fasta = None
        self._temp_file = None
        self._temp_index = None

        # Handle gzipped files
        actual_file = reference_file
        if reference_file.endswith('.gz'):
            logger.warning(
                f"Detected gzipped FASTA file: {reference_file}. "
                "Note: Using gzipped FASTA may require re-indexing on each use. "
                f"Consider using bgzip (pysam-compatible) for better performance."
            )
            actual_file = self._extract_gzip(reference_file)

        # Open reference file
        try:
            self.fasta = pysam.FastaFile(actual_file)
            logger.info(f"Opened reference file: {reference_file}")
            logger.info(f"Contigs: {list(self.fasta.references)}")
        except Exception as e:
            logger.error(f"Failed to open reference file {reference_file}: {e}")
            raise

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def _extract_gzip(self, gz_file: str) -> str:
        """
        Extract gzipped FASTA to a temporary file.

        Args:
            gz_file: Path to gzipped FASTA file

        Returns:
            Path to temporary extracted file
        """
        # Create temporary directory and file
        temp_dir = tempfile.mkdtemp()
        self._temp_file = os.path.join(temp_dir, Path(gz_file).stem)

        # Extract gzipped content
        with gzip.open(gz_file, 'rt') as f_in:
            with open(self._temp_file, 'w') as f_out:
                chunk_size = 8192
                while True:
                    chunk = f_in.read(chunk_size)
                    if not chunk:
                        break
                    f_out.write(chunk)

        logger.info(f"Extracted gzipped FASTA to temporary file: {self._temp_file}")
        return self._temp_file

    def close(self):
        """Close the reference file and cleanup temporary files."""
        if self.fasta is not None:
            self.fasta.close()
            self.fasta = None
            logger.debug(f"Closed reference file: {self.reference_file}")

        # Cleanup temporary extracted file if needed
        if self._temp_file and os.path.exists(self._temp_file):
            try:
                # Remove both the file and the temp directory
                temp_dir = os.path.dirname(self._temp_file)
                os.unlink(self._temp_file)
                os.rmdir(temp_dir)
                logger.debug(f"Removed temporary directory: {temp_dir}")
            except Exception as e:
                logger.warning(f"Failed to remove temporary file: {e}")
            self._temp_file = None

    def get_contigs(self) -> List[str]:
        """Get list of contig names."""
        return list(self.fasta.references)

    def get_contig_length(self, contig: str) -> int:
        """Get length of a specific contig."""
        return self.fasta.get_reference_length(contig)

    def extract_sequence(self, region: GenomicRegion) -> str:
        """
        Extract sequence from reference genome.

        Args:
            region: Genomic region to extract

        Returns:
            Sequence string (forward strand)
        """
        # Check bounds
        contig_length = self.get_contig_length(region.chrom)
        if region.start < 0 or region.end > contig_length:
            logger.error(f"Region {region} outside contig bounds (0-{contig_length})")
            raise ValueError("Region outside contig bounds")

        # Extract sequence (returns forward strand)
        try:
            seq = self.fasta.fetch(region.chrom, region.start, region.end)
        except Exception as e:
            logger.error(f"Failed to extract sequence for {region}: {e}")
            raise

        logger.debug(f"Extracted {len(seq)} bp from {region}")
        return seq

    def extract_sequence_bed(self, bed_line: str) -> Tuple[str, Optional[str]]:
        """
        Extract sequence from a BED line.

        Args:
            bed_line: BED format line (chrom\tstart\tend[\tname][\t...][\tstrand])

        Returns:
            Tuple of (sequence, name) where name may be None
        """
        fields = bed_line.strip().split('\t')

        if len(fields) < 3:
            raise ValueError(f"Invalid BED line (need at least 3 fields): {bed_line}")

        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        name = fields[3] if len(fields) > 3 else None
        strand = fields[5] if len(fields) > 5 else "+"

        region = GenomicRegion(chrom=chrom, start=start, end=end, name=name, strand=strand)
        seq = self.extract_sequence(region)

        logger.info(f"Extracted {len(seq)} bp for BED entry: {chrom}:{start}-{end}")
        return seq, name

    def extract_region_with_context(self, chrom: str, start: int, end: int,
                                    flank: int = 50, strand: str = "+") -> str:
        """
        Extract sequence with flanking context.

        Args:
            chrom: Chromosome/contig name
            start: Start position (0-based)
            end: End position (exclusive)
            flank: Number of flanking bases to include
            strand: Strand orientation

        Returns:
            Sequence string including flanking regions
        """
        contig_length = self.get_contig_length(chrom)

        # Calculate bounds with flanking
        start_with_flank = max(0, start - flank)
        end_with_flank = min(contig_length, end + flank)

        # Extract extended sequence
        region = GenomicRegion(
            chrom=chrom,
            start=start_with_flank,
            end=end_with_flank,
            strand=strand
        )
        seq = self.extract_sequence(region)

        logger.debug(f"Extracted {len(seq)} bp (with {flank} bp flank) for {chrom}:{start}-{end}")
        return seq

    def extract_regions_from_bed(self, bed_file: str) -> Iterator[Tuple[str, GenomicRegion]]:
        """
        Extract sequences from all regions in a BED file.

        Args:
            bed_file: Path to BED file

        Yields:
            Tuple of (sequence_string, GenomicRegion)
        """
        with open(bed_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue

                try:
                    seq, name = self.extract_sequence_bed(line)

                    fields = line.strip().split('\t')
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    strand = fields[5] if len(fields) > 5 else "+"

                    region = GenomicRegion(
                        chrom=chrom, start=start, end=end,
                        name=name, strand=strand
                    )

                    yield seq, region

                except Exception as e:
                    logger.error(f"Error on line {line_num} of {bed_file}: {e}")
                    raise

    def extract_transcript_sequence(self, chrom: str, exons: List[Tuple[int, int]],
                                   strand: str = "+") -> str:
        """
        Extract spliced transcript sequence from exon coordinates.

        Args:
            chrom: Chromosome name
            exons: List of (start, end) tuples for exons
            strand: Strand orientation

        Returns:
            Spliced transcript sequence
        """
        transcript_seq = ""

        for start, end in exons:
            region = GenomicRegion(chrom=chrom, start=start, end=end, strand=strand)
            exon_seq = self.extract_sequence(region)
            transcript_seq += exon_seq

        logger.info(f"Extracted transcript sequence: {len(transcript_seq)} bp from {len(exons)} exons")
        return transcript_seq


class BamSequenceExtractor:
    """
    Extract sequences and information from BAM alignment files.

    Handles both basecalled sequences from BAM and reference coordinates.
    """

    def __init__(self, bam_file: str):
        """
        Initialize BAM sequence extractor.

        Args:
            bam_file: Path to BAM file
        """
        self.bam_file = bam_file
        self.bam = None

        try:
            self.bam = pysam.AlignmentFile(bam_file, 'rb')
            logger.info(f"Opened BAM file: {bam_file}")
        except Exception as e:
            logger.error(f"Failed to open BAM file {bam_file}: {e}")
            raise

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def close(self):
        """Close the BAM file."""
        if self.bam is not None:
            self.bam.close()
            self.bam = None
            logger.debug(f"Closed BAM file: {self.bam_file}")

    def get_sequence_from_read(self, read_id: str) -> Optional[str]:
        """
        Extract basecalled sequence for a specific read.

        Args:
            read_id: Read identifier

        Returns:
            Basecalled sequence string or None if not found
        """
        # Note: BAM files don't have random access by read_id
        # Need to iterate through file or use index
        logger.warning("get_sequence_from_read requires iterating through BAM - slow for large files")

        for record in self.bam:
            if record.query_name == read_id:
                if record.query_sequence:
                    seq = record.query_sequence
                    logger.debug(f"Found sequence for read {read_id}: {len(seq)} bp")
                    return seq
                else:
                    logger.debug(f"Read {read_id} has no sequence")
                    return None

        logger.warning(f"Read {read_id} not found in BAM")
        return None

    def extract_sequences_in_region(self, chrom: str, start: int, end: int,
                                   min_mapq: int = 20) -> Iterator[Tuple[str, str, bool]]:
        """
        Extract sequences from all reads overlapping a region.

        Args:
            chrom: Chromosome/contig name
            start: Start position (0-based)
            end: End position (exclusive)
            min_mapq: Minimum mapping quality

        Yields:
            Tuple of (read_id, sequence, is_reverse)
        """
        try:
            for record in self.bam.fetch(chrom, start, end):
                if record.mapping_quality < min_mapq:
                    continue

                if record.is_unmapped:
                    continue

                if not record.query_sequence:
                    continue

                read_id = record.query_name
                seq = record.query_sequence
                is_reverse = record.is_reverse

                logger.debug(f"Extracted {len(seq)} bp from read {read_id}")
                yield read_id, seq, is_reverse

        except ValueError as e:
            logger.error(f"Error fetching region {chrom}:{start}-{end}: {e}")
            raise

    def get_aligned_pairs(self, read_id: str) -> Optional[List[Tuple[int, int]]]:
        """
        Get aligned pairs (reference_pos, query_pos) for a read.

        Args:
            read_id: Read identifier

        Returns:
            List of (ref_pos, query_pos) tuples or None
        """
        for record in self.bam:
            if record.query_name == read_id:
                if record.has_tag('MD'):
                    pairs = []
                    for qp, rp, _ in record.get_aligned_pairs():
                        if qp is not None and rp is not None:
                            pairs.append((rp, qp))
                    return pairs
                break

        return None

    def extract_soft_clipped_sequences(self, min_clip_len: int = 10) -> Iterator[Tuple[str, str, str]]:
        """
        Extract soft-clipped sequences from alignments.

        Args:
            min_clip_len: Minimum soft clip length to report

        Yields:
            Tuple of (read_id, clipped_seq, position: "start" or "end")
        """
        for record in self.bam:
            if record.is_unmapped:
                continue

            read_id = record.query_name
            cigar = record.cigar

            if not cigar:
                continue

            # Check for soft clips at beginning
            if cigar[0][0] == 4:  # 4 = soft clip
                clip_len = cigar[0][1]
                if clip_len >= min_clip_len:
                    clipped_seq = record.query_sequence[:clip_len]
                    yield read_id, clipped_seq, "start"

            # Check for soft clips at end
            if cigar[-1][0] == 4:  # 4 = soft clip
                clip_len = cigar[-1][1]
                if clip_len >= min_clip_len:
                    clipped_seq = record.query_sequence[-clip_len:]
                    yield read_id, clipped_seq, "end"


def reverse_complement(seq: str) -> str:
    """
    Compute reverse complement of DNA/RNA sequence.

    Args:
        seq: DNA/RNA sequence string

    Returns:
        Reverse complement string
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
                  'N': 'N', 'n': 'n'}

    return "".join(complement.get(base, base) for base in reversed(seq))


def extract_sequence_from_bed_and_fasta(bed_file: str, fasta_file: str,
                                        output_file: Optional[str] = None) -> Dict[str, str]:
    """
    Extract sequences from BED file using reference FASTA.

    Args:
        bed_file: Path to BED file with regions
        fasta_file: Path to reference FASTA file
        output_file: Optional output file path (FASTA format)

    Returns:
        Dictionary of region_name -> sequence
    """
    sequences = {}

    with ReferenceExtractor(fasta_file) as extractor:
        for seq, region in extractor.extract_regions_from_bed(bed_file):
            name = region.name or f"{region.chrom}:{region.start}-{region.end}"
            sequences[name] = seq

            logger.info(f"Extracted {len(seq)} bp for {name}")

    # Write to output file if requested
    if output_file:
        with open(output_file, 'w') as f:
            for name, seq in sequences.items():
                f.write(f">{name}\n{seq}\n")
        logger.info(f"Wrote {len(sequences)} sequences to {output_file}")

    return sequences


def get_sequence_at_position(chrom: str, pos: int, fasta_file: str,
                             flank: int = 50) -> Tuple[str, GenomicRegion]:
    """
    Get sequence at a specific genomic position with flanking context.

    Args:
        chrom: Chromosome name
        pos: Position (0-based)
        fasta_file: Path to reference FASTA
        flank: Number of flanking bases

    Returns:
        Tuple of (sequence, GenomicRegion)
    """
    with ReferenceExtractor(fasta_file) as extractor:
        start = pos - flank
        end = pos + flank + 1

        region = GenomicRegion(chrom=chrom, start=start, end=end, strand="+")
        seq = extractor.extract_sequence(region)

        logger.info(f"Extracted {len(seq)} bp at {chrom}:{pos} (±{flank} bp)")

        return seq, region
