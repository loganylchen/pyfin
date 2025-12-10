"""
Region separator for splitting BAM files by gene regions.

Separates reads from a BAM file into isolated gene regions based on GTF annotation,
particularly useful for RNA modification detection workflows. Also includes functionality
to detect and filter out fusion-like reads.
"""

import os
import pysam
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Iterator
import logging
from dataclasses import dataclass

from .io_bam import BamParser, ReadAlignment
from .io_gtf import GtfParser, TranscriptFeature

logger = logging.getLogger(__name__)


@dataclass
class GeneRegion:
    """Represents a gene region with its reads."""
    gene_id: str
    gene_name: Optional[str]
    chrom: str
    start: int
    end: int
    strand: str
    transcripts: List[TranscriptFeature]
    reads: List[ReadAlignment]
    read_count: int
    is_valid: bool


class FusionReadDetector:
    """
    Detects fusion-like reads that may indicate structural variants or
    alignment artifacts.
    """

    def __init__(
        self,
        min_alignment_blocks: int = 2,
        max_block_gap: int = 100000,
        min_block_size: int = 50,
        max_trimmed_bases: int = 20
    ):
        """
        Initialize fusion read detector.

        Args:
            min_alignment_blocks: Minimum number of alignment blocks to
                consider as potential fusion (default: 2)
            max_block_gap: Maximum allowed gap between alignment blocks in
                the same read (default: 100kb)
            min_block_size: Minimum size for an alignment block to be valid
            max_trimmed_bases: Maximum number of bases that can be trimmed at
                ends and still be considered a valid block
        """
        self.min_alignment_blocks = min_alignment_blocks
        self.max_block_gap = max_block_gap
        self.min_block_size = min_block_size
        self.max_trimmed_bases = max_trimmed_bases

    def is_fusion_like(self, alignment: pysam.AlignedSegment) -> bool:
        """
        Check if a read is fusion-like.

        Args:
            alignment: pysam alignment object

        Returns:
            True if read is fusion-like, False otherwise
        """
        # Unmapped reads or supplementary alignments are candidates
        if alignment.is_unmapped or alignment.is_supplementary:
            return True

        # Check for multiple alignment blocks (potential fusions)
        if hasattr(alignment, 'get_blocks'):
            blocks = alignment.get_blocks()
            if blocks:
                # Check number of blocks
                if len(blocks) >= self.min_alignment_blocks:
                    # Check gaps between blocks
                    for i in range(len(blocks) - 1):
                        prev_end = blocks[i][1]
                        next_start = blocks[i + 1][0]
                        gap = next_start - prev_end

                        # If gap is too large, this could be a fusion
                        if gap > self.max_block_gap:
                            return True

                # Check if any block is too small to be reliable
                for block_start, block_end in blocks:
                    block_size = block_end - block_start
                    if block_size < self.min_block_size:
                        return True

        # Check for abnormal clipping patterns
        cigartuples = alignment.cigartuples
        if cigartuples:
            total_length = sum(length for op, length in cigartuples)

            # Check soft clipping at both ends (chimeric reads often have this)
            left_clip = cigartuples[0][1] if cigartuples[0][0] in [4, 5] else 0
            right_clip = cigartuples[-1][1] if cigartuples[-1][0] in [4, 5] else 0

            # If more than max_trimmed_bases are clipped on both ends
            if left_clip > self.max_trimmed_bases and right_clip > self.max_trimmed_bases:
                if left_clip + right_clip > total_length * 0.5:  # More than 50% clipped
                    return True

        # Check mapping quality (fusion reads often have low mapq)
        if alignment.mapping_quality < 10:
            # But only flag as fusion if it has other suspicious features
            if hasattr(alignment, 'get_blocks') and alignment.get_blocks():
                if len(alignment.get_blocks()) >= 2:
                    return True

        return False


class RegionSeparator:
    """
    Separates reads from a BAM file into isolated gene regions based on GTF annotation.

    This class provides methods to split a BAM file into smaller BAM files containing
    only reads that overlap specific gene regions. It also filters out fusion-like reads
    that may indicate structural variants or alignment artifacts.
    """

    def __init__(
        self,
        bam_file: str,
        gtf_file: str,
        output_dir: str,
        detect_fusions: bool = True,
        fusion_detector: Optional[FusionReadDetector] = None
    ):
        """
        Initialize region separator.

        Args:
            bam_file: Path to input BAM file (must be indexed)
            gtf_file: Path to GTF annotation file
            output_dir: Directory to store separated BAM files
            detect_fusions: Whether to detect and filter fusion-like reads
            fusion_detector: Optional custom FusionReadDetector instance
        """
        self.bam_file = bam_file
        self.gtf_file = gtf_file
        self.output_dir = Path(output_dir)
        self.detect_fusions = detect_fusions

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Initialize parsers
        self.bam_parser = BamParser(self.bam_file)
        self.gtf_parser = GtfParser(self.gtf_file)

        # Initialize fusion detector
        if detect_fusions:
            self.fusion_detector = fusion_detector or FusionReadDetector()
        else:
            self.fusion_detector = None

        # Validate input files
        self._validate_inputs()

        # Statistics
        self.stats = {
            'total_reads': 0,
            'fusion_reads': 0,
            'reads_per_gene': {},
            'genes_with_reads': 0,
            'total_genes': 0
        }

    def _validate_inputs(self):
        """Validate that input files exist and are properly formatted."""
        # Check BAM file
        if not Path(self.bam_file).exists():
            raise FileNotFoundError(f"BAM file not found: {self.bam_file}")

        # Check BAM index
        bai_file = f"{self.bam_file}.bai"
        if not Path(bai_file).exists():
            raise FileNotFoundError(f"BAM index not found: {bai_file}. Please index with 'samtools index'")

        # Check GTF file
        if not Path(self.gtf_file).exists():
            raise FileNotFoundError(f"GTF file not found: {self.gtf_file}")

        logger.info(f"Validated inputs: BAM={self.bam_file}, GTF={self.gtf_file}")

    def extract_gene_regions(self, min_reads: int = 1) -> Iterator[GeneRegion]:
        """
        Extract reads for each gene region defined in GTF.

        Args:
            min_reads: Minimum number of reads required to keep a gene region

        Yields:
            GeneRegion objects containing reads for each gene
        """
        logger.info("Extracting gene regions...")

        # Get all genes from GTF
        genes = self.gtf_parser.get_all_genes()
        self.stats['total_genes'] = len(genes)

        logger.info(f"Found {len(genes)} genes in GTF")

        # Open BAM file with pysam for direct access
        with pysam.AlignmentFile(self.bam_file, 'rb') as bam:
            for gene_id, transcripts in genes.items():
                if not transcripts:
                    continue

                # Get gene boundaries
                gene_name = transcripts[0].attrs.get('gene_name')
                chrom = transcripts[0].chrom
                strand = transcripts[0].strand

                # Calculate gene boundaries
                gene_start = min(t.start for t in transcripts)
                gene_end = max(t.end for t in transcripts)

                # Add padding to capture reads that extend beyond annotated regions
                pad_start = max(0, gene_start - 1000)
                pad_end = gene_end + 1000

                # Fetch all reads in the region
                reads = []
                fusion_count = 0

                try:
                    for record in bam.fetch(chrom, pad_start, pad_end):
                        self.stats['total_reads'] += 1

                        # Check if this is a fusion-like read
                        if self.detect_fusions and self.fusion_detector.is_fusion_like(record):
                            self.stats['fusion_reads'] += 1
                            fusion_count += 1
                            continue

                        # Check if read actually overlaps the gene region
                        if record.reference_end <= gene_start or record.reference_start >= gene_end:
                            continue

                        # Convert to ReadAlignment
                        read_alignment = self._pysam_to_read_alignment(record)
                        reads.append(read_alignment)

                except ValueError:
                    # Chromosome not in BAM
                    logger.warning(f"Chromosome {chrom} not found in BAM, skipping gene {gene_id}")
                    continue

                if len(reads) >= min_reads:
                    self.stats['reads_per_gene'][gene_id] = len(reads)
                    self.stats['genes_with_reads'] += 1

                    yield GeneRegion(
                        gene_id=gene_id,
                        gene_name=gene_name,
                        chrom=chrom,
                        start=gene_start,
                        end=gene_end,
                        strand=strand,
                        transcripts=transcripts,
                        reads=reads,
                        read_count=len(reads),
                        is_valid=True
                    )
                else:
                    logger.debug(f"Gene {gene_id} has only {len(reads)} reads (< {min_reads}), skipping")

    def write_region_bams(self, min_reads: int = 1, prefix: str = "gene"):
        """
        Write extracted reads to separate BAM files for each gene region.

        Args:
            min_reads: Minimum number of reads required to write a region BAM
            prefix: Prefix for output BAM filenames

        Returns:
            List of paths to created BAM files
        """
        logger.info("Writing region BAM files...")

        created_files = []

        # Open input BAM to copy header
        with pysam.AlignmentFile(self.bam_file, 'rb') as input_bam:
            header = input_bam.header

            # Process each gene region
            for region in self.extract_gene_regions(min_reads=min_reads):
                # Create output filename
                safe_gene_id = region.gene_id.replace('/', '_').replace(':', '_')
                if region.gene_name:
                    filename = f"{prefix}_{safe_gene_id}_{region.gene_name}.bam"
                else:
                    filename = f"{prefix}_{safe_gene_id}.bam"

                output_path = self.output_dir / filename

                # Write reads to BAM
                with pysam.AlignmentFile(str(output_path), 'wb', header=header) as output_bam:
                    # Need to fetch reads again for BAM writing
                    with pysam.AlignmentFile(self.bam_file, 'rb') as bam:
                        try:
                            pad_start = max(0, region.start - 1000)
                            pad_end = region.end + 1000

                            for record in bam.fetch(region.chrom, pad_start, pad_end):
                                # Skip fusion reads
                                if self.detect_fusions and self.fusion_detector.is_fusion_like(record):
                                    continue

                                # Check overlap
                                if record.reference_end <= region.start or record.reference_start >= region.end:
                                    continue

                                # Write the record
                                output_bam.write(record)
                        except ValueError:
                            continue

                # Index the output BAM
                pysam.index(str(output_path))

                created_files.append(str(output_path))
                logger.info(f"Created: {output_path} ({region.read_count} reads)")

        logger.info(f"Created {len(created_files)} region BAM files")
        return created_files

    def write_region_list(self, output_file: str, min_reads: int = 1):
        """
        Write a list of regions with read counts to a file.

        Args:
            output_file: Path to output file
            min_reads: Minimum reads threshold
        """
        with open(output_file, 'w') as f:
            f.write("gene_id\tgene_name\tchrom\tstart\tend\tstrand\tread_count\ttranscript_count\n")

            for region in self.extract_gene_regions(min_reads=min_reads):
                f.write(f"{region.gene_id}\t{region.gene_name or 'NA'}\t"
                        f"{region.chrom}\t{region.start}\t{region.end}\t"
                        f"{region.strand}\t{region.read_count}\t"
                        f"{len(region.transcripts)}\n")

        logger.info(f"Region list written to: {output_file}")

    def get_statistics(self) -> Dict:
        """Return statistics about the separation process."""
        return self.stats.copy()

    def print_statistics(self):
        """Print statistics summary."""
        stats = self.get_statistics()

        print("========== Region Separation Statistics ==========")
        print(f"Total genes in GTF: {stats['total_genes']}")
        print(f"Genes with reads (≥1): {stats['genes_with_reads']}")
        print(f"Total reads processed: {stats['total_reads']}")
        print(f"Fusion-like reads filtered: {stats['fusion_reads']}")

        if stats['genes_with_reads'] > 0:
            avg_reads = sum(stats['reads_per_gene'].values()) / stats['genes_with_reads']
            print(f"Average reads per gene: {avg_reads:.1f}")

            # Top 5 genes by read count
            sorted_genes = sorted(stats['reads_per_gene'].items(),
                                  key=lambda x: x[1], reverse=True)[:5]
            print("\nTop 5 genes by read count:")
            for gene_id, read_count in sorted_genes:
                print(f"  {gene_id}: {read_count} reads")

        print("=" * 50)

    def _pysam_to_read_alignment(self, record: pysam.AlignedSegment) -> ReadAlignment:
        """Convert pysam AlignedSegment to ReadAlignment object."""
        alignment_score = None
        if record.has_tag('AS'):
            alignment_score = record.get_tag('AS')

        query_qualities = None
        if record.query_qualities is not None:
            query_qualities = list(record.query_qualities)

        return ReadAlignment(
            read_id=record.query_name,
            contig=record.reference_name,
            ref_start=record.reference_start,
            ref_end=record.reference_end,
            ref_length=record.reference_length,
            query_length=record.query_length,
            mapq=record.mapping_quality,
            is_reverse=record.is_reverse,
            is_secondary=record.is_secondary,
            is_supplementary=record.is_supplementary,
            alignment_score=alignment_score,
            query_sequence=record.query_sequence,
            query_qualities=query_qualities,
            aligned_pairs=None
        )


def separate_regions(
    bam_file: str,
    gtf_file: str,
    output_dir: str,
    min_reads: int = 1,
    detect_fusions: bool = True,
    write_bams: bool = True,
    region_list: Optional[str] = None,
    prefix: str = "gene"
) -> List[str]:
    """
    Convenience function to separate BAM file by gene regions.

    Args:
        bam_file: Path to input BAM file
        gtf_file: Path to GTF annotation file
        output_dir: Output directory for region files
        min_reads: Minimum reads required per region
        detect_fusions: Whether to filter fusion-like reads
        write_bams: Whether to write individual BAM files
        region_list: Optional path to write region list TSV
        prefix: Prefix for output BAM filenames

    Returns:
        List of paths to created BAM files (if write_bams=True)
    """
    separator = RegionSeparator(
        bam_file=bam_file,
        gtf_file=gtf_file,
        output_dir=output_dir,
        detect_fusions=detect_fusions
    )

    # Write region list if requested
    if region_list:
        separator.write_region_list(region_list, min_reads=min_reads)

    # Write BAM files
    created_files = []
    if write_bams:
        created_files = separator.write_region_bams(min_reads=min_reads, prefix=prefix)

    # Print statistics
    separator.print_statistics()

    return created_files
