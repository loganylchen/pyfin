"""
Interval Manager for generating isolated genomic intervals with strand separation

This module provides functions to generate non-overlapping genomic intervals
by analyzing both GTF annotation and BAM read alignments. Fusion reads are
identified and skipped during interval generation.

Key features:
- Strand-separated intervals (no mixing of + and -)
- No gene/transcript IDs stored in intervals
- Read counts tracked per interval
- Only essential info: chrom, start, end, strand, read_count
"""

from typing import List, Dict, Tuple, Optional, Iterator, Set, Any, Union
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass

import pysam

from .io_bam import BamReader
from .io_gtf import GTFReader
from ..utils.log_config import get_package_logger

logger = get_package_logger(__name__)


@dataclass
class GenomicInterval:
    """Non-overlapping genomic interval with strand separation"""
    chrom: str
    start: int  # 0-based
    end: int
    strand: Optional[str] = None  # '+' or '-', None for unstranded
    read_count: int = 0  # Number of reads mapped to this interval
    interval_id: Optional[str] = None

    @property
    def interval_tuple(self) -> Tuple[str, int, int]:
        """Get interval as (chrom, start, end) tuple"""
        return (self.chrom, self.start, self.end)

    @property
    def region_string(self) -> str:
        """Get region string for BAM/GTF extraction"""
        return f"{self.chrom}:{self.start}-{self.end}"

    def overlaps(self, other: 'GenomicInterval') -> bool:
        """Check if this interval overlaps with another"""
        return (self.chrom == other.chrom and
                self.start < other.end and
                self.end > other.start)

    def merge(self, other: 'GenomicInterval') -> 'GenomicInterval':
        """Merge with another interval (must overlap AND same strand)"""
        if self.strand != other.strand:
            raise ValueError("Cannot merge intervals on different strands")

        return GenomicInterval(
            chrom=self.chrom,
            start=min(self.start, other.start),
            end=max(self.end, other.end),
            strand=self.strand,
            read_count=self.read_count + other.read_count
        )


def intervals_to_bed(
    intervals: List[GenomicInterval],
    output_path: Union[str, Path, None] = None) -> Optional[str]:
    """
    Export GenomicInterval objects to BED format.
    
    Args:
        intervals: List of GenomicInterval objects
        output_path: Path to output file, file object, or None to return string
    Returns:
        BED formatted string if output_path is None, otherwise None
    """
    if not intervals:
        if output_path is None:
            return ""
        return None
    
    # Sort intervals by chrom, then start position
    bed_lines = []
    
    # Add UCSC BED header if requested
    
    for i, interval in enumerate(intervals):
        # Required BED3 columns: chrom, start, end
        chrom = interval.chrom
        start = interval.start  # Already 0-based
        end = interval.end      # Already exclusive
        
        # Optional columns
        name = interval.interval_id if interval.interval_id else f"interval_{i}"
        
      
        score = interval.read_count if interval.read_count else 0
        
        # Strand handling
        strand = interval.strand if interval.strand else "."
        
        # Build BED fields
        bed_fields = [chrom, str(start), str(end),name, str(score), strand]
    
        
        # Join fields with tabs (BED format requirement)
        bed_line = "\t".join(bed_fields)
        bed_lines.append(bed_line)
    
    bed_content = "\n".join(bed_lines)
    
    # Output to file, file object, or return string
    if output_path is None:
        return bed_content
    elif hasattr(output_path, 'write'):
        # output_path is a file-like object
        output_path.write(bed_content)
        return None
    else:
        # output_path is a string/Path
        with open(output_path, 'w') as f:
            f.write(bed_content)
        return None


def is_fusion_read(read_dict: Dict) -> bool:
    """
    Check if a read is a fusion candidate based on soft-clip length

    Args:
        read_dict: Dictionary representation of the read from BamReader

    Returns:
        True if read is a fusion candidate
    """
    SOFT_CLIP_THRESHOLD = 50

    if not read_dict.get('cigartuples'):
        return False

    # Check for long soft-clips at read ends
    cigartuples = read_dict['cigartuples']
    if not cigartuples:
        return False

    # Check first and last CIGAR operations for soft-clips
    first_op, first_len = cigartuples[0]
    last_op, last_len = cigartuples[-1]

    # Check if soft-clips exceed threshold at either end
    if first_op == 4 and first_len >= SOFT_CLIP_THRESHOLD:
        return True
    if last_op == 4 and last_len >= SOFT_CLIP_THRESHOLD:
        return True

    return False


def extract_strand_from_read(read_dict: Dict) -> Optional[str]:
    """
    Extract strand information from a read alignment

    For some BAM files, strand may not be directly available.
    This function can be extended based on specific needs.

    Args:
        read_dict: Dictionary representation of the read

    Returns:
        Strand ('+' or '-') or None if not available
    """
    # Placeholder - can be extended based on specific BAM format
    # Some BAM files store strand in flags or tags
    return '+' if read_dict.get('is_forward') else '-'


def generate_intervals_from_reads(bam_path: str, max_reads: Optional[int] = None) -> Tuple[List[Dict], Set[str]]:
    """
    Generate intervals from BAM read alignments

    Args:
        bam_path: Path to BAM file
        max_reads: Optional limit on number of reads to process

    Returns:
        Tuple of (read_alignments, fusion_read_ids)
    """
    logger.info("Processing BAM reads for interval generation...")

    fusion_read_ids = set()
    read_alignments = []

    with BamReader(bam_path) as bam_reader:
        bam_reader.open()

        read_count = 0
        for alignment in bam_reader._alignment_file.fetch(until_eof=True):
            if max_reads is not None and read_count >= max_reads:
                break

            read_dict = bam_reader.alignment_to_dict(alignment)

            # Check for fusion reads
            if is_fusion_read(read_dict):
                fusion_read_ids.add(read_dict['query_name'])
                read_count += 1
                if read_count % 10000 == 0:
                    logger.info(f"Processed {read_count} reads...")
                continue

            # Keep the read alignment info
            read_alignments.append(read_dict)
            read_count += 1

            if read_count % 10000 == 0:
                logger.info(f"Processed {read_count} reads...")

    logger.info(f"Collected {len(read_alignments)} read alignments, {len(fusion_read_ids)} fusion reads")
    return read_alignments, fusion_read_ids


def cluster_intervals(read_alignments: List[Dict], max_gap: int = 0) -> List[GenomicInterval]:
    """
    Cluster overlapping or nearby intervals by chromosome and strand

    Args:
        read_alignments: List of read alignment dictionaries
        max_gap: Maximum gap between intervals to merge (default: 0bp)

    Returns:
        List of clustered GenomicInterval objects (strand-separated)
    """
    if not read_alignments:
        return []

    logger.info(f"Clustering {len(read_alignments)} read alignments with max_gap={max_gap}...")

    # Group by chromosome and strand
    chrom_strand_intervals = defaultdict(list)
    for read_dict in read_alignments:
        chrom = read_dict.get('reference_name')
        strand = extract_strand_from_read(read_dict) or None
        start = read_dict.get('reference_start')
        end = read_dict.get('reference_end')
        read_id = read_dict.get('query_name')

        if chrom and start is not None and end is not None:
            key = (chrom, strand)  # Separate by strand
            chrom_strand_intervals[key].append((start, end, read_id))

    clustered = []
    interval_counter = 0

    for (chrom, strand), intervals in sorted(chrom_strand_intervals.items()):
        if not intervals:
            continue

        # Count reads per position for proper read counting
        position_counts = defaultdict(set)  # position -> set of read_ids
        for start, end, read_id in intervals:
            # Mark all positions in the interval as covered by this read
            for pos in range(start, end):
                position_counts[pos].add(read_id)

        # Sort intervals by start position
        intervals_sorted = sorted(intervals, key=lambda x: x[0])

        current_cluster = None
        current_reads = set()

        for start, end, read_id in intervals_sorted:
            if current_cluster is None:
                # Start new cluster
                current_cluster = (start, end)
                current_reads = {read_id}
            elif start <= current_cluster[1] + max_gap:
                # Merge with current cluster
                current_cluster = (
                    current_cluster[0],
                    max(current_cluster[1], end)
                )
                current_reads.add(read_id)
            else:
                # Save current cluster and start new one
                clustered.append(GenomicInterval(
                    chrom=chrom,
                    start=current_cluster[0],
                    end=current_cluster[1],
                    strand=strand,
                    read_count=len(current_reads),
                    interval_id=f"interval_{interval_counter:06d}"
                ))
                interval_counter += 1
                current_cluster = (start, end)
                current_reads = {read_id}

        # Save the last cluster
        if current_cluster:
            clustered.append(GenomicInterval(
                chrom=chrom,
                start=current_cluster[0],
                end=current_cluster[1],
                strand=strand,
                read_count=len(current_reads),
                interval_id=f"interval_{interval_counter:06d}"
            ))
            interval_counter += 1

    logger.info(f"Created {len(clustered)} clustered intervals (strand-separated)")
    return clustered


def generate_intervals_from_gtf(gtf_path: str) -> List[GenomicInterval]:
    """
    Generate intervals from GTF annotation (for merging with read intervals)

    Args:
        gtf_path: Path to GTF file

    Returns:
        List of genomic intervals from transcripts
    """
    logger.info("Generating intervals from GTF annotation...")

    with GTFReader(gtf_path) as gtf_reader:
        gtf_reader.open()
        gtf_reader.parse()

        intervals = []
        for transcript in gtf_reader.iterate_transcripts():
            interval = GenomicInterval(
                chrom=transcript.chrom,
                start=transcript.start,
                end=transcript.end,
                strand=transcript.strand,
                interval_id=f"gtf_{transcript.transcript_id}"
            )
            intervals.append(interval)

    logger.info(f"Generated {len(intervals)} intervals from GTF")
    return intervals


def merge_gtf_and_read_intervals(gtf_intervals: List[GenomicInterval],
                                read_intervals: List[GenomicInterval],
                                max_gap: int = 50) -> List[GenomicInterval]:
    """
    Merge GTF annotation intervals with read coverage intervals

    Args:
        gtf_intervals: Intervals from GTF
        read_intervals: Intervals from BAM reads (clustered)
        max_gap: Maximum gap for merging intervals

    Returns:
        List of merged, non-overlapping GenomicInterval objects (strand-separated)
    """
    logger.info("Merging GTF intervals with read coverage...")

    # Combine all intervals
    all_intervals = gtf_intervals + read_intervals

    # Group by chromosome and strand
    chrom_strand_intervals = defaultdict(list)
    for interval in all_intervals:
        key = (interval.chrom, interval.strand)
        chrom_strand_intervals[key].append(interval)

    # Merge intervals on each chromosome and strand
    merged = []

    # Sort with proper handling of None strand values
    # Convert None to a sortable value (use 'z' so it sorts after '+' and '-')
    for (chrom, strand), intervals in sorted(
        chrom_strand_intervals.items(),
        key=lambda x: (x[0][0], 'z' if x[0][1] is None else x[0][1])
    ):
        if not intervals:
            continue

        # Sort by start position
        intervals_sorted = sorted(intervals, key=lambda x: x.start)

        merged_chrom = []
        current = intervals_sorted[0]

        for interval in intervals_sorted[1:]:
            # Check if intervals should be merged (same strand check already done)
            if (current.overlaps(interval) or interval.start <= current.end + max_gap):
                # Merge intervals
                try:
                    current = current.merge(interval)
                except ValueError:
                    # Different strand (shouldn't happen with our grouping)
                    merged_chrom.append(current)
                    current = interval
            else:
                # No overlap, save current
                merged_chrom.append(current)
                current = interval

        # Save the last interval
        merged_chrom.append(current)
        merged.extend(merged_chrom)

    logger.info(f"Generated {len(merged)} final intervals (strand-separated)")
    return merged


def generate_isolated_intervals(bam_path: str,
                                gtf_path: Optional[str] = None,
                                max_gap: int =0,
                                max_reads: Optional[int] = None, tmp_dir: Optional[str] = None) -> Dict[str, Any]:
    """
    Generate isolated, non-overlapping genomic intervals with strand separation

    This is the main function that generates intervals by analyzing both
    GTF annotation and BAM read alignments. Fusion reads are identified
    and skipped during interval generation.

    Args:
        bam_path: Path to BAM file (required)
        gtf_path: Path to GTF file (optional, improves annotation)
        max_gap: Maximum gap between intervals to merge (default: 0)
        max_reads: Optional limit on number of reads to process

    Returns:
        Dictionary with:
        - 'intervals': List of GenomicInterval objects (strand-separated)
        - 'fusion_read_ids': Set of fusion read IDs
        - 'num_reads_processed': Number of reads processed
    """
    logger.info("=" * 60)
    logger.info("Starting interval generation")
    logger.info(f"BAM: {bam_path}")
    logger.info(f"GTF: {gtf_path}")
    logger.info(f"Max gap: {max_gap}bp")
    logger.info("=" * 60)

    # Process BAM reads
    read_alignments, fusion_read_ids = generate_intervals_from_reads(bam_path, max_reads)

    # Cluster read alignments
    read_intervals = cluster_intervals(read_alignments, max_gap)
    # logger.debug(read_intervals)
    if tmp_dir:
        Path(tmp_dir).mkdir(exist_ok=True)
        tmp_read_intervals = Path(tmp_dir) / "read_intervals.bed"
        # tmp_read_intervals.mkdir(exist_ok=True)
        intervals_to_bed(read_intervals,tmp_read_intervals)
    # Generate intervals from GTF if provided
    gtf_intervals = []
    if gtf_path:
        gtf_intervals = generate_intervals_from_gtf(gtf_path)
    if tmp_dir:
        tmp_gtf_intervals = Path(tmp_dir) / "gtf_intervals.bed"
        # tmp_gtf_intervals.mkdir(exist_ok=True)
        intervals_to_bed(gtf_intervals,tmp_gtf_intervals)
    # Merge intervals
    final_intervals = merge_gtf_and_read_intervals(
        gtf_intervals,
        read_intervals,
        max_gap=max_gap
    )
    # logger.debug(final_intervals)
    # Sort intervals by chromosome, strand, position
    strand_order = {'+': 0, '-': 1, None: 2}
    final_intervals.sort(key=lambda x: (x.chrom, strand_order.get(x.strand, 2), x.start))
    if tmp_dir:
        tmp_final_intervals = Path(tmp_dir) / "final_intervals.bed"
        # tmp_final_intervals.mkdir(exist_ok=True)
        intervals_to_bed(final_intervals,tmp_final_intervals)
    # Log summary
    logger.info("=" * 60)
    logger.info("Interval generation completed")
    logger.info(f"Total intervals: {len(final_intervals)}")
    logger.info(f"Fusion reads detected: {len(fusion_read_ids)}")

    # Log intervals by strand
    strand_counts = defaultdict(int)
    for interval in final_intervals:
        strand_counts[interval.strand] += 1

    
    logger.info("Intervals by strand:")
    logger.debug(strand_counts.items())
    for strand, count in sorted(strand_counts.items()):
        logger.info(f"  {strand if strand else 'unstranded'}: {count}")

    logger.info("=" * 60)

    # Log first 5 intervals
    logger.info("First 5 intervals:")
    for i, interval in enumerate(final_intervals[:5], 1):
        logger.info(f"  {i}. {interval.region_string} [strand: {interval.strand}] "
                   f"(reads: {interval.read_count})")

    if len(final_intervals) > 5:
        logger.info(f"  ... and {len(final_intervals) - 5} more intervals")

    return {
        'intervals': final_intervals,
        'fusion_read_ids': fusion_read_ids,
        'num_reads_processed': len(read_alignments) + len(fusion_read_ids)
    }


def extract_reads_for_interval(bam_path: str,
                               interval: GenomicInterval,
                               fusion_read_ids: Optional[Set[str]] = None) -> List[Dict]:
    """
    Extract reads for a specific interval from BAM file

    Args:
        bam_path: Path to BAM file
        interval: GenomicInterval to extract reads for
        fusion_read_ids: Optional set of fusion read IDs to skip

    Returns:
        List of read dictionaries
    """
    if fusion_read_ids is None:
        fusion_read_ids = set()

    reads = []

    with BamReader(bam_path) as bam_reader:
        bam_reader.open()

        # Fetch reads from the interval
        for alignment in bam_reader._alignment_file.fetch(
            interval.chrom,
            interval.start,
            interval.end
        ):
            read_dict = bam_reader.alignment_to_dict(alignment)
            read_id = read_dict['query_name']

            # Skip fusion reads
            if read_id in fusion_read_ids or is_fusion_read(read_dict):
                continue

            # Check if read actually overlaps
            read_start = read_dict.get('reference_start')
            read_end = read_dict.get('reference_end')

            if (read_start is not None and read_end is not None and
                read_start < interval.end and read_end > interval.start):
                reads.append(read_dict)

    logger.debug(f"Extracted {len(reads)} reads for interval {interval.region_string}")
    return reads


def extract_annotation_for_interval(gtf_path: str,
                                    interval: GenomicInterval) -> List[Dict]:
    """
    Extract annotation features for a specific interval from GTF file

    Args:
        gtf_path: Path to GTF file
        interval: GenomicInterval to extract annotation for

    Returns:
        List of annotation features
    """
    features = []

    with GTFReader(gtf_path) as gtf_reader:
        gtf_reader.open()
        gtf_reader.parse()

        for transcript in gtf_reader.iterate_transcripts():
            # Check if transcript overlaps interval
            if (transcript.chrom == interval.chrom and
                transcript.start < interval.end and transcript.end > interval.start):
                features.append({
                    'transcript_id': transcript.transcript_id,
                    'gene_id': getattr(transcript, 'gene_id', None),
                    'start': transcript.start,
                    'end': transcript.end,
                    'strand': transcript.strand
                })

    logger.debug(f"Extracted {len(features)} features for interval {interval.region_string}")
    return features
