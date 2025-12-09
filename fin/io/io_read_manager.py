"""
Read Subset Manager for BAM/GTF/FASTA integration

This module provides classes for organizing nanopore Direct RNA-seq reads into logical subsets
for efficient processing: transcript-specific reads, fusion candidates, and unannotated reads.
Reads are categorized without creating intermediate files until explicitly requested.
"""

from typing import List, Dict, Optional, Iterator, Tuple, Set, Any, Generator
from pathlib import Path
from collections import defaultdict
import pysam
import tempfile
import os

from .io_bam import BamReader
from .io_gtf import GTFReader

# Import package logger
from ..utils.log_config import get_package_logger

logger = get_package_logger(__name__)


class ReadSubset:
    """Represents a subset of reads for a transcriptomic interval"""

    def __init__(self, subset_id: str, chrom: str, start: int, end: int,
                 strand: Optional[str] = None, transcript_ids: Optional[List[str]] = None,
                 read_ids: Optional[Set[str]] = None, is_fusion: bool = False,
                 fusion_breakpoints: Optional[List[Tuple[int, int]]] = None):
        """
        Initialize read subset

        Args:
            subset_id: Unique subset identifier
            chrom: Chromosome name
            start: Start position (0-based)
            end: End position
            strand: Strand information (+, -, or None)
            transcript_ids: List of transcript IDs this subset represents
            read_ids: Set of read IDs in this subset
            is_fusion: Whether this is a fusion candidate subset
            fusion_breakpoints: List of fusion breakpoint coordinates
        """
        self.subset_id = subset_id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.transcript_ids = transcript_ids or []
        self.read_ids = read_ids or set()
        self.is_fusion = is_fusion
        self.fusion_breakpoints = fusion_breakpoints or []
        self.bam_reader = None  # Will be set by manager

    @property
    def interval(self) -> Tuple[str, int, int]:
        """Get interval as (chrom, start, end) tuple"""
        return (self.chrom, self.start, self.end)

    @property
    def num_reads(self) -> int:
        """Get number of reads in this subset"""
        return len(self.read_ids)

    @property
    def region_string(self) -> str:
        """Get region string for BAM fetching"""
        return f"{self.chrom}:{self.start}-{self.end}"

    def __str__(self) -> str:
        if self.is_fusion:
            tx_info = f"({len(self.fusion_breakpoints)} breakpoints)"
        else:
            tx_info = f"({len(self.transcript_ids)} transcripts)"
        return (f"ReadSubset({self.subset_id}, {self.region_string}, "
                f"{self.num_reads} reads, {tx_info}")


class ReadDataBundle:
    """
    Bundle containing reads and associated data for processing

    This class holds the actual read data, sequence, and annotation information
    for a subset of reads, ready for downstream analysis.
    """

    def __init__(self, subset: ReadSubset, reads: List[Dict[str, Any]],
                 reference_sequences: Optional[Dict[str, str]] = None):
        """
        Initialize read data bundle

        Args:
            subset: The ReadSubset this bundle corresponds to
            reads: List of read dictionaries from BamReader
            reference_sequences: Dictionary of reference sequences for the region
        """
        self.subset = subset
        self.reads = reads
        self.reference_sequences = reference_sequences or {}

    @property
    def subset_id(self) -> str:
        """Get subset ID"""
        return self.subset.subset_id

    @property
    def num_reads(self) -> int:
        """Get number of reads"""
        return len(self.reads)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary representation"""
        return {
            'subset_id': self.subset_id,
            'interval': self.subset.interval,
            'num_reads': self.num_reads,
            'is_fusion': self.subset.is_fusion,
            'transcript_ids': self.subset.transcript_ids,
            'fusion_breakpoints': self.subset.fusion_breakpoints,
            'reads': self.reads,
        }


class ReadSubsetManager:
    """
    Manager for organizing reads into structured subsets

    This class analyzes a BAM file in the context of GTF annotation to group reads
    by transcriptomic intervals. It provides iteration over subsets and can generate
    actual BAM files for each subset when needed.

    Features:
    - Groups reads by transcriptomic intervals from GTF
    - Identifies fusion candidate reads (multi-mapping, chimeric)
    - Separates unannotated reads into their own subset
    - Allows lazy generation of subset BAM files
    - Ensures reads are unique to one subset (except fusion candidates)
    """

    def __init__(self, bam_path: str, gtf_path: str, fasta_path: Optional[str] = None):
        """
        Initialize manager

        Args:
            bam_path: Path to BAM file
            gtf_path: Path to GTF annotation file
            fasta_path: Optional path to reference FASTA file
        """
        self.bam_path = Path(bam_path)
        self.gtf_path = Path(gtf_path)
        self.fasta_path = Path(fasta_path) if fasta_path else None

        # Readers
        self.bam_reader = BamReader(str(self.bam_path))
        self.gtf_reader = GTFReader(str(self.gtf_path))

        # Data structures
        self.subsets = []  # List of all ReadSubset objects
        self.fusion_subset = None  # Special subset for fusion candidates
        self.unannotated_subset = None  # Special subset for unannotated reads
        self.read_id_to_subset = {}  # Map read IDs to subset IDs

        # Caches
        self._fasta_sequences = {}  # Cache for FASTA sequences
        self._reads_by_id = {}  # Cache for read objects

        logger.info(f"Initialized ReadSubsetManager with BAM: {bam_path}, "
                   f"GTF: {gtf_path}, FASTA: {fasta_path}")

    def __enter__(self):
        """Context manager entry"""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()

    def open(self):
        """Open all files"""
        self.bam_reader.open()
        self.gtf_reader.open()
        self.gtf_reader.parse()  # Parse GTF to build transcript structure

        # Load FASTA if available
        if self.fasta_path and self.fasta_path.exists():
            self._load_fasta()

    def close(self):
        """Close all files"""
        self.bam_reader.close()
        self.gtf_reader.close()

    def _load_fasta(self):
        """Load reference sequences from FASTA"""
        try:
            from .io_fasta import FASTAReader
            with FASTAReader(str(self.fasta_path)) as fasta_reader:
                for record in fasta_reader.iterate_records():
                    self._fasta_sequences[record.id] = record.sequence
            logger.info(f"Loaded {len(self._fasta_sequences)} reference sequences")
        except Exception as e:
            logger.warning(f"Failed to load FASTA file: {e}")

    def analyze_reads(self, max_reads: Optional[int] = None):
        """
        Analyze all reads and categorize them into subsets

        This method iterates through all reads in the BAM file and groups them
        based on their alignment to annotated transcripts or identifies them
        as fusion candidates/unannotated reads.

        Args:
            max_reads: Optional limit on number of reads to process
        """
        logger.info("Starting read analysis and categorization...")

        # Initialize special subsets
        self.fusion_subset = ReadSubset(
            subset_id="FUSION_CANDIDATES",
            chrom="multi",
            start=0,
            end=0,
            is_fusion=True
        )
        self.fusion_subset.bam_reader = self.bam_reader

        self.unannotated_subset = ReadSubset(
            subset_id="UNANNOTATED",
            chrom="unannotated",
            start=0,
            end=0
        )
        self.unannotated_subset.bam_reader = self.bam_reader

        read_count = 0
        # First pass: collect all read alignments and their intervals
        read_intervals = defaultdict(list)  # read_id -> list of intervals

        logger.info("First pass: collecting read intervals from BAM...")
        for alignment in self.bam_reader._alignment_file.fetch(until_eof=True):
            if max_reads is not None and read_count >= max_reads:
                break

            read_dict = self.bam_reader.alignment_to_dict(alignment)
            read_id = read_dict['query_name']

            # Store read in cache
            self._reads_by_id[read_id] = alignment

            chrom = read_dict.get('reference_name')
            start = read_dict.get('reference_start', 0)
            end = read_dict.get('reference_end', 0)

            if chrom and start and end:
                read_intervals[read_id].append((chrom, start, end))

            read_count += 1
            if read_count % 10000 == 0:
                logger.info(f"Processed {read_count} reads in first pass...")

        # Build transcript subsets based on BOTH GTF annotation and read coverage
        self._build_transcript_subsets_with_reads(read_intervals)

        # Second pass: categorize reads
        logger.info("Second pass: categorizing reads...")
        categorized_count = 0
        for alignment in self.bam_reader._alignment_file.fetch(until_eof=True):
            if max_reads is not None and categorized_count >= max_reads:
                break

            read_dict = self.bam_reader.alignment_to_dict(alignment)
            read_id = read_dict['query_name']

            # Skip if read_id not in initial pass (shouldn't happen)
            if read_id not in read_intervals:
                continue

            # Check for fusion candidates first
            if self._is_fusion_candidate(alignment, read_dict):
                self.fusion_subset.read_ids.add(read_id)
                self.read_id_to_subset[read_id] = self.fusion_subset.subset_id
                categorized_count += 1
                continue

            # Try to assign to transcript-based subset
            assigned = self._assign_to_transcript_subset(alignment, read_dict)

            if not assigned:
                # Assign to unannotated subset
                self.unannotated_subset.read_ids.add(read_id)
                self.read_id_to_subset[read_id] = self.unannotated_subset.subset_id

            categorized_count += 1

            if categorized_count % 10000 == 0:
                logger.info(f"Categorized {categorized_count} reads...")

        # Filter out empty subsets and add special subsets
        self.subsets = [s for s in self.subsets if s.num_reads > 0]
        self.subsets.append(self.fusion_subset)
        self.subsets.append(self.unannotated_subset)

        logger.info(f"Completed analysis: {len(self.subsets)} subsets, "
                   f"{categorized_count} total reads categorized")

        # Sort subsets by chromosome and position
        self.subsets.sort(key=lambda s: (s.chrom, s.start))

    def _build_transcript_subsets_with_reads(self, read_intervals: Dict[str, List[Tuple[str, int, int]]]):
        """
        Build transcript-based subsets considering both GTF annotation and BAM read coverage

        This ensures subsets are created where there are actually reads aligning, not just where
        there are annotated transcripts.

        Args:
            read_intervals: Dictionary mapping read_id to list of (chrom, start, end) intervals
        """
        logger.info("Building transcript subsets based on annotation and read coverage...")

        # Group reads by chromosome and position (clustering)
        chrom_clusters = defaultdict(list)

        for read_id, intervals in read_intervals.items():
            for chrom, start, end in intervals:
                chrom_clusters[chrom].append((start, end, read_id))

        # Process each chromosome separately
        subset_counter = 0
        for chrom in sorted(chrom_clusters.keys()):
            if chrom == 'unmapped':
                continue

            # Sort intervals by start position
            intervals = sorted(chrom_clusters[chrom], key=lambda x: x[0])

            # Merge overlapping intervals to create clusters
            # Using hierarchical clustering approach
            clusters = []
            current_cluster = []

            for start, end, read_id in intervals:
                if not current_cluster:
                    current_cluster = [(start, end, read_id)]
                else:
                    # Check if this interval overlaps or is close to current cluster
                    cluster_start = min(s for s, e, r in current_cluster)
                    cluster_end = max(e for s, e, r in current_cluster)

                    # Merge if there's significant overlap or they're close (< 1000bp)
                    if start <= cluster_end + 1000:
                        current_cluster.append((start, end, read_id))
                    else:
                        # Save current cluster and start new one
                        clusters.append(current_cluster)
                        current_cluster = [(start, end, read_id)]

            if current_cluster:
                clusters.append(current_cluster)

            # For each cluster, check if it overlaps with transcripts
            logger.info(f"Processing {len(clusters)} clusters on {chrom}...")

            for cluster in clusters:
                cluster_start = min(s for s, e, r in cluster)
                cluster_end = max(e for s, e, r in cluster)

                # Find transcripts that overlap this cluster
                overlapping_transcripts = []
                for transcript in self.gtf_reader.iterate_transcripts():
                    if transcript.chrom == chrom and \
                       transcript.start < cluster_end and transcript.end > cluster_start:
                        overlapping_transcripts.append(transcript)

                if overlapping_transcripts:
                    # This cluster has annotation - create transcript subset
                    merged_start = min(cluster_start, min(t.start for t in overlapping_transcripts))
                    merged_end = max(cluster_end, max(t.end for t in overlapping_transcripts))

                    # Determine strand (use majority strand from transcripts)
                    strands = [t.strand for t in overlapping_transcripts if t.strand]
                    strand = max(set(strands), key=strands.count) if strands else None

                    subset_id = f"TX_{subset_counter:06d}_{chrom}_{cluster_start}_{cluster_end}"
                    subset = ReadSubset(
                        subset_id=subset_id,
                        chrom=chrom,
                        start=merged_start,
                        end=merged_end,
                        strand=strand,
                        transcript_ids=[t.transcript_id for t in overlapping_transcripts]
                    )
                    subset.bam_reader = self.bam_reader
                    self.subsets.append(subset)
                    subset_counter += 1

                # Note: clusters without transcript overlap will be handled as unannotated

        logger.info(f"Created {subset_counter} transcript-based subsets")

    def _build_transcript_subsets(self):
        """Build transcript-based subsets from GTF"""
        subset_counter = 0

        for transcript in self.gtf_reader.iterate_transcripts():
            # Group transcripts that overlap heavily into the same subset
            # to avoid having too many tiny subsets
            merged = False
            for subset in self.subsets:
                if self._should_merge_transcripts(transcript, subset):
                    subset.transcript_ids.append(transcript.transcript_id)
                    # Update interval to cover both
                    subset.start = min(subset.start, transcript.start)
                    subset.end = max(subset.end, transcript.end)
                    merged = True
                    break

            if not merged:
                subset_id = f"TX_{subset_counter:06d}_{transcript.chrom}_{transcript.start}_{transcript.end}"
                subset = ReadSubset(
                    subset_id=subset_id,
                    chrom=transcript.chrom,
                    start=transcript.start,
                    end=transcript.end,
                    strand=transcript.strand,
                    transcript_ids=[transcript.transcript_id]
                )
                subset.bam_reader = self.bam_reader
                self.subsets.append(subset)
                subset_counter += 1

    def _should_merge_transcripts(self, transcript: Any, subset: ReadSubset) -> bool:
        """
        Determine if a transcript should be merged into an existing subset

        Args:
            transcript: GTF transcript
            subset: Existing subset

        Returns:
            True if should merge
        """
        if transcript.chrom != subset.chrom:
            return False

        # Don't merge across strands
        if subset.strand and transcript.strand and subset.strand != transcript.strand:
            return False

        # Check if overlap is significant (e.g., > 50% of transcript)
        overlap_start = max(subset.start, transcript.start)
        overlap_end = min(subset.end, transcript.end)

        if overlap_start < overlap_end:
            transcript_length = transcript.end - transcript.start
            overlap_length = overlap_end - overlap_start

            # Merge if overlap > 50% or gap < 5kb
            gap = max(0, subset.start - transcript.end, transcript.start - subset.end)
            return (overlap_length > transcript_length * 0.5 or
                    (gap < 5000 and overlap_length > 0))

        return False

    def _is_fusion_candidate(self, alignment: pysam.AlignedSegment, read_dict: Dict[str, Any]) -> bool:
        """
        Check if a read is a fusion candidate based on soft-clip length

        Fusion reads typically have long soft-clipped sequences at their ends that didn't
        align to the reference, indicating they span a fusion junction.

        Args:
            alignment: Pysam alignment object
            read_dict: Dictionary representation of the read

        Returns:
            True if fusion candidate
        """
        # Default threshold for soft-clip length (can be made configurable)
        SOFT_CLIP_THRESHOLD = 50  # bases

        if not read_dict.get('cigartuples'):
            return False

        # Check for long soft-clips at read ends
        cigartuples = read_dict['cigartuples']

        # Soft-clip operations: 4 (soft clipping)
        # Check first and last CIGAR operations for soft-clips
        first_op, first_len = cigartuples[0]
        last_op, last_len = cigartuples[-1]

        # Check if soft-clips exceed threshold at either end
        if first_op == 4 and first_len >= SOFT_CLIP_THRESHOLD:
            return True
        if last_op == 4 and last_len >= SOFT_CLIP_THRESHOLD:
            return True

        return False

    def _assign_to_transcript_subset(self, alignment: pysam.AlignedSegment,
                                     read_dict: Dict[str, Any]) -> bool:
        """
        Assign read to appropriate transcript-based subset

        Args:
            alignment: Pysam alignment object
            read_dict: Dictionary representation of the read

        Returns:
            True if assigned, False if unannotated
        """
        chrom = read_dict.get('reference_name')
        if not chrom:
            return False

        start = read_dict.get('reference_start', 0)
        end = read_dict.get('reference_end', 0)

        if not start or not end:
            return False

        # Find overlapping transcript subsets
        overlapping_subsets = []
        for subset in self.subsets:
            if subset.chrom == chrom and subset.start < end and subset.end > start:
                # Calculate overlap
                overlap_start = max(subset.start, start)
                overlap_end = min(subset.end, end)
                overlap_length = overlap_end - overlap_start

                if overlap_length > 0:
                    overlapping_subsets.append((subset, overlap_length))

        if not overlapping_subsets:
            return False

        # Assign to subset with maximum overlap
        best_subset, _ = max(overlapping_subsets, key=lambda x: x[1])
        best_subset.read_ids.add(read_dict['query_name'])
        self.read_id_to_subset[read_dict['query_name']] = best_subset.subset_id

        return True

    def iterate_subsets(self) -> Generator[ReadSubset, None, None]:
        """
        Iterator over all read subsets

        Yields:
            ReadSubset objects
        """
        yield from self.subsets

    def get_subset(self, subset_id: str) -> Optional[ReadSubset]:
        """
        Get a specific subset by ID

        Args:
            subset_id: Subset ID

        Returns:
            ReadSubset object or None
        """
        for subset in self.subsets:
            if subset.subset_id == subset_id:
                return subset
        return None

    def get_data_bundle(self, subset: ReadSubset) -> ReadDataBundle:
        """
        Get full data bundle with reads and sequences for a subset

        Args:
            subset: ReadSubset object

        Returns:
            ReadDataBundle object
        """
        # Fetch reads for this subset
        reads = []
        for read_id in subset.read_ids:
            if read_id in self._reads_by_id:
                read_dict = self.bam_reader.alignment_to_dict(self._reads_by_id[read_id])
                reads.append(read_dict)

        # Get reference sequences if FASTA is available
        reference_sequences = {}
        if self.fasta_path and subset.chrom in self._fasta_sequences:
            seq_start = max(0, subset.start - 1000)  # Add some flanking sequence
            seq_end = min(len(self._fasta_sequences[subset.chrom]), subset.end + 1000)
            reference_sequences[subset.chrom] = self._fasta_sequences[subset.chrom][seq_start:seq_end]

        return ReadDataBundle(subset, reads, reference_sequences)

    def write_subset_bam(self, subset: ReadSubset, output_path: Optional[str] = None) -> str:
        """
        Write reads from a subset to a new BAM file

        Args:
            subset: ReadSubset object
            output_path: Output BAM file path (auto-generated if None)

        Returns:
            Path to created BAM file
        """
        if not output_path:
            output_path = f"{subset.subset_id}.bam"

        logger.info(f"Writing {subset.num_reads} reads to {output_path}")

        try:
            # Create output file with same header
            with pysam.AlignmentFile(str(self.bam_path), 'rb') as input_bam:
                with pysam.AlignmentFile(str(output_path), 'wb', template=input_bam) as output_bam:
                    for read_id in subset.read_ids:
                        if read_id in self._reads_by_id:
                            output_bam.write(self._reads_by_id[read_id])

            # Index the output BAM
            pysam.index(str(output_path))

            logger.info(f"Successfully created BAM file: {output_path}")
            return str(output_path)

        except Exception as e:
            logger.error(f"Failed to write subset BAM: {e}")
            raise

    def get_summary(self) -> Dict[str, Any]:
        """
        Get summary statistics

        Returns:
            Dictionary with summary statistics
        """
        return {
            'num_subsets': len(self.subsets),
            'num_fusion_reads': self.fusion_subset.num_reads if self.fusion_subset else 0,
            'num_unannotated_reads': self.unannotated_subset.num_reads if self.unannotated_subset else 0,
            'subset_sizes': [s.num_reads for s in self.subsets],
            'subsets': [{
                'id': s.subset_id,
                'chrom': s.chrom,
                'start': s.start,
                'end': s.end,
                'num_reads': s.num_reads,
                'num_transcripts': len(s.transcript_ids),
                'is_fusion': s.is_fusion
            } for s in self.subsets]
        }

    def clear_read_cache(self):
        """Clear the read cache to free memory"""
        self._reads_by_id.clear()


def create_subset_manager(bam_path: str, gtf_path: str, fasta_path: Optional[str] = None,
                          max_reads: Optional[int] = None) -> ReadSubsetManager:
    """
    Convenience function to create and initialize a ReadSubsetManager

    Args:
        bam_path: Path to BAM file
        gtf_path: Path to GTF annotation file
        fasta_path: Optional path to FASTA file
        max_reads: Optional limit on number of reads to process

    Returns:
        Initialized ReadSubsetManager
    """
    manager = ReadSubsetManager(bam_path, gtf_path, fasta_path)
    with manager:
        manager.analyze_reads(max_reads=max_reads)
    return manager
