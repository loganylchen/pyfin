"""
BAM Clustering Pipeline - Extract reads and transcripts by intervals and cluster by 3' end positions
"""

from typing import List, Dict, Set, Optional, Tuple, Any
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import pickle

from ..io.io_bam import BamReader
from ..io.io_fasta import FASTAReader
from ..io.io_gtf import GTFReader
from ..io.interval_manager import GenomicInterval, extract_reads_for_interval
from ..utils.log_config import get_package_logger

logger = get_package_logger(__name__)


@dataclass
class ClusteredData:
    """Data container for a single cluster"""
    cluster_id: str
    chrom: str
    strand: str
    three_prime_positions: List[int]
    read_sequences: Dict[str, str]
    read_qualities: Dict[str, List[int]]
    read_cigars: Dict[str, List[Tuple[int, str]]]
    transcript_sequences: Dict[str, str]
    transcript_ids: List[str]
    interval: GenomicInterval


@dataclass
class ThreePrimePosition:
    """Represents a 3' end position"""
    position: int
    is_read: bool
    read_id: Optional[str] = None
    transcript_id: Optional[str] = None
    cigar: Optional[List[Tuple[int, str]]] = None
    strand: str = '+'


class BamClusteringPipeline:
    """
    Pipeline for extracting reads and transcripts by intervals, clustering by 3' end positions,
    and preparing sequences for downstream analysis.
    """

    def __init__(self, bam_path: str, fasta_path: str, gtf_path: Optional[str] = None):
        """
        Initialize pipeline

        Args:
            bam_path: Path to BAM file
            fasta_path: Path to reference FASTA file
            gtf_path: Optional path to GTF annotation file
        """
        self.bam_path = bam_path
        self.fasta_path = fasta_path
        self.gtf_path = gtf_path
        self.fusion_read_ids: Set[str] = set()

    def process_intervals(self, intervals: List[GenomicInterval],
                         max_reads_per_interval: Optional[int] = None,
                         distance_threshold: int = 100,
                         fusion_read_ids: Optional[Set[str]] = None,
                         output_path: Optional[str] = None) -> List[ClusteredData]:
        """
        Process all intervals: extract reads, cluster by 3' ends, prepare sequences

        Args:
            intervals: List of GenomicInterval objects
            max_reads_per_interval: Maximum reads to extract per interval
            distance_threshold: Distance threshold for 3' end clustering
            fusion_read_ids: Set of fusion read IDs to skip
            output_path: Optional path to save results

        Returns:
            List of ClusteredData objects, one per cluster
        """
        if fusion_read_ids is None:
            fusion_read_ids = set()
        if self.fusion_read_ids:
            fusion_read_ids.update(self.fusion_read_ids)

        all_clusters = []

        logger.info(f"Processing {len(intervals)} intervals")

        for interval_idx, interval in enumerate(intervals):
            logger.info(f"Processing interval {interval_idx + 1}/{len(intervals)}: {interval.region_string}")

            # Extract reads for this interval
            reads = self.extract_reads(interval, max_reads_per_interval, fusion_read_ids)

            # Extract transcripts from GTF if available
            transcripts = self.extract_transcripts(interval)

            # Get 3' end positions for reads and transcripts
            three_prime_positions = self.get_three_prime_positions(reads, transcripts, interval)

            # Cluster by 3' end positions
            clusters = self.cluster_three_prime_positions(three_prime_positions, distance_threshold)

            # Extract sequences for each cluster
            for cluster_idx, cluster_positions in enumerate(clusters):
                cluster_data = self.prepare_cluster_data(
                    cluster_positions, reads, transcripts, interval,
                    cluster_id=f"{interval.interval_id}_cluster_{cluster_idx}"
                )
                all_clusters.append(cluster_data)

        logger.info(f"Generated {len(all_clusters)} total clusters from {len(intervals)} intervals")

        # Save results if output path provided
        if output_path:
            self.save_clusters(all_clusters, output_path)

        return all_clusters

    def extract_reads(self, interval: GenomicInterval,
                     max_reads: Optional[int] = None,
                     fusion_read_ids: Optional[Set[str]] = None) -> List[Dict[str, Any]]:
        """
        Extract reads for a specific interval

        Args:
            interval: GenomicInterval to extract reads for
            max_reads: Maximum number of reads to extract
            fusion_read_ids: Set of fusion read IDs to skip

        Returns:
            List of read dictionaries
        """
        logger.debug(f"Extracting reads for interval: {interval.region_string}")

        # Use existing function from interval_manager
        reads = extract_reads_for_interval(self.bam_path, interval, fusion_read_ids)

        # Apply max_reads limit if specified
        if max_reads and len(reads) > max_reads:
            reads = reads[:max_reads]
            logger.info(f"Limited reads to {max_reads} out of {len(reads)} total")

        logger.info(f"Extracted {len(reads)} reads for interval {interval.region_string}")
        return reads

    def extract_transcripts(self, interval: GenomicInterval) -> List[Any]:
        """
        Extract transcripts that overlap with the interval

        Args:
            interval: GenomicInterval to extract transcripts for

        Returns:
            List of GTFTranscript objects
        """
        if not self.gtf_path:
            logger.debug("No GTF path provided, skipping transcript extraction")
            return []

        logger.debug(f"Extracting transcripts for interval: {interval.region_string}")

        transcripts = []
        with GTFReader(self.gtf_path) as gtf_reader:
            gtf_reader.open()
            gtf_reader.parse()

            for transcript in gtf_reader.iterate_transcripts():
                # Check overlap with interval
                if (transcript.chrom == interval.chrom and
                    transcript.start < interval.end and
                    transcript.end > interval.start):
                    transcripts.append(transcript)

        logger.info(f"Extracted {len(transcripts)} transcripts for interval {interval.region_string}")
        return transcripts

    def get_three_prime_positions(self, reads: List[Dict], transcripts: List,
                                interval: GenomicInterval) -> List[ThreePrimePosition]:
        """
        Get 3' end positions for all reads and transcripts

        Args:
            reads: List of read dictionaries
            transcripts: List of GTFTranscript objects
            interval: Interval we're working on

        Returns:
            List of ThreePrimePosition objects
        """
        three_prime_positions = []

        # Process reads
        for read in reads:
            strand = read.get('is_forward', True)
            strand_char = '+' if strand else '-'

            # Get read sequence
            read_seq = read.get('query_sequence', '')
            if not read_seq:
                logger.debug(f"No sequence for read {read['query_name']}, skipping")
                continue

            # Calculate 3' end position based on strand
            # For forward (+) reads: 3' end is reference_end
            # For reverse (-) reads: 3' end is reference_start
            if strand:  # Forward strand
                three_prime_pos = read.get('reference_end', 0)
            else:  # Reverse strand
                three_prime_pos = read.get('reference_start', 0)

            three_prime_positions.append(ThreePrimePosition(
                position=three_prime_pos,
                is_read=True,
                read_id=read['query_name'],
                strand=strand_char,
                cigar=read.get('cigartuples')
            ))

        # Process transcripts
        for transcript in transcripts:
            # Calculate 3' end position based on strand
            # For + strand: 3' end is transcript end
            # For - strand: 3' end is transcript start
            if transcript.strand == '+':
                three_prime_pos = transcript.end
            else:
                three_prime_pos = transcript.start

            three_prime_positions.append(ThreePrimePosition(
                position=three_prime_pos,
                is_read=False,
                transcript_id=transcript.transcript_id,
                strand=transcript.strand
            ))

        logger.info(f"Collected {len(three_prime_positions)} 3' end positions (reads + transcripts)")
        return three_prime_positions

    def cluster_three_prime_positions(self,
                                    three_prime_positions: List[ThreePrimePosition],
                                    distance_threshold: int = 100) -> List[List[ThreePrimePosition]]:
        """
        Cluster 3' end positions by genomic proximity

        Args:
            three_prime_positions: List of ThreePrimePosition objects
            distance_threshold: Maximum distance to cluster positions together

        Returns:
            List of clusters, where each cluster is a list of ThreePrimePosition objects
        """
        if not three_prime_positions:
            return []

        logger.debug(f"Clustering {len(three_prime_positions)} 3' end positions with threshold {distance_threshold}")

        # Sort by position
        sorted_positions = sorted(three_prime_positions, key=lambda x: x.position)

        clusters = []
        current_cluster = [sorted_positions[0]]

        for i in range(1, len(sorted_positions)):
            curr_pos = sorted_positions[i]
            prev_pos = sorted_positions[i - 1]

            # Check if within distance threshold
            if curr_pos.position - prev_pos.position <= distance_threshold:
                # Add to current cluster
                current_cluster.append(curr_pos)
            else:
                # Close current cluster and start new one
                if current_cluster:
                    clusters.append(current_cluster)
                current_cluster = [curr_pos]

        # Add the last cluster
        if current_cluster:
            clusters.append(current_cluster)

        logger.info(f"Created {len(clusters)} clusters from {len(three_prime_positions)} positions")

        # Log cluster sizes
        cluster_sizes = [len(c) for c in clusters]
        logger.debug(f"Cluster sizes: min={min(cluster_sizes) if cluster_sizes else 0}, "
                    f"max={max(cluster_sizes) if cluster_sizes else 0}, "
                    f"mean={sum(cluster_sizes)/len(cluster_sizes) if cluster_sizes else 0:.1f}")

        return clusters

    def prepare_cluster_data(self,
                           cluster_positions: List[ThreePrimePosition],
                           reads: List[Dict],
                           transcripts: List,
                           interval: GenomicInterval,
                           cluster_id: str) -> ClusteredData:
        """
        Prepare all data for a single cluster

        Args:
            cluster_positions: List of ThreePrimePosition objects in this cluster
            reads: List of all read dictionaries for the interval
            transcripts: List of all GTFTranscript objects for the interval
            interval: The interval this cluster belongs to
            cluster_id: Unique ID for this cluster

        Returns:
            ClusteredData object
        """
        # Extract sequences for reads in this cluster
        read_sequences = {}
        read_qualities = {}
        read_cigars = {}
        transcript_sequences = {}
        transcript_ids = []

        chrom = interval.chrom

        # Find reads in this cluster
        cluster_read_ids = {pos.read_id for pos in cluster_positions if pos.is_read and pos.read_id}

        logger.debug(f"Preparing data for cluster {cluster_id} with {len(cluster_read_ids)} reads")

        # Extract read sequences
        for read in reads:
            read_id = read['query_name']
            if read_id in cluster_read_ids:
                # Get read sequence
                seq = read.get('query_sequence')
                if seq:
                    read_sequences[read_id] = seq

                # Get base qualities
                qualities = read.get('query_qualities')
                if qualities:
                    read_qualities[read_id] = qualities

                # Get CIGAR
                cigar = read.get('cigartuples')
                if cigar:
                    read_cigars[read_id] = cigar

        # Extract transcript sequences
        cluster_transcript_ids = {pos.transcript_id for pos in cluster_positions
                                 if not pos.is_read and pos.transcript_id}

        if self.gtf_path and cluster_transcript_ids and transcripts:
            with FASTAReader(self.fasta_path) as fasta_reader:
                for transcript in transcripts:
                    if transcript.transcript_id in cluster_transcript_ids:
                        transcript_ids.append(transcript.transcript_id)

                        # Get spliced sequence
                        chrom_seq = fasta_reader.get_sequence(chrom)
                        if chrom_seq:
                            transcript_seq = transcript.get_spliced_sequence(chrom_seq)
                            transcript_sequences[transcript.transcript_id] = transcript_seq

        # Collect all 3' end positions
        three_prime_positions = [pos.position for pos in cluster_positions]

        cluster_data = ClusteredData(
            cluster_id=cluster_id,
            chrom=chrom,
            strand=interval.strand or '+',
            three_prime_positions=three_prime_positions,
            read_sequences=read_sequences,
            read_qualities=read_qualities,
            read_cigars=read_cigars,
            transcript_sequences=transcript_sequences,
            transcript_ids=transcript_ids,
            interval=interval
        )

        logger.debug(f"Prepared cluster {cluster_id}: {len(read_sequences)} reads, "
                    f"{len(transcript_sequences)} transcripts")

        return cluster_data

    def save_clusters(self, clusters: List[ClusteredData], output_path: str):
        """
        Save clusters to file

        Args:
            clusters: List of ClusteredData objects
            output_path: Path to save results
        """
        output_path = Path(output_path)

        # Create output directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Pickle the clusters
        with open(output_path, 'wb') as f:
            pickle.dump(clusters, f)

        logger.info(f"Saved {len(clusters)} clusters to {output_path}")

        # Also save summary stats
        summary_path = output_path.with_suffix('.summary.txt')
        with open(summary_path, 'w') as f:
            f.write(f"BAM Clustering Pipeline Summary\n")
            f.write(f"================================\n\n")
            f.write(f"Total clusters: {len(clusters)}\n\n")
            f.write(f"Cluster details:\n")

            for cluster in clusters:
                f.write(f"  {cluster.cluster_id}:\n")
                f.write(f"    Reads: {len(cluster.read_sequences)}\n")
                f.write(f"    Transcripts: {len(cluster.transcript_sequences)}\n")
                f.write(f"    Strand: {cluster.strand}\n")
                f.write(f"    3' positions: {sorted(cluster.three_prime_positions)}\n\n")

        logger.info(f"Saved summary to {summary_path}")
