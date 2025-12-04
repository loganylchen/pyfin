"""
BAM Clustering Pipeline - Extract reads and transcripts by intervals and cluster by 3' end positions
"""

from typing import List, Dict, Set, Optional, Tuple, Any, Generator
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import pickle
import pysam
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
    sequences: List[str]
    ids: List[str]






class ThreePrimePositionClustering:
    """
    Pipeline for extracting reads and transcripts by intervals, clustering by 3' end positions,
    and preparing sequences for downstream analysis.
    """

    def __init__(self, bam_path: str, fasta_path: str):
        """
        Initialize pipeline

        Args:
            bam_path: Path to BAM file
            fasta_path: Path to reference FASTA file
        """
        self.bam_path = Path(bam_path)
        self.fasta_path = Path(fasta_path)
        self.open()
        

    def open(self):
        """Open the alignment file"""
        try:
            mode = 'rb' if self.bam_path.suffix.lower() == '.bam' else 'r'
            self._alignment_file = BamReader(self.bam_path)
            self._is_open = True
            
            
            self._reference_file = pysam.FastaFile(str(self.fasta_path))
            logger.info(f"Opened reference file: {self.fasta_path}")
            
        except Exception as e:
            logger.error(f"Failed to open alignment/reference file: {e}")
            raise

    def close(self):
        """Close the alignment file"""
        if self._alignment_file is not None:
            self._alignment_file.close()
            self._alignment_file = None
        if self._reference_file is not None:
            self._reference_file.close()
            self._reference_file = None
            logger.info(f"Closed alignment file: {self.fasta_path}")
        self._is_open = False
            



    def iter_interval(self, intervals: List[GenomicInterval],
                         max_reads_per_interval: Optional[int] = None,
                         min_reads_per_interval: int = 10, 
                         distance_threshold: int = 24) -> Generator[Tuple, None, None]:
        """
        Process all intervals: extract reads, cluster by 3' ends, prepare sequences

        Args:
            intervals: List of GenomicInterval objects
            max_reads_per_interval: Maximum reads to extract per interval
            distance_threshold: Distance threshold for 3' end clustering

        Yields:
            ClusteredData object
        """


        all_clusters = []

        logger.info(f"Processing {len(intervals)} intervals")

        for interval_idx, interval in enumerate(intervals):
            logger.info(f"Processing interval {interval_idx + 1}/{len(intervals)}: {interval.region_string}")
            if interval.read_count <= min_reads_per_interval:
                logger.debug(f"Skipped interval {interval_idx + 1} due to low read count {interval.read_count}")
                continue
            # Cluster by 3' end positions
            clusters = self.cluster_three_prime_positions(interval.attrs,interval.three_prime_pos)
            read_sequences = self.extract_read_sequence(interval,max_reads=max_reads_per_interval)
            reference_sequences = self.extract_ref_sequence(interval)
            yield (read_sequences,reference_sequences,clusters)
            # Extract sequences for each cluster
            # for cluster_idx, cluster_dict in enumerate(clusters):
            #     cluster_data = self.prepare_cluster_data(
            #         cluster_dict['ids'],
            #         cluster_dict['3_positions'],
            #         read_sequences,
            #         reference_sequences,
            #         interval.strand
            #     )
                

        # logger.info(f"Generated {len(all_clusters)} total clusters from {len(intervals)} intervals")

  

    def extract_read_sequence(self, interval: GenomicInterval,
                     max_reads: Optional[int] = None) -> Dict[str, str]:
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
        reads = [i for i in interval.attrs if not i.startswith('gtf_')]

        # Apply max_reads limit if specified
        if max_reads and len(reads) > max_reads:
            reads = reads[:max_reads]
            logger.info(f"Limited reads to {max_reads} out of {len(reads)} total")
        read_sequences = {read.query_name:self._alignment_file.get_reference_sequence(read) for read in self._alignment_file.fetch(region=interval.region_string) if self._alignment_file.get_reference_sequence(read)}
        logger.info(f"Extracted {len(reads)} reads for interval {interval.region_string}")
        return read_sequences

    def extract_ref_sequence(self, interval: GenomicInterval) -> Dict[str, str]:
        """
        Extract transcripts that overlap with the interval

        Args:
            interval: GenomicInterval to extract transcripts for

        Returns:
            List of GTFTranscript objects
        """
        transcripts = [i.replace('gtf_','') for i in interval.attrs if i.startswith('gtf_')]
        transcripts_sequence = {transcript:self._reference_file.fetch(transcript) for transcript in transcripts}
        logger.info(f"Extracted {len(transcripts)} transcripts for interval {interval.region_string}")
        return transcripts_sequence


    def cluster_three_prime_positions(self,
                                    ids: List[str],
                                    three_prime_positions: List[int],
                                    distance_threshold: int = 24) -> List[Dict[str,List]]:
        """
        Cluster 3' end positions by genomic proximity

        Args:
            three_prime_positions: List of position
            distance_threshold: Maximum distance to cluster positions together

        Returns:
            List of clusters, where each cluster is a list of ThreePrimePosition objects
        """
        if not three_prime_positions:
            return []

        logger.debug(f"Clustering {len(three_prime_positions)} 3' end positions with threshold {distance_threshold}")

        # Sort by position
        sorted_indices = sorted(range(len(three_prime_positions)),key=lambda i: three_prime_positions[i],)
    
        # Step 2: Reorder positions and names using the sorted indices
        sorted_positions = [three_prime_positions[i] for i in sorted_indices]
        sorted_ids = [ids[i] for i in sorted_indices]
    

        clusters = []
        cluster_ids = []
        current_cluster = [sorted_positions[0]]
        current_id = [sorted_ids[0]]
        for i in range(1, len(sorted_positions)):
            curr_pos = sorted_positions[i]
            prev_pos = sorted_positions[i - 1]
            curr_id = sorted_ids[i]

            # Check if within distance threshold
            if curr_pos - prev_pos <= distance_threshold:
                # Add to current cluster
                current_cluster.append(curr_pos)
                current_id.append(curr_id)
            else:
                # Close current cluster and start new one
                if current_cluster:
                    clusters.append(current_cluster)
                    cluster_ids.append(current_id)
                current_cluster = [curr_pos]
                current_id = [curr_id]

        # Add the last cluster
        if current_cluster:
            clusters.append(current_cluster)
            cluster_ids.append(current_id)

        logger.info(f"Created {len(clusters)} clusters from {len(three_prime_positions)} positions")

        # Log cluster sizes
        cluster_sizes = [len(c) for c in clusters]
        logger.debug(f"Cluster sizes: min={min(cluster_sizes) if cluster_sizes else 0}, "
                    f"max={max(cluster_sizes) if cluster_sizes else 0}, "
                    f"mean={sum(cluster_sizes)/len(cluster_sizes) if cluster_sizes else 0:.1f}")

        return [{'ids':i,'3_positions':c} for i,c in zip(cluster_ids,clusters)]

    def prepare_cluster_data(self,
                           ids: List[str],
                           three_prime_positions: List[int],
                           read_seqs: Dict[str,str],
                           ref_seqs: Dict[str,str],
                           strand: str) -> ClusteredData:

        # Prepare the read seq set 
        read_seq_list = [(j, read_seqs[j],three_prime_positions[i]) for i,j in enumerate(ids) if not j.startswith('gtf_')]
        ref_seq_list = [(j.replace('gtf_',''),ref_seqs[j.replace('gtf_','')],three_prime_positions[i]) for i,j in enumerate(ids) if j.startswith('gtf_')]
        contained_read = set()
        
        for read_i in read_seq_list:
      
            read_id, read_seq, read_3_end = read_i
            for ref_i in ref_seq_list:
                ref_id, ref_seq, ref_3_end = ref_i
                if strand == '+':
                    end_dif = read_3_end- ref_3_end
                else:
                    end_dif =  ref_3_end - read_3_end
                if end_dif >= 0:
                    if read_seq[:-end_dif] in ref_seq:
                        contained_read.add(read_id)
                        break
                else:
                    if read_seq in ref_seq[:end_dif]:
                        contained_read.add(read_id)
                        break
                    
        potential_novels = sorted([i for i in read_seq_list if i[0] not in contained_read],key=lambda x: len(x[1]))
        for i, read_i in enumerate(potential_novels):
            read_id, read_seq, read_3_end = read_i
            for other_i in potential_novels[i+1:]:
                if ref_id == read_id:
                    continue
                ref_id, ref_seq, ref_3_end = other_i
                if strand == '+':
                    end_dif = read_3_end- ref_3_end
                else:
                    end_dif =  ref_3_end - read_3_end
                if end_dif >= 0:
                    if read_seq[:-end_dif] in ref_seq:
                        contained_read.add(read_id)
                        break
                else:
                    if read_seq in ref_seq[:end_dif]:
                        contained_read.add(read_id)
                        break
                        
                        
    
    #     # Extract sequences for reads in this cluster
    #     read_sequences = {}
    #     read_qualities = {}
    #     read_cigars = {}
    #     transcript_sequences = {}
    #     transcript_ids = []

    #     chrom = interval.chrom

    #     # Find reads in this cluster
    #     cluster_read_ids = {pos.read_id for pos in cluster_positions if pos.is_read and pos.read_id}

    #     logger.debug(f"Preparing data for cluster {cluster_id} with {len(cluster_read_ids)} reads")

    #     # Extract read sequences
    #     for read in reads:
    #         read_id = read['query_name']
    #         if read_id in cluster_read_ids:
    #             # Get read sequence
    #             seq = read.get('query_sequence')
    #             if seq:
    #                 read_sequences[read_id] = seq

    #             # Get base qualities
    #             qualities = read.get('query_qualities')
    #             if qualities:
    #                 read_qualities[read_id] = qualities

    #             # Get CIGAR
    #             cigar = read.get('cigartuples')
    #             if cigar:
    #                 read_cigars[read_id] = cigar

    #     # Extract transcript sequences
    #     cluster_transcript_ids = {pos.transcript_id for pos in cluster_positions
    #                              if not pos.is_read and pos.transcript_id}

    #     if self.gtf_path and cluster_transcript_ids and transcripts:
    #         with FASTAReader(self.fasta_path) as fasta_reader:
    #             for transcript in transcripts:
    #                 if transcript.transcript_id in cluster_transcript_ids:
    #                     transcript_ids.append(transcript.transcript_id)

    #                     # Get spliced sequence
    #                     chrom_seq = fasta_reader.get_sequence(chrom)
    #                     if chrom_seq:
    #                         transcript_seq = transcript.get_spliced_sequence(chrom_seq)
    #                         transcript_sequences[transcript.transcript_id] = transcript_seq

    #     # Collect all 3' end positions
    #     three_prime_positions = [pos.position for pos in cluster_positions]

    #     cluster_data = ClusteredData(
    #         cluster_id=cluster_id,
    #         chrom=chrom,
    #         strand=interval.strand or '+',
    #         three_prime_positions=three_prime_positions,
    #         read_sequences=read_sequences,
    #         read_qualities=read_qualities,
    #         read_cigars=read_cigars,
    #         transcript_sequences=transcript_sequences,
    #         transcript_ids=transcript_ids,
    #         interval=interval
    #     )

    #     logger.debug(f"Prepared cluster {cluster_id}: {len(read_sequences)} reads, "
    #                 f"{len(transcript_sequences)} transcripts")

    #     return cluster_data


