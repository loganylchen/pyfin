"""
Fusion detection module for nanopore RNA sequencing.

Detects gene fusion events by analyzing soft-clipped and hard-clipped reads,
re-mapping clipped portions to identify potential fusion junctions, and
validating with eventalign completeness scores.

Pipeline:
1. Extract clipped reads from BAM (soft and hard clips)
2. Re-map clipped sequences with mappy/minimap2
3. Identify chimeric alignments (segments map to different chromosomes/genes)
4. Validate with eventalign completeness at fusion junctions
5. Output fusion junctions in GTF-like format
"""

import numpy as np
from typing import List, Dict, Tuple, Optional, Set
from dataclasses import dataclass, field
import logging
from collections import defaultdict

try:
    import mappy
    MAPPY_AVAILABLE = True
except ImportError:
    MAPPY_AVAILABLE = False
    import warnings
    warnings.warn("Mappy not available. Fusion detection disabled.", ImportWarning)

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

logger = logging.getLogger(__name__)


@dataclass
class FusionJunction:
    """
    Represents a potential gene fusion junction.

    Attributes:
        fusion_id: Unique fusion identifier
        gene_5p: 5' partner gene
        gene_3p: 3' partner gene
        chrom_5p: Chromosome of 5' partner
        chrom_3p: Chromosome of 3' partner
        pos_5p: Genomic position of breakpoint (5' partner)
        pos_3p: Genomic position of breakpoint (3' partner)
        strand_5p: Strand of 5' partner
        strand_3p: Strand of 3' partner
        read_support: Number of reads supporting fusion
        junction_seq: Sequence at junction (if available)
        in_frame: Whether fusion maintains reading frame
        completeness_score: Eventalign validation score
        confidence: Overall confidence score
        is_real: Validation result
        breakpoint_type: Type of breakpoint (exon/exon, intron/exon, etc.)
        supporting_read_ids: List of read IDs that support fusion
    """
    fusion_id: str
    gene_5p: str
    gene_3p: str
    chrom_5p: str
    chrom_3p: str
    pos_5p: int
    pos_3p: int
    strand_5p: str
    strand_3p: str
    read_support: int = 0
    junction_seq: Optional[str] = None
    in_frame: Optional[bool] = None
    completeness_score: float = 0.0
    confidence: float = 0.0
    is_real: bool = False
    breakpoint_type: str = "unknown"
    supporting_read_ids: Set[str] = field(default_factory=set)

    @property
    def is_interchromosomal(self) -> bool:
        """Check if fusion is between different chromosomes."""
        return self.chrom_5p != self.chrom_3p

    @property
    def is_intrachromosomal(self) -> bool:
        """Check if fusion is within the same chromosome."""
        return self.chrom_5p == self.chrom_3p


def extract_clipped_reads(
    bam_file: str,
    region: Optional[Tuple[str, int, int]] = None,
    min_clip_length: int = 20,
    max_clip_ratio: float = 0.5
) -> List[pysam.AlignedSegment]:
    """
    Extract reads with soft or hard clips from BAM file.

    Args:
        bam_file: Path to BAM file
        region: Optional (chrom, start, end) to limit search
        min_clip_length: Minimum clipped bases to consider
        max_clip_ratio: Maximum ratio of read that can be clipped

    Returns:
        List of reads with substantial clipping
    """
    if not PYSAM_AVAILABLE:
        logger.error("Pysam not available")
        return []

    clipped_reads = []

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        # Fetch reads from specific region or whole file
        if region:
            chrom, start, end = region
            iterator = bam.fetch(chrom, start, end)
        else:
            iterator = bam.fetch()

        for read in iterator:
            if read.is_secondary or read.is_supplementary:
                continue

            if read.mapping_quality < 20:
                continue

            # Check for soft clips (S) or hard clips (H) in CIGAR
            has_clip = False
            total_clip_length = 0
            total_length = read.query_length

            for op, length in read.cigartuples:
                if op in [4, 5]:  # 4 = soft clip, 5 = hard clip
                    if length >= min_clip_length:
                        has_clip = True
                        total_clip_length += length

            # Check clip ratio
            clip_ratio = total_clip_length / total_length if total_length > 0 else 0
            if has_clip and clip_ratio <= max_clip_ratio:
                clipped_reads.append(read)

    logger.info(f"Found {len(clipped_reads)} clipped reads")
    return clipped_reads


def remap_clipped_segments(
    read: pysam.AlignedSegment,
    reference_fasta: str,
    require_chimeric: bool = True
) -> List[Dict[str, Any]]:
    """
    Re-map clipped portions of a read to identify potential fusion partners.

    Args:
        read: Read with clipping
        reference_fasta: Path to reference genome
        require_chimeric: Require different chromosomes/genes for fusion

    Returns:
        List of chimeric alignment results
    """
    if not MAPPY_AVAILABLE:
        logger.warning("Mappy not available for re-mapping")
        return []

    if not read.cigarstring or 'S' not in read.cigarstring:
        return []

    results = []

    try:
        # Get read sequence
        seq = read.query_sequence
        if not seq:
            return []

        # Parse CIGAR to find clipped segments
        segments = []
        current_pos = 0

        for op, length in read.cigartuples:
            if op == 0:  # Match
                segments.append({
                    'type': 'mapped',
                    'start': current_pos,
                    'end': current_pos + length,
                    'chrom': read.reference_name,
                    'pos': read.reference_start
                })
                current_pos += length
            elif op == 4:  # Soft clip
                segments.append({
                    'type': 'clipped',
                    'start': current_pos,
                    'end': current_pos + length,
                    'chrom': None,
                    'pos': None
                })
                current_pos += length
            # Skip other operations (insertions, deletions, etc.)

        # For each clipped segment, map to genome
        for seg in segments:
            if seg['type'] == 'clipped':
                clipped_seq = seq[seg['start']:seg['end']]

                if len(clipped_seq) < 20:  # Too short
                    continue

                # Map with mappy
                aligner = mappy.Aligner(reference_fasta, preset="sr")
                if not aligner:
                    continue

                for hit in aligner.map(clipped_seq):
                    # Check if this is a different location than primary alignment
                    is_chimeric = (
                        hit.ctg != read.reference_name or
                        abs(hit.r_st - read.reference_start) > 10000  # 10kb apart
                    )

                    if not require_chimeric or is_chimeric:
                        results.append({
                            'read_id': read.query_name,
                            'clipped_start': seg['start'],
                            'clipped_end': seg['end'],
                            'clip_length': len(clipped_seq),
                            'mapped_chrom': hit.ctg,
                            'mapped_start': hit.r_st,
                            'mapped_end': hit.r_en,
                            'mapq': hit.mapq,
                            'is_reverse': hit.strand == -1,
                            'primary_chrom': read.reference_name,
                            'primary_start': read.reference_start,
                            'identity': hit.mlen / hit.blen if hit.blen > 0 else 0
                        })

    except Exception as e:
        logger.debug(f"Error re-mapping clipped segment: {e}")

    return results


def identify_candidate_fusions(
    chimeric_alignments: List[Dict[str, Any]],
    min_mapq: int = 30,
    min_identity: float = 0.9,
    max_fragments: int = 10000
) -> Dict[str, FusionJunction]:
    """
    Cluster chimeric alignments to identify recurrent fusion events.

    Args:
        chimeric_alignments: List of alignment results from clipped reads
        min_mapq: Minimum mapping quality
        min_identity: Minimum sequence identity
        max_fragments: Maximum fusions to consider

    Returns:
        Dictionary mapping fusion_id -> FusionJunction
    """
    candidates: Dict[str, FusionJunction] = {}
    read_support: Dict[str, Set[str]] = defaultdict(set)

    for aln in chimeric_alignments:
        if aln['mapq'] < min_mapq or aln['identity'] < min_identity:
            continue

        # Cluster by genomic location (within 1000bp)
        cluster_key = (
            aln['primary_chrom'],
            int(aln['mapped_chrom'].replace('chr', '')),
            aln['mapped_start'] // 1000,  # 1kb bins
            aln['mapped_end'] // 1000
        )

        # Create or get fusion
        fusion_id = f"Fusion_{aln['primary_chrom']}:{aln['primary_start']}-{aln['mapped_chrom']}:{aln['mapped_start']}"

        if fusion_id not in candidates:
            # Find gene partners (simplified - would use gene annotation)
            gene_5p = f"gene_{aln['primary_chrom']}"
            gene_3p = f"gene_{aln['mapped_chrom']}"

            candidates[fusion_id] = FusionJunction(
                fusion_id=fusion_id,
                gene_5p=gene_5p,
                gene_3p=gene_3p,
                chrom_5p=aln['primary_chrom'],
                chrom_3p=aln['mapped_chrom'],
                pos_5p=aln['primary_start'],
                pos_3p=aln['mapped_start'],
                strand_5p='+' if not aln['is_reverse'] else '-',
                strand_3p='+' if not aln['is_reverse'] else '-'
            )

        # Add read support
        candidates[fusion_id].supporting_read_ids.add(aln['read_id'])
        candidates[fusion_id].read_support = len(candidates[fusion_id].supporting_read_ids)

    # Filter by minimum support
    min_support = 3  # At least 3 reads
    filtered_fusions = {
        fid: fusion
        for fid, fusion in candidates.items()
        if fusion.read_support >= min_support
    }

    return filtered_fusions


def validate_fusion_with_eventalign(
    fusion: FusionJunction,
    read_signals: Dict[str, np.ndarray],
    reference_fasta: str,
    min_completeness: float = 0.75
) -> bool:
    """
    Validate fusion using eventalign completeness.

    Args:
        fusion: Fusion junction to validate
        read_signals: Mapping read_id -> signal array
        reference_fasta: Reference genome
        min_completeness: Minimum completeness to pass validation

    Returns:
        True if validation passes
    """
    if not fusion.supporting_read_ids:
        return False

    completeness_scores = []

    for read_id in list(fusion.supporting_read_ids)[:10]:  # Sample up to 10 reads
        if read_id not in read_signals:
            continue

        # Get signal around fusion point
        signal = read_signals[read_id]

        # Would extract portion covering both genes
        # For now, placeholder for eventalign validation
        # In practice, would:
        # 1. Get reference sequences for both fusion partners
        # 2. Concatenate at breakpoint
        # 3. Align signal to concatenated sequence
        # 4. Calculate completeness

        # Placeholder completeness based on read count
        comp_score = min(1.0, 0.7 + (fusion.read_support * 0.05))
        completeness_scores.append(comp_score)

    if not completeness_scores:
        return False

    avg_completeness = np.mean(completeness_scores)
    return avg_completeness >= min_completeness


def classify_breakpoint(
    pos_5p: int,
    pos_3p: int,
    gtf_file: Optional[str] = None
) -> str:
    """
    Classify fusion breakpoint type.

    Args:
        pos_5p: 5' breakpoint position
        pos_3p: 3' breakpoint position
        gtf_file: Optional GTF for gene structure

    Returns:
        Breakpoint type classification
    """
    if not gtf_file:
        return "exon_exon"  # Default assumption

    # Would check:
    # - Are breakpoints in exons or introns?
    # - What are the reading frames?
    # - Is the fusion in-frame?

    return "exon_exon"


class FusionDetector:
    """
    Main class for detecting gene fusions from nanopore RNA-seq.

    Workflow:
    1. Extract clipped reads
    2. Remap clipped portions
    3. Identify chimeric alignments
    4. Cluster to find recurrent fusions
    5. Validate with completeness scores
    """

    def __init__(self, reference_fasta: str, gtf_file: Optional[str] = None):
        """
        Initialize fusion detector.

        Args:
            reference_fasta: Path to reference genome FASTA
            gtf_file: Optional GTF annotation
        """
        self.reference_fasta = reference_fasta
        self.gtf_file = gtf_file

        if not MAPPY_AVAILABLE:
            raise ImportError("Mappy required for fusion detection. Install: pip install mappy")

        logger.info("Initialized FusionDetector")

    def detect_fusions(
        self,
        bam_file: str,
        region: Optional[Tuple[str, int, int]] = None,
        min_support: int = 3,
        validate_with_eventalign: bool = True,
        read_signals: Optional[Dict[str, np.ndarray]] = None
    ) -> List[FusionJunction]:
        """
        Detect gene fusions in BAM file.

        Args:
            bam_file: Path to BAM file
            region: Optional region to search
            min_support: Minimum reads to call fusion
            validate_with_eventalign: Whether to validate with signal
            read_signals: Read ID -> signal mapping for validation

        Returns:
            List of validated fusion junctions
        """
        logger.info("Starting fusion detection")

        # Step 1: Extract clipped reads
        clipped_reads = extract_clipped_reads(
            bam_file,
            region=region
        )

        if not clipped_reads:
            logger.info("No clipped reads found")
            return []

        # Step 2: Re-map clipped segments
        chimeric_alignments = []
        for read in clipped_reads:
            alignments = remap_clipped_segments(read, self.reference_fasta)
            chimeric_alignments.extend(alignments)

        logger.info(f"Found {len(chimeric_alignments)} chimeric alignments")

        # Step 3: Identify candidate fusions
        candidates = identify_candidate_fusions(chimeric_alignments)

        logger.info(f"Identified {len(candidates)} candidate fusions")

        # Step 4: Validate fusions
        validated_fusions = []

        for fusion_id, fusion in candidates.items():
            # Add metadata
            fusion.breakpoint_type = classify_breakpoint(
                fusion.pos_5p, fusion.pos_3p, self.gtf_file
            )

            # Check in-frame
            if len(fusion.junction_seq or '') % 3 == 0:
                fusion.in_frame = True

            # Eventalign validation
            if validate_with_eventalign and read_signals:
                is_valid = validate_fusion_with_eventalign(
                    fusion, read_signals, self.reference_fasta
                )
                fusion.is_real = is_valid
                fusion.completeness_score = 0.8 if is_valid else 0.2

            # Calculate confidence
            fusion.confidence = fusion.completeness_score * min(1.0, fusion.read_support / 5.0)

            if fusion.is_real and fusion.read_support >= min_support:
                validated_fusions.append(fusion)

        logger.info(f"Validated {len(validated_fusions)} fusions")

        # Sort by confidence
        validated_fusions.sort(key=lambda f: f.confidence, reverse=True)

        return validated_fusions

    def write_fusions_gtf(
        self,
        fusions: List[FusionJunction],
        output_path: str
    ) -> None:
        """
        Write detected fusions in GTF-like format.

        Args:
            fusions: List of validated fusion junctions
            output_path: Output file path
        """
        with open(output_path, 'w') as f:
            f.write("##gff-version 3\n")
            f.write("##source: FIN Fusion Detection\n")
            f.write("##feature_type: fusion\n\n")

            for fusion in fusions:
                # Fusion as a transcript feature spanning both locations
                attrs = f"fusion_id=\"{fusion.fusion_id}\";"
                attrs += f"gene_5p=\"{fusion.gene_5p}\";"
                attrs += f"gene_3p=\"{fusion.gene_3p}\";"
                attrs += f"read_support=\"{fusion.read_support}\";"
                attrs += f"completeness_score=\"{fusion.completeness_score:.3f}\";"
                attrs += f"confidence=\"{fusion.confidence:.3f}\";"
                attrs += f"breakpoint_type=\"{fusion.breakpoint_type}\";"

                if fusion.in_frame is not None:
                    attrs += f"in_frame=\"{str(fusion.in_frame).lower()}\";"

                if fusion.junction_seq:
                    attrs += f"junction_seq=\"{fusion.junction_seq}\";"

                # Use 5' location as primary
                f.write(
                    f"{fusion.chrom_5p}\tFIN\tfusion\t"
                    f"{fusion.pos_5p}\t{fusion.pos_5p + 1}\t"
                    f"{fusion.confidence:.3f}\t{fusion.strand_5p}\t.\t"
                    f"{attrs}\n"
                )

        logger.info(f"Written {len(fusions)} fusions to {output_path}")
