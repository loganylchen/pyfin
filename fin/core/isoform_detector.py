"""
Main pipeline for isoform detection from nanopore direct RNA-seq.

Integrates all components:
- GTF parsing
- Sequence extraction
- Event alignment
- Completeness scoring
- DTW clustering
- Integration and validation
- Fusion detection (optional)
"""

import numpy as np
from typing import List, Tuple, Optional, Dict, Any
import logging
from dataclasses import dataclass

from fin.core.isoform import IsoformSequence, ValidatedIsoform, IsoformEvidence
from fin.core.completeness import calculate_completeness
from fin.core.integration_matrix import validate_isoforms, cluster_reads_by_signal, build_completeness_matrix
from fin.core.fusion_detector import FusionDetector, FusionJunction
from fin.io.io_gtf import GtfParser, TranscriptFeature
from fin.io.io_manager import IOManager, ReadData
from fin.io.sequence_extractor import ReferenceExtractor
from fin.core.f5c_wrapper import F5CWrapper
from fin.utils.config import PipelineConfig

logger = logging.getLogger(__name__)


class IsoformDetector:
    """
    Main pipeline for isoform and fusion detection from nanopore direct RNA-seq.

    Processes a single gene region through the complete pipeline:
    1. Extract isoforms (annotated + read-derived)
    2. Align signals to isoforms via eventalign
    3. Calculate completeness scores
    4. Cluster reads by signal similarity
    5. Integrate metrics to validate isoforms
    6. Optional: Detect gene fusions from clipped reads
    """

    def __init__(
        self,
        config: PipelineConfig,
        reference_fasta: str,
        signal_dir: Optional[str] = None,
        use_cuda: bool = True,
        enable_fusion_detection: bool = False
    ):
        """
        Initialize detector.

        Args:
            config: PipelineConfiguration object
            reference_fasta: Path to reference genome FASTA
            signal_dir: Directory containing signal files (FAST5/POD5/SLOW5)
            use_cuda: Enable GPU acceleration
            enable_fusion_detection: Whether to detect gene fusions
        """
        self.config = config
        self.use_cuda = use_cuda
        self.enable_fusion_detection = enable_fusion_detection

        # Initialize components
        self.seq_extractor = ReferenceExtractor(reference_fasta)
        self.io_manager = IOManager() if signal_dir else None
        self.f5c = F5CWrapper()
        self.fusion_detector = FusionDetector(reference_fasta) if enable_fusion_detection else None

        logger.info("IsoformDetector initialized")

    def process_gene_region(
        self,
        chrom: str,
        start: int,
        end: int,
        bam_file: str,
        gtf_file: Optional[str] = None,
        detect_fusions: bool = False
    ) -> Tuple[List[ValidatedIsoform], Dict[str, str], Optional[List[FusionJunction]]]:
        """
        Process a single gene region.

        Args:
            chrom: Chromosome name
            start: Start coordinate (0-based)
            end: End coordinate
            bam_file: Path to BAM file
            gtf_file: Optional GTF file with annotations
            detect_fusions: Whether to detect gene fusions

        Returns:
            Tuple of:
            - List of ValidatedIsoform objects
            - Dictionary mapping read_id -> assigned_isoform_id
            - List of FusionJunction (if fusion detection enabled)
        """
        # Step 1: Extract isoform sequences
        isoform_sequences = self._extract_isoform_sequences(
            chrom, start, end, bam_file, gtf_file
        )

        logger.info(f"Extracted {len(isoform_sequences)} isoform sequences")

        # Step 2: Extract reads and their signals
        reads_data = self._extract_reads_in_region(chrom, start, end, bam_file)

        if not reads_data:
            logger.warning(f"No reads found in region {chrom}:{start}-{end}")
            return [], {}, None

        logger.info(f"Extracted {len(reads_data)} reads")

        # Step 3: Eventalign each read to each isoform
        logger.info(f"Aligning {len(reads_data)} reads to {len(isoform_sequences)} isoforms...")
        completeness_scores = self._align_all_reads_to_isoforms(
            reads_data, isoform_sequences
        )

        # Step 4: Build completeness matrix
        completeness_matrix = build_completeness_matrix(
            reads_data, isoform_sequences, completeness_scores
        )

        # Step 5: Cluster reads by signal similarity (DTW)
        read_signals = [read.signal for read in reads_data if read.signal is not None]
        dtw_labels = cluster_reads_by_signal(
            read_signals,
            self.config.algorithm.dtw_similarity_threshold,
            use_cuda=self.use_cuda
        ) if read_signals else np.arange(len(reads_data))

        # Step 6: Integrate metrics and validate isoforms
        validated_isoforms = self._validate_isoforms(
            isoform_sequences, completeness_matrix, dtw_labels, reads_data
        )

        logger.info(f"Validated {len(validated_isoforms)} isoforms")

        # Step 7: Assign reads to isoforms
        read_assignments = self._assign_reads_to_isoforms(
            completeness_matrix, validated_isoforms
        )

        logger.info(f"Assigned {len(read_assignments)} reads to isoforms")

        # Optional: Fusion detection
        detected_fusions = None
        if detect_fusions and self.fusion_detector:
            logger.info("Detecting gene fusions...")
            detected_fusions = self.fusion_detector.detect_fusions(
                bam_file=bam_file,
                region=(chrom, start, end),
                min_support=self.config.algorithm.min_read_support,
                validate_with_eventalign=True,
                read_signals={read.read_id: read.signal for read in reads_data if read.signal is not None}
            )
            logger.info(f"Detected {len(detected_fusions)} fusions")

        return validated_isoforms, read_assignments, detected_fusions

    def _extract_isoform_sequences(
        self,
        chrom: str,
        start: int,
        end: int,
        bam_file: str,
        gtf_file: Optional[str]
    ) -> List[IsoformSequence]:
        """
        Extract isoform sequences from annotations and reads.

        Returns both annotated (from GTF) and read-derived isoforms.
        """
        isoforms = []

        # Extract annotated isoforms from GTF
        if gtf_file:
            try:
                gtf_parser = GtfParser(gtf_file)
                annotated_isoforms = gtf_parser.get_transcripts_in_region(
                    chrom, start, end
                )

                # Convert to IsoformSequence objects
                for transcript in annotated_isoforms:
                    seq = self.seq_extractor.extract_transcript_sequence(
                        transcript.chrom,
                        transcript.exons,
                        transcript.strand
                    )

                    isoform = IsoformSequence(
                        isoform_id=transcript.transcript_id,
                        gene_id=transcript.gene_id,
                        chrom=transcript.chrom,
                        strand=transcript.strand,
                        exons=transcript.exons,
                        sequence=seq,
                        source="annotation",
                        metadata=transcript.attrs
                    )
                    isoforms.append(isoform)

                logger.info(f"Extracted {len(annotated_isoforms)} annotated isoforms")

            except Exception as e:
                logger.warning(f"Failed to parse GTF: {e}")

        # Extract read-derived isoforms from BAM
        # TODO: Implement read-derived isoform extraction
        # For now, only using annotated isoforms

        return isoforms

    def _extract_reads_in_region(
        self,
        chrom: str,
        start: int,
        end: int,
        bam_file: str
    ) -> List[ReadData]:
        """Extract all reads overlapping the region."""
        io_manager = IOManager()
        io_manager.add_alignment_file(bam_file)

        # Get reads in region
        reads = io_manager.get_reads_in_region(chrom, start, end)

        logger.info(f"Found {len(reads)} reads in region")

        # If signal files available, attach signals
        if self.io_manager:
            reads_with_signal = 0
            for read in reads:
                try:
                    signal = self.io_manager.get_signal_by_read_id(read.read_id)
                    read.signal = signal
                    reads_with_signal += 1
                except KeyError:
                    logger.debug(f"No signal found for read {read.read_id}")

            logger.info(f"Attached signals to {reads_with_signal} reads")

        return reads

    def _align_all_reads_to_isoforms(
        self,
        reads_data: List[ReadData],
        isoforms: List[IsoformSequence]
    ) -> Dict[Tuple[str, str], Any]:
        """
        Align each read to each isoform via eventalign.

        Returns dictionary mapping (read_id, isoform_id) -> completeness_result
        """
        from fin.core.completeness import calculate_completeness
        from fin.core.eventalign import AlignedEvent

        completeness_scores = {}
        total_alignments = 0
        successful_alignments = 0

        for read_idx, read in enumerate(reads_data):
            if read.signal is None:
                continue

            # Detect events in signal
            try:
                events = self.f5c.detect_events(read.signal)
            except Exception as e:
                logger.debug(f"Event detection failed for {read.read_id}: {e}")
                continue

            # Align to each isoform
            for isoform in isoforms:
                try:
                    # Align events to isoform sequence
                    aligned = self.f5c.align_to_sequence(
                        events, isoform.sequence, is_rna=True
                    )

                    # Calculate completeness
                    try:
                        comp_result = calculate_completeness(
                            aligned, isoform.sequence,
                            min_event_duration=self.config.eventalign.min_event_duration
                        )
                        completeness_scores[(read.read_id, isoform.isoform_id)] = comp_result
                        successful_alignments += 1
                    except Exception as e:
                        logger.debug(f"Completeness calculation failed: {e}")
                        # Use zeros for failed calculations
                        from fin.core.completeness import CompletenessResult
                        completeness_scores[(read.read_id, isoform.isoform_id)] = CompletenessResult()

                    total_alignments += 1

                except Exception as e:
                    logger.debug(f"Alignment failed for {read.read_id}-{isoform.isoform_id}: {e}")
                    # Use zeros for failed alignments
                    from fin.core.completeness import CompletenessResult
                    completeness_scores[(read.read_id, isoform.isoform_id)] = CompletenessResult()

        logger.info(f"Eventalign completed: {successful_alignments}/{total_alignments} successful")

        return completeness_scores

    def _validate_isoforms(
        self,
        isoforms: List[IsoformSequence],
        completeness_matrix: np.ndarray,
        dtw_labels: np.ndarray,
        reads: List[ReadData]
    ) -> List[ValidatedIsoform]:
        """Validate which isoforms are real."""
        from fin.core.integration_matrix import validate_isoforms

        validated_indices = validate_isoforms(
            completeness_matrix,
            dtw_labels,
            self.config.algorithm.min_read_support,
            self.config.algorithm.min_completeness_threshold,
            self.config.algorithm.cluster_consistency_threshold
        )

        validated_isoforms = []
        for isoform_idx, read_indices in validated_indices.items():
            isoform = isoforms[isoform_idx]

            # Calculate support metrics
            read_ids = [reads[ri].read_id for ri in read_indices if ri < len(reads)]
            completeness_scores = [
                completeness_matrix[ri, isoform_idx] for ri in read_indices if ri < len(reads)
            ]

            validated = ValidatedIsoform(
                isoform_id=isoform.isoform_id,
                gene_id=isoform.gene_id,
                chrom=isoform.chrom,
                strand=isoform.strand,
                exons=isoform.exons,
                sequence=isoform.sequence,
                read_support=len(read_indices),
                avg_completeness=np.mean(completeness_scores) if completeness_scores else 0.0,
                confidence_score=np.mean(completeness_scores) * len(read_indices) if completeness_scores else 0.0,
                supporting_reads=read_ids,
                is_novel=(isoform.source == "read")
            )
            validated_isoforms.append(validated)

            logger.debug(
                f"Validated {isoform.isoform_id}: {len(read_indices)} reads, "
                f"completeness={np.mean(completeness_scores):.3f}"
            )

        return validated_isoforms

    def _assign_reads_to_isoforms(
        self,
        completeness_matrix: np.ndarray,
        validated_isoforms: List[ValidatedIsoform]
    ) -> Dict[str, str]:
        """
        Assign each read to the most likely isoform.

        Returns:
            Dict mapping read_id -> isoform_id
        """
        if not validated_isoforms:
            return {}

        read_assignments = {}

        # Get mapping from isoform_id -> matrix column index
        isoform_columns = {
            valid.isoform_id: i for i, valid in enumerate(validated_isoforms)
        }

        # Find best isoform for each read
        for read_idx in range(completeness_matrix.shape[0]):
            best_isoform = None
            best_score = 0.0

            for valid_iso in validated_isoforms:
                if valid_iso.isoform_id not in isoform_columns:
                    continue

                col_idx = isoform_columns[valid_iso.isoform_id]
                score = completeness_matrix[read_idx, col_idx]

                if score > best_score and score >= self.config.algorithm.min_completeness_threshold:
                    best_score = score
                    best_isoform = valid_iso.isoform_id

            if best_isoform:
                read_assignments[f"read_{read_idx}"] = best_isoform

        return read_assignments

    def quantify_isoforms(
        self,
        validated_isoforms: List[ValidatedIsoform],
        completeness_scores: Dict[Tuple[str, str], Any]
    ) -> None:
        """
        Quantify isoform expression levels.

        Calculates:
        - estimated_count: Estimated number of reads
        - tpm: Transcripts Per Million
        - fpkm: Fragments Per Kilobase Million
        - isoform_fraction: Fraction of gene's expression

        Args:
            validated_isoforms: List of validated isoforms (modified in-place)
            completeness_scores: Completeness results for scoring
        """
        # TODO: Implement EM quantification
        # For now, basic weighting

        for isoform in validated_isoforms:
            # Basic count estimate (raw reads weighted by completeness)
            total_completeness = sum(
                scores.completeness_score for scores in completeness_scores.values()
                if scores  # Check if exists
            )

            if total_completeness > 0:
                self_count = sum(
                    scores.completeness_score
                    for (read_id, iso_id), scores in completeness_scores.items()
                    if iso_id == isoform.isoform_id and scores
                )
                isoform.estimated_count = self_count
            else:
                isoform.estimated_count = float(isoform.read_support)

            # Calculate TPM (simplified)
            if isoform.length > 0:
                isoform.fpkm = isoform.estimated_count / (isoform.length / 1000.0)
                isoform.tpm = isoform.fpkm

            # Isoform fraction of gene
            # Would need all isoforms of same gene for proper calculation
            isoform.isoform_fraction = min(1.0, isoform.read_support / max(1, isoform.read_support))

    def run_complete_pipeline(
        self,
        genes: List[str],
        bam_file: str,
        gtf_file: str,
        detect_fusions: bool = False
    ) -> Dict[str, Dict[str, Any]]:
        """
        Run complete pipeline on multiple genes.

        Args:
            genes: List of gene names or regions
            bam_file: Path to BAM file
            gtf_file: Path to GTF annotations
            detect_fusions: Whether to detect fusions

        Returns:
            Dictionary mapping gene -> results dict
        """
        all_results = {}

        for gene in genes:
            logger.info(f"Processing gene: {gene}")

            # Parse region
            if ':' in gene:
                chrom, coords = gene.split(':', 1)
                start, end = map(int, coords.split('-', 1))
            else:
                # Look up in GTF
                from fin.io.io_gtf import GtfParser
                parser = GtfParser(gtf_file)
                gene_transcripts = parser.get_all_genes().get(gene, [])

                if not gene_transcripts:
                    logger.warning(f"Gene {gene} not found in GTF")
                    continue

                chrom = gene_transcripts[0].chrom
                start = min(min(e[0] for e in t.exons) for t in gene_transcripts)
                end = max(max(e[1] for e in t.exons) for t in gene_transcripts)

            # Run detection
            isoforms, assignments, fusions = self.process_gene_region(
                chrom=chrom,
                start=start,
                end=end,
                bam_file=bam_file,
                gtf_file=gtf_file,
                detect_fusions=detect_fusions
            )

            all_results[gene] = {
                'validated_isoforms': isoforms,
                'read_assignments': assignments,
                'detected_fusions': fusions,
                'chrom': chrom,
                'start': start,
                'end': end
            }

        return all_results
