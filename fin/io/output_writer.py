"""
Output writers for isoform detection results.

Provides functionality to write:
1. GTF format: Validated isoforms with coordinates and metadata
2. BED12 format: Read-to-isoform assignments with mapping information
"""

from typing import List, Dict, Any, Tuple, Optional
import logging
import pysam

from fin.core.isoform import ValidatedIsoform

logger = logging.getLogger(__name__)


def write_gtf(
    validated_isoforms: List[ValidatedIsoform],
    output_path: str,
    source: str = "fin",
    write_version: bool = True
) -> None:
    """
    Write validated isoforms in GTF format.

    GTF format specification:
    chrom\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes

    Features written:
    - gene: Represents the gene locus
    - transcript: Represents the transcript/isoform
    - exon: Individual exons within the transcript

    Args:
        validated_isoforms: List of ValidatedIsoform objects
        output_path: Path to output GTF file
        source: Source field (default: "fin")
        write_version: Whether to write GTF version header

    Examples:
        >>> isoforms = detector.process_gene_region(chrom, start, end, bam, gtf)
        >>> write_gtf(isoforms, "output.isoforms.gtf")
    """
    if not validated_isoforms:
        logger.warning("No validated isoforms to write")
        return

    # Group isoforms by gene
    gene_isoforms: Dict[str, List[ValidatedIsoform]] = {}
    for isoform in validated_isoforms:
        if isoform.gene_id not in gene_isoforms:
            gene_isoforms[isoform.gene_id] = []
        gene_isoforms[isoform.gene_id].append(isoform)

    with open(output_path, 'w') as f:
        # Write GTF version header
        if write_version:
            f.write("##gff-version 3\n")
            f.write("##source: fin (nanopore isoform detection)\n")
            f.write(f"##isoforms: {len(validated_isoforms)}\n\n")

        # Write genes and their isoforms
        for gene_id, isoforms in gene_isoforms.items():
            # Calculate gene bounds (spanning all isoforms)
            all_starts = [min(e[0] for e in iso.exons) for iso in isoforms]
            all_ends = [max(e[1] for e in iso.exons) for iso in isoforms]

            gene_start = min(all_starts)
            gene_end = max(all_ends)
            gene_chrom = isoforms[0].chrom
            gene_strand = isoforms[0].strand

            # Calculate average score for gene
            avg_gene_score = sum(iso.confidence_score for iso in isoforms) / len(isoforms)

            # 1. Write gene feature
            gene_attrs = f"gene_id=\"{gene_id}\";"
            gene_attrs += f"gene_name=\"{gene_id}\";"
            gene_attrs += f"num_isoforms=\"{len(isoforms)}\";"
            gene_attrs += f"avg_confidence=\"{avg_gene_score:.3f}\";"

            f.write(
                f"{gene_chrom}\t{source}\tgene\t{gene_start}\t{gene_end}\t",
                f"{avg_gene_score:.3f}\t{gene_strand}\t.\t{gene_attrs}\n"
            )

            # Write each isoform
            for isoform in isoforms:
                # 2. Write transcript feature
                trans_attrs = f"gene_id=\"{isoform.gene_id}\";"
                trans_attrs += f"transcript_id=\"{isoform.isoform_id}\";"
                trans_attrs += f"read_support=\"{isoform.read_support}\";"
                trans_attrs += f"avg_completeness=\"{isoform.avg_completeness:.3f}\";"
                trans_attrs += f"is_novel=\"{str(isoform.is_novel).lower()}\";"

                # Add quantification metrics
                if hasattr(isoform, 'estimated_count') and isoform.estimated_count is not None:
                    trans_attrs += f"estimated_count=\"{isoform.estimated_count:.2f}\";"
                if hasattr(isoform, 'tpm') and isoform.tpm is not None:
                    trans_attrs += f"TPM=\"{isoform.tpm:.2f}\";"
                if hasattr(isoform, 'fpkm') and isoform.fpkm is not None:
                    trans_attrs += f"FPKM=\"{isoform.fpkm:.2f}\";"
                if hasattr(isoform, 'isoform_fraction') and isoform.isoform_fraction is not None:
                    trans_attrs += f"isoform_fraction=\"{isoform.isoform_fraction:.3f}\";"

                # Add metadata fields
                if isoform.metadata:
                    for key, value in isoform.metadata.items():
                        if key not in ['gene_id', 'transcript_id']:
                            trans_attrs += f"{key}=\"{value}\";"


                trans_start = min(e[0] for e in isoform.exons)
                trans_end = max(e[1] for e in isoform.exons)

                f.write(
                    f"{isoform.chrom}\t{source}\ttranscript\t",
                    f"{trans_start}\t{trans_end}\t",
                    f"{isoform.confidence_score:.3f}\t",
                    f"{isoform.strand}\t.\t{trans_attrs}\n"
                )

                # 3. Write exon features
                for i, (exon_start, exon_end) in enumerate(isoform.exons, 1):
                    exon_attrs = f"gene_id=\"{isoform.gene_id}\";"
                    exon_attrs += f"transcript_id=\"{isoform.isoform_id}\";"
                    exon_attrs += f"exon_number=\"{i}\";"
                    exon_attrs += f"exon_type=\"{isoform.metadata.get('exon_type', 'exon')}\";"

                    f.write(
                        f"{isoform.chrom}\t{source}\texon\t",
                        f"{exon_start}\t{exon_end}\t",
                        f"{isoform.confidence_score:.3f}\t",
                        f"{isoform.strand}\t.\t{exon_attrs}\n"
                    )

    logger.info(f"Written {len(validated_isoforms)} validated isoforms to {output_path}")


def write_bed12(
    read_assignments: Dict[str, str],
    bam_file: str,
    validated_isoforms: List[ValidatedIsoform],
    output_path: str,
    min_mapq: int = 10
) -> None:
    """
    Write read-to-isoform assignments in BED12 format.

    BED12 fields:
    chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd,
    itemRgb, blockCount, blockSizes, blockStarts

    The 'thick' fields are used to indicate the isoform coordinates.

    Args:
        read_assignments: Dictionary mapping read_id -> isoform_id
        bam_file: Path to BAM file (for read alignment details)
        validated_isoforms: List of ValidatedIsoform objects
        output_path: Path to output BED12 file
        min_mapq: Minimum mapping quality to include

    Examples:
        >>> validated_isoforms, assignments = detector.process_gene_region(...)
        >>> write_bed12(assignments, "alignments.bam", validated_isoforms, "output.reads.bed12")
    """
    if not read_assignments:
        logger.warning("No read assignments to write")
        return

    # Create lookup for isoform information
    isoform_lookup: Dict[str, ValidatedIsoform] = {
        isoform.isoform_id: isoform for isoform in validated_isoforms
    }

    # Track statistics
    total_reads = 0
    assigned_reads = 0
    unassigned_reads = 0

    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        with open(output_path, 'w') as f:
            # Write BED header
            f.write("track name=FNIsoforms description=\"FIN Isoform Assignments\" useScore=1\n")

            # Process each read in BAM file
            for read in bam.fetch():
                if read.is_secondary or read.is_supplementary:
                    continue

                if read.mapping_quality < min_mapq:
                    continue

                total_reads += 1
                read_id = read.query_name

                # Check if read is assigned to an isoform
                if read_id not in read_assignments:
                    unassigned_reads += 1
                    continue

                isoform_id = read_assignments[read_id]

                if isoform_id not in isoform_lookup:
                    logger.warning(f"Isoform {isoform_id} not found for read {read_id}")
                    unassigned_reads += 1
                    continue

                isoform = isoform_lookup[isoform_id]
                assigned_reads += 1

                # Get alignment blocks (exons)
                blocks = read.get_blocks()
                if not blocks:
                    continue

                # BED coordinate system is 0-based, half-open
                chrom_start = read.reference_start
                chrom_end = read.reference_end

                # Read name with isoform info
                name = f"{read_id}|isoform={isoform_id}|gene={isoform.gene_id}"

                # Score: high for good alignments, low for poor ones
                score = min(1000, max(0, int(read.mapping_quality * 10)))

                # Strand
                strand = '-' if read.is_reverse else '+'

                # Use thick fields to show isoform coordinates
                thick_start = min(e[0] for e in isoform.exons)
                thick_end = max(e[1] for e in isoform.exons)

                # Item RGB color: green for novel, blue for annotated
                if isoform.is_novel:
                    item_rgb = "0,255,0"  # Green for novel
                else:
                    item_rgb = "0,0,255"  # Blue for annotated

                # Block information
                block_count = len(blocks)
                block_sizes = [str(end - start) for start, end in blocks]
                block_starts = [str(start - chrom_start) for start, end in blocks]

                # Write BED12 line
                f.write(
                    f"{read.reference_name}\t"
                    f"{chrom_start}\t"
                    f"{chrom_end}\t"
                    f"{name}\t"
                    f"{score}\t"
                    f"{strand}\t"
                    f"{thick_start}\t"
                    f"{thick_end}\t"
                    f"{item_rgb}\t"
                    f"{block_count}\t"
                    f"{','.join(block_sizes)}\t"
                    f"{','.join(block_starts)}\n"
                )

    logger.info(
        f"Written {assigned_reads} assigned reads to {output_path} "
        f"({unassigned_reads} unassigned out of {total_reads} total)"
    )


def write_read_assignments_tsv(
    read_assignments: Dict[str, str],
    completeness_results: Dict[Tuple[str, str], 'CompletenessResult'],
    output_path: str
) -> None:
    """
    Write read assignments as TSV file with details.

    More detailed than BED12, includes completeness scores.

    Columns:
    - read_id
    - assigned_isoform
    - completeness_score
    - num_aligned_events
    - alignment_quality

    Args:
        read_assignments: Mapping read_id -> isoform_id
        completeness_results: Mapping (read_id, isoform_id) -> CompletenessResult
        output_path: Path to output TSV file
    """
    with open(output_path, 'w') as f:
        # Write header
        f.write("read_id\tassigned_isoform\tcompleteness_score\tnum_aligned_events\talignment_quality\n")

        # Write assignments
        for read_id, isoform_id in read_assignments.items():
            # Get completeness result
            comp_result = completeness_results.get((read_id, isoform_id))

            if comp_result:
                comp_score = f"{comp_result.completeness_score:.3f}"
                num_events = str(comp_result.num_aligned)
                quality = f"{comp_result.mean_event_quality:.3f}"
            else:
                comp_score = "NA"
                num_events = "NA"
                quality = "NA"

            f.write(f"{read_id}\t{isoform_id}\t{comp_score}\t{num_events}\t{quality}\n")

    logger.info(f"Written read assignments TSV to {output_path}")


def write_metrics_summary(
    validated_isoforms: List[ValidatedIsoform],
    completeness_matrix: 'np.ndarray',
    output_path: str
) -> None:
    """
    Write summary metrics for all isoforms.

    Args:
        validated_isoforms: List of ValidatedIsoform objects
        completeness_matrix: N x M matrix of completeness scores
        output_path: Path to output summary file
    """
    import json

    summary = {
        "total_isoforms": len(validated_isoforms),
        "novel_isoforms": sum(1 for iso in validated_isoforms if iso.is_novel),
        "annotated_isoforms": sum(1 for iso in validated_isoforms if not iso.is_novel),
        "isoforms_by_gene": {},
        "completeness_statistics": {}
    }

    # Group by gene
    for isoform in validated_isoforms:
        gene = isoform.gene_id
        if gene not in summary["isoforms_by_gene"]:
            summary["isoforms_by_gene"][gene] = []

        summary["isoforms_by_gene"][gene].append({
            "isoform_id": isoform.isoform_id,
            "read_support": isoform.read_support,
            "avg_completeness": f"{isoform.avg_completeness:.3f}",
            "confidence": f"{isoform.confidence_score:.3f}",
            "is_novel": isoform.is_novel,
            "length": isoform.length,
            "num_exons": isoform.num_exons
        })

    # Calculate matrix statistics if available
    if completeness_matrix is not None and completeness_matrix.size > 0:
        flat_scores = completeness_matrix.flatten()
        flat_scores = flat_scores[flat_scores > 0]  # Non-zero only

        if len(flat_scores) > 0:
            summary["completeness_statistics"] = {
                "mean": float(flat_scores.mean()),
                "median": float(np.median(flat_scores)),
                "std": float(flat_scores.std()),
                "min": float(flat_scores.min()),
                "max": float(flat_scores.max()),
                "num_comparisons": len(flat_scores)
            }

    # Write JSON
    with open(output_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info(f"Written metrics summary to {output_path}")


def write_fusions_gtf(
    fusions: List[Any],
    output_path: str
) -> None:
    """
    Write fusion junctions in GTF-like format.

    This is a wrapper function that provides GTF output for fusion junctions.
    For proper fusion detection, use fin.core.fusion_detector.FusionDetector
    which handles read extraction, re-mapping, and validation.

    Args:
        fusions: List of fusion junction objects
        output_path: Path to output GTF file
    """
    if not fusions:
        logger.warning("No fusions to write")
        return

    with open(output_path, 'w') as f:
        f.write("##gff-version 3\n")
        f.write("##source: FIN Fusion Detection\n")
        f.write("##feature_type: fusion\n\n")

        for fusion in fusions:
            # Basic GTF format for fusions
            # Extract attributes with safe defaults
            fusion_id = getattr(fusion, 'fusion_id', 'unknown')
            gene_5p = getattr(fusion, 'gene_5p', 'unknown')
            gene_3p = getattr(fusion, 'gene_3p', 'unknown')
            read_support = getattr(fusion, 'read_support', 0)

            attrs = f"fusion_id=\"{fusion_id}\";"
            attrs += f"gene_5p=\"{gene_5p}\";"
            attrs += f"gene_3p=\"{gene_3p}\";"
            attrs += f"read_support=\"{read_support}\""

            # Position information with safe defaults
            chrom = getattr(fusion, 'chrom_5p', 'chr1')
            start = getattr(fusion, 'pos_5p', 1)
            end = getattr(fusion, 'pos_3p', 100)
            score = getattr(fusion, 'confidence', 0.0)
            strand = getattr(fusion, 'strand', '+')

            # Format: chrom\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes
            f.write(f"{chrom}\tfin\tfusion\t{start}\t{end}\t{score:.2f}\t{strand}\t.\t{attrs}\n")

    logger.info(f"Written {len(fusions)} fusions to {output_path}")
