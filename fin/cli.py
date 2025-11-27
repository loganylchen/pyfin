"""
Command-line interface for FIN isoform detection pipeline.

Provides a CLI tool for detecting and validating RNA isoforms from
nanopore direct RNA sequencing data using signal-level analysis
(eventalign) and DTW clustering.
"""

import click
import logging
import sys
from pathlib import Path
from typing import Tuple, Optional, Dict, List
import numpy as np

from fin.core.isoform_detector import IsoformDetector
from fin.utils.config import PipelineConfig, ConfigManager
from fin.io.output_writer import write_gtf, write_bed12, write_read_assignments_tsv
from fin.core.isoform import ValidatedIsoform
from fin.io.io_gtf import GtfParser

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def validate_input_files(
    bam_file: str,
    reference: str,
    gtf_file: Optional[str] = None,
    signal_dir: Optional[str] = None
) -> bool:
    """
    Validate that all input files exist and are readable.

    Args:
        bam_file: Path to BAM file
        reference: Path to reference FASTA
        gtf_file: Path to GTF file (optional)
        signal_dir: Path to signal directory (optional)

    Returns:
        True if all files are valid
    """
    errors = []

    if not Path(bam_file).exists():
        errors.append(f"BAM file not found: {bam_file}")

    if not Path(reference).exists():
        errors.append(f"Reference file not found: {reference}")

    if gtf_file and not Path(gtf_file).exists():
        errors.append(f"GTF file not found: {gtf_file}")

    if signal_dir and not Path(signal_dir).is_dir():
        errors.append(f"Signal directory not found: {signal_dir}")

    if errors:
        for error in errors:
            logger.error(error)
        return False

    return True


def parse_gene_region(gene: str, gtf_file: Optional[str] = None) -> Tuple[str, int, int]:
    """
    Parse gene region from string.

    Supports two formats:
    1. "chr:start-end" (e.g., "chr17:7668420-7687960")
    2. Gene name (e.g., "TP53") - requires GTF file

    Args:
        gene: Gene string
        gtf_file: Path to GTF file (required for gene name format)

    Returns:
        Tuple of (chrom, start, end)

    Raises:
        ValueError: If region cannot be parsed
    """
    if ':' in gene:
        # Region format: chr:start-end
        chrom, coords = gene.split(':', 1)
        if '-' not in coords:
            raise ValueError(f"Invalid region format: {gene}. Expected chr:start-end")

        start_str, end_str = coords.split('-', 1)
        try:
            start = int(start_str.replace(',', ''))
            end = int(end_str.replace(',', ''))
        except ValueError:
            raise ValueError(f"Invalid coordinates in: {gene}")

        if start >= end:
            raise ValueError(f"Start >= end in: {gene}")

        logger.info(f"Parsed region: {chrom}:{start}-{end}")
        return chrom, start, end

    else:
        # Gene name format
        if not gtf_file:
            raise ValueError(f"GTF file required to parse gene name: {gene}")

        logger.info(f"Looking up gene '{gene}' in GTF...")

        try:
            parser = GtfParser(gtf_file)
            gene_transcripts = parser.get_all_genes().get(gene)

            if not gene_transcripts:
                raise ValueError(f"Gene '{gene}' not found in GTF")

            # Get bounds of all transcripts
            chrom = gene_transcripts[0].chrom
            start = min(min(e[0] for e in t.exons) for t in gene_transcripts)
            end = max(max(e[1] for e in t.exons) for t in gene_transcripts)

            logger.info(f"Found gene '{gene}' at {chrom}:{start}-{end}")
            return chrom, start, end

        except Exception as e:
            raise ValueError(f"Error looking up gene '{gene}': {e}")


@click.command()
@click.option(
    '--config', '-c',
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help='Path to configuration YAML file'
)
@click.option(
    '--gene', '-g',
    type=str,
    required=True,
    help='Gene name or region (e.g., "TP53" or "chr17:7668420-7687960")'
)
@click.option(
    '--bam', '-b',
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help='Path to BAM alignment file'
)
@click.option(
    '--reference', '-r',
    type=click.Path(exists=True, path_type=Path),
    required=True,
    help='Path to reference genome FASTA file'
)
@click.option(
    '--gtf',
    type=click.Path(exists=True, path_type=Path),
    help='Path to GTF annotation file (optional but recommended)'
)
@click.option(
    '--signal-dir', '-s',
    type=click.Path(exists=True, path_type=Path),
    help='Directory containing nanopore signal files (FAST5/POD5/SLOW5)'
)
@click.option(
    '--output-dir', '-o',
    type=click.Path(path_type=Path),
    default=Path('./outputs'),
    help='Output directory (default: ./outputs)'
)
@click.option(
    '--prefix', '-p',
    type=str,
    default='fin',
    help='Output file prefix (default: fin)'
)
@click.option(
    '--verbose', '-v',
    is_flag=True,
    help='Enable verbose logging'
)
def main(
    config: Path,
    gene: str,
    bam: Path,
    reference: Path,
    gtf: Optional[Path],
    signal_dir: Optional[Path],
    output_dir: Path,
    prefix: str,
    verbose: bool
):
    """
    FIN: Isoform Detection Pipeline for Nanopore Direct RNA Sequencing

    Detects and validates RNA isoforms by integrating signal-level analysis
    (eventalign completeness) with sequence-based validation (DTW clustering).

    Input:
    - BAM file with read alignments
    - GTF file with reference annotations (optional)
    - Signal files (FAST5/POD5/SLOW5) for signal analysis
    - Reference genome FASTA

    Output:
    - GTF file with validated isoforms
    - BED12 file with read-to-isoform assignments

    Examples:

    # Run on specific region with all data:
    fin --config config.yaml \
        --gene chr17:7668420-7687960 \
        --bam alignments.bam \
        --reference reference.fasta \
        --gtf annotations.gtf \
        --signal-dir signals/

    # Run on gene name (requires GTF):
    fin --config config.yaml \
        --gene TP53 \
        --bam alignments.bam \
        --reference reference.fasta \
        --gtf annotations.gtf

    # Run without signal files (sequence only):
    fin --config config.yaml \
        --gene chr17:7668420-7687960 \
        --bam alignments.bam \
        --reference reference.fasta
    """
    # Configure logging
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logger.info("=" * 60)
    logger.info("FIN: Nanopore RNA Isoform Detection Pipeline")
    logger.info("=" * 60)

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Output directory: {output_dir}")

    # Validate input files
    logger.info("Validating input files...")
    if not validate_input_files(
        str(bam),
        str(reference),
        str(gtf) if gtf else None,
        str(signal_dir) if signal_dir else None
    ):
        logger.error("Input validation failed")
        sys.exit(1)

    # Load configuration
    try:
        logger.info(f"Loading configuration: {config}")
        pipeline_config = ConfigManager.load_config(str(config))

        # Validate configuration
        is_valid, message = ConfigManager.validate_config(pipeline_config)
        if not is_valid:
            logger.error(f"Configuration validation failed: {message}")
            sys.exit(1)

        logger.info("Configuration loaded and validated")

        # Log key parameters
        logger.info(f"  min_completeness_threshold: {pipeline_config.algorithm.min_completeness_threshold}")
        logger.info(f"  dtw_similarity_threshold: {pipeline_config.algorithm.dtw_similarity_threshold}")
        logger.info(f"  min_read_support: {pipeline_config.algorithm.min_read_support}")
        logger.info(f"  use_cuda: {pipeline_config.dtw.use_cuda}")

    except Exception as e:
        logger.error(f"Failed to load configuration: {e}")
        sys.exit(1)

    # Parse gene region
    try:
        logger.info(f"Parsing gene/region: {gene}")
        chrom, start, end = parse_gene_region(gene, str(gtf) if gtf else None)
        logger.info(f"Processing: {chrom}:{start:,}-{end:,} ({end - start:,} bp)")
    except Exception as e:
        logger.error(f"Failed to parse gene/region: {e}")
        sys.exit(1)

    # Initialize detector
    try:
        logger.info("Initializing isoform detector...")
        detector = IsoformDetector(
            config=pipeline_config,
            reference_fasta=str(reference),
            signal_dir=str(signal_dir) if signal_dir else None,
            use_cuda=pipeline_config.dtw.use_cuda
        )
        logger.info("Detector initialized successfully")
    except Exception as e:
        logger.error(f"Failed to initialize detector: {e}")
        sys.exit(1)

    # Run detection
    try:
        logger.info("Running isoform detection...")
        logger.info("-" * 60)

        validated_isoforms, read_assignments = detector.process_gene_region(
            chrom=chrom,
            start=start,
            end=end,
            bam_file=str(bam),
            gtf_file=str(gtf) if gtf else None
        )

        logger.info("-" * 60)
        logger.info("Detection completed successfully")

    except Exception as e:
        logger.error(f"Error during detection: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    # Calculate statistics
    total_isoforms = len(validated_isoforms)
    novel_isoforms = sum(1 for iso in validated_isoforms if iso.is_novel)
    annotated_isoforms = total_isoforms - novel_isoforms
    total_assigned_reads = len(read_assignments)

    # Print results summary
    logger.info("\n" + "=" * 60)
    logger.info("RESULTS SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Processed region: {chrom}:{start:,}-{end:,}")
    logger.info(f"Validated isoforms: {total_isoforms}")
    logger.info(f"  - Annotated: {annotated_isoforms}")
    logger.info(f"  - Novel: {novel_isoforms}")
    logger.info(f"Assigned reads: {total_assigned_reads}")

    if total_isoforms > 0:
        avg_support = np.mean([iso.read_support for iso in validated_isoforms])
        avg_completeness = np.mean([iso.avg_completeness for iso in validated_isoforms])
        logger.info(f"Average read support: {avg_support:.1f}")
        logger.info(f"Average completeness: {avg_completeness:.3f}")

    # Generate output file paths
    gtf_output = output_dir / f"{prefix}.isoforms.gtf"
    bed_output = output_dir / f"{prefix}.reads.bed12"
    tsv_output = output_dir / f"{prefix}.assignments.tsv"

    # Write outputs
    try:
        logger.info("\nWriting outputs...")

        if total_isoforms > 0:
            logger.info(f"  - GTF: {gtf_output}")
            write_gtf(
                validated_isoforms,
                str(gtf_output),
                source="fin"
            )

            if read_assignments:
                logger.info(f"  - BED12: {bed_output}")
                write_bed12(
                    read_assignments,
                    str(bam),
                    validated_isoforms,
                    str(bed_output)
                )

                logger.info(f"  - TSV: {tsv_output}")
                write_read_assignments_tsv(
                    read_assignments,
                    {},
                    str(tsv_output)
                )

        logger.info("Outputs written successfully")

    except Exception as e:
        logger.error(f"Error writing outputs: {e}")
        sys.exit(1)

    # Print completion message
    logger.info("\n" + "=" * 60)
    logger.info("✓ FIN pipeline completed successfully!")
    logger.info("=" * 60)

    sys.exit(0)


if __name__ == '__main__':
    main()
