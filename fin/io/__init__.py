"""
I/O module for fin.

This module handles reading and writing of various file formats:
- Bedfiles and BAM files (mapping results)
- Fast5/Blow5/Pod5 files (nanopore raw signals)
- GTF/GFF annotation files
- Output formats for RNA modification and isoform information
"""

from .io_bam import BamParser, parse_bed_file, ReadAlignment
from .io_signal import Fast5Parser, Pod5Parser, ReadSignal, get_signal_parser
from .sequence_extractor import (
    ReferenceExtractor, BamSequenceExtractor, GenomicRegion,
    reverse_complement, extract_sequence_from_bed_and_fasta,
    get_sequence_at_position
)
from .io_gtf import GtfParser, TranscriptFeature, parse_gtf_attributes
from .output_writer import (
    write_gtf, write_bed12, write_read_assignments_tsv, write_metrics_summary
)

__all__ = [
    "io_bam",
    "io_signal",
    "io_output",
    "BamParser",
    "ReadAlignment",
    "Fast5Parser",
    "Pod5Parser",
    "ReadSignal",
    "get_signal_parser",
    "ReferenceExtractor",
    "BamSequenceExtractor",
    "GenomicRegion",
    "reverse_complement",
    "extract_sequence_from_bed_and_fasta",
    "get_sequence_at_position",
]
