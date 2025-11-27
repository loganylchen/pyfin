"""
Unified I/O Manager for handling all input file types.

Coordinates reading from multiple file formats and provides a unified interface
for retrieving read alignments and corresponding signal data.
"""

import os
from typing import Dict, List, Optional, Iterator, Tuple, Union
import logging
from pathlib import Path

from .io_bam import BamParser, ReadAlignment
from .io_signal import (
    Fast5Parser, Pod5Parser, ReadSignal,
    get_signal_parser, is_fast5_file, is_pod5_file
)
from ..utils.sequences import reverse_complement

logger = logging.getLogger(__name__)


class InputConfig:
    """Configuration for input files."""

    def __init__(
        self,
        bam_file: str,
        signal_dir: Optional[str] = None,
        signal_files: Optional[List[str]] = None
    ):
        """
        Initialize input configuration.

        Args:
            bam_file: Path to BAM alignment file
            signal_dir: Directory containing signal files (fast5/pod5/blow5)
            signal_files: List of explicit signal file paths
        """
        self.bam_file = bam_file
        self.signal_dir = signal_dir
        self.signal_files = signal_files or []

        # Auto-discover signal files if directory provided
        if signal_dir and not signal_files:
            self._discover_signal_files()

    def _discover_signal_files(self):
        """Auto-discover signal files in directory."""
        if not self.signal_dir:
            return

        supported_extensions = {
            '.fast5', '.pod5', '.slow5', '.blow5'
        }

        for ext in supported_extensions:
            self.signal_files.extend(
                str(p) for p in Path(self.signal_dir).rglob(f"*{ext}")
            )

        logger.info(f"Discovered {len(self.signal_files)} signal files in {self.signal_dir}")


class ReadData:
    """
    Container for all data related to a single read.

    Combines alignment information, raw signal, and sequence data.
    """

    def __init__(
        self,
        read_id: str,
        alignment: Optional[ReadAlignment] = None,
        signal: Optional[ReadSignal] = None,
        sequence: Optional[str] = None
    ):
        """
        Initialize read data container.

        Args:
            read_id: Unique read identifier
            alignment: Read alignment information
            signal: Raw nanopore signal
            sequence: Basecalled sequence
        """
        self.read_id = read_id
        self.alignment = alignment
        self.signal = signal
        self.sequence = sequence

        # Derived properties
        self.is_valid = alignment is not None and signal is not None
        self.contig = alignment.contig if alignment else None
        self.ref_start = alignment.ref_start if alignment else None
        self.ref_end = alignment.ref_end if alignment else None
        self.direction = alignment.is_reverse if alignment else None

    def get_reference_coords(self, contig: str, start: int, end: int):
        """
        Get reference coordinates for this read.

        Args:
            contig: Reference contig name
            start: Start position
            end: End position

        Returns:
            Tuple of (ref_start, ref_end) or None if not aligned
        """
        if not self.alignment:
            return None

        if self.alignment.contig != contig:
            return None

        return self.alignment.ref_start, self.alignment.ref_end

    def get_orientation_strand(self, is_rna: bool = False) -> str:
        """
        Get orientation strand considering RNA coordinate issues.

        Corectly handles RNA 5'→3' sequences with 3'→5' nanopore signals.

        Args:
            is_rna: Whether this is RNA data

        Returns:
            Strand orientation string
        """
        if not self.alignment:
            return "UNKNOWN"

        if is_rna:
            # For RNA, BAM strand indicates original sequence orientation
            # Signal is generated in reverse (3'→5')
            return "-" if self.alignment.is_reverse else "+"
        else:
            # For DNA, straightforward
            return "-" if self.alignment.is_reverse else "+"

    def get_sequence_for_signal(self, is_rna: bool = False) -> Optional[str]:
        """
        Get basecalled sequence matching signal orientation.

        For RNA, if the read is on '+' strand in BAM, the signal is generated
        from 3'→5' and the sequence needs to be reverse complemented.

        Args:
            is_rna: Whether this is RNA data

        Returns:
            Sequence in signal coordinate system
        """
        if not self.sequence:
            return None

        seq = self.sequence

        # Apply strand correction for signal matching
        if is_rna and self.alignment:
            if self.alignment.is_reverse:
                # Already reverse, no change needed
                pass
            else:
                # Forward strand in BAM means signal is reverse complement
                seq = reverse_complement(seq)

        return seq


class IOManager:
    """
    Unified I/O manager for handling all input files.

    Manages multiple BAM files (native and IVT) and corresponding signal files,
    providing a unified interface to retrieve read data.
    """

    def __init__(
        self,
        native_input: Union[InputConfig, str],
        ivt_input: Optional[Union[InputConfig, str]] = None,
        is_rna: bool = True
    ):
        """
        Initialize I/O manager.

        Args:
            native_input: Input config for native RNA (or path string)
            ivt_input: Input config for IVT RNA (or path string, or None)
            is_rna: Whether this is RNA data (affects strand handling)
        """
        # Convert string paths to InputConfig
        if isinstance(native_input, str):
            native_input = InputConfig(bam_file=native_input)

        if isinstance(ivt_input, str):
            ivt_input = InputConfig(bam_file=ivt_input)

        self.native_input = native_input
        self.ivt_input = ivt_input
        self.is_rna = is_rna

        # Parsers
        self.native_bam_parser = None
        self.ivt_bam_parser = None

        # Signal file indices
        self.signal_file_index = {}  # read_id -> signal_file_path

        # Build indices
        self._build_signal_index()

    def __enter__(self):
        """Context manager entry."""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit."""
        self.close()

    def _build_signal_index(self):
        """Build index mapping read_id to signal file."""
        logger.info("Building signal file index...")

        all_signal_files = []

        # Collect from native input
        if self.native_input:
            all_signal_files.extend(self.native_input.signal_files)

        # Collect from IVT input
        if self.ivt_input:
            all_signal_files.extend(self.ivt_input.signal_files)

        # Index all signal files
        total_reads = 0
        for signal_file in all_signal_files:
            try:
                parser_class = get_signal_parser(signal_file)

                with parser_class(signal_file) as parser:
                    read_ids = parser.get_read_ids()

                    for read_id in read_ids:
                        self.signal_file_index[read_id] = signal_file
                        total_reads += 1

                logger.debug(f"Indexed {len(read_ids)} reads from {signal_file}")

            except Exception as e:
                logger.error(f"Failed to index {signal_file}: {e}")

        logger.info(f"Signal file index built with {total_reads} reads from {len(all_signal_files)} files")

    def open(self):
        """Open BAM files."""
        # Open native BAM
        if self.native_input:
            self.native_bam_parser = BamParser(self.native_input.bam_file)
            self.native_bam_parser.open()

        # Open IVT BAM
        if self.ivt_input:
            self.ivt_bam_parser = BamParser(self.ivt_input.bam_file)
            self.ivt_bam_parser.open()

    def close(self):
        """Close all open files."""
        if self.native_bam_parser:
            self.native_bam_parser.close()

        if self.ivt_bam_parser:
            self.ivt_bam_parser.close()

    def get_read_signal(self, read_id: str) -> Optional[ReadSignal]:
        """
        Get raw signal for a specific read.

        Args:
            read_id: Read identifier

        Returns:
            ReadSignal object or None if not found
        """
        if read_id not in self.signal_file_index:
            logger.warning(f"No signal file for read {read_id}")
            return None

        signal_file = self.signal_file_index[read_id]
        parser_class = get_signal_parser(signal_file)

        with parser_class(signal_file) as parser:
            return parser.get_read_signal(read_id)

    def iterate_reads(
        self,
        sample_type: str = "native"
    ) -> Iterator[ReadData]:
        """
        Iterate over all reads, combining alignment and signal.

        Args:
            sample_type: "native" or "ivt"

        Yields:
            ReadData objects
        """
        # Select parser based on sample type
        if sample_type == "native":
            bam_parser = self.native_bam_parser
        elif sample_type == "ivt" and self.ivt_bam_parser:
            bam_parser = self.ivt_bam_parser
        else:
            raise ValueError(f"Invalid sample type: {sample_type}")

        if bam_parser is None:
            logger.warning(f"No {sample_type} BAM parser available")
            return

        # Iterate alignments
        for alignment in bam_parser.parse_alignments():
            # Get signal for this read
            signal = self.get_read_signal(alignment.read_id)

            # Skip if no signal
            if signal is None:
                logger.debug(f"No signal for read {alignment.read_id}, skipping")
                continue

            # Get sequence if available
            sequence = alignment.query_sequence

            # Apply strand correction for RNA if needed
            if self.is_rna and sequence and alignment.is_reverse:
                # For RNA on reverse strand, signal and sequence are aligned
                pass
            elif self.is_rna and sequence and not alignment.is_reverse:
                # For RNA on forward strand, signal is reverse complement
                pass

            yield ReadData(
                read_id=alignment.read_id,
                alignment=alignment,
                signal=signal,
                sequence=sequence
            )

    def get_region_reads(self, contig: str, start: int, end: int, sample_type: str = "native") -> List[ReadData]:
        """
        Get all reads covering a specific genomic region.

        Args:
            contig: Contig/chromosome name
            start: Start position (0-based)
            end: End position (exclusive)
            sample_type: "native" or "ivt"

        Returns:
            List of ReadData objects
        """
        bam_parser = self.native_bam_parser if sample_type == "native" else self.ivt_bam_parser

        if bam_parser is None:
            return []

        reads = []
        for alignment in bam_parser.fetch_region(contig, start, end):
            signal = self.get_read_signal(alignment.read_id)

            if signal:
                reads.append(ReadData(
                    read_id=alignment.read_id,
                    alignment=alignment,
                    signal=signal,
                    sequence=alignment.query_sequence
                ))

        return reads

    def get_comparison_pairs(self, contig: str, start: int, end: int) -> List[Tuple[ReadData, List[ReadData]]]:
        """
        Get native reads and corresponding IVT reads for comparison.

        For each native read covering a region, find IVT reads mapped to the same region.

        Args:
            contig: Contig/chromosome
            start: Start position
            end: End position

        Returns:
            List of (native_read, [ivt_reads...]) tuples
        """
        pairs = []

        # Get native reads
        native_reads = self.get_region_reads(contig, start, end, sample_type="native")

        # Get IVT reads for the same region
        ivt_reads = self.get_region_reads(contig, start, end, sample_type="ivt")

        # Group by position
        for native_read in native_reads:
            # Find IVT reads that overlap this read
            overlapping_ivt = []
            native_start, native_end = native_read.ref_start, native_read.ref_end

            for ivt_read in ivt_reads:
                # Check overlap
                if (ivt_read.ref_start < native_end) and (ivt_read.ref_end > native_start):
                    overlapping_ivt.append(ivt_read)

            if overlapping_ivt:
                pairs.append((native_read, overlapping_ivt))

        logger.info(f"Found {len(pairs)} native reads with {len(ivt_reads)} IVT reads")

        return pairs
