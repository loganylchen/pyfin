"""
FASTQ file format parser for basecalled reads
"""

from typing import Iterator, Tuple, Optional
from pathlib import Path
import logging


logger = logging.getLogger(__name__)


class FastqReader:
    """Reader for FASTQ files containing basecalled reads"""

    def __init__(self, file_path: str):
        """
        Initialize FASTQ reader

        Args:
            file_path: Path to the FASTQ file
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"FASTQ file not found: {file_path}")

        self._file_handle = None
        self._is_open = False

    def __enter__(self):
        """Context manager entry"""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()

    def open(self):
        """Open the FASTQ file"""
        try:
            self._file_handle = open(self.file_path, 'r')
            self._is_open = True
            logger.info(f"Opened FASTQ file: {self.file_path}")
        except Exception as e:
            logger.error(f"Failed to open FASTQ file {self.file_path}: {e}")
            raise

    def close(self):
        """Close the FASTQ file"""
        if self._file_handle is not None:
            self._file_handle.close()
            self._file_handle = None
            self._is_open = False
            logger.info(f"Closed FASTQ file: {self.file_path}")

    @property
    def is_open(self) -> bool:
        """Check if file is open"""
        return self._is_open

    def iter_reads(self) -> Iterator[Tuple[str, str, str]]:
        """
        Iterate through all reads in the FASTQ file

        Yields:
            Tuple of (read_id, sequence, quality_scores)
        """
        if not self._is_open:
            raise RuntimeError("FASTQ file not opened. Call open() first.")

        line_number = 0
        while True:
            # Read header line (starts with @)
            header = self._file_handle.readline()
            line_number += 1
            if not header:
                # End of file
                break
            if not header.startswith('@'):
                logger.warning(f"Invalid FASTQ format at line {line_number}: expected header starting with '@'")
                continue

            # Read sequence line
            sequence = self._file_handle.readline().strip()
            line_number += 1
            if not sequence:
                logger.warning(f"Incomplete FASTQ record at line {line_number}")
                break

            # Read plus line (starts with +)
            plus = self._file_handle.readline()
            line_number += 1
            if not plus.startswith('+'):
                logger.warning(f"Invalid FASTQ format at line {line_number}: expected '+'")
                continue

            # Read quality line
            quality = self._file_handle.readline().strip()
            line_number += 1
            if not quality:
                logger.warning(f"Incomplete FASTQ record at line {line_number}")
                break

            # Extract read ID from header (remove @ and take first part)
            read_id = header[1:].split()[0]

            yield read_id, sequence, quality

    def get_read_by_id(self, target_id: str) -> Optional[Tuple[str, str, str]]:
        """
        Get a specific read by ID

        Args:
            target_id: Read ID to search for

        Returns:
            Tuple of (read_id, sequence, quality) or None if not found
        """
        if not self._is_open:
            raise RuntimeError("FASTQ file not opened. Call open() first.")

        # Store current position
        current_pos = self._file_handle.tell()

        # Reset to beginning of file
        self._file_handle.seek(0)

        for read_id, sequence, quality in self.iter_reads():
            if read_id == target_id:
                # Restore position and return
                self._file_handle.seek(current_pos)
                return read_id, sequence, quality

        # Not found - restore position
        self._file_handle.seek(current_pos)
        return None

    def get_all_reads(self) -> dict:
        """
        Load all reads into a dictionary

        Returns:
            Dictionary mapping read_id to (sequence, quality) tuples
        """
        reads = {}
        for read_id, sequence, quality in self.iter_reads():
            reads[read_id] = (sequence, quality)
        return reads

    def count_reads(self) -> int:
        """
        Count the number of reads in the file

        Returns:
            Number of reads
        """
        if not self._is_open:
            raise RuntimeError("FASTQ file not opened. Call open() first.")

        # Store current position
        current_pos = self._file_handle.tell()

        # Reset to beginning of file
        self._file_handle.seek(0)

        count = sum(1 for _ in self.iter_reads())

        # Restore position
        self._file_handle.seek(current_pos)

        return count

    def get_file_stats(self) -> dict:
        """
        Get file statistics

        Returns:
            Dictionary with file information
        """
        if not self._is_open:
            raise RuntimeError("FASTQ file not opened. Call open() first.")

        stats = {
            'filename': str(self.file_path),
            'num_reads': self.count_reads(),
        }

        # Get read length statistics from first 1000 reads
        lengths = []
        for i, (_, sequence, _) in enumerate(self.iter_reads()):
            if i >= 1000:
                break
            lengths.append(len(sequence))

        if lengths:
            stats['min_length'] = min(lengths)
            stats['max_length'] = max(lengths)
            stats['mean_length'] = sum(lengths) / len(lengths)

        return stats
