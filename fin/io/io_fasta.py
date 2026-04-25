"""
FASTA file format parser
Supports reading and writing FASTA files with sequence records
"""

from typing import List, Dict, Optional, Generator, Union
from pathlib import Path
import logging

try:
    import pysam
except ImportError:
    raise ImportError(
        "pysam is required for FASTA/FASTQ support. "
        "Install it with: pip install pysam"
    )


logger = logging.getLogger(__name__)


class FASTARecord:
    """Represents a single FASTA record"""

    def __init__(self, header: str, sequence: str):
        """
        Initialize FASTA record

        Args:
            header: FASTA header line (without the '>')
            sequence: DNA/RNA/protein sequence
        """
        self.header = header
        self.sequence = sequence
        self.id = self._extract_id()
        self.description = self._extract_description()

    def _extract_id(self) -> str:
        """Extract sequence ID from header (first word)"""
        return self.header.split()[0] if self.header else ""

    def _extract_description(self) -> str:
        """Extract description from header (everything after first word)"""
        parts = self.header.split(maxsplit=1)
        return parts[1] if len(parts) > 1 else ""

    @property
    def length(self) -> int:
        """Get sequence length"""
        return len(self.sequence)

    @property
    def gc_content(self) -> float:
        """Calculate GC content (0-100)"""
        if not self.sequence:
            return 0.0
        gc_count = self.sequence.upper().count('G') + self.sequence.upper().count('C')
        return (gc_count / len(self.sequence)) * 100

    @property
    def sequence_upper(self) -> str:
        """Get sequence in uppercase"""
        return self.sequence.upper()

    @property
    def sequence_lower(self) -> str:
        """Get sequence in lowercase"""
        return self.sequence.lower()

    def __str__(self) -> str:
        """String representation (FASTA format)"""
        return f">{self.header}\n{self._format_sequence()}"

    def _format_sequence(self, width: int = 60) -> str:
        """Format sequence with line breaks"""
        return '\n'.join(self.sequence[i:i+width] for i in range(0, len(self.sequence), width))

    def subsequence(self, start: int, end: int) -> str:
        """
        Extract subsequence (0-based, end exclusive)

        Args:
            start: Start position
            end: End position

        Returns:
            Subsequence string
        """
        return self.sequence[start:end]

    def complement(self) -> str:
        """
        Get complement of DNA sequence

        Returns:
            Complement sequence

        Raises:
            ValueError: If sequence contains invalid DNA characters
        """
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                         'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
        try:
            return ''.join(complement_map[base] for base in self.sequence)
        except KeyError as e:
            raise ValueError(f"Invalid DNA base: {e}")

    def reverse_complement(self) -> str:
        """
        Get reverse complement of DNA sequence

        Returns:
            Reverse complement sequence
        """
        return self.complement()[::-1]

    def count_nucleotides(self) -> Dict[str, int]:
        """
        Count occurrences of each nucleotide

        Returns:
            Dictionary mapping nucleotides to counts
        """
        counts = {}
        for base in self.sequence.upper():
            counts[base] = counts.get(base, 0) + 1
        return counts


class FASTAReader:
    """Reader for FASTA files using pysam"""

    def __init__(self, file_path: str):
        """
        Initialize FASTA reader

        Args:
            file_path: Path to the FASTA file
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"FASTA file not found: {file_path}")

        self._fasta_file = None
        self._record_count = None

    def __enter__(self):
        """Context manager entry"""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()

    def open(self):
        """Open the FASTA file"""
        try:
            self._fasta_file = pysam.FastxFile(str(self.file_path))
            logger.info(f"Opened FASTA file: {self.file_path}")
        except Exception as e:
            logger.error(f"Failed to open FASTA file {self.file_path}: {e}")
            raise

    def close(self):
        """Close the FASTA file"""
        if self._fasta_file is not None:
            self._fasta_file.close()
            self._fasta_file = None
            logger.info(f"Closed FASTA file: {self.file_path}")

    def _parse_records(self) -> Generator[FASTARecord, None, None]:
        """
        Generator to parse records from file

        Yields:
            FASTARecord objects
        """
        if self._fasta_file is None:
            raise RuntimeError("FASTA file not opened. Call open() first.")

        try:
            for entry in self._fasta_file:
                yield FASTARecord(entry.name, entry.sequence)
        except Exception as e:
            logger.error(f"Error reading FASTA file: {e}")
            raise

    def get_records(self) -> List[FASTARecord]:
        """
        Get all records from the file

        Returns:
            List of FASTARecord objects
        """
        records = []
        # Create a new FastxFile iterator since we can't reuse it
        try:
            with pysam.FastxFile(str(self.file_path)) as fasta_file:
                for entry in fasta_file:
                    records.append(FASTARecord(entry.name, entry.sequence))
        except Exception as e:
            logger.error(f"Error reading FASTA file: {e}")
            raise
        return records

    def iterate_records(self) -> Generator[FASTARecord, None, None]:
        """
        Generator to iterate through all records

        Yields:
            FASTARecord objects
        """
        yield from self._parse_records()

    def get_record_by_id(self, record_id: str) -> Optional[FASTARecord]:
        """
        Get a specific record by ID

        Args:
            record_id: Sequence ID to search for

        Returns:
            FASTARecord object or None if not found
        """
        for record in self.iterate_records():
            if record.id == record_id:
                return record
        return None

    def get_records_by_ids(self, record_ids: List[str]) -> List[FASTARecord]:
        """
        Get multiple records by IDs

        Args:
            record_ids: List of sequence IDs to search for

        Returns:
            List of FASTARecord objects
        """
        id_set = set(record_ids)
        records = []
        for record in self.iterate_records():
            if record.id in id_set:
                records.append(record)
        return records

    def get_record_count(self) -> int:
        """
        Get total number of records

        Returns:
            Number of records
        """
        if self._record_count is None:
            self._record_count = sum(1 for _ in self.iterate_records())
        return self._record_count

    def get_sequence_ids(self) -> List[str]:
        """
        Get list of sequence IDs

        Returns:
            List of sequence IDs
        """
        return [record.id for record in self.iterate_records()]

    def get_total_length(self) -> int:
        """
        Get total length of all sequences

        Returns:
            Total number of bases across all records
        """
        return sum(record.length for record in self.iterate_records())

    def get_summary(self) -> Dict[str, Union[int, float]]:
        """
        Get summary statistics for the FASTA file

        Returns:
            Dictionary with summary statistics
        """
        # Need to reopen file since pysam.FastxFile doesn't support seeking
        records = self.get_records()
        if not records:
            return {"num_records": 0, "total_length": 0, "avg_length": 0, "min_length": 0, "max_length": 0}

        lengths = [record.length for record in records]
        return {
            "num_records": len(records),
            "total_length": sum(lengths),
            "avg_length": sum(lengths) / len(lengths),
            "min_length": min(lengths),
            "max_length": max(lengths)
        }

    @staticmethod
    def write_records(records: List[FASTARecord], output_path: str, mode: str = 'w'):
        """
        Write records to a FASTA file

        Args:
            records: List of FASTARecord objects
            output_path: Output file path
            mode: Write mode ('w' for write, 'a' for append)
        """
        try:
            with open(output_path, mode) as f:
                for record in records:
                    f.write(str(record) + '\n')

            logger.info(f"Wrote {len(records)} records to {output_path}")

        except Exception as e:
            logger.error(f"Failed to write FASTA file {output_path}: {e}")
            raise

    @staticmethod
    def create_record(sequence_id: str, sequence: str, description: str = "") -> FASTARecord:
        """
        Create a new FASTA record

        Args:
            sequence_id: Sequence ID
            sequence: Sequence string
            description: Optional description

        Returns:
            FASTARecord object
        """
        header = sequence_id if not description else f"{sequence_id} {description}"
        return FASTARecord(header, sequence)
