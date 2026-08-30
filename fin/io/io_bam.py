"""
BAM/SAM file format parser using pysam
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Iterator, Tuple, Any, Generator
from pathlib import Path
import logging

try:
    import pysam
except ImportError:
    raise ImportError(
        "pysam is required for BAM/SAM support. "
        "Install it with: pip install pysam"
    )


logger = logging.getLogger(__name__)


@dataclass
class RegionReadBatch:
    """Region fetch result with an explicit completeness contract.

    ``complete`` is True only when the region iterator was exhausted without
    an error and without an early ``max_reads`` stop. Consumers that derive
    HARD-NEGATIVE evidence ("no read supports X") must require
    ``complete=True``; a partial batch can only support positive statements.
    """

    reads: List[Dict[str, Any]] = field(default_factory=list)
    complete: bool = True
    error: Optional[str] = None


class BamReader:
    """Reader for BAM/SAM alignment files"""

    def __init__(self, file_path: str, index: bool = True):
        """
        Initialize BAM/SAM reader

        Args:
            file_path: Path to the BAM/SAM file
            index: Whether to load the index for random access
        """
        self.file_path = Path(file_path)
        if not self.file_path.exists():
            raise FileNotFoundError(f"BAM/SAM file not found: {file_path}")

        self._alignment_file = None
        self._index = index
        self._is_open = False

    def __enter__(self):
        """Context manager entry"""
        self.open()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        self.close()

    def open(self):
        """Open the alignment file"""
        try:
            mode = 'rb' if self.file_path.suffix.lower() == '.bam' else 'r'
            self._alignment_file = pysam.AlignmentFile(str(self.file_path), mode, index_filename=str(self.file_path) + '.bai' if self._index else None)
            self._is_open = True

            logger.info(f"Opened alignment file: {self.file_path}")
        except Exception as e:
            logger.error(f"Failed to open alignment file {self.file_path}: {e}")
            raise

    def close(self):
        """Close the alignment file"""
        if self._alignment_file is not None:
            self._alignment_file.close()
            self._alignment_file = None
            self._is_open = False
            logger.info(f"Closed alignment file: {self.file_path}")

    @property
    def is_open(self) -> bool:
        """Check if file is open"""
        return self._is_open

    @property
    def header(self) -> Dict[str, Any]:
        """
        Get file header

        Returns:
            Dictionary containing header information
        """
        if self._alignment_file is None:
            raise RuntimeError("Alignment file not opened. Call open() first.")
        return dict(self._alignment_file.header)

    @property
    def references(self) -> List[Tuple[str, int]]:
        """
        Get list of reference sequences

        Returns:
            List of tuples (reference_name, reference_length)
        """
        if self._alignment_file is None:
            raise RuntimeError("Alignment file not opened. Call open() first.")

        return [(self._alignment_file.get_reference_name(i), self._alignment_file.get_reference_length(i))
                for i in range(self._alignment_file.nreferences)]

    def count(self, region: Optional[str] = None, read_callback: str = 'all') -> int:
        """
        Count alignments

        Args:
            region: Region to count (e.g., 'chr1:1000-2000')
            read_callback: Read filter (all, mapped, unmapped)

        Returns:
            Number of alignments
        """
        if self._alignment_file is None:
            raise RuntimeError("Alignment file not opened. Call open() first.")

        try:
            return self._alignment_file.count(region=region, read_callback=read_callback)
        except Exception as e:
            logger.error(f"Failed to count alignments: {e}")
            return 0

    def fetch(self, reference: Optional[str] = None, start: Optional[int] = None,
              end: Optional[int] = None, region: Optional[str] = None,
              until_eof: bool = False) -> Iterator[pysam.AlignedSegment]:
        """
        Fetch alignments from a region

        Args:
            reference: Reference name
            start: Start position (0-based)
            end: End position
            region: Region string (e.g., 'chr1:1000-2000')
            until_eof: Continue to end of file

        Returns:
            Iterator over AlignedSegment objects
        """
        if self._alignment_file is None:
            raise RuntimeError("Alignment file not opened. Call open() first.")

        try:
            # Pass by keyword: pysam.fetch's 5th positional parameter is `tid`,
            # not `until_eof`. Passing until_eof positionally lands it in the tid
            # slot, and tid=False (==0) overrides `contig`, so every region fetch
            # silently returns reads from reference 0 regardless of `reference`.
            return self._alignment_file.fetch(
                contig=reference, start=start, stop=end,
                region=region, until_eof=until_eof,
            )
        except Exception as e:
            logger.error(f"Failed to fetch alignments: {e}")
            return iter([])

    def get_read(self, read_id: str) -> Optional[pysam.AlignedSegment]:
        """
        Get a specific read by name (requires indexed BAM)

        Args:
            read_id: Read name/ID

        Returns:
            AlignedSegment object or None if not found
        """
        if self._alignment_file is None:
            raise RuntimeError("Alignment file not opened. Call open() first.")

        try:
            # BAM files can be queried by read name if indexed properly
            # Note: pysam doesn't directly support fetch by name, so we iterate
            # For large files, this is inefficient - better to use fetch with region
            for read in self._alignment_file.fetch(until_eof=True):
                if read.query_name == read_id:
                    return read
            return None
        except Exception as e:
            logger.error(f"Failed to get read {read_id}: {e}")
            return None

    def complement(self,seq) -> str:
        
        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                         'N': 'N', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}
        try:
            return ''.join(complement_map[base] for base in seq)
        except KeyError as e:
            raise ValueError(f"Invalid DNA base: {e}")

    def reverse_complement(self,seq) -> str:
        if seq:
            return self.complement(seq)[::-1]
        else:
            return None

    def get_reference_sequence(self, alignment:pysam.AlignedSegment) -> str:
        if alignment.is_secondary:
            return None
        try:
            forward_ref_seq = alignment.get_reference_sequence().upper()
        except ValueError:
            # MD tag not present in BAM (minimap2 without --MD).
            # Downstream consumers default to "" so we can return None safely.
            return None
        return forward_ref_seq if not alignment.is_reverse else self.reverse_complement(forward_ref_seq)
    
    def alignment_to_dict(self, alignment: pysam.AlignedSegment) -> Dict[str, Any]:
        """
        Convert AlignedSegment to dictionary

        Args:
            alignment: AlignedSegment object

        Returns:
            Dictionary representation
        """
        if alignment is None:
            return {}

        return {
            'query_name': alignment.query_name,
            'flag': alignment.flag,
            'is_reverse': alignment.is_reverse,
            'is_secondary': alignment.is_secondary,
            'is_duplicate': alignment.is_duplicate,
            'is_supplementary': alignment.is_supplementary,
            'cigartuples': alignment.cigartuples,
            'reference_name': alignment.reference_name,
            'reference_id': alignment.reference_id,
            'reference_start': alignment.reference_start,
            'reference_end': alignment.reference_end,
            'query_alignment_start': alignment.query_alignment_start,
            'query_alignment_end': alignment.query_alignment_end,
            'query_alignment_length': alignment.query_alignment_length,
            'mapping_quality': alignment.mapping_quality,
            'template_length': alignment.template_length,
            'reference_sequence': self.get_reference_sequence(alignment),
            'is_forward': alignment.is_forward,
            'is_mapped': not alignment.is_unmapped,
            'query_length': alignment.query_length,
            'query_sequence': alignment.query_sequence,
        }

    def get_reads_in_region(self, region: str, max_reads: Optional[int] = None) -> List[Dict[str, Any]]:
        """
        Get all reads in a specific region

        Args:
            region: Region string (e.g., 'chr1:1000-2000')
            max_reads: Maximum number of reads to return

        Returns:
            List of read dictionaries. Legacy contract: a mid-iteration fetch
            error is swallowed and whatever was accumulated is returned, so
            the list may silently be INCOMPLETE. Evidence consumers that
            treat absence as zero support must use
            :meth:`get_reads_in_region_checked` instead.
        """
        return self.get_reads_in_region_checked(region, max_reads).reads

    def get_reads_in_region_checked(
        self, region: str, max_reads: Optional[int] = None
    ) -> RegionReadBatch:
        """Fetch region reads with an explicit completeness verdict.

        Returns a :class:`RegionReadBatch`. ``complete`` is False when the
        underlying iterator raised mid-fetch (``error`` records it) or when
        ``max_reads`` stopped the scan early — in both cases the batch must
        not be used as exhaustive (hard-negative) evidence.
        """
        if self._alignment_file is None:
            raise RuntimeError("Alignment file not opened. Call open() first.")

        batch = RegionReadBatch()
        try:
            for i, alignment in enumerate(self._alignment_file.fetch(region=region)):
                if max_reads is not None and i >= max_reads:
                    batch.complete = False
                    break
                batch.reads.append(self.alignment_to_dict(alignment))
        except Exception as e:
            logger.error(f"Failed to get reads in region {region}: {e}")
            batch.complete = False
            batch.error = f"{type(e).__name__}: {e}"
        return batch

    def iterate_alignments(self) -> Generator[Dict[str, Any], None, None]:
        """
        Generator to iterate through all alignments

        Yields:
            Alignment dictionary
        """
        if self._alignment_file is None:
            raise RuntimeError("Alignment file not opened. Call open() first.")

        try:
            for alignment in self._alignment_file.fetch(until_eof=True):
                yield self.alignment_to_dict(alignment)
        except Exception as e:
            logger.error(f"Failed to iterate alignments: {e}")

    def get_coverage(self, region: str) -> List[int]:
        """
        Get coverage array for a region

        Args:
            region: Region string (e.g., 'chr1:1000-2000')

        Returns:
            List of coverage values
        """
        if self._alignment_file is None:
            raise RuntimeError("Alignment file not opened. Call open() first.")

        try:
            ref_name, start_end = region.split(':')
            start, end = map(int, start_end.split('-'))

            coverage = [0] * (end - start)
            for alignment in self._alignment_file.fetch(ref_name, start, end):
                if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                    continue

                # Calculate overlap with region
                read_start = max(alignment.reference_start, start)
                read_end = min(alignment.reference_end, end)

                if read_start < read_end:
                    for pos in range(read_start - start, read_end - start):
                        if 0 <= pos < len(coverage):
                            coverage[pos] += 1

            return coverage

        except Exception as e:
            logger.error(f"Failed to get coverage for {region}: {e}")
            return []

    def get_mate(self, alignment: pysam.AlignedSegment) -> Optional[pysam.AlignedSegment]:
        """
        Get mate alignment for a paired-end read (requires indexed BAM)

        Args:
            alignment: AlignedSegment object

        Returns:
            Mate AlignedSegment or None if not found
        """
        if self._alignment_file is None:
            raise RuntimeError("Alignment file not opened. Call open() first.")

        if not alignment.is_proper_pair:
            return None

        try:
            mate = None
            for read in self._alignment_file.fetch(alignment.next_reference_name,
                                                  alignment.next_reference_start,
                                                  alignment.next_reference_start + 1):
                if (read.query_name == alignment.query_name and
                    ((read.is_read1 and alignment.is_read2) or (read.is_read2 and alignment.is_read1))):
                    mate = read
                    break
            return mate
        except Exception as e:
            logger.error(f"Failed to get mate: {e}")
            return None

    def get_file_stats(self) -> Dict[str, Any]:
        """
        Get file statistics

        Returns:
            Dictionary with file statistics
        """
        if self._alignment_file is None:
            raise RuntimeError("Alignment file not opened. Call open() first.")

        stats = {
            'filename': str(self.file_path),
            'file_format': 'BAM' if self.file_path.suffix.lower() == '.bam' else 'SAM',
            'nreferences': self._alignment_file.nreferences,
            'references': self.references,
            'header': dict(self._alignment_file.header),
        }

        # Count various types of reads
        try:
            stats['total_reads'] = self.count()
            stats['mapped_reads'] = self.count(read_callback='mapped')
            stats['unmapped_reads'] = self.count(read_callback='unmapped')
            stats['mapped_pct'] = (stats['mapped_reads'] / stats['total_reads'] * 100) if stats['total_reads'] > 0 else 0
        except Exception:
            stats['total_reads'] = 0
            stats['mapped_reads'] = 0
            stats['unmapped_reads'] = 0
            stats['mapped_pct'] = 0

        return stats

    @staticmethod
    def sort_and_index(bam_path: str, output_path: Optional[str] = None) -> bool:
        """
        Sort and index a BAM file

        Args:
            bam_path: Path to input BAM file
            output_path: Path to output sorted BAM file (default: input + '.sorted.bam')

        Returns:
            True if successful, False otherwise
        """
        try:
            input_path = Path(bam_path)
            output = Path(output_path) if output_path else input_path.with_suffix('.sorted.bam')

            # Sort the BAM file
            pysam.sort('-o', str(output), str(input_path))

            # Index the sorted BAM file
            pysam.index(str(output))

            logger.info(f"Successfully sorted and indexed BAM: {output}")
            return True

        except Exception as e:
            logger.error(f"Failed to sort and index BAM: {e}")
            return False
