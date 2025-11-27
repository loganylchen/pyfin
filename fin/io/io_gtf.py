"""
GTF/GFF annotation file parser.

Provides functionality to parse GTF/GFF files and extract transcript
structures (isoforms) for use in isoform detection pipeline.

Handles:
- GTF format (http://mblab.wustl.edu/GTF2.html)
- GFF3 format
- Exon grouping by transcript_id or Parent attributes
- Transcript coordinate extraction
- Strand-aware coordinate handling
"""

from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Iterator, Optional, Set
import re
import logging

try:
    import pysam
    PYSAM_AVAILABLE = True
except ImportError:
    PYSAM_AVAILABLE = False

logger = logging.getLogger(__name__)


@dataclass
class TranscriptFeature:
    """
    Represents a transcript feature parsed from GTF/GFF.

    Attributes:
        transcript_id: Unique transcript identifier (e.g., ENST00000000001)
        gene_id: Gene identifier (e.g., ENSG00000000001)
        chrom: Chromosome name
        strand: Strand (+, -, or .)
        exons: List of exon (start, end) coordinates in genomic space
        attrs: Dictionary of all attributes from the GTF
        start: Minimum start coordinate (transcript start)
        end: Maximum end coordinate (transcript end)
    """
    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: List[Tuple[int, int]] = field(default_factory=list)
    attrs: Dict[str, str] = field(default_factory=dict)
    start: Optional[int] = None
    end: Optional[int] = None

    def __post_init__(self):
        """Calculate transcript bounds after initialization."""
        if self.exons:
            self.start = min(e[0] for e in self.exons)
            self.end = max(e[1] for e in self.exons)

    @property
    def length(self) -> int:
        """Return total length of transcript (sum of exon lengths)."""
        return sum(end - start for start, end in self.exons)

    @property
    def introns(self) -> List[Tuple[int, int]]:
        """Return list of intron coordinates."""
        if len(self.exons) < 2:
            return []

        introns = []
        for i in range(len(self.exons) - 1):
            _, prev_end = self.exons[i]
            next_start, _ = self.exons[i + 1]
            introns.append((prev_end, next_start))
        return introns

    def get_exon_number(self, genomic_pos: int) -> Optional[int]:
        """
        Get exon number (1-indexed) containing the genomic position.

        Args:
            genomic_pos: Genomic position (0-based)

        Returns:
            Exon number if position is within a transcr exon, None otherwise
        """
        for i, (start, end) in enumerate(self.exons, 1):
            if start <= genomic_pos < end:
                return i
        return None


def parse_gtf_attributes(attribute_string: str) -> Dict[str, str]:
    """
    Parse GTF/GFF attribute string.

    GTF format: "key1 \"value1\"; key2 \"value2\";"
    GFF3 format: "key1=value1;key2=value2"

    Args:
        attribute_string: Attribute field from GTF/GFF

    Returns:
        Dictionary mapping keys to values
    """
    if not attribute_string or attribute_string == '.':
        return {}

    attrs = {}

    # Try GFF3 format first (key=value)
    if '=' in attribute_string or ';' in attribute_string:
        for field in attribute_string.rstrip(';').split(';'):
            field = field.strip()
            if '=' in field:
                key, value = field.split('=', 1)
                attrs[key] = value
            elif ' ' in field:
                # GTF format
                parts = field.split(' ', 1)
                if len(parts) == 2:
                    key = parts[0]
                    value = parts[1].strip('"')
                    attrs[key] = value

    if not attrs:
        # Legacy GTF parsing with regex
        gtf_pattern = r'(\w+)\s+"([^"]+)"'
        matches = re.finditer(gtf_pattern, attribute_string)
        for match in matches:
            key, value = match.groups()
            attrs[key] = value

    return attrs


class GtfParser:
    """
    Parser for GTF/GFF annotation files.

    Supports both tabix-indexed and non-indexed GTF files.
    """

    def __init__(self, filename: str):
        """
        Initialize GTF parser.

        Args:
            filename: Path to GTF or GFF file
        """
        self.filename = filename
        self.transcripts: Dict[str, TranscriptFeature] = {}

        # Try to open as tabix indexed
        self.use_tabix = False
        if PYSAM_AVAILABLE:
            try:
                self.tabix_file = pysam.TabixFile(filename)
                self.use_tabix = True
            except Exception as e:
                logger.info(f"Could not open {filename} as tabix: {e}. Using sequential access.")
                self.use_tabix = False

    def parse(self) -> Iterator[TranscriptFeature]:
        """
        Parse entire GTF file and yield transcript features.

        This method builds transcript models by aggregating exons.

        Yields:
            TranscriptFeature objects with exons grouped
        """
        temp_transcripts: Dict[str, TranscriptFeature] = {}

        with open(self.filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                try:
                    fields = line.split('\t')
                    if len(fields) < 9:
                        logger.warning(f"Skipping malformed line {line_num}: {line[:50]}")
                        continue

                    chrom, source, feature, start, end, score, strand, frame, attrs_str = fields[:9]

                    # Skip if not an exon
                    if feature != 'exon':
                        continue

                    start, end = int(start), int(end)

                    # Parse attributes
                    attrs = parse_gtf_attributes(attrs_str)

                    # Get required IDs
                    transcript_id = attrs.get('transcript_id', '')
                    if not transcript_id:
                        continue  # Skip transcripts without ID

                    gene_id = attrs.get('gene_id', transcript_id)

                    # Initialize transcript if needed
                    if transcript_id not in temp_transcripts:
                        temp_transcripts[transcript_id] = TranscriptFeature(
                            transcript_id=transcript_id,
                            gene_id=gene_id,
                            chrom=chrom,
                            strand=strand,
                            attrs=attrs
                        )

                    # Add exon
                    temp_transcripts[transcript_id].exons.append((start, end))

                except Exception as e:
                    logger.warning(f"Error parsing line {line_num}: {e}")
                    continue

        # Sort exons by coordinate and yield transcripts
        for transcript in temp_transcripts.values():
            # Sort exons based on strand
            if transcript.strand == '-':
                transcript.exons.sort(reverse=True)
            else:
                transcript.exons.sort()

            # Check for overlapping or malformed exons
            clean_exons = self._validate_exons(transcript.exons)
            if clean_exons:
                transcript.exons = clean_exons
                yield transcript

    def _validate_exons(self, exons: List[Tuple[int, int]]) -> Optional[List[Tuple[int, int]]]:
        """
        Validate and clean exon list.

        Args:
            exons: List of (start, end) exon coordinates

        Returns:
            Cleaned exon list, or None if invalid
        """
        if not exons:
            return None

        clean_exons = []
        for i, (start, end) in enumerate(exons):
            if start >= end:
                logger.warning(f"Invalid exon {i}: start={start} >= end={end}")
                continue
            clean_exons.append((start, end))

        # Check for overlapping exons
        if len(clean_exons) > 1:
            for i in range(len(clean_exons) - 1):
                _, prev_end = clean_exons[i]
                next_start, _ = clean_exons[i + 1]
                if prev_end > next_start:
                    logger.warning(f"Overlapping exons detected: {clean_exons}")
                    return None

        return clean_exons if clean_exons else None

    def get_transcripts_in_region(self, chrom: str, start: int, end: int) -> List[TranscriptFeature]:
        """
        Get all transcripts overlapping a genomic region.

        Args:
            chrom: Chromosome name
            start: Start coordinate (0-based)
            end: End coordinate

        Returns:
            List of TranscriptFeature objects overlapping the region
        """
        transcripts = []

        if self.use_tabix and PYSAM_AVAILABLE:
            # Use tabix for fast region queries
            try:
                # Tabix is 0-based, inclusive-start, exclusive-end
                for line in self.tabix_file.fetch(chrom, start, end):
                    fields = line.split('\t')
                    if len(fields) < 9:
                        continue

                    if fields[2] != 'exon':
                        continue

                    transcript_start, transcript_end = int(fields[3]), int(fields[4])

                    # Check overlap
                    if (transcript_start <= end and transcript_end >= start):
                        # Parse this transcript
                        attrs = parse_gtf_attributes(fields[8])
                        transcript_id = attrs.get('transcript_id', '')
                        gene_id = attrs.get('gene_id', transcript_id)

                        if not transcript_id:
                            continue

                        # Build transcript object
                        transcript = TranscriptFeature(
                            transcript_id=transcript_id,
                            gene_id=gene_id,
                            chrom=chrom,
                            strand=fields[6],
                            attrs=attrs
                        )
                        transcript.exons.append((transcript_start, transcript_end))

                        transcripts.append(transcript)

                return transcripts

            except Exception as e:
                logger.warning(f"Tabix fetch failed: {e}")

        # Fallback: sequential parsing
        for transcript in self.parse():
            if transcript.chrom != chrom:
                continue

            transcript_start = min(e[0] for e in transcript.exons)
            transcript_end = max(e[1] for e in transcript.exons)

            if transcript_start <= end and transcript_end >= start:
                transcripts.append(transcript)

        return transcripts

    def get_all_genes(self) -> Dict[str, List[TranscriptFeature]]:
        """
        Get all genes and their transcripts from the GTF file.

        Returns:
            Dictionary mapping gene_id to list of TranscriptFeature objects
        """
        genes: Dict[str, List[TranscriptFeature]] = {}

        for transcript in self.parse():
            if transcript.gene_id not in genes:
                genes[transcript.gene_id] = []
            genes[transcript.gene_id].append(transcript)

        return genes

    def get_transcript_by_id(self, transcript_id: str) -> Optional[TranscriptFeature]:
        """
        Get a specific transcript by ID.

        Args:
            transcript_id: Transcript identifier

        Returns:
            TranscriptFeature object or None if not found
        """
        for transcript in self.parse():
            if transcript.transcript_id == transcript_id:
                return transcript
        return None


def create_gtf_index(gtf_file: str) -> Optional[str]:
    """
    Create tabix index for GTF file if it doesn't exist.

    Args:
        gtf_file: Path to GTF file

    Returns:
        Path to index file if successful, None otherwise
    """
    if not PYSAM_AVAILABLE:
        logger.warning("pysam not available for GTF indexing")
        return None

    try:
        # Sort GTF if needed
        import subprocess
        sorted_gtf = gtf_file.replace('.gtf', '.sorted.gtf')

        if not Path(sorted_gtf).exists():
            logger.info(f"Sorting GTF file: {gtf_file}")
            subprocess.run([
                'sort', '-k1,1', '-k4,4n', gtf_file
            ], stdout=open(sorted_gtf, 'w'), check=True)

        # Index with tabix
        index_file = sorted_gtf + '.tbi'
        pysam.tabix_index(sorted_gtf, preset='gff')
        logger.info(f"Created GTF index: {index_file}")

        return index_file

    except Exception as e:
        logger.error(f"Failed to create GTF index: {e}")
        return None


# Example usage and test
if __name__ == "__main__":
    # Simple test
    test_attrs = 'gene_id "ENSG00000141510"; transcript_id "ENST00000269305";'
    parsed = parse_gtf_attributes(test_attrs)
    print("Test attribute parsing:", parsed)
