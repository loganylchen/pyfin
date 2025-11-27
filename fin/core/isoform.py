"""
Isoform data structures for nanopore RNA isoform detection.

Defines data models for representing transcript sequences,
isform evidence, and validated isoforms with support metrics.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Any, Tuple, Optional, Set
import numpy as np


@dataclass
class IsoformSequence:
    """>
    Represents an isoform sequence (whether from annotation or derived from reads).

    This data structure holds both the genomic coordinates and the actual
    nucleotide sequence of a transcript isoform.

    Attributes:
        isoform_id: Unique identifier
                     - For annotated: transcript_id from GTF
                     - For read-derived: read_id or generated ID
        gene_id: Associated gene identifier
        chrom: Chromosome name
        strand: Strand orientation (+, -, or .)
        exons: List of exon (start, end) coordinates in genomic space
        sequence: Actual nucleotide sequence (from reference or read)
        source: Origin of the isoform ("annotation" or "read")
        metadata: Additional information (e.g., from GTF attributes)
    """
    isoform_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: List[Tuple[int, int]] = field(default_factory=list)
    sequence: str = ""
    source: str = "annotation"  # "annotation" or "read"
    metadata: Dict[str, Any] = field(default_factory=dict)

    @property
    def length(self) -> int:
        """Return total transcript length (sum of exon lengths)."""
        return sum(end - start for start, end in self.exons)

    @property
    def start(self) -> Optional[int]:
        """Return minimum start coordinate."""
        return min(e[0] for e in self.exons) if self.exons else None

    @property
    def end(self) -> Optional[int]:
        """Return maximum end coordinate."""
        return max(e[1] for e in self.exons) if self.exons else None

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

    @property
    def num_exons(self) -> int:
        """Return number of exons."""
        return len(self.exons)

    @property
    def num_introns(self) -> int:
        """Return number of introns."""
        return len(self.introns)

    def get_exon_number(self, genomic_pos: int) -> Optional[int]:
        """
        Get exon number (1-indexed) containing the genomic position.

        Args:
            genomic_pos: Genomic position (0-based)

        Returns:
            Exon number if position is within transcript exon, None otherwise
        """
        for i, (start, end) in enumerate(self.exons, 1):
            if start <= genomic_pos < end:
                return i
        return None

    def get_exon_bounds(self, exon_number: int) -> Optional[Tuple[int, int]]:
        """
        Get genomic coordinates of a specific exon.

        Args:
            exon_number: 1-indexed exon number

        Returns:
            (start, end) tuple or None if exon doesn't exist
        """
        if 1 <= exon_number <= len(self.exons):
            return self.exons[exon_number - 1]
        return None

    def is_pos_in_exon(self, genomic_pos: int) -> bool:
        """Check if a genomic position is within any exon."""
        return self.get_exon_number(genomic_pos) is not None

    def is_pos_in_intron(self, genomic_pos: int) -> bool:
        """Check if a genomic position is within any intron."""
        for start, end in self.introns:
            if start <= genomic_pos < end:
                return True
        return False

    def __repr__(self) -> str:
        return (f"IsoformSequence({self.isoform_id}, gene={self.gene_id}, "
                f"{self.num_exons} exons, length={self.length})")


@dataclass
class IsoformEvidence:
    """
    Evidence supporting an isoform based on eventalign results.

    This stores the results of aligning a read's signal to an isoform
    sequence, including completeness scores and alignment metrics.

    Attributes:
        isoform_id: ID of the isoform
        read_id: ID of the read
        completeness_score: Proportion of isoform covered by aligned events (0-1)
        num_aligned_events: Number of events that aligned to isoform
        alignment_score: Mean alignment score from eventalign
        signal_length: Length of the read signal in samples
        aligned_positions: Set of isoform positions with alignments
    """
    isoform_id: str
    read_id: str
    completeness_score: float = 0.0
    num_aligned_events: int = 0
    alignment_score: float = 0.0
    signal_length: int = 0
    aligned_positions: Optional[Set[int]] = None

    @property
    def is_high_quality(self) -> bool:
        """Check if this is high-quality evidence."""
        return (self.completeness_score >= 0.8 and
                self.num_aligned_events >= 50 and
                self.alignment_score >= 0.8)

    def __repr__(self) -> str:
        return (f"IsoformEvidence({self.read_id} -> {self.isoform_id}, "
                f"completeness={self.completeness_score:.3f})")


@dataclass
class ValidatedIsoform:
    """
    Final validated isoform with support metrics.

    This represents an isoform that has been validated through the
    pipeline by combining eventalign completeness and DTW clustering.

    Attributes:
        isoform_id: Unique isoform identifier
        gene_id: Associated gene
        chrom: Chromosome
        strand: Strand orientation
        exons: Exon coordinates
        sequence: Nucleotide sequence
        read_support: Number of supporting reads
        avg_completeness: Mean completeness across supporting reads
        confidence_score: Overall confidence (combination of support and quality)
        supporting_reads: List of read IDs that support this isoform
        is_novel: True if not present in reference annotation
    """
    isoform_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: List[Tuple[int, int]]
    sequence: str
    read_support: int = 0
    avg_completeness: float = 0.0
    confidence_score: float = 0.0
    supporting_reads: List[str] = field(default_factory=list)
    is_novel: bool = False
    metadata: Dict[str, Any] = field(default_factory=dict)
    estimated_count: Optional[float] = None  # EM algorithm estimated read count
    tpm: Optional[float] = None  # Transcripts Per Million (expression level)
    fpkm: Optional[float] = None  # Fragments Per Kilobase Million (expression level)
    isoform_fraction: Optional[float] = None  # Fraction of gene's expression

    def __post_init__(self):
        """Initialize optional quantification fields if not set."""
        if self.estimated_count is None:
            self.estimated_count = float(self.read_support)
        if self.isoform_fraction is None and self.read_support > 0:
            self.isoform_fraction = 1.0
        if self.tpm is None:
            self.tpm = 0.0
        if self.fpkm is None:
            self.fpkm = 0.0

    @property
    def length(self) -> int:
        """Return total length."""
        return sum(end - start for start, end in self.exons)

    @property
    def start(self) -> Optional[int]:
        """Return minimum start coordinate."""
        return min(e[0] for e in self.exons) if self.exons else None

    @property
    def end(self) -> Optional[int]:
        """Return maximum end coordinate."""
        return max(e[1] for e in self.exons) if self.exons else None

    @property
    def num_exons(self) -> int:
        """Return number of exons."""
        return len(self.exons)

    @property
    def has_strong_support(self) -> bool:
        """Check if isoform has strong read support."""
        return (self.read_support >= 10 and
                self.avg_completeness >= 0.85 and
                self.confidence_score >= 8.0)

    @property
    def is_well_supported(self) -> bool:
        """Check if isoform meets minimum support criteria."""
        return self.read_support >= 5 and self.avg_completeness >= 0.80

    def __repr__(self) -> str:
        support_type = "novel" if self.is_novel else "annotated"
        return (f"ValidatedIsoform({self.isoform_id}, type={support_type}, "
                f"{self.read_support} reads, "
                f"completeness={self.avg_completeness:.3f})")

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "isoform_id": self.isoform_id,
            "gene_id": self.gene_id,
            "chrom": self.chrom,
            "strand": self.strand,
            "exons": self.exons,
            "length": self.length,
            "read_support": self.read_support,
            "avg_completeness": self.avg_completeness,
            "confidence_score": self.confidence_score,
            "is_novel": self.is_novel,
            "metadata": self.metadata
        }


# Convenience functions
def group_isoforms_by_gene(isoforms: List[IsoformSequence]) -> Dict[str, List[IsoformSequence]]:
    """Group isoforms by gene_id."""
    genes: Dict[str, List[IsoformSequence]] = {}
    for isoform in isoforms:
        if isoform.gene_id not in genes:
            genes[isoform.gene_id] = []
        genes[isoform.gene_id].append(isoform)
    return genes


def filter_isoforms_by_length(
    isoforms: List[IsoformSequence],
    min_length: int = 200,
    max_length: Optional[int] = None
) -> List[IsoformSequence]:
    """Filter isoforms by length."""
    filtered = []
    for isoform in isoforms:
        if isoform.length < min_length:
            continue
        if max_length and isoform.length > max_length:
            continue
        filtered.append(isoform)
    return filtered


def find_novel_isoforms(
    validated_isoforms: List[ValidatedIsoform]
) -> List[ValidatedIsoform]:
    """Extract only novel isoforms from validated list."""
    return [vi for vi in validated_isoforms if vi.is_novel]


def sort_isoforms_by_support(
    validated_isoforms: List[ValidatedIsoform],
    descending: bool = True
) -> List[ValidatedIsoform]:
    """Sort isoforms by read support and completeness."""
    return sorted(
        validated_isoforms,
        key=lambda vi: (vi.read_support, vi.avg_completeness),
        reverse=descending
    )
