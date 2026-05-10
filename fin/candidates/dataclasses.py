"""Data classes for transcript candidate discovery."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Set, Tuple

from fin.io.interval_manager import GenomicInterval


@dataclass(frozen=True)
class IntronChain:
    """Hashable representation of a transcript's intron structure.

    Each intron is a (start, end) tuple in genomic coordinates.
    """

    introns: Tuple[Tuple[int, int], ...]

    @property
    def num_introns(self) -> int:
        return len(self.introns)

    @property
    def is_single_exon(self) -> bool:
        return len(self.introns) == 0


@dataclass
class TranscriptCandidate:
    """A candidate transcript (from GTF annotation, novel discovery, or fusion).

    Required fields (positional):
        candidate_id, intron_chain, three_prime_pos, sequence, source,
        supporting_read_ids, chrom, strand, start, end.

    Optional fusion fields (populated only when source == "fusion"):
        fusion_junction: Tuple of (left_genomic_pos, right_genomic_pos) marking
            where the two parental transcripts are joined.
        breakpoint_left: (chrom, pos, strand) for the left fusion partner.
        breakpoint_right: (chrom, pos, strand) for the right fusion partner.
    """

    candidate_id: str
    intron_chain: IntronChain
    three_prime_pos: int
    sequence: str
    source: str  # "gtf", "novel", or "fusion"
    supporting_read_ids: Set[str]
    chrom: str
    strand: str
    start: int
    end: int

    # Fusion-specific (populated only when source == "fusion")
    fusion_junction: Optional[Tuple[int, int]] = None
    breakpoint_left: Optional[Tuple[str, int, str]] = None   # (chrom, pos, strand)
    breakpoint_right: Optional[Tuple[str, int, str]] = None

    @property
    def length(self) -> int:
        # Fusion candidates set start=posA on chromA and end=posB on chromB.
        # Subtracting them is meaningless (and can be negative) so we fall
        # back to the assembled sequence length, which is the only thing
        # downstream code actually uses for fusion candidates.
        if self.source == "fusion":
            return len(self.sequence) if self.sequence else 0
        return self.end - self.start

    @property
    def num_supporting_reads(self) -> int:
        return len(self.supporting_read_ids)


@dataclass
class CandidateSet:
    """Collection of transcript candidates for a genomic interval."""

    interval: GenomicInterval
    candidates: List[TranscriptCandidate]
    read_ids: Set[str]

    @property
    def num_candidates(self) -> int:
        return len(self.candidates)

    def get_candidate(self, candidate_id: str) -> Optional[TranscriptCandidate]:
        for c in self.candidates:
            if c.candidate_id == candidate_id:
                return c
        return None

    def candidate_ids(self) -> List[str]:
        return [c.candidate_id for c in self.candidates]
