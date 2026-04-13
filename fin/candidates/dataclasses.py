"""Data classes for transcript candidate discovery."""

from __future__ import annotations

from dataclasses import dataclass, field
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
    """A candidate transcript (from GTF annotation or novel discovery)."""

    candidate_id: str
    intron_chain: IntronChain
    three_prime_pos: int
    sequence: str
    source: str  # "gtf" or "novel"
    supporting_read_ids: Set[str]
    chrom: str
    strand: str
    start: int
    end: int

    @property
    def length(self) -> int:
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
