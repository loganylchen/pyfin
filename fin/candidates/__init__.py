"""Transcript candidate discovery modules."""

from .dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from .discovery import discover_candidates
from .intron_chains import (
    extract_intron_chain,
    group_reads_by_three_prime_and_intron_chain,
    gtf_transcript_to_intron_chain,
    pick_representative_read,
)

__all__ = [
    "IntronChain",
    "TranscriptCandidate",
    "CandidateSet",
    "discover_candidates",
    "extract_intron_chain",
    "gtf_transcript_to_intron_chain",
    "group_reads_by_three_prime_and_intron_chain",
    "pick_representative_read",
]
