"""Scoring modules: eventalign parsing, signal DTW, composite scoring (krill-only)."""

from .eventalign_parser import (
    ReadCandidateScore,
    build_distance_matrix,
    parse_eventalign_tsv,
)
from .signal_dtw import compute_read_to_read_dtw, extract_signal_segments
from .composite import (
    CompositeScore,
    compute_coherence_score,
    compute_discrimination_score,
    compute_combined_score,
    score_candidates_composite,
)

__all__ = [
    "ReadCandidateScore",
    "parse_eventalign_tsv",
    "build_distance_matrix",
    "extract_signal_segments",
    "compute_read_to_read_dtw",
    "CompositeScore",
    "compute_coherence_score",
    "compute_discrimination_score",
    "compute_combined_score",
    "score_candidates_composite",
]
