"""Scoring modules: eventalign parsing, signal DTW (krill-only)."""

from .eventalign_parser import (
    ReadCandidateScore,
    build_distance_matrix,
    parse_eventalign_tsv,
)
from .signal_dtw import compute_read_to_read_dtw, extract_signal_segments

__all__ = [
    "ReadCandidateScore",
    "parse_eventalign_tsv",
    "build_distance_matrix",
    "extract_signal_segments",
    "compute_read_to_read_dtw",
]
