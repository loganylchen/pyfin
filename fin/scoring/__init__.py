"""Scoring modules: external tools, eventalign parsing, signal DTW."""

from .eventalign_parser import (
    ReadCandidateScore,
    build_distance_matrix,
    parse_eventalign_tsv,
)
from .external_tools import ExternalToolPaths, ExternalToolRunner
from .signal_dtw import compute_read_to_read_dtw, extract_signal_segments

__all__ = [
    "ExternalToolPaths",
    "ExternalToolRunner",
    "ReadCandidateScore",
    "parse_eventalign_tsv",
    "build_distance_matrix",
    "extract_signal_segments",
    "compute_read_to_read_dtw",
]
