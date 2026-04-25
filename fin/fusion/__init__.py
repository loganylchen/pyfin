"""Fusion transcript detection and candidate building.

Exports
-------
Breakpoint
    Dataclass representing a split-read fusion breakpoint.
parse_sa_tags
    Parse SA tags from a BAM file into raw Breakpoints.
cluster_breakpoints
    Cluster and filter raw Breakpoints.
build_fusion_candidates
    Build TranscriptCandidate objects from clustered Breakpoints (requires
    ``fin.fusion.candidates`` which is implemented separately; gracefully
    resolves to *None* if that module is not yet available).
"""

from __future__ import annotations

from fin.fusion.breakpoints import Breakpoint, cluster_breakpoints, parse_sa_tags

try:
    from fin.fusion.candidates import build_fusion_candidates
except ImportError:
    build_fusion_candidates = None  # type: ignore[assignment]

__all__ = [
    "Breakpoint",
    "parse_sa_tags",
    "cluster_breakpoints",
    "build_fusion_candidates",
]
