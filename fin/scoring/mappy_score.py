"""Reconstruct an alignment score (AS) from a mappy hit + reject large indels.

minimap2 map-ont defaults: -A2 -B4 -O4,24 -E2,1 (two-piece affine gap).
AS = A*matches - B*mismatches - sum_gaps min(O1+k*E1, O2+k*E2)

A single indel longer than ``M1_MAX_INDEL_BP`` (default 50) means the read is
treated as not mapped to that candidate (``score_hit`` returns ``None``), so a
structural exon-difference DEL/INS does not masquerade as a good alignment.
"""
from __future__ import annotations

import os
from typing import Optional

_A, _B = 2, 4
_O1, _E1 = 4, 2
_O2, _E2 = 24, 1


def _max_indel_bp() -> int:
    try:
        return int(os.environ.get("M1_MAX_INDEL_BP", "50"))
    except ValueError:
        return 50


def _gap_cost(k: int) -> int:
    return min(_O1 + k * _E1, _O2 + k * _E2)


def score_hit(hit, max_indel_bp: Optional[int] = None) -> Optional[float]:
    """Reconstructed map-ont AS for a mappy hit, or None if a single indel > cap.

    Args:
        hit: a mappy alignment hit (must expose ``.cigar``, ``.NM``, ``.mlen``).
        max_indel_bp: per-indel length cap; if any single I/D op exceeds it the
            hit is rejected (None). Defaults to ``M1_MAX_INDEL_BP`` env (50).

    Returns:
        float AS, or None if rejected / cigar unavailable.
    """
    if max_indel_bp is None:
        max_indel_bp = _max_indel_bp()
    cigar = getattr(hit, "cigar", None)
    if not cigar:
        return None
    indel_bases = 0
    gap_cost = 0
    for ln, op in cigar:
        if op == 1 or op == 2:  # insertion / deletion
            if ln > max_indel_bp:
                return None
            indel_bases += ln
            gap_cost += _gap_cost(ln)
    mismatches = max(0, hit.NM - indel_bases)
    return float(_A * hit.mlen - _B * mismatches - gap_cost)
