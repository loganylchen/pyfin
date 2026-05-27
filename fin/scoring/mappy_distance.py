"""M1: read×candidate distance via mappy (minimap2) alignment scores.

Provides ``compute_mappy_distance`` which aligns each read sequence against
every candidate transcript using ``mappy.Aligner`` (preset ``map-ont``) and
returns a non-negative distance matrix suitable for ``em_with_coherence``:

    d[i, j] = max_score_across_candidates - score(read_i, cand_j)

where ``score`` is the best primary-alignment score for the read against
candidate j (``hit.score`` when available, falling back to ``hit.mlen``).
Reads with no hit against a candidate get the row max (worst), so they are
not penalised relative to candidates they failed against.

The transform to a *distance* keeps the same monotonic ordering as the raw
mappy AS (higher AS = closer), guarantees ``d >= 0`` (required by EM), and
is invariant to the absolute AS scale per read.

Public function:
    compute_mappy_distance(read_seqs, candidates, read_ids) -> np.ndarray
"""
from __future__ import annotations

import logging
from typing import Dict, List

import numpy as np

from fin.candidates.dataclasses import TranscriptCandidate
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit

logger = logging.getLogger(__name__)


def compute_mappy_distance(
    read_seqs: Dict[str, str],
    candidates: List[TranscriptCandidate],
    read_ids: List[str],
) -> np.ndarray:
    """Return (n_reads, n_cands) distance matrix from mappy alignments.

    Args:
        read_seqs: Mapping read_id -> sequence (already strand-corrected).
        candidates: Ordered list of TranscriptCandidate (axis 1).
        read_ids: Ordered read IDs (axis 0). Reads not in ``read_seqs`` get
            a row of zeros (i.e. uninformative; will be reweighted in EM).

    Returns:
        np.ndarray float32 shape (len(read_ids), len(candidates)). All
        entries are non-negative.
    """
    import mappy

    n_r = len(read_ids)
    n_c = len(candidates)
    if n_r == 0 or n_c == 0:
        return np.zeros((n_r, n_c), dtype=np.float32)

    # Build per-candidate aligners once.
    aligners: List[mappy.Aligner | None] = []
    for c in candidates:
        if not c.sequence:
            aligners.append(None)
            continue
        a = mappy.Aligner(seq=c.sequence, preset=get_m1_preset())
        aligners.append(a if a else None)

    # Raw AS matrix; missing entries get -inf so we can compute per-row max.
    raw = np.full((n_r, n_c), -np.inf, dtype=np.float64)
    for i, rid in enumerate(read_ids):
        seq = read_seqs.get(rid, "")
        if not seq:
            continue
        for j, aln in enumerate(aligners):
            if aln is None:
                continue
            best = None
            for hit in aln.map(seq):
                s = score_hit(hit)
                if s is None:
                    continue
                if best is None or s > best:
                    best = s
            if best is not None:
                raw[i, j] = float(best)

    # Per-row distance = row_max - score; missing entries (-inf) get the
    # full row range so they are no better than any actual hit.
    dist = np.zeros((n_r, n_c), dtype=np.float32)
    for i in range(n_r):
        row = raw[i]
        finite_mask = np.isfinite(row)
        if not finite_mask.any():
            # No hits at all → leave row at zeros (uninformative).
            continue
        row_max = float(row[finite_mask].max())
        # Replace -inf with the row's worst (most negative) finite score
        # minus 1, so missing aligns convert to the largest distance.
        row_min = float(row[finite_mask].min())
        clean = np.where(finite_mask, row, row_min - 1.0)
        dist[i] = (row_max - clean).astype(np.float32)
    return dist
