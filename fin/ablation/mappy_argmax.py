"""R1 bypass: mappy alignment-score argmax assignment (no signal, no EM)."""

from __future__ import annotations

import logging
from typing import Dict, Iterable, List, Optional, Tuple

from fin.candidates.dataclasses import TranscriptCandidate

logger = logging.getLogger(__name__)


def mappy_argmax_assignment(
    reads: Iterable[Tuple[str, str]],
    candidates: List[TranscriptCandidate],
) -> Dict[str, str]:
    """Align each read against every candidate and assign by best mappy AS.

    Args:
        reads: Iterable of (read_id, read_sequence) tuples.
        candidates: List of TranscriptCandidate; each must have a non-empty
            ``sequence`` and ``candidate_id``.

    Returns:
        Dict mapping read_id -> best candidate_id (reads with no hit dropped).
    """
    import mappy

    aligners: List[Tuple[str, mappy.Aligner]] = []
    for cand in candidates:
        if not cand.sequence:
            logger.warning(
                "mappy_argmax_assignment: candidate %s has empty sequence, skipping",
                cand.candidate_id,
            )
            continue
        aln = mappy.Aligner(seq=cand.sequence, preset="map-ont")
        if not aln:
            logger.warning(
                "mappy_argmax_assignment: failed to build aligner for %s",
                cand.candidate_id,
            )
            continue
        aligners.append((cand.candidate_id, aln))

    assignment: Dict[str, str] = {}
    if not aligners:
        return assignment

    dropped = 0
    for read_id, seq in reads:
        if not seq:
            dropped += 1
            continue
        best_cid: Optional[str] = None
        best_score: float = float("-inf")
        for cid, aln in aligners:
            for hit in aln.map(seq):
                # mappy hit has integer alignment score `.mlen` and `.NM`; the
                # AS-equivalent is `hit.mlen - hit.NM` weighted by primary,
                # but mappy exposes a direct alignment score via `.score` in
                # newer versions. Fall back to mlen on older versions.
                hit_score = getattr(hit, "score", None)
                if hit_score is None:
                    hit_score = hit.mlen
                if hit_score > best_score:
                    best_score = hit_score
                    best_cid = cid
        if best_cid is not None:
            assignment[read_id] = best_cid
        else:
            dropped += 1

    if dropped:
        logger.info("mappy_argmax_assignment: dropped %d reads with no hit", dropped)
    return assignment


def per_tx_counts_from_argmax(
    assignment: Dict[str, str],
    candidates: List[TranscriptCandidate],
) -> Dict[str, float]:
    """Aggregate per-candidate counts from a read->candidate assignment.

    Args:
        assignment: read_id -> candidate_id.
        candidates: All candidates (used so unassigned ones show count=0).

    Returns:
        Dict candidate_id -> count (float, but always integer-valued for R1).
    """
    counts: Dict[str, float] = {c.candidate_id: 0.0 for c in candidates}
    for cid in assignment.values():
        if cid in counts:
            counts[cid] += 1.0
        else:
            # Unknown candidate id (shouldn't happen if assigners came from
            # the same candidate list); still tally for completeness.
            counts[cid] = counts.get(cid, 0.0) + 1.0
    return counts
