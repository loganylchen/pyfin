"""R1 bypass: mappy alignment-only assignment (no signal, no EM).

Two variants are provided:
 - ``mappy_argmax_assignment``: each read → single best-AS candidate (hard).
 - ``mappy_multimap_responsibilities``: each read distributes AS-weighted
   responsibility across all hit candidates (soft; salmon/NanoCount-style
   alignment-only baseline). This is the default R1 path because hard argmax
   discards alignment uncertainty and unfairly handicaps the baseline.
"""

from __future__ import annotations

import logging
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np

from fin.candidates.dataclasses import TranscriptCandidate
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit

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
        aln = mappy.Aligner(seq=cand.sequence, preset=get_m1_preset())
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
                # Reconstructed map-ont AS; None when a single indel exceeds the
                # cap (structural exon difference → treated as not mapped).
                hit_score = score_hit(hit)
                if hit_score is None:
                    continue
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


def mappy_multimap_responsibilities(
    reads: Iterable[Tuple[str, str]],
    candidates: List[TranscriptCandidate],
) -> Tuple[np.ndarray, List[str]]:
    """Multi-mapping AS-weighted soft assignment (R1 default baseline).

    For each read, aligns to every candidate, takes the max alignment score per
    candidate, and distributes responsibility proportional to AS across all
    candidates the read hits. This is the standard "alignment-only" baseline
    used by salmon/NanoCount: r_ic = AS_ic / Σ_c' AS_ic'.

    Args:
        reads: Iterable of (read_id, read_sequence) tuples.
        candidates: List of TranscriptCandidate. Output columns align to this
            list's order; candidates with empty sequence get all-zero columns.

    Returns:
        (R, kept_read_ids) where R[i, j] is the AS-weighted responsibility of
        read i (kept_read_ids[i]) for candidate j (candidates[j].candidate_id).
        Rows sum to 1.0. Reads with no hit anywhere are dropped.
    """
    import mappy

    aligners: List[Optional["mappy.Aligner"]] = []
    for cand in candidates:
        if not cand.sequence:
            aligners.append(None)
            continue
        aln = mappy.Aligner(seq=cand.sequence, preset=get_m1_preset())
        if not aln:
            logger.warning(
                "mappy_multimap_responsibilities: failed to build aligner for %s",
                cand.candidate_id,
            )
            aligners.append(None)
            continue
        aligners.append(aln)

    # Optional env-var knobs (R1 tuning sweep):
    #   MAPPY_R1_T (float)       softmax temperature; if unset → linear normalize
    #   MAPPY_R1_MIN_AS (float)  drop hits with AS < this value before weighting
    import os as _os
    _T = _os.environ.get("MAPPY_R1_T")
    _T = float(_T) if _T else None
    _MIN_AS = float(_os.environ.get("MAPPY_R1_MIN_AS", "0") or "0")

    n_cands = len(candidates)
    rows: List[np.ndarray] = []
    kept_read_ids: List[str] = []
    dropped = 0

    for read_id, seq in reads:
        if not seq:
            dropped += 1
            continue
        row = np.zeros(n_cands, dtype=np.float32)
        for j, aln in enumerate(aligners):
            if aln is None:
                continue
            best = 0.0
            for hit in aln.map(seq):
                s = score_hit(hit)
                if s is None:
                    continue
                if s > best:
                    best = float(s)
            if best >= _MIN_AS:
                row[j] = best
        if _T is not None:
            # softmax(AS / T) over candidates with AS > 0
            mask = row > 0
            if not mask.any():
                dropped += 1
                continue
            logits = row.copy()
            logits[~mask] = -np.inf
            logits = logits / _T
            logits -= logits[mask].max()
            ex = np.zeros_like(row)
            ex[mask] = np.exp(logits[mask])
            row = ex / ex.sum()
        else:
            total = float(row.sum())
            if total <= 0.0:
                dropped += 1
                continue
            row /= total
        rows.append(row)
        kept_read_ids.append(read_id)

    if dropped:
        logger.info(
            "mappy_multimap_responsibilities: dropped %d reads with no hit",
            dropped,
        )
    if not rows:
        return np.zeros((0, n_cands), dtype=np.float32), []
    return np.stack(rows, axis=0), kept_read_ids


def per_tx_counts_from_responsibilities(
    R: np.ndarray,
    candidates: List[TranscriptCandidate],
) -> Dict[str, float]:
    """Sum responsibilities per candidate.

    Args:
        R: (n_reads, n_candidates) responsibility matrix; columns aligned to
            ``candidates`` order.
        candidates: Candidate list (defines output keys).

    Returns:
        Dict candidate_id -> fractional count (float).
    """
    if R.size == 0:
        return {c.candidate_id: 0.0 for c in candidates}
    col_sums = R.sum(axis=0)
    return {c.candidate_id: float(col_sums[j]) for j, c in enumerate(candidates)}
