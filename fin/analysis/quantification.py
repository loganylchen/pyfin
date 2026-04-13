"""Probability-weighted transcript quantification from EM assignments."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np

from fin.candidates.dataclasses import TranscriptCandidate

logger = logging.getLogger(__name__)


@dataclass
class QuantResult:
    """Quantification result for a single transcript candidate."""

    candidate_id: str
    abundance: float  # probability-weighted read count
    confidence: float  # mean assignment probability for assigned reads
    num_assigned_reads: int  # number of hard-assigned reads
    source: str  # "gtf" or "novel"


def quantify_transcripts(
    R: np.ndarray,
    hard_assignments: np.ndarray,
    candidates: List[TranscriptCandidate],
    read_ids: List[str],
) -> List[QuantResult]:
    """Compute probability-weighted transcript abundance.

    Args:
        R: Soft assignment matrix of shape (n_reads, n_candidates).
            R[i, j] = P(read_i belongs to candidate_j).
        hard_assignments: Array of shape (n_reads,) with hard cluster indices.
        candidates: List of TranscriptCandidate objects (column order matches R).
        read_ids: List of read IDs (row order matches R).

    Returns:
        List of QuantResult, one per candidate.
    """
    n_reads, n_cands = R.shape

    results = []
    for j, cand in enumerate(candidates):
        # Probability-weighted count
        abundance = float(np.sum(R[:, j]))

        # Hard-assigned reads for this candidate
        assigned_mask = hard_assignments == j
        num_assigned = int(np.sum(assigned_mask))

        # Mean confidence of assigned reads
        if num_assigned > 0:
            confidence = float(np.mean(R[assigned_mask, j]))
        else:
            confidence = 0.0

        results.append(
            QuantResult(
                candidate_id=cand.candidate_id,
                abundance=abundance,
                confidence=confidence,
                num_assigned_reads=num_assigned,
                source=cand.source,
            )
        )

    return results


def aggregate_across_intervals(
    interval_results: List[List[QuantResult]],
) -> Dict[str, QuantResult]:
    """Aggregate quantification results across multiple intervals.

    For the same candidate_id appearing in multiple intervals, sums
    abundance and num_assigned_reads, and averages confidence.

    Args:
        interval_results: List of per-interval QuantResult lists.

    Returns:
        Dict mapping candidate_id -> aggregated QuantResult.
    """
    agg: Dict[str, dict] = {}

    for results in interval_results:
        for qr in results:
            if qr.candidate_id not in agg:
                agg[qr.candidate_id] = {
                    "abundance": 0.0,
                    "confidence_sum": 0.0,
                    "confidence_count": 0,
                    "num_assigned_reads": 0,
                    "source": qr.source,
                }
            a = agg[qr.candidate_id]
            a["abundance"] += qr.abundance
            a["num_assigned_reads"] += qr.num_assigned_reads
            if qr.num_assigned_reads > 0:
                a["confidence_sum"] += qr.confidence * qr.num_assigned_reads
                a["confidence_count"] += qr.num_assigned_reads

    result = {}
    for cid, a in agg.items():
        confidence = (
            a["confidence_sum"] / a["confidence_count"]
            if a["confidence_count"] > 0
            else 0.0
        )
        result[cid] = QuantResult(
            candidate_id=cid,
            abundance=a["abundance"],
            confidence=confidence,
            num_assigned_reads=a["num_assigned_reads"],
            source=a["source"],
        )

    return result
