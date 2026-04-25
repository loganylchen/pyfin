"""Probability-weighted transcript quantification from EM assignments."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

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
    chrom: str = ""
    strand: str = "."
    start: int = 0
    end: int = 0
    exons: Tuple[Tuple[int, int], ...] = ()  # 0-based
    gene_id: str = ""
    coherence_score: float = 0.0
    discrimination_score: float = 0.0
    combined_score: float = 0.0
    breakpoint_left: Optional[Tuple[str, int, str]] = None
    breakpoint_right: Optional[Tuple[str, int, str]] = None


def _exons_from_candidate(c: TranscriptCandidate) -> Tuple[Tuple[int, int], ...]:
    """Derive exon coordinates from a candidate's intron chain and boundaries."""
    introns = c.intron_chain.introns
    if not introns:
        return ((c.start, c.end),)
    exons = [(c.start, introns[0][0])]
    for i in range(len(introns) - 1):
        exons.append((introns[i][1], introns[i + 1][0]))
    exons.append((introns[-1][1], c.end))
    return tuple(exons)


def quantify_transcripts(
    R: np.ndarray,
    hard_assignments: np.ndarray,
    candidates: List[TranscriptCandidate],
    read_ids: List[str],
) -> List[QuantResult]:
    """Compute probability-weighted transcript abundance (vectorized).

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

    # Vectorized abundance: sum of soft assignments per candidate
    abundance = R.sum(axis=0)  # shape (n_cands,)

    # Vectorized hard assignment counts and confidence
    results = []
    for j, cand in enumerate(candidates):
        assigned_mask = hard_assignments == j
        num_assigned = int(assigned_mask.sum())
        confidence = float(R[assigned_mask, j].mean()) if num_assigned > 0 else 0.0

        results.append(
            QuantResult(
                candidate_id=cand.candidate_id,
                abundance=float(abundance[j]),
                confidence=confidence,
                num_assigned_reads=num_assigned,
                source=cand.source,
                chrom=cand.chrom,
                strand=cand.strand,
                start=cand.start,
                end=cand.end,
                exons=_exons_from_candidate(cand),
            )
        )

    return results


def compute_tpm(
    results: Dict[str, QuantResult],
    transcript_lengths: Dict[str, int],
) -> Dict[str, float]:
    """Compute Transcripts Per Million (TPM) from abundance and transcript lengths.

    TPM = (abundance / length_kb) / sum(abundance / length_kb) * 1e6

    Args:
        results: Dict mapping candidate_id -> QuantResult.
        transcript_lengths: Dict mapping candidate_id -> spliced exon length in bp.

    Returns:
        Dict mapping candidate_id -> TPM value.
    """
    # Compute RPK (reads per kilobase) for each transcript
    rpk = {}
    for cid, qr in results.items():
        length = transcript_lengths.get(cid, 0)
        if length > 0 and qr.abundance > 0:
            rpk[cid] = qr.abundance / (length / 1000.0)
        else:
            rpk[cid] = 0.0

    # Normalize so sum = 1e6
    total_rpk = sum(rpk.values())
    if total_rpk > 0:
        scaling = 1e6 / total_rpk
        return {cid: val * scaling for cid, val in rpk.items()}
    else:
        return {cid: 0.0 for cid in results}


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
                    "chrom": qr.chrom,
                    "strand": qr.strand,
                    "start": qr.start,
                    "end": qr.end,
                    "exons": qr.exons,
                    "gene_id": qr.gene_id,
                    "coherence_sum": 0.0,
                    "discrimination_sum": 0.0,
                    "combined_sum": 0.0,
                    "score_weight": 0.0,
                    "breakpoint_left": qr.breakpoint_left,
                    "breakpoint_right": qr.breakpoint_right,
                }
            a = agg[qr.candidate_id]
            a["abundance"] += qr.abundance
            a["num_assigned_reads"] += qr.num_assigned_reads
            if qr.num_assigned_reads > 0:
                a["confidence_sum"] += qr.confidence * qr.num_assigned_reads
                a["confidence_count"] += qr.num_assigned_reads
            # Read-count-weighted score aggregation; fall back to unit weight
            # for intervals with zero assigned reads so we still record the score.
            score_weight = float(qr.num_assigned_reads) if qr.num_assigned_reads > 0 else 1.0
            a["coherence_sum"] += qr.coherence_score * score_weight
            a["discrimination_sum"] += qr.discrimination_score * score_weight
            a["combined_sum"] += qr.combined_score * score_weight
            a["score_weight"] += score_weight
            if a["breakpoint_left"] is None:
                a["breakpoint_left"] = qr.breakpoint_left
            if a["breakpoint_right"] is None:
                a["breakpoint_right"] = qr.breakpoint_right

    result = {}
    for cid, a in agg.items():
        confidence = (
            a["confidence_sum"] / a["confidence_count"]
            if a["confidence_count"] > 0
            else 0.0
        )
        weight = a["score_weight"]
        if weight > 0.0:
            coherence = a["coherence_sum"] / weight
            discrimination = a["discrimination_sum"] / weight
            combined = a["combined_sum"] / weight
        else:
            coherence = 0.0
            discrimination = 0.0
            combined = 0.0
        result[cid] = QuantResult(
            candidate_id=cid,
            abundance=a["abundance"],
            confidence=confidence,
            num_assigned_reads=a["num_assigned_reads"],
            source=a["source"],
            chrom=a["chrom"],
            strand=a["strand"],
            start=a["start"],
            end=a["end"],
            exons=a["exons"],
            gene_id=a["gene_id"],
            coherence_score=coherence,
            discrimination_score=discrimination,
            combined_score=combined,
            breakpoint_left=a["breakpoint_left"],
            breakpoint_right=a["breakpoint_right"],
        )

    return result
