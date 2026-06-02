"""Probability-weighted transcript quantification from EM assignments."""

from __future__ import annotations

import logging
from collections import defaultdict
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
    # P0-6: read IDs hard-assigned to this candidate. Used by
    # aggregate_across_intervals to dedup reads that span multiple intervals.
    assigned_read_ids: Tuple[str, ...] = ()
    # T8: maximum EM responsibility any single read assigned to this transcript.
    # max_R = R[:, j].max() after EM. Used for FP-by-EM analysis.
    max_R: float = 0.0


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
    read_ids_arr = np.asarray(read_ids)
    results = []
    for j, cand in enumerate(candidates):
        assigned_mask = hard_assignments == j
        num_assigned = int(assigned_mask.sum())
        confidence = float(R[assigned_mask, j].mean()) if num_assigned > 0 else 0.0
        assigned_ids = (
            tuple(read_ids_arr[assigned_mask].tolist())
            if num_assigned > 0
            else ()
        )

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
                assigned_read_ids=assigned_ids,
            )
        )

    return results


def _qr_overlap(a: QuantResult, b: QuantResult) -> bool:
    """1-D genomic interval overlap on (chrom, strand, [start, end))."""
    if a.chrom != b.chrom or a.strand != b.strand:
        return False
    return a.start < b.end and b.start < a.end


def isoform_fraction_drops(
    results: Dict[str, QuantResult],
    min_fraction: float,
) -> set:
    """Return candidate_ids to drop by the minimum-isoform-fraction heuristic.

    For every NOVEL multi-exon transcript C, compute its locus-relative
    abundance::

        relabund(C) = C.abundance / max(abundance over the NOVEL multi-exon
                       transcripts that OVERLAP C, C included)

    and drop C when ``relabund(C) < min_fraction``. This is the standard
    Cufflinks ``--min-isoform-fraction`` / StringTie ``-f`` minor-isoform
    suppression rule: the low-fraction tail of a locus is enriched for
    incompletely-spliced precursors (pre-mRNA), RT/template-switching
    artifacts, and assembly noise rather than genuine isoforms.

    GTF-passthrough (``source="gtf"``), fusion (``source="fusion"``), and
    single-exon (mono) candidates are EXEMPT and never dropped. ``min_fraction
    <= 0`` disables the filter (returns an empty set).

    SIRV WARNING: synthetic SIRV has no genuine low-abundance isoform tail, so
    its F1-optimal threshold (~0.4) is 4-40x more aggressive than real tools
    and would gut recall on a real transcriptome. Keep this conservative
    (literature default 0.01) and re-tune on real data.
    """
    if min_fraction <= 0.0:
        return set()
    # Bucket by (chrom, strand): only candidates on the same chromosome and
    # strand can share a locus, so the pairwise overlap scan stays within a
    # bucket instead of running all-vs-all across the whole genome. (Within a
    # locus-dense single chrom/strand it is still O(n^2); a sweep-line could be
    # added if a real run ever makes that the bottleneck.)
    buckets: Dict[Tuple[str, str], List[QuantResult]] = defaultdict(list)
    for qr in results.values():
        if qr.source == "novel" and len(qr.exons) >= 2:
            buckets[(qr.chrom, qr.strand)].append(qr)
    drops: set = set()
    for bucket in buckets.values():
        for qr in bucket:
            locus_max = qr.abundance
            for other in bucket:
                if other.candidate_id != qr.candidate_id and _qr_overlap(qr, other):
                    if other.abundance > locus_max:
                        locus_max = other.abundance
            # The locus-dominant isoform (relabund == 1.0) is NEVER dropped,
            # even if a caller mis-configures min_fraction > 1.0.
            if qr.abundance >= locus_max:
                continue
            if locus_max > 0 and (qr.abundance / locus_max) < min_fraction:
                drops.add(qr.candidate_id)
    return drops


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
                    "abundance_sum": 0.0,
                    "abundance_weight": 0.0,  # for dedup-weighted average
                    "num_assigned_total": 0,  # legacy fallback when no read IDs
                    "confidence_sum": 0.0,
                    "confidence_count": 0,
                    "assigned_read_ids": set(),  # P0-6: union across intervals
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
                    "max_R": 0.0,
                }
            a = agg[qr.candidate_id]
            # Preserve max EM responsibility across intervals (a candidate's
            # max_R is the highest single-read responsibility it ever sees).
            if qr.max_R > a["max_R"]:
                a["max_R"] = qr.max_R
            # P0-6: track per-interval abundance and dedup read IDs across
            # intervals. Naive sum double-counts reads that overlap interval
            # boundaries; we instead average per-read abundance across the
            # intervals that observe each read.
            a["abundance_sum"] += qr.abundance
            a["abundance_weight"] += max(qr.num_assigned_reads, 1)
            a["num_assigned_total"] += qr.num_assigned_reads
            a["assigned_read_ids"].update(qr.assigned_read_ids)
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
        # P0-6: dedup read counts and rescale abundance by the dedup ratio.
        # If the same read is hard-assigned to this candidate in two
        # intervals, naive sum reports 2x the abundance. We rescale by
        # n_unique / n_observations so abundance per unique read stays
        # consistent. When callers don't populate `assigned_read_ids`
        # (legacy fixtures, intervals with no hard assignments), fall
        # back to the legacy sum-based behavior.
        unique_ids = a["assigned_read_ids"]
        num_assigned_unique = len(unique_ids)
        num_assigned_total = a["num_assigned_total"]
        ab_weight = a["abundance_weight"]
        if num_assigned_unique > 0 and ab_weight > 0:
            # avg abundance per assigned-read observation, then scale by
            # unique read count so candidates only seen in one interval are
            # unaffected.
            abundance = a["abundance_sum"] * (num_assigned_unique / ab_weight)
            num_assigned_out = num_assigned_unique
        else:
            # Legacy / no-read-id fallback: keep naive sum.
            abundance = a["abundance_sum"]
            num_assigned_out = num_assigned_total
        result[cid] = QuantResult(
            candidate_id=cid,
            abundance=abundance,
            confidence=confidence,
            num_assigned_reads=num_assigned_out,
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
            assigned_read_ids=tuple(sorted(unique_ids)),
            max_R=a["max_R"],
        )

    return result
