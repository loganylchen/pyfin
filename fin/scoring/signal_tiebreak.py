"""Signal-based tiebreaker for ambiguous M1 EM assignments.

When M1 assigns equal distance (d=0) to multiple candidates for a read,
use diff-region signal quality to break the tie. Only applies to reads
where R_best < ambig_threshold — confident reads are untouched.

The key insight: don't compare global signal fit (which is diluted by shared
exons). Instead, compare signal fit ONLY at genomic regions where the tied
candidates DIFFER (unique exons, different junctions). This focuses signal
information on the most discriminative positions.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

from fin.candidates.dataclasses import TranscriptCandidate
from fin.scoring.eventalign_parser import ReadCandidateScore

logger = logging.getLogger(__name__)


def _get_exons(candidate: TranscriptCandidate) -> List[Tuple[int, int]]:
    """Derive exon intervals from intron_chain + start/end."""
    introns = candidate.intron_chain.introns if candidate.intron_chain else ()
    if not introns:
        return [(candidate.start, candidate.end)]
    exons = []
    # Before first intron
    exons.append((candidate.start, introns[0][0]))
    # Between introns
    for i in range(len(introns) - 1):
        exons.append((introns[i][1], introns[i + 1][0]))
    # After last intron
    exons.append((introns[-1][1], candidate.end))
    return [(s, e) for s, e in exons if e > s]


def _exon_set(candidate: TranscriptCandidate) -> Set[Tuple[int, int]]:
    """Return set of (start, end) exon intervals for a candidate."""
    return set(_get_exons(candidate))


def _unique_genomic_regions(
    cand_a: TranscriptCandidate,
    cand_b: TranscriptCandidate,
) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    """Find genomic regions unique to each candidate.

    Returns (unique_to_a, unique_to_b) as lists of (start, end) intervals.
    """
    exons_a = _exon_set(cand_a)
    exons_b = _exon_set(cand_b)

    # Exons in A but not in B
    unique_a = sorted(exons_a - exons_b)
    unique_b = sorted(exons_b - exons_a)

    return unique_a, unique_b


def _genomic_to_cdna_positions(
    exons: List[Tuple[int, int]],
    genomic_regions: List[Tuple[int, int]],
) -> Set[int]:
    """Map genomic intervals to cDNA (0-based) positions.

    Each exon contributes sequential cDNA positions. A genomic position
    within exon i at offset k maps to sum(exon_lengths[:i]) + k.
    """
    sorted_exons = sorted(exons)
    cdna_positions = set()
    cdna_offset = 0

    for ex_start, ex_end in sorted_exons:
        ex_len = ex_end - ex_start
        for g_start, g_end in genomic_regions:
            # Overlap between this exon and the genomic region
            ov_start = max(ex_start, g_start)
            ov_end = min(ex_end, g_end)
            if ov_start < ov_end:
                for pos in range(ov_start - ex_start, ov_end - ex_start):
                    cdna_positions.add(cdna_offset + pos)
        cdna_offset += ex_len

    return cdna_positions


def _diff_region_nll(
    score: ReadCandidateScore,
    diff_cdna_positions: Set[int],
) -> Tuple[float, int]:
    """Compute mean NLL at diff-region cDNA positions.

    Returns (mean_nll, n_events_at_diff). If no events at diff positions,
    returns (penalty, 0) — the read doesn't cover the diff region at all,
    which is strong evidence against this candidate.
    """
    if not diff_cdna_positions:
        return 0.0, 0

    lls = []
    for pos in diff_cdna_positions:
        if pos in score._position_ll:
            lls.append(score._position_ll[pos])

    if not lls:
        # No events at diff positions → read doesn't cover this candidate's
        # unique region → strong evidence this is the wrong candidate
        return 10.0, 0  # high penalty

    mean_nll = -sum(lls) / len(lls)  # negative LL → positive distance
    return mean_nll, len(lls)


def signal_tiebreak(
    R: np.ndarray,
    read_ids: List[str],
    candidates: List[TranscriptCandidate],
    scores_by_pair: Dict[Tuple[str, str], ReadCandidateScore],
    ambig_threshold: float = 0.90,
) -> np.ndarray:
    """Update R matrix by resolving M1 ties using diff-region signal quality.

    For each read where R_best < ambig_threshold:
      1. Find tied candidates (those with R[i,j] > 0.1 * R_best)
      2. Find genomic regions unique to each tied candidate
      3. Compute signal NLL at those unique regions
      4. Assign read to candidate with lowest diff-region NLL

    Args:
        R: (n_reads, n_cands) responsibility matrix from M1 EM.
        read_ids: Read names corresponding to R rows.
        candidates: Ordered candidate list corresponding to R columns.
        scores_by_pair: (read_name, candidate_id) → ReadCandidateScore with
            _position_ll populated.
        ambig_threshold: Reads with R_best below this are considered ambiguous.

    Returns:
        Updated R matrix (copy). Confident reads unchanged.
    """
    R_new = R.copy()
    n_reads, n_cands = R.shape

    if n_cands < 2:
        return R_new

    # Pre-compute candidate exons
    cand_exons = []
    for c in candidates:
        cand_exons.append(_get_exons(c))

    n_tiebroken = 0
    n_ambiguous = 0

    for i in range(n_reads):
        best_R = R[i].max()
        if best_R >= ambig_threshold:
            continue

        n_ambiguous += 1
        rid = read_ids[i]

        # Find tied candidates: those with significant R share
        tied_mask = R[i] > 0.1 * best_R
        tied_indices = np.where(tied_mask)[0]
        if len(tied_indices) < 2:
            continue

        # For each pair of tied candidates, find diff regions and score
        # Strategy: for each tied candidate, sum up its "unique region NLL"
        # across all pairwise comparisons with other tied candidates
        cand_diff_scores = {}  # cand_idx → (total_diff_nll, total_diff_events)

        for j in tied_indices:
            cid = candidates[j].candidate_id
            score = scores_by_pair.get((rid, cid))
            if score is None or not score._position_ll:
                cand_diff_scores[j] = (10.0, 0)  # no data → penalty
                continue

            # Aggregate unique-region NLL across all other tied candidates
            total_nll = 0.0
            total_events = 0
            n_comparisons = 0

            for k in tied_indices:
                if k == j:
                    continue
                unique_j, _ = _unique_genomic_regions(candidates[j], candidates[k])
                if not unique_j:
                    continue

                # Map unique genomic regions to cDNA positions on candidate j
                diff_cdna = _genomic_to_cdna_positions(cand_exons[j], unique_j)
                if not diff_cdna:
                    continue

                nll, n_ev = _diff_region_nll(score, diff_cdna)
                total_nll += nll
                total_events += n_ev
                n_comparisons += 1

            if n_comparisons > 0:
                cand_diff_scores[j] = (total_nll / n_comparisons, total_events)
            else:
                # No unique regions found (all candidates identical)
                cand_diff_scores[j] = (0.0, 0)

        if not cand_diff_scores:
            continue

        # Pick winner: lowest diff-region NLL (best signal match at unique regions)
        # Tie-break by number of events (more events = more evidence)
        scored = [(j, nll, n_ev) for j, (nll, n_ev) in cand_diff_scores.items()]

        # Sort by NLL (lower=better), then by -events (more=better)
        scored.sort(key=lambda x: (x[1], -x[2]))
        winner_j = scored[0][0]
        winner_nll = scored[0][1]
        second_nll = scored[1][1] if len(scored) > 1 else winner_nll

        # Only tiebreak if there's meaningful difference
        if winner_nll < second_nll - 0.01:
            # Assign full responsibility to winner
            R_new[i, :] = 0.0
            R_new[i, winner_j] = 1.0
            n_tiebroken += 1

    logger.info(
        "Signal tiebreak: %d ambiguous reads, %d tiebroken (%.1f%%)",
        n_ambiguous, n_tiebroken,
        100 * n_tiebroken / max(n_ambiguous, 1),
    )
    return R_new
