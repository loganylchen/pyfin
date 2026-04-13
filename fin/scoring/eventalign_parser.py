"""Parse f5c eventalign TSV output into scoring matrices."""

from __future__ import annotations

import csv
import logging
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

logger = logging.getLogger(__name__)


@dataclass
class ReadCandidateScore:
    """Aggregated eventalign score for one (read, candidate) pair."""

    read_name: str
    candidate_id: str
    total_log_likelihood: float = 0.0
    coverage: float = 0.0
    event_rmse: float = 0.0
    signal_start_idx: int = 0
    signal_end_idx: int = 0
    num_events: int = 0
    candidate_length: int = 0

    # Internal accumulators
    _positions: Set[int] = field(default_factory=set, repr=False)
    _squared_errors: List[float] = field(default_factory=list, repr=False)

    def finalize(self, candidate_length: int):
        """Compute final metrics from accumulated events."""
        self.candidate_length = candidate_length
        if candidate_length > 0:
            self.coverage = len(self._positions) / candidate_length
        if self._squared_errors:
            self.event_rmse = math.sqrt(
                sum(self._squared_errors) / len(self._squared_errors)
            )


def parse_eventalign_tsv(
    tsv_path: str,
    candidate_lengths: Optional[Dict[str, int]] = None,
) -> List[ReadCandidateScore]:
    """Parse f5c eventalign 16-column TSV into per-(read, candidate) scores.

    Expected columns (tab-separated):
        contig, position, reference_kmer, read_name, strand, event_index,
        event_level_mean, event_stdv, event_length, model_kmer, model_mean,
        model_stdv, standardized_level, start_idx, end_idx, samples

    Args:
        tsv_path: Path to eventalign TSV file.
        candidate_lengths: Optional dict of candidate_id -> sequence length
            for computing coverage.

    Returns:
        List of ReadCandidateScore, one per (read, candidate) pair.
    """
    # Accumulate by (read_name, contig)
    accum: Dict[Tuple[str, str], ReadCandidateScore] = {}

    with open(tsv_path, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)

        # Detect header row
        if header and header[0].lower().startswith("contig"):
            pass  # skip header
        else:
            # No header, process first line
            if header:
                _process_row(header, accum)

        for row in reader:
            if not row or len(row) < 15:
                continue
            _process_row(row, accum)

    # Finalize scores
    results = []
    for score in accum.values():
        cand_len = 0
        if candidate_lengths and score.candidate_id in candidate_lengths:
            cand_len = candidate_lengths[score.candidate_id]
        score.finalize(cand_len)
        results.append(score)

    logger.info(
        "Parsed %d (read, candidate) pairs from %s", len(results), tsv_path
    )
    return results


def _process_row(
    row: List[str],
    accum: Dict[Tuple[str, str], ReadCandidateScore],
):
    """Process a single eventalign TSV row."""
    try:
        contig = row[0]
        position = int(row[1])
        read_name = row[3]
        event_level_mean = float(row[6])
        model_mean = float(row[10])
        model_stdv = float(row[11])
        start_idx = int(row[13])
        end_idx = int(row[14])
    except (IndexError, ValueError):
        return

    key = (read_name, contig)
    if key not in accum:
        accum[key] = ReadCandidateScore(
            read_name=read_name,
            candidate_id=contig,
            signal_start_idx=start_idx,
            signal_end_idx=end_idx,
        )

    score = accum[key]
    score.num_events += 1
    score._positions.add(position)

    # Log likelihood: log P(event | model) assuming Gaussian
    if model_stdv > 0:
        z = (event_level_mean - model_mean) / model_stdv
        ll = -0.5 * (z * z + math.log(2 * math.pi) + 2 * math.log(model_stdv))
        score.total_log_likelihood += ll

    # Squared error for RMSE
    score._squared_errors.append((event_level_mean - model_mean) ** 2)

    # Update signal index bounds
    score.signal_start_idx = min(score.signal_start_idx, start_idx)
    score.signal_end_idx = max(score.signal_end_idx, end_idx)


def build_distance_matrix(
    scores: List[ReadCandidateScore],
    read_ids: List[str],
    candidate_ids: List[str],
) -> np.ndarray:
    """Build a distance matrix from eventalign scores.

    The matrix is compatible with em_with_coherence() input format.
    Distance = -total_log_likelihood (lower is better for EM).

    Args:
        scores: Parsed ReadCandidateScore list.
        read_ids: Ordered list of read IDs (rows).
        candidate_ids: Ordered list of candidate IDs (columns).

    Returns:
        np.ndarray of shape (n_reads, n_candidates).
        Missing pairs get a large distance value.
    """
    read_idx = {rid: i for i, rid in enumerate(read_ids)}
    cand_idx = {cid: i for i, cid in enumerate(candidate_ids)}

    n_reads = len(read_ids)
    n_cands = len(candidate_ids)

    # Initialize with large distance
    dist = np.full((n_reads, n_cands), fill_value=1e6, dtype=np.float64)

    for s in scores:
        ri = read_idx.get(s.read_name)
        ci = cand_idx.get(s.candidate_id)
        if ri is not None and ci is not None:
            # Use negative log-likelihood as distance
            dist[ri, ci] = -s.total_log_likelihood

    return dist
