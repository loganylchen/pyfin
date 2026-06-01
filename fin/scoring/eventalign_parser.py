"""Parse f5c eventalign TSV output into scoring matrices."""

from __future__ import annotations

import csv
import logging
import math
from dataclasses import dataclass, field
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

    # AC7-pre: per-event records — (cDNA reference_position, signal_start_idx,
    # signal_end_idx). Opt-in via parse_eventalign_tsv(collect_events=True)
    # to protect real-data runs from a 2-3 GB memory regression.
    events: List[Tuple[int, int, int]] = field(default_factory=list)

    # Gap-based distance: weighted gap runs with context quality.
    gap_distance: float = 0.0
    # Outlier-based distance: outlier_fraction + max_bad_run_normalized.
    outlier_distance: float = 0.0

    # Internal accumulators
    _positions: Set[int] = field(default_factory=set, repr=False)
    _squared_errors: List[float] = field(default_factory=list, repr=False)
    _position_ll: Dict[int, float] = field(default_factory=dict, repr=False)
    _outlier_count: int = field(default=0, repr=False)
    _cur_bad_run: int = field(default=0, repr=False)
    _max_bad_run: int = field(default=0, repr=False)

    def finalize(self, candidate_length: int, gap_context_n: int = 5):
        """Compute final metrics from accumulated events."""
        self.candidate_length = candidate_length
        if candidate_length > 0:
            self.coverage = len(self._positions) / candidate_length
        if self._squared_errors:
            self.event_rmse = math.sqrt(
                sum(self._squared_errors) / len(self._squared_errors)
            )
        self.gap_distance = self._compute_gap_distance(gap_context_n)
        # Outlier distance: fraction + normalized max run.
        if self.num_events > 0:
            frac = self._outlier_count / self.num_events
            run_norm = self._max_bad_run / self.num_events
            self.outlier_distance = frac + run_norm
        else:
            self.outlier_distance = 0.0

    def _compute_gap_distance(self, context_n: int = 5) -> float:
        """Weighted gap-run distance with flanking-context confidence.

        For each run of consecutive unaligned positions within the mapped span:
          weight = gap_length * context_confidence
        where context_confidence = exp(mean_LL of flanking events), clamped.

        Returns sum of weights / mapped_span. Higher = more structural gaps.
        """
        if len(self._positions) < 2:
            return 0.0

        sorted_pos = sorted(self._positions)
        min_pos, max_pos = sorted_pos[0], sorted_pos[-1]
        mapped_span = max_pos - min_pos + 1
        if mapped_span <= 1:
            return 0.0

        # Find gap runs (consecutive positions without events)
        covered = self._positions
        gap_runs: List[Tuple[int, int]] = []  # (start, end_exclusive)
        i = min_pos
        while i <= max_pos:
            if i not in covered:
                gap_start = i
                while i <= max_pos and i not in covered:
                    i += 1
                gap_runs.append((gap_start, i))
            else:
                i += 1

        if not gap_runs:
            return 0.0

        total_score = 0.0
        for gap_start, gap_end in gap_runs:
            gap_len = gap_end - gap_start

            # Flanking context: LL of N events before and after gap
            left_lls = [
                self._position_ll[p]
                for p in range(gap_start - context_n, gap_start)
                if p in self._position_ll
            ]
            right_lls = [
                self._position_ll[p]
                for p in range(gap_end, gap_end + context_n)
                if p in self._position_ll
            ]
            flanking = left_lls + right_lls

            if flanking:
                mean_ll = sum(flanking) / len(flanking)
                # exp(mean_ll) ∈ (0,1): good fit → high confidence gap is real
                ctx_conf = math.exp(max(mean_ll, -10.0))
            else:
                # Terminal gap — low confidence
                ctx_conf = 0.01

            total_score += gap_len * ctx_conf

        return total_score / mapped_span


def parse_eventalign_tsv(
    tsv_path: str,
    candidate_lengths: Optional[Dict[str, int]] = None,
    collect_events: bool = False,
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
        collect_events: When True, append per-event (position, signal_start_idx,
            signal_end_idx) tuples to ReadCandidateScore.events. Off by default
            for memory safety (real-data runs can produce 5-20k events/read).
            Only the diff-region DTW m4 path requires this.

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
                _process_row(header, accum, collect_events=collect_events)

        for row in reader:
            if not row or len(row) < 15:
                continue
            _process_row(row, accum, collect_events=collect_events)

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
    collect_events: bool = False,
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
        # Track per-position LL for gap context computation
        score._position_ll[position] = ll
        # Outlier tracking (|z| > 2.0)
        if abs(z) > 2.0:
            score._outlier_count += 1
            score._cur_bad_run += 1
            if score._cur_bad_run > score._max_bad_run:
                score._max_bad_run = score._cur_bad_run
        else:
            score._cur_bad_run = 0

    # Squared error for RMSE
    score._squared_errors.append((event_level_mean - model_mean) ** 2)

    # Update signal index bounds
    score.signal_start_idx = min(score.signal_start_idx, start_idx)
    score.signal_end_idx = max(score.signal_end_idx, end_idx)

    # AC7-pre: opt-in per-event record collection for diff-region DTW.
    if collect_events:
        score.events.append((position, start_idx, end_idx))


def build_distance_matrix(
    scores: List[ReadCandidateScore],
    read_ids: List[str],
    candidate_ids: List[str],
    metric: str = "nll",
) -> np.ndarray:
    """Build a distance matrix from eventalign scores.

    The matrix is compatible with em_with_coherence() input format.

    Args:
        scores: Parsed ReadCandidateScore list.
        read_ids: Ordered list of read IDs (rows).
        candidate_ids: Ordered list of candidate IDs (columns).
        metric: Distance metric. ``"nll"`` = mean negative log-likelihood
            per event (default, legacy). ``"gap"`` = weighted gap-run
            distance with flanking-context confidence. ``"outlier"`` =
            z-score outlier fraction + max consecutive bad run (normalized).

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
            if metric == "gap":
                dist[ri, ci] = s.gap_distance
            elif metric == "outlier":
                dist[ri, ci] = s.outlier_distance
            else:
                # Mean negative log-likelihood per event (legacy).
                n_events = max(int(s.num_events), 1)
                dist[ri, ci] = -s.total_log_likelihood / n_events

    return dist
