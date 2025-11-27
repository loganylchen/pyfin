"""
Completeness scoring for eventalign results.

Provides functionality to calculate alignment completeness scores
for read-isoform alignments, measuring how much of the isoform
sequence is covered by aligned events.
"""

import numpy as np
from typing import List, Dict, Any, Set, Tuple, Optional
from dataclasses import dataclass
import logging

from fin.core.eventalign import AlignedEvent

logger = logging.getLogger(__name__)


@dataclass
class CompletenessResult:
    """
    Result of completeness calculation.

    Attributes:
        completeness_score: Proportion of isoform covered by events (0.0-1.0)
        aligned_positions: Set of isoform positions that have alignments
        total_events: Total number of events considered
        num_aligned: Number of events that aligned
        mean_event_quality: Average quality of aligned events
        coverage_distribution: Position-by-position coverage counts
        percentile_90_coverage: 90th percentile of coverage
    """
    completeness_score: float = 0.0
    aligned_positions: Set[int] = None
    total_events: int = 0
    num_aligned: int = 0
    mean_event_quality: float = 0.0
    coverage_distribution: Optional[Dict[int, int]] = None
    percentile_90_coverage: float = 0.0

    def __post_init__(self):
        """Initialize set after creation."""
        if self.aligned_positions is None:
            self.aligned_positions = set()
        if self.coverage_distribution is None:
            self.coverage_distribution = {}

    @property
    def is_high_completeness(self) -> bool:
        """Check if this is high completeness (>0.8)."""
        return self.completeness_score >= 0.8

    @property
    def alignment_proportion(self) -> float:
        """Return proportion of events that aligned."""
        if self.total_events == 0:
            return 0.0
        return self.num_aligned / self.total_events

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "completeness_score": self.completeness_score,
            "aligned_bases": len(self.aligned_positions),
            "total_events": self.total_events,
            "num_aligned": self.num_aligned,
            "mean_event_quality": self.mean_event_quality,
            "alignment_proportion": self.alignment_proportion
        }


def calculate_completeness(
    aligned_events: List[AlignedEvent],
    isoform_seq: str,
    min_event_duration: int = 1,
) -> CompletenessResult:
    """
    Calculate eventalign completeness score for read-isoform alignment.

    Completeness measures the proportion of isoform sequence that is
    covered by aligned events. Higher completeness indicates better
    alignment quality and alignment to the entire isoform.

    Args:
        aligned_events: List of AlignedEvent objects from eventalign
        isoform_seq: Full isoform nucleotide sequence
        min_event_duration: Minimum event duration (in samples) to include

    Returns:
        CompletenessResult with scores and metrics

    Examples:
        >>> aligned_events = align_to_sequence(events, isoform_seq)
        >>> result = calculate_completeness(aligned_events, isoform_seq)
        >>> print(f"Completeness: {result.completeness_score:.3f}")
        >>> print(f"Aligned positions: {len(result.aligned_positions)}")
    """
    if not isoform_seq:
        logger.warning("Empty isoform sequence provided")
        return CompletenessResult()

    if not aligned_events:
        logger.warning("No aligned events provided")
        return CompletenessResult(
            aligned_positions=set(),
            coverage_distribution={}
        )

    result = CompletenessResult()
    result.total_events = len(aligned_events)

    # Track aligned positions with coverage
    # Position 0-based indexing into isoform sequence
    coverage = np.zeros(len(isoform_seq), dtype=int)
    quality_scores = []

    # Filter events and map to positions
    for event in aligned_events:
        # Skip if event doesn't meet quality criteria
        if event.kmer_length != len(event.reference_kmer):
            continue

        if event.reference_kmer is None or len(event.reference_kmer) == 0:
            continue

        # Check event duration filter
        if hasattr(event, 'duration') and event.duration < min_event_duration:
            continue

        try:
            # Find position in isoform sequence
            # reference_position is the position in the reference (isoform)
            # We align based on the kmer start position
            pos = event.reference_position

            if pos is None or pos < 0:
                continue

            # Add coverage for each base in the kmer
            kmer_length = len(event.reference_kmer)
            end_pos = min(pos + kmer_length, len(isoform_seq))

            if pos < len(coverage):
                coverage[pos:end_pos] += 1

                # Track aligned base positions
                for i in range(pos, end_pos):
                    if i < len(isoform_seq):
                        result.aligned_positions.add(i)

                result.num_aligned += 1

                # Track quality
                if hasattr(event, 'mean') or hasattr(event, 'stdv'):
                    quality_scores.append(1.0)

        except (AttributeError, TypeError) as e:
            logger.debug(f"Error processing event: {e}")
            continue

    # Calculate completeness score
    if len(isoform_seq) > 0:
        result.completeness_score = len(result.aligned_positions) / len(isoform_seq)

    # Store coverage distribution
    unique_coverage = np.unique(coverage)
    for cov in unique_coverage:
        count = np.sum(coverage == cov)
        result.coverage_distribution[int(cov)] = int(count)

    # Calculate 90th percentile coverage
    if len(coverage) > 0:
        result.percentile_90_coverage = float(np.percentile(coverage, 90))

    # Calculate mean quality
    if quality_scores:
        result.mean_event_quality = np.mean(quality_scores)

    return result


def filter_events_by_quality(
    aligned_events: List[AlignedEvent],
    min_mean_signal: Optional[float] = None,
    max_stdv: Optional[float] = None,
    min_duration: Optional[int] = None
) -> List[AlignedEvent]:
    """
    Filter aligned events by quality metrics.

    Args:
        aligned_events: List of AlignedEvent objects
        min_mean_signal: Minimum mean signal value
        max_stdv: Maximum standard deviation
        min_duration: Minimum duration in samples

    Returns:
        Filtered list of events
    """
    filtered = []

    for event in aligned_events:
        # Check mean signal
        if min_mean_signal is not None:
            if not hasattr(event, 'mean') or event.mean < min_mean_signal:
                continue

        # Check std dev
        if max_stdv is not None:
            if hasattr(event, 'stdv') and event.stdv > max_stdv:
                continue

        # Check duration
        if min_duration is not None:
            if hasattr(event, 'duration') and event.duration < min_duration:
                continue

        filtered.append(event)

    return filtered


def calculate_per_exon_completeness(
    aligned_events: List[AlignedEvent],
    isoform_seq: str,
    exon_bounds: List[Tuple[int, int]]
) -> Dict[int, float]:
    """
    Calculate completeness separately for each exon.

    Args:
        aligned_events: AlignedEvent list
        isoform_seq: Full isoform sequence
        exon_bounds: List of (start, end) for each exon

    Returns:
        Dictionary mapping exon_number (1-indexed) -> completeness score
    """
    if not aligned_events or not isoform_seq:
        return {}

    per_exon_scores = {}

    for exon_num, (exon_start, exon_end) in enumerate(exon_bounds, 1):
        # Create subset sequence for this exon
        if exon_start >= len(isoform_seq) or exon_end > len(isoform_seq):
            continue

        exon_seq = isoform_seq[exon_start:exon_end]

        # Find events that align to this exon
        exon_events = []
        for event in aligned_events:
            try:
                pos = event.reference_position
                if pos is None:
                    continue

                kmer_length = len(event.reference_kmer) if event.reference_kmer else 0
                event_end = pos + kmer_length

                # Check if event overlaps with exon
                if pos < exon_end and event_end > exon_start:
                    exon_events.append(event)
            except (AttributeError, TypeError):
                continue

        # Calculate completeness for this exon
        if exon_events:
            result = calculate_completeness(exon_events, exon_seq)
            per_exon_scores[exon_num] = result.completeness_score

    return per_exon_scores


def normalize_completeness_by_length(
    completeness_score: float,
    isoform_length: int,
    read_length: int
) -> float:
    """
    Normalize completeness score accounting for length differences.

    Prevents penalization of longer isoforms where only a subregion
    is expected to align.

    Args:
        completeness_score: Raw completeness (0-1)
        isoform_length: Length of isoform sequence
        read_length: Estimated length of read based on signal

    Returns:
        Normalized completeness score
    """
    if isoform_length == 0:
        return 0.0

    # If read is much shorter than isoform, normalize based on read coverage
    coverage_ratio = read_length / isoform_length if isoform_length > 0 else 0.0

    if coverage_ratio < 0.5:
        # For short reads aligning to long isoforms,
        # use completeness relative to read coverage
        expected_completeness = coverage_ratio
        if expected_completeness > 0:
            return min(1.0, completeness_score / expected_completeness)

    return completeness_score


def summarize_completeness_scores(
    completeness_scores: List[CompletenessResult]
) -> Dict[str, float]:
    """
    Summarize multiple completeness scores.

    Args:
        completeness_scores: List of CompletenessResult objects

    Returns:
        Dictionary with summary statistics
    """
    if not completeness_scores:
        return {}

    scores = [cs.completeness_score for cs in completeness_scores]
    aligned_bases = [len(cs.aligned_positions) for cs in completeness_scores]
    num_events = [cs.num_aligned for cs in completeness_scores]

    return {
        "mean_completeness": float(np.mean(scores)),
        "median_completeness": float(np.median(scores)),
        "std_completeness": float(np.std(scores)),
        "min_completeness": float(np.min(scores)),
        "max_completeness": float(np.max(scores)),
        "mean_aligned_bases": float(np.mean(aligned_bases)),
        "mean_events_per_read": float(np.mean(num_events)),
        "total_reads": len(completeness_scores)
    }


def compare_stranded_completeness(
    forward_completeness: CompletenessResult,
    reverse_completeness: CompletenessResult,
) -> Tuple[str, CompletenessResult]:
    """
    Compare completeness on both strands and return best.

    RNA from nanopore is sequenced 3'->5', causing strand orientation
    issues. This function helps determine the correct orientation.

    Args:
        forward_completeness: Completeness on forward strand
        reverse_completeness: Completeness on reverse strand

    Returns:
        Tuple of (strand, best CompletenessResult)
    """
    if forward_completeness.completeness_score >= reverse_completeness.completeness_score:
        return "+", forward_completeness
    else:
        return "-", reverse_completeness
