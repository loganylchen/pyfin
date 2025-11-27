"""
Event alignment module for aligning nanopore signals to sequences.

This module implements adaptive banded event alignment algorithm similar to f5c and nanopolish.

Based on the ABEA (Adaptive Banded Event Alignment) algorithm from:
https://github.com/hasindu2008/f5c
"""

import numpy as np
from typing import List, Tuple, Optional, Dict, Any
from dataclasses import dataclass
import logging

logger = logging.getLogger(__name__)


@dataclass
class Event:
    """Nanopore signal event with statistical features."""
    mean: float
    stdv: float
    start: int
    length: int
    samples: Optional[np.ndarray] = None

    def __post_init__(self):
        if self.samples is None:
            self.samples = np.array([])


@dataclass
class AlignedEvent:
    """An event aligned to a specific k-mer position."""
    ref_position: int
    ref_kmer: str
    event_idx: int
    event_mean: float
    event_stdv: float
    event_length: int
    strand: str
    alignment_score: float
    model_mean: float
    model_stdv: float
    state: int  # 0=MATCH, 1=STAY, 2=SKIP


@dataclass
class KmerModel:
    """Pore model for k-mer current levels."""
    k: int = 5
    levels: Optional[Dict[str, float]] = None
    variances: Optional[Dict[str, float]] = None
    name: str = "r9.4_5mer"

    def __post_init__(self):
        if self.levels is None:
            self.levels = {}
        if self.variances is None:
            self.variances = {}

    def get_level(self, kmer: str) -> float:
        """Get expected current level for a k-mer."""
        if kmer in self.levels:
            return self.levels[kmer]

        # Fallback: simple model based on GC content
        gc_content = (kmer.count('G') + kmer.count('C')) / len(kmer)
        return 80.0 + gc_content * 20.0

    def get_stdv(self, kmer: str) -> float:
        """Get expected standard deviation for a k-mer."""
        if kmer in self.variances:
            return self.variances[kmer]
        return 2.0


class AdaptiveBandedAligner:
    """
    Adaptive Banded Event Aligner (ABEA).

    Aligns nanopore signal events to reference sequence using banded dynamic programming.
    Based on the algorithm from f5c (https://github.com/hasindu2008/f5c).
    """

    def __init__(
        self,
        band_width: int = 100,
        kmer_size: int = 5,
        lp_step: float = 0.0,      # Log probability of step transition
        lp_stay: float = -0.5,     # Log probability of stay transition
        lp_skip: float = -3.0,     # Log probability of skip transition
    ):
        """
        Initialize the aligner.

        Args:
            band_width: Width of the band in DP matrix
            kmer_size: Size of k-mers (typically 5 for DNA, 6 for RNA004)
            lp_step: Log probability of normal progression (step)
            lp_stay: Log probability of event splitting (stay)
            lp_skip: Log probability of skipping (very rare)
        """
        self.band_width = band_width
        self.kmer_size = kmer_size
        self.lp_step = lp_step
        self.lp_stay = lp_stay
        self.lp_skip = lp_skip

        self.model = KmerModel(k=kmer_size)

    @staticmethod
    def log_normal_pdf(x: float, mean: float, stdv: float) -> float:
        """
        Log probability density of normal distribution.

        Args:
            x: Observed value
            mean: Mean of distribution
            stdv: Standard deviation

        Returns:
            Log probability density
        """
        if stdv <= 0:
            return -np.inf

        # Constant term: -0.5 * log(2 * pi)
        log_const = -0.9189385332046727

        # Calculate log probability
        diff = x - mean
        log_prob = log_const - np.log(stdv) - (diff * diff) / (2 * stdv * stdv)

        return log_prob

    @staticmethod
    def coordinate_to_band(e_idx: int, k_idx: int, band_lower_left: List[Tuple[int, int]]) -> int:
        """
        Convert (event_idx, kmer_idx) coordinate to band index.

        Args:
            e_idx: Event index
            k_idx: Kmer index
            band_lower_left: Mapping from band to coordinates

        Returns:
            Band index if within band, -1 otherwise
        """
        # Find band that contains this coordinate
        for band_idx, (be, bk) in enumerate(band_lower_left):
            # Check if (e_idx, k_idx) is in this band
            if be <= e_idx < be + len(band_lower_left) and \
               k_idx == bk + (e_idx - be):
                return band_idx
        return -1

    def generate_kmers(self, sequence: str) -> List[str]:
        """
        Generate overlapping kmers from sequence.

        Args:
            sequence: Reference sequence

        Returns:
            List of kmers
        """
        kmers = []
        for i in range(len(sequence) - self.kmer_size + 1):
            kmer = sequence[i:i + self.kmer_size]
            kmers.append(kmer)
        return kmers

    def _calculate_suzuki_kasahara(
        self,
        score_diag: float,
        score_up: float,
        score_left: float
    ) -> Tuple[float, int]:
        """
        Suzuki-Kasahara band movement rule.

        Determines direction based on scores from previous band.

        Args:
            score_diag: Score from diagonal (step)
            score_up: Score from above (stay)
            score_left: Score from left (skip)

        Returns:
            Tuple of (max_score, direction)
            Direction: 0=diag/step, 1=up/stay, 2=left/skip
        """
        if score_diag >= score_up and score_diag >= score_left:
            return score_diag, 0
        elif score_up >= score_left:
            return score_up, 1
        else:
            return score_left, 2

    def align(
        self,
        events: List[Event],
        sequence: str,
        is_rna: bool = False
    ) -> List[AlignedEvent]:
        """
        Align events to sequence using adaptive banded DP.

        Args:
            events: List of detected events
            sequence: Reference sequence
            is_rna: Whether this is RNA (requires reverse complement and signal flipping)

        Returns:
            List of aligned events
        """
        if not events:
            logger.warning("No events provided for alignment")
            return []

        if len(sequence) < self.kmer_size:
            logger.warning(f"Sequence too short: {len(sequence)} < {self.kmer_size}")
            return []

        n_events = len(events)

        # Handle RNA: reverse complement sequence and reverse events
        if is_rna:
            sequence = self._reverse_complement(sequence)
            events = events[::-1]  # Reverse event order

        # Generate kmers from sequence
        kmers = self.generate_kmers(sequence)
        n_kmers = len(kmers)

        logger.debug(f"Aligning {n_events} events to {n_kmers} kmers")

        # Initialize adaptive band
        band_lower_left = self._init_band(n_events, n_kmers)
        n_bands = len(band_lower_left)

        # Initialize DP matrices
        band_scores = np.full((n_bands, self.band_width), -np.inf, dtype=float)
        trace_table = np.zeros((n_bands, self.band_width), dtype=np.int8)

        # Fill DP band matrix
        self._fill_dp_matrix(
            events, kmers, band_lower_left,
            band_scores, trace_table
        )

        # Backtrace to get alignment
        aligned_events = self._backtrace(
            events, kmers, sequence,
            band_lower_left, band_scores, trace_table,
            is_rna
        )

        return aligned_events

    def _reverse_complement(self, seq: str) -> str:
        """Compute reverse complement of sequence."""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                      'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
        return ''.join(complement.get(base, base) for base in reversed(seq))

    def _init_band(self, n_events: int, n_kmers: int) -> List[Tuple[int, int]]:
        """
        Initialize adaptive band lower-left coordinates.

        The band follows a diagonal path from (0, 0) to (n_events, n_kmers)
        with adaptive width to handle length variations.

        Args:
            n_events: Number of events
            n_kmers: Number of kmers

        Returns:
            List of (event_idx, kmer_idx) for each band
        """
        band_lower_left = [(0, 0)]

        slope = n_kmers / n_events if n_events > 0 else 1.0

        for i in range(1, min(n_events, n_kmers) + self.band_width):
            # Start near diagonal: kmer_idx ≈ slope * event_idx
            e_idx = min(i, n_events - 1)
            k_idx = min(int(slope * e_idx), n_kmers - 1)
            band_lower_left.append((e_idx, k_idx))

        # Add one more for edge case
        band_lower_left.append((n_events - 1, n_kmers - 1))

        return band_lower_left

    def _fill_dp_matrix(
        self,
        events: List[Event],
        kmers: List[str],
        band_lower_left: List[Tuple[int, int]],
        band_scores: np.ndarray,
        trace_table: np.ndarray
    ) -> None:
        """
        Fill the DP matrix using banded alignment.

        Args:
            events: List of events
            kmers: List of kmers
            band_lower_left: Band coordinates
            band_scores: DP score matrix
            trace_table: Traceback matrix
        """
        n_bands = band_scores.shape[0]

        # Initialize first band
        band_scores[0, 0] = 0.0

        for band_idx in range(1, n_bands):
            lower_left = band_lower_left[band_idx]
            e_idx = lower_left[0]
            k_idx = lower_left[1]

            # Calculate band boundaries
            max_row = min(self.band_width, len(events) - e_idx)

            for row in range(max_row):
                curr_e_idx = e_idx + row
                curr_k_idx = k_idx + row

                if curr_e_idx >= len(events) or curr_k_idx >= len(kmers):
                    continue

                # Get current event and kmer
                event = events[curr_e_idx]
                kmer = kmers[curr_k_idx]

                # Calculate emission probability
                model_mean = self.model.get_level(kmer)
                model_stdv = self.model.get_stdv(kmer)
                lp_emission = self.log_normal_pdf(event.mean, model_mean, model_stdv)

                # Get scores from previous band positions
                score_diag = -np.inf
                score_up = -np.inf
                score_left = -np.inf

                if band_idx >= 2 and row >= 1:
                    score_diag = band_scores[band_idx - 2, row - 1]
                if band_idx >= 1:
                    score_up = band_scores[band_idx - 1, row]
                if band_idx >= 1 and row >= 1:
                    score_left = band_scores[band_idx - 1, row - 1]

                # Calculate transition scores
                score_step = score_diag + lp_emission + self.lp_step if score_diag > -np.inf else -np.inf
                score_stay = score_up + lp_emission + self.lp_stay if score_up > -np.inf else -np.inf
                score_skip = score_left + lp_emission + self.lp_skip if score_left > -np.inf else -np.inf

                # Choose best path
                if score_step >= score_stay and score_step >= score_skip:
                    band_scores[band_idx, row] = score_step
                    trace_table[band_idx, row] = 0  # DIAG/STEP
                elif score_stay >= score_skip:
                    band_scores[band_idx, row] = score_stay
                    trace_table[band_idx, row] = 1  # UP/STAY
                else:
                    band_scores[band_idx, row] = score_skip
                    trace_table[band_idx, row] = 2  # LEFT/SKIP

    def _backtrace(
        self,
        events: List[Event],
        kmers: List[str],
        sequence: str,
        band_lower_left: List[Tuple[int, int]],
        band_scores: np.ndarray,
        trace_table: np.ndarray,
        is_rna: bool
    ) -> List[AlignedEvent]:
        """
        Backtrace through DP matrix to extract alignment.

        Args:
            events: List of events
            kmers: List of kmers
            sequence: Reference sequence
            band_lower_left: Band coordinates
            band_scores: DP score matrix
            trace_table: Traceback matrix
            is_rna: Whether this is RNA

        Returns:
            List of aligned events
        """
        aligned_events = []

        # Start from best score at end
        n_bands = band_scores.shape[0]
        band_idx = n_bands - 1

        # Find best position in last band
        best_row = np.argmax(band_scores[band_idx, :])
        row = best_row

        lower_left = band_lower_left[band_idx]
        e_idx = lower_left[0] + row
        k_idx = lower_left[1] + row

        total_score = band_scores[band_idx, row]

        logger.debug(f"Backtrace from band {band_idx}, row {row}, score {total_score}")

        while band_idx > 0 and row >= 0:
            lower_left = band_lower_left[band_idx]
            e_idx = lower_left[0] + min(row, len(events) - lower_left[0] - 1)
            k_idx = lower_left[1] + min(row, len(kmers) - lower_left[1] - 1)

            if e_idx >= len(events) or k_idx >= len(kmers):
                break

            direction = trace_table[band_idx, row]

            if direction == 0:  # DIAG/STEP
                # Create aligned event
                kmer = kmers[k_idx]
                event = events[e_idx]

                aligned_event = AlignedEvent(
                    ref_position=k_idx,
                    ref_kmer=kmer,
                    event_idx=e_idx,
                    event_mean=event.mean,
                    event_stdv=event.stdv,
                    event_length=event.length,
                    strand="-" if is_rna else "+",
                    alignment_score=band_scores[band_idx, row],
                    model_mean=self.model.get_level(kmer),
                    model_stdv=self.model.get_stdv(kmer),
                    state=0
                )
                aligned_events.append(aligned_event)

                # Move diagonally
                band_idx -= 2
                row -= 1

            elif direction == 1:  # UP/STAY
                # Stay: same kmer, previous event
                band_idx -= 1

            elif direction == 2:  # LEFT/SKIP
                # Skip: same event, previous kmer
                band_idx -= 1

            else:
                break

        # Reverse to get forward order
        aligned_events.reverse()

        return aligned_events
