"""Tests for transcript quantification."""

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.analysis.quantification import (
    QuantResult,
    aggregate_across_intervals,
    quantify_transcripts,
)
from fin.candidates.dataclasses import IntronChain, TranscriptCandidate


def _make_candidate(cid, source="gtf"):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=()),
        three_prime_pos=1000,
        sequence="ACGT",
        source=source,
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=0,
        end=1000,
    )


class TestQuantifyTranscripts:
    """Tests for probability-weighted quantification."""

    def test_uniform_assignment(self):
        """With uniform R, abundance should be equal."""
        R = np.array([[0.5, 0.5], [0.5, 0.5]])
        hard = np.array([0, 1])
        candidates = [_make_candidate("tx1"), _make_candidate("tx2")]
        read_ids = ["r1", "r2"]

        results = quantify_transcripts(R, hard, candidates, read_ids)

        assert len(results) == 2
        assert results[0].abundance == 1.0
        assert results[1].abundance == 1.0

    def test_hard_assignment(self):
        """With hard assignment R, abundance should reflect assignment."""
        R = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 0.0]])
        hard = np.array([0, 1, 0])
        candidates = [_make_candidate("tx1"), _make_candidate("tx2")]
        read_ids = ["r1", "r2", "r3"]

        results = quantify_transcripts(R, hard, candidates, read_ids)

        assert results[0].abundance == 2.0
        assert results[0].num_assigned_reads == 2
        assert results[0].confidence == 1.0
        assert results[1].abundance == 1.0
        assert results[1].num_assigned_reads == 1

    def test_soft_assignment(self):
        """With soft R, abundance is probability-weighted."""
        R = np.array([[0.8, 0.2], [0.3, 0.7]])
        hard = np.array([0, 1])
        candidates = [_make_candidate("tx1"), _make_candidate("tx2")]
        read_ids = ["r1", "r2"]

        results = quantify_transcripts(R, hard, candidates, read_ids)

        np.testing.assert_almost_equal(results[0].abundance, 1.1)  # 0.8 + 0.3
        np.testing.assert_almost_equal(results[1].abundance, 0.9)  # 0.2 + 0.7
        assert results[0].confidence == 0.8  # only r1 assigned to tx1
        assert results[1].confidence == 0.7  # only r2 assigned to tx2

    def test_no_assigned_reads(self):
        """Candidate with no hard-assigned reads should have 0 confidence."""
        R = np.array([[0.1, 0.9]])
        hard = np.array([1])
        candidates = [_make_candidate("tx1"), _make_candidate("tx2")]
        read_ids = ["r1"]

        results = quantify_transcripts(R, hard, candidates, read_ids)

        assert results[0].num_assigned_reads == 0
        assert results[0].confidence == 0.0


class TestAggregateAcrossIntervals:
    """Tests for cross-interval aggregation."""

    def test_same_candidate_aggregated(self):
        """Same candidate in two intervals should sum abundance."""
        interval1 = [
            QuantResult("tx1", abundance=3.0, confidence=0.9, num_assigned_reads=3, source="gtf"),
        ]
        interval2 = [
            QuantResult("tx1", abundance=2.0, confidence=0.8, num_assigned_reads=2, source="gtf"),
        ]

        agg = aggregate_across_intervals([interval1, interval2])

        assert "tx1" in agg
        assert agg["tx1"].abundance == 5.0
        assert agg["tx1"].num_assigned_reads == 5
        # Weighted average confidence: (0.9*3 + 0.8*2) / 5 = 4.3/5 = 0.86
        np.testing.assert_almost_equal(agg["tx1"].confidence, 0.86)

    def test_different_candidates(self):
        """Different candidates should stay separate."""
        interval1 = [
            QuantResult("tx1", abundance=3.0, confidence=0.9, num_assigned_reads=3, source="gtf"),
            QuantResult("tx2", abundance=1.0, confidence=0.7, num_assigned_reads=1, source="novel"),
        ]

        agg = aggregate_across_intervals([interval1])

        assert len(agg) == 2
        assert agg["tx1"].abundance == 3.0
        assert agg["tx2"].abundance == 1.0
        assert agg["tx2"].source == "novel"
