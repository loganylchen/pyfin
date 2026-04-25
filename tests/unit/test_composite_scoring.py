"""Unit tests for composite scoring functions (US-008, AC-06/07/08)."""

from __future__ import annotations

import numpy as np

from fin.scoring.composite import (
    CompositeScore,
    compute_coherence_score,
    compute_combined_score,
    compute_discrimination_score,
    score_candidates_composite,
)
from fin.candidates.dataclasses import IntronChain, TranscriptCandidate


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_candidate(cid: str) -> TranscriptCandidate:
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=()),
        three_prime_pos=1000,
        sequence="ACGT",
        source="novel",
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=0,
        end=100,
    )


# ---------------------------------------------------------------------------
# AC-06  coherence
# ---------------------------------------------------------------------------


class TestCoherenceScore:
    def test_coherence_tight_cluster(self):
        """Reads tightly clustered (small distances relative to global median) → coherence > 0.8."""
        # Build a matrix where candidate 0's reads (0,1,2) are very close together
        # but the global median is much larger (set by cross-cluster distances).
        n = 6
        dist_rr = np.ones((n, n)) * 5.0  # baseline large distances
        np.fill_diagonal(dist_rr, 0.0)
        # Within-cluster distances for reads 0-2 are tiny
        for i in range(3):
            for j in range(3):
                if i != j:
                    dist_rr[i, j] = 0.01

        # R: reads 0-2 fully assigned to candidate 0
        R = np.zeros((n, 2))
        R[:3, 0] = 1.0
        R[3:, 1] = 1.0

        score = compute_coherence_score(dist_rr, R, candidate_idx=0)
        assert score > 0.8, f"Expected > 0.8 for tight cluster, got {score}"

    def test_coherence_loose_cluster(self):
        """Reads loosely clustered (large within-cluster distances) → coherence < 0.4."""
        # Build a matrix where candidate 0's reads have large within-cluster distances
        # and the global median is small (set by tight cross-cluster distances).
        n = 6
        dist_rr = np.ones((n, n)) * 0.01  # baseline small distances
        np.fill_diagonal(dist_rr, 0.0)
        # Within-cluster distances for reads 0-2 are huge
        for i in range(3):
            for j in range(3):
                if i != j:
                    dist_rr[i, j] = 20.0

        R = np.zeros((n, 2))
        R[:3, 0] = 1.0
        R[3:, 1] = 1.0

        score = compute_coherence_score(dist_rr, R, candidate_idx=0)
        assert score < 0.4, f"Expected < 0.4 for loose cluster, got {score}"

    def test_coherence_empty_cluster(self):
        """No reads assigned to candidate → returns 0.0."""
        n = 4
        dist_rr = np.ones((n, n))
        np.fill_diagonal(dist_rr, 0.0)

        # All weight on candidate 1, none on 0
        R = np.zeros((n, 2))
        R[:, 1] = 1.0

        score = compute_coherence_score(dist_rr, R, candidate_idx=0)
        assert score == 0.0, f"Expected 0.0 for empty cluster, got {score}"


# ---------------------------------------------------------------------------
# AC-07  discrimination
# ---------------------------------------------------------------------------


class TestDiscriminationScore:
    def test_discrimination_clear_winner(self):
        """Candidate 0 much closer than all others → discrimination > 0.7."""
        n_reads, n_tx = 5, 3
        dist = np.full((n_reads, n_tx), 2.0)
        dist[:, 0] = 0.1  # candidate 0 is clearly best

        score = compute_discrimination_score(dist, candidate_idx=0)
        assert score > 0.7, f"Expected > 0.7 for clear winner, got {score}"

    def test_discrimination_tied(self):
        """All candidates equidistant → discrimination ≈ 0.5 (within 0.1)."""
        n_reads, n_tx = 10, 3
        dist = np.ones((n_reads, n_tx))  # all equal

        score = compute_discrimination_score(dist, candidate_idx=0)
        assert abs(score - 0.5) < 0.1, f"Expected ≈ 0.5 for tied, got {score}"

    def test_discrimination_single_candidate(self):
        """Only 1 candidate (n_tx == 1) → returns exactly 0.5."""
        dist = np.array([[0.3], [0.7], [0.1]])  # shape (3, 1)

        score = compute_discrimination_score(dist, candidate_idx=0)
        assert score == 0.5, f"Expected exactly 0.5 for single candidate, got {score}"


# ---------------------------------------------------------------------------
# AC-08  combined
# ---------------------------------------------------------------------------


class TestCombinedScore:
    def test_combined_score_formula(self):
        """combined(0.8, 0.9, alpha=0.5) ≈ sqrt(0.8 * 0.9) ≈ 0.8485."""
        expected = np.sqrt(0.8 * 0.9)
        result = compute_combined_score(0.8, 0.9, alpha=0.5)
        assert abs(result - expected) < 1e-4, f"Expected {expected:.6f}, got {result:.6f}"

    def test_combined_alpha_extremes(self):
        """alpha=1.0 → equals coherence; alpha=0.0 → equals discrimination."""
        coh, dis = 0.7, 0.4

        r_coh = compute_combined_score(coh, dis, alpha=1.0)
        assert abs(r_coh - coh) < 1e-9, f"alpha=1.0 should equal coherence, got {r_coh}"

        r_dis = compute_combined_score(coh, dis, alpha=0.0)
        assert abs(r_dis - dis) < 1e-9, f"alpha=0.0 should equal discrimination, got {r_dis}"

    def test_combined_bounds(self):
        """Random inputs in [0, 1] → output in [0, 1]."""
        rng = np.random.default_rng(42)
        for _ in range(200):
            c = float(rng.uniform(0, 1))
            d = float(rng.uniform(0, 1))
            a = float(rng.uniform(0, 1))
            result = compute_combined_score(c, d, alpha=a)
            assert 0.0 <= result <= 1.0, (
                f"combined({c:.3f}, {d:.3f}, alpha={a:.3f}) = {result} out of [0,1]"
            )


# ---------------------------------------------------------------------------
# score_candidates_composite
# ---------------------------------------------------------------------------


class TestScoreCandidatesComposite:
    def _make_arrays(self, n_reads: int = 6, n_tx: int = 3):
        rng = np.random.default_rng(0)
        dist_rr = rng.uniform(0.1, 1.0, (n_reads, n_reads))
        dist_rr = (dist_rr + dist_rr.T) / 2.0
        np.fill_diagonal(dist_rr, 0.0)

        dist_rt = rng.uniform(0.1, 2.0, (n_reads, n_tx))

        R = rng.dirichlet(np.ones(n_tx), size=n_reads)

        return dist_rr, dist_rt, R

    def test_score_candidates_composite_basic(self):
        """3 candidates → 3 CompositeScore objects with matching candidate_ids."""
        candidates = [_make_candidate(f"tx{i}") for i in range(3)]
        dist_rr, dist_rt, R = self._make_arrays(n_reads=6, n_tx=3)

        results = score_candidates_composite(candidates, dist_rt, dist_rr, R)

        assert len(results) == 3
        for i, cs in enumerate(results):
            assert isinstance(cs, CompositeScore)
            assert cs.candidate_id == f"tx{i}"
            assert 0.0 <= cs.coherence <= 1.0
            assert 0.0 <= cs.discrimination <= 1.0
            assert 0.0 <= cs.combined <= 1.0

    def test_score_candidates_composite_use_gpu_kwarg_accepted(self):
        """use_gpu=True must not raise even when cupy is absent."""
        candidates = [_make_candidate("tx0")]
        dist_rr = np.array([[0.0]])
        dist_rt = np.array([[0.5]])
        R = np.array([[1.0]])

        # Should not raise regardless of cupy availability
        results = score_candidates_composite(
            candidates, dist_rt, dist_rr, R, use_gpu=True
        )
        assert len(results) == 1
        assert results[0].candidate_id == "tx0"
