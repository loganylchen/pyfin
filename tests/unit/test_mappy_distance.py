"""Tests for fin.scoring.mappy_distance.compute_mappy_distance."""

from __future__ import annotations

import numpy as np
import pytest

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.scoring.mappy_distance import compute_mappy_distance

pytest.importorskip("mappy")


def _cand(cid: str, seq: str) -> TranscriptCandidate:
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=()),
        three_prime_pos=0,
        sequence=seq,
        source="gtf",
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=0,
        end=len(seq),
    )


# Long enough sequences for minimap2 map-ont preset to align.
SEQ_A = "ACGT" * 100  # 400bp
SEQ_B = "TGCA" * 100  # 400bp, distinct content


class TestShapeAndTypes:
    def test_empty_reads(self):
        d = compute_mappy_distance({}, [_cand("c1", SEQ_A)], [])
        assert d.shape == (0, 1)
        assert d.dtype == np.float32

    def test_empty_candidates(self):
        d = compute_mappy_distance({"r1": SEQ_A}, [], ["r1"])
        assert d.shape == (1, 0)

    def test_both_empty(self):
        d = compute_mappy_distance({}, [], [])
        assert d.shape == (0, 0)


class TestDistanceProperties:
    def test_distance_non_negative(self):
        cands = [_cand("c1", SEQ_A), _cand("c2", SEQ_B)]
        seqs = {"r1": SEQ_A, "r2": SEQ_B}
        d = compute_mappy_distance(seqs, cands, ["r1", "r2"])
        assert (d >= 0).all()

    def test_self_match_is_min_per_row(self):
        # A read identical to candidate j should have d[i,j] == row min.
        cands = [_cand("c1", SEQ_A), _cand("c2", SEQ_B)]
        seqs = {"r1": SEQ_A, "r2": SEQ_B}
        d = compute_mappy_distance(seqs, cands, ["r1", "r2"])
        # r1 most similar to c1: d[0,0] should be ≤ d[0,1]
        assert d[0, 0] <= d[0, 1]
        # r2 most similar to c2
        assert d[1, 1] <= d[1, 0]

    def test_missing_read_seq_yields_zero_row(self):
        cands = [_cand("c1", SEQ_A)]
        d = compute_mappy_distance({}, cands, ["r1"])
        np.testing.assert_array_equal(d[0], np.zeros(1, dtype=np.float32))

    def test_empty_candidate_sequence_yields_no_hits(self):
        cands = [_cand("c_empty", ""), _cand("c2", SEQ_A)]
        seqs = {"r1": SEQ_A}
        d = compute_mappy_distance(seqs, cands, ["r1"])
        # Aligner for empty cand is None → that column treated as missing.
        # Aligner for c2 hits → d[0,1] is row max (=0 here since only one finite).
        assert d.shape == (1, 2)
        assert d[0, 1] == 0.0  # only hit, equals row_max - row_max
        # Missing column gets (row_max - (row_min - 1)) = 1.0 distance penalty
        assert d[0, 0] > d[0, 1]

    def test_dtype_float32(self):
        d = compute_mappy_distance(
            {"r1": SEQ_A}, [_cand("c1", SEQ_A)], ["r1"]
        )
        assert d.dtype == np.float32
