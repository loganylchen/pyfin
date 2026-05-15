"""Tests for fin/scoring/diff_region_dtw.py.

Covers:
  - extract_diff_regions: basic diff detection, single candidate, merged regions
  - genomic_region_to_cdna: plus-strand, minus-strand, intronic overlap
  - cdna_region_to_signal_range: event lookup
  - compute_diff_region_m4: zero matrix on no diff regions (AC5), NaN handling
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval
from fin.scoring.diff_region_dtw import (
    cdna_region_to_signal_range,
    cluster_candidates_by_chain,
    compute_class_partitioned_m4,
    compute_diff_region_m4,
    extract_diff_regions,
    genomic_region_to_cdna,
)
from fin.scoring.eventalign_parser import ReadCandidateScore


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_candidate(
    cid: str,
    chrom: str,
    start: int,
    end: int,
    strand: str = "+",
    introns: Tuple[Tuple[int, int], ...] = (),
    seq: str = "",
) -> TranscriptCandidate:
    interval = GenomicInterval(chrom=chrom, start=start, end=end, strand=strand)
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=introns),
        three_prime_pos=end,
        sequence=seq or ("A" * (end - start)),
        source="novel",
        supporting_read_ids=frozenset(),
        chrom=chrom,
        strand=strand,
        start=start,
        end=end,
    )


def _make_score(
    read_name: str,
    cid: str,
    events: List[Tuple[int, int, int]],
) -> ReadCandidateScore:
    s = ReadCandidateScore(read_name=read_name, candidate_id=cid)
    s.events = events
    return s


# ---------------------------------------------------------------------------
# extract_diff_regions
# ---------------------------------------------------------------------------

class TestExtractDiffRegions:
    def test_fewer_than_two_candidates_returns_empty(self):
        c = _make_candidate("c1", "chr1", 100, 200)
        assert extract_diff_regions([c]) == []
        assert extract_diff_regions([]) == []

    def test_identical_single_exon_no_diff(self):
        """Two single-exon candidates with same span → no diff bases."""
        c1 = _make_candidate("c1", "chr1", 100, 300)
        c2 = _make_candidate("c2", "chr1", 100, 300)
        assert extract_diff_regions([c1, c2]) == []

    def test_full_intron_padded_window(self):
        """One cand single-exon, other has 100bp intron → window is the full
        unioned intron (150-250) padded ±pad_bp=6 → (144, 256)."""
        c1 = _make_candidate("c1", "chr1", 100, 300)
        c2 = _make_candidate("c2", "chr1", 100, 300, introns=((150, 250),))
        assert extract_diff_regions([c1, c2]) == [(144, 256)]

    def test_small_wobble_unioned_to_full_intron(self):
        """Donor wobble: c1 intron 150-250, c2 intron 153-250 → raw cluster
        is 3bp at 150-153, but unioning all overlapping introns expands to
        the full intron (150, 250), padded ±6 → (144, 256)."""
        c1 = _make_candidate("c1", "chr1", 100, 300, introns=((150, 250),))
        c2 = _make_candidate("c2", "chr1", 100, 300, introns=((153, 250),))
        regions = extract_diff_regions([c1, c2])
        assert regions == [(144, 256)]

    def test_pad_bp_parameter(self):
        """Custom pad_bp adjusts the flanking exonic context."""
        c1 = _make_candidate("c1", "chr1", 100, 300)
        c2 = _make_candidate("c2", "chr1", 100, 300, introns=((150, 250),))
        # pad_bp=0 → exact intron span
        assert extract_diff_regions([c1, c2], pad_bp=0) == [(150, 250)]
        # pad_bp=10 → padded ±10
        assert extract_diff_regions([c1, c2], pad_bp=10) == [(140, 260)]

    def test_pad_clipped_to_candidate_span(self):
        """Pad is clipped to the union of candidate spans."""
        c1 = _make_candidate("c1", "chr1", 100, 300)
        c2 = _make_candidate("c2", "chr1", 100, 300, introns=((102, 298),))
        # Intron union 102-298, pad ±6 would yield (96, 304) but clipped to
        # candidate-span union (100, 300).
        assert extract_diff_regions([c1, c2]) == [(100, 300)]

    def test_multiple_small_wobbles(self):
        """Two independent donor/acceptor wobbles → two separate windows."""
        c1 = _make_candidate("c1", "chr1", 0, 600, introns=((100, 200), (400, 500)))
        c2 = _make_candidate("c2", "chr1", 0, 600, introns=((104, 200), (400, 506)))
        regions = extract_diff_regions([c1, c2])
        # Cluster 1: raw 100-104, union introns (100,200)∪(104,200)=(100,200),
        # pad ±6 → (94, 206).
        # Cluster 2: raw 500-506, union introns (400,500)∪(400,506)=(400,506),
        # pad ±6 → (394, 512).
        assert regions == [(94, 206), (394, 512)]

    def test_adjacent_padded_windows_merge(self):
        """Two intron diffs whose padded windows touch should merge."""
        # Two introns 100-110 and 120-130; pad ±6 → (94, 116) and (114, 136)
        # → merge to (94, 136).
        c1 = _make_candidate("c1", "chr1", 0, 200, introns=((100, 110), (120, 130)))
        c2 = _make_candidate("c2", "chr1", 0, 200)
        regions = extract_diff_regions([c1, c2])
        assert regions == [(94, 136)]

    def test_no_diff_if_all_candidates_have_same_intron(self):
        """Both candidates have the same intron → no exon/intron split."""
        c1 = _make_candidate("c1", "chr1", 0, 500, introns=((100, 200),))
        c2 = _make_candidate("c2", "chr1", 0, 500, introns=((100, 200),))
        assert extract_diff_regions([c1, c2]) == []


# ---------------------------------------------------------------------------
# genomic_region_to_cdna
# ---------------------------------------------------------------------------

class TestGenomicRegionToCdna:
    def test_single_exon_full_overlap(self):
        """Single exon 100-300, region 100-300 → cDNA 0-200."""
        c = _make_candidate("c1", "chr1", 100, 300)
        result = genomic_region_to_cdna(c, (100, 300))
        assert result == (0, 200)

    def test_single_exon_partial_overlap(self):
        """Single exon 100-300, region 150-250 → cDNA 50-150."""
        c = _make_candidate("c1", "chr1", 100, 300)
        result = genomic_region_to_cdna(c, (150, 250))
        assert result == (50, 150)

    def test_no_overlap_returns_none(self):
        c = _make_candidate("c1", "chr1", 100, 300)
        assert genomic_region_to_cdna(c, (400, 500)) is None

    def test_intronic_region_returns_none(self):
        """Region entirely within intron → no exonic bases → None."""
        c = _make_candidate("c1", "chr1", 0, 500, introns=((100, 400),))
        assert genomic_region_to_cdna(c, (150, 350)) is None

    def test_two_exon_plus_strand(self):
        """Exons: 0-100, 200-300; region 0-100 → cDNA 0-100."""
        c = _make_candidate("c1", "chr1", 0, 300, introns=((100, 200),))
        result = genomic_region_to_cdna(c, (0, 100))
        assert result == (0, 100)

    def test_two_exon_plus_strand_second_exon(self):
        """Exons: 0-100, 200-300; region 200-300 → cDNA 100-200."""
        c = _make_candidate("c1", "chr1", 0, 300, introns=((100, 200),))
        result = genomic_region_to_cdna(c, (200, 300))
        assert result == (100, 200)

    def test_minus_strand_single_exon(self):
        """Minus strand single exon 0-100; region 0-100.
        cDNA is reverse-complemented: spliced_len=100, c_start_rc=100-100=0, c_end_rc=100-0=100."""
        c = _make_candidate("c1", "chr1", 0, 100, strand="-")
        result = genomic_region_to_cdna(c, (0, 100))
        assert result == (0, 100)

    def test_minus_strand_partial(self):
        """Minus strand single exon 0-100; region 60-100.
        Genomic c_start=60, c_end=100. RC: c_start_rc=100-100=0, c_end_rc=100-60=40."""
        c = _make_candidate("c1", "chr1", 0, 100, strand="-")
        result = genomic_region_to_cdna(c, (60, 100))
        assert result == (0, 40)


# ---------------------------------------------------------------------------
# cdna_region_to_signal_range
# ---------------------------------------------------------------------------

class TestCdnaRegionToSignalRange:
    def test_no_events_returns_none(self):
        s = ReadCandidateScore(read_name="r1", candidate_id="c1")
        assert cdna_region_to_signal_range(s, (0, 100)) is None

    def test_events_in_range(self):
        s = _make_score("r1", "c1", [(10, 100, 200), (20, 200, 300), (50, 300, 400)])
        result = cdna_region_to_signal_range(s, (0, 60))
        assert result == (100, 400)

    def test_events_outside_range_returns_none(self):
        s = _make_score("r1", "c1", [(100, 0, 50)])
        assert cdna_region_to_signal_range(s, (0, 50)) is None

    def test_boundary_exclusive_end(self):
        """pos == c_end is excluded (half-open interval)."""
        s = _make_score("r1", "c1", [(50, 0, 100)])
        assert cdna_region_to_signal_range(s, (0, 50)) is None

    def test_single_event_in_range(self):
        s = _make_score("r1", "c1", [(5, 10, 20)])
        result = cdna_region_to_signal_range(s, (0, 10))
        assert result == (10, 20)

    def test_zero_width_signal_range_returns_none(self):
        """sig_lo == sig_hi → None."""
        s = _make_score("r1", "c1", [(5, 50, 50)])
        assert cdna_region_to_signal_range(s, (0, 10)) is None


# ---------------------------------------------------------------------------
# compute_diff_region_m4
# ---------------------------------------------------------------------------

class TestComputeDiffRegionM4:
    def test_empty_read_ids_returns_empty_matrix(self):
        m4 = compute_diff_region_m4(
            read_ids=[],
            candidates=[],
            scores_by_pair={},
            signal_reader=None,
            interval_start=0,
            interval_end=1000,
        )
        assert m4.shape == (0, 0)

    def test_no_diff_regions_returns_zero_matrix(self):
        """AC5: when no diff regions, must return zeros (not NaN)."""
        # Two identical single-exon candidates → no diff regions.
        c1 = _make_candidate("c1", "chr1", 0, 500)
        c2 = _make_candidate("c2", "chr1", 0, 500)
        read_ids = ["r1", "r2"]
        m4 = compute_diff_region_m4(
            read_ids=read_ids,
            candidates=[c1, c2],
            scores_by_pair={},
            signal_reader=None,
            interval_start=0,
            interval_end=500,
        )
        assert m4.shape == (2, 2)
        assert np.all(m4 == 0.0), f"Expected all zeros, got:\n{m4}"

    def test_diagonal_is_zero(self):
        """Self-distance must always be 0."""
        c1 = _make_candidate("c1", "chr1", 0, 500)
        c2 = _make_candidate("c2", "chr1", 0, 500, introns=((100, 400),))
        # No signal reader → all segments will be empty → off-diagonal NaN.
        m4 = compute_diff_region_m4(
            read_ids=["r1", "r2"],
            candidates=[c1, c2],
            scores_by_pair={},
            signal_reader=None,
            interval_start=0,
            interval_end=500,
        )
        assert m4[0, 0] == 0.0
        assert m4[1, 1] == 0.0

    def test_single_read_returns_1x1_zero(self):
        c1 = _make_candidate("c1", "chr1", 0, 200)
        c2 = _make_candidate("c2", "chr1", 0, 200, introns=((50, 150),))
        m4 = compute_diff_region_m4(
            read_ids=["r1"],
            candidates=[c1, c2],
            scores_by_pair={},
            signal_reader=None,
            interval_start=0,
            interval_end=200,
        )
        assert m4.shape == (1, 1)
        assert m4[0, 0] == 0.0

    def test_any_cand_projection_when_host_region_is_intronic(self):
        """Regression: when a read's max-LL host puts the diff region inside
        an intron (so genomic→cDNA projection fails), we must fall back to
        ANY other candidate that DOES project the region. Without this fix,
        every off-diagonal cell stays NaN on intervals where the wrong
        candidate happens to win the LL race.

        Setup:
          c1 (host)  intron (150,160)  → diff region (150,153) is intronic
          c2 (alt)   intron (153,160)  → diff region (150,153) is exonic

        Both reads carry events for c2 only in the diff-region range; c1's
        higher total_log_likelihood would make it the "host" under the old
        single-host code path, which would then drop the region for being
        intronic. The any-cand path must instead project through c2 and
        compute a finite distance.
        """
        c1 = _make_candidate("c1", "chr1", 0, 300, introns=((150, 160),))
        c2 = _make_candidate("c2", "chr1", 0, 300, introns=((153, 160),))
        # Sanity: padded window covers full unioned intron (150, 160) ± 6.
        assert extract_diff_regions([c1, c2]) == [(144, 166)]

        # c1 would be the host (higher LL), but its cDNA frame skips bases
        # 150-160 (intron). Only c2 has the region in an exon and only c2
        # carries events covering the cDNA bases mapped from genomic 150-153.
        s_r1_c1 = ReadCandidateScore(read_name="r1", candidate_id="c1")
        s_r1_c1.total_log_likelihood = 100.0
        s_r1_c1.events = []  # no events for c1
        s_r1_c2 = ReadCandidateScore(read_name="r1", candidate_id="c2")
        s_r1_c2.total_log_likelihood = 10.0
        # genomic 150-153 → c2 cDNA 150-153 (first exon 0-153). Events at
        # cDNA pos 150, 151, 152 land in the region with signal samples 0..30.
        s_r1_c2.events = [(150, 0, 10), (151, 10, 20), (152, 20, 30)]

        s_r2_c1 = ReadCandidateScore(read_name="r2", candidate_id="c1")
        s_r2_c1.total_log_likelihood = 100.0
        s_r2_c1.events = []
        s_r2_c2 = ReadCandidateScore(read_name="r2", candidate_id="c2")
        s_r2_c2.total_log_likelihood = 10.0
        s_r2_c2.events = [(150, 5, 15), (151, 15, 25), (152, 25, 35)]

        scores = {
            ("r1", "c1"): s_r1_c1, ("r1", "c2"): s_r1_c2,
            ("r2", "c1"): s_r2_c1, ("r2", "c2"): s_r2_c2,
        }

        class _FakeSignalReader:
            def __init__(self, signals):
                self._signals = signals

            def get_picoamp_signal(self, rid):
                sig = self._signals.get(rid)
                return (sig, None) if sig is not None else None

        # Two distinct signals long enough to slice [0..35]. DTW on these
        # will give a finite (non-NaN) distance.
        rng = np.random.default_rng(0)
        sig_r1 = rng.normal(0.0, 1.0, 100).astype(np.float32)
        sig_r2 = rng.normal(0.5, 1.0, 100).astype(np.float32)
        reader = _FakeSignalReader({"r1": sig_r1, "r2": sig_r2})

        m4 = compute_diff_region_m4(
            read_ids=["r1", "r2"],
            candidates=[c1, c2],
            scores_by_pair=scores,
            signal_reader=reader,
            interval_start=0,
            interval_end=300,
            use_gpu=False,
        )
        assert m4.shape == (2, 2)
        # Diagonal is always 0.
        assert m4[0, 0] == 0.0 and m4[1, 1] == 0.0
        # Off-diagonal MUST be finite (this is the bug-fix assertion).
        assert np.isfinite(m4[0, 1]), f"any-cand path failed: m4={m4}"
        assert np.isfinite(m4[1, 0])
        # Symmetric.
        assert m4[0, 1] == m4[1, 0]

    def test_cross_class_entries_are_zero(self):
        """Reads assigned to different classes must have M[i,j]=0 even if
        signal DTW would otherwise produce a distance. Cross-class info is
        handled by the chain partition, not by signal DTW."""
        # Two distinct classes (different intron counts).
        c_a = _make_candidate("cA", "chr1", 0, 500)  # class 0: single exon
        c_b = _make_candidate("cB", "chr1", 0, 500, introns=((100, 200),))  # class 1
        # Add a wobble partner to class 1 so it has intra-class diff regions.
        c_b2 = _make_candidate("cB2", "chr1", 0, 500, introns=((103, 200),))

        s_rA = ReadCandidateScore(read_name="rA", candidate_id="cA")
        s_rA.total_log_likelihood = 100.0
        s_rA.events = []
        s_rB = ReadCandidateScore(read_name="rB", candidate_id="cB")
        s_rB.total_log_likelihood = 100.0
        s_rB.events = []
        scores = {("rA", "cA"): s_rA, ("rB", "cB"): s_rB}

        m4 = compute_diff_region_m4(
            read_ids=["rA", "rB"],
            candidates=[c_a, c_b, c_b2],
            scores_by_pair=scores,
            signal_reader=None,
            interval_start=0,
            interval_end=500,
            use_gpu=False,
        )
        # rA in class 0, rB in class 1 → cross-class pair = 0.
        assert m4[0, 1] == 0.0
        assert m4[1, 0] == 0.0

    def test_singleton_class_is_zero(self):
        """A class with a single candidate has no intra-class disagreement
        and must yield zero off-diagonal entries for its reads."""
        c1 = _make_candidate("c1", "chr1", 0, 500)  # singleton class
        s_r1 = ReadCandidateScore(read_name="r1", candidate_id="c1")
        s_r1.total_log_likelihood = 100.0
        s_r2 = ReadCandidateScore(read_name="r2", candidate_id="c1")
        s_r2.total_log_likelihood = 100.0
        scores = {("r1", "c1"): s_r1, ("r2", "c1"): s_r2}

        m4 = compute_diff_region_m4(
            read_ids=["r1", "r2"],
            candidates=[c1],
            scores_by_pair=scores,
            signal_reader=None,
            interval_start=0,
            interval_end=500,
            use_gpu=False,
        )
        assert m4.shape == (2, 2)
        assert np.all(m4 == 0.0)


# ---------------------------------------------------------------------------
# cluster_candidates_by_chain
# ---------------------------------------------------------------------------

class TestClusterCandidatesByChain:
    def test_empty(self):
        assert cluster_candidates_by_chain([]) == []

    def test_all_single_exon_one_class(self):
        c1 = _make_candidate("c1", "chr1", 0, 100)
        c2 = _make_candidate("c2", "chr1", 0, 100)
        c3 = _make_candidate("c3", "chr1", 0, 100)
        classes = cluster_candidates_by_chain([c1, c2, c3])
        assert classes == [[0, 1, 2]]

    def test_single_exon_vs_spliced_separate(self):
        c1 = _make_candidate("c1", "chr1", 0, 300)
        c2 = _make_candidate("c2", "chr1", 0, 300, introns=((100, 200),))
        classes = cluster_candidates_by_chain([c1, c2])
        # Different intron counts → different classes.
        assert sorted(map(sorted, classes)) == [[0], [1]]

    def test_wobble_groups_within_tolerance(self):
        """±6bp wobble at donor and acceptor → same class at chain_wobble=6."""
        c1 = _make_candidate("c1", "chr1", 0, 500, introns=((100, 300),))
        c2 = _make_candidate("c2", "chr1", 0, 500, introns=((104, 296),))
        c3 = _make_candidate("c3", "chr1", 0, 500, introns=((120, 320),))  # outside
        classes = cluster_candidates_by_chain([c1, c2, c3], chain_wobble=6)
        # c1 & c2 cluster; c3 is its own class.
        flat = sorted(map(sorted, classes))
        assert [0, 1] in flat
        assert [2] in flat

    def test_zero_wobble_strict(self):
        c1 = _make_candidate("c1", "chr1", 0, 500, introns=((100, 300),))
        c2 = _make_candidate("c2", "chr1", 0, 500, introns=((101, 300),))
        classes = cluster_candidates_by_chain([c1, c2], chain_wobble=0)
        assert sorted(map(sorted, classes)) == [[0], [1]]


# ---------------------------------------------------------------------------
# compute_class_partitioned_m4
# ---------------------------------------------------------------------------

class _FakeSignalReader:
    """Returns a deterministic synthetic signal per read.

    The signal shape is unique per read so DTW distances are finite and
    non-zero between different reads (avoiding the zero-distance shortcut
    that would mask matrix structure).
    """

    def __init__(self, n_per_read: int = 64, seed: int = 0):
        self.n = n_per_read
        self.seed = seed

    def _make(self, read_id: str) -> np.ndarray:
        # Stable hash from the read_id so signals are reproducible.
        rng = np.random.default_rng(abs(hash((read_id, self.seed))) % (2**32))
        base = rng.standard_normal(self.n).astype(np.float32) * 5.0 + 100.0
        return base

    def get_picoamp_signal(self, read_id):
        return self._make(read_id), {}

    def get_calibrated_signal(self, read_id):
        return self._make(read_id), {}


class TestComputeClassPartitionedM4:
    def test_empty_reads(self):
        c1 = _make_candidate("c1", "chr1", 0, 200)
        m4 = compute_class_partitioned_m4(
            read_ids=[],
            candidates=[c1],
            cand_ids=["c1"],
            dist_read_cand=np.zeros((0, 1)),
            signal_reader=None,
        )
        assert m4.shape == (0, 0)

    def test_empty_candidates(self):
        m4 = compute_class_partitioned_m4(
            read_ids=["r1", "r2"],
            candidates=[],
            cand_ids=[],
            dist_read_cand=np.zeros((2, 0)),
            signal_reader=None,
        )
        assert m4.shape == (2, 2)
        assert np.all(m4 == 0.0)

    def test_single_class_all_intra(self):
        """All candidates same chain → one class → all off-diag finite & equal-
        fill (max-of-intra trivially = intra distances themselves; no inter)."""
        c1 = _make_candidate("c1", "chr1", 0, 500, introns=((100, 300),))
        c2 = _make_candidate("c2", "chr1", 0, 500, introns=((103, 300),))  # wobble
        read_ids = ["r1", "r2", "r3"]
        dist = np.array([
            [0.0, 1.0],   # r1 prefers c1
            [1.0, 0.0],   # r2 prefers c2
            [0.5, 0.5],   # r3 either
        ], dtype=np.float32)
        m4 = compute_class_partitioned_m4(
            read_ids=read_ids,
            candidates=[c1, c2],
            cand_ids=["c1", "c2"],
            dist_read_cand=dist,
            signal_reader=_FakeSignalReader(),
            use_gpu=False,
        )
        assert m4.shape == (3, 3)
        # Diagonal zero.
        assert m4[0, 0] == 0.0 and m4[1, 1] == 0.0 and m4[2, 2] == 0.0
        # Symmetric.
        assert np.allclose(m4, m4.T)
        # All off-diagonal are intra-class DTW results (finite, > 0 for
        # distinct synthetic signals).
        off = m4[~np.eye(3, dtype=bool)]
        assert np.all(np.isfinite(off))

    def test_multi_class_inter_is_constant_max(self):
        """Two classes; inter-class cells must equal max of intra-class
        distances (default inter_class_fill='max')."""
        # Class A: single-exon
        c1 = _make_candidate("c1", "chr1", 0, 500)
        c2 = _make_candidate("c2", "chr1", 0, 500)
        # Class B: spliced
        c3 = _make_candidate("c3", "chr1", 0, 500, introns=((100, 300),))
        c4 = _make_candidate("c4", "chr1", 0, 500, introns=((100, 300),))
        read_ids = ["rA1", "rA2", "rB1", "rB2"]
        # rA* prefer c1/c2 (class 0); rB* prefer c3/c4 (class 1).
        dist = np.array([
            [0.0, 0.1, 5.0, 5.1],
            [0.2, 0.0, 5.0, 5.2],
            [5.0, 5.0, 0.0, 0.1],
            [5.0, 5.0, 0.2, 0.0],
        ], dtype=np.float32)
        m4 = compute_class_partitioned_m4(
            read_ids=read_ids,
            candidates=[c1, c2, c3, c4],
            cand_ids=["c1", "c2", "c3", "c4"],
            dist_read_cand=dist,
            signal_reader=_FakeSignalReader(),
            use_gpu=False,
            inter_class_fill="max",
        )
        # Intra-class distances: m4[0,1] (class A), m4[2,3] (class B).
        intra = [m4[0, 1], m4[2, 3]]
        expected_fill = float(max(intra))
        # Inter-class cells:
        for i, j in [(0, 2), (0, 3), (1, 2), (1, 3)]:
            assert m4[i, j] == pytest.approx(expected_fill, rel=1e-6)
            assert m4[j, i] == pytest.approx(expected_fill, rel=1e-6)
        # Symmetric + zero diagonal.
        assert np.allclose(m4, m4.T)
        for k in range(4):
            assert m4[k, k] == 0.0

    def test_multi_class_inter_is_user_constant(self):
        c1 = _make_candidate("c1", "chr1", 0, 500)
        c2 = _make_candidate("c2", "chr1", 0, 500, introns=((100, 300),))
        read_ids = ["r1", "r2"]
        dist = np.array([[0.0, 5.0], [5.0, 0.0]], dtype=np.float32)
        m4 = compute_class_partitioned_m4(
            read_ids=read_ids,
            candidates=[c1, c2],
            cand_ids=["c1", "c2"],
            dist_read_cand=dist,
            signal_reader=_FakeSignalReader(),
            use_gpu=False,
            inter_class_fill=42.0,
        )
        assert m4[0, 1] == pytest.approx(42.0)
        assert m4[1, 0] == pytest.approx(42.0)
        assert m4[0, 0] == 0.0 and m4[1, 1] == 0.0

    def test_unassigned_read_treated_as_inter(self):
        """Read with all-NaN distance row → class -1 → inter-class with all
        other reads (fill constant)."""
        c1 = _make_candidate("c1", "chr1", 0, 500)
        c2 = _make_candidate("c2", "chr1", 0, 500)
        read_ids = ["r1", "r2", "r3"]
        dist = np.array([
            [0.0, 0.1],
            [0.2, 0.0],
            [np.nan, np.nan],   # r3 unassigned
        ], dtype=np.float32)
        m4 = compute_class_partitioned_m4(
            read_ids=read_ids,
            candidates=[c1, c2],
            cand_ids=["c1", "c2"],
            dist_read_cand=dist,
            signal_reader=_FakeSignalReader(),
            use_gpu=False,
            inter_class_fill=7.0,
        )
        # r1 vs r2 are same-class → intra DTW result (finite).
        assert np.isfinite(m4[0, 1])
        # r3 vs anyone is "inter" because r3 is unassigned (class -1).
        assert m4[0, 2] == pytest.approx(7.0)
        assert m4[1, 2] == pytest.approx(7.0)
        assert m4[2, 0] == pytest.approx(7.0)
        assert m4[2, 1] == pytest.approx(7.0)
        # Diagonal zero.
        assert m4[2, 2] == 0.0

    def test_no_signal_reader_returns_constant_fill_for_inter(self):
        """No reader → no intra distances computed → inter fill defaults to 0
        (max over empty list)."""
        c1 = _make_candidate("c1", "chr1", 0, 500)
        c2 = _make_candidate("c2", "chr1", 0, 500, introns=((100, 300),))
        read_ids = ["r1", "r2"]
        dist = np.array([[0.0, 5.0], [5.0, 0.0]], dtype=np.float32)
        m4 = compute_class_partitioned_m4(
            read_ids=read_ids,
            candidates=[c1, c2],
            cand_ids=["c1", "c2"],
            dist_read_cand=dist,
            signal_reader=None,
            use_gpu=False,
        )
        assert m4.shape == (2, 2)
        assert np.all(m4 == 0.0)
