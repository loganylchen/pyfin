"""Tests for the production ``_quant_m2_em`` pure tie-break junction-NLL seed.

``_quant_m2_em`` (the production default quant engine) seeds ``em_with_coherence``
with a PURE tie-break junction-NLL distance: M1/AS only selects each read's
best-AS tie set, and the per-event junction NLL is the sole graded distance over
that tie set. The d_tx skeleton rules are:

  * cells outside a read's tie set            -> MISSING (1e6)
  * a unique best-AS read (tie of one)        -> 0.0 on its single candidate
  * a >=2 tie with NO scorable NLL            -> flat 0.0 across the tie (1/K split)
  * a scored tie cell                         -> (nll - lo)        [then /sigma2]
  * an unscorable cell inside a scored tie    -> (hi - lo) + PAD   [then /sigma2]

where sigma2 = max(median of the positive scored cells, 1e-3).

M3 read×read DTW coherence is opt-in (``m3_coherence``): default OFF means
``build_m3_coherence`` is never called and the EM runs with beta=0.

These tests mock the I/O boundaries (mappy multimap, the krill ``_tie_nll`` pass,
``em_with_coherence``, ``quantify_transcripts``) so only the seed math is exercised.
"""

from __future__ import annotations

from unittest.mock import patch

import numpy as np
import pytest

from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner
from fin.scoring.m2_junction_nll import MISSING

# Read ids (rows) and candidates (cols) used by every scenario below.
READ_IDS = ["r0", "r1", "r2", "r3"]
N_R = len(READ_IDS)
N_C = 3

# Per-read tie sets + scored junction-NLLs injected via the mocked _tie_nll.
#   r0: unique best-AS (tie of one) on candidate 1.
#   r1: 2-way tie {0,1}, both scored (nll 2.0 / 5.0).
#   r2: 3-way tie {0,1,2}, cands 0 & 2 scored (1.0 / 4.0), cand 1 unscorable.
#   r3: 2-way tie {0,1}, nothing scorable -> flat split.
TIES = {0: [1], 1: [0, 1], 2: [0, 1, 2], 3: [0, 1]}
NLLS = {1: {0: 2.0, 1: 5.0}, 2: {0: 1.0, 2: 4.0}}


def _cand(cid):
    seq = "ACGT" * 25
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=()),
        three_prime_pos=0,
        sequence=seq,
        source="novel",
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=0,
        end=len(seq),
    )


def _candidate_set():
    cands = [_cand(f"c{j}") for j in range(N_C)]
    read_seqs = {rid: "ACGT" * 25 for rid in READ_IDS}
    return CandidateSet(
        interval=GenomicInterval(chrom="chr1", start=0, end=100, strand="+"),
        candidates=cands,
        read_ids=set(read_seqs),
        read_sequences=read_seqs,
    )


def _interval():
    return GenomicInterval(chrom="chr1", start=0, end=100, strand="+")


class _FakeAligner:
    """mappy.Aligner stub: yields one dummy hit per read so the raw AS pass marks
    every read as kept (score_hit is patched to a positive constant in _run). The
    actual AS values are irrelevant because _tie_nll is mocked to drive ties."""

    def __init__(self, *a, **k):
        pass

    def map(self, *a, **k):
        return iter((object(),))


def _fake_tie_nll(self, kept_read_ids, read_seqs, cand_list, aligners, raw):
    n_ties = sum(1 for t in TIES.values() if len(t) >= 2)
    n_refined = len(NLLS)
    # cover_by_read is ignored on the gate-OFF path (these skeleton tests).
    return dict(NLLS), dict(TIES), n_ties, n_refined, {}


def _run(cfg, captured, m3_matrix=None, tie_nll_fn=_fake_tie_nll):
    """Run _quant_m2_em with all I/O boundaries mocked; capture EM inputs."""

    def _fake_em(**kwargs):
        captured.update(kwargs)
        n_r = kwargs["dist_read_to_tx"].shape[0]
        n_c = kwargs["dist_read_to_tx"].shape[1]
        return (
            np.zeros((n_r, n_c), dtype=np.float64),
            np.zeros(n_r, dtype=np.int64),
            None,
        )

    cs = _candidate_set()
    runner = PipelineRunner(cfg)

    with patch("mappy.Aligner", _FakeAligner), patch(
        "fin.scoring.mappy_score.score_hit", return_value=1.0,
    ), patch.object(PipelineRunner, "_tie_nll", tie_nll_fn), patch.object(
        PipelineRunner, "_eff_lengths", return_value=None
    ), patch(
        "fin.pipeline.runner.em_with_coherence", side_effect=_fake_em
    ), patch(
        "fin.pipeline.runner.quantify_transcripts", return_value=[]
    ), patch(
        "fin.scoring.m3_junction_coherence.build_m3_coherence",
        return_value=(m3_matrix if m3_matrix is not None
                      else np.zeros((N_R, N_R), dtype=np.float32)),
    ) as mock_m3:
        runner._quant_m2_em(cs, list(READ_IDS), _interval())
    return mock_m3


class TestPureSeedSkeleton:
    """The d_tx skeleton math (M3 off)."""

    def setup_method(self):
        self.captured = {}
        cfg = PipelineConfig(
            bam_path="/tmp/x.bam", quant_mode="m2_em", use_gpu=False,
            m2_diff_cover_gate=False,
        )
        self.mock_m3 = _run(cfg, self.captured)
        self.d_tx = np.asarray(self.captured["dist_read_to_tx"])
        # sigma2 = median of positive scored cells [3, 3, 4] = 3.0
        self.sigma2 = 3.0

    def test_outside_tie_is_missing(self):
        # r0 only ties candidate 1; 0 and 2 stay MISSING.
        assert self.d_tx[0, 0] == pytest.approx(MISSING)
        assert self.d_tx[0, 2] == pytest.approx(MISSING)
        # r1 ties {0,1}; candidate 2 stays MISSING.
        assert self.d_tx[1, 2] == pytest.approx(MISSING)

    def test_unique_best_is_zero(self):
        assert self.d_tx[0, 1] == pytest.approx(0.0)

    def test_scored_tie_cell_is_nll_minus_lo_over_sigma2(self):
        # r1: lo=2, hi=5 -> cell0=(2-2)=0, cell1=(5-2)=3 -> 3/3 = 1.0
        assert self.d_tx[1, 0] == pytest.approx(0.0)
        assert self.d_tx[1, 1] == pytest.approx(3.0 / self.sigma2)

    def test_unscorable_cell_in_scored_tie_is_hi_minus_lo_plus_pad(self):
        # r2: lo=1, hi=4 -> cell0=0, cell2=(4-1)=3 -> 3/3=1.0;
        # cell1 unscorable -> (hi-lo)+PAD = 3+1 = 4 -> 4/3.
        assert self.d_tx[2, 0] == pytest.approx(0.0)
        assert self.d_tx[2, 2] == pytest.approx(3.0 / self.sigma2)
        assert self.d_tx[2, 1] == pytest.approx(4.0 / self.sigma2)

    def test_unscorable_tie_is_flat_zero(self):
        # r3: no scorable NLL -> flat 0.0 across the tie (1/K split).
        assert self.d_tx[3, 0] == pytest.approx(0.0)
        assert self.d_tx[3, 1] == pytest.approx(0.0)

    def test_m3_off_means_no_dtw_and_beta_zero(self):
        self.mock_m3.assert_not_called()
        assert self.captured["beta"] == 0.0
        d_rr = np.asarray(self.captured["dist_read_to_read"])
        assert d_rr.shape == (N_R, N_R)
        assert not d_rr.any()
        assert self.captured["sigma"] == 1.0


# --- Diff-region coverage gate (gate ON): proportional redistribution --------
# Scenario 1 (build a covered ratio, then a fuzzy read follows it):
#   r0: tie{0,1} nlls 1/5 -> margin 4 >= 0.5, COVERED -> hard to 0; vote[0]+=1
#   r1: tie{0,1} nlls 1/5 -> margin 4 >= 0.5, COVERED -> hard to 0; vote[0]+=1
#   r2: tie{0,1} nlls 5/1 -> margin 4 >= 0.5, COVERED -> hard to 1; vote[1]+=1
#   => covered_vote = [2, 1, 0]
#   r3: tie{0,1} NO nlls -> fuzzy -> redistribute by [2,1] -> p=[2/3,1/3]
ON_TIES = {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1]}
ON_NLLS = {0: {0: 1.0, 1: 5.0}, 1: {0: 1.0, 1: 5.0}, 2: {0: 5.0, 1: 1.0}}
ON_COVER = {0: True, 1: True, 2: True}  # r3 absent -> fallback covered=True (unused: no nlls)


def _fake_tie_nll_on(self, kept_read_ids, read_seqs, cand_list, aligners, raw):
    return dict(ON_NLLS), dict(ON_TIES), 4, len(ON_NLLS), dict(ON_COVER)


# Scenario 2 (the prior is fed ONLY by covered+distinguishing reads):
#   r0: tie{0,1} nlls 2.0/2.05 -> margin 0.05 < 0.5, COVERED -> FLAT (no vote)
#   r1: tie{0,1} nlls 1/5 -> margin 4 >= 0.5, NOT covered -> hard to 0 (no vote)
#   r2: tie{0,1} NO nlls -> fuzzy; covered_vote over {0,1} == 0 -> FLAT fallback
#   r3: tie{0,1} NO nlls -> fuzzy -> FLAT fallback
ON2_TIES = {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1]}
ON2_NLLS = {0: {0: 2.0, 1: 2.05}, 1: {0: 1.0, 1: 5.0}}
ON2_COVER = {0: True, 1: False}


def _fake_tie_nll_on2(self, kept_read_ids, read_seqs, cand_list, aligners, raw):
    return dict(ON2_NLLS), dict(ON2_TIES), 4, len(ON2_NLLS), dict(ON2_COVER)


class TestDiffCoverGate:
    """Gate-ON decision matrix with proportional redistribution (no read drop)."""

    def setup_method(self):
        self.captured = {}
        cfg = PipelineConfig(
            bam_path="/tmp/x.bam", quant_mode="m2_em", use_gpu=False,
            m2_diff_cover_gate=True, m2_diff_cover_margin=0.5,
        )
        _run(cfg, self.captured, tie_nll_fn=_fake_tie_nll_on)
        self.d_tx = np.asarray(self.captured["dist_read_to_tx"])

    def test_no_reads_dropped(self):
        # Every read keeps its row -> recall-safe (no row subsetting).
        assert self.d_tx.shape == (N_R, N_C)

    def test_hard_assign_winner(self):
        # r0: margin 4.0 >= 0.5 -> full mass to candidate 0; others MISSING.
        assert self.d_tx[0, 0] == pytest.approx(0.0)
        assert self.d_tx[0, 1] == pytest.approx(MISSING)
        assert self.d_tx[0, 2] == pytest.approx(MISSING)
        # r2: hard to candidate 1.
        assert self.d_tx[2, 1] == pytest.approx(0.0)
        assert self.d_tx[2, 0] == pytest.approx(MISSING)

    def test_fuzzy_read_follows_covered_ratio(self):
        # r3: no signal -> redistribute by covered_vote [2,1] over tie {0,1}.
        # d_tx = -log(p): p0=2/3 -> 0.405465, p1=1/3 -> 1.098612; cand2 MISSING.
        assert self.d_tx[3, 0] == pytest.approx(-np.log(2 / 3))
        assert self.d_tx[3, 1] == pytest.approx(-np.log(1 / 3))
        assert self.d_tx[3, 2] == pytest.approx(MISSING)
        # softmax(-d_tx) over the tie reproduces the 2:1 ratio.
        w = np.exp(-(self.d_tx[3, :2] - self.d_tx[3, :2].min()))
        p = w / w.sum()
        assert p[0] == pytest.approx(2 / 3)
        assert p[1] == pytest.approx(1 / 3)


class TestDiffCoverGatePriorSource:
    """Only covered+distinguishing reads feed the redistribution prior."""

    def setup_method(self):
        self.captured = {}
        cfg = PipelineConfig(
            bam_path="/tmp/x.bam", quant_mode="m2_em", use_gpu=False,
            m2_diff_cover_gate=True, m2_diff_cover_margin=0.5,
        )
        _run(cfg, self.captured, tie_nll_fn=_fake_tie_nll_on2)
        self.d_tx = np.asarray(self.captured["dist_read_to_tx"])

    def test_covered_indistinguishable_flat_when_no_prior(self):
        # r0: covered but margin 0.05 < 0.5 -> fuzzy; the tie has no covered prior
        # (no covered+distinguishing read), so it falls back to a flat 1/K split.
        assert self.d_tx[0, 0] == pytest.approx(0.0)
        assert self.d_tx[0, 1] == pytest.approx(0.0)
        assert self.d_tx[0, 2] == pytest.approx(MISSING)

    def test_not_covered_distinguishing_still_hard(self):
        # r1: not covered but margin 4 >= 0.5 -> hard to candidate 0.
        assert self.d_tx[1, 0] == pytest.approx(0.0)
        assert self.d_tx[1, 1] == pytest.approx(MISSING)

    def test_fuzzy_falls_back_to_flat_when_no_covered_votes(self):
        # Neither r0 (indistinguishable) nor r1 (not covered) feeds covered_vote,
        # so r2/r3 fuzzy reads see an all-zero prior -> flat 1/K fallback.
        for i in (2, 3):
            assert self.d_tx[i, 0] == pytest.approx(0.0)
            assert self.d_tx[i, 1] == pytest.approx(0.0)
            assert self.d_tx[i, 2] == pytest.approx(MISSING)

    def test_no_reads_dropped(self):
        assert self.d_tx.shape == (N_R, N_C)


# Scenario 4 (covered+indistinguishable now FOLLOWS the covered prior, not flat):
#   r0: tie{0,1} nlls 1.0/5.0 -> margin 4 >= 0.5, COVERED -> hard to 0; vote[0]+=1
#   r1: tie{0,1} nlls 1.0/5.0 -> margin 4 >= 0.5, COVERED -> hard to 0; vote[0]+=1
#   r2: tie{0,1} nlls 2.0/2.05 -> margin 0.05 < 0.5, COVERED -> fuzzy; prior
#       [vote0=2, vote1=0] -> p=[1,0] -> hard-follows to candidate 0.
#   r3: tie{0,1} nlls 2.0/2.05 -> margin 0.05 < 0.5, NOT covered -> fuzzy; same prior.
ON4_TIES = {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1]}
ON4_NLLS = {0: {0: 1.0, 1: 5.0}, 1: {0: 1.0, 1: 5.0},
            2: {0: 2.0, 1: 2.05}, 3: {0: 2.0, 1: 2.05}}
ON4_COVER = {0: True, 1: True, 2: True, 3: False}


def _fake_tie_nll_on4(self, kept_read_ids, read_seqs, cand_list, aligners, raw):
    return dict(ON4_NLLS), dict(ON4_TIES), 4, len(ON4_NLLS), dict(ON4_COVER)


class TestDiffCoverGateCoveredIndistFollowsPrior:
    """covered+indistinguishable reads now follow the covered-read ratio (the
    margin is the SOLE decider of hard-assign vs ratio-follow)."""

    def setup_method(self):
        self.captured = {}
        cfg = PipelineConfig(
            bam_path="/tmp/x.bam", quant_mode="m2_em", use_gpu=False,
            m2_diff_cover_gate=True, m2_diff_cover_margin=0.5,
        )
        _run(cfg, self.captured, tie_nll_fn=_fake_tie_nll_on4)
        self.d_tx = np.asarray(self.captured["dist_read_to_tx"])

    def test_covered_indistinguishable_follows_prior(self):
        # r2 covered+indist -> prior [2,0] => p=[1,0]: -log(1)=0 on cand0, MISSING cand1.
        assert self.d_tx[2, 0] == pytest.approx(0.0)
        assert self.d_tx[2, 1] == pytest.approx(MISSING)

    def test_not_covered_indistinguishable_follows_same_prior(self):
        assert self.d_tx[3, 0] == pytest.approx(0.0)
        assert self.d_tx[3, 1] == pytest.approx(MISSING)


# Scenario 3 (fallback path: cover_by_read empty -> coverage uncomputable):
#   r0: tie{0,1} nlls 1/5 -> margin 4 distinguishes -> hard to 0, but coverage is
#       UNKNOWN (not in cover map) -> must NOT seed the prior.
#   r1..r3: tie{0,1} NO nlls -> fuzzy; no candidate earned prior votes -> flat.
ON3_TIES = {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1]}
ON3_NLLS = {0: {0: 1.0, 1: 5.0}}
ON3_COVER: dict = {}  # empty -> mimics the per-read _tie_nll fallback path


def _fake_tie_nll_on3(self, kept_read_ids, read_seqs, cand_list, aligners, raw):
    return dict(ON3_NLLS), dict(ON3_TIES), 4, len(ON3_NLLS), dict(ON3_COVER)


class TestDiffCoverGateFallbackCoverage:
    """When coverage is uncomputable (empty cover map, the per-read fallback),
    distinguishing reads still hard-assign but must NOT seed the prior, so fuzzy
    reads fall back to the flat 1/K split (recall-safe)."""

    def setup_method(self):
        self.captured = {}
        cfg = PipelineConfig(
            bam_path="/tmp/x.bam", quant_mode="m2_em", use_gpu=False,
            m2_diff_cover_gate=True, m2_diff_cover_margin=0.5,
        )
        _run(cfg, self.captured, tie_nll_fn=_fake_tie_nll_on3)
        self.d_tx = np.asarray(self.captured["dist_read_to_tx"])

    def test_distinguishing_read_still_hard_assigns(self):
        assert self.d_tx[0, 0] == pytest.approx(0.0)
        assert self.d_tx[0, 1] == pytest.approx(MISSING)

    def test_fuzzy_reads_flat_not_biased_by_unverified_coverage(self):
        # r0's coverage is unknown -> no prior -> r1/r2/r3 fuzzy -> flat 1/K.
        for i in (1, 2, 3):
            assert self.d_tx[i, 0] == pytest.approx(0.0)
            assert self.d_tx[i, 1] == pytest.approx(0.0)
            assert self.d_tx[i, 2] == pytest.approx(MISSING)


class TestTieNllKrillUnavailable:
    """When krill has no aligner, _tie_nll still exposes AS tie sets on the gate-ON
    path (so the gate's unique-best/drop logic is well-defined), but leaves them
    empty on the gate-OFF path (legacy byte-identical behavior)."""

    def _call(self, gate):
        cands = [_cand(f"c{j}") for j in range(N_C)]
        kept = ["r0", "r1"]
        read_seqs = {"r0": "ACGT" * 25, "r1": "ACGT" * 25}
        # r0 unique-best on candidate 1; r1 ties candidates {0, 2}.
        raw = np.array(
            [[-np.inf, 5.0, -np.inf], [3.0, -np.inf, 3.0]], dtype=np.float64
        )
        cfg = PipelineConfig(
            bam_path="/tmp/x.bam", quant_mode="m2_em", use_gpu=False,
            m2_diff_cover_gate=gate,
        )
        runner = PipelineRunner(cfg)
        with patch(
            "fin.scoring.krill_aligner.make_krill_aligner",
            return_value=(None, False),
        ):
            return runner._tie_nll(kept, read_seqs, cands, [None] * N_C, raw)

    def test_gate_on_populates_ties_from_raw(self):
        nlls, ties, n_ties, n_ref, cover = self._call(gate=True)
        assert ties == {0: [1], 1: [0, 2]}
        assert nlls == {} and cover == {}

    def test_gate_off_leaves_ties_empty(self):
        nlls, ties, n_ties, n_ref, cover = self._call(gate=False)
        assert ties == {}
        assert nlls == {} and cover == {}


class TestM3Gate:
    """m3_coherence=True turns on the DTW term at beta=em_beta."""

    def test_m3_on_calls_build_and_uses_em_beta(self):
        captured = {}
        # Non-zero DTW so the sigma3 normalization path runs.
        m3 = np.full((N_R, N_R), 2.0, dtype=np.float32)
        np.fill_diagonal(m3, 0.0)
        cfg = PipelineConfig(
            bam_path="/tmp/x.bam", quant_mode="m2_em", use_gpu=False,
            m3_coherence=True, em_beta=1.0, m2_diff_cover_gate=False,
        )
        mock_m3 = _run(cfg, captured, m3_matrix=m3)
        mock_m3.assert_called_once()
        assert captured["beta"] == 1.0
        d_rr = np.asarray(captured["dist_read_to_read"])
        assert d_rr.shape == (N_R, N_R)
        # sigma3 = median of non-zero DTW (all 2.0) -> normalized to 1.0.
        assert d_rr.max() == pytest.approx(1.0)
