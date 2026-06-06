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
    """mappy.Aligner stub: builds, maps to nothing (raw stays -inf; _tie_nll is
    mocked so the raw AS matrix is unused for tie detection here)."""

    def __init__(self, *a, **k):
        pass

    def map(self, *a, **k):
        return iter(())


def _fake_tie_nll(self, kept_read_ids, read_seqs, cand_list, aligners, raw):
    n_ties = sum(1 for t in TIES.values() if len(t) >= 2)
    n_refined = len(NLLS)
    return dict(NLLS), dict(TIES), n_ties, n_refined


def _run(cfg, captured, m3_matrix=None):
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
        "fin.ablation.mappy_argmax.mappy_multimap_responsibilities",
        return_value=(np.ones((N_R, N_C)), list(READ_IDS)),
    ), patch.object(PipelineRunner, "_tie_nll", _fake_tie_nll), patch.object(
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


class TestM3Gate:
    """m3_coherence=True turns on the DTW term at beta=em_beta."""

    def test_m3_on_calls_build_and_uses_em_beta(self):
        captured = {}
        # Non-zero DTW so the sigma3 normalization path runs.
        m3 = np.full((N_R, N_R), 2.0, dtype=np.float32)
        np.fill_diagonal(m3, 0.0)
        cfg = PipelineConfig(
            bam_path="/tmp/x.bam", quant_mode="m2_em", use_gpu=False,
            m3_coherence=True, em_beta=1.0,
        )
        mock_m3 = _run(cfg, captured, m3_matrix=m3)
        mock_m3.assert_called_once()
        assert captured["beta"] == 1.0
        d_rr = np.asarray(captured["dist_read_to_read"])
        assert d_rr.shape == (N_R, N_R)
        # sigma3 = median of non-zero DTW (all 2.0) -> normalized to 1.0.
        assert d_rr.max() == pytest.approx(1.0)
