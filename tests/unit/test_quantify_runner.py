"""Tests for QuantifyRunner config wiring (US-012)."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from fin.pipeline.config import PipelineConfig
from fin.pipeline.quantify_runner import QuantifyRunner, SampleInput


def _make_runner(config: PipelineConfig) -> QuantifyRunner:
    """Helper: build a QuantifyRunner with a dummy config and no real files."""
    return QuantifyRunner(
        gtf_path="/dev/null",
        genome_fasta_path="/dev/null",
        samples=[],
        output_dir="/tmp/pyfin_test",
        config=config,
    )


def _make_process_interval_mocks():
    """Return mock objects needed to exercise _process_interval without I/O."""
    interval = MagicMock()
    interval.region_string = "chr1:1-1000"
    interval.chrom = "chr1"

    sample = SampleInput(
        name="s1",
        bam_path="/dev/null",
        fastq_path="/dev/null",
        signal_path="/dev/null",
    )

    tool_runner = MagicMock()
    tool_runner.score_candidates.return_value = []

    signal_reader = MagicMock()
    sample_work_dir = MagicMock()

    return interval, sample, tool_runner, signal_reader, sample_work_dir


# ---------------------------------------------------------------------------
# AC: use_gpu forwarded to em_with_coherence
# ---------------------------------------------------------------------------

@patch("fin.pipeline.quantify_runner.em_with_coherence")
@patch("fin.pipeline.quantify_runner.discover_gtf_only")
@patch("fin.pipeline.quantify_runner.build_distance_matrix")
@patch("fin.pipeline.quantify_runner.parse_eventalign_tsv")
@patch("fin.pipeline.quantify_runner.extract_signal_segments")
@patch("fin.pipeline.quantify_runner.compute_read_to_read_dtw")
@patch("fin.pipeline.quantify_runner.quantify_transcripts")
def test_quantify_runner_forwards_use_gpu(
    mock_quant,
    mock_dtw,
    mock_segments,
    mock_parse,
    mock_build_dm,
    mock_discover,
    mock_em,
):
    """em_with_coherence must receive use_gpu from PipelineConfig."""
    n_reads, n_tx = 3, 2
    mock_em.return_value = (
        np.ones((n_reads, n_tx)) / n_tx,
        np.zeros(n_reads, dtype=int),
        [],
    )
    mock_quant.return_value = []
    mock_dtw.return_value = np.zeros((n_reads, n_reads))
    mock_build_dm.return_value = np.zeros((n_reads, n_tx))
    mock_segments.return_value = {}

    candidate_set = MagicMock()
    candidate_set.num_candidates = n_tx
    candidate_set.read_ids = ["r0", "r1", "r2"]
    candidate_set.candidate_ids.return_value = ["t0", "t1"]
    candidate_set.candidates = [MagicMock(), MagicMock()]
    mock_discover.return_value = candidate_set
    mock_parse.return_value = []

    for gpu_val in (True, False):
        mock_em.reset_mock()
        cfg = PipelineConfig(bam_path="/dev/null", use_gpu=gpu_val, use_prior=False)
        runner = _make_runner(cfg)
        runner._genome_fasta = {"chr1": "ACGT" * 250}

        interval, sample, tool_runner, signal_reader, sample_work_dir = (
            _make_process_interval_mocks()
        )
        sample_work_dir.__truediv__ = lambda self, other: MagicMock(
            __truediv__=lambda s, o: MagicMock(mkdir=MagicMock(), exists=MagicMock(return_value=False))
        )

        with patch("pathlib.Path.mkdir"):
            runner._process_interval(
                interval, sample, tool_runner, signal_reader, sample_work_dir
            )

        call_kwargs = mock_em.call_args.kwargs
        assert call_kwargs["use_gpu"] == gpu_val, (
            f"Expected use_gpu={gpu_val}, got {call_kwargs.get('use_gpu')}"
        )


# ---------------------------------------------------------------------------
# AC: score_alpha read from config
# ---------------------------------------------------------------------------

def test_quantify_runner_reads_score_alpha():
    """Config score_alpha should be accessible on the runner via self.config."""
    cfg = PipelineConfig(bam_path="/dev/null", score_alpha=0.75)
    runner = _make_runner(cfg)
    assert runner.config.score_alpha == 0.75


# ---------------------------------------------------------------------------
# AC: use_prior=False passes prior_weights=None
# ---------------------------------------------------------------------------

@patch("fin.pipeline.quantify_runner.em_with_coherence")
@patch("fin.pipeline.quantify_runner.discover_gtf_only")
@patch("fin.pipeline.quantify_runner.build_distance_matrix")
@patch("fin.pipeline.quantify_runner.parse_eventalign_tsv")
@patch("fin.pipeline.quantify_runner.extract_signal_segments")
@patch("fin.pipeline.quantify_runner.compute_read_to_read_dtw")
@patch("fin.pipeline.quantify_runner.quantify_transcripts")
def test_quantify_runner_use_prior_false_passes_none(
    mock_quant,
    mock_dtw,
    mock_segments,
    mock_parse,
    mock_build_dm,
    mock_discover,
    mock_em,
):
    """When use_prior=False, em_with_coherence must be called with prior_weights=None."""
    n_reads, n_tx = 2, 2
    mock_em.return_value = (
        np.ones((n_reads, n_tx)) / n_tx,
        np.zeros(n_reads, dtype=int),
        [],
    )
    mock_quant.return_value = []
    mock_dtw.return_value = np.zeros((n_reads, n_reads))
    mock_build_dm.return_value = np.zeros((n_reads, n_tx))
    mock_segments.return_value = {}

    candidate_set = MagicMock()
    candidate_set.num_candidates = n_tx
    candidate_set.read_ids = ["r0", "r1"]
    candidate_set.candidate_ids.return_value = ["t0", "t1"]
    candidate_set.candidates = [MagicMock(), MagicMock()]
    mock_discover.return_value = candidate_set
    mock_parse.return_value = []

    cfg = PipelineConfig(bam_path="/dev/null", use_prior=False)
    runner = _make_runner(cfg)
    runner._genome_fasta = {"chr1": "ACGT" * 250}

    interval, sample, tool_runner, signal_reader, sample_work_dir = (
        _make_process_interval_mocks()
    )

    with patch("pathlib.Path.mkdir"):
        runner._process_interval(
            interval, sample, tool_runner, signal_reader, sample_work_dir
        )

    call_kwargs = mock_em.call_args.kwargs
    assert call_kwargs["prior_weights"] is None


# ---------------------------------------------------------------------------
# AC: default config has use_prior=True
# ---------------------------------------------------------------------------

def test_quantify_runner_config_defaults():
    """Default PipelineConfig must have use_prior=True."""
    cfg = PipelineConfig(bam_path="/dev/null")
    assert cfg.use_prior is True
    runner = _make_runner(cfg)
    assert runner.config is cfg
    assert runner.config.use_prior is True


# ---------------------------------------------------------------------------
# AC: backward compat — use_prior=False gives same em call shape as pre-change
# ---------------------------------------------------------------------------

@patch("fin.pipeline.quantify_runner.em_with_coherence")
@patch("fin.pipeline.quantify_runner.discover_gtf_only")
@patch("fin.pipeline.quantify_runner.build_distance_matrix")
@patch("fin.pipeline.quantify_runner.parse_eventalign_tsv")
@patch("fin.pipeline.quantify_runner.extract_signal_segments")
@patch("fin.pipeline.quantify_runner.compute_read_to_read_dtw")
@patch("fin.pipeline.quantify_runner.quantify_transcripts")
def test_quantify_runner_backward_compat_no_prior(
    mock_quant,
    mock_dtw,
    mock_segments,
    mock_parse,
    mock_build_dm,
    mock_discover,
    mock_em,
):
    """With use_prior=False, em_with_coherence is called with prior_weights=None,
    which is bit-identical to the pre-change behavior (verified by US-009 tests)."""
    n_reads, n_tx = 4, 3
    mock_em.return_value = (
        np.ones((n_reads, n_tx)) / n_tx,
        np.zeros(n_reads, dtype=int),
        [],
    )
    mock_quant.return_value = []
    mock_dtw.return_value = np.zeros((n_reads, n_reads))
    mock_build_dm.return_value = np.zeros((n_reads, n_tx))
    mock_segments.return_value = {}

    candidate_set = MagicMock()
    candidate_set.num_candidates = n_tx
    candidate_set.read_ids = ["r0", "r1", "r2", "r3"]
    candidate_set.candidate_ids.return_value = ["t0", "t1", "t2"]
    candidate_set.candidates = [MagicMock(), MagicMock(), MagicMock()]
    mock_discover.return_value = candidate_set
    mock_parse.return_value = []

    cfg = PipelineConfig(bam_path="/dev/null", use_prior=False, use_gpu=False)
    runner = _make_runner(cfg)
    runner._genome_fasta = {"chr1": "ACGT" * 250}

    interval, sample, tool_runner, signal_reader, sample_work_dir = (
        _make_process_interval_mocks()
    )

    with patch("pathlib.Path.mkdir"):
        runner._process_interval(
            interval, sample, tool_runner, signal_reader, sample_work_dir
        )

    call_kwargs = mock_em.call_args.kwargs
    # prior_weights=None is equivalent to pre-change (no prior_weights kwarg)
    assert call_kwargs["prior_weights"] is None
    assert call_kwargs["use_gpu"] is False


# ---------------------------------------------------------------------------
# AC: use_prior=True derives a valid probability vector and forwards to EM
# ---------------------------------------------------------------------------

@patch("fin.pipeline.quantify_runner.em_with_coherence")
@patch("fin.pipeline.quantify_runner.discover_gtf_only")
@patch("fin.pipeline.quantify_runner.build_distance_matrix")
@patch("fin.pipeline.quantify_runner.parse_eventalign_tsv")
@patch("fin.pipeline.quantify_runner.extract_signal_segments")
@patch("fin.pipeline.quantify_runner.compute_read_to_read_dtw")
@patch("fin.pipeline.quantify_runner.quantify_transcripts")
def test_quantify_runner_use_prior_true_yields_valid_weights(
    mock_quant,
    mock_dtw,
    mock_segments,
    mock_parse,
    mock_build_dm,
    mock_discover,
    mock_em,
):
    """When use_prior=True, em_with_coherence receives a valid probability vector."""
    n_reads, n_tx = 4, 3
    mock_em.return_value = (
        np.ones((n_reads, n_tx)) / n_tx,
        np.zeros(n_reads, dtype=int),
        [],
    )
    mock_quant.return_value = []
    # Distinct distances to make composite scores non-degenerate.
    rng = np.random.default_rng(0)
    mock_dtw.return_value = np.abs(rng.standard_normal((n_reads, n_reads)))
    mock_build_dm.return_value = np.abs(rng.standard_normal((n_reads, n_tx)))
    mock_segments.return_value = {}

    candidate_set = MagicMock()
    candidate_set.num_candidates = n_tx
    candidate_set.read_ids = ["r0", "r1", "r2", "r3"]
    candidate_set.candidate_ids.return_value = ["t0", "t1", "t2"]
    candidate_set.candidates = [MagicMock(candidate_id=f"t{i}") for i in range(n_tx)]
    mock_discover.return_value = candidate_set
    mock_parse.return_value = []

    cfg = PipelineConfig(
        bam_path="/dev/null",
        use_prior=True,
        prior_weight_cap=10.0,
        use_gpu=False,
    )
    runner = _make_runner(cfg)
    runner._genome_fasta = {"chr1": "ACGT" * 250}

    interval, sample, tool_runner, signal_reader, sample_work_dir = (
        _make_process_interval_mocks()
    )

    with patch("pathlib.Path.mkdir"):
        runner._process_interval(
            interval, sample, tool_runner, signal_reader, sample_work_dir
        )

    pw = mock_em.call_args.kwargs["prior_weights"]
    assert pw is not None
    assert pw.shape == (n_tx,)
    assert np.all(pw > 0)
    assert pw.sum() == pytest.approx(1.0)
    # Cap=10 enforces an upper bound of uniform*cap on each weight.
    uniform = 1.0 / n_tx
    assert np.all(pw <= uniform * 10.0 + 1e-9)


# ---------------------------------------------------------------------------
# AC: composite scores are written back onto each QuantResult
# ---------------------------------------------------------------------------

@patch("fin.pipeline.quantify_runner.em_with_coherence")
@patch("fin.pipeline.quantify_runner.discover_gtf_only")
@patch("fin.pipeline.quantify_runner.build_distance_matrix")
@patch("fin.pipeline.quantify_runner.parse_eventalign_tsv")
@patch("fin.pipeline.quantify_runner.extract_signal_segments")
@patch("fin.pipeline.quantify_runner.compute_read_to_read_dtw")
@patch("fin.pipeline.quantify_runner.quantify_transcripts")
def test_quantify_runner_populates_quant_score_fields(
    mock_quant,
    mock_dtw,
    mock_segments,
    mock_parse,
    mock_build_dm,
    mock_discover,
    mock_em,
):
    """coherence_score / discrimination_score / combined_score must be set on QuantResult."""
    from fin.analysis.quantification import QuantResult

    n_reads, n_tx = 3, 2
    mock_em.return_value = (
        np.ones((n_reads, n_tx)) / n_tx,
        np.zeros(n_reads, dtype=int),
        [],
    )
    qr0 = QuantResult(
        candidate_id="t0", abundance=1.0, confidence=0.9,
        num_assigned_reads=2, source="gtf",
    )
    qr1 = QuantResult(
        candidate_id="t1", abundance=0.0, confidence=0.0,
        num_assigned_reads=0, source="gtf",
    )
    mock_quant.return_value = [qr0, qr1]
    rng = np.random.default_rng(1)
    mock_dtw.return_value = np.abs(rng.standard_normal((n_reads, n_reads)))
    mock_build_dm.return_value = np.abs(rng.standard_normal((n_reads, n_tx)))
    mock_segments.return_value = {}

    candidate_set = MagicMock()
    candidate_set.num_candidates = n_tx
    candidate_set.read_ids = ["r0", "r1", "r2"]
    candidate_set.candidate_ids.return_value = ["t0", "t1"]
    candidate_set.candidates = [
        MagicMock(candidate_id="t0"),
        MagicMock(candidate_id="t1"),
    ]
    mock_discover.return_value = candidate_set
    mock_parse.return_value = []

    cfg = PipelineConfig(bam_path="/dev/null", use_prior=False, use_gpu=False)
    runner = _make_runner(cfg)
    runner._genome_fasta = {"chr1": "ACGT" * 250}

    interval, sample, tool_runner, signal_reader, sample_work_dir = (
        _make_process_interval_mocks()
    )

    with patch("pathlib.Path.mkdir"):
        result = runner._process_interval(
            interval, sample, tool_runner, signal_reader, sample_work_dir
        )

    assert result is not None
    quant, _assignment = result
    assert quant is mock_quant.return_value
    # populate_quant_scores must have stamped the score fields on both QRs.
    for qr in quant:
        assert 0.0 <= qr.coherence_score <= 1.0
        assert 0.0 <= qr.discrimination_score <= 1.0
        assert 0.0 <= qr.combined_score <= 1.0


# ---------------------------------------------------------------------------
# AC: DTW subsampling honors max_reads_per_interval_for_dtw
# ---------------------------------------------------------------------------

@patch("fin.pipeline.quantify_runner.em_with_coherence")
@patch("fin.pipeline.quantify_runner.discover_gtf_only")
@patch("fin.pipeline.quantify_runner.build_distance_matrix")
@patch("fin.pipeline.quantify_runner.parse_eventalign_tsv")
@patch("fin.pipeline.quantify_runner.extract_signal_segments")
@patch("fin.pipeline.quantify_runner.compute_read_to_read_dtw")
@patch("fin.pipeline.quantify_runner.quantify_transcripts")
def test_quantify_runner_dtw_subsampling(
    mock_quant,
    mock_dtw,
    mock_segments,
    mock_parse,
    mock_build_dm,
    mock_discover,
    mock_em,
):
    """When read count exceeds max_reads_per_interval_for_dtw, reads are subsampled."""
    n_reads_total, n_tx = 50, 2
    cap = 10
    mock_em.return_value = (
        np.ones((cap, n_tx)) / n_tx,
        np.zeros(cap, dtype=int),
        [],
    )
    mock_quant.return_value = []
    mock_dtw.return_value = np.zeros((cap, cap))
    mock_build_dm.return_value = np.zeros((cap, n_tx))
    mock_segments.return_value = {}

    candidate_set = MagicMock()
    candidate_set.num_candidates = n_tx
    candidate_set.read_ids = [f"r{i:03d}" for i in range(n_reads_total)]
    candidate_set.candidate_ids.return_value = ["t0", "t1"]
    candidate_set.candidates = [MagicMock(candidate_id=f"t{i}") for i in range(n_tx)]
    mock_discover.return_value = candidate_set
    mock_parse.return_value = []

    cfg = PipelineConfig(
        bam_path="/dev/null",
        use_prior=False,
        use_gpu=False,
        max_reads_per_interval_for_dtw=cap,
    )
    runner = _make_runner(cfg)
    runner._genome_fasta = {"chr1": "ACGT" * 250}

    interval, sample, tool_runner, signal_reader, sample_work_dir = (
        _make_process_interval_mocks()
    )

    with patch("pathlib.Path.mkdir"):
        runner._process_interval(
            interval, sample, tool_runner, signal_reader, sample_work_dir
        )

    # build_distance_matrix and compute_read_to_read_dtw must have been called
    # with the subsampled read list, not all 50 reads.
    bdm_args = mock_build_dm.call_args
    sub_reads = bdm_args.args[1] if len(bdm_args.args) > 1 else bdm_args.kwargs.get("read_ids")
    assert len(sub_reads) == cap

    dtw_args = mock_dtw.call_args
    sub_reads_dtw = dtw_args.args[1] if len(dtw_args.args) > 1 else dtw_args.kwargs.get("read_ids")
    assert len(sub_reads_dtw) == cap


# ---------------------------------------------------------------------------
# AC: PipelineConfig field values override constructor params
# ---------------------------------------------------------------------------

def test_quantify_runner_config_overrides_ctor_params():
    """When config is provided, its values must win over constructor params."""
    cfg = PipelineConfig(
        bam_path="/dev/null",
        signal_format="pod5",
        use_gpu=False,
        em_sigma=2.5,
        em_beta=0.7,
        em_max_iter=42,
        em_tol=1e-3,
        f5c_path="/usr/local/bin/f5c",
    )
    runner = QuantifyRunner(
        gtf_path="/dev/null",
        genome_fasta_path="/dev/null",
        samples=[],
        output_dir="/tmp/pyfin_test",
        # Constructor values that should be overridden by config:
        signal_format="slow5",
        use_gpu=True,
        em_sigma=1.0,
        em_beta=0.5,
        em_max_iter=1000,
        em_tol=1e-4,
        f5c_path="f5c",
        config=cfg,
    )
    assert runner.signal_format == "pod5"
    assert runner.use_gpu is False
    assert runner.em_sigma == 2.5
    assert runner.em_beta == 0.7
    assert runner.em_max_iter == 42
    assert runner.em_tol == 1e-3
    assert runner.f5c_path == "/usr/local/bin/f5c"


def test_quantify_runner_no_config_uses_ctor_params():
    """When config is None, constructor parameters drive runner behavior."""
    runner = QuantifyRunner(
        gtf_path="/dev/null",
        genome_fasta_path="/dev/null",
        samples=[],
        output_dir="/tmp/pyfin_test",
        signal_format="pod5",
        use_gpu=False,
        em_sigma=3.3,
        config=None,
    )
    assert runner.signal_format == "pod5"
    assert runner.use_gpu is False
    assert runner.em_sigma == 3.3
