"""Tests for QuantifyRunner config wiring (US-012), krill-only signal path."""

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
    """Return the args _process_interval needs (krill path), without I/O."""
    interval = MagicMock()
    interval.region_string = "chr1:1-1000"
    interval.chrom = "chr1"

    sample = SampleInput(
        name="s1",
        bam_path="/dev/null",
        fastq_path="/dev/null",
        signal_path="/dev/null",
    )

    sample_work_dir = MagicMock()
    cached_reads: list = []

    return interval, sample, sample_work_dir, cached_reads


def _wire_krill_mocks(
    mock_discover,
    mock_mm,
    mock_m2,
    mock_m3,
    mock_em,
    mock_quant,
    n_reads: int,
    n_tx: int,
    *,
    seed: int = 0,
    quant_return=None,
):
    """Configure the krill-path mocks for a single interval.

    The shape-coupled mocks use side_effects so DTW subsampling (which changes
    the read axis between the mappy floor and the krill matrices) stays
    internally consistent.
    """
    read_ids = [f"r{i:03d}" for i in range(n_reads)]
    cand_ids = [f"t{i}" for i in range(n_tx)]

    candidate_set = MagicMock()
    candidate_set.num_candidates = n_tx
    candidate_set.read_ids = read_ids
    candidate_set.candidate_ids.return_value = cand_ids
    # sequence="" so the per-candidate alignment BAM loop is skipped.
    candidate_set.candidates = [
        MagicMock(candidate_id=cid, sequence="") for cid in cand_ids
    ]
    candidate_set.read_sequences = {rid: "ACGT" for rid in read_ids}
    mock_discover.return_value = candidate_set

    # mappy structural floor: keep all reads, uniform responsibilities.
    mock_mm.return_value = (np.ones((n_reads, n_tx)) / n_tx, list(read_ids))

    rng = np.random.default_rng(seed)

    def _m2(dtw_ids, *args, **kwargs):
        # Moderate distances in (0, 0.5): below the 0.999 no-data sentinel.
        return 0.1 + 0.3 * np.abs(rng.standard_normal((len(dtw_ids), n_tx)))

    mock_m2.side_effect = _m2
    mock_m3.side_effect = lambda dtw_ids, *a, **k: np.zeros(
        (len(dtw_ids), len(dtw_ids)), dtype=np.float32
    )

    def _em(**kwargs):
        d = kwargs["dist_read_to_tx"]
        nr, nt = d.shape
        return (np.ones((nr, nt)) / nt, np.zeros(nr, dtype=int), [])

    mock_em.side_effect = _em
    mock_quant.return_value = [] if quant_return is None else quant_return
    return candidate_set


# Common patch stack for the krill _process_interval path. Decorators apply
# bottom-up, so the wrapped test receives:
#   mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant, mock_align
def _krill_patches(fn):
    # Applied innermost-first: the first patch applied maps to the first arg.
    fn = patch("fin.pipeline.quantify_runner.discover_gtf_only")(fn)
    fn = patch("fin.pipeline.quantify_runner.mappy_multimap_responsibilities")(fn)
    fn = patch("fin.pipeline.quantify_runner._build_m2_krill")(fn)
    fn = patch("fin.pipeline.quantify_runner.build_m3_coherence")(fn)
    fn = patch("fin.pipeline.quantify_runner.em_with_coherence")(fn)
    fn = patch("fin.pipeline.quantify_runner.quantify_transcripts")(fn)
    fn = patch("fin.pipeline.quantify_runner.align_reads_with_mappy")(fn)
    return fn


# ---------------------------------------------------------------------------
# AC: use_gpu forwarded to em_with_coherence
# ---------------------------------------------------------------------------

@_krill_patches
def test_quantify_runner_forwards_use_gpu(
    mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant, mock_align
):
    """em_with_coherence must receive use_gpu from PipelineConfig."""
    for gpu_val in (True, False):
        mock_em.reset_mock()
        _wire_krill_mocks(
            mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant,
            n_reads=3, n_tx=2,
        )
        cfg = PipelineConfig(
            bam_path="/dev/null", use_gpu=gpu_val, use_prior=False,
            persist_R_matrix=False,
        )
        runner = _make_runner(cfg)
        runner._genome_fasta = {"chr1": "ACGT" * 250}

        interval, sample, sample_work_dir, cached_reads = (
            _make_process_interval_mocks()
        )
        runner._process_interval(interval, sample, sample_work_dir, cached_reads)

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

@_krill_patches
def test_quantify_runner_use_prior_false_passes_none(
    mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant, mock_align
):
    """When use_prior=False, em_with_coherence must get prior_weights=None."""
    _wire_krill_mocks(
        mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant,
        n_reads=2, n_tx=2,
    )
    cfg = PipelineConfig(
        bam_path="/dev/null", use_prior=False, persist_R_matrix=False,
    )
    runner = _make_runner(cfg)
    runner._genome_fasta = {"chr1": "ACGT" * 250}

    interval, sample, sample_work_dir, cached_reads = _make_process_interval_mocks()
    runner._process_interval(interval, sample, sample_work_dir, cached_reads)

    assert mock_em.call_args.kwargs["prior_weights"] is None


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

@_krill_patches
def test_quantify_runner_backward_compat_no_prior(
    mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant, mock_align
):
    """use_prior=False -> prior_weights=None and use_gpu forwarded."""
    _wire_krill_mocks(
        mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant,
        n_reads=4, n_tx=3,
    )
    cfg = PipelineConfig(
        bam_path="/dev/null", use_prior=False, use_gpu=False,
        persist_R_matrix=False,
    )
    runner = _make_runner(cfg)
    runner._genome_fasta = {"chr1": "ACGT" * 250}

    interval, sample, sample_work_dir, cached_reads = _make_process_interval_mocks()
    runner._process_interval(interval, sample, sample_work_dir, cached_reads)

    call_kwargs = mock_em.call_args.kwargs
    assert call_kwargs["prior_weights"] is None
    assert call_kwargs["use_gpu"] is False


# ---------------------------------------------------------------------------
# AC: use_prior=True derives a valid probability vector and forwards to EM
# ---------------------------------------------------------------------------

@_krill_patches
def test_quantify_runner_use_prior_true_yields_valid_weights(
    mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant, mock_align
):
    """When use_prior=True, em_with_coherence receives a valid probability vector."""
    n_tx = 3
    _wire_krill_mocks(
        mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant,
        n_reads=4, n_tx=n_tx, seed=1,
    )
    cfg = PipelineConfig(
        bam_path="/dev/null",
        use_prior=True,
        prior_weight_cap=10.0,
        use_gpu=False,
        persist_R_matrix=False,
    )
    runner = _make_runner(cfg)
    runner._genome_fasta = {"chr1": "ACGT" * 250}

    interval, sample, sample_work_dir, cached_reads = _make_process_interval_mocks()
    runner._process_interval(interval, sample, sample_work_dir, cached_reads)

    pw = mock_em.call_args.kwargs["prior_weights"]
    assert pw is not None
    assert pw.shape == (n_tx,)
    assert np.all(pw > 0)
    assert pw.sum() == pytest.approx(1.0)
    uniform = 1.0 / n_tx
    assert np.all(pw <= uniform * 10.0 + 1e-9)


# ---------------------------------------------------------------------------
# AC: composite scores are written back onto each QuantResult
# ---------------------------------------------------------------------------

@_krill_patches
def test_quantify_runner_populates_quant_score_fields(
    mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant, mock_align
):
    """coherence/discrimination/combined_score must be set on each QuantResult."""
    from fin.analysis.quantification import QuantResult

    qr0 = QuantResult(
        candidate_id="t0", abundance=1.0, confidence=0.9,
        num_assigned_reads=2, source="gtf",
    )
    qr1 = QuantResult(
        candidate_id="t1", abundance=0.0, confidence=0.0,
        num_assigned_reads=0, source="gtf",
    )
    _wire_krill_mocks(
        mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant,
        n_reads=3, n_tx=2, seed=1, quant_return=[qr0, qr1],
    )

    cfg = PipelineConfig(
        bam_path="/dev/null", use_prior=False, use_gpu=False,
        persist_R_matrix=False,
    )
    runner = _make_runner(cfg)
    runner._genome_fasta = {"chr1": "ACGT" * 250}

    interval, sample, sample_work_dir, cached_reads = _make_process_interval_mocks()
    result = runner._process_interval(
        interval, sample, sample_work_dir, cached_reads
    )

    assert result is not None
    quant, _assignment = result
    assert quant is mock_quant.return_value
    for qr in quant:
        assert 0.0 <= qr.coherence_score <= 1.0
        assert 0.0 <= qr.discrimination_score <= 1.0
        assert 0.0 <= qr.combined_score <= 1.0


# ---------------------------------------------------------------------------
# AC: DTW subsampling honors max_reads_per_interval_for_dtw
# ---------------------------------------------------------------------------

@_krill_patches
def test_quantify_runner_dtw_subsampling(
    mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant, mock_align
):
    """When read count exceeds the cap, the krill matrices use subsampled reads."""
    cap = 10
    _wire_krill_mocks(
        mock_discover, mock_mm, mock_m2, mock_m3, mock_em, mock_quant,
        n_reads=50, n_tx=2,
    )
    cfg = PipelineConfig(
        bam_path="/dev/null",
        use_prior=False,
        use_gpu=False,
        max_reads_per_interval_for_dtw=cap,
        persist_R_matrix=False,
    )
    runner = _make_runner(cfg)
    runner._genome_fasta = {"chr1": "ACGT" * 250}

    interval, sample, sample_work_dir, cached_reads = _make_process_interval_mocks()
    runner._process_interval(interval, sample, sample_work_dir, cached_reads)

    # _build_m2_krill and build_m3_coherence must receive the subsampled list.
    m2_reads = mock_m2.call_args.args[0]
    assert len(m2_reads) == cap
    m3_reads = mock_m3.call_args.args[0]
    assert len(m3_reads) == cap


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
        config=cfg,
    )
    assert runner.signal_format == "pod5"
    assert runner.use_gpu is False
    assert runner.em_sigma == 2.5
    assert runner.em_beta == 0.7
    assert runner.em_max_iter == 42
    assert runner.em_tol == 1e-3


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
