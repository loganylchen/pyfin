"""Integration tests for PipelineRunner.process_interval() (US-011).

Exercises the six-phase pipeline with heavy mocking — no real BAM/f5c/CUDA.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from fin.candidates.dataclasses import (
    CandidateSet,
    IntronChain,
    TranscriptCandidate,
)
from fin.io.interval_manager import GenomicInterval
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _make_candidate(candidate_id: str, source: str = "gtf") -> TranscriptCandidate:
    return TranscriptCandidate(
        candidate_id=candidate_id,
        intron_chain=IntronChain(introns=()),
        three_prime_pos=100,
        sequence="ACGT" * 25,
        source=source,
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=0,
        end=100,
    )


def _make_interval() -> GenomicInterval:
    return GenomicInterval(chrom="chr1", start=0, end=1000, strand="+")


def _make_config(**overrides) -> PipelineConfig:
    cfg = PipelineConfig(
        bam_path="/dev/null",
        work_dir="/tmp/pyfin_integration_test",
        use_gpu=False,
        use_prior=False,
        fusion_enabled=False,
        em_max_iter=2,
        em_tol=1e-2,
    )
    for k, v in overrides.items():
        setattr(cfg, k, v)
    return cfg


def _make_runner(config: PipelineConfig, num_reads: int = 10) -> PipelineRunner:
    runner = PipelineRunner(config)

    # Stub tool_runner: score_candidates returns an empty TSV list — we'll
    # intercept parse_eventalign_tsv so it doesn't matter what paths are here.
    tool_runner = MagicMock()
    tool_runner.score_candidates.return_value = []
    runner._tool_runner = tool_runner

    # No real genome/gtf/signal needed.
    runner._gtf_reader = None
    runner._genome_fasta = {"chr1": "A" * 2000}
    runner._signal_reader = MagicMock()
    return runner


# ---------------------------------------------------------------------------
# Helpers to patch the pipeline phases
# ---------------------------------------------------------------------------


class _PhaseMocks:
    """Container for mocks used during process_interval patching."""

    def __init__(self, em_mock, composite_mock, discover_mock,
                 parse_eventalign_mock, build_dist_mock, extract_seg_mock,
                 dtw_mock, quantify_mock):
        self.em = em_mock
        self.composite = composite_mock
        self.discover = discover_mock
        self.parse_eventalign = parse_eventalign_mock
        self.build_dist = build_dist_mock
        self.extract_seg = extract_seg_mock
        self.dtw = dtw_mock
        self.quantify = quantify_mock


def _patch_phases(num_reads: int = 10, num_cands: int = 2, candidates=None):
    """Return a context manager stack that patches every external phase."""
    if candidates is None:
        candidates = [_make_candidate(f"tx{i}") for i in range(num_cands)]

    read_ids = [f"read_{i}" for i in range(num_reads)]
    candidate_set = CandidateSet(
        interval=_make_interval(),
        candidates=candidates,
        read_ids=set(read_ids),
    )

    dist_tx = np.random.RandomState(0).rand(num_reads, num_cands)
    dist_rr = np.zeros((num_reads, num_reads))
    # Make it symmetric, non-negative.
    dist_rr = dist_rr + 0.1
    np.fill_diagonal(dist_rr, 0.0)

    R = np.full((num_reads, num_cands), 1.0 / num_cands)
    hard = np.zeros(num_reads, dtype=int)

    from fin.analysis.quantification import QuantResult
    from fin.scoring.composite import CompositeScore

    quant_results = [
        QuantResult(
            candidate_id=c.candidate_id,
            abundance=float(num_reads) / num_cands,
            confidence=1.0 / num_cands,
            num_assigned_reads=num_reads // num_cands,
            source=c.source,
            chrom=c.chrom,
            strand=c.strand,
            start=c.start,
            end=c.end,
            exons=((c.start, c.end),),
        )
        for c in candidates
    ]
    composite_scores = [
        CompositeScore(
            candidate_id=c.candidate_id,
            coherence=0.8,
            discrimination=0.6,
            combined=0.7,
        )
        for c in candidates
    ]

    patches = {
        "discover_candidates": patch(
            "fin.pipeline.runner.discover_candidates",
            return_value=candidate_set,
        ),
        "parse_eventalign_tsv": patch(
            "fin.pipeline.runner.parse_eventalign_tsv",
            return_value=[],
        ),
        "build_distance_matrix": patch(
            "fin.pipeline.runner.build_distance_matrix",
            return_value=dist_tx,
        ),
        "extract_signal_segments": patch(
            "fin.pipeline.runner.extract_signal_segments",
            return_value={},
        ),
        "compute_read_to_read_dtw": patch(
            "fin.pipeline.runner.compute_read_to_read_dtw",
            return_value=dist_rr,
        ),
        "score_candidates_composite": patch(
            "fin.pipeline.runner.score_candidates_composite",
            return_value=composite_scores,
        ),
        "em_with_coherence": patch(
            "fin.pipeline.runner.em_with_coherence",
            return_value=(R, hard, [0.0]),
        ),
        "quantify_transcripts": patch(
            "fin.pipeline.runner.quantify_transcripts",
            return_value=quant_results,
        ),
    }
    return patches, candidate_set, quant_results, composite_scores


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_process_interval_populates_scores(tmp_path):
    """AC-14: QuantResult objects have score fields populated from Phase 4."""
    cfg = _make_config(work_dir=str(tmp_path))
    runner = _make_runner(cfg)

    patches, candidate_set, quant_results, composite_scores = _patch_phases(
        num_reads=10, num_cands=2
    )

    with patches["discover_candidates"], patches["parse_eventalign_tsv"], \
         patches["build_distance_matrix"], patches["extract_signal_segments"], \
         patches["compute_read_to_read_dtw"], patches["score_candidates_composite"], \
         patches["em_with_coherence"], patches["quantify_transcripts"]:
        results = runner.process_interval(_make_interval())

    assert results is not None
    assert len(results) == 2
    for qr in results:
        assert qr.coherence_score == pytest.approx(0.8)
        assert qr.discrimination_score == pytest.approx(0.6)
        assert qr.combined_score == pytest.approx(0.7)


def test_process_interval_fusion_enabled(tmp_path):
    """AC-16: fusion_enabled=True merges fusion candidates into the set."""
    cfg = _make_config(work_dir=str(tmp_path), fusion_enabled=True)
    runner = _make_runner(cfg)

    fusion_cand = _make_candidate("fusion_abc", source="fusion")
    gtf_cand = _make_candidate("tx0", source="gtf")

    # Discovery returns only a GTF candidate; fusion phase adds the fusion.
    read_ids = {f"read_{i}" for i in range(5)}
    base_set = CandidateSet(
        interval=_make_interval(),
        candidates=[gtf_cand],
        read_ids=read_ids,
    )

    from fin.fusion.breakpoints import Breakpoint

    fake_bp = Breakpoint(
        chromA="chr1", posA=100, strandA="+",
        chromB="chr2", posB=200, strandB="+",
        support_count=3,
        supporting_read_ids={"read_0", "read_1", "read_2"},
    )

    patches, _, quant_results, composite_scores = _patch_phases(
        num_reads=5, num_cands=2,
        candidates=[gtf_cand, fusion_cand],
    )
    # Override discover to return just the GTF; fusion merges in.
    patches["discover_candidates"] = patch(
        "fin.pipeline.runner.discover_candidates",
        return_value=base_set,
    )

    with patches["discover_candidates"], patches["parse_eventalign_tsv"], \
         patches["build_distance_matrix"], patches["extract_signal_segments"], \
         patches["compute_read_to_read_dtw"], patches["score_candidates_composite"], \
         patches["em_with_coherence"], patches["quantify_transcripts"], \
         patch("fin.fusion.parse_sa_tags", return_value=[fake_bp]) as psa, \
         patch("fin.fusion.cluster_breakpoints", return_value=[fake_bp]), \
         patch("fin.fusion.build_fusion_candidates", return_value=[fusion_cand]):
        results = runner.process_interval(_make_interval())

    psa.assert_called_once()
    # The fusion candidate_id should appear among the returned results.
    ids = {qr.candidate_id for qr in results}
    assert "fusion_abc" in ids


def test_process_interval_fusion_disabled_by_default(tmp_path):
    """fusion_enabled=False -> parse_sa_tags NOT called."""
    cfg = _make_config(work_dir=str(tmp_path), fusion_enabled=False)
    runner = _make_runner(cfg)

    patches, _, _, _ = _patch_phases(num_reads=5, num_cands=2)

    with patches["discover_candidates"], patches["parse_eventalign_tsv"], \
         patches["build_distance_matrix"], patches["extract_signal_segments"], \
         patches["compute_read_to_read_dtw"], patches["score_candidates_composite"], \
         patches["em_with_coherence"], patches["quantify_transcripts"], \
         patch("fin.fusion.parse_sa_tags") as psa:
        runner.process_interval(_make_interval())

    psa.assert_not_called()


def test_process_interval_use_gpu_threaded(tmp_path):
    """Both EM and composite scoring receive the config.use_gpu flag."""
    cfg = _make_config(work_dir=str(tmp_path), use_gpu=False)
    runner = _make_runner(cfg)

    patches, _, _, _ = _patch_phases(num_reads=6, num_cands=2)

    with patches["discover_candidates"], patches["parse_eventalign_tsv"], \
         patches["build_distance_matrix"], patches["extract_signal_segments"], \
         patches["compute_read_to_read_dtw"] as dtw_m, \
         patches["score_candidates_composite"] as comp_m, \
         patches["em_with_coherence"] as em_m, \
         patches["quantify_transcripts"]:
        runner.process_interval(_make_interval())

    _, em_kwargs = em_m.call_args
    assert em_kwargs.get("use_gpu") is False

    _, comp_kwargs = comp_m.call_args
    assert comp_kwargs.get("use_gpu") is False

    _, dtw_kwargs = dtw_m.call_args
    assert dtw_kwargs.get("use_gpu") is False


def test_process_interval_prior_passed(tmp_path):
    """use_prior=True -> EM gets a non-None prior_weights summing to 1."""
    cfg = _make_config(work_dir=str(tmp_path), use_prior=True)
    runner = _make_runner(cfg)

    patches, _, _, _ = _patch_phases(num_reads=6, num_cands=3)

    with patches["discover_candidates"], patches["parse_eventalign_tsv"], \
         patches["build_distance_matrix"], patches["extract_signal_segments"], \
         patches["compute_read_to_read_dtw"], \
         patches["score_candidates_composite"], \
         patches["em_with_coherence"] as em_m, \
         patches["quantify_transcripts"]:
        runner.process_interval(_make_interval())

    _, kwargs = em_m.call_args
    prior = kwargs.get("prior_weights")
    assert prior is not None
    assert prior.shape == (3,)
    assert prior.sum() == pytest.approx(1.0, abs=1e-6)


def test_process_interval_prior_disabled(tmp_path):
    """use_prior=False -> EM gets prior_weights=None."""
    cfg = _make_config(work_dir=str(tmp_path), use_prior=False)
    runner = _make_runner(cfg)

    patches, _, _, _ = _patch_phases(num_reads=6, num_cands=3)

    with patches["discover_candidates"], patches["parse_eventalign_tsv"], \
         patches["build_distance_matrix"], patches["extract_signal_segments"], \
         patches["compute_read_to_read_dtw"], \
         patches["score_candidates_composite"], \
         patches["em_with_coherence"] as em_m, \
         patches["quantify_transcripts"]:
        runner.process_interval(_make_interval())

    _, kwargs = em_m.call_args
    assert kwargs.get("prior_weights") is None


def test_process_interval_read_ids_aligned_assertion(tmp_path):
    """Mismatched DTW / read-tx matrix row counts -> AssertionError."""
    cfg = _make_config(work_dir=str(tmp_path))
    runner = _make_runner(cfg)

    patches, _, _, _ = _patch_phases(num_reads=5, num_cands=2)
    # Break the shape contract: dist_read_to_read has wrong row count.
    patches["compute_read_to_read_dtw"] = patch(
        "fin.pipeline.runner.compute_read_to_read_dtw",
        return_value=np.zeros((3, 3)),  # mismatched with dist_read_to_tx (5x2)
    )

    with patches["discover_candidates"], patches["parse_eventalign_tsv"], \
         patches["build_distance_matrix"], patches["extract_signal_segments"], \
         patches["compute_read_to_read_dtw"], \
         patches["score_candidates_composite"], \
         patches["em_with_coherence"], patches["quantify_transcripts"]:
        with pytest.raises(AssertionError):
            runner.process_interval(_make_interval())


def test_process_interval_dtw_subsampling_triggered(tmp_path):
    """With max_reads_per_interval_for_dtw=5 and 20 reads, DTW gets <=5 reads."""
    cfg = _make_config(
        work_dir=str(tmp_path),
        max_reads_per_interval_for_dtw=5,
    )
    runner = _make_runner(cfg)

    candidates = [_make_candidate(f"tx{i}") for i in range(2)]
    read_ids = {f"read_{i}" for i in range(20)}
    cs = CandidateSet(
        interval=_make_interval(),
        candidates=candidates,
        read_ids=read_ids,
    )

    # Distance matrices sized to subsampled read count = 5.
    dist_tx = np.zeros((5, 2)) + 0.1
    dist_rr = np.zeros((5, 5)) + 0.1
    np.fill_diagonal(dist_rr, 0.0)
    R = np.full((5, 2), 0.5)
    hard = np.zeros(5, dtype=int)

    from fin.analysis.quantification import QuantResult
    from fin.scoring.composite import CompositeScore

    qr = [
        QuantResult(
            candidate_id=c.candidate_id, abundance=2.5, confidence=0.5,
            num_assigned_reads=2, source=c.source, chrom=c.chrom,
            strand=c.strand, start=c.start, end=c.end, exons=((0, 100),),
        )
        for c in candidates
    ]
    css = [
        CompositeScore(candidate_id=c.candidate_id, coherence=0.5,
                       discrimination=0.5, combined=0.5)
        for c in candidates
    ]

    with patch("fin.pipeline.runner.discover_candidates", return_value=cs), \
         patch("fin.pipeline.runner.parse_eventalign_tsv", return_value=[]), \
         patch("fin.pipeline.runner.build_distance_matrix",
               return_value=dist_tx) as bd_m, \
         patch("fin.pipeline.runner.extract_signal_segments", return_value={}), \
         patch("fin.pipeline.runner.compute_read_to_read_dtw",
               return_value=dist_rr) as dtw_m, \
         patch("fin.pipeline.runner.score_candidates_composite", return_value=css), \
         patch("fin.pipeline.runner.em_with_coherence",
               return_value=(R, hard, [0.0])), \
         patch("fin.pipeline.runner.quantify_transcripts", return_value=qr):
        runner.process_interval(_make_interval())

    # build_distance_matrix and DTW both receive the subsampled read list.
    bd_args, _ = bd_m.call_args
    dtw_args, dtw_kwargs = dtw_m.call_args
    # build_distance_matrix(scores, read_ids, candidate_ids)
    assert len(bd_args[1]) <= 5
    # compute_read_to_read_dtw(segments, read_ids, use_gpu=...)
    assert len(dtw_args[1]) <= 5
