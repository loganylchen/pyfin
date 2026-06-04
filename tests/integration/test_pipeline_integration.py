"""Integration tests for PipelineRunner.process_interval() (US-011).

Exercises the krill-only quant_mode='m2_em' EM path with heavy mocking — no
real BAM/krill/CUDA. The signal matrices (M2 read×tx distance, M3 read×read
coherence) and EM are stubbed; the test asserts the orchestration contract
(dispatch, fusion augmentation, use_gpu threading).
"""

from __future__ import annotations

from unittest.mock import patch

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
        signal_path="/dev/null",
        work_dir="/tmp/pyfin_integration_test",
        use_gpu=False,
        use_prior=False,
        fusion_enabled=False,
        quant_mode="m2_em",
        m4_source="diff_region",
        em_max_iter=2,
        em_tol=1e-2,
        # Avoid the post-quant fulllen BAM fetch (no real BAM here).
        min_fulllen_fraction=0.0,
    )
    for k, v in overrides.items():
        setattr(cfg, k, v)
    return cfg


def _make_runner(config: PipelineConfig) -> PipelineRunner:
    runner = PipelineRunner(config)
    # No real genome/gtf/signal needed.
    runner._gtf_reader = None
    runner._genome_fasta = {"chr1": "A" * 2000}
    runner._signal_reader = None
    return runner


# ---------------------------------------------------------------------------
# Helpers to patch the krill m2_em phases
# ---------------------------------------------------------------------------


def _patch_m2_phases(num_reads: int = 10, num_cands: int = 2, candidates=None):
    """Return a dict of patchers for the m2_em path plus the stub quant results.

    The signal matrices and EM are mocked at their source modules because
    ``_quant_m2_em`` imports them function-locally.
    """
    if candidates is None:
        candidates = [_make_candidate(f"tx{i}") for i in range(num_cands)]
    num_cands = len(candidates)

    read_ids = [f"read_{i}" for i in range(num_reads)]
    candidate_set = CandidateSet(
        interval=_make_interval(),
        candidates=candidates,
        read_ids=set(read_ids),
    )

    R_mm = np.full((num_reads, num_cands), 1.0 / num_cands)
    # M2 read×tx distances in (0, 0.5): below the 0.999 no-data sentinel.
    dist_tx = 0.1 + 0.3 * np.random.RandomState(0).rand(num_reads, num_cands)
    dist_rr = np.zeros((num_reads, num_reads), dtype=np.float32)
    R = np.full((num_reads, num_cands), 1.0 / num_cands)
    hard = np.zeros(num_reads, dtype=int)

    from fin.analysis.quantification import QuantResult

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

    patches = {
        "discover_candidates": patch(
            "fin.pipeline.runner.discover_candidates",
            return_value=candidate_set,
        ),
        "mappy_multimap_responsibilities": patch(
            "fin.ablation.mappy_argmax.mappy_multimap_responsibilities",
            return_value=(R_mm, list(read_ids)),
        ),
        "_build_m2_krill": patch(
            "fin.scoring.krill_tiebreak._build_m2_krill",
            return_value=dist_tx,
        ),
        "build_m3_coherence": patch(
            "fin.scoring.m3_junction_coherence.build_m3_coherence",
            return_value=dist_rr,
        ),
        "em_with_coherence": patch(
            "fin.pipeline.runner.em_with_coherence",
            return_value=(R, hard, []),
        ),
        "quantify_transcripts": patch(
            "fin.pipeline.runner.quantify_transcripts",
            return_value=quant_results,
        ),
    }
    return patches, candidate_set, quant_results


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_process_interval_m2_em_smoke(tmp_path):
    """m2_em dispatch returns the quantify_transcripts results with max_R set."""
    cfg = _make_config(work_dir=str(tmp_path))
    runner = _make_runner(cfg)

    patches, _candidate_set, quant_results = _patch_m2_phases(
        num_reads=10, num_cands=2
    )

    with patches["discover_candidates"], patches["mappy_multimap_responsibilities"], \
         patches["_build_m2_krill"], patches["build_m3_coherence"], \
         patches["em_with_coherence"], patches["quantify_transcripts"]:
        results = runner.process_interval(_make_interval())

    assert results is not None
    assert len(results) == 2
    for qr in results:
        # max_R is stamped from the (uniform) EM responsibility matrix.
        assert qr.max_R == pytest.approx(0.5)


def test_process_interval_use_gpu_threaded(tmp_path):
    """em_with_coherence receives the config.use_gpu flag."""
    cfg = _make_config(work_dir=str(tmp_path), use_gpu=False)
    runner = _make_runner(cfg)

    patches, _, _ = _patch_m2_phases(num_reads=6, num_cands=2)

    with patches["discover_candidates"], patches["mappy_multimap_responsibilities"], \
         patches["_build_m2_krill"], patches["build_m3_coherence"], \
         patches["em_with_coherence"] as em_m, patches["quantify_transcripts"]:
        runner.process_interval(_make_interval())

    _, em_kwargs = em_m.call_args
    assert em_kwargs.get("use_gpu") is False


def test_process_interval_m4_none_disables_coherence(tmp_path):
    """m4_source='none' -> EM beta=0 and build_m3_coherence is not called."""
    cfg = _make_config(work_dir=str(tmp_path), m4_source="none")
    runner = _make_runner(cfg)

    patches, _, _ = _patch_m2_phases(num_reads=6, num_cands=2)

    with patches["discover_candidates"], patches["mappy_multimap_responsibilities"], \
         patches["_build_m2_krill"], patches["build_m3_coherence"] as m3_m, \
         patches["em_with_coherence"] as em_m, patches["quantify_transcripts"]:
        runner.process_interval(_make_interval())

    m3_m.assert_not_called()
    _, em_kwargs = em_m.call_args
    assert em_kwargs.get("beta") == 0.0


def test_process_interval_fusion_enabled(tmp_path):
    """AC-16: fusion_enabled=True merges fusion candidates into the set."""
    cfg = _make_config(work_dir=str(tmp_path), fusion_enabled=True)
    runner = _make_runner(cfg)

    fusion_cand = _make_candidate("fusion_abc", source="fusion")
    gtf_cand = _make_candidate("tx0", source="gtf")

    # Discovery returns only a GTF candidate; the fusion phase adds the fusion.
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

    patches, _, _ = _patch_m2_phases(
        num_reads=5, num_cands=2, candidates=[gtf_cand, fusion_cand],
    )
    # Override discover to return just the GTF; fusion merges in.
    patches["discover_candidates"] = patch(
        "fin.pipeline.runner.discover_candidates",
        return_value=base_set,
    )

    with patches["discover_candidates"], patches["mappy_multimap_responsibilities"], \
         patches["_build_m2_krill"], patches["build_m3_coherence"], \
         patches["em_with_coherence"], patches["quantify_transcripts"], \
         patch("fin.fusion.parse_sa_tags", return_value=[fake_bp]) as psa, \
         patch("fin.fusion.cluster_breakpoints", return_value=[fake_bp]), \
         patch("fin.fusion.build_fusion_candidates", return_value=[fusion_cand]):
        results = runner.process_interval(_make_interval())

    psa.assert_called_once()
    ids = {qr.candidate_id for qr in results}
    assert "fusion_abc" in ids


def test_process_interval_fusion_disabled_by_default(tmp_path):
    """fusion_enabled=False -> parse_sa_tags NOT called."""
    cfg = _make_config(work_dir=str(tmp_path), fusion_enabled=False)
    runner = _make_runner(cfg)

    patches, _, _ = _patch_m2_phases(num_reads=5, num_cands=2)

    with patches["discover_candidates"], patches["mappy_multimap_responsibilities"], \
         patches["_build_m2_krill"], patches["build_m3_coherence"], \
         patches["em_with_coherence"], patches["quantify_transcripts"], \
         patch("fin.fusion.parse_sa_tags") as psa:
        runner.process_interval(_make_interval())

    psa.assert_not_called()
