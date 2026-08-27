"""Integration tests for PipelineRunner.process_interval() (US-011).

Exercises the krill-only quant_mode='m2_em' EM path with heavy mocking — no
real BAM/krill/CUDA. The pure tie-break junction-NLL seed (``_tie_nll``), the
per-candidate mappy aligners and EM are stubbed; the test asserts the
orchestration contract, including dispatch, fusion augmentation, GPU threading,
and the production zero-coherence EM call.
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


class _FakeHit:
    """Minimal mappy hit: a clean 100-match alignment so the real score_hit
    returns a positive AS (2*100), marking the read as kept in the raw pass."""

    cigar = [(100, 0)]  # 100 matches, no indels
    NM = 0
    mlen = 100


class _FakeAligner:
    """mappy.Aligner stub: yields one clean hit per read so the raw AS pass marks
    every read as kept. The tie sets still come from the mocked ``_tie_nll``, so
    the raw matrix values beyond ">0" are unused."""

    def __init__(self, *a, **k):
        pass

    def map(self, *a, **k):
        return iter((_FakeHit(),))


def _fake_tie_nll(self, kept_read_ids, read_seqs, cand_list, aligners, raw):
    """Give every read a unique-best tie on candidate 0 (no scored NLLs)."""
    ties = {i: [0] for i in range(len(kept_read_ids))}
    return {}, ties, 0, 0, {}


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
        read_sequences={rid: "ACGT" * 25 for rid in read_ids},
    )

    R_mm = np.full((num_reads, num_cands), 1.0 / num_cands)
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
        "mappy_aligner": patch("mappy.Aligner", _FakeAligner),
        "tie_nll": patch.object(PipelineRunner, "_tie_nll", _fake_tie_nll),
        "eff_lengths": patch.object(
            PipelineRunner, "_eff_lengths", return_value=None
        ),
        "em_with_coherence": patch(
            "fin.pipeline.assignment.em_with_coherence",
            return_value=(R, hard, []),
        ),
        "quantify_transcripts": patch(
            "fin.pipeline.assignment.quantify_transcripts",
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
         patches["mappy_aligner"], patches["tie_nll"], patches["eff_lengths"], \
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
         patches["mappy_aligner"], patches["tie_nll"], patches["eff_lengths"], \
         patches["em_with_coherence"] as em_m, patches["quantify_transcripts"]:
        runner.process_interval(_make_interval())

    _, em_kwargs = em_m.call_args
    assert em_kwargs.get("use_gpu") is False


def test_process_interval_disables_read_coherence(tmp_path):
    """Production M2 calls the generic EM engine with zero coherence."""
    cfg = _make_config(work_dir=str(tmp_path))
    runner = _make_runner(cfg)

    patches, _, _ = _patch_m2_phases(num_reads=6, num_cands=2)

    with patches["discover_candidates"], patches["mappy_multimap_responsibilities"], \
         patches["mappy_aligner"], patches["tie_nll"], patches["eff_lengths"], \
         patches["em_with_coherence"] as em_m, patches["quantify_transcripts"]:
        runner.process_interval(_make_interval())

    _, em_kwargs = em_m.call_args
    assert em_kwargs.get("beta") == 0.0
    assert not np.asarray(em_kwargs["dist_read_to_read"]).any()


def _fake_bam_reader_with_reads(read_dicts):
    """Return a patch target factory: BamReader(...) as ctx -> .get_reads_in_region."""
    from unittest.mock import MagicMock

    reader = MagicMock()
    reader.get_reads_in_region.return_value = read_dicts
    ctx = MagicMock()
    ctx.__enter__.return_value = reader
    ctx.__exit__.return_value = False
    factory = MagicMock(return_value=ctx)
    return factory


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
        read_sequences={rid: "ACGT" * 25 for rid in read_ids},
    )

    patches, _, _ = _patch_m2_phases(
        num_reads=5, num_cands=2, candidates=[gtf_cand, fusion_cand],
    )
    # Override discover to return just the GTF; fusion merges in.
    patches["discover_candidates"] = patch(
        "fin.pipeline.runner.discover_candidates",
        return_value=base_set,
    )

    fake_bam = _fake_bam_reader_with_reads([{"query_name": "read_0"}])

    with patches["discover_candidates"], patches["mappy_multimap_responsibilities"], \
         patches["mappy_aligner"], patches["tie_nll"], patches["eff_lengths"], \
         patches["em_with_coherence"], patches["quantify_transcripts"], \
         patch.object(runner, "_get_fusion_genome_aligner", return_value=object()), \
         patch("fin.io.io_bam.BamReader", fake_bam), \
         patch("fin.fusion.detect_fusion_candidates",
               return_value=[fusion_cand]) as det:
        results = runner.process_interval(_make_interval())

    det.assert_called_once()
    ids = {qr.candidate_id for qr in results}
    assert "fusion_abc" in ids


def test_process_interval_fusion_disabled_by_default(tmp_path):
    """fusion_enabled=False -> fusion detection NOT invoked."""
    cfg = _make_config(work_dir=str(tmp_path), fusion_enabled=False)
    runner = _make_runner(cfg)

    patches, _, _ = _patch_m2_phases(num_reads=5, num_cands=2)

    with patches["discover_candidates"], patches["mappy_multimap_responsibilities"], \
         patches["mappy_aligner"], patches["tie_nll"], patches["eff_lengths"], \
         patches["em_with_coherence"], patches["quantify_transcripts"], \
         patch("fin.fusion.detect_fusion_candidates") as det:
        runner.process_interval(_make_interval())

    det.assert_not_called()
