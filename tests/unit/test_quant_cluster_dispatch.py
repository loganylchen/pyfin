"""Tests for the within-cluster quant (``_quant_cluster``, quant_mode="cluster").

``_quant_cluster`` scores each generation cluster (``CandidateSet.clusters``) in
isolation: a read is aligned only against that cluster's members (M1), straddling
ties go to the summed-LLR M2 tiebreak, and containment/truncation ties fall to the
cluster's main peak. These tests run signal-free (``signal_path=""`` -> no krill ->
``m2_resolve`` is None), so every tie is deterministic. The BAM span fetch is
monkeypatched so the test needs no real BAM file.
"""

from __future__ import annotations

import random

import mappy
import pytest

from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval
from fin.pipeline import runner as runner_mod
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner


def _rnd(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


# A distinctive 5' prefix + a shared 3' body, each long/complex enough to map. The
# "full" member carries both; the "short" (truncation) member carries only the
# shared body. Verified mappy AS: a read of the shared body alone scores 300 on
# BOTH members (a containment tie), while a read of the full sequence scores 600 on
# the full member vs 300 on the short one (uniquely-best -> its independent anchor).
_RNG = random.Random(42)
PREFIX = _rnd(150, _RNG)
SHARED = _rnd(150, _RNG)
FULL = PREFIX + SHARED
OTHER = _rnd(300, random.Random(7))


def _cand(cid, seq, source="novel"):
    # Empty intron chain -> _differing_coords is empty -> no tie ever "straddles"
    # a differing junction -> every tie is a containment tie (main-peak routing),
    # which is exactly the deterministic path this test pins.
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=()),
        three_prime_pos=0,
        sequence=seq,
        source=source,
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=0,
        end=len(seq),
    )


def _interval():
    return GenomicInterval(chrom="chr1", start=0, end=200, strand="+")


class _FakeRead:
    def __init__(self, name, start, end):
        self.query_name = name
        self.reference_start = start
        self.reference_end = end
        self.is_unmapped = False
        self.is_secondary = False
        self.is_supplementary = False


class _FakeBam:
    def __init__(self, reads):
        self._reads = reads

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def fetch(self, *a, **k):
        return iter(self._reads)


def _patch_bam(monkeypatch, reads):
    monkeypatch.setattr(
        runner_mod.pysam, "AlignmentFile", lambda *a, **k: _FakeBam(reads)
    )


def _runner():
    # signal_path="" -> krill leg skipped -> m2_resolve stays None (fine: no
    # straddling ties occur here anyway).
    cfg = PipelineConfig(bam_path="x", quant_mode="cluster")
    return PipelineRunner(cfg)


@pytest.mark.skipif(mappy is None, reason="mappy required for cluster quant")
class TestQuantClusterDispatch:
    def test_returns_one_result_per_candidate(self, monkeypatch):
        cands = [_cand("c_full", FULL), _cand("c_short", SHARED)]
        # Cluster reads = union of members' supporting_read_ids (scoped by read_ids).
        cands[0].supporting_read_ids.update({"anchor", "shared"})
        cands[1].supporting_read_ids.add("shared")
        cs = CandidateSet(
            interval=_interval(),
            candidates=cands,
            read_ids={"anchor", "shared"},
            read_sequences={"anchor": FULL, "shared": SHARED},
            clusters=[["c_full", "c_short"]],
        )
        _patch_bam(monkeypatch, [
            _FakeRead("anchor", 0, len(FULL)),
            _FakeRead("shared", len(PREFIX), len(FULL)),
        ])
        runner = _runner()
        results = runner._quant_cluster(cs, ["anchor", "shared"], _interval())
        assert len(results) == 2
        assert {r.candidate_id for r in results} == {"c_full", "c_short"}

    def test_containment_read_em_concentrates_on_anchor(self, monkeypatch):
        # c_full has an independent anchor read (aligns uniquely-best to it). The
        # "shared" read ties on both members and does not straddle a differing
        # junction -> ambiguous (1/k); EM then pulls it onto c_full (the abundance
        # the anchor establishes), starving the unsupported c_short. No 0.5/0.5 split.
        cands = [_cand("c_full", FULL), _cand("c_short", SHARED)]
        cands[0].supporting_read_ids.update({"anchor", "shared"})
        cands[1].supporting_read_ids.add("shared")
        cs = CandidateSet(
            interval=_interval(),
            candidates=cands,
            read_ids={"anchor", "shared"},
            read_sequences={"anchor": FULL, "shared": SHARED},
            clusters=[["c_full", "c_short"]],
        )
        _patch_bam(monkeypatch, [
            _FakeRead("anchor", 0, len(FULL)),
            _FakeRead("shared", len(PREFIX), len(FULL)),
        ])
        runner = _runner()
        out = {r.candidate_id: r
               for r in runner._quant_cluster(cs, ["anchor", "shared"], _interval())}
        assert out["c_full"].abundance > 1.99      # anchor + EM-pulled shared read
        assert out["c_short"].abundance < 0.01     # starved (EM abundance ~0)

    def test_clusters_none_falls_back_to_single_cluster(self, monkeypatch):
        # clusters=None -> the whole candidate list is treated as ONE cluster and
        # the mode still runs (unscoped), emitting a QuantResult per candidate.
        cands = [_cand("c_full", FULL), _cand("c_other", OTHER)]
        cands[0].supporting_read_ids.add("r1")
        cs = CandidateSet(
            interval=_interval(),
            candidates=cands,
            read_ids={"r1"},
            read_sequences={"r1": FULL},
            clusters=None,
        )
        _patch_bam(monkeypatch, [_FakeRead("r1", 0, len(FULL))])
        runner = _runner()
        out = {r.candidate_id: r
               for r in runner._quant_cluster(cs, ["r1"], _interval())}
        assert set(out) == {"c_full", "c_other"}
        assert out["c_full"].abundance == pytest.approx(1.0)
        assert out["c_other"].abundance == pytest.approx(0.0)

    def test_unreferenced_candidate_scored_via_extra_cluster(self, monkeypatch):
        # c_other is NOT listed in clusters (e.g. a post-discovery/fusion candidate).
        # It must still be scored (grouped into the one extra cluster), not emitted
        # with silent zero abundance.
        cands = [_cand("c_full", FULL), _cand("c_other", OTHER)]
        cands[0].supporting_read_ids.add("r1")
        cands[1].supporting_read_ids.add("r2")
        cs = CandidateSet(
            interval=_interval(),
            candidates=cands,
            read_ids={"r1", "r2"},
            read_sequences={"r1": FULL, "r2": OTHER},
            clusters=[["c_full"]],          # c_other deliberately unreferenced
        )
        _patch_bam(monkeypatch, [
            _FakeRead("r1", 0, len(FULL)),
            _FakeRead("r2", 0, len(OTHER)),
        ])
        out = {r.candidate_id: r
               for r in _runner()._quant_cluster(cs, ["r1", "r2"], _interval())}
        assert out["c_full"].abundance == pytest.approx(1.0)
        assert out["c_other"].abundance == pytest.approx(1.0)

    def test_clusters_threaded_from_candidate_set(self, monkeypatch):
        # Two disjoint single-member clusters: each read is scored ONLY within its
        # own cluster, so cross-cluster candidates never compete for a read.
        cands = [_cand("c_full", FULL), _cand("c_other", OTHER)]
        cs = CandidateSet(
            interval=_interval(),
            candidates=cands,
            read_ids={"r1", "r2"},
            read_sequences={"r1": FULL, "r2": OTHER},
            clusters=[["c_full"], ["c_other"]],
        )
        cands[0].supporting_read_ids.add("r1")
        cands[1].supporting_read_ids.add("r2")
        _patch_bam(monkeypatch, [
            _FakeRead("r1", 0, len(FULL)),
            _FakeRead("r2", 0, len(OTHER)),
        ])
        runner = _runner()
        out = {r.candidate_id: r
               for r in runner._quant_cluster(cs, ["r1", "r2"], _interval())}
        assert out["c_full"].abundance == pytest.approx(1.0)
        assert out["c_other"].abundance == pytest.approx(1.0)
