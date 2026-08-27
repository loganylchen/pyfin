"""Regression tests for per-read M2 eventalign input validation."""
from __future__ import annotations

import sys
import types

import numpy as np

import fin.pipeline.assignment as assignment
import fin.scoring.krill_aligner as krill_aligner
import fin.scoring.m2_junction_nll as m2
import fin.scoring.mappy_score as mappy_score
from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.pipeline.config import PipelineConfig


def _candidate(cid: str, donor: int, sequence: str) -> TranscriptCandidate:
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(
            introns=((200, 300), (donor, 500)),
        ),
        three_prime_pos=600,
        sequence=sequence,
        source="novel",
        supporting_read_ids=set(),
        chrom="chrT",
        strand="+",
        start=100,
        end=600,
    )


class _Hit:
    def __init__(self, start: int, end: int):
        self.r_st = start
        self.r_en = end


class _Aligner:
    def __init__(self, invalid_slice_for_c_read: bool = False):
        self.invalid_slice_for_c_read = invalid_slice_for_c_read

    def map(self, read_sequence: str):
        if self.invalid_slice_for_c_read and read_sequence.startswith("C"):
            return [_Hit(90, 111)]
        return [_Hit(0, 50)]


def test_invalid_candidate_slice_abstains_without_poisoning_batch(monkeypatch):
    payloads = []
    fake_krill = types.ModuleType("krill")

    def align_reads_variants(_signal, reads_variants, **_kwargs):
        payloads.append(reads_variants)
        assert set(reads_variants) == {"valid"}
        assert all(
            set(sequence) <= set("ACGTUacgtu")
            for variants in reads_variants.values()
            for sequence in variants.values()
        )
        return {
            rid: [
                {"variant_label": label, "status": 0}
                for label in variants
            ]
            for rid, variants in reads_variants.items()
        }

    fake_krill.align_reads_variants = align_reads_variants
    monkeypatch.setitem(sys.modules, "krill", fake_krill)
    monkeypatch.setattr(
        krill_aligner,
        "make_krill_aligner",
        lambda *_args, **_kwargs: (object(), False),
    )
    monkeypatch.setattr(krill_aligner, "krill_thread_count", lambda: 1)
    monkeypatch.setattr(mappy_score, "score_hit", lambda _hit: 10.0)
    monkeypatch.setattr(
        m2,
        "_mean_nll_in_window",
        lambda result, _candidate, _windows, reduce="mean": (
            (1.0 if result["variant_label"] == "a" else 2.0),
            5,
        ),
    )

    seq_a = "A" * 300
    seq_b = "A" * 100 + "N" + "A" * 199
    candidates = [
        _candidate("a", 400, seq_a),
        _candidate("b", 410, seq_b),
    ]
    config = PipelineConfig(
        bam_path="/tmp/input.bam",
        signal_path="/tmp/input.blow5",
        m2_metric="summed_llr",
        m2_summed_llr_flank=4,
        m2_diff_cover_gate=False,
        use_gpu=False,
    )

    nlls, ties, n_ties, n_refined, coverage = assignment.tie_nll(
        config,
        ["valid", "invalid"],
        {"valid": "AAAA", "invalid": "CCCC"},
        candidates,
        [_Aligner(), _Aligner(invalid_slice_for_c_read=True)],
        np.array([[10.0, 10.0], [10.0, 10.0]]),
    )

    assert len(payloads) == 1
    assert ties == {0: [0, 1], 1: [0, 1]}
    assert nlls == {0: {0: 1.0, 1: 2.0}}
    assert n_ties == 2
    assert n_refined == 1
    assert coverage == {}
