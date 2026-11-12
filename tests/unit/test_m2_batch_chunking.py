"""The M2 batched eventalign is chunked to bound memory; chunking must not enter the result.

The backend here is a fake whose answer for a (read, candidate) pair depends only
on that pair -- which is the property of krill's CPU eventalign that the claim
rests on. What these tests pin is everything on the PyFIN side: the payload
handed to the backend, the reduction, the coverage gate, the covered-vote prior,
the EM responsibilities and the support gate.
"""
from __future__ import annotations

import random
import sys
import types
import zlib

import numpy as np
import pytest

import fin.pipeline.assignment as assignment
import fin.scoring.krill_aligner as krill_aligner
import fin.scoring.m2_junction_nll as m2
import fin.scoring.mappy_score as mappy_score
from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.pipeline.config import PipelineConfig

N_READS = 12
LABELS = "abc"
DONORS = {"a": 400, "b": 410, "c": 420}
SEQS = {"a": "ACGT" * 100, "b": "CCGT" * 100, "c": "GCGT" * 100}
SHIFT = {"a": 0, "b": 3, "c": 7}
WHOLE = 10**9


def _read_ids():
    return [f"r{k:02d}" for k in range(N_READS)]


def _read_seq(k: int) -> str:
    return "A" * (5 + k)


def _tie(k: int):
    """The read's exact best-AS tie, by construction."""
    return {0: "abc", 1: "ab", 2: "bc", 3: "a"}[k % 4]


def _bounds(k: int, label: str):
    """Where the read's best hit lies on the candidate: the independent oracle."""
    start = len(_read_seq(k)) + SHIFT[label]
    return start, start + 40 + len(_read_seq(k))


def _expected_payload():
    """What the backend must be handed for the whole interval, from first principles."""
    variants, starts = {}, {}
    for k, rid in enumerate(_read_ids()):
        if len(_tie(k)) < 2:
            continue
        variants[rid] = {c: SEQS[c][slice(*_bounds(k, c))] for c in _tie(k)}
        starts[rid] = {c: _bounds(k, c)[0] for c in _tie(k)}
    return variants, starts


def _expected_nll(k: int, label: str) -> float:
    start, end = _bounds(k, label)
    content = (zlib.crc32(SEQS[label][start:end].encode()) % 7) * 1e-3
    if k % 3 == 0:  # the signal separates the structures: margin 4 >> threshold 1
        return {"a": 1.0, "b": 5.0, "c": 9.0}[label] + content
    return 2.0 + content  # indistinguishable: margin < threshold


def _candidate(cid: str) -> TranscriptCandidate:
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=((200, 300), (DONORS[cid], 500))),
        three_prime_pos=600,
        sequence=SEQS[cid],
        source="novel",
        supporting_read_ids=set(),
        chrom="chrT",
        strand="+",
        start=100,
        end=600,
    )


class _Hit:
    def __init__(self, label: str, k: int):
        self.label, self.k = label, k
        self.r_st, self.r_en = _bounds(k, label)


class _Aligner:
    def __init__(self, seq=None, preset=None):
        self.label = next(c for c in LABELS if SEQS[c] == seq)

    def map(self, read_sequence: str):
        return [_Hit(self.label, len(read_sequence) - 5)]


def _config(batch_pairs, gate, metric="summed_llr"):
    return PipelineConfig(
        bam_path="/tmp/input.bam",
        signal_path="/tmp/input.blow5",
        m2_metric=metric,
        m2_summed_llr_flank=4,
        m2_summed_llr_margin=1.0,
        m2_diff_cover_gate=gate,
        use_gpu=False,
        m2_batch_pairs=batch_pairs,
    )


class _Backend:
    """A fake krill that checks every payload against the oracle."""

    def __init__(self, *, fail_call=None, none_call=None, drop_read=None, shuffle=False,
                 batch_api=True):
        self.calls, self.starts, self.singular = [], [], []
        self.fail_call, self.none_call = fail_call, none_call
        self.drop_read, self.shuffle = drop_read, shuffle
        self.module = types.ModuleType("krill")
        if batch_api:
            self.module.align_reads_variants = self.align_reads_variants

    def align_reads_variants(self, _signal, reads_variants, *, start, **_kwargs):
        index = len(self.calls)
        self.calls.append({rid: dict(v) for rid, v in reads_variants.items()})
        self.starts.append({rid: dict(v) for rid, v in start.items()})
        exp_variants, exp_starts = _expected_payload()
        for rid, variants in reads_variants.items():
            assert list(variants.items()) == list(exp_variants[rid].items())
            assert list(start[rid].items()) == list(exp_starts[rid].items())
        assert list(start) == list(reads_variants)
        if index == self.fail_call:
            raise RuntimeError("simulated krill failure")
        if index == self.none_call:
            return None
        rids = [rid for rid in reads_variants if rid != self.drop_read]
        if self.shuffle:
            random.Random(index).shuffle(rids)
        return {
            rid: [
                {"variant_label": label, "status": 0, "k": int(rid[1:])}
                for label in reads_variants[rid]
            ]
            for rid in rids
        }

    def read_cand_mean_nll(self, rid, _seq, cand, _windows, *_a, **_k):
        self.singular.append((rid, cand.candidate_id))
        return 1000.0 + ord(cand.candidate_id), 5


def _install(monkeypatch, backend):
    monkeypatch.setitem(sys.modules, "krill", backend.module)
    fake_mappy = types.ModuleType("mappy")
    fake_mappy.Aligner = _Aligner
    monkeypatch.setitem(sys.modules, "mappy", fake_mappy)
    monkeypatch.setattr(
        krill_aligner, "make_krill_aligner", lambda *_a, **_k: (object(), False)
    )
    monkeypatch.setattr(krill_aligner, "krill_thread_count", lambda: 1)
    monkeypatch.setattr(
        mappy_score, "score_hit", lambda hit: 10.0 if hit.label in _tie(hit.k) else 5.0
    )
    monkeypatch.setattr(
        m2, "_mean_nll_in_window",
        lambda result, cand, _windows, reduce="mean": (
            _expected_nll(result["k"], cand.candidate_id), 5),
    )
    monkeypatch.setattr(
        m2, "_mean_nll_in_gset",
        lambda result, cand, _gset: (_expected_nll(result["k"], cand.candidate_id), 5),
    )
    # even reads straddle every disagreement span, odd reads sit inside one
    monkeypatch.setattr(
        m2, "event_genomic_positions",
        lambda result, _cand: [0, 1000] if result["k"] % 2 == 0 else [450],
    )
    monkeypatch.setattr(m2, "read_cand_mean_nll", backend.read_cand_mean_nll)


def _tie_nll(monkeypatch, batch_pairs, gate=False, metric="summed_llr", **backend_kwargs):
    backend = _Backend(**backend_kwargs)
    _install(monkeypatch, backend)
    candidates = [_candidate(c) for c in LABELS]
    raw = np.array(
        [[10.0 if c in _tie(k) else 5.0 for c in LABELS] for k in range(N_READS)]
    )
    out = assignment.tie_nll(
        _config(batch_pairs, gate, metric),
        _read_ids(),
        {rid: _read_seq(k) for k, rid in enumerate(_read_ids())},
        candidates,
        [_Aligner(seq=SEQS[c]) for c in LABELS],
        raw,
    )
    return out, backend


def _assign(monkeypatch, batch_pairs, **backend_kwargs):
    backend = _Backend(**backend_kwargs)
    _install(monkeypatch, backend)
    config = _config(batch_pairs, gate=True)
    candidate_set = types.SimpleNamespace(
        candidates=[_candidate(c) for c in LABELS],
        read_sequences={rid: _read_seq(k) for k, rid in enumerate(_read_ids())},
    )
    quant = assignment.Assigner(config).assign(
        candidate_set,
        _read_ids(),
        tie_nll_fn=lambda *a: assignment.tie_nll(config, *a),
        eff_lengths_fn=lambda _cands: None,
    )
    return quant, backend


def _same_quant(a, b):
    assert np.array_equal(a.d_tx, b.d_tx)
    assert np.array_equal(a.R, b.R)
    assert np.array_equal(a.hard_assignments, b.hard_assignments)
    assert np.array_equal(a.structural_d_tx_raw, b.structural_d_tx_raw)
    assert a.nlls_by_read == b.nlls_by_read
    assert a.ties_by_read == b.ties_by_read
    assert a.cover_by_read == b.cover_by_read
    assert (a.n_ties, a.n_refined) == (b.n_ties, b.n_refined)
    drops = [
        m2.support_gate_drops(q.cand_list, q.ties_by_read, q.nlls_by_read, tie_ok=True)
        for q in (a, b)
    ]
    assert drops[0] == drops[1]


BOUNDS = [1, 2, 3, 5, 8, 64]


def test_one_call_carries_the_whole_interval_exactly(monkeypatch):
    (nlls, ties, n_ties, n_refined, cover), backend = _tie_nll(monkeypatch, WHOLE)
    exp_variants, exp_starts = _expected_payload()

    assert len(backend.calls) == 1
    assert list(backend.calls[0]) == list(exp_variants)  # read order
    assert backend.calls[0] == exp_variants and backend.starts[0] == exp_starts
    assert n_ties == n_refined == len(exp_variants) == 9
    assert cover == {}
    for k in range(N_READS):
        tie = [LABELS.index(c) for c in _tie(k)]
        assert ties[k] == tie
        if len(tie) < 2:
            assert k not in nlls
        else:
            assert nlls[k] == {j: _expected_nll(k, LABELS[j]) for j in tie}


@pytest.mark.parametrize("batch_pairs", BOUNDS)
def test_chunks_partition_the_payload_and_reproduce_the_result(monkeypatch, batch_pairs):
    whole, _ = _tie_nll(monkeypatch, WHOLE)
    chunked, backend = _tie_nll(monkeypatch, batch_pairs)
    exp_variants, exp_starts = _expected_payload()

    assert chunked == whole
    assert backend.singular == []
    merged_variants, merged_starts = {}, {}
    for call, start in zip(backend.calls, backend.starts):
        assert not set(call) & set(merged_variants)
        merged_variants.update(call)
        merged_starts.update(start)
    assert list(merged_variants) == list(exp_variants)
    assert merged_variants == exp_variants and merged_starts == exp_starts


@pytest.mark.parametrize("metric", ["mean", "sqrt_count_mean_llr"])
@pytest.mark.parametrize("batch_pairs", [1, 5, 64])
def test_the_other_reducers_are_independent_of_the_bound_too(monkeypatch, metric, batch_pairs):
    whole, whole_backend = _tie_nll(monkeypatch, WHOLE, metric=metric)
    chunked, backend = _tie_nll(monkeypatch, batch_pairs, metric=metric)
    assert len(whole_backend.calls) == 1 and whole[3] == 9
    assert chunked == whole
    assert sum(len(v) for call in backend.calls for v in call.values()) == 21


@pytest.mark.parametrize("batch_pairs", BOUNDS)
def test_a_call_exceeds_the_bound_by_less_than_one_read(monkeypatch, batch_pairs):
    _, backend = _tie_nll(monkeypatch, batch_pairs)
    sizes = [sum(len(v) for v in call.values()) for call in backend.calls]
    largest_tie = 3
    assert all(size < batch_pairs + largest_tie for size in sizes)
    assert all(size >= batch_pairs for size in sizes[:-1])  # only the last may be short
    assert sum(sizes) == 21  # 3 reads x 3 + 6 reads x 2


@pytest.mark.parametrize("batch_pairs", BOUNDS)
@pytest.mark.parametrize("shuffle", [False, True])
def test_gated_assignment_is_independent_of_the_bound(monkeypatch, batch_pairs, shuffle):
    whole, _ = _assign(monkeypatch, WHOLE)
    chunked, _ = _assign(monkeypatch, batch_pairs, shuffle=shuffle)
    _same_quant(whole, chunked)


def test_gated_assignment_matches_first_principles(monkeypatch):
    quant, _ = _assign(monkeypatch, 4)
    a, b, c = 0, 1, 2
    missing = np.float32(assignment.MISSING)
    d = quant.d_tx

    # even reads are covered, odd reads are not -- for every refined read
    assert quant.cover_by_read == {k: k % 2 == 0 for k in range(N_READS) if k % 4 != 3}
    # read 0: tie abc, separated, covered -> one-hot on a, and a vote for a
    assert list(d[0]) == [0.0, missing, missing]
    # read 6: tie bc, separated, covered -> one-hot on b, and a vote for b
    assert list(d[6]) == [missing, 0.0, missing]
    # read 9: tie ab, separated but NOT covered -> one-hot on a, no vote
    assert list(d[9]) == [0.0, missing, missing]
    # read 3: sole best a, never a tie
    assert list(d[3]) == [0.0, missing, missing]
    # the prior is therefore a:1, b:1, c:0. Ambiguous reads follow it:
    half = np.float32(-np.log(0.5))
    assert list(d[1]) == [half, half, missing]        # tie ab -> 1/2, 1/2
    assert list(d[2]) == [missing, 0.0, missing]      # tie bc -> all on b
    assert list(d[4]) == [half, half, missing]        # tie abc -> c has no vote
    assert quant.hard_assignments[0] == a and quant.hard_assignments[6] == b
    assert quant.hard_assignments[2] == b
    assert c not in set(quant.hard_assignments.tolist())


@pytest.mark.parametrize("mode", ["raise", "none"])
def test_a_failed_chunk_falls_back_per_read_without_touching_the_others(monkeypatch, mode):
    whole, _ = _tie_nll(monkeypatch, WHOLE)
    kwargs = {"fail_call": 1} if mode == "raise" else {"none_call": 1}
    (nlls, ties, n_ties, n_refined, _cover), backend = _tie_nll(monkeypatch, 5, **kwargs)

    failed = list(backend.calls[1])
    assert failed == ["r02", "r04"]
    assert sorted({rid for rid, _ in backend.singular}) == failed
    assert (ties, n_ties, n_refined) == (whole[1], whole[2], whole[3])
    for k in whole[0]:
        if f"r{k:02d}" in failed:
            assert nlls[k] == {LABELS.index(c): 1000.0 + ord(c) for c in _tie(k)}
        else:
            assert nlls[k] == whole[0][k]


def test_a_failed_chunk_keeps_coverage_elsewhere_and_never_seeds_the_prior(monkeypatch):
    whole, _ = _assign(monkeypatch, WHOLE)
    # bound 1 puts every read in its own call; the third call is read 2 alone
    quant, backend = _assign(monkeypatch, 1, fail_call=2)
    missing = np.float32(assignment.MISSING)
    half = np.float32(-np.log(0.5))

    assert list(backend.calls[2]) == ["r02"]
    assert {rid for rid, _ in backend.singular} == {"r02"}
    # coverage is known only for reads whose call was batch-aligned
    assert quant.cover_by_read == {k: v for k, v in whole.cover_by_read.items() if k != 2}
    # the fallback scores separate b from c (margin 1 >= threshold 1), so read 2
    # is hard-assigned to b -- but with no coverage it must not vote for b
    assert list(quant.d_tx[2]) == [missing, 0.0, missing]
    # the prior is therefore still a:1, b:1 (reads 0 and 6). Had read 2 voted it
    # would be a:1, b:2, and these ambiguous reads would split 1/3 : 2/3
    assert list(quant.d_tx[1]) == [half, half, missing]
    assert list(quant.d_tx[4]) == [half, half, missing]
    for k in range(N_READS):
        if k != 2:
            assert np.array_equal(quant.d_tx[k], whole.d_tx[k])


def test_a_read_missing_from_one_chunk_abstains_alone(monkeypatch):
    whole, _ = _tie_nll(monkeypatch, WHOLE)
    (nlls, _ties, _n_ties, n_refined, _cover), backend = _tie_nll(
        monkeypatch, 5, drop_read="r05"
    )
    assert backend.singular == []
    assert n_refined == whole[3] - 1
    assert nlls == {k: v for k, v in whole[0].items() if k != 5}


def test_a_backend_without_the_batch_api_scores_every_read_singly(monkeypatch):
    (nlls, _ties, n_ties, n_refined, cover), backend = _tie_nll(
        monkeypatch, 5, batch_api=False
    )
    assert backend.calls == []
    assert n_ties == n_refined == 9 and cover == {}
    assert [rid for rid, _ in backend.singular] == [
        rid for k, rid in enumerate(_read_ids()) for _c in _tie(k) if len(_tie(k)) > 1
    ]
    assert all(
        nlls[k] == {LABELS.index(c): 1000.0 + ord(c) for c in _tie(k)} for k in nlls
    )


def test_pair_bounded_chunks_keeps_reads_whole_and_ordered():
    per_read = {"x": {"a": 0, "b": 1}, "y": {"a": 0}, "z": {"a": 0, "b": 1, "c": 2}, "w": {"a": 0}}
    assert list(assignment._pair_bounded_chunks(per_read, 3)) == [["x", "y"], ["z"], ["w"]]
    assert list(assignment._pair_bounded_chunks(per_read, 1)) == [["x"], ["y"], ["z"], ["w"]]
    assert list(assignment._pair_bounded_chunks(per_read, 100)) == [["x", "y", "z", "w"]]
    assert list(assignment._pair_bounded_chunks({}, 5)) == []


def test_batch_bound_must_be_positive(tmp_path):
    bam, signal = tmp_path / "input.bam", tmp_path / "input.blow5"
    bam.touch()
    signal.touch()
    PipelineConfig(bam_path=str(bam), signal_path=str(signal), m2_batch_pairs=1).validate()
    config = PipelineConfig(bam_path=str(bam), signal_path=str(signal), m2_batch_pairs=0)
    with pytest.raises(ValueError, match="m2_batch_pairs"):
        config.validate()
