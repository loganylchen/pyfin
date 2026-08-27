"""Unit tests for ``m2_junction_nll.m2_resolve_tie`` (Stage M2-0).

The production tie-resolution entry point. Early-return paths (no window, <2
candidates) need no signal deps. The decision-aggregation logic is exercised by
monkeypatching ``read_cand_mean_nll`` (the krill leg is validated end-to-end
in the Docker ablation), so these tests run on a host without krill.
"""
from __future__ import annotations

import fin.scoring.m2_junction_nll as m2
from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.pipeline.assignment import _m2_score_region
from fin.pipeline.config import PipelineConfig


def _cand(cid, introns, start, end, strand="+", seq=None):
    # Default sequence length = spliced transcript length so _tx2genome maps
    # every in-window position.
    if seq is None:
        introns_s = sorted(introns)
        pos = start
        tx = 0
        for s, e in introns_s:
            if s > pos:
                tx += s - pos
            pos = max(pos, e)
        if end > pos:
            tx += end - pos
        seq = "A" * tx
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=end if strand == "+" else start,
        sequence=seq,
        source="novel",
        supporting_read_ids=set(),
        chrom="chrT",
        strand=strand,
        start=start,
        end=end,
    )


def _wobble_pair():
    """Two same-class candidates: 2nd intron donor shifted 400 -> 410.

    internal_diff_regions == [(398, 412)] -> class_junction_window_set non-empty.
    """
    a = _cand("a", [(200, 300), (400, 500)], 100, 600)
    b = _cand("b", [(200, 300), (410, 500)], 100, 600)
    return a, b


def test_assignment_summed_metric_uses_tight_windows():
    a, b = _wobble_pair()
    cfg = PipelineConfig(
        bam_path="/tmp/x.bam", m2_metric="summed_llr", m2_summed_llr_flank=6
    )
    gset, windows = _m2_score_region(cfg, [a, b])
    assert gset is None
    assert windows == m2.diff_junction_windows([a, b], flank=6)
    assert windows


def test_assignment_summed_metric_abstains_on_containment():
    short = _cand("short", [(400, 500)], 300, 600)
    long = _cand("long", [(200, 300), (400, 500)], 100, 600)
    cfg = PipelineConfig(bam_path="/tmp/x.bam", m2_metric="summed_llr")
    assert _m2_score_region(cfg, [short, long]) == (None, [])


def test_assignment_off_metric_has_no_signal_window():
    a, b = _wobble_pair()
    cfg = PipelineConfig(bam_path="/tmp/x.bam", m2_metric="off")
    assert _m2_score_region(cfg, [a, b]) == (None, [])


# --- early returns (no signal deps) ---------------------------------------


def test_fewer_than_two_candidates_returns_none():
    a, _ = _wobble_pair()
    assert m2.m2_resolve_tie("r0", "ACGT", [a], "/x.blow5") == (None, 0.0)


def test_empty_read_seq_returns_none():
    a, b = _wobble_pair()
    assert m2.m2_resolve_tie("r0", "", [a, b], "/x.blow5") == (None, 0.0)


def test_no_discrimination_window_returns_none():
    # Identical intron chains, only terminal start/end differ -> no internal
    # diff -> empty window -> M2 abstains (caller keeps its split default).
    a = _cand("a", [(200, 300), (400, 500)], 100, 600)
    d = _cand("d", [(200, 300), (400, 500)], 120, 650)
    assert m2.m2_resolve_tie("r0", "ACGT", [a, d], "/x.blow5") == (None, 0.0)


# --- decision aggregation (monkeypatched scorer) --------------------------


def _patch_scores(monkeypatch, scores):
    """scores: dict candidate_id -> (nll, n_events)."""

    def fake(read_id, read_seq, cand, windows, krill_aligner, mappy_aligner,
             sig_path, pore, gset=None, use_gpu=True, num_thread=0, reduce="mean"):
        return scores.get(cand.candidate_id, (float("nan"), 0))

    monkeypatch.setattr(m2, "read_cand_mean_nll", fake)


def test_summed_llr_metric_uses_tight_window_and_picks_lower(monkeypatch):
    # metric="summed_llr" builds diff_junction_windows (non-empty for a wobble pair)
    # and returns the lower-score candidate; margin is the summed gap.
    a, b = _wobble_pair()
    _patch_scores(monkeypatch, {"a": (3.0, 8), "b": (9.0, 8)})
    idx, margin = m2.m2_resolve_tie(
        "r0", "ACGT", [a, b], "/x.blow5",
        krill_aligner=object(), mappy_aligners=[object(), object()],
        metric="summed_llr",
    )
    assert idx == 0
    assert margin == 6.0


def test_summed_llr_unanimous_returns_none(monkeypatch):
    # No differing junction boundary -> diff_junction_windows empty -> abstain.
    a = _cand("a", [(200, 300), (400, 500)], 100, 600)
    d = _cand("d", [(200, 300), (400, 500)], 120, 650)
    _patch_scores(monkeypatch, {"a": (3.0, 8), "d": (9.0, 8)})
    assert m2.m2_resolve_tie(
        "r0", "ACGT", [a, d], "/x.blow5",
        krill_aligner=object(), mappy_aligners=[object(), object()],
        metric="summed_llr",
    ) == (None, 0.0)


def test_two_way_picks_lower_nll_with_margin(monkeypatch):
    a, b = _wobble_pair()
    _patch_scores(monkeypatch, {"a": (1.0, 5), "b": (2.5, 5)})
    idx, margin = m2.m2_resolve_tie(
        "r0", "ACGT", [a, b], "/x.blow5",
        krill_aligner=object(), mappy_aligners=[object(), object()],
    )
    assert idx == 0
    assert margin == 1.5


def test_winner_is_global_min_regardless_of_order(monkeypatch):
    a, b = _wobble_pair()
    _patch_scores(monkeypatch, {"a": (3.0, 5), "b": (0.5, 5)})
    idx, margin = m2.m2_resolve_tie(
        "r0", "ACGT", [a, b], "/x.blow5",
        krill_aligner=object(), mappy_aligners=[object(), object()],
    )
    assert idx == 1
    assert margin == 2.5


def test_single_scored_candidate_abstains(monkeypatch):
    # Technical scoring absence for b is not biological evidence against b.
    a, b = _wobble_pair()
    _patch_scores(monkeypatch, {"a": (1.2, 7), "b": (float("nan"), 0)})
    assert m2.m2_resolve_tie(
        "r0", "ACGT", [a, b], "/x.blow5",
        krill_aligner=object(), mappy_aligners=[object(), object()],
    ) == (None, 0.0)


def test_no_candidate_covered_returns_none(monkeypatch):
    a, b = _wobble_pair()
    _patch_scores(monkeypatch, {"a": (float("nan"), 0), "b": (float("nan"), 0)})
    out = m2.m2_resolve_tie(
        "r0", "ACGT", [a, b], "/x.blow5",
        krill_aligner=object(), mappy_aligners=[object(), object()],
    )
    assert out == (None, 0.0)


def test_none_mappy_aligner_causes_abstention(monkeypatch):
    a, b = _wobble_pair()
    _patch_scores(monkeypatch, {"a": (1.0, 5), "b": (0.1, 5)})
    assert m2.m2_resolve_tie(
        "r0", "ACGT", [a, b], "/x.blow5",
        krill_aligner=object(), mappy_aligners=[object(), None],
    ) == (None, 0.0)


# --- return_scored: score-gated fallback split (mk1) ----------------------


def test_return_scored_both_scored_best_first(monkeypatch):
    # Both candidates score; scored_idxs is sorted best (lowest NLL) first.
    a, b = _wobble_pair()
    _patch_scores(monkeypatch, {"a": (2.5, 5), "b": (1.0, 5)})
    idx, margin, scored = m2.m2_resolve_tie(
        "r0", "ACGT", [a, b], "/x.blow5",
        krill_aligner=object(), mappy_aligners=[object(), object()],
        return_scored=True,
    )
    assert idx == 1
    assert margin == 1.5
    assert scored == [1, 0]


def test_return_scored_mixed_preserves_full_m1_fallback(monkeypatch):
    a, b = _wobble_pair()
    _patch_scores(monkeypatch, {"a": (1.2, 7), "b": (float("nan"), 0)})
    assert m2.m2_resolve_tie(
        "r0", "ACGT", [a, b], "/x.blow5",
        krill_aligner=object(), mappy_aligners=[object(), object()],
        return_scored=True,
    ) == (None, 0.0, [])


def test_return_scored_none_scored_empty_list(monkeypatch):
    # No candidate scores -> (None, 0.0, []) so caller keeps its full split (M1).
    a, b = _wobble_pair()
    _patch_scores(monkeypatch, {"a": (float("nan"), 0), "b": (float("nan"), 0)})
    out = m2.m2_resolve_tie(
        "r0", "ACGT", [a, b], "/x.blow5",
        krill_aligner=object(), mappy_aligners=[object(), object()],
        return_scored=True,
    )
    assert out == (None, 0.0, [])


def test_return_scored_early_return_has_empty_list():
    # Early returns (no window) also yield an empty scored list under the flag.
    a, b = _wobble_pair()
    assert m2.m2_resolve_tie(
        "r0", "", [a, b], "/x.blow5", return_scored=True
    ) == (None, 0.0, [])
