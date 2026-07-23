"""Unit tests for per-cluster read assignment policy (fin/pipeline/cluster_quant.py).

Pure logic: mappy AS rows and the M2 resolver are injected, so no krill/mappy needed.
"""
from __future__ import annotations

from fin.pipeline.cluster_quant import (
    assign_cluster_reads,
    _differing_coords,
    _straddles_difference,
)

FULL = ((100, 200), (300, 400))      # 2-intron full-length
TRUNC = ((100, 200),)                # 3'-truncation: exact prefix sub-chain
WOBBLE = ((100, 204), (300, 400))    # first-intron acceptor wobbled 200->204


# --- helpers ---------------------------------------------------------------------
def test_differing_coords_wobble():
    # (300,400) is shared -> excluded; only the wobbled first intron's coords remain.
    assert _differing_coords([FULL, WOBBLE], [0, 1]) == {100, 200, 204}


def test_differing_coords_unanimous_empty():
    assert _differing_coords([FULL, FULL], [0, 1]) == set()


def test_straddles_true_when_span_covers_differing_intron():
    # FULL vs TRUNC differ by intron (300,400); a read spanning it straddles.
    assert _straddles_difference((50, 450), [FULL, TRUNC], [0, 1]) is True


def test_straddles_false_for_containment_read():
    # read only covers the shared (100,200) region -> does not cross (300,400).
    assert _straddles_difference((90, 260), [FULL, TRUNC], [0, 1]) is False


# --- unique-best anchors ---------------------------------------------------------
def test_unique_best_assigns_100pct_and_survives():
    a = assign_cluster_reads(
        [FULL, TRUNC],
        as_rows={"r": {0: 1000.0}},
        spans={"r": (100, 400)},
    )
    assert a.weights["r"] == {0: 1.0}
    assert a.n_unique == 1
    assert 0 in a.survivors and 1 not in a.survivors


# --- containment tie -> ambiguous 1/k -> EM concentrates on the anchored member --
def test_containment_tie_em_concentrates_on_anchor():
    # 'anchor' uniquely supports FULL; 'cont' is an ambiguous containment tie (1/k).
    # EM redistributes cont toward FULL (the abundance the anchor establishes); the
    # unsupported TRUNC starves rather than being hard-routed.
    a = assign_cluster_reads(
        [FULL, TRUNC],
        as_rows={"anchor": {0: 1000.0}, "cont": {0: 500.0, 1: 500.0}},
        spans={"anchor": (100, 400), "cont": (95, 260)},
    )
    assert a.weights["anchor"] == {0: 1.0}
    assert a.weights["cont"][0] > 0.99          # EM pulls cont onto FULL
    assert a.weights["cont"][1] < 0.01
    assert a.n_containment == 1
    assert a.survivors == {0}                    # TRUNC (EM abundance ~0) starved


def test_truncation_with_own_support_survives():
    # TRUNC has its own unique read (a genuine short isoform) -> both survive.
    a = assign_cluster_reads(
        [FULL, TRUNC],
        as_rows={"f": {0: 1000.0}, "t": {1: 900.0}},
        spans={"f": (100, 400), "t": (100, 200)},
        min_support=1.0,
    )
    assert a.survivors == {0, 1}


# --- straddling tie -> M2 (strict best, no margin threshold) ----------------------
def test_straddling_tie_m2_strict_best_assigns_winner():
    calls = []
    def m2(rid, tie):
        calls.append((rid, tuple(tie)))
        return (0, 3.5)                        # strict winner (margin > 0)
    a = assign_cluster_reads(
        [FULL, WOBBLE],
        as_rows={"r": {0: 1000.0, 1: 1000.0}},
        spans={"r": (50, 450)},                # straddles the wobbled intron
        m2_resolve=m2,
    )
    assert calls == [("r", (0, 1))]
    assert a.weights["r"] == {0: 1.0}
    assert a.n_m2_called == 1 and a.n_m2_confident == 1


def test_straddling_tie_tiny_margin_still_wins():
    # No margin threshold: any strictly-positive margin makes the read certain.
    def m2(rid, tie):
        return (1, 0.05)                       # strict winner = member 1
    a = assign_cluster_reads(
        [FULL, WOBBLE],
        as_rows={"r": {0: 1000.0, 1: 1000.0}},
        spans={"r": (50, 450)},
        m2_resolve=m2,
    )
    assert a.weights["r"] == {1: 1.0}
    assert a.n_m2_confident == 1


def test_straddling_tie_nll_tie_stays_ambiguous():
    # margin == 0 (an NLL tie) is NOT a strict best -> read stays ambiguous (1/k) ->
    # EM then follows the anchor.
    def m2(rid, tie):
        return (0, 0.0)
    a = assign_cluster_reads(
        [FULL, WOBBLE],
        as_rows={"anchor": {0: 1000.0}, "r": {0: 1000.0, 1: 1000.0}},
        spans={"anchor": (50, 450), "r": (50, 450)},
        m2_resolve=m2,
    )
    assert a.n_deferred == 1 and a.n_m2_confident == 0
    assert a.weights["r"][0] > a.weights["r"][1]   # EM pulls toward the anchor


def test_ambiguous_even_split_when_no_anchor():
    def m2(rid, tie):
        return None                            # abstain -> ambiguous
    a = assign_cluster_reads(
        [FULL, WOBBLE],
        as_rows={"r": {0: 1000.0, 1: 1000.0}},
        spans={"r": (50, 450)},
        m2_resolve=m2,
    )
    assert a.weights["r"] == {0: 0.5, 1: 0.5}   # symmetric, no evidence
    assert a.n_deferred == 1


def test_no_m2_resolver_leaves_straddling_ambiguous():
    a = assign_cluster_reads(
        [FULL, WOBBLE],
        as_rows={"r": {0: 1000.0, 1: 1000.0}},
        spans={"r": (50, 450)},
    )
    assert a.n_m2_called == 0 and a.n_deferred == 1


def test_m1_tie_margin_makes_near_ties_ambiguous():
    # best beats 2nd by only 10 AS (< margin 20) -> NOT unique-best -> the two are
    # TIED (ambiguous), so a near-tie wobble read does not anchor the shadow.
    a = assign_cluster_reads(
        [FULL, WOBBLE],
        as_rows={"r": {0: 2000.0, 1: 1990.0}},
        spans={"r": (50, 450)},          # straddles -> but no m2 -> ambiguous 1/k
        m1_tie_margin=20.0,
    )
    assert a.n_unique == 0
    assert a.weights["r"] == {0: 0.5, 1: 0.5}


def test_m1_tie_margin_keeps_clear_best_unique():
    # gap 40 > margin 20 -> member 0 is a genuine unique-best anchor.
    a = assign_cluster_reads(
        [FULL, WOBBLE],
        as_rows={"r": {0: 2000.0, 1: 1960.0}},
        spans={"r": (50, 450)},
        m1_tie_margin=20.0,
    )
    assert a.n_unique == 1
    assert a.weights["r"] == {0: 1.0}


def test_min_support_threshold_drops_weak_member():
    a = assign_cluster_reads(
        [FULL, TRUNC],
        as_rows={"f1": {0: 1000.0}, "f2": {0: 1000.0}, "t": {1: 900.0}},
        spans={"f1": (100, 400), "f2": (100, 400), "t": (100, 200)},
        min_support=2.0,
    )
    assert a.survivors == {0}                  # member 1 has only 1 read < 2
