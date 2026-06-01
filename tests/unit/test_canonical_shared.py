"""Unit tests for the shared canonical splice primitives.

Covers ``fin.candidates.canonical``: ``parse_motifs``, ``canonical_intron_alts``
and ``chain_all_canonical`` (the candidate-independent pure chain version).
These are the single source of truth imported by both the production pipeline
and the ablation harness, so the harness-side behaviour (exercised by
``tests/unit/test_canonical_search.py``) and these must stay in lock-step.
"""
from __future__ import annotations

import pytest

from fin.candidates.canonical import (
    MOTIFS_EXTENDED,
    MOTIFS_STRICT,
    canonical_intron_alts,
    chain_all_canonical,
    parse_motifs,
)

GTAG = frozenset({("GT", "AG")})


def _genome_single() -> str:
    """40-bp genome: canonical GT at 10 and AG at both 16 and 18."""
    g = ["C"] * 40
    g[10], g[11] = "G", "T"        # donor GT, intron start 10
    g[16], g[17] = "A", "G"        # acceptor AG -> ne=18
    g[18], g[19] = "A", "G"        # acceptor AG -> ne=20
    return "".join(g)


# --------------------------------------------------------------------------- #
# parse_motifs
# --------------------------------------------------------------------------- #
def test_parse_motifs_default_empty():
    assert parse_motifs([]) == GTAG
    assert parse_motifs(None) == GTAG


def test_parse_motifs_multi_and_uppercase():
    got = parse_motifs(["gt-ag", "GC-AG", "AT-AC"])
    assert got == MOTIFS_EXTENDED


def test_parse_motifs_invalid_token():
    with pytest.raises(ValueError):
        parse_motifs(["GTAG"])
    with pytest.raises(ValueError):
        parse_motifs(["G-AG"])
    with pytest.raises(ValueError):
        parse_motifs(["GT-A-AG"])


def test_module_constants():
    assert MOTIFS_STRICT == GTAG
    assert ("GC", "AG") in MOTIFS_EXTENDED and ("AT", "AC") in MOTIFS_EXTENDED


# --------------------------------------------------------------------------- #
# canonical_intron_alts
# --------------------------------------------------------------------------- #
def test_alts_plus_sorted_by_distance():
    g = _genome_single()
    alts = canonical_intron_alts(10, 20, "+", g, GTAG, 2)
    assert alts == [(10, 20), (10, 18)]


def test_alts_none_when_no_motif():
    assert canonical_intron_alts(5, 15, "+", "C" * 40, GTAG, 2) == []


def test_alts_respects_window():
    g = _genome_single()
    assert canonical_intron_alts(10, 20, "+", g, GTAG, 0) == [(10, 20)]


def test_alts_minus_strand():
    g = list("C" * 40)
    g[10], g[11] = "C", "T"        # s=10 -> acceptor revcomp("CT")="AG"
    g[18], g[19] = "A", "C"        # e=20 -> donor revcomp("AC")="GT"
    alts = canonical_intron_alts(10, 20, "-", "".join(g), GTAG, 2)
    assert (10, 20) in alts


# --------------------------------------------------------------------------- #
# chain_all_canonical
# --------------------------------------------------------------------------- #
def test_chain_empty_passes():
    assert chain_all_canonical((), _genome_single(), "+") is True


def test_chain_canonical_plus():
    assert chain_all_canonical([(10, 20)], _genome_single(), "+", GTAG) is True


def test_chain_noncanonical_rejected():
    assert chain_all_canonical([(5, 15)], _genome_single(), "+", GTAG) is False


def test_chain_out_of_bounds_rejected():
    assert chain_all_canonical([(0, 1)], _genome_single(), "+", GTAG) is False
    assert chain_all_canonical([(38, 60)], _genome_single(), "+", GTAG) is False


def test_chain_empty_genome_rejected():
    assert chain_all_canonical([(10, 20)], "", "+", GTAG) is False


def test_chain_extended_motif():
    # GC-AG intron at (10, 20): donor GC, acceptor AG.
    g = ["C"] * 40
    g[10], g[11] = "G", "C"        # donor GC
    g[18], g[19] = "A", "G"        # acceptor AG
    genome = "".join(g)
    assert chain_all_canonical([(10, 20)], genome, "+", GTAG) is False
    assert chain_all_canonical([(10, 20)], genome, "+", MOTIFS_EXTENDED) is True
