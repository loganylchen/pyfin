"""Tests for the motif-aware Stage C expansion (expand_canonical_chain_alternatives_v2).

The v2 expansion enumerates PAIRED canonical (s', e') alternatives near each
read-derived junction via the shared ``canonical.canonical_intron_alts`` grid,
fixing the legacy donor/acceptor-independent Cartesian bug. The original chain
is ALWAYS retained.
"""

from __future__ import annotations

from fin.candidates.canonical import MOTIFS_EXTENDED, MOTIFS_STRICT
from fin.candidates.dataclasses import IntronChain
from fin.candidates.intron_chains import expand_canonical_chain_alternatives_v2

# Genome layout (0-based), mirrors test_canonical_gate:
#   0-4   "AAAAA"
#   5-6   "GT"      <- canonical donor
#   7-12  "CCCCCC"
#   13-14 "AG"      <- canonical acceptor
#   15-19 "AAAAA"
# Canonical + strand intron = (5, 15): genome[5:7]="GT", genome[13:15]="AG".
GENOME = "AAAAA" + "GT" + "CCCCCC" + "AG" + "AAAAA"
CANON = (5, 15)
# Off-by-one donor: genome[6:8]="TC" (non-canonical), acceptor still AG.
WOBBLE = (6, 15)


def _rc(qname, introns):
    return ({"query_name": qname}, IntronChain(introns=tuple(introns)))


class TestExpandV2:
    def test_search_bp_zero_returns_empty(self):
        out = expand_canonical_chain_alternatives_v2(
            [_rc("r1", [WOBBLE])], GENOME, "+", MOTIFS_EXTENDED, search_bp=0
        )
        assert out == {}

    def test_empty_genome_returns_empty(self):
        out = expand_canonical_chain_alternatives_v2(
            [_rc("r1", [WOBBLE])], "", "+", MOTIFS_EXTENDED, search_bp=4
        )
        assert out == {}

    def test_original_chain_always_first(self):
        out = expand_canonical_chain_alternatives_v2(
            [_rc("r1", [WOBBLE])], GENOME, "+", MOTIFS_EXTENDED, search_bp=4
        )
        chains = out["r1"]
        assert chains[0].introns == (WOBBLE,)

    def test_wobble_recovers_canonical_alternative(self):
        out = expand_canonical_chain_alternatives_v2(
            [_rc("r1", [WOBBLE])], GENOME, "+", MOTIFS_EXTENDED, search_bp=4
        )
        intron_sets = {c.introns for c in out["r1"]}
        assert (WOBBLE,) in intron_sets
        assert (CANON,) in intron_sets

    def test_mono_read_keeps_only_original(self):
        out = expand_canonical_chain_alternatives_v2(
            [_rc("mono", [])], GENOME, "+", MOTIFS_EXTENDED, search_bp=4
        )
        assert out["mono"] == [IntronChain(introns=())]

    def test_already_canonical_retained(self):
        # A read already on the canonical junction keeps it (and may list itself
        # once only — no duplicate).
        out = expand_canonical_chain_alternatives_v2(
            [_rc("r1", [CANON])], GENOME, "+", MOTIFS_STRICT, search_bp=4
        )
        intron_sets = [c.introns for c in out["r1"]]
        assert (CANON,) in intron_sets
        assert intron_sets.count((CANON,)) == 1

    def test_paired_grid_no_illegal_pairings(self):
        # With only GT-AG strict, the single canonical alt near WOBBLE is CANON.
        # The legacy independent-donor/acceptor Cartesian could emit spurious
        # pairings; the paired grid must not.
        out = expand_canonical_chain_alternatives_v2(
            [_rc("r1", [WOBBLE])], GENOME, "+", MOTIFS_STRICT, search_bp=4
        )
        intron_sets = {c.introns for c in out["r1"]}
        assert intron_sets == {(WOBBLE,), (CANON,)}

    def test_max_chains_cap_keeps_only_original(self):
        # Cap of 1 forces "no expansion": only the original chain survives.
        out = expand_canonical_chain_alternatives_v2(
            [_rc("r1", [WOBBLE])], GENOME, "+", MOTIFS_EXTENDED,
            search_bp=4, max_chains_per_read=1,
        )
        assert [c.introns for c in out["r1"]] == [(WOBBLE,)]

    def test_minus_strand_canonical_alt(self):
        # - strand intron (s,e): donor=revcomp(genome[e-2:e]),
        # acceptor=revcomp(genome[s:s+2]). Canonical GT-AG needs
        # genome[e-2:e]="AC" and genome[s:s+2]="CT".
        genome = "AAAAA" + "CT" + "CCCCCC" + "AC" + "AAAAA"
        canon = (5, 15)
        wobble = (5, 14)  # acceptor genome[12:14]="CA" -> revcomp "TG" non-canon
        out = expand_canonical_chain_alternatives_v2(
            [_rc("m", [wobble])], genome, "-", MOTIFS_STRICT, search_bp=4
        )
        intron_sets = {c.introns for c in out["m"]}
        assert (wobble,) in intron_sets
        assert (canon,) in intron_sets

    def test_missing_query_name_skipped(self):
        out = expand_canonical_chain_alternatives_v2(
            [({"query_name": None}, IntronChain(introns=(WOBBLE,)))],
            GENOME, "+", MOTIFS_EXTENDED, search_bp=4,
        )
        assert out == {}
