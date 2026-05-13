"""Tests for ``snap_junctions`` in fin.candidates.intron_chains.

The snap step clusters donor (intron start, + strand) and acceptor (intron
end, + strand) positions across reads in an interval, then rewrites each
read's chain to per-cluster consensus positions. Canonical GT-AG motifs
are preferred when available.
"""

from __future__ import annotations

from fin.candidates.dataclasses import IntronChain
from fin.candidates.intron_chains import (
    _cluster_positions,
    _has_canonical_motif,
    snap_junctions,
)


def _make(name: str, introns):
    """Build a minimal (read_dict, IntronChain) tuple."""
    return ({"query_name": name}, IntronChain(introns=tuple(introns)))


def test_empty_input_returns_empty():
    assert snap_junctions([], genome_seq="", strand="+", snap_bp=6) == []


def test_snap_disabled_when_snap_bp_zero():
    rc = [_make("r1", [(100, 200)])]
    out = snap_junctions(rc, genome_seq="A" * 300, strand="+",
                          snap_bp=0, min_support=2)
    assert out[0][1].introns == ((100, 200),)


def test_two_reads_no_motif_snap_to_mode():
    # donors at 100 / 106 (within 6 bp), acceptors at 200 / 200.
    # Genome has no GT/AG motifs anywhere -> fall back to mode/median.
    genome = "N" * 300
    rc = [
        _make("r1", [(100, 200)]),
        _make("r2", [(106, 200)]),
    ]
    out = snap_junctions(rc, genome_seq=genome, strand="+",
                          snap_bp=6, min_support=2)
    chains = [c.introns for _, c in out]
    # Both donors snap to the same position; mode tie -> median (n//2 index,
    # matching the convention used by 3'-end consensus elsewhere).
    assert chains[0] == chains[1]
    assert chains[0][0][0] in (100, 106)
    assert chains[0][0][1] == 200


def test_canonical_motif_wins_over_mode():
    # Donors at 100,100,103. Mode = 100. But genome has GT at position 103,
    # so canonical motif should snap all three to 103.
    genome = bytearray(b"N" * 300)
    genome[103] = ord("G")
    genome[104] = ord("T")
    # Acceptor motif AG at intron end e=200 -> genome[198:200]="AG"
    genome[198] = ord("A")
    genome[199] = ord("G")
    rc = [
        _make("r1", [(100, 200)]),
        _make("r2", [(100, 200)]),
        _make("r3", [(103, 200)]),
    ]
    out = snap_junctions(rc, genome_seq=genome.decode(), strand="+",
                          snap_bp=6, min_support=2)
    chains = [c.introns for _, c in out]
    # All three donors should snap to 103 because of the GT motif.
    for ch in chains:
        assert ch[0][0] == 103
        assert ch[0][1] == 200


def test_single_read_below_min_support_not_snapped():
    genome = "N" * 300
    rc = [_make("r1", [(100, 200)])]
    out = snap_junctions(rc, genome_seq=genome, strand="+",
                          snap_bp=6, min_support=2)
    # Cluster size 1 < min_support 2 -> unchanged.
    assert out[0][1].introns == ((100, 200),)


def test_degenerate_intron_dropped_after_snap():
    # Two introns close together. Snap collapses donor of intron2 onto its
    # own acceptor (or near it), creating s >= e -> dropped.
    genome = "N" * 500
    rc = [
        _make("r1", [(100, 150), (152, 200)]),
        _make("r2", [(100, 150), (152, 200)]),
        # An extra read with donor at 148 inside intron2's donor cluster
        # forces the donor consensus near 150 == intron2's acceptor.
        _make("r3", [(100, 150), (148, 200)]),
        _make("r4", [(100, 150), (150, 200)]),
    ]
    out = snap_junctions(rc, genome_seq=genome, strand="+",
                          snap_bp=6, min_support=2)
    for _, chain in out:
        for s, e in chain.introns:
            assert e > s, "degenerate intron should have been dropped"


def test_minus_strand_motif_check():
    # On - strand: donor = AC at genome[e-2:e], acceptor = CT at genome[s:s+2].
    # Build a genome with AC at intron end 200 and CT at intron start 100.
    g = bytearray(b"N" * 400)
    # Acceptor (- strand) CT at start positions 100, 100, 103: place CT at 103.
    g[103] = ord("C")
    g[104] = ord("T")
    # Donor (- strand) AC at end positions 200, 200, 197: place AC at 197 (=> e-2:e for e=197).
    g[195] = ord("A")
    g[196] = ord("C")
    rc = [
        _make("r1", [(100, 200)]),
        _make("r2", [(100, 200)]),
        _make("r3", [(103, 197)]),
    ]
    out = snap_junctions(rc, genome_seq=g.decode(), strand="-",
                          snap_bp=6, min_support=2)
    chains = [c.introns for _, c in out]
    # Acceptor (s on - strand) should snap to 103 (CT motif).
    # Donor (e on - strand) should snap to 197 (AC motif).
    for ch in chains:
        s, e = ch[0]
        assert s == 103
        assert e == 197


def test_cluster_positions_basic():
    # 3 cluster: [100,103,106], [120], [200,201]
    clusters = _cluster_positions([100, 106, 103, 120, 200, 201], gap=6)
    assert [sorted(c) for c in clusters] == [[100, 103, 106], [120], [200, 201]]


def test_cluster_positions_empty():
    assert _cluster_positions([], gap=6) == []


def test_has_canonical_motif_plus_strand():
    # genome[100:102] = "GT" -> donor at pos 100 on + strand is canonical.
    g = "N" * 100 + "GT" + "N" * 98 + "AG" + "N" * 100
    assert _has_canonical_motif(100, "donor", "+", g) is True
    # Acceptor on + strand uses genome[e-2:e]; e=202 -> genome[200:202]="AG".
    assert _has_canonical_motif(202, "acceptor", "+", g) is True
    # Wrong position is not canonical.
    assert _has_canonical_motif(50, "donor", "+", g) is False


def test_has_canonical_motif_out_of_bounds():
    assert _has_canonical_motif(-5, "donor", "+", "ACGT") is False
    assert _has_canonical_motif(100, "donor", "+", "ACGT") is False
