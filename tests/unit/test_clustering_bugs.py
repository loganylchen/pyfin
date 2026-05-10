"""Regression tests for fin/analysis/clustering.py bug fixes.

Covers the novel-detection inner loop (formerly using a stale `ref_id` from
the outer loop) in `ThreePrimePositionClustering.prepare_cluster_data`.
"""

from __future__ import annotations

from fin.analysis.clustering import ThreePrimePositionClustering


def _make_clusterer():
    """Build a clusterer instance without opening BAM/FASTA."""
    return ThreePrimePositionClustering.__new__(ThreePrimePositionClustering)


def test_novel_detection_collapses_contained_read_plus_strand():
    """A shorter novel read fully contained in a longer novel read should be
    marked contained. Tests the post-fix inner loop variable scoping.
    """
    clusterer = _make_clusterer()

    # Two novel reads on '+' strand: short read is a 3'-truncated version
    # of long read (long extends 5 bp further at the 3' end).
    read_seqs = {
        "long_read": "ACGTACGTACGT",   # length 12
        "short_read": "ACGTACG",       # length 7, prefix of long_read
    }
    ref_seqs = {}  # No GTF refs -> falls into novel detection loop

    ids = ["long_read", "short_read"]
    # 3' end positions on '+' strand: long_read extends past short_read
    three_prime_positions = [112, 107]

    result = clusterer.prepare_cluster_data(
        ids=ids,
        three_prime_positions=three_prime_positions,
        read_seqs=read_seqs,
        ref_seqs=ref_seqs,
        strand="+",
    )

    # short_read is contained in long_read -> only long_read survives as ref
    assert "long_read" in result.ref_seqs
    assert "short_read" not in result.ref_seqs


def test_novel_detection_keeps_distinct_reads():
    """Two unrelated novel reads must both survive."""
    clusterer = _make_clusterer()

    read_seqs = {
        "read_a": "AAAAAAAAAAAAA",
        "read_b": "GGGGGGGGGGGGG",
    }
    ref_seqs = {}
    ids = ["read_a", "read_b"]
    three_prime_positions = [100, 100]

    result = clusterer.prepare_cluster_data(
        ids=ids,
        three_prime_positions=three_prime_positions,
        read_seqs=read_seqs,
        ref_seqs=ref_seqs,
        strand="+",
    )

    assert "read_a" in result.ref_seqs
    assert "read_b" in result.ref_seqs


def test_novel_detection_minus_strand():
    """Same containment logic should hold on '-' strand with inverted end_dif."""
    clusterer = _make_clusterer()

    read_seqs = {
        "long_read": "ACGTACGTACGT",
        "short_read": "ACGTACG",
    }
    ref_seqs = {}
    ids = ["long_read", "short_read"]
    # On '-' strand, larger position = upstream; long extends further 3'
    three_prime_positions = [100, 105]

    result = clusterer.prepare_cluster_data(
        ids=ids,
        three_prime_positions=three_prime_positions,
        read_seqs=read_seqs,
        ref_seqs=ref_seqs,
        strand="-",
    )

    assert "long_read" in result.ref_seqs
    assert "short_read" not in result.ref_seqs
