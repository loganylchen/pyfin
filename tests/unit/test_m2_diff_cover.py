"""Unit tests for the diff-region coverage gate helpers.

``wobble_diff_spans`` emits a genomic span bracketing each wobbling junction in a
tie (donor->acceptor), ``read_straddles`` tests whether a read's events traverse
such a span, and ``event_genomic_positions`` projects eventalign positions back to
genomic. These back the m2_em diff-region coverage gate (config.m2_diff_cover_gate).
"""
from __future__ import annotations

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.scoring.m2_junction_nll import (
    event_genomic_positions,
    read_straddles,
    structural_wobble_clusters,
    wobble_diff_spans,
)


def _cand(introns, start=100, end=700, strand="+"):
    seq = "ACGT" * 200
    return TranscriptCandidate(
        candidate_id="c",
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=end,
        sequence=seq,
        source="novel",
        supporting_read_ids=set(),
        chrom="chr1",
        strand=strand,
        start=start,
        end=end,
    )


class TestWobbleDiffSpans:
    def test_donor_wobble(self):
        # introns differ at the donor (200 vs 203); same acceptor 300.
        a = _cand([(200, 300)])
        b = _cand([(203, 300)])
        assert wobble_diff_spans([a, b], flank=2) == [(198, 302)]

    def test_acceptor_wobble(self):
        a = _cand([(200, 300)])
        b = _cand([(200, 297)])
        assert wobble_diff_spans([a, b], flank=2) == [(198, 302)]

    def test_multi_junction_only_wobbling_emitted(self):
        # junction1 wobbles (200/203), junction2 unanimous (500,600) -> skipped.
        a = _cand([(200, 300), (500, 600)])
        b = _cand([(203, 300), (500, 600)])
        assert wobble_diff_spans([a, b], flank=2) == [(198, 302)]

    def test_unanimous_is_empty(self):
        a = _cand([(200, 300)])
        b = _cand([(200, 300)])
        assert wobble_diff_spans([a, b], flank=2) == []

    def test_exon_skip_counts_as_wobble(self):
        # b lacks the second intron entirely -> that cluster is not unanimous.
        a = _cand([(200, 300), (500, 600)])
        b = _cand([(200, 300)])
        assert wobble_diff_spans([a, b], flank=2) == [(498, 602)]

    def test_single_candidate_is_empty(self):
        assert wobble_diff_spans([_cand([(200, 300)])], flank=2) == []


class TestReadStraddles:
    def test_both_sides_straddles(self):
        assert read_straddles([195, 305], 198, 302) is True

    def test_acceptor_side_only_does_not(self):
        assert read_straddles([305, 308], 198, 302) is False

    def test_donor_side_only_does_not(self):
        assert read_straddles([195, 199], 198, 302) is False

    def test_empty_events_do_not(self):
        assert read_straddles([], 198, 302) is False


class TestEventGenomicPositions:
    def test_single_exon_projection(self):
        # single-exon (+) candidate: tx offset p -> genomic start+p.
        c = _cand([], start=100, end=200)
        res = {"position": [0, 5, 10]}
        assert event_genomic_positions(res, c) == [100, 105, 110]

    def test_missing_position_key(self):
        c = _cand([], start=100, end=200)
        assert event_genomic_positions({}, c) == []


class TestStructuralWobbleClusters:
    """Structural clustering: same intron count + every junction within +-bp,
    all sources; the basis for the post-EM abundance-ratio shadow drop."""

    def _members(self, clusters):
        return sorted(sorted(v) for v in clusters.values() if len(v) > 1)

    def test_wobble_within_bp_clusters(self):
        # three +-2bp donor wobbles of one true junction -> one cluster.
        cs = [
            _cand([(200, 300)]), _cand([(202, 300)]), _cand([(198, 300)]),
        ]
        cl = structural_wobble_clusters(cs, bp=6)
        assert self._members(cl) == [[0, 1, 2]]

    def test_beyond_bp_separate(self):
        cs = [_cand([(200, 300)]), _cand([(215, 300)])]  # 15bp > 6
        assert self._members(structural_wobble_clusters(cs, bp=6)) == []

    def test_different_intron_count_separate(self):
        cs = [_cand([(200, 300)]), _cand([(200, 300), (400, 500)])]
        assert self._members(structural_wobble_clusters(cs, bp=6)) == []

    def test_opposite_strand_separate(self):
        cs = [_cand([(200, 300)], strand="+"), _cand([(200, 300)], strand="-")]
        assert self._members(structural_wobble_clusters(cs, bp=6)) == []

    def test_mono_exon_never_clustered(self):
        cs = [_cand([], start=0, end=500), _cand([], start=1, end=500)]
        assert self._members(structural_wobble_clusters(cs, bp=6)) == []

    def test_gtf_and_novel_cluster_together(self):
        # a GTF anchor + its +-bp novel shadow land in one cluster (source-agnostic).
        gtf = _cand([(200, 300)])
        gtf.source = "gtf"
        nov = _cand([(203, 300)])  # default source novel
        cl = structural_wobble_clusters([gtf, nov], bp=6)
        assert self._members(cl) == [[0, 1]]
