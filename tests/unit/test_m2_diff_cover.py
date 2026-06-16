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
