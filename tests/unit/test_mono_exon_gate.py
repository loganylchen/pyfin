"""Unit tests for Lever 3 — mono-exon (single-exon) read-support gate.

``mono_exon_drops`` drops a NOVEL single-exon candidate whose hard read count
< min_reads OR genomic length < min_len. Multi-exon, GTF and fusion are exempt;
both thresholds 0 disables (byte-identical).
"""
from __future__ import annotations

from fin.analysis.quantification import QuantResult, mono_exon_drops


def _qr(cid, *, source="novel", reads=10, start=0, end=1000, n_exons=1):
    exons = ((start, end),) if n_exons == 1 else tuple(
        (start + i * 100, start + i * 100 + 50) for i in range(n_exons)
    )
    return QuantResult(
        candidate_id=cid, abundance=float(reads), confidence=1.0,
        num_assigned_reads=reads, source=source, chrom="chr1", strand="+",
        start=start, end=end, exons=exons,
    )


def _drop(results, **kw):
    return mono_exon_drops({q.candidate_id: q for q in results}, **kw)


class TestDisabled:
    def test_both_thresholds_zero_is_noop(self):
        r = [_qr("A", reads=1, start=0, end=10)]  # would fail any gate
        assert _drop(r, min_reads=0, min_len=0) == set()


class TestReadSupport:
    def test_low_read_mono_dropped(self):
        r = [_qr("lo", reads=2), _qr("hi", reads=20)]
        assert _drop(r, min_reads=5, min_len=0) == {"lo"}

    def test_exactly_threshold_kept(self):
        # >= min_reads survives (drop only when strictly below)
        r = [_qr("eq", reads=5)]
        assert _drop(r, min_reads=5, min_len=0) == set()


class TestLength:
    def test_short_mono_dropped(self):
        r = [_qr("short", start=0, end=120), _qr("long", start=0, end=5000)]
        assert _drop(r, min_reads=0, min_len=200) == {"short"}

    def test_either_threshold_fails_drops(self):
        # long but low-read -> dropped by reads; short but high-read -> by len
        r = [_qr("lowread", reads=1, start=0, end=5000),
             _qr("shortlen", reads=50, start=0, end=100)]
        assert _drop(r, min_reads=5, min_len=200) == {"lowread", "shortlen"}


class TestExemptions:
    def test_multi_exon_novel_exempt(self):
        r = [_qr("multi", reads=1, start=0, end=80, n_exons=3)]
        assert _drop(r, min_reads=5, min_len=200) == set()

    def test_gtf_mono_exempt(self):
        r = [_qr("g", source="gtf", reads=1, start=0, end=50)]
        assert _drop(r, min_reads=5, min_len=200) == set()

    def test_fusion_mono_exempt(self):
        r = [_qr("f", source="fusion", reads=1, start=0, end=50)]
        assert _drop(r, min_reads=5, min_len=200) == set()

    def test_high_support_long_real_mono_kept(self):
        # a real intronless gene: many reads, long -> survives
        r = [_qr("histone", reads=200, start=0, end=900)]
        assert _drop(r, min_reads=5, min_len=200) == set()
