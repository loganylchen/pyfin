"""Tests for benchmarks/em_subset_evaluate.py chain matching."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def _load_evaluator():
    """Load em_subset_evaluate.py as a module (it lives under benchmarks/)."""
    p = Path(__file__).resolve().parents[2] / "benchmarks" / "em_subset_evaluate.py"
    spec = importlib.util.spec_from_file_location("em_subset_evaluate", p)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["em_subset_evaluate"] = mod  # so @dataclass can find __module__
    spec.loader.exec_module(mod)
    return mod


ev = _load_evaluator()


def _write_gtf(path: Path, records: list[tuple]) -> None:
    """records: list of (chrom, feat, start_1based, end_1based, strand, tid).

    Single-exon writes both transcript + exon. Multi-exon callers handle order.
    """
    with open(path, "w") as f:
        for chrom, feat, s1, e1, strand, tid in records:
            attrs = f'transcript_id "{tid}"; gene_id "{tid}";'
            f.write(f"{chrom}\tt\t{feat}\t{s1}\t{e1}\t.\t{strand}\t.\t{attrs}\n")


class TestGtfParsing:
    def test_single_exon(self, tmp_path):
        gtf = tmp_path / "x.gtf"
        _write_gtf(gtf, [
            ("chr1", "exon", 11, 20, "+", "t1"),
        ])
        txs = ev._load_gtf(str(gtf))
        assert "t1" in txs
        # GTF 1-based inclusive 11..20 → 0-based half-open (10, 20)
        assert txs["t1"].exons == [(10, 20)]
        assert txs["t1"].intron_chain == ()
        assert txs["t1"].span == (10, 20)

    def test_multi_exon_intron_chain(self, tmp_path):
        gtf = tmp_path / "x.gtf"
        _write_gtf(gtf, [
            ("chr1", "exon", 1, 100, "+", "t2"),
            ("chr1", "exon", 201, 300, "+", "t2"),
            ("chr1", "exon", 401, 500, "+", "t2"),
        ])
        txs = ev._load_gtf(str(gtf))
        # introns are (e1_end, e2_start) in 0-based half-open
        assert txs["t2"].intron_chain == ((100, 200), (300, 400))


class TestPerfectMatch:
    def test_identical_gtfs_perfect_score(self, tmp_path):
        ref = tmp_path / "ref.gtf"
        _write_gtf(ref, [
            ("chr1", "exon", 1, 100, "+", "g1"),
            ("chr1", "exon", 201, 300, "+", "g1"),
        ])
        gt_txs = ev._load_gtf(str(ref))
        gt_multi, gt_mono = ev._build_chain_keys(gt_txs)
        cr = ev._evaluate_gtf(str(ref), gt_multi, gt_mono, "m2+m3", "with-gtf")
        assert cr.tp == 1
        assert cr.fp == 0
        assert cr.fn == 0
        assert cr.precision == 1.0
        assert cr.recall == 1.0
        assert cr.f1 == 1.0


class TestMixedMatch:
    def test_partial_overlap(self, tmp_path):
        ref = tmp_path / "ref.gtf"
        pred = tmp_path / "pred.gtf"
        _write_gtf(ref, [
            ("chr1", "exon", 1, 100, "+", "g1"),
            ("chr1", "exon", 201, 300, "+", "g1"),
            ("chr1", "exon", 1, 100, "+", "g2"),
            ("chr1", "exon", 401, 500, "+", "g2"),
        ])
        # pred recovers g1 but invents an off-by-10 chain for g2
        _write_gtf(pred, [
            ("chr1", "exon", 1, 100, "+", "p1"),
            ("chr1", "exon", 201, 300, "+", "p1"),
            ("chr1", "exon", 1, 100, "+", "p2"),
            ("chr1", "exon", 411, 500, "+", "p2"),  # shifted donor
        ])
        gt_multi, gt_mono = ev._build_chain_keys(ev._load_gtf(str(ref)))
        cr = ev._evaluate_gtf(str(pred), gt_multi, gt_mono, "all", "no-gtf")
        assert cr.tp == 1
        assert cr.fp == 1
        assert cr.fn == 1
        assert cr.precision == 0.5
        assert cr.recall == 0.5


class TestStrandSensitive:
    def test_opposite_strand_is_fp(self, tmp_path):
        ref = tmp_path / "ref.gtf"
        pred = tmp_path / "pred.gtf"
        _write_gtf(ref, [
            ("chr1", "exon", 1, 100, "+", "g1"),
            ("chr1", "exon", 201, 300, "+", "g1"),
        ])
        _write_gtf(pred, [
            ("chr1", "exon", 1, 100, "-", "p1"),
            ("chr1", "exon", 201, 300, "-", "p1"),
        ])
        gt_multi, gt_mono = ev._build_chain_keys(ev._load_gtf(str(ref)))
        cr = ev._evaluate_gtf(str(pred), gt_multi, gt_mono, "m1", "with-gtf")
        # Same intron chain on the opposite strand should NOT match.
        assert cr.tp == 0
        assert cr.fp == 1
        assert cr.fn == 1
