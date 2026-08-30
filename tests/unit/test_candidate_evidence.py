"""Observable-feature evidence table: correctness and missing-data sentinels."""
from collections import Counter

from fin.analysis.candidate_evidence import (
    EVIDENCE_COLUMNS,
    compute_candidate_evidence,
    write_candidate_evidence,
)
from fin.analysis.quantification import QuantResult
from fin.candidates.canonical import parse_motifs


def _qr(cid, exons, *, abundance, reads, family=None, strand="+",
        chrom="chr1", assigned=(), source="novel"):
    exons = tuple(exons)
    return QuantResult(
        candidate_id=cid, abundance=float(abundance), confidence=0.9,
        num_assigned_reads=reads, source=source, chrom=chrom, strand=strand,
        start=exons[0][0], end=exons[-1][1], exons=exons, family_id=family,
        assigned_read_ids=tuple(assigned),
    )


def _rows_by_id(rows):
    return {r["candidate_id"]: r for r in rows}


def test_family_share_rank_and_containment_geometry():
    parent = _qr("parent", [(0, 100), (200, 300), (400, 500)],
                 abundance=30.0, reads=30, family="F")
    child = _qr("child", [(210, 300), (400, 480)],
                abundance=10.0, reads=10, family="F")
    other = _qr("other", [(1000, 1100)], abundance=5.0, reads=5)
    rows = _rows_by_id(compute_candidate_evidence(
        {"parent": parent, "child": child, "other": other}))

    assert rows["parent"]["family_size"] == 2
    assert rows["parent"]["family_share"] == 0.75
    assert rows["parent"]["family_rank"] == 1
    assert rows["child"]["family_rank"] == 2
    assert rows["child"]["family_dominant_share"] == 0.75
    # child's chain (300,400) is a strict subchain of parent's ((100,200),(300,400))
    assert rows["child"]["is_subchain_of_sibling"] == 1
    assert rows["parent"]["is_superchain_of_sibling"] == 1
    assert rows["parent"]["is_subchain_of_sibling"] == 0
    # solo mono candidate: own family, share 1
    assert rows["other"]["family_size"] == 1
    assert rows["other"]["family_share"] == 1.0
    assert rows["other"]["is_mono"] == 1


def test_junction_support_features_and_missing_sentinels():
    qr = _qr("x", [(0, 100), (200, 300), (400, 500)], abundance=8.0, reads=8)
    support = {("chr1", "+"): Counter({(100, 200): 10, (300, 400): 1})}
    row = _rows_by_id(compute_candidate_evidence(
        {"x": qr}, junction_support=support))["x"]
    assert row["weakest_junction_support"] == 1.0
    assert row["median_junction_support"] == 10.0
    assert row["n_junctions_below3"] == 1

    missing = _rows_by_id(compute_candidate_evidence({"x": qr}))["x"]
    assert missing["weakest_junction_support"] == -1.0
    assert missing["n_junctions_below3"] == -1

    mono = _qr("m", [(0, 100)], abundance=5.0, reads=5)
    mono_row = _rows_by_id(compute_candidate_evidence(
        {"m": mono}, junction_support=support))["m"]
    assert mono_row["weakest_junction_support"] == -1.0


def test_canonical_fraction_uses_strand_aware_motifs():
    # + strand intron (4, 16): donor GT at 4..6, acceptor AG at 14..16
    genome = {"chr1": "AAAA" + "GT" + "TTTTTTTT" + "AG" + "AAAA"}
    qr = _qr("c", [(0, 4), (16, 20)], abundance=4.0, reads=4)
    row = _rows_by_id(compute_candidate_evidence(
        {"c": qr}, genome=genome, canonical_motifs=parse_motifs(["GT-AG"])))["c"]
    assert row["canonical_fraction"] == 1.0

    bad = {"chr1": "AAAA" + "CC" + "TTTTTTTT" + "CC" + "AAAA"}
    row2 = _rows_by_id(compute_candidate_evidence(
        {"c": qr}, genome=bad, canonical_motifs=parse_motifs(["GT-AG"])))["c"]
    assert row2["canonical_fraction"] == 0.0

    row3 = _rows_by_id(compute_candidate_evidence({"c": qr}))["c"]
    assert row3["canonical_fraction"] == -1.0


def test_end_support_fractions_are_strand_aware():
    qr = _qr("e", [(1000, 1200), (1400, 1600)], abundance=3.0, reads=3,
             strand="-", assigned=("r1", "r2", "r3"))
    read_ends = {
        "r1": (1000, 1600),   # both ends exact
        "r2": (1450, 1600),   # 5' (=end for '-') ok, 3' (=start) far
        "r3": (1000, 1300),   # 3' ok, 5' far
    }
    row = _rows_by_id(compute_candidate_evidence(
        {"e": qr}, read_ends=read_ends, end_window_bp=25))["e"]
    assert row["n_end_reads"] == 3
    assert abs(row["end5_support_frac"] - 2 / 3) < 1e-3
    assert abs(row["end3_support_frac"] - 2 / 3) < 1e-3
    assert abs(row["fulllen_frac"] - 1 / 3) < 1e-3

    no_spans = _rows_by_id(compute_candidate_evidence({"e": qr}))["e"]
    assert no_spans["end5_support_frac"] == -1.0


def test_write_is_atomic_and_schema_stable(tmp_path):
    qr = _qr("w", [(0, 100)], abundance=2.0, reads=2)
    rows = compute_candidate_evidence({"w": qr})
    out = tmp_path / "candidate_evidence.tsv"
    write_candidate_evidence(out, rows)
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == list(EVIDENCE_COLUMNS)
    assert len(lines) == 2
    assert not (tmp_path / "candidate_evidence.tsv.tmp").exists()
