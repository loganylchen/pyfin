"""Catalog completion: opt-in emission of unreached annotation entries."""
import dataclasses

import pytest

from fin.analysis.catalog_completion import (
    ADAPTIVE_CATALOG_FRACTION,
    CATALOG_SOURCE,
    complete_from_catalog,
    structure_key,
)
from fin.analysis.quantification import QuantResult
from fin.pipeline.config import PipelineConfig


def _qr(cid, exons, *, chrom="chr1", strand="+", source="novel", ab=5.0):
    exons = tuple(exons)
    return QuantResult(
        candidate_id=cid, abundance=ab, confidence=1.0,
        num_assigned_reads=int(ab), source=source, chrom=chrom, strand=strand,
        start=exons[0][0], end=exons[-1][1], exons=exons,
    )


class _Tx:
    def __init__(self, tid, gid, chrom, strand, exons):
        self.transcript_id = tid
        self.gene_id = gid
        self.chrom = chrom
        self.strand = strand
        self.exons = list(exons)


class _Reader:
    def __init__(self, txs):
        self.transcripts = {t.transcript_id: t for t in txs}


# --------------------------------------------------------------- structure key
def test_structure_key_separates_chrom_and_strand():
    """Identical intron patterns on different contigs/strands are different."""
    ex = ((100, 200), (300, 400))
    a = structure_key("chr1", "+", ex)
    b = structure_key("chr2", "+", ex)
    c = structure_key("chr1", "-", ex)
    assert a != b and a != c and b != c


def test_structure_key_is_intron_based_for_multi_exon():
    """Terminal coordinates may move without changing identity."""
    a = structure_key("chr1", "+", ((100, 200), (300, 400)))
    b = structure_key("chr1", "+", ((150, 200), (300, 380)))
    assert a == b


def test_structure_key_mono_exon_uses_exact_span():
    a = structure_key("chr1", "+", ((100, 200),))
    b = structure_key("chr1", "+", ((100, 201),))
    assert a != b
    assert a[2] == "MONO"


# ------------------------------------------------------------------- off mode
def test_off_mode_returns_the_same_object_untouched():
    agg = {"x": _qr("x", [(100, 200), (300, 400)])}
    out, audit = complete_from_catalog(agg, None, mode="off")
    assert out is agg
    assert audit["added"] == 0
    assert audit["catalog_mode"] == "off"


# -------------------------------------------------------------- complete mode
def test_complete_adds_only_structures_not_already_emitted():
    agg = {"asm": _qr("asm", [(100, 200), (300, 400)])}
    reader = _Reader([
        # same intron chain as the assembled model -> must NOT be added
        _Tx("same", "g1", "chr1", "+", [(90, 200), (300, 410)]),
        # different chain -> must be added
        _Tx("novelcat", "g2", "chr1", "+", [(1000, 1100), (1200, 1300)]),
    ])
    out, audit = complete_from_catalog(agg, reader, mode="complete")
    assert audit["already_assembled"] == 1
    assert audit["added"] == 1
    assert "novelcat" in out
    assert "same" not in out
    assert len(out) == 2


def test_added_entries_are_zero_abundance_and_marked_catalog():
    agg = {"asm": _qr("asm", [(100, 200), (300, 400)])}
    reader = _Reader([_Tx("cat", "g", "chr1", "+", [(9000, 9100), (9200, 9300)])])
    out, _ = complete_from_catalog(agg, reader, mode="complete")
    added = out["cat"]
    assert added.source == CATALOG_SOURCE
    assert added.abundance == 0.0
    assert added.num_assigned_reads == 0
    assert added.exons == ((9000, 9100), (9200, 9300))


def test_completion_never_mutates_or_drops_assembled_models():
    asm = _qr("asm", [(100, 200), (300, 400)], ab=7.5)
    agg = {"asm": asm}
    reader = _Reader([_Tx("cat", "g", "chr1", "+", [(9000, 9100), (9200, 9300)])])
    out, _ = complete_from_catalog(agg, reader, mode="complete")
    assert out["asm"] is asm
    assert out["asm"].abundance == 7.5
    assert set(agg) <= set(out)


def test_completion_is_deterministic():
    agg = {"asm": _qr("asm", [(100, 200), (300, 400)])}
    txs = [_Tx(f"t{i}", "g", "chr1", "+", [(i * 1000, i * 1000 + 50),
                                           (i * 1000 + 100, i * 1000 + 150)])
           for i in range(1, 8)]
    a, aa = complete_from_catalog(agg, _Reader(txs), mode="complete")
    b, ba = complete_from_catalog(agg, _Reader(list(reversed(txs))),
                                  mode="complete")
    assert list(a) == list(b)
    assert aa["added_ids"] == ba["added_ids"]


def test_adaptive_keeps_novel_models_when_catalog_is_sparse():
    # One catalog-derived structure versus three read-derived structures gives
    # fraction 0.25, safely below the frozen threshold. Even overlapping novel
    # models must remain available to recover annotation omissions.
    gtf = _qr("gtf", [(100, 200), (300, 400)], source="gtf")
    agg = {
        "gtf": gtf,
        "n1": _qr("n1", [(110, 210), (310, 410)]),
        "n2": _qr("n2", [(120, 220), (320, 420)]),
        "n3": _qr("n3", [(130, 230), (330, 430)]),
    }
    reader = _Reader([_Tx("cat", "g", "chr1", "+",
                              [(100, 200), (300, 400)])])
    out, audit = complete_from_catalog(agg, reader, mode="adaptive")
    assert set(out) == set(agg)
    assert audit["catalog_fraction"] == pytest.approx(0.25)
    assert audit["adaptive_suppression_applied"] is False
    assert audit["novel_suppressed"] == 0


def test_adaptive_suppresses_only_after_catalog_fraction_gate():
    # Nine catalog structures versus one novel gives fraction 0.90, above the
    # 0.87 threshold frozen before holdout evaluation.
    txs = [
        _Tx(f"cat{i}", "g", "chr1", "+",
            [(i * 1000 + 100, i * 1000 + 200),
             (i * 1000 + 300, i * 1000 + 400)])
        for i in range(9)
    ]
    agg = {
        "gtf": _qr("gtf", txs[0].exons, source="gtf"),
        # Overlaps cat1, which is absent before completion, so the transfer
        # target is a newly added source="catalog" entry.
        "artifact": _qr("artifact", [(1110, 1210), (1320, 1420)]),
    }
    out, audit = complete_from_catalog(agg, _Reader(txs), mode="adaptive")
    assert ADAPTIVE_CATALOG_FRACTION == 0.87
    assert audit["catalog_fraction"] == pytest.approx(0.9)
    assert audit["adaptive_suppression_applied"] is True
    assert audit["novel_suppressed"] == 1
    assert audit["novel_suppressed_effective"] == 1
    assert audit["transferred_models"] == 1
    assert audit["transfers_to_completion_entries"] == 1
    assert audit["transfers"][0]["to_source"] == CATALOG_SOURCE
    assert audit["restored_models"] == 0
    assert audit["mass_before"] == pytest.approx(10.0)
    assert audit["mass_after"] == pytest.approx(10.0)
    assert audit["mass_delta"] == pytest.approx(0.0, abs=1e-12)
    assert "artifact" not in out
    assert len(out) == 9
    assert sum(q.abundance for q in out.values()) == pytest.approx(10.0)


def test_mono_exon_catalog_entry_is_added_when_span_differs():
    agg = {"asm": _qr("asm", [(100, 200)])}
    reader = _Reader([_Tx("m", "g", "chr1", "+", [(100, 250)])])
    out, audit = complete_from_catalog(agg, reader, mode="complete")
    assert audit["added"] == 1 and "m" in out


def test_read_extent_sibling_covers_its_mono_parent():
    """A surviving ``X.ext`` (mono, span re-estimated from the reads) IS the
    annotated model X: completion must not add X's rejected span back."""
    sib = _qr("m1.ext", [(100, 260)], source="gtf")
    sib.extent_of = "m1"
    out, audit = complete_from_catalog(
        {"m1.ext": sib}, _Reader([_Tx("m1", "g", "chr1", "+", [(100, 200)])]), mode="complete")
    assert set(out) == {"m1.ext"}
    assert audit["already_assembled"] == 1 and audit["catalog_models_for_fraction"] == 1


def test_transcript_without_exons_is_skipped_not_crashed():
    agg = {"asm": _qr("asm", [(100, 200), (300, 400)])}
    reader = _Reader([_Tx("empty", "g", "chr1", "+", [])])
    out, audit = complete_from_catalog(agg, reader, mode="complete")
    assert audit["skipped_no_exons"] == 1
    assert audit["added"] == 0
    assert len(out) == 1


def test_id_collision_is_namespaced_not_overwritten():
    agg = {"dup": _qr("dup", [(100, 200), (300, 400)])}
    reader = _Reader([_Tx("dup", "g", "chr2", "+", [(10, 20), (30, 40)])])
    out, audit = complete_from_catalog(agg, reader, mode="complete")
    assert audit["added"] == 1
    assert out["dup"].chrom == "chr1"          # original untouched
    assert out["catalog_dup"].chrom == "chr2"  # completion namespaced


# ----------------------------------------------------------------- contract
def test_complete_without_annotation_fails_loudly():
    """An explicitly requested mode must not silently degrade."""
    with pytest.raises(ValueError, match="requires an annotation"):
        complete_from_catalog({}, None, mode="complete")


def test_unknown_mode_rejected():
    with pytest.raises(ValueError, match="unknown catalog_mode"):
        complete_from_catalog({}, None, mode="bogus")


def test_config_rejects_catalog_mode_without_gtf(tmp_path):
    # validate() checks that the inputs exist before anything else
    bam, genome = tmp_path / "input.bam", tmp_path / "genome.fa"
    bam.touch()
    genome.touch()
    cfg = PipelineConfig.from_profile(
        "sirv", bam_path=str(bam), genome_fasta_path=str(genome))
    bad = dataclasses.replace(cfg, catalog_mode="complete", gtf_path="")
    with pytest.raises(ValueError, match="catalog_mode requires --gtf"):
        bad.validate()


def test_no_profile_enables_catalog_mode():
    for name in ("custom", "real-drna", "real-drna-precision", "sirv"):
        cfg = PipelineConfig.from_profile(name, bam_path="x.bam")
        assert cfg.catalog_mode == "off", name
