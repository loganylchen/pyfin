"""Mechanical paired-run abundance-refit audit tests."""
from __future__ import annotations

import json

import pytest

from experiments.prod_validation.abundance_refit_audit import audit


HEADER = (
    "candidate_id\tgene_id\tchrom\tstrand\tstart\tend\tsource\t"
    "abundance\tconfidence\tcoherence_score\tdiscrimination_score\t"
    "combined_score\tnum_reads\ttpm\tmax_R\n"
)


def _write_run(path, *, abundance, chrom="chr1", diagnostics=False):
    path.mkdir(parents=True)
    (path / "assembly.gtf").write_text(
        f'{chrom}\tpyfin\tgene\t1\t100\t.\t+\t.\tgene_id "g";\n'
        f'{chrom}\tpyfin\ttranscript\t1\t100\t.\t+\t.\tgene_id "g"; '
        f'transcript_id "t"; abundance "{abundance:.4f}"; confidence "0.5000"; '
        'num_reads "1"; transcript_source "novel"; coherence_score "0.0000"; '
        'discrimination_score "0.0000"; combined_score "0.0000";\n'
        f'{chrom}\tpyfin\texon\t1\t100\t.\t+\t.\tgene_id "g"; '
        'transcript_id "t"; exon_number "1";\n'
    )
    (path / "scores.tsv").write_text(
        HEADER
        + f"t\tg\t{chrom}\t+\t0\t100\tnovel\t{abundance:.4f}\t0.5000\t"
        "0.0000\t0.0000\t0.0000\t1\t1000000.0000\t0.5000\n"
    )
    if diagnostics:
        (path / "abundance_refit.json").write_text(json.dumps({
            "effective": True,
            "mass_balance_error": 0.0,
            "max_per_read_conservation_error": 0.0,
            "selection_orphaned_reads": 0,
            "unassigned_mass": 0.0,
            "rescued_mass": 0.5,
        }))


def test_audit_accepts_numeric_only_changes_and_legacy_identity(tmp_path):
    on = tmp_path / "on"
    off = tmp_path / "off"
    legacy = tmp_path / "legacy"
    _write_run(on, abundance=1.0, diagnostics=True)
    _write_run(off, abundance=0.5)
    _write_run(legacy, abundance=0.5)
    before = (on / "abundance_refit.json").read_bytes()
    report = audit(on, off, legacy)
    assert report["gtf_structural_identity"] is True
    assert report["legacy_gtf_byte_identity"] is True
    assert report["legacy_tsv_byte_identity"] is True
    assert report["tsv_quant_columns_changed"]["abundance"] == 1
    assert report["structural_identity"] == "passed"

    # The audited artifact must never be rewritten: stamping it would make an
    # audited run differ byte-wise from an unaudited one.
    assert (on / "abundance_refit.json").read_bytes() == before
    verdict = json.loads((on / "refit_identity.json").read_text())
    assert verdict["structural_identity"] == "passed"


def test_audit_accepts_same_coordinate_gtf_record_reordering(tmp_path):
    on = tmp_path / "on"
    off = tmp_path / "off"
    _write_run(on, abundance=1.0, diagnostics=True)
    _write_run(off, abundance=0.5)

    on_block = on.joinpath("assembly.gtf").read_text().splitlines()
    off_block = off.joinpath("assembly.gtf").read_text().splitlines()
    second = [
        line.replace('gene_id "g"', 'gene_id "h"').replace(
            'transcript_id "t"', 'transcript_id "u"'
        )
        for line in off_block
    ]
    on.joinpath("assembly.gtf").write_text(
        "\n".join(second + on_block) + "\n"
    )
    off.joinpath("assembly.gtf").write_text(
        "\n".join(off_block + second) + "\n"
    )
    second_row = (
        "u\th\tchr1\t+\t0\t100\tnovel\t0.5000\t0.5000\t"
        "0.0000\t0.0000\t0.0000\t1\t1000000.0000\t0.5000\n"
    )
    for path in (on, off):
        path.joinpath("scores.tsv").write_text(
            path.joinpath("scores.tsv").read_text() + second_row
        )

    report = audit(on, off)
    assert report["gtf_record_multiset_identity"] is True
    assert report["gtf_line_order_identity"] is False
    assert report["gtf_line_order_differences"] > 0


def test_audit_rejects_structural_change(tmp_path):
    on = tmp_path / "on"
    off = tmp_path / "off"
    _write_run(on, abundance=1.0, chrom="chr2", diagnostics=True)
    _write_run(off, abundance=0.5, chrom="chr1")
    with pytest.raises(AssertionError, match="structurally"):
        audit(on, off)
