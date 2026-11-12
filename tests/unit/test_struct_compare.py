"""The regression gate's comparator must ignore minted-id spelling and nothing else."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_PATH = Path(__file__).resolve().parents[2] / "experiments" / "prod_validation" / "v1_struct_compare.py"
_spec = importlib.util.spec_from_file_location("v1_struct_compare", _PATH)
cmp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cmp)


def _gtf(tmp_path, name, transcripts, gene_spans=None):
    """transcripts: (gene_id, transcript_id, chrom, strand, exons, abundance[, exon_chrom])."""
    lines, genes = [], {}
    for gid, tid, chrom, strand, exons, abundance, *rest in transcripts:
        exon_chrom = rest[0] if rest else chrom
        start, end = exons[0][0], exons[-1][1]
        lo, hi = genes.get(gid, (chrom, strand, start, end))[2:]
        genes[gid] = (chrom, strand, min(lo, start), max(hi, end))
        source = "novel" if tid.startswith("novel_") else "gtf"
        lines.append(
            f'{chrom}\tpyfin\ttranscript\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gid}"; '
            f'transcript_id "{tid}"; abundance "{abundance}"; transcript_source "{source}";'
        )
        for n, (s, e) in enumerate(exons, 1):
            lines.append(
                f'{exon_chrom}\tpyfin\texon\t{s}\t{e}\t.\t{strand}\t.\tgene_id "{gid}"; '
                f'transcript_id "{tid}"; exon_number "{n}";'
            )
    for gid, (chrom, strand, start, end) in genes.items():
        lines.append(f'{chrom}\tpyfin\tgene\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gid}";')
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return str(path)


A, B, C, D = ([(100, 200), (300, 400)], [(100, 200), (310, 400)],
              [(100, 200), (320, 400)], [(100, 200), (330, 400)])


def _same(tmp_path, left, right):
    return cmp.records(_gtf(tmp_path, "l.gtf", left)) == cmp.records(_gtf(tmp_path, "r.gtf", right))


def test_respelling_minted_ids_is_not_a_difference(tmp_path):
    left = [("novel_aa", "novel_a1", "chr1", "+", A, "3.0"), ("novel_aa", "novel_a2", "chr1", "+", B, "2.0"),
            ("ENSG1", "ENST1", "chr1", "+", C, "5.0"), ("ENSG1", "novel_a3", "chr1", "+", D, "1.0")]
    right = [("novel_ff", "novel_09", "chr1", "+", A, "3.0"), ("novel_ff", "novel_08", "chr1", "+", B, "2.0"),
             ("ENSG1", "ENST1", "chr1", "+", C, "5.0"), ("ENSG1", "novel_07", "chr1", "+", D, "1.0")]
    assert _same(tmp_path, left, right)


def test_regrouping_the_same_transcripts_into_other_novel_genes_is_a_difference(tmp_path):
    left = [("novel_01", "novel_a1", "chr1", "+", A, "1.0"), ("novel_01", "novel_a2", "chr1", "+", B, "1.0"),
            ("novel_02", "novel_a3", "chr1", "+", C, "1.0"), ("novel_02", "novel_a4", "chr1", "+", D, "1.0")]
    right = [("novel_01", "novel_a1", "chr1", "+", A, "1.0"), ("novel_01", "novel_a3", "chr1", "+", C, "1.0"),
             ("novel_02", "novel_a2", "chr1", "+", B, "1.0"), ("novel_02", "novel_a4", "chr1", "+", D, "1.0")]
    assert not _same(tmp_path, left, right)


def test_moving_a_novel_transcript_into_an_annotated_gene_is_a_difference(tmp_path):
    left = [("ENSG1", "ENST1", "chr1", "+", A, "1.0"), ("novel_01", "novel_a1", "chr1", "+", B, "1.0")]
    right = [("ENSG1", "ENST1", "chr1", "+", A, "1.0"), ("ENSG1", "novel_a1", "chr1", "+", B, "1.0")]
    assert not _same(tmp_path, left, right)


@pytest.mark.parametrize("field, value", [
    (2, "chr2"), (3, "-"), (4, [(100, 200), (301, 400)]), (5, "3.5"), (1, "ENST2"),
])
def test_any_structural_or_attribute_change_is_a_difference(tmp_path, field, value):
    left = [("ENSG1", "ENST1", "chr1", "+", A, "3.0")]
    changed = list(left[0])
    changed[field] = value
    assert not _same(tmp_path, left, [tuple(changed)])


def test_an_exon_on_another_chromosome_is_a_difference(tmp_path):
    left = [("ENSG1", "ENST1", "chr1", "+", A, "3.0")]
    right = [("ENSG1", "ENST1", "chr1", "+", A, "3.0", "chr2")]
    assert not _same(tmp_path, left, right)


@pytest.mark.parametrize("old, new", [
    ('exon_number "2";', 'exon_number "7";'),
    ('\tgene\t100\t400\t.\t+\t.\tgene_id "ENSG1";', '\tgene\t100\t400\t.\t+\t.\tgene_id "ENSG1"; gene_name "X";'),
    ('\tgene\t100\t400\t', '\tgene\t100\t401\t'),
])
def test_gene_and_exon_metadata_are_compared(tmp_path, old, new):
    left = _gtf(tmp_path, "l.gtf", [("ENSG1", "ENST1", "chr1", "+", A, "3.0")])
    text = Path(left).read_text()
    assert old in text
    right = tmp_path / "r.gtf"
    right.write_text(text.replace(old, new))
    assert cmp.records(left) != cmp.records(str(right))


def test_a_duplicated_model_is_a_difference(tmp_path):
    one = [("novel_01", "novel_a1", "chr1", "+", A, "1.0")]
    two = one + [("novel_02", "novel_a2", "chr1", "+", A, "1.0")]
    assert not _same(tmp_path, one, two)


def test_malformed_input_is_reported_not_compared(tmp_path, capsys, monkeypatch):
    good = _gtf(tmp_path, "good.gtf", [("ENSG1", "ENST1", "chr1", "+", A, "3.0")])
    bad = tmp_path / "bad.gtf"
    bad.write_text(Path(good).read_text().replace('exon\t100\t200\t.\t+\t.\tgene_id "ENSG1"',
                                                  'exon\t100\t200\t.\t+\t.\tgene_id "ENSG9"'))
    monkeypatch.setattr("sys.argv", ["cmp", good, str(bad)])
    assert cmp.main() == 2
    assert "MALFORMED" in capsys.readouterr().out
