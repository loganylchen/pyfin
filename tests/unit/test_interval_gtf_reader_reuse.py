"""Lending the parsed annotation to interval generation must change nothing.

The runner parses the annotation in `setup()` and then called
`generate_intervals_from_gtf(path)`, which parsed it AGAIN -- twice per process.
Interval workers spawn rather than fork, so each worker repeated both parses: a
cell with `--threads 2` parsed a 1.5 GB GENCODE six times. Measured on a paired
probe whose only variable was annotation size (40,163 vs 252,894 transcripts on
an identical BAM and identical loci), parsing was 235 s of a 359 s excess.

Interval order is load-bearing beyond the intervals themselves: workers receive
integer INDICES into a list they rebuild for themselves, so an order that
differed between parent and worker would silently mean different loci.

Gates below follow Codex's review: ownership and lifecycle, then direct
differential equivalence over every dataclass field in order.
"""
import dataclasses
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from fin.io.interval_manager import generate_intervals_from_gtf  # noqa: E402
from fin.io.io_gtf import GTFReader  # noqa: E402

GTF = """chr1\tsrc\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
chr1\tsrc\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "t1";
chr1\tsrc\texon\t150\t260\t.\t-\t.\tgene_id "g2"; transcript_id "t2";
chr19\tsrc\texon\t500\t900\t.\t-\t.\tgene_id "g3"; transcript_id "t3";
chr19\tsrc\texon\t950\t1200\t.\t+\t.\tgene_id "g4"; transcript_id "t4";
chr21\tsrc\texon\t10\t80\t.\t+\t.\tgene_id "g5"; transcript_id "t5";
"""


@pytest.fixture
def annotation(tmp_path):
    path = tmp_path / "annotation.gtf"
    path.write_text(GTF)
    return path


def _parsed(path):
    reader = GTFReader(str(path))
    reader.open()
    reader.parse()
    return reader


def _fields(intervals):
    return [dataclasses.astuple(i) for i in intervals]


def test_a_borrowed_reader_reproduces_the_parsed_result_field_for_field(annotation):
    from_path = generate_intervals_from_gtf(str(annotation))
    reader = _parsed(annotation)
    borrowed = generate_intervals_from_gtf(str(annotation), reader=reader)
    assert _fields(borrowed) == _fields(from_path)
    # order, not just membership: workers index into this list by integer
    assert [i.interval_id for i in borrowed] == [i.interval_id for i in from_path]


def test_the_borrowed_reader_is_not_closed_or_mutated(annotation):
    reader = _parsed(annotation)
    before = {tid: (tx.chrom, tx.start, tx.end, tx.strand)
              for tid, tx in reader.transcripts.items()}
    generate_intervals_from_gtf(str(annotation), reader=reader)
    assert reader._parsed is True
    assert reader._file is not None, "a borrowed reader must not be closed"
    assert {tid: (tx.chrom, tx.start, tx.end, tx.strand)
            for tid, tx in reader.transcripts.items()} == before
    # and it is still usable afterwards
    assert _fields(generate_intervals_from_gtf(str(annotation), reader=reader))


def test_a_reader_is_never_reparsed_because_parse_appends(annotation):
    """`GTFReader.parse` rewinds the file but does NOT clear its dictionaries,
    so calling it twice appends duplicate features. Borrowing must never call
    it. Demonstrated here so the reason the contract exists is on record."""
    reader = _parsed(annotation)
    exons_before = len(reader.transcripts["t1"].exons)
    generate_intervals_from_gtf(str(annotation), reader=reader)
    assert len(reader.transcripts["t1"].exons) == exons_before


def test_an_unparsed_reader_is_refused(annotation):
    reader = GTFReader(str(annotation))
    reader.open()
    with pytest.raises(ValueError, match="already be parsed"):
        generate_intervals_from_gtf(str(annotation), reader=reader)


def test_a_reader_for_a_different_annotation_is_refused(annotation, tmp_path):
    other = tmp_path / "other.gtf"
    other.write_text(GTF)
    with pytest.raises(ValueError, match="not"):
        generate_intervals_from_gtf(str(other), reader=_parsed(annotation))


def test_an_empty_but_parsed_annotation_is_not_an_unparsed_one(tmp_path):
    empty = tmp_path / "empty.gtf"
    empty.write_text("# no features\n")
    reader = _parsed(empty)
    assert generate_intervals_from_gtf(str(empty), reader=reader) == []
    assert generate_intervals_from_gtf(str(empty)) == []


def test_generate_isolated_intervals_refuses_a_reader_without_a_path(annotation):
    from fin.io.interval_manager import generate_isolated_intervals

    with pytest.raises(ValueError, match="without a gtf_path"):
        generate_isolated_intervals("/nonexistent.bam", gtf_path=None,
                                    gtf_reader=_parsed(annotation))
