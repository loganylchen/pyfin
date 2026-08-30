"""Medium 9: hard-negative junction evidence requires a complete region scan.

`get_reads_in_region` historically swallowed a mid-iteration fetch error and
returned the partial accumulation, so support gates read undercounts as
biological zero. The checked API makes completeness explicit, and
`compute_observed_junctions` discards partial evidence so gates abstain.
"""
from types import SimpleNamespace

import pytest

import fin.io.io_bam as io_bam
import fin.pipeline.evidence as evidence
from fin.io.io_bam import BamReader, RegionReadBatch
from fin.io.interval_manager import GenomicInterval


def _reader(fetch_factory):
    reader = BamReader.__new__(BamReader)
    reader._alignment_file = SimpleNamespace(fetch=fetch_factory)
    reader.alignment_to_dict = lambda alignment: alignment
    return reader


def _interval():
    return GenomicInterval(chrom="chr1", start=0, end=1000, strand="+")


READ = {
    "is_secondary": False,
    "is_supplementary": False,
    "is_mapped": True,
    "is_reverse": False,
    "reference_start": 10,
    # 20M 100N 20M -> one intron at (30, 130)
    "cigartuples": ((0, 20), (3, 100), (0, 20)),
}


def test_checked_fetch_marks_mid_iteration_error_incomplete():
    def fetch(**_kw):
        yield READ
        raise OSError("truncated bgzf block")

    batch = _reader(fetch).get_reads_in_region_checked("chr1:1-1000")
    assert batch.reads == [READ]
    assert batch.complete is False
    assert "OSError" in batch.error


def test_legacy_fetch_keeps_partial_list_contract():
    def fetch(**_kw):
        yield READ
        raise OSError("boom")

    assert _reader(fetch).get_reads_in_region("chr1:1-1000") == [READ]


def test_checked_fetch_max_reads_truncation_is_incomplete():
    def fetch(**_kw):
        yield READ
        yield READ

    batch = _reader(fetch).get_reads_in_region_checked("chr1:1-1000", max_reads=1)
    assert len(batch.reads) == 1
    assert batch.complete is False
    assert batch.error is None


def test_checked_fetch_exhausted_iterator_is_complete():
    def fetch(**_kw):
        yield READ

    batch = _reader(fetch).get_reads_in_region_checked("chr1:1-1000")
    assert batch.complete is True
    assert batch.error is None
    assert len(batch.reads) == 1


def test_unopened_reader_raises():
    reader = BamReader.__new__(BamReader)
    reader._alignment_file = None
    with pytest.raises(RuntimeError):
        reader.get_reads_in_region_checked("chr1:1-1000")


class _FakeBam:
    def __init__(self, batch):
        self._batch = batch

    def __call__(self, _path):
        return self

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def get_reads_in_region_checked(self, region, max_reads=None):
        return self._batch


def _patched_observed(monkeypatch, tmp_path, batch):
    bam = tmp_path / "input.bam"
    bam.write_bytes(b"x")
    monkeypatch.setattr(io_bam, "BamReader", _FakeBam(batch))
    return evidence.compute_observed_junctions(str(bam), _interval())


def test_observed_junctions_abstain_on_incomplete_fetch(monkeypatch, tmp_path):
    partial = RegionReadBatch(reads=[READ], complete=False, error="OSError: x")
    assert _patched_observed(monkeypatch, tmp_path, partial) is None


def test_observed_junctions_built_from_complete_fetch(monkeypatch, tmp_path):
    full = RegionReadBatch(reads=[READ], complete=True)
    observed = _patched_observed(monkeypatch, tmp_path, full)
    assert observed is not None
    assert observed["+"][(30, 130)] == 1


def test_observed_junctions_missing_bam_is_none():
    assert evidence.compute_observed_junctions("/nonexistent.bam", _interval()) is None


def test_interval_bundle_memoizes_single_computation(monkeypatch, tmp_path):
    calls = []
    full = RegionReadBatch(reads=[READ], complete=True)

    class CountingBam(_FakeBam):
        def get_reads_in_region_checked(self, region, max_reads=None):
            calls.append(region)
            return self._batch

    bam = tmp_path / "input.bam"
    bam.write_bytes(b"x")
    monkeypatch.setattr(io_bam, "BamReader", CountingBam(full))
    bundle = evidence.IntervalBundle(interval=_interval(), bam_path=str(bam))
    first = bundle.observed_junctions
    second = bundle.observed_junctions
    assert first is second
    assert len(calls) == 1
