"""Lazy genome mapping: byte equality with the eager loader, cache, pickling."""
import pickle

import pytest

from fin.io.io_fasta import FASTAReader
from fin.io.lazy_genome import LazyGenomeFasta, open_genome


@pytest.fixture()
def fasta(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_text(
        ">chr1 primary assembly\n"
        "ACGTacgtNNNACGT\n"
        "ggttAA\n"
        ">chr2\n"
        "TTTTccccGGGG\n"
        ">chrM circular\n"
        "ACacACac\n"
    )
    return path


def _eager(path):
    seqs = {}
    with FASTAReader(str(path)) as reader:
        for record in reader.iterate_records():
            seqs[record.id] = record.sequence
    return seqs


def test_lazy_equals_eager_byte_for_byte(fasta):
    eager = _eager(fasta)
    lazy = LazyGenomeFasta(str(fasta))
    assert set(lazy) == set(eager)
    for chrom, seq in eager.items():
        assert lazy[chrom] == seq
    assert len(lazy) == len(eager)


def test_mapping_protocol_and_missing_key(fasta):
    lazy = LazyGenomeFasta(str(fasta))
    assert "chr1" in lazy
    assert "chrX" not in lazy
    assert lazy.get("chrX") is None
    assert bool(lazy) is True
    with pytest.raises(KeyError):
        lazy["chrX"]


def test_cache_is_bounded(fasta):
    lazy = LazyGenomeFasta(str(fasta), cache_chroms=1)
    _ = lazy["chr1"]
    _ = lazy["chr2"]
    assert len(lazy._cache) == 1
    assert "chr2" in lazy._cache


def test_pickles_to_path_only(fasta):
    lazy = LazyGenomeFasta(str(fasta), cache_chroms=3)
    _ = lazy["chr1"]  # populate cache + handle
    clone = pickle.loads(pickle.dumps(lazy))
    assert clone._handle is None and not clone._cache
    assert clone["chr2"] == lazy["chr2"]


def test_open_genome_fallback_eager(tmp_path):
    bad = tmp_path / "notfasta.fa"
    bad.write_text("this is not fasta\n")
    # pysam cannot index it -> eager fallback -> plain dict (empty or parsed)
    result = open_genome(str(bad), lazy=True)
    assert isinstance(result, dict)


def test_open_genome_eager_mode(fasta):
    result = open_genome(str(fasta), lazy=False)
    assert isinstance(result, dict)
    assert result == _eager(fasta)


def test_cli_exposes_genome_flags():
    """The documented escape hatch must actually exist on the CLI."""
    from click.testing import CliRunner
    from fin.cli import main

    out = CliRunner().invoke(main, ["--help"]).output
    assert "--lazy-genome" in out and "--no-lazy-genome" in out
    assert "--genome-cache-chroms" in out


def test_runner_cleanup_closes_lazy_genome(fasta):
    from fin.pipeline.runner import PipelineRunner

    runner = PipelineRunner.__new__(PipelineRunner)
    runner._gtf_reader = None
    runner._signal_reader = None
    runner._genome_fasta = LazyGenomeFasta(str(fasta))
    _ = runner._genome_fasta["chr1"]
    assert runner._genome_fasta._handle is not None
    runner.cleanup()
    assert runner._genome_fasta._handle is None
