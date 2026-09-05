"""The junction registry is collected on every run and recorded as a sidecar."""
from fin.analysis.junction_support import merge_support_shards, write_support_shard


class _Iv:
    def __init__(self, chrom, start, end, strand):
        self.chrom, self.start, self.end, self.strand = chrom, start, end, strand
        self.region_string = f"{chrom}:{start}-{end}"


def test_shards_merge_sorted_and_shard_dir_removed(tmp_path):
    write_support_shard(tmp_path, _Iv("chr2", 0, 10, "+"), {(5, 8): (3, 1)})
    write_support_shard(tmp_path, _Iv("chr1", 0, 10, "-"), {(7, 9): (2, 2), (1, 4): (9, 5)})
    out = merge_support_shards(tmp_path)
    lines = open(out).read().strip().split("\n")
    assert lines[0].split("\t")[0] == "chrom"
    assert [l.split("\t")[:5] for l in lines[1:]] == [
        ["chr1", "-", "chr1:0-10", "1", "4"],
        ["chr1", "-", "chr1:0-10", "7", "9"],
        ["chr2", "+", "chr2:0-10", "5", "8"],
    ]
    assert not (tmp_path / "_jsupport_shards").exists()


def test_no_shards_still_writes_the_header(tmp_path):
    # An empty registry must still leave the sidecar: a missing file cannot be
    # told apart from "the recording step never ran".
    from fin.analysis.junction_support import TSV_FIELDS

    out = merge_support_shards(tmp_path)
    assert out is not None
    lines = (tmp_path / "junction_support.tsv").read_text().splitlines()
    assert lines == ["\t".join(TSV_FIELDS)]


def test_discovery_always_collects_evidence():
    """collect_chain_evidence no longer gates the scan: depths and junction_support are
    populated even when no opt-in gate is on."""
    import inspect

    from fin.candidates import discovery

    src = inspect.getsource(discovery._chain_cluster_candidates)
    assert "if collect_chain_evidence and" not in src
    assert "if collect_chain_evidence else" not in src
