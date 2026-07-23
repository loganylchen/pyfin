"""Mono (single-exon) GTF matching in chain-cluster discovery.

Single-exon transcripts all share the empty intron chain, so matching them to GTF by
chain equality collapses every single-exon read onto the first mono GTF in an interval.
_best_overlap_gtf matches by genomic overlap instead.
"""
from fin.candidates.discovery import _best_overlap_gtf


class _GC:
    def __init__(self, cid):
        self.candidate_id = cid


# two single-exon GTF transcripts at different loci in one interval
MONO = [(_GC("X"), 1000, 1200), (_GC("Y"), 5000, 5200)]


def test_overlap_picks_correct_locus_not_first():
    # reads over Y must match Y, not the first-inserted X (the collapse bug)
    assert _best_overlap_gtf(MONO, 5010, 5190).candidate_id == "Y"
    assert _best_overlap_gtf(MONO, 1010, 1190).candidate_id == "X"


def test_no_overlap_returns_none_for_novel():
    assert _best_overlap_gtf(MONO, 9000, 9100) is None


def test_empty_gtf_returns_none_p00_safe():
    # no GTF (p00): nothing to match -> None -> novel candidate (unchanged behaviour)
    assert _best_overlap_gtf([], 1010, 1190) is None


def test_partial_overlap_prefers_larger_overlap():
    # a read straddling both leans to the transcript it overlaps more
    assert _best_overlap_gtf(MONO, 1150, 5050).candidate_id == "X"   # 50bp X vs 50bp Y tie -> X first on tie
    assert _best_overlap_gtf(MONO, 1190, 5150).candidate_id == "Y"   # 10bp X vs 150bp Y -> Y


from fin.candidates.discovery import _spatial_read_clusters


def test_spatial_clusters_split_distinct_loci():
    # two single-exon loci: reads r1,r2 near 1000; r3,r4 near 5000 -> two groups
    spans = {"r1": (1000, 1100), "r2": (1050, 1150), "r3": (5000, 5100), "r4": (5050, 5150)}
    groups = _spatial_read_clusters(set(spans), spans)
    assert len(groups) == 2
    assert {frozenset(g) for g in groups} == {frozenset(["r1", "r2"]), frozenset(["r3", "r4"])}


def test_spatial_clusters_merge_overlapping():
    spans = {"a": (1000, 1200), "b": (1150, 1300), "c": (1250, 1400)}  # chained overlap
    groups = _spatial_read_clusters(set(spans), spans)
    assert len(groups) == 1


def test_spatial_clusters_drop_spanless_reads():
    spans = {"a": (1000, 1100)}
    groups = _spatial_read_clusters({"a", "nospan"}, spans)
    assert groups == [["a"]]
