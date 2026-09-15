"""The genome cache's bounds.

`LazyGenomeFasta` returns a WHOLE chromosome per `__getitem__`, so its LRU is
the difference between one fetch per chromosome and one fetch per access. The
count bound of 2 was chosen for the assembly pipeline, which sweeps one locus
at a time. The fusion adjudication walks candidate molecules in read-id order,
so the chromosomes it needs arrive in effectively random order, and measured on
a whole-genome run under that default it took 137,042 hits, 499 misses, **497
evictions**, and fetched 66.6 Gb from a 3.1 Gb genome -- 21 times over, 246 s.

A size bound is the one a caller can pick without knowing the reference.
"""
import pytest

from fin.io.lazy_genome import LazyGenomeFasta, open_genome

CHROMS = {"c1": "A" * 400, "c2": "C" * 400, "c3": "G" * 400, "big": "T" * 5000}


@pytest.fixture
def reference(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in CHROMS.items()))
    return path


class _CountingHandle:
    """pysam's FastaFile is a C extension whose methods cannot be rebound, so
    count fetches by standing in front of the handle rather than patching it."""

    def __init__(self, handle, fetched):
        self._handle = handle
        self._fetched = fetched
        self.references = handle.references

    def fetch(self, chrom, *args, **kwargs):
        self._fetched.append(chrom)
        return self._handle.fetch(chrom, *args, **kwargs)

    def close(self):
        self._handle.close()


def _counting(genome):
    fetched = []
    genome._handle = _CountingHandle(genome._fasta(), fetched)
    return fetched


def test_the_count_bound_thrashes_on_a_rotating_access_pattern(reference):
    """The behaviour that motivated the change, stated as a test rather than a
    claim: three chromosomes visited in rotation, a cache of two, and every
    single access is a miss."""
    genome = LazyGenomeFasta(str(reference), cache_chroms=2)
    fetched = _counting(genome)
    for _ in range(4):
        for name in ("c1", "c2", "c3"):
            assert genome[name] == CHROMS[name]
    assert len(fetched) == 12  # every access


def test_a_size_bound_holds_the_same_working_set(reference):
    genome = LazyGenomeFasta(str(reference), cache_chroms=None, cache_bases=2000)
    fetched = _counting(genome)
    for _ in range(4):
        for name in ("c1", "c2", "c3"):
            assert genome[name] == CHROMS[name]
    assert len(fetched) == 3  # one per chromosome


def test_a_budget_smaller_than_the_working_set_thrashes_exactly_as_before(reference):
    """The honest limit of this change: bounding by SIZE does not conjure
    capacity. A budget of 900 bases holds two 400-base chromosomes, and three
    visited in rotation still miss on every access -- the same pathology as the
    count bound, for the same reason. What the size bound buys is that a caller
    can choose a budget correctly without knowing the reference; it is not a
    better eviction policy."""
    genome = LazyGenomeFasta(str(reference), cache_chroms=None, cache_bases=900)
    fetched = _counting(genome)
    for _ in range(3):
        for name in ("c1", "c2", "c3"):
            assert genome[name] == CHROMS[name]
    assert len(fetched) == 9
    assert len(genome._cache) == 2  # it did evict; it just could not help


def test_an_entry_larger_than_the_whole_budget_is_still_kept(reference):
    """Evicting it on insert would make every access a miss -- worse than no
    cache. One entry always survives."""
    genome = LazyGenomeFasta(str(reference), cache_chroms=None, cache_bases=100)
    fetched = _counting(genome)
    for _ in range(3):
        assert genome["big"] == CHROMS["big"]
    assert len(fetched) == 1
    assert list(genome._cache) == ["big"]


def test_both_bounds_apply_when_both_are_given(reference):
    genome = LazyGenomeFasta(str(reference), cache_chroms=2, cache_bases=10**9)
    for name in ("c1", "c2", "c3"):
        genome[name]
    assert list(genome._cache) == ["c2", "c3"]


def test_a_cache_with_no_bound_at_all_is_refused(reference):
    with pytest.raises(ValueError, match="bound"):
        LazyGenomeFasta(str(reference), cache_chroms=None, cache_bases=None)


def test_eviction_accounting_does_not_drift(reference):
    genome = LazyGenomeFasta(str(reference), cache_chroms=None, cache_bases=900)
    for _ in range(5):
        for name in ("c1", "c2", "c3", "c1"):
            genome[name]
    assert genome._cached_bases == sum(len(genome._cache[name]) for name in genome._cache)


def test_the_pickle_form_carries_both_bounds(reference):
    genome = LazyGenomeFasta(str(reference), cache_chroms=None, cache_bases=4096)
    factory, args = genome.__reduce__()
    rebuilt = factory(*args)
    assert (rebuilt._cache_chroms, rebuilt._cache_bases) == (None, 4096)
    assert rebuilt["c1"] == CHROMS["c1"]


def test_open_genome_passes_the_size_bound_through(reference):
    genome = open_genome(str(reference), cache_chroms=None, cache_bases=4096)
    assert genome._cache_bases == 4096 and genome._cache_chroms is None
    assert genome["c2"] == CHROMS["c2"]


def test_the_default_is_unchanged_for_every_other_caller(reference):
    genome = open_genome(str(reference))
    assert (genome._cache_chroms, genome._cache_bases) == (2, None)
