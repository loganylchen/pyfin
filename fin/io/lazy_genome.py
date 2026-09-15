"""Lazily-fetched genome FASTA mapping with a bounded per-process cache.

Memory attribution (`experiments/prod_validation/memory_attribution_profile.py`)
showed the dominant pipeline memory cost is per-worker duplicated state led by
the eagerly-loaded whole-genome dict (~3.1 GB x N spawn workers), while the
responsibility ledgers are ~0.1 GB. This mapping keeps the dict contract that
every consumer already uses (``chrom in genome``, ``genome[chrom]``,
``genome.get(chrom)``) but fetches whole-chromosome sequences on demand
through pysam's indexed FASTA reader and retains only a small LRU of
chromosomes per process.

Correctness contract: for an indexed FASTA, ``LazyGenomeFasta(path)[chrom]``
must equal the eager ``FASTAReader`` dict value byte-for-byte (same key = the
header's first word, same raw sequence case). The pipeline's byte-identity
acceptance runs are the enforcement.

Picklable for spawn workers: only the path and cache size cross the process
boundary; each worker opens its own file handle lazily.
"""
from __future__ import annotations

import logging
from collections import OrderedDict
from collections.abc import Mapping
from pathlib import Path
from typing import Iterator, Optional

logger = logging.getLogger(__name__)


class LazyGenomeFasta(Mapping):
    """dict-like read-only genome access backed by an indexed FASTA."""

    def __init__(self, path: str, cache_chroms: int = 2, cache_bases: Optional[int] = None):
        """``cache_chroms`` bounds the cache by COUNT, ``cache_bases`` by SIZE.

        A count bound alone cannot be set safely without knowing the genome: 2
        is right for a pipeline that walks one chromosome at a time and
        catastrophic for one whose access order is effectively random across
        the assembly. Measured on a whole-library fusion run against GRCh38:
        137,042 hits, 499 misses, **497 evictions** -- nearly every miss
        evicting the entry the next one wanted -- and 66.6 Gb fetched from a
        3.1 Gb genome, 21 times over, costing 246 s.

        A size bound is the one a caller can choose without knowing the
        reference: it means the same thing on a spike-in FASTA (everything
        fits) and on a 30 Gb genome (bounded). Pass ``cache_chroms=None`` to
        bound by size alone. An entry larger than the whole budget is still
        kept -- evicting it on insert would turn every access into a miss,
        which is worse than no cache at all.
        """
        self._path = str(path)
        self._cache_chroms = None if cache_chroms is None else max(1, int(cache_chroms))
        self._cache_bases = None if cache_bases is None else max(1, int(cache_bases))
        if self._cache_chroms is None and self._cache_bases is None:
            raise ValueError("genome cache needs a bound on count, size, or both")
        self._handle = None
        self._names: Optional[tuple] = None
        self._cached_bases = 0
        self._cache: "OrderedDict[str, str]" = OrderedDict()

    # -- lazy pysam handle -------------------------------------------------
    def _fasta(self):
        if self._handle is None:
            import pysam

            self._handle = pysam.FastaFile(self._path)
            self._names = tuple(self._handle.references)
        return self._handle

    # -- Mapping protocol --------------------------------------------------
    def __getitem__(self, chrom: str) -> str:
        cached = self._cache.get(chrom)
        if cached is not None:
            self._cache.move_to_end(chrom)
            return cached
        handle = self._fasta()
        if chrom not in self._names_set():
            raise KeyError(chrom)
        seq = handle.fetch(chrom)
        self._cache[chrom] = seq
        self._cache.move_to_end(chrom)
        self._cached_bases += len(seq)
        while len(self._cache) > 1 and self._over_budget():
            _, evicted = self._cache.popitem(last=False)
            self._cached_bases -= len(evicted)
        return seq

    def _over_budget(self):
        if self._cache_chroms is not None and len(self._cache) > self._cache_chroms:
            return True
        return self._cache_bases is not None and self._cached_bases > self._cache_bases

    def _names_set(self):
        if self._names is None:
            self._fasta()
        return set(self._names)

    def __contains__(self, chrom) -> bool:  # Mapping default would fetch
        return chrom in self._names_set()

    def __iter__(self) -> Iterator[str]:
        self._fasta()
        return iter(self._names)

    def __len__(self) -> int:
        self._fasta()
        return len(self._names)

    def __bool__(self) -> bool:  # empty-genome check without loading data
        try:
            return len(self) > 0
        except Exception:
            return False

    # -- pickling for spawn workers ---------------------------------------
    def __reduce__(self):
        return (LazyGenomeFasta, (self._path, self._cache_chroms, self._cache_bases))

    def close(self) -> None:
        if self._handle is not None:
            try:
                self._handle.close()
            finally:
                self._handle = None


def open_genome(path: str, *, lazy: bool = True, cache_chroms: Optional[int] = 2,
                cache_bases: Optional[int] = None):
    """Open a genome FASTA as a mapping.

    ``lazy=True`` (default) uses :class:`LazyGenomeFasta` when the FASTA can
    be indexed (pysam builds a .fai automatically when the directory is
    writable); on any failure it falls back to the eager whole-genome dict so
    behavior degrades to the historical path rather than erroring.
    """
    if lazy:
        try:
            lg = LazyGenomeFasta(path, cache_chroms=cache_chroms, cache_bases=cache_bases)
            lg._fasta()  # force index open so failures fall back here
            return lg
        except Exception as exc:
            logger.warning(
                "lazy genome unavailable for %s (%s); falling back to eager "
                "load", path, exc,
            )
    from fin.io.io_fasta import FASTAReader

    seqs = {}
    with FASTAReader(path) as reader:
        for record in reader.iterate_records():
            seqs[record.id] = record.sequence
    return seqs


__all__ = ["LazyGenomeFasta", "open_genome"]
