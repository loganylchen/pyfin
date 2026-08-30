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

    def __init__(self, path: str, cache_chroms: int = 2):
        self._path = str(path)
        self._cache_chroms = max(1, int(cache_chroms))
        self._handle = None
        self._names: Optional[tuple] = None
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
        while len(self._cache) > self._cache_chroms:
            self._cache.popitem(last=False)
        return seq

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
        return (LazyGenomeFasta, (self._path, self._cache_chroms))

    def close(self) -> None:
        if self._handle is not None:
            try:
                self._handle.close()
            finally:
                self._handle = None


def open_genome(path: str, *, lazy: bool = True, cache_chroms: int = 2):
    """Open a genome FASTA as a mapping.

    ``lazy=True`` (default) uses :class:`LazyGenomeFasta` when the FASTA can
    be indexed (pysam builds a .fai automatically when the directory is
    writable); on any failure it falls back to the eager whole-genome dict so
    behavior degrades to the historical path rather than erroring.
    """
    if lazy:
        try:
            lg = LazyGenomeFasta(path, cache_chroms=cache_chroms)
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
