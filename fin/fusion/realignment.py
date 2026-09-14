"""Independent genome realignment for clipped arms and competing explanations."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Tuple

import logging

from fin.fusion.records import AlignmentSegment, deduplicate_segments, segment_from_cigar

logger = logging.getLogger(__name__)


_POOL_ALIGNER = None


def _map_one(item):
    """Worker entry: the aligner and its index are inherited through fork."""
    read_id, sequence, segments = item
    return _POOL_ALIGNER.map_molecule(read_id, sequence, segments)


@dataclass(frozen=True)
class GenomeRealignment:
    read_id: str
    status: str
    reason: str
    segments: Tuple[AlignmentSegment, ...]
    searched_windows: Tuple[Tuple[int, int], ...]


def segment_from_hit(read_id, hit, window_start, window_length, read_length):
    """Convert PAF-style query intervals without reversing the molecule order."""
    if (window_start < 0 or window_length <= 0 or window_start + window_length > read_length
            or not 0 <= hit.q_st < hit.q_en <= window_length or hit.strand not in (-1, 1)):
        raise ValueError("invalid realignment query coordinates")
    leading, trailing = hit.q_st, window_length - hit.q_en
    if hit.strand == -1:
        leading, trailing = trailing, leading
    cigar = (f"{leading}S" if leading else "") + hit.cigar_str + (f"{trailing}S" if trailing else "")
    segment = segment_from_cigar(
        read_id, hit.ctg, hit.r_st, "+" if hit.strand == 1 else "-", cigar, hit.mapq,
        edit_distance=hit.NM, source="fusion_genome_realign", secondary=not hit.is_primary,
    )
    if segment.read_length != window_length or segment.blocks[-1].reference_end > hit.ctg_len:
        raise ValueError("realignment CIGAR disagrees with query or reference length")
    if (min(block.reference_start for block in segment.blocks) != hit.r_st
            or max(block.reference_end for block in segment.blocks) != hit.r_en):
        raise ValueError("realignment CIGAR disagrees with reference interval")
    blocks = tuple(replace(block, query_start=block.query_start + window_start,
                           query_end=block.query_end + window_start) for block in segment.blocks)
    return replace(segment, read_length=read_length, blocks=blocks)


class FusionGenomeAligner:
    """Own a splice-aware index and expose mapping limits as unresolved evidence."""

    def __init__(self, reference, *, max_hits=32, max_windows=64, max_segments=128,
                 min_anchor=25, min_identity=0.7, threads=1):
        import mappy

        if (min(max_hits, max_windows, max_segments, min_anchor, threads) <= 0
                or not 0 < min_identity <= 1):
            raise ValueError("invalid fusion genome alignment limits")
        self.max_hits = max_hits
        self.max_windows = max_windows
        self.max_segments = max_segments
        self.min_anchor = min_anchor
        self.min_identity = min_identity
        self.aligner = mappy.Aligner(str(reference), preset="splice", best_n=max_hits,
                                    n_threads=threads)
        if not self.aligner:
            raise ValueError("cannot build fusion genome index")

    def map_many(self, molecules, workers=1):
        """Realign a batch, optionally across forked workers sharing the index.

        mappy releases no GIL for `map`, so threads do not help; `n_threads`
        only speeds the index build. Forking after the index exists gives each
        worker the whole index copy-on-write. Results come back keyed by read
        id and are consumed in the caller's deterministic order, so the number
        of workers cannot change the analysis.

        Operational limits, deliberately not hidden: not thread-safe or
        reentrant (`_POOL_ALIGNER` is process-global, one batch at a time);
        fork-based, so Linux only; an ordinary worker exception propagates and
        the pool is cleaned up, but abrupt native worker death would leave the
        blocking map awaiting a task -- crash supervision is separate work.
        Workers touch only the inherited mappy index, never the genome
        mapping; a worker that needs the genome must reopen it itself.
        """
        if int(workers) < 1:
            raise ValueError("realignment workers must be a positive integer")
        items = list(molecules)
        if workers <= 1 or len(items) < 2:
            return {read_id: self.map_molecule(read_id, sequence, segments)
                    for read_id, sequence, segments in items}
        import multiprocessing

        context = multiprocessing.get_context("fork")
        chunk = max(1, min(64, (len(items) + workers - 1) // workers))
        logger.info("fusion realignment: %d molecules across %d forked workers, chunk %d",
                    len(items), workers, chunk)
        global _POOL_ALIGNER
        _POOL_ALIGNER = self
        try:
            with context.Pool(processes=workers) as pool:
                results = pool.map(_map_one, items, chunksize=chunk)
        finally:
            _POOL_ALIGNER = None
        return {result.read_id: result for result in results}

    def map_molecule(self, read_id, sequence, existing_segments=()):
        """Search whole molecule, existing arms and independently clipped tails.

        Low-MAPQ alternatives remain visible: a partner cannot be declared unique
        by discarding its competing genomic mappings. Mapping is heuristic; a
        completed search is not a proof that no genomic alternative exists.
        """
        if not sequence:
            return GenomeRealignment(read_id, "unresolved", "missing_sequence", (), ())
        sequence = sequence.upper().replace("U", "T")
        length = len(sequence)
        existing_segments = tuple(existing_segments)
        if any(segment.read_id != read_id or segment.read_length != length for segment in existing_segments):
            raise ValueError("existing segments disagree with realigned molecule")
        pending = {(0, length)}
        pending.update((segment.query_start, segment.query_end) for segment in existing_segments)
        for segment in existing_segments:
            for start, end in ((0, segment.query_start), (segment.query_end, length)):
                if end - start >= self.min_anchor:
                    pending.add((start, end))
        searched = set()
        evidence = []
        while pending:
            if len(pending | searched) > self.max_windows:
                return GenomeRealignment(read_id, "unresolved", "window_limit", (), tuple(sorted(searched)))
            start, end = min(pending)
            pending.remove((start, end))
            searched.add((start, end))
            hits = []
            try:
                for hit in self.aligner.map(sequence[start:end]):
                    hits.append(hit)
                    if len(hits) >= self.max_hits:
                        return GenomeRealignment(read_id, "unresolved", "hit_limit", (), tuple(sorted(searched)))
            except (ValueError, RuntimeError) as error:
                return GenomeRealignment(read_id, "unresolved", f"mapping_failure:{type(error).__name__}",
                                         (), tuple(sorted(searched)))
            for hit in hits:
                if hit.blen <= 0 or hit.mlen / hit.blen < self.min_identity:
                    continue
                try:
                    segment = segment_from_hit(read_id, hit, start, end - start, length)
                except ValueError:
                    return GenomeRealignment(read_id, "unresolved", "invalid_mapping", (), tuple(sorted(searched)))
                if segment.aligned_bases < self.min_anchor:
                    continue
                evidence.append(segment)
                if len(deduplicate_segments(evidence)) > self.max_segments:
                    return GenomeRealignment(read_id, "unresolved", "segment_limit", (), tuple(sorted(searched)))
                windows = ((start, segment.query_start), (segment.query_end, end))
                for window in windows:
                    if window[1] - window[0] >= self.min_anchor and window not in searched:
                        pending.add(window)
        return GenomeRealignment(read_id, "complete", "genome_search_completed",
                                 deduplicate_segments(evidence), tuple(sorted(searched)))
