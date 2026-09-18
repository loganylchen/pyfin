"""Bounded full-molecule path search; no cross-read fusion-path assembly."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple

from fin.fusion.records import AlignmentSegment, deduplicate_segments


@dataclass(frozen=True)
class ReadPath:
    read_id: str
    segments: Tuple[AlignmentSegment, ...]

    def __post_init__(self):
        if not self.segments or any(segment.read_id != self.read_id for segment in self.segments):
            raise ValueError("read path must contain segments from one molecule")
        if len({segment.read_length for segment in self.segments}) != 1:
            raise ValueError("read path segments disagree on full query length")
        for previous, following in zip(self.segments, self.segments[1:]):
            if (following.query_start <= previous.query_start
                    or following.query_end <= previous.query_end):
                raise ValueError("read path must advance along the full forward query")

    @property
    def query_covered_bases(self) -> int:
        intervals = sorted((block.query_start, block.query_end)
                           for segment in self.segments for block in segment.blocks)
        total = 0
        covered_end = 0
        for start, end in intervals:
            total += max(0, end - max(covered_end, start))
            covered_end = max(covered_end, end)
        return total

    @property
    def coverage_fraction(self) -> float:
        return self.query_covered_bases / self.segments[0].read_length

    def effective_coverage(self, domain) -> float:
        """Coverage over the molecule-level domain ``(lead, trail)``: the
        terminal bases the genome search could not place anywhere are excluded
        from the denominator (see :func:`coverage_domain`)."""
        lead, trail = domain
        length = self.segments[0].read_length
        lo, hi = lead, length - trail
        if hi <= lo:
            return 0.0
        intervals = sorted((max(lo, block.query_start), min(hi, block.query_end))
                           for segment in self.segments for block in segment.blocks)
        total = 0
        covered_end = lo
        for start, end in intervals:
            if end <= start:
                continue
            total += max(0, end - max(covered_end, start))
            covered_end = max(covered_end, end)
        return total / (hi - lo)

    @property
    def junctions(self) -> tuple:
        return tuple(
            ((previous.chrom, previous.exit_position, previous.strand),
             (following.chrom, following.entry_position, following.strand))
            for previous, following in zip(self.segments, self.segments[1:])
        )

    @property
    def geometry_key(self) -> tuple:
        return tuple(segment.geometry_key for segment in self.segments)


@dataclass(frozen=True)
class PathSearchResult:
    read_id: str
    status: str
    reason: str
    paths: Tuple[ReadPath, ...]
    segment_count: int


def coverage_domain(segments: Iterable[AlignmentSegment], read_length: int, *,
                    max_end_bases: int = 100, max_total_bases: int = 120, max_fraction: float = 0.2):
    """Terminal query bases no detailed alignment could place: ``(lead, trail)``.

    Direct-RNA molecules carry non-genomic sequence at the ends (adapter and
    poly(A) carry-through at the 3' end, a few bases at the 5' end). Judged on
    the UNION of every detailed mapping of the molecule -- low-MAPQ, normal and
    paralogous alternatives included -- so an independently placed terminal arm
    is never excluded. Only a whole eligible overhang is excluded, never a
    slice of a longer one, and never anything internal. The bounds are
    proposed acceptance settings, not literature constants.
    """
    if read_length <= 0 or min(max_end_bases, max_total_bases) < 0 or not 0 <= max_fraction <= 1:
        raise ValueError("invalid coverage-domain bounds")
    placed = [(block.query_start, block.query_end) for segment in segments
              if _detailed(segment) for block in segment.blocks]
    if not placed:
        return 0, 0
    lead = min(start for start, _ in placed)
    trail = read_length - max(end for _, end in placed)
    lead = lead if 0 < lead <= max_end_bases else 0
    trail = trail if 0 < trail <= max_end_bases else 0
    if lead + trail > max_total_bases or lead + trail > max_fraction * read_length:
        # too much unexplained sequence to wave through: keep the full read
        # as the denominator and let the coverage rule speak
        return 0, 0
    return lead, trail


def _detailed(segment):
    return any(source != "sa" for source in segment.sources)


def find_read_paths(
    segments: Iterable[AlignmentSegment],
    *,
    min_mapping_quality: int = 20,
    min_anchor_bases: int = 25,
    max_query_gap: int = 29,
    max_query_overlap: int = 10,
    min_read_coverage: float = 0.8,
    coverage_tolerance_bases: int = 5,
    max_segments: int = 64,
    max_paths: int = 64,
    domain=(0, 0),
    max_expansions: int = 20000,
) -> PathSearchResult:
    """Keep every coverage-qualified path or abstain when search cannot be completed.

    A unique path is only a mapping result, not a confident fusion call. Normal
    whole-read mappings remain competitors and gene/path validation comes later.
    Coverage is judged over ``domain`` (see :func:`coverage_domain`), and NO
    path is discarded for merely covering a few bases fewer than the best one:
    covered bases are not a defensible substitute for the sequence adjudication
    that follows, and a slightly longer wrong geometry must not erase the right
    hypothesis. ``coverage_tolerance_bases`` is retained for signature
    compatibility and no longer prunes.

    EVERY directed path through the eligible segments is a hypothesis, not only
    the maximal ones. The search used to start at graph roots and emit only
    leaves, so an alignment that could be extended never appeared on its own:
    chr1 950M50S alone returned the chr1 path (coverage 0.95), but with a chr2
    950S50M tail available only chr1->chr2 came back, and the still-qualifying
    normal explanation was gone before the adjudication that exists to weigh
    it. A walk now starts at every node and every prefix is a candidate; each
    index-path has one first node and one extension sequence, so nothing is
    produced twice. Whether the shorter explanation is COMPETITIVE is decided
    downstream by the stitched-template score, which already charges for every
    unaligned query base -- enumeration stays neutral.

    Two budgets, kept distinct because a small frontier is not a small search:
    ``max_paths`` bounds the QUALIFYING paths retained (the number adjudication
    must weigh), and ``max_expansions`` bounds the walk itself, which a graph of
    many mutually-exclusive placements can make exponential while every path
    fails coverage. Either exhausted means abstain with nothing, never a
    truncated set presented as complete.
    """
    if not 0 <= min_mapping_quality <= 254 or min_anchor_bases <= 0:
        raise ValueError("invalid mapping quality or anchor threshold")
    if (min(max_query_gap, max_query_overlap, coverage_tolerance_bases) < 0
            or min(max_segments, max_paths, max_expansions) <= 0
            or not 0 < min_read_coverage <= 1):
        raise ValueError("invalid path search bounds")
    if min(domain) < 0:
        raise ValueError("invalid coverage domain")
    supplied = tuple(segments)
    read_ids = {segment.read_id for segment in supplied}
    if len(read_ids) > 1:
        raise ValueError("path search cannot combine different molecules")
    read_id = next(iter(read_ids), "")
    if len({segment.read_length for segment in supplied}) > 1:
        return PathSearchResult(read_id, "unresolved", "inconsistent_read_length", (), len(supplied))
    unique = deduplicate_segments(supplied)
    if len(unique) > max_segments:
        return PathSearchResult(read_id, "unresolved", "segment_limit", (), len(unique))
    eligible = sorted(
        (segment for segment in unique
         if segment.effective_mapping_quality >= min_mapping_quality
         and segment.aligned_bases >= min_anchor_bases),
        key=lambda segment: (segment.query_start, segment.query_end, segment.geometry_key),
    )
    if not eligible:
        return PathSearchResult(read_id, "unresolved", "no_eligible_alignment", (), len(unique))
    successors = [[] for _ in eligible]
    for previous_index, previous in enumerate(eligible):
        for following_index in range(previous_index + 1, len(eligible)):
            following = eligible[following_index]
            gap = following.query_start - previous.query_end
            if (following.query_start > previous.query_start
                    and following.query_end > previous.query_end
                    and -max_query_overlap <= gap <= max_query_gap):
                successors[previous_index].append(following_index)
    covered = []
    expansions = 0
    # Starts are taken one at a time, not queued together: a chain of N
    # segments has N starts, and counting them as pending work would abstain
    # on a molecule whose only qualifying path is the full chain.
    for start in range(len(eligible)):
        stack = [(start,)]
        while stack:
            indices = stack.pop()
            expansions += 1
            if expansions > max_expansions:
                return PathSearchResult(read_id, "unresolved", "search_budget", (), len(unique))
            path = ReadPath(read_id, tuple(eligible[index] for index in indices))
            if path.effective_coverage(domain) >= min_read_coverage:
                covered.append(path)
                if len(covered) > max_paths:
                    return PathSearchResult(read_id, "unresolved", "path_limit", (), len(unique))
            stack.extend(indices + (index,) for index in reversed(successors[indices[-1]]))
    if not covered:
        return PathSearchResult(read_id, "unresolved", "insufficient_query_coverage", (), len(unique))
    competitive = tuple(sorted(covered, key=lambda path: path.geometry_key))
    status = "unique" if len(competitive) == 1 else "ambiguous"
    reason = "one_query_path" if status == "unique" else "competing_query_paths"
    return PathSearchResult(read_id, status, reason, competitive, len(unique))
