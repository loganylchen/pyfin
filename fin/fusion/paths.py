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
) -> PathSearchResult:
    """Keep all competitive paths or abstain when search cannot be completed.

    A unique path is only a mapping result, not a confident fusion call. Normal
    whole-read mappings remain competitors and gene/path validation comes later.
    """
    if not 0 <= min_mapping_quality <= 254 or min_anchor_bases <= 0:
        raise ValueError("invalid mapping quality or anchor threshold")
    if (min(max_query_gap, max_query_overlap, coverage_tolerance_bases) < 0
            or min(max_segments, max_paths) <= 0 or not 0 < min_read_coverage <= 1):
        raise ValueError("invalid path search bounds")
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
    has_predecessor = set()
    for previous_index, previous in enumerate(eligible):
        for following_index in range(previous_index + 1, len(eligible)):
            following = eligible[following_index]
            gap = following.query_start - previous.query_end
            if (following.query_start > previous.query_start
                    and following.query_end > previous.query_end
                    and -max_query_overlap <= gap <= max_query_gap):
                successors[previous_index].append(following_index)
                has_predecessor.add(following_index)
    pending = [(index,) for index in reversed(range(len(eligible))) if index not in has_predecessor]
    completed = []
    while pending:
        indices = pending.pop()
        following = successors[indices[-1]]
        if following:
            pending.extend(indices + (index,) for index in reversed(following))
        else:
            completed.append(ReadPath(read_id, tuple(eligible[index] for index in indices)))
        if len(pending) + len(completed) > max_paths:
            return PathSearchResult(read_id, "unresolved", "path_limit", (), len(unique))
    covered = [path for path in completed if path.coverage_fraction >= min_read_coverage]
    if not covered:
        return PathSearchResult(read_id, "unresolved", "insufficient_query_coverage", (), len(unique))
    best_coverage = max(path.query_covered_bases for path in covered)
    competitive = tuple(sorted(
        (path for path in covered
         if path.query_covered_bases >= best_coverage - coverage_tolerance_bases),
        key=lambda path: path.geometry_key,
    ))
    status = "unique" if len(competitive) == 1 else "ambiguous"
    reason = "one_query_path" if status == "unique" else "competing_query_paths"
    return PathSearchResult(read_id, status, reason, competitive, len(unique))
