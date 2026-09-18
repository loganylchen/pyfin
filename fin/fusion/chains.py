"""Deterministic whole-chain support and bounded-diameter event grouping."""

from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Iterable, Tuple

from fin.fusion.annotation import AnnotatedReadPath


def topology_key(path):
    return tuple((arm.chrom, arm.strand, arm.gene_ids, len(arm.blocks)) for arm in path.arms)


def coordinate_key(path):
    coordinates = []
    for index, arm in enumerate(path.arms):
        if index > 0:
            coordinates.append(arm.entry_position)
        for previous, following in zip(arm.blocks, arm.blocks[1:]):
            coordinates.extend((previous.reference_end, following.reference_start)
                               if arm.strand == "+" else
                               (previous.reference_start, following.reference_end))
        if index < len(path.arms) - 1:
            coordinates.append(arm.exit_position)
    return tuple(coordinates)


def _distance(first, second):
    return max((abs(left - right) for left, right in zip(first, second)), default=0)


@dataclass(frozen=True)
class FusionChain:
    chain_id: str
    representative: AnnotatedReadPath
    full_support_read_ids: Tuple[str, ...]
    partial_support_read_ids: Tuple[str, ...]
    edge_support_read_ids: Tuple[Tuple[str, ...], ...]
    raw_coordinate_modes: Tuple[Tuple[Tuple[int, ...], int], ...]
    status: str
    coordinate_status: str


def _is_subpath(partial, complete, tolerance):
    count = len(partial.arms)
    if count >= len(complete.arms) or count < 2:
        return None
    matching = []
    for offset in range(len(complete.arms) - count + 1):
        candidate = AnnotatedReadPath(
            complete.read_id, complete.arms[offset:offset + count],
            complete.event_class, complete.identity_status,
        )
        if (topology_key(partial) == topology_key(candidate)
                and _distance(coordinate_key(partial), coordinate_key(candidate)) <= tolerance):
            matching.append(offset)
    return matching[0] if len(matching) == 1 else None


def build_fusion_chains(
    paths: Iterable[AnnotatedReadPath], *, min_full_support=2, breakpoint_tolerance=6,
) -> Tuple[FusionChain, ...]:
    """Only complete observed paths seed chains; partial edge evidence never does."""
    if min_full_support <= 0 or breakpoint_tolerance < 0:
        raise ValueError("invalid fusion support or breakpoint tolerance")
    by_read = defaultdict(dict)
    for path in paths:
        key = topology_key(path), coordinate_key(path)
        previous = by_read[path.read_id].get(key)
        geometry = lambda candidate: tuple(
            tuple((block.reference_start, block.reference_end,
                   block.query_start, block.query_end, block.aligned_bases)
                  for block in arm.blocks) for arm in candidate.arms
        )
        if previous is None or geometry(path) < geometry(previous):
            by_read[path.read_id][key] = path
    unique_paths = [next(iter(alternatives.values())) for alternatives in by_read.values()
                    if len(alternatives) == 1]
    eligible = [path for path in unique_paths
                if len(path.arms) >= 2 and path.identity_status != "ambiguous"]
    grouped = defaultdict(list)
    for path in eligible:
        grouped[topology_key(path)].append(path)
    clusters = []
    for topology in sorted(grouped):
        ordered = sorted(grouped[topology], key=lambda path: (coordinate_key(path), path.read_id))
        local_clusters = []
        for path in ordered:
            coordinates = coordinate_key(path)
            compatible = [members for members in local_clusters
                          if all(_distance(coordinates, coordinate_key(member)) <= breakpoint_tolerance
                                 for member in members)]
            if compatible:
                chosen = min(compatible, key=lambda members: (
                    max(_distance(coordinates, coordinate_key(member)) for member in members),
                    coordinate_key(members[0]), members[0].read_id,
                ))
                chosen.append(path)
            else:
                local_clusters.append([path])
        clusters.extend(local_clusters)
    output = []
    for members in clusters:
        modes = Counter(coordinate_key(path) for path in members)
        best_count = max(modes.values())
        best_modes = sorted(mode for mode, count in modes.items() if count == best_count)
        representative = min((path for path in members if coordinate_key(path) == best_modes[0]),
                             key=lambda path: (tuple(arm.exons for arm in path.arms), path.read_id))
        full_support = tuple(sorted({path.read_id for path in members}))
        edge_support = [set(full_support) for _ in representative.junctions]
        partial_support = set()
        for partial in eligible:
            if partial.read_id in full_support:
                continue
            offset = _is_subpath(partial, representative, breakpoint_tolerance)
            if offset is not None:
                partial_support.add(partial.read_id)
                for edge_index in range(offset, offset + len(partial.arms) - 1):
                    edge_support[edge_index].add(partial.read_id)
        status = ("insufficient_full_chain_support" if len(full_support) < min_full_support
                  else "readthrough_or_splicing" if representative.event_class == "readthrough_or_splicing_candidate"
                  # a circular reference's seam is continuity, never a fusion,
                  # however many molecules cross it; the support is kept
                  else "circular_reference_origin" if representative.event_class == "circular_reference_origin_spanning"
                  else "supported_candidate")
        coordinate_status = "exact_observed" if len(modes) == 1 else "unresolved_nearby_modes"
        payload = json.dumps((topology_key(representative), sorted(modes)), separators=(",", ":"))
        chain_id = "fusion_chain_" + hashlib.sha256(payload.encode()).hexdigest()[:20]
        output.append(FusionChain(
            chain_id, representative, full_support, tuple(sorted(partial_support)),
            tuple(tuple(sorted(read_ids)) for read_ids in edge_support),
            tuple(sorted(modes.items())), status, coordinate_status,
        ))
    return tuple(sorted(output, key=lambda chain: chain.chain_id))
