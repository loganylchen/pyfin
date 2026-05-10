"""Fusion breakpoint parsing and clustering from BAM SA tags."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

import pysam


@dataclass
class Breakpoint:
    """A split-read fusion breakpoint inferred from a supplementary alignment (SA) tag.

    Attributes:
        chromA: Reference name for the primary alignment (left partner).
        posA: 0-based reference start of the primary alignment.
        strandA: Strand of the primary alignment ('+' or '-').
        chromB: Reference name for the supplementary alignment (right partner).
        posB: 0-based reference start of the supplementary alignment.
        strandB: Strand of the supplementary alignment ('+' or '-').
        support_count: Number of reads supporting this breakpoint.
        supporting_read_ids: Set of query names that support this breakpoint.
    """

    chromA: str
    posA: int
    strandA: str
    chromB: str
    posB: int
    strandB: str
    support_count: int = 1
    supporting_read_ids: Set[str] = field(default_factory=set)


def _parse_region(region: str) -> Tuple[str, Optional[int], Optional[int]]:
    """Parse a region string like 'chr1:1000-2000' into (chrom, start, end)."""
    if ":" in region:
        chrom, coords = region.split(":", 1)
        if "-" in coords:
            start_str, end_str = coords.split("-", 1)
            return chrom, int(start_str), int(end_str)
        return chrom, int(coords), None
    return region, None, None


def parse_sa_tags(
    bam_path: str,
    region: Optional[str] = None,
    min_mapq: int = 10,
) -> List[Breakpoint]:
    """Parse supplementary alignment (SA) tags from a BAM file to extract fusion breakpoints.

    Each SA tag entry is a semicolon-separated list of supplementary alignments in the
    format: ``rname,pos,strand,CIGAR,mapQ,NM;...``

    For each primary read that carries an SA tag, one Breakpoint is emitted per SA
    entry (provided both the primary and supplementary mapping qualities meet
    *min_mapq*).

    Args:
        bam_path: Path to the indexed BAM file.
        region: Optional region string (e.g. ``'chr1:1000-5000'``) passed to
            ``pysam.AlignmentFile.fetch``.  If *None*, all reads are iterated.
        min_mapq: Minimum mapping quality for both the primary and supplementary
            alignment.  Reads below this threshold are skipped.

    Returns:
        A list of :class:`Breakpoint` objects, one per qualifying SA entry.
    """
    breakpoints: List[Breakpoint] = []

    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        if region is not None:
            chrom, start, end = _parse_region(region)
            fetch_iter = bam.fetch(chrom, start, end)
        else:
            fetch_iter = bam.fetch()

        for read in fetch_iter:
            # Skip reads without SA tag or below primary MAPQ threshold
            if not read.has_tag("SA"):
                continue
            if read.mapping_quality is None or read.mapping_quality < min_mapq:
                continue

            strandA = "-" if read.is_reverse else "+"
            chromA = read.reference_name
            # P1-α: posA is the breakpoint coordinate on the primary side.
            # For forward-strand reads the alignment runs left-to-right and the
            # breakpoint sits at reference_end; for reverse-strand reads the
            # alignment runs right-to-left and the breakpoint sits at
            # reference_start.
            if read.is_reverse:
                posA = read.reference_start
            else:
                posA = read.reference_end if read.reference_end is not None else read.reference_start

            sa_string: str = read.get_tag("SA")
            # SA tag ends with ';' so the last element after split is empty
            for entry in sa_string.rstrip(";").split(";"):
                entry = entry.strip()
                if not entry:
                    continue
                parts = entry.split(",")
                if len(parts) < 6:
                    continue
                chromB = parts[0]
                posB = int(parts[1]) - 1  # SA pos is 1-based; convert to 0-based
                strandB = parts[2]
                # parts[3] = CIGAR, parts[4] = mapQ, parts[5] = NM
                try:
                    sa_mapq = int(parts[4])
                except ValueError:
                    continue
                if sa_mapq < min_mapq:
                    continue

                bp = Breakpoint(
                    chromA=chromA,
                    posA=posA,
                    strandA=strandA,
                    chromB=chromB,
                    posB=posB,
                    strandB=strandB,
                    support_count=1,
                    supporting_read_ids={read.query_name},
                )
                breakpoints.append(bp)
    finally:
        bam.close()

    return breakpoints


def cluster_breakpoints(
    breakpoints: List[Breakpoint],
    max_dist: int = 500,
    min_support: int = 2,
) -> List[Breakpoint]:
    """Cluster nearby breakpoints and filter by minimum support.

    Two breakpoints are merged into the same cluster (single-linkage) when ALL
    of the following hold for *some* pair within the cluster:

    * ``chromA``, ``chromB``, ``strandA``, and ``strandB`` are identical.
    * ``|posA_i - posA_j| <= max_dist``
    * ``|posB_i - posB_j| <= max_dist``

    P1-β: clustering is implemented via union-find on edges that satisfy the
    proximity predicate, so the result is *order-independent* — feeding the
    same multiset of breakpoints in any permutation yields the same clusters.
    The merged representative uses the rounded mean of all member positions
    for *posA* and *posB*.

    Args:
        breakpoints: Raw breakpoints, typically from :func:`parse_sa_tags`.
        max_dist: Maximum coordinate distance on each side for two breakpoints
            to be considered the same event.
        min_support: Clusters whose ``support_count`` falls below this value
            after merging are discarded.

    Returns:
        Merged and filtered :class:`Breakpoint` list (deterministic order).
    """
    if not breakpoints:
        return []

    n = len(breakpoints)
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx == ry:
            return
        # Lower index becomes parent → deterministic regardless of call order.
        if rx < ry:
            parent[ry] = rx
        else:
            parent[rx] = ry

    # Group indices by (chromA, chromB, strandA, strandB) signature.
    groups: Dict[Tuple[str, str, str, str], List[int]] = defaultdict(list)
    for i, bp in enumerate(breakpoints):
        key = (bp.chromA, bp.chromB, bp.strandA, bp.strandB)
        groups[key].append(i)

    # Within each group, sort by posA so we can short-circuit the inner loop
    # once posA difference exceeds max_dist. Output is identical to a full
    # O(g²) scan because we only skip pairs that cannot be neighbours.
    for indices in groups.values():
        indices.sort(key=lambda i: (breakpoints[i].posA, breakpoints[i].posB, i))
        for a in range(len(indices)):
            ia = indices[a]
            bp_a = breakpoints[ia]
            for b in range(a + 1, len(indices)):
                ib = indices[b]
                bp_b = breakpoints[ib]
                if bp_b.posA - bp_a.posA > max_dist:
                    break  # sorted: no further j can be closer on posA
                if abs(bp_b.posB - bp_a.posB) <= max_dist:
                    union(ia, ib)

    # Collect connected components.
    clusters_by_root: Dict[int, List[Breakpoint]] = defaultdict(list)
    for i, bp in enumerate(breakpoints):
        clusters_by_root[find(i)].append(bp)

    # Deterministic output ordering: sort clusters by signature + min positions.
    cluster_lists = list(clusters_by_root.values())
    cluster_lists.sort(
        key=lambda cl: (
            cl[0].chromA,
            cl[0].chromB,
            cl[0].strandA,
            cl[0].strandB,
            min(m.posA for m in cl),
            min(m.posB for m in cl),
        )
    )

    merged: List[Breakpoint] = []
    for cluster in cluster_lists:
        total_support = sum(m.support_count for m in cluster)
        if total_support < min_support:
            continue
        mean_posA = round(sum(m.posA for m in cluster) / len(cluster))
        mean_posB = round(sum(m.posB for m in cluster) / len(cluster))
        all_read_ids: Set[str] = set()
        for m in cluster:
            all_read_ids |= m.supporting_read_ids
        rep = cluster[0]
        merged.append(
            Breakpoint(
                chromA=rep.chromA,
                posA=mean_posA,
                strandA=rep.strandA,
                chromB=rep.chromB,
                posB=mean_posB,
                strandB=rep.strandB,
                support_count=total_support,
                supporting_read_ids=all_read_ids,
            )
        )

    return merged
