"""Fusion breakpoint parsing and clustering from BAM SA tags."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Set, Tuple

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
            posA = read.reference_start

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

    Two breakpoints are merged into the same cluster when ALL of the following
    hold:

    * ``chromA``, ``chromB``, ``strandA``, and ``strandB`` are identical.
    * ``|posA_i - posA_j| <= max_dist``
    * ``|posB_i - posB_j| <= max_dist``

    Clustering is performed with a simple greedy O(N²) pass (adequate for the
    low breakpoint counts expected in practice).  The merged representative uses
    the rounded mean of all member positions for *posA* and *posB*.

    Args:
        breakpoints: Raw breakpoints, typically from :func:`parse_sa_tags`.
        max_dist: Maximum coordinate distance on each side for two breakpoints
            to be considered the same event.
        min_support: Clusters whose ``support_count`` falls below this value
            after merging are discarded.

    Returns:
        Merged and filtered :class:`Breakpoint` list.
    """
    if not breakpoints:
        return []

    # Each element in `clusters` is a list of Breakpoints assigned to that cluster
    clusters: List[List[Breakpoint]] = []

    for bp in breakpoints:
        placed = False
        for cluster in clusters:
            rep = cluster[0]
            # Must share chrom/strand signature
            if (
                rep.chromA != bp.chromA
                or rep.chromB != bp.chromB
                or rep.strandA != bp.strandA
                or rep.strandB != bp.strandB
            ):
                continue
            # Check proximity against every member already in the cluster
            compatible = all(
                abs(m.posA - bp.posA) <= max_dist and abs(m.posB - bp.posB) <= max_dist
                for m in cluster
            )
            if compatible:
                cluster.append(bp)
                placed = True
                break
        if not placed:
            clusters.append([bp])

    merged: List[Breakpoint] = []
    for cluster in clusters:
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
