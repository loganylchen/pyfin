"""Intron chain extraction and read grouping for candidate discovery."""

from __future__ import annotations

import logging
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

from fin.candidates.dataclasses import IntronChain
from fin.io.interval_manager import extract_three_prime_pos

logger = logging.getLogger(__name__)

# pysam CIGAR operation codes
_CIGAR_N = 3  # skipped region (intron)


def extract_intron_chain(
    cigartuples: List[Tuple[int, int]], reference_start: int
) -> IntronChain:
    """Extract intron chain from CIGAR tuples.

    Walks the CIGAR string collecting N (skip) operations as introns.

    Args:
        cigartuples: List of (operation, length) from pysam alignment.
        reference_start: 0-based reference start position.

    Returns:
        IntronChain with sorted intron intervals.
    """
    introns = []
    ref_pos = reference_start

    for op, length in cigartuples:
        if op == _CIGAR_N:
            introns.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op in (0, 2, 7, 8):  # M, D, =, X consume reference
            ref_pos += length
        # S(4), I(1), H(5), P(6) do not consume reference

    return IntronChain(introns=tuple(introns))


def gtf_transcript_to_intron_chain(transcript) -> IntronChain:
    """Convert a GTFTranscript to an IntronChain via gaps between sorted exons.

    Args:
        transcript: GTFTranscript instance with .exons attribute.

    Returns:
        IntronChain derived from gaps between exons.
    """
    exons = sorted(transcript.exons, key=lambda e: e[0])
    introns = []
    for i in range(len(exons) - 1):
        intron_start = exons[i][1]
        intron_end = exons[i + 1][0]
        if intron_end > intron_start:
            introns.append((intron_start, intron_end))
    return IntronChain(introns=tuple(introns))


def _cluster_positions(positions: List[int], gap: int) -> List[List[int]]:
    """Greedy 1D clustering: positions within ``gap`` bp join the active cluster.

    Args:
        positions: Unsorted list of integer coordinates. Duplicates are kept.
        gap: Maximum allowed distance from the previous (sorted) position to
            stay in the same cluster.

    Returns:
        List of clusters; each cluster is a list of positions in sorted order
        (duplicates preserved so callers can compute support counts).
    """
    if not positions:
        return []
    sorted_pos = sorted(positions)
    clusters: List[List[int]] = [[sorted_pos[0]]]
    for p in sorted_pos[1:]:
        if p - clusters[-1][-1] <= gap:
            clusters[-1].append(p)
        else:
            clusters.append([p])
    return clusters


def _has_canonical_motif(
    pos: int,
    kind: str,
    strand: str,
    genome_seq: str,
) -> bool:
    """Check whether a donor/acceptor position carries the canonical splice motif.

    Splice site conventions (0-based, half-open intron [s, e)):
        + strand: donor=GT at genome[s:s+2], acceptor=AG at genome[e-2:e].
        - strand: donor=AC at genome[e-2:e], acceptor=CT at genome[s:s+2].

    Args:
        pos: For a donor on + strand or acceptor on - strand, this is the
            intron start ``s``. For an acceptor on + strand or donor on -
            strand, this is the intron end ``e``.
        kind: "donor" or "acceptor".
        strand: "+" or "-".
        genome_seq: Full chromosome sequence string.

    Returns:
        True if the 2 bp motif at the expected offset matches.
    """
    if not genome_seq:
        return False
    n = len(genome_seq)
    if strand == "+":
        if kind == "donor":
            if pos < 0 or pos + 2 > n:
                return False
            return genome_seq[pos:pos + 2].upper() == "GT"
        # acceptor on + strand: pos is intron end e, motif at [e-2, e)
        if pos - 2 < 0 or pos > n:
            return False
        return genome_seq[pos - 2:pos].upper() == "AG"
    # '-' strand
    if kind == "donor":
        # donor on - strand: pos is intron end e, motif at [e-2, e) == "AC"
        if pos - 2 < 0 or pos > n:
            return False
        return genome_seq[pos - 2:pos].upper() == "AC"
    # acceptor on - strand: pos is intron start s, motif at [s, s+2) == "CT"
    if pos < 0 or pos + 2 > n:
        return False
    return genome_seq[pos:pos + 2].upper() == "CT"


def _pick_consensus(
    cluster_positions: List[int],
    kind: str,
    strand: str,
    genome_seq: str,
) -> int:
    """Pick a consensus position for a cluster of donor or acceptor sites.

    Priority order (per the junction-snap algorithm):
      1. Most-supported position whose flanking motif is canonical.
      2. Most-supported position overall (mode).
      3. Median (tie-break — applied implicitly because Counter.most_common is
         stable; we additionally pick the median when ties remain).

    Args:
        cluster_positions: Positions inside one cluster (duplicates allowed,
            their counts encode read support).
        kind: "donor" or "acceptor".
        strand: "+" or "-".
        genome_seq: Full chromosome sequence string.

    Returns:
        Chosen consensus coordinate.
    """
    counts = Counter(cluster_positions)
    # Group positions by their support count, descending.
    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))

    # 1. Highest-support canonical-motif position.
    top_support = ordered[0][1]
    top_tier = [p for p, c in ordered if c == top_support]
    canonical_top = [p for p in top_tier if _has_canonical_motif(p, kind, strand, genome_seq)]
    if canonical_top:
        return min(canonical_top, key=lambda p: (abs(p - _median(top_tier)), p))

    # Look further down support tiers for any canonical motif before falling
    # back to mode/median; this lets a slightly-less-supported canonical site
    # beat noisy non-canonical positions at higher support.
    for p, _c in ordered:
        if _has_canonical_motif(p, kind, strand, genome_seq):
            # Only override the mode if its support is at least half the mode's.
            if counts[p] * 2 >= top_support:
                return p
            break

    # 2. Mode (most-supported overall). Among ties, take the median.
    if len(top_tier) == 1:
        return top_tier[0]
    return _median(top_tier)


def _median(values: List[int]) -> int:
    """Return the lower median of ``values`` (matches existing 3'-end logic)."""
    s = sorted(values)
    return s[len(s) // 2]


def snap_junctions(
    read_chains: List[Tuple[Dict, IntronChain]],
    genome_seq: str,
    strand: str,
    snap_bp: int = 6,
    min_support: int = 2,
) -> List[Tuple[Dict, IntronChain]]:
    """Snap noisy splice-site coordinates onto per-cluster consensus positions.

    Algorithm:
      1. Pool every read's donor (intron start) and acceptor (intron end)
         positions independently.
      2. Greedy 1D-cluster each pool with ``snap_bp`` as the join gap.
      3. Pick a consensus per cluster, preferring canonical splice motifs
         (see :func:`_pick_consensus`).
      4. Rewrite each read's chain to the consensus coordinates. Clusters
         with fewer than ``min_support`` reads keep their original positions
         to avoid over-snapping rare junctions.
      5. Drop any rewritten intron that becomes degenerate (``s >= e``).

    The function does not mutate the input read dicts; it returns a new list
    of ``(read_dict, IntronChain)`` tuples with replaced (or unchanged) chains.

    Args:
        read_chains: Per-read ``(read_dict, IntronChain)`` pairs.
        genome_seq: Full chromosome sequence string used for motif lookup.
        strand: "+" or "-".
        snap_bp: Maximum within-cluster gap. Must be > 0 to do any work.
        min_support: Minimum number of reads in a cluster before it is
            snapped to its consensus.

    Returns:
        A list parallel to ``read_chains`` with (possibly) rewritten chains.
    """
    if snap_bp <= 0 or not read_chains:
        return list(read_chains)

    donors: List[int] = []
    acceptors: List[int] = []
    for _rd, chain in read_chains:
        for s, e in chain.introns:
            donor_pos = s if strand == "+" else e
            acceptor_pos = e if strand == "+" else s
            donors.append(donor_pos)
            acceptors.append(acceptor_pos)

    donor_clusters = _cluster_positions(donors, snap_bp)
    acceptor_clusters = _cluster_positions(acceptors, snap_bp)

    donor_map: Dict[int, int] = {}
    for cluster in donor_clusters:
        if len(cluster) < min_support:
            for p in cluster:
                donor_map[p] = p
            continue
        consensus = _pick_consensus(cluster, "donor", strand, genome_seq)
        for p in cluster:
            donor_map[p] = consensus

    acceptor_map: Dict[int, int] = {}
    for cluster in acceptor_clusters:
        if len(cluster) < min_support:
            for p in cluster:
                acceptor_map[p] = p
            continue
        consensus = _pick_consensus(cluster, "acceptor", strand, genome_seq)
        for p in cluster:
            acceptor_map[p] = consensus

    snapped: List[Tuple[Dict, IntronChain]] = []
    for rd, chain in read_chains:
        new_introns: List[Tuple[int, int]] = []
        for s, e in chain.introns:
            if strand == "+":
                new_s = donor_map.get(s, s)
                new_e = acceptor_map.get(e, e)
            else:
                new_s = acceptor_map.get(s, s)
                new_e = donor_map.get(e, e)
            if new_e > new_s:
                new_introns.append((new_s, new_e))
            # else: degenerate after snap, drop the intron
        snapped.append((rd, IntronChain(introns=tuple(new_introns))))

    return snapped


def _find_canonical_motif_positions(
    genome_seq: str,
    center: int,
    kind: str,
    strand: str,
    search_bp: int,
) -> List[int]:
    """All positions in [center-search_bp, center+search_bp] whose 2 bp motif
    matches the canonical donor (GT/AC) or acceptor (AG/CT) for ``strand``.

    The motif convention is the same as :func:`_has_canonical_motif`.
    """
    out: List[int] = []
    for delta in range(-search_bp, search_bp + 1):
        pos = center + delta
        if _has_canonical_motif(pos, kind, strand, genome_seq):
            out.append(pos)
    return out


def expand_canonical_chain_alternatives(
    read_chains: List[Tuple[Dict, IntronChain]],
    genome_seq: str,
    strand: str,
    search_bp: int = 4,
    max_chains_per_read: int = 16,
) -> Dict[str, List[IntronChain]]:
    """For each read, build the set of intron-chain alternatives obtained by
    moving each junction's donor and/or acceptor onto any canonical GT-AG
    motif within ``search_bp`` of the read's CIGAR-derived position.

    The original chain is always retained (so reads whose alignments already
    sit on the true position do not lose strict matches). Junctions without
    any nearby canonical motif also fall back to the original position.

    The total number of chains per read is bounded by ``max_chains_per_read``:
    if the per-intron Cartesian product would exceed the cap, the read keeps
    only its original chain (no expansion).

    Args:
        read_chains: Per-read ``(read_dict, IntronChain)`` pairs.
        genome_seq: Full chromosome sequence string for motif lookup.
        strand: "+" or "-".
        search_bp: Window radius (bp) around each junction position.
        max_chains_per_read: Hard cap on the per-read alternative count.

    Returns:
        Mapping ``query_name -> [IntronChain, ...]`` (each list non-empty;
        contains at least the original chain).
    """
    import itertools

    out: Dict[str, List[IntronChain]] = {}
    if search_bp <= 0:
        return out

    for rd, chain in read_chains:
        qname = rd.get("query_name")
        if qname is None:
            continue
        if not chain.introns:
            out[qname] = [chain]
            continue

        per_intron_alts: List[List[Tuple[int, int]]] = []
        for s, e in chain.introns:
            donor_center = s if strand == "+" else e
            acceptor_center = e if strand == "+" else s

            donor_alts = _find_canonical_motif_positions(
                genome_seq, donor_center, "donor", strand, search_bp
            )
            acceptor_alts = _find_canonical_motif_positions(
                genome_seq, acceptor_center, "acceptor", strand, search_bp
            )
            # Always include the original positions so reads whose CIGAR is
            # already correct do not lose their strict match.
            if donor_center not in donor_alts:
                donor_alts.append(donor_center)
            if acceptor_center not in acceptor_alts:
                acceptor_alts.append(acceptor_center)

            intron_alts: List[Tuple[int, int]] = []
            seen = set()
            for d in donor_alts:
                for a in acceptor_alts:
                    ns, ne = (d, a) if strand == "+" else (a, d)
                    if ne > ns and (ns, ne) not in seen:
                        seen.add((ns, ne))
                        intron_alts.append((ns, ne))
            if not intron_alts:
                intron_alts = [(s, e)]
            per_intron_alts.append(intron_alts)

        # Quick size check before exploding the Cartesian product.
        total = 1
        for ia in per_intron_alts:
            total *= len(ia)
            if total > max_chains_per_read:
                break
        if total > max_chains_per_read:
            out[qname] = [chain]
            continue

        chains_for_read: List[IntronChain] = []
        seen_chains = set()
        for combo in itertools.product(*per_intron_alts):
            key = tuple(combo)
            if key in seen_chains:
                continue
            seen_chains.add(key)
            chains_for_read.append(IntronChain(introns=key))
        if not chains_for_read:
            chains_for_read = [chain]
        out[qname] = chains_for_read

    return out


def group_reads_by_three_prime_and_intron_chain(
    read_dicts: List[Dict],
    strand: str,
    threshold: int = 24,
    chain_override: Optional[Dict[str, IntronChain]] = None,
    multi_chain_override: Optional[Dict[str, List[IntronChain]]] = None,
) -> Dict[Tuple[int, IntronChain], List[Dict]]:
    """Group reads by clustered 3' end position and intron chain.

    1. Extract 3' position and intron chain per read.
    2. Sort by 3' position.
    3. Cluster 3' positions within *threshold* bp.
    4. Sub-group each cluster by intron chain.
    5. Use consensus (median) 3' position per cluster.

    Args:
        read_dicts: List of alignment dicts from BamReader.alignment_to_dict().
        strand: '+' or '-'.
        threshold: Maximum distance between 3' ends in a cluster.
        chain_override: Optional mapping ``query_name -> IntronChain`` that
            replaces the CIGAR-derived chain (used for junction-snap output).
        multi_chain_override: Optional mapping ``query_name -> [IntronChain, ...]``
            that produces multiple grouping entries per read (used for the
            canonical-motif alternative expansion). Takes precedence over
            ``chain_override`` when both are provided for the same read.

    Returns:
        Dict mapping (consensus_3prime, IntronChain) -> list of read dicts.
    """
    # Annotate each read with its 3' pos and intron chain
    annotated: List[Tuple[int, IntronChain, Dict]] = []
    for rd in read_dicts:
        three_prime = extract_three_prime_pos(rd)
        if three_prime is None:
            continue
        cigar = rd.get("cigartuples")
        ref_start = rd.get("reference_start")
        if cigar is None or ref_start is None:
            continue
        qname = rd.get("query_name")
        chains_for_this_read: Optional[List[IntronChain]] = None
        if multi_chain_override is not None and qname is not None \
                and qname in multi_chain_override:
            chains_for_this_read = multi_chain_override[qname]
        elif chain_override is not None and qname is not None \
                and qname in chain_override:
            chains_for_this_read = [chain_override[qname]]
        if chains_for_this_read is None:
            chains_for_this_read = [extract_intron_chain(cigar, ref_start)]
        for chain in chains_for_this_read:
            annotated.append((three_prime, chain, rd))

    if not annotated:
        return {}

    # Sort by 3' position
    annotated.sort(key=lambda x: x[0])

    # Cluster by 3' proximity
    clusters: List[List[Tuple[int, IntronChain, Dict]]] = []
    current_cluster: List[Tuple[int, IntronChain, Dict]] = [annotated[0]]

    for item in annotated[1:]:
        if item[0] - current_cluster[-1][0] <= threshold:
            current_cluster.append(item)
        else:
            clusters.append(current_cluster)
            current_cluster = [item]
    clusters.append(current_cluster)

    # Sub-group each cluster by intron chain, compute consensus 3'
    groups: Dict[Tuple[int, IntronChain], List[Dict]] = {}
    for cluster in clusters:
        # Consensus 3' = median of cluster positions
        positions = [x[0] for x in cluster]
        consensus_3prime = sorted(positions)[len(positions) // 2]

        by_chain: Dict[IntronChain, List[Dict]] = defaultdict(list)
        for _, chain, rd in cluster:
            by_chain[chain].append(rd)

        for chain, reads in by_chain.items():
            groups[(consensus_3prime, chain)] = reads

    return groups


def pick_representative_read(read_dicts: List[Dict]) -> Dict:
    """Pick the longest read as the representative for a group.

    Args:
        read_dicts: List of alignment dicts.

    Returns:
        The read dict with the largest query_length.
    """
    return max(read_dicts, key=lambda rd: rd.get("query_length", 0))
