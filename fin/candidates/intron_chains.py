"""Intron chain extraction and read grouping for candidate discovery."""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from fin.candidates.dataclasses import IntronChain
from fin.io.interval_manager import extract_three_prime_pos, extract_strand_from_read

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


def group_reads_by_three_prime_and_intron_chain(
    read_dicts: List[Dict],
    strand: str,
    threshold: int = 24,
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
        chain = extract_intron_chain(cigar, ref_start)
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
