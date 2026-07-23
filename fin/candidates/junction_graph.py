"""De-novo intron graph: assemble full-length transcript intron chains from
partial/truncated read chains.

Measured motivation (p00, no annotation): pyfin's dominant de-novo error is
truncation — candidates that are a contiguous sub-chain of a real transcript
(gffcompare class 'c', ~16% of output), because dRNA reads are 3'-biased so a
read covers only the 3' junctions. pyfin currently emits each read's exact chain
as a candidate, so a truncated read -> a truncated candidate.

This module pools all reads' junctions, clusters wobbled junctions to a
read-count consensus, builds a directed read-adjacency graph (nodes = consensus
junctions ordered by genomic coord; edges = "these two junctions are consecutive
in >= min_edge_reads reads"), and EXTENDS each read's chain through
UNAMBIGUOUSLY-supported edges to a maximal chain. Extension stops at genuine
branch points (a node with >1 supported successor), so it lengthens truncated
chains WITHOUT fabricating wrong junction combinations.

NOT annotation-snap: consensus coordinates and edges come only from reads.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from typing import Dict, List, Set, Tuple

Junction = Tuple[int, int]
Chain = Tuple[Junction, ...]


def cluster_junctions(chains: List[Chain], tol: int) -> Dict[Junction, Junction]:
    """Map every observed junction to a consensus junction. Greedy by read count:
    the most-observed junction seeds a centroid; a junction within ``tol`` bp on
    BOTH donor and acceptor maps to the first matching centroid. Deterministic
    (count desc, then coordinate)."""
    counts: Counter = Counter()
    for ch in chains:
        for j in ch:
            counts[j] += 1
    consensus: Dict[Junction, Junction] = {}
    centroids: List[Junction] = []
    for j in sorted(counts, key=lambda x: (-counts[x], x)):
        hit = None
        for cj in centroids:
            if abs(j[0] - cj[0]) <= tol and abs(j[1] - cj[1]) <= tol:
                hit = cj
                break
        if hit is None:
            centroids.append(j)
            consensus[j] = j
        else:
            consensus[j] = hit
    return consensus


def build_graph(cons_chains: List[Chain]) -> Tuple[Counter, Counter]:
    """Node + edge read-support counts from consensus-mapped chains.
    Edge (a, b) = a and b are consecutive junctions in a read's chain."""
    node_counts: Counter = Counter()
    edge_counts: Counter = Counter()
    for ch in cons_chains:
        for j in ch:
            node_counts[j] += 1
        for a, b in zip(ch, ch[1:]):
            edge_counts[(a, b)] += 1
    return node_counts, edge_counts


def anchored_5p_junctions(
    cons_chains: List[Chain],
    five_prime_pos: List[int],
    strand: str,
    tss_tol: int,
    min_tss_reads: int,
    tss_frac: float,
) -> Set[Junction]:
    """Return the set of 5'-terminal consensus junctions that sit at a real TSS.

    dRNA is 3'-biased: a truncated read of a LONG transcript stops at a RANDOM 5'
    position, but a COMPLETE read of a genuine short isoform stops at that
    isoform's real TSS. So among the reads whose 5'-most junction is ``J`` (they
    begin at ``J``), if a strong fraction pile their read-5'-ends within
    ``tss_tol`` bp of one position, ``J`` marks a real transcription start and a
    chain ending there must NOT be extended further 5'-ward (that would erase the
    short isoform). Empirically (p00) such peaks carry 40-90% of a locus's reads;
    degradation-only starts scatter and never cross ``tss_frac``.

    5'-most junction is ``chain[0]`` on '+' (5' = small coord) and ``chain[-1]``
    on '-' (5' = large coord).
    """
    by_term: Dict[Junction, List[int]] = defaultdict(list)
    for cc, p in zip(cons_chains, five_prime_pos):
        if not cc or p is None:
            continue
        term = cc[0] if strand == "+" else cc[-1]
        by_term[term].append(p)
    anchored: Set[Junction] = set()
    for term, positions in by_term.items():
        if len(positions) < min_tss_reads:
            continue
        best = max(sum(1 for q in positions if abs(q - p) <= tss_tol)
                   for p in positions)
        if best >= min_tss_reads and best >= tss_frac * len(positions):
            anchored.add(term)
    return anchored


def assemble_chains(
    chains: List[Chain],
    tol: int,
    min_edge_reads: int,
    five_prime_pos: List[int] = None,
    strand: str = "+",
    tss_tol: int = 20,
    min_tss_reads: int = 3,
    tss_frac: float = 0.4,
) -> Tuple[List[Chain], Dict[int, Chain]]:
    """Assemble maximal chains by unambiguous extension.

    Returns ``(cons_chains, extended_by_index)`` where ``cons_chains[i]`` is input
    chain ``i`` mapped to consensus junctions, and ``extended_by_index[i]`` is the
    maximal chain that chain ``i`` extends into (following only edges from a node
    that has exactly ONE supported successor / predecessor). Two input reads whose
    (possibly truncated) chains lie on the same maximal path map to the SAME
    extended chain, so downstream grouping by the extended chain merges them.

    5'-TSS brake: if ``five_prime_pos`` (per-input-chain read-5'-end genomic
    positions) is given, a chain whose 5'-terminal junction is a real TSS
    (:func:`anchored_5p_junctions`) is NOT extended in the 5' direction, so
    genuine short isoforms contained inside a longer transcript survive instead of
    being merged away. Without ``five_prime_pos`` behaviour is byte-identical to
    the un-braked assembly (default-off contract).
    """
    cons = cluster_junctions(chains, tol)
    cons_chains: List[Chain] = [tuple(cons[j] for j in ch) for ch in chains]
    _, edge_counts = build_graph(cons_chains)

    succ: Dict[Junction, List[Junction]] = defaultdict(list)
    pred: Dict[Junction, List[Junction]] = defaultdict(list)
    for (a, b), n in edge_counts.items():
        if n >= min_edge_reads:
            succ[a].append(b)
            pred[b].append(a)

    anchored: Set[Junction] = set()
    if five_prime_pos is not None:
        anchored = anchored_5p_junctions(
            cons_chains, five_prime_pos, strand, tss_tol, min_tss_reads, tss_frac
        )

    def unique_succ(j: Junction):
        s = succ.get(j)
        return s[0] if s and len(s) == 1 else None

    def unique_pred(j: Junction):
        p = pred.get(j)
        return p[0] if p and len(p) == 1 else None

    cache: Dict[Chain, Chain] = {}

    def extend(cc: Chain) -> Chain:
        if cc in cache:
            return cache[cc]
        if not cc:
            cache[cc] = cc
            return cc
        chain = list(cc)
        seen: Set[Junction] = set(chain)
        # 5' terminus of THIS read's chain (fixed; 3'-ward extension won't move it)
        five_term = cc[0] if strand == "+" else cc[-1]
        braked = five_term in anchored
        # right extension is 5'-ward on '-' strand, 3'-ward on '+'
        block_right = braked and strand == "-"
        # left extension is 5'-ward on '+' strand, 3'-ward on '-'
        block_left = braked and strand == "+"
        while not block_right:
            nxt = unique_succ(chain[-1])
            if nxt is None or nxt in seen:
                break
            chain.append(nxt)
            seen.add(nxt)
        while not block_left:
            prv = unique_pred(chain[0])
            if prv is None or prv in seen:
                break
            chain.insert(0, prv)
            seen.add(prv)
        out = tuple(chain)
        cache[cc] = out
        return out

    extended_by_index: Dict[int, Chain] = {
        i: extend(cc) for i, cc in enumerate(cons_chains)
    }
    return cons_chains, extended_by_index
