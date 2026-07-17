"""Generation-side intron-chain clustering.

Cluster-first candidate generation. Pool all reads, group by their EXACT (raw,
un-snapped) intron chain, then:

  1. COLLAPSE exact sub-chains  -- a chain that is a strict contiguous sub-chain of a
     longer chain with IDENTICAL coordinates folds into it (reads pooled, chain gone).
  2. CLUSTER the survivors      -- union-find groups the remaining distinct chains into
     one cluster when related by any of:
       * WOBBLE      : same intron count, every junction within ``wobble_bp``,
       * CASSETTE    : K vs K-1 chains differing by one small (< ``cassette_max_exon_bp``)
                       skipped exon (minimap2's small-exon-skip artefact),
       * CONTAINMENT : a shorter chain that is a contiguous sub-chain of a longer one
                       within ``wobble_bp`` per junction (the non-exact leftovers that
                       step 1 did not fold).

A cluster KEEPS all its member candidates (distinct chains + pooled reads); only exact
sub-chains disappear. NO consensus snapping (coordinates are never rewritten). The
cluster is the downstream unit: reads attach to a cluster instead of splitting across
near-duplicates, and which members are genuine isoforms is decided later by the
per-cluster computation. The caller runs with canonical_search_bp = 0 (no shadows).

Mono-exon chains never join a structural cluster; each unique mono chain is a singleton.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple

from fin.candidates.dataclasses import IntronChain

Chain = Tuple[Tuple[int, int], ...]


@dataclass
class GenCandidate:
    """One candidate inside a cluster: a distinct intron chain + its pooled reads."""
    chain: IntronChain
    read_ids: Set[str] = field(default_factory=set)

    @property
    def num_introns(self) -> int:
        return len(self.chain.introns)


@dataclass
class ChainCluster:
    """A structural family of candidate chains sharing their reads' assignment.

    ``members`` are the distinct surviving candidates (exact sub-chains already folded
    away). ``read_ids`` is the union of all members' reads. ``representative`` is the
    longest member (convenience; NOT a snapped consensus -- original coordinates)."""
    members: List[GenCandidate]
    read_ids: Set[str] = field(default_factory=set)

    @property
    def representative(self) -> GenCandidate:
        return max(self.members,
                   key=lambda m: (m.num_introns, len(m.read_ids), m.chain.introns))


# ---------------------------------------------------------------------------
# chain-relation predicates (single interval/strand: no chrom/strand bucketing)
# ---------------------------------------------------------------------------
def _exact_subchain(short: Chain, long_: Chain) -> bool:
    """``short`` is a strict EXACT (identical coords) contiguous sub-chain of ``long_``."""
    n, m = len(long_), len(short)
    if m == 0 or m >= n:
        return False
    return any(long_[off:off + m] == short for off in range(n - m + 1))


def _wobble(a: Chain, b: Chain, bp: int) -> bool:
    """Same intron count, every donor/acceptor within ``bp``."""
    if len(a) != len(b):
        return False
    return all(abs(s1 - s2) <= bp and abs(e1 - e2) <= bp
               for (s1, e1), (s2, e2) in zip(a, b))


def _contains(short: Chain, long_: Chain, bp: int) -> bool:
    """``short`` is a strictly-shorter contiguous sub-chain of ``long_`` within ``bp``."""
    n, m = len(long_), len(short)
    if m == 0 or m >= n:
        return False
    for off in range(n - m + 1):
        if all(abs(short[k][0] - long_[off + k][0]) <= bp and
               abs(short[k][1] - long_[off + k][1]) <= bp for k in range(m)):
            return True
    return False


def _cassette(a: Chain, b: Chain, bp: int, max_exon_bp: int) -> bool:
    """``a`` (K introns) is a cassette-skip sibling of ``b`` (K-1 introns): at one
    position ``i`` the two introns ``a[i], a[i+1]`` replace ``b[i]``'s span (outer
    donor/acceptor within bp), the skipped exon between them is < ``max_exon_bp``, and
    every other intron matches within bp."""
    if max_exon_bp <= 0 or len(a) != len(b) + 1:
        return False
    for i in range(len(a) - 1):
        if any(abs(a[k][0] - b[k][0]) > bp or abs(a[k][1] - b[k][1]) > bp
               for k in range(i)):
            continue
        if abs(a[i][0] - b[i][0]) > bp or abs(a[i + 1][1] - b[i][1]) > bp:
            continue
        rest = len(b) - i - 1
        if any(abs(a[i + 2 + k][0] - b[i + 1 + k][0]) > bp or
               abs(a[i + 2 + k][1] - b[i + 1 + k][1]) > bp for k in range(rest)):
            continue
        cas = a[i + 1][0] - a[i][1]
        if 0 < cas < max_exon_bp:
            return True
    return False


def _related(a: Chain, b: Chain, wobble_bp: int, cassette_max_exon_bp: int) -> bool:
    """True iff ``a`` and ``b`` should share a cluster (wobble / cassette / containment)."""
    if _wobble(a, b, wobble_bp):
        return True
    lo, hi = (a, b) if len(a) < len(b) else (b, a)
    if len(lo) < len(hi):
        if _cassette(hi, lo, wobble_bp, cassette_max_exon_bp):
            return True
        if _contains(lo, hi, wobble_bp):
            return True
    return False


def cluster_read_chains(
    read_chains: List[Tuple[Dict, IntronChain]],
    *,
    wobble_bp: int = 6,
    cassette_max_exon_bp: int = 70,
) -> List[ChainCluster]:
    """Build generation-side clusters from per-read chains (NO snap, 3' ignored).

    Args:
        read_chains: per-read ``(read_dict, IntronChain)`` (read_dict needs
            "query_name"; CIGAR-derived chains, used as-is).
        wobble_bp: per-junction tolerance for wobble / cassette / containment joins.
        cassette_max_exon_bp: max skipped-exon length for a cassette join (0 off).

    Returns:
        List of :class:`ChainCluster`. EVERY exact contiguous sub-chain is folded
        (unconditionally, no read-support guard) into its longest exact container --
        generation does pure structural collapse and does NOT decide whether a folded
        sub-chain is a real short isoform; that judgement is a SEPARATE downstream
        recovery step (e.g. read-end pileup), not done here. The survivors are unioned
        into clusters by wobble/cassette/containment; each unique mono chain is its own
        singleton cluster.
    """
    # 1. group reads by EXACT chain (no 3', no snap).
    by_chain: Dict[Chain, Set[str]] = {}
    for rd, chain in read_chains:
        rid = rd.get("query_name")
        if rid is None:
            continue
        by_chain.setdefault(tuple(chain.introns), set()).add(rid)

    # 2. collapse EVERY exact contiguous sub-chain into its longest exact container,
    #    UNCONDITIONALLY (no read-support guard). Whether a folded sub-chain is a
    #    genuine short isoform is decided later by a dedicated recovery step, never
    #    guessed here from read counts.
    chains_sorted = sorted(by_chain, key=lambda c: (-len(c), -len(by_chain[c]), c))
    kept: List[Chain] = []
    reads: Dict[Chain, Set[str]] = {}
    for c in chains_sorted:
        container = next((k for k in kept if _exact_subchain(c, k)), None)
        if container is None:
            kept.append(c)
            reads[c] = set(by_chain[c])
        else:
            reads[container] |= by_chain[c]

    # 3. union-find the survivors by wobble / cassette / containment (multi-exon only).
    multi = [c for c in kept if len(c) >= 1]
    mono = [c for c in kept if len(c) == 0]
    parent = {c: c for c in multi}

    def find(x: Chain) -> Chain:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(len(multi)):
        for j in range(i + 1, len(multi)):
            a, b = multi[i], multi[j]
            if find(a) != find(b) and _related(a, b, wobble_bp, cassette_max_exon_bp):
                parent[find(a)] = find(b)

    comps: Dict[Chain, List[Chain]] = {}
    for c in multi:
        comps.setdefault(find(c), []).append(c)

    clusters: List[ChainCluster] = []
    for comp in comps.values():
        members = [GenCandidate(IntronChain(introns=c), set(reads[c])) for c in comp]
        union: Set[str] = set()
        for c in comp:
            union |= reads[c]
        clusters.append(ChainCluster(members=members, read_ids=union))

    for c in mono:  # mono singletons
        clusters.append(ChainCluster(
            members=[GenCandidate(IntronChain(introns=c), set(reads[c]))],
            read_ids=set(reads[c])))

    return clusters
