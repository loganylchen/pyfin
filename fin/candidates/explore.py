"""Explore step: family-local intron-graph path enumeration (additive candidate expansion).

Runs AFTER clustering (fin.candidates.chain_cluster.cluster_families) and BEFORE any fold. It is a
pure RECALL step: it never drops or rewrites anything -- it ADDS assembled full-length candidates that
no single read produced, by tiling read-phased junction edges across reads within one family. The
downstream re-mapping / selection does the TP/FP discrimination, so generating many candidates here is
acceptable by design.

Why it exists: a read's exact chain is only the introns THAT read spanned. dRNA 5' degradation means a
long transcript's reads are each partial (read A = [I1,I2], read B = [I2,I3]), so the full chain
[I1,I2,I3] is on NO single read and is therefore absent from the exact-chain candidate set. The union
of read-phased edges (I1->I2 from A, I2->I3 from B) DOES contain the full path; walking it reconstructs
the missing candidate. This is the dominant de-novo recall gap (measured ~74% of missed full-length-
supported truth are this cross-read-truncation case).

Design (agreed with the user + Codex gpt-5.6-sol):
- NODES = every distinct OBSERVED intron coordinate in the family (wobble variants kept as SEPARATE
  nodes -- NO consensus collapse, NO coordinate rewrite).
- EDGES = ordered pairs (A, B) that appear CONSECUTIVELY in some variant chain (i.e. co-observed on the
  reads that produced that variant). Read-phased only: a path's every adjacency is read-backed. Genomic
  order + positive exon between introns make the graph a DAG (no cycles, no same-slot double-use).
- ENUMERATE maximal source->sink paths with a deterministic BOTTLENECK-FIRST beam (bounded), keeping
  N=1 edges (precision is a downstream gate, default 2). Emit only chains that are NOT already a variant
  (baseline exact variants are kept by the caller). Cap the assembled additions per family.

CAVEAT (Codex): read-phased edges prevent LOCAL fabrication but NOT global recombination at a shared
node -- reads A-B-X-D and A-C-X-E also yield A-B-X-E / A-C-X-D whose every edge is observed but whose
whole chain is not. That is accepted for an additive, downstream-filtered generation step; the
per-path ``min_edge_support`` is exposed so downstream can weight it. Richer recombination/handoff
provenance is a follow-up.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple

Chain = Tuple[Tuple[int, int], ...]


@dataclass
class ExploredPath:
    """One assembled candidate chain (not observed whole on any single read).

    ``min_edge_support`` = the bottleneck (minimum over the path's edges of the number of reads that
    co-observed that adjacency) -- the key quality signal for downstream weighting. ``reads`` = union
    of the reads of every variant that is a contiguous sub-chain of this path (the reads that tile it).
    """
    chain: Chain
    min_edge_support: int
    n_edges: int
    reads: Set[str] = field(default_factory=set)


def _build_graph(family):
    """Return (node_support, edge_support, succ, pred) for a family's read-derived variants.

    A read supports every intron in its variant chain and every consecutive pair. Reads are disjoint
    across variants (a read produces exactly one exact chain), so summing per-variant read counts is a
    true union count. Only genomically valid edges (positive exon between the two introns) are kept, so
    the graph is a DAG in genomic order.
    """
    node_support: Counter = Counter()
    edge_support: Counter = Counter()
    succ: Dict[Tuple[int, int], Set[Tuple[int, int]]] = defaultdict(set)
    pred: Dict[Tuple[int, int], Set[Tuple[int, int]]] = defaultdict(set)
    for v in family.variants:
        r = len(family.variant_reads.get(v, ()) or ())
        for intr in v:
            node_support[intr] += r
        for a, b in zip(v, v[1:]):
            if b[0] <= a[1]:            # require a real exon between introns (genomic order)
                continue
            edge_support[(a, b)] += r
            succ[a].add(b)
            pred[b].add(a)
    return node_support, edge_support, succ, pred


def explore_family_paths(
    family,
    *,
    beam_width: int = 256,
    max_new: int = 512,
) -> List[ExploredPath]:
    """Enumerate assembled full-length candidates for one :class:`ChainFamily`.

    Builds the read-phased intron DAG, beam-enumerates maximal source->sink paths (deterministic,
    bottleneck-first), and returns the paths that are NOT already a variant of the family (the caller
    keeps all exact variants separately). Additive and reversible: with the graph empty or single-node
    families, returns []. Never mutates the family.
    """
    if max_new <= 0:
        return []
    beam_width = max(1, beam_width)
    node_support, edge_support, succ, pred = _build_graph(family)
    if not node_support:
        return []
    variants_set = set(family.variants)
    nodes = set(node_support)
    sources = sorted(n for n in nodes if not pred.get(n))
    if not sources:
        return []  # a pure cycle is impossible (DAG); no source means nothing to walk

    def _edge_scores(path: Tuple) -> Tuple[int, int]:
        if len(path) < 2:
            return (0, 0)
        e = [edge_support[(path[i], path[i + 1])] for i in range(len(path) - 1)]
        return (min(e), sum(e))            # bottleneck first, then cumulative

    completed: List[Tuple] = []
    # Termination: the DAG is finite with no node revisits, so the walk ends after at most
    # len(nodes) layers. The initial frontier is every source S (bounded by the family's variant
    # count) and is NOT truncated -- one-node paths all score (0, 0), so a beam cut here would
    # deterministically drop a high-support source before it is ever explored (Codex gpt-5.6-sol).
    # Work: the first layer is O(S x max_outdegree) transient and may retain O(S) completions;
    # every later layer prunes the frontier to beam_width, so it is O(beam_width x max_outdegree)
    # transient, giving total `completed` retention O(S + beam_width x depth). We also deliberately
    # do NOT early-exit on a completion budget: a path's bottleneck-first rank is only known once it
    # reaches its sink, so cutting early on a count of completions can drop the STRONGEST assembled
    # chain in favour of shallow, weak ones that reached a sink first. Rank + the max_new cap are
    # applied once at the end instead.
    frontier: List[Tuple] = [(s,) for s in sources]
    while frontier:
        nxt: List[Tuple] = []
        for path in frontier:
            outs = sorted(succ.get(path[-1], ()))
            if not outs:
                completed.append(path)     # reached a sink -> a maximal path
                continue
            for b in outs:
                if b in path:              # never revisit a node (DAG guarantees, guard anyway)
                    continue
                nxt.append(path + (b,))
        if not nxt:
            break
        # deterministic prune: stable sort by coords, then by (bottleneck, cumulative) desc.
        nxt.sort()
        nxt.sort(key=_edge_scores, reverse=True)
        frontier = nxt[:beam_width]

    # dedup + keep only NEW multi-intron assembled chains, rank by support then coords.
    seen: Set[Chain] = set()
    out: List[ExploredPath] = []
    for path in sorted(completed, key=lambda p: (-_edge_scores(p)[0], -_edge_scores(p)[1], p)):
        chain: Chain = tuple(path)
        if len(chain) < 2 or chain in variants_set or chain in seen:
            continue
        seen.add(chain)
        mn, _ = _edge_scores(path)
        reads: Set[str] = set()
        for v in family.variants:                     # reads of variants that tile this path
            if _is_subchain(v, chain):
                reads |= (family.variant_reads.get(v, set()) or set())
        out.append(ExploredPath(chain=chain, min_edge_support=mn,
                                n_edges=len(chain) - 1, reads=reads))
        if len(out) >= max_new:
            break
    return out


def _is_subchain(short: Chain, long_: Chain) -> bool:
    """True iff ``short`` is a contiguous exact sub-sequence of ``long_`` (or equal)."""
    n, m = len(long_), len(short)
    if m == 0 or m > n:
        return False
    return any(long_[o:o + m] == short for o in range(n - m + 1))
