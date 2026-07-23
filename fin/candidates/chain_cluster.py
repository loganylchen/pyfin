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
    """One candidate inside a cluster: a distinct intron chain + its pooled reads.

    ``folded`` holds the SHADOW sub-chains folded into this candidate (each an exact
    contiguous sub-chain, with its OWN reads). They are dormant provenance for the
    downstream short-isoform recovery step: their truncation boundary is a proposed
    TSS/polyA site and their reads are the direct endpoint evidence. ``read_ids``
    already INCLUDES every shadow's reads (pooled), so shadows never add read mass."""
    chain: IntronChain
    read_ids: Set[str] = field(default_factory=set)
    folded: List["GenCandidate"] = field(default_factory=list)

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


def _subchain_offset(short: Chain, long_: Chain) -> int:
    """Offset of the first exact ``short`` occurrence in ``long_`` (-1 if none). The
    introns of ``long_`` outside ``[off, off+len(short))`` are the container's EXTRA
    introns relative to ``short`` -- positions where ``short`` splices differently."""
    n, m = len(long_), len(short)
    if m == 0 or m >= n:
        return -1
    for off in range(n - m + 1):
        if long_[off:off + m] == short:
            return off
    return -1


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


def _monoexon_in_exon(s: int, e: int, chain: Chain, lo: int, hi: int) -> bool:
    """``[s, e)`` (a single-exon read) lies wholly inside ONE exon of a multi-exon
    candidate whose exonic footprint is ``[lo, hi)`` and whose introns are ``chain``
    (each ``[d, a)`` half-open). True iff the read is within the footprint AND overlaps
    no intron -- i.e. it is a fragment of that candidate's exonic sequence, not an
    intronic feature (a different gene) and not something that crosses a junction."""
    if s < lo or e > hi:
        return False
    for d, a in chain:                       # overlap of [s,e) and intron [d,a)
        if s < a and d < e:
            return False
    return True


def cluster_read_chains(
    read_chains: List[Tuple[Dict, IntronChain]],
    *,
    wobble_bp: int = 6,
    cassette_max_exon_bp: int = 70,
    fold_monoexon_contained: bool = False,
    fold_span_guard: bool = False,
) -> List[ChainCluster]:
    """Build generation-side clusters from per-read chains (NO snap, 3' ignored).

    Args:
        read_chains: per-read ``(read_dict, IntronChain)`` (read_dict needs
            "query_name"; CIGAR-derived chains, used as-is).
        wobble_bp: per-junction tolerance for wobble / cassette / containment joins.
        cassette_max_exon_bp: max skipped-exon length for a cassette join (0 off).
        fold_monoexon_contained: when True, a single-exon (chain-less) read whose
            aligned span lies wholly inside one exon of a multi-exon candidate is folded
            into that candidate (its read joins the candidate's support) instead of
            emitting a standalone mono candidate. This absorbs 5'/3' degradation
            fragments of a spliced transcript. Reads inside an INTRON (a different gene)
            or not contained in any multi-exon candidate stay as mono candidates. Off by
            default (byte-identical legacy behaviour).
        fold_span_guard: when True, a read only folds into an exact-sub-chain container
            if its aligned span does NOT run exonically across one of the container's
            EXTRA introns (introns present in the container but absent from the read's
            chain). A read that spans such an intron contradicts it (retained-intron or
            alternative isoform) and is kept as its own candidate instead of being
            absorbed. Off by default (byte-identical legacy behaviour).

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

    # Read spans (0-based [start, end)) -- needed by the span-guarded fold and the
    # mono-exon fold. Built once, only when a span-aware step is enabled.
    read_span: Dict[str, Tuple[int, int]] = {}
    if fold_span_guard or fold_monoexon_contained:
        for rd, _chain in read_chains:
            rid = rd.get("query_name")
            rs = rd.get("reference_start")
            re_ = rd.get("reference_end")
            if rid is not None and rs is not None and re_ is not None:
                read_span[rid] = (rs, re_)

    # 2. collapse every exact contiguous sub-chain into its longest exact container.
    #    Default: UNCONDITIONAL (no read-support guard). With fold_span_guard, a read
    #    folds ONLY if its aligned span does not run exonically across one of the
    #    container's EXTRA introns (retained-intron / alternative isoform reads are kept
    #    as their own candidate instead of being absorbed). Whether a *consistent*
    #    folded sub-chain is a genuine short isoform is still a SEPARATE downstream
    #    recovery step (read-end pileup), never guessed here from read counts.
    chains_sorted = sorted(by_chain, key=lambda c: (-len(c), -len(by_chain[c]), c))
    kept: List[Chain] = []
    reads: Dict[Chain, Set[str]] = {}
    folded_of: Dict[Chain, List[Tuple[Chain, Set[str]]]] = {}  # container -> shadows
    for c in chains_sorted:
        container = next((k for k in kept if _exact_subchain(c, k)), None)
        if container is None:
            kept.append(c)
            reads[c] = set(by_chain[c])
            continue
        if not fold_span_guard:
            reads[container] |= by_chain[c]
            folded_of.setdefault(container, []).append((c, set(by_chain[c])))
            continue
        # span guard: split c's reads into those consistent with the container (fold)
        # and those exonic across one of the container's extra introns (keep as c).
        off = _subchain_offset(c, container)
        extra = container[:off] + container[off + len(c):]
        fold_r: Set[str] = set()
        keep_r: Set[str] = set()
        for rid in by_chain[c]:
            sp = read_span.get(rid)
            if sp is not None and any(sp[0] <= d and sp[1] >= a for d, a in extra):
                keep_r.add(rid)            # exonic across an extra intron -> distinct
            else:
                fold_r.add(rid)
        if fold_r:
            reads[container] |= fold_r
            folded_of.setdefault(container, []).append((c, fold_r))
        if keep_r:
            kept.append(c)
            reads[c] = keep_r

    # 3. union-find the survivors by wobble / cassette / containment (multi-exon only).
    multi = [c for c in kept if len(c) >= 1]
    mono = [c for c in kept if len(c) == 0]

    # 3b. (optional) fold single-exon reads that are wholly inside a multi-exon
    #     candidate's exon into that candidate -- 5'/3' degradation fragments of a
    #     spliced transcript, otherwise emitted as standalone mono FPs. Per-read (the
    #     mono chain pools every chain-less read), deterministic container choice.
    if fold_monoexon_contained and mono and multi:
        multi_span: Dict[Chain, Tuple[int, int]] = {}
        for c in multi:
            pts = [read_span[r] for r in reads[c] if r in read_span]
            if pts:
                multi_span[c] = (min(p[0] for p in pts), max(p[1] for p in pts))
        containers = [c for c in multi if c in multi_span]
        for mc in mono:
            remaining: Set[str] = set()
            for rid in reads[mc]:
                sp = read_span.get(rid)
                if sp is None:
                    remaining.add(rid)
                    continue
                s, e = sp
                hits = [c for c in containers
                        if _monoexon_in_exon(s, e, c, *multi_span[c])]
                if not hits:
                    remaining.add(rid)
                    continue
                best = max(hits, key=lambda c: (len(reads[c]), len(c), c))
                reads[best].add(rid)
                folded_of.setdefault(best, []).append(((), {rid}))
            reads[mc] = remaining
        mono = [c for c in mono if reads[c]]

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

    def _mk(c: Chain) -> GenCandidate:
        """Member candidate with its folded exact-sub-chain shadows attached."""
        return GenCandidate(
            IntronChain(introns=c), set(reads[c]),
            folded=[GenCandidate(IntronChain(introns=fc), set(fr))
                    for fc, fr in folded_of.get(c, [])])

    clusters: List[ChainCluster] = []
    for comp in comps.values():
        members = [_mk(c) for c in comp]
        union: Set[str] = set()
        for c in comp:
            union |= reads[c]
        clusters.append(ChainCluster(members=members, read_ids=union))

    for c in mono:  # mono singletons
        clusters.append(ChainCluster(members=[_mk(c)], read_ids=set(reads[c])))

    return clusters
