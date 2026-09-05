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
from typing import Dict, List, Optional, Set, Tuple

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


@dataclass
class ChainFamily:
    """A clustering-only family of related MULTI-EXON structural variants.

    ``variants`` are the distinct exact intron chains (structures) grouped together by
    wobble / cassette / containment -- NO coordinates are merged or snapped, and every
    variant is KEPT (no fold). ``read_pool`` is the union of all reads whose exact chain is
    in this family: once a family forms, reads belong to the FAMILY, not to an individual
    variant -- which variant a read truly supports is re-decided later, by evidence, in the
    per-cluster assignment. ``variant_reads`` is the per-variant discovery provenance (which
    exact chain each read produced), kept only for choosing a representative and for the
    later path-completion step; it is NOT the final read ownership.

    ``gtf_members`` records which of this family's ``variants`` are (also) annotated transcripts:
    each is ``(gtf_id, intron_chain)``. An annotated chain is an ORDINARY variant of the family -- it
    is a key in ``variants``, it takes part in the union-find exactly like a read-derived chain (it can
    join and can bridge), and it takes part in :func:`collapse` as a container or as a foldee. It
    carries no discovery reads of its own: ``variant_reads`` for a pure-GTF chain is empty and it is
    never a ``read_pool`` source; a read chain that EXACTLY equals an annotated chain is simply the
    same variant key. The ``source`` tag lives on the emitted candidate; whether annotation deserves
    any leniency is decided at final filtering, never here.
    """
    variants: List[Chain]
    read_pool: Set[str] = field(default_factory=set)
    variant_reads: Dict[Chain, Set[str]] = field(default_factory=dict)
    gtf_members: List[Tuple[str, Chain]] = field(default_factory=list)
    # Annotation-only chains (no reads of their own) that the reads do not
    # CONTRADICT. Only these may act as fold containers: absorbing a read
    # sub-chain asserts "those reads are fragments of this annotated transcript",
    # which a corrupted annotation must not be allowed to assert.
    uncontradicted_gtf: Set[Chain] = field(default_factory=set)

    @property
    def representative(self) -> Optional[Chain]:
        """Longest read variant (tie-broken by discovery support then coords), or None for a
        GTF-only family.

        A convenience pick, NOT a snapped consensus -- coordinates are never rewritten.
        """
        if not self.variants:
            return None
        return max(self.variants,
                   key=lambda c: (len(c), len(self.variant_reads.get(c, ())), c))


@dataclass
class ClusterResult:
    """Output of the clustering-only step: multi-exon families + the interval mono bucket.

    ``mono_reads`` is the holding bucket of every single-exon (chain-less) read in the
    interval -- deliberately UNPROCESSED here (not folded, not locus-split). Whether each is
    a genuine intronless transcript or a fragment of a surviving multi-exon candidate is a
    later decision, made after multi-exon resolution (see the clustering redesign).
    """
    families: List[ChainFamily]
    mono_reads: Set[str] = field(default_factory=set)


# ---------------------------------------------------------------------------
# chain-relation predicates (single interval/strand: no chrom/strand bucketing)
# ---------------------------------------------------------------------------

# Reference-consuming CIGAR footprint of one read: (match_blocks, deletion_blocks,
# skip_blocks), each a tuple of half-open genomic [start, end) blocks. Insertions
# consume no reference and never break reference continuity, so they are not kept.
ReadAlignment = Tuple[Tuple[Tuple[int, int], ...],
                      Tuple[Tuple[int, int], ...],
                      Tuple[Tuple[int, int], ...]]


def read_alignment_from_cigar(cigar, ref_start: int) -> ReadAlignment:
    """Build a :data:`ReadAlignment` from pysam-style ``cigartuples``."""
    pos = int(ref_start)
    match: List[Tuple[int, int]] = []
    dels: List[Tuple[int, int]] = []
    skips: List[Tuple[int, int]] = []
    for op, ln in cigar:
        if op in (0, 7, 8):          # M/=/X
            match.append((pos, pos + ln))
            pos += ln
        elif op == 2:                # D
            dels.append((pos, pos + ln))
            pos += ln
        elif op == 3:                # N
            skips.append((pos, pos + ln))
            pos += ln
        # I/S/H/P consume no reference
    return tuple(match), tuple(dels), tuple(skips)


def _overlap_len(blocks, lo: int, hi: int) -> int:
    return sum(max(0, min(e, hi) - max(s, lo)) for s, e in blocks)


def _cigar_exonically_covers(
    aln: ReadAlignment, donor: int, acceptor: int, *,
    flank_bp: int = 10, max_deletion_bp: int = 3, max_total_deletion_bp: int = 6,
) -> bool:
    """True iff the read's CIGAR footprint reads exonically THROUGH the intron
    ``[donor, acceptor)`` -- positive retained-intron / alternative-isoform evidence,
    as opposed to a 5'-degraded fragment that merely starts downstream of it.

    Window ``W = [donor - flank_bp, acceptor + flank_bp)``. Requires, inside W:
    the reference footprint (M/=/X/D) covers all of W (a reference_start/_end span
    is NOT enough: a big D or a soft clip leaves holes); no N; no single D longer
    than ``max_deletion_bp`` and total D at most ``max_total_deletion_bp``; each
    exonic flank has >= ``flank_bp - 2`` matched bases; the intron body has at
    least ``(acceptor - donor) - max_total_deletion_bp`` matched bases. Terminal
    clips are not a rule of their own -- they simply produce no match block.
    """
    match, dels, skips = aln
    w_lo, w_hi = donor - flank_bp, acceptor + flank_bp
    if w_hi <= w_lo:
        return False
    if any(s < w_hi and e > w_lo for s, e in skips):
        return False
    d_in = [(max(s, w_lo), min(e, w_hi)) for s, e in dels if s < w_hi and e > w_lo]
    if any(e - s > max_deletion_bp for s, e in d_in):
        return False
    if sum(e - s for s, e in d_in) > max_total_deletion_bp:
        return False
    footprint = _overlap_len(match, w_lo, w_hi) + sum(e - s for s, e in d_in)
    if footprint < w_hi - w_lo:
        return False
    need_flank = max(0, flank_bp - 2)
    if _overlap_len(match, w_lo, donor) < need_flank:
        return False
    if _overlap_len(match, acceptor, w_hi) < need_flank:
        return False
    body_need = max(0, (acceptor - donor) - max_total_deletion_bp)
    return _overlap_len(match, donor, acceptor) >= body_need


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


def _family_sort_key(fam: ChainFamily):
    """Total, deterministic ordering key over a family's chains (read variants + GTF members).

    Every family has >=1 chain, so ``min`` is well-defined even for a GTF-only family (empty
    ``variants``). The extra components make the order total under any (near-impossible) min tie.
    """
    chains = list(fam.variants) + [c for _, c in fam.gtf_members]
    return (min(chains), tuple(fam.variants), tuple(sorted(g for g, _ in fam.gtf_members)))


def cluster_families(
    read_chains: List[Tuple[Dict, IntronChain]],
    *,
    wobble_bp: int = 6,
    cassette_max_exon_bp: int = 70,
    gtf_variants: Optional[List[Tuple[str, Chain]]] = None,
    observed_junctions: Optional[Set[Tuple[int, int]]] = None,
) -> ClusterResult:
    """Clustering ONLY: group an interval's reads into multi-exon structural families.

    The clean, single-responsibility clustering step (see the clustering redesign). Unlike
    :func:`cluster_read_chains`, it does NOT fold exact sub-chains, does NOT fold or spatially
    split mono reads, and never rewrites coordinates. It ONLY:

      1. groups reads by their EXACT intron chain (3' ignored, no snap),
      2. buckets every single-exon (chain-less) read into ``mono_reads`` UNTOUCHED,
      3. union-finds the distinct multi-exon chains into families by wobble / cassette /
         containment (single-linkage; see :func:`_related`),
      4. returns each family's variant chains + a POOLED read set (reads belong to the
         family, not a member) + the per-variant discovery provenance,
      5. (optional) enters each annotated multi-exon chain as an ORDINARY zero-read variant key,
         so annotation joins a family, collapses and is emitted exactly like a read-derived chain
         (only its ``source`` tag differs; leniency, if any, is a final-filter decision). A read
         chain equal to an annotated chain IS that key. Annotation-only chains are attached AFTER
         the read-only union-find and never bridge two read-supported components.

    Deterministic: families and their variants are coordinate-sorted; ``mono_reads`` is a
    plain set. Fold, path-completion, and mono resolution are SEPARATE later per-cluster
    steps -- this function makes no structural-collapse or read-ownership decision.

    Args:
        read_chains: per-read ``(read_dict, IntronChain)`` (read_dict needs "query_name";
            CIGAR-derived chains, used as-is).
        wobble_bp: per-junction tolerance for wobble / cassette / containment joins.
        cassette_max_exon_bp: max skipped-exon length for a cassette join (0 disables).
        gtf_variants: optional annotation hypotheses as ``(gtf_id, intron_chain)`` pairs. An
            annotated chain that a read also produced is simply that key. An annotation-only chain
            becomes a zero-read variant key; whether it may merge two read-derived components is
            decided by ``observed_junctions`` (see below). They never contribute to a
            ``read_pool``. Mono (empty-chain) GTF entries are ignored here (handled by the later
            mono finalizer).
        observed_junctions: the interval's set of ``(donor, acceptor)`` junctions actually seen in
            read alignments. An annotation-only chain may join the union-find (and hence bridge)
            only when no junction of it is CONTRADICTED by the reads -- annotation may reshape the
            topology where the data is silent, but not where the data disagrees. None disables
            bridging entirely (the conservative choice for callers with no registry).

    Returns:
        A :class:`ClusterResult` (multi-exon families + the interval mono bucket).
    """
    # 1. group reads by EXACT chain; single-exon (chain-less) reads -> the mono bucket.
    #    Annotated multi-exon chains enter the SAME key space as zero-read variants
    #    (a read chain equal to an annotated chain is the same key). Mono (empty-chain)
    #    GTF entries are skipped here -- the mono finalizer matches them by overlap.
    by_chain: Dict[Chain, Set[str]] = {}
    gtf_ids_by_chain: Dict[Chain, List[str]] = {}
    mono_reads: Set[str] = set()
    for rd, chain in read_chains:
        rid = rd.get("query_name")
        if rid is None:
            continue
        introns = tuple(chain.introns)
        if not introns:
            mono_reads.add(rid)
        else:
            by_chain.setdefault(introns, set()).add(rid)
    read_chain_keys = set(by_chain)
    for gid, gchain in sorted(gtf_variants or (), key=lambda gv: (gv[1], gv[0])):
        gchain = tuple(gchain)
        if not gchain:
            continue
        by_chain.setdefault(gchain, set())
        gtf_ids_by_chain.setdefault(gchain, []).append(gid)
    # Annotated chains that are NOT also read chains: they are hypotheses, and a
    # hypothesis with zero reads must not change the connectivity of the graph the
    # DATA produced (see below).
    gtf_only_keys = set(gtf_ids_by_chain) - read_chain_keys

    # 2. union-find the distinct multi-exon chains by wobble / cassette / containment.
    #    Sorted vertex order -> deterministic; the resulting partition (connected
    #    components of the symmetric _related graph) is order-independent regardless.
    multi = sorted(by_chain)
    parent: Dict[Chain, Chain] = {c: c for c in multi}

    def find(x: Chain) -> Chain:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    # Pass 1 operates on the chains the data does not DISAGREE with: every read
    # chain, plus the annotation-only chains none of whose junctions is
    # contradicted by a read junction at the same splice site. Those are ordinary
    # vertices -- they may join and bridge like any read chain.
    def _uncontradicted(j: Tuple[int, int]) -> bool:
        """The reads do not DISAGREE with this annotated junction.

        Three cases, and only the middle one blocks:
          * the junction is observed in a read           -> the data confirms it;
          * the reads put a DIFFERENT junction within the
            alignment tolerance of one of its splice sites -> the data contradicts it;
          * no read junction is anywhere near it         -> the data is silent,
            which is not disagreement (coverage is finite, and a real junction of
            a lowly expressed isoform is routinely unobserved).
        This is what separates a truncated annotation (every junction real, so
        nothing contradicts it) from a merged or jittered one (its fabricated
        junction sits right beside the true junction the reads did produce).
        """
        if j in observed_junctions:
            return True
        d, a = j
        for od, oa in observed_junctions:
            if (od, oa) != j and (abs(od - d) <= wobble_bp
                                  or abs(oa - a) <= wobble_bp):
                # The reads put a DIFFERENT junction within the alignment
                # tolerance of one of this junction's splice sites. The window is
                # the project's existing junction-equivalence tolerance (a
                # property of alignment precision), NOT a constant chosen against
                # a benchmark's corruption magnitude.
                return False
        return True

    supported_gtf = {
        g for g in gtf_only_keys
        if observed_junctions is not None and all(_uncontradicted(j) for j in g)
    }
    unsupported_gtf = gtf_only_keys - supported_gtf
    read_multi = [c for c in multi if c not in unsupported_gtf]
    for i in range(len(read_multi)):
        for j in range(i + 1, len(read_multi)):
            a, b = read_multi[i], read_multi[j]
            if find(a) != find(b) and _related(a, b, wobble_bp, cassette_max_exon_bp):
                parent[find(a)] = find(b)

    # Pass 2: attach each annotation-only chain to AT MOST ONE existing component.
    # It is an ordinary variant of the family it joins (a legal fold container, and
    # emitted through the same path as a read chain), but it can never BRIDGE two
    # read-supported components: single-linkage would otherwise let one unsupported
    # -- possibly corrupted -- annotation merge two families and silently rewrite
    # the collapse containers, the EM competition scope and the isoform-fraction
    # denominator (Codex gpt-5.6-sol 2026-09-04: "annotation as a prior is fine;
    # annotation changing the topology the data produced is not").
    # An annotation-only chain has no reads of its own, so whether it may change
    # the read-derived partition is decided by the DATA, junction by junction:
    #
    #   every junction of the annotated chain is observed in some read  ->  it may
    #   join ALL related components, bridging them if need be;
    #   otherwise                                                       ->  it may
    #   join at most ONE (and none when several are equally related).
    #
    # A chain whose junctions are NOT all observed gets the conservative rule: it
    # may join at most one related component, and none when several are equally
    # related (nothing shows which it belongs to, and an arbitrary tie-break must
    # not decide a fold container).
    #
    # Rationale, and it is measured rather than assumed. Letting annotation bridge
    # unconditionally is a large recall win where the annotation is trustworthy
    # (Sequin expression contract: c_trunc3p +4.8 F1, full/p90/p99/skip/spurious
    # +3.4) because the annotated chain pools reads that the data alone splits.
    # But it is exactly how a CORRUPTED annotation does damage: a fabricated
    # merge junction welds two real families together (c_merge50 -2.4 F1, Pr
    # -7.9; c_jitter50 -1.3). The junction-support test separates the two cases at
    # the source: a truncated annotation's junctions are all real, so it bridges;
    # a merged or jittered annotation contains a junction no read ever produced,
    # so it cannot. `observed_junctions` is the interval's raw read junction
    # registry, which discovery collects on every run.
    for g in sorted(unsupported_gtf):
        roots = {find(c) for c in read_multi
                 if _related(g, c, wobble_bp, cassette_max_exon_bp)}
        if len(roots) == 1:
            parent[find(g)] = roots.pop()
    # Annotation-only chains never join each other directly either: two unsupported
    # annotations forming their own family would recreate the bridge one level down.

    comps: Dict[Chain, List[Chain]] = {}
    for c in multi:
        comps.setdefault(find(c), []).append(c)

    # 3. build families: variants coord-sorted, read pool = union of members' reads,
    #    per-variant discovery reads kept as provenance (not final ownership).
    families: List[ChainFamily] = []
    for comp in comps.values():
        variants = sorted(comp)
        pool: Set[str] = set()
        for c in comp:
            pool |= by_chain[c]
        families.append(ChainFamily(
            variants=variants,
            read_pool=pool,
            variant_reads={c: set(by_chain[c]) for c in comp},
            gtf_members=[(gid, c) for c in variants for gid in gtf_ids_by_chain.get(c, ())],
            uncontradicted_gtf={c for c in variants if c in supported_gtf},
        ))

    families.sort(key=_family_sort_key)

    return ClusterResult(families=families, mono_reads=mono_reads)


def collapse(family: ChainFamily, *, span_guard: bool = False,
             read_alignments: Optional[Dict[str, ReadAlignment]] = None) -> ChainCluster:
    """Collapse a single (post-explore) family into a :class:`ChainCluster`, folding EXACT
    contiguous sub-chains into their longest exact container.

    The explicit, single-family FOLD step of the clustering redesign (replaces the fold that
    the old :func:`cluster_read_chains` did inline). Design (locked with the user + Codex
    gpt-5.6-sol xhigh):

    - **fold-all** (philosophy A): every variant that is an EXACT contiguous sub-chain of a
      longer variant folds into its LONGEST exact container. dRNA 5' degradation ladders (same 3'
      end, successively shorter 5') are exactly these exact sub-chains; whether a folded sub-chain
      is a genuine short isoform is decided LATER by the post-EM 5'-TSS read-end recovery -- never
      guessed here from read counts.
    - **span_guard** (optional, the CIGAR-aware successor of ``cluster_read_chains``'s
      ``fold_span_guard``): when set, a sub-chain read folds into its container ONLY if its
      CIGAR footprint does NOT read exonically THROUGH one of the container's EXTRA introns
      (introns the container splices but the sub-chain lacks; see
      :func:`_cigar_exonically_covers`). A read that genuinely reads through such an intron is
      positive evidence for a retained-intron / alternative isoform, NOT 5' degradation, so it
      is KEPT as its own candidate (the sub-chain survives with only its contradicting reads).
      5' degradation reads start DOWNSTREAM of the missing introns and fold; so do reads whose
      reference span merely brackets the intron via a large deletion, an N, or a soft clip.
      Needs ``read_alignments`` (query_name -> :data:`ReadAlignment`); without it, span_guard
      is a no-op (folds unconditionally).
    - **exact only**: wobble / cassette / non-exact-containment siblings are NOT folded (they
      are not exact sub-chains); they stay as separate members and compete in EM, where the
      abundance-based selection decides (Codex's pre-EM=canonicalisation / post-EM=biological
      boundary).
    - **absorb-only**: a fold target is always an existing longer variant; no chain is
      synthesised (synthesis is :func:`explore`'s job).
    - reads are pooled: a survivor's ``read_ids`` already includes every folded shadow's reads
      (shadows add no mass); folded shadows are kept as dormant provenance (``GenCandidate.folded``)
      that feeds the downstream 5'-TSS recovery.

    The PER-MEMBER fold/shadow structure matches the old inline fold, so the discovery emission
    loop, ``shadow_map`` and TSS recovery consume it unchanged. EM SCOPING may differ, though:
    :func:`cluster_families` groups by single-linkage over ``_related`` (wobble/cassette/containment
    within ``wobble_bp``), so an exact sub-chain can BRIDGE two otherwise-unrelated maximal chains
    into one family (e.g. ``(I1,I2)``–``(I2,)``–``(I2,I3)``). The old path folded that sub-chain away
    BEFORE clustering, so it produced two clusters; here the two maximal chains stay co-scoped in one
    cluster. This is an INTENDED consequence of doing clustering once (in ``cluster_families``) and
    fold second -- a behavior change validated by metrics, not a byte-identical reproduction. collapse
    does NOT re-cluster the survivors.

    Precondition (guaranteed by ``cluster_families`` output): ``variant_reads`` PARTITIONS
    ``read_pool`` (each read is in exactly one variant's set), so no read is lost or double-counted
    and ``cluster.read_ids == family.read_pool``. A variant missing from ``variant_reads`` contributes
    no reads (defensive; e.g. a future explore-assembled variant whose read ownership the adapter must
    register before calling collapse). An annotated chain that READS also produced is an
    ordinary key and may contain; an annotation-ONLY key (zero reads) may neither fold nor be
    folded -- it must not delete a read-supported sub-chain before any evidence is weighed.
    Annotated chains are never foldees themselves: the fold is a read-level degradation judgement that a
    genuine 5'-truncated annotated isoform can never win, so annotated keys always stay members
    and EM + the final GTF gates decide their existence. Deterministic (coordinate-sorted); never
    mutates ``family``.
    """
    vr = family.variant_reads
    alns = read_alignments or {}
    annotated = {c for _, c in family.gtf_members}
    # Annotation-only = annotated AND carrying no reads of its own. Such a key may
    # not ABSORB a read-derived sub-chain: a hypothesis with zero read support must
    # not delete a read-supported one before any evidence is weighed. (An annotated
    # chain that reads also produced is an ordinary read key and may contain.)
    # An annotation-only chain (annotated, no reads of its own) may absorb a read
    # sub-chain ONLY when the reads do not contradict it. Measured on the Sequin
    # expression contract: allowing it unconditionally is +3.4 Sn / -1.1 Pr
    # overall -- a large recall win wherever the annotation is trustworthy
    # (c_trunc3p +4.6 F1, full/p90/p99/skip/spurious +2.9) because the reads of a
    # 5'-degraded fragment get pooled onto the right model, but a clear loss
    # exactly where the annotation is corrupted (c_merge50 -1.7 F1 with Pr -6.6,
    # c_jitter50 -1.3 with Pr -3.8). The contradiction test separates the two at
    # the source, so neither case needs a threshold.
    annotation_only = {c for c in annotated
                       if not vr.get(c) and c not in family.uncontradicted_gtf}
    # longest-first: `next(exact-subchain of a kept chain)` then picks the LONGEST container.
    chains_sorted = sorted(family.variants, key=lambda c: (-len(c), -len(vr.get(c, ())), c))
    kept: List[Chain] = []
    reads: Dict[Chain, Set[str]] = {}
    folded_of: Dict[Chain, List[Tuple[Chain, Set[str]]]] = {}  # container -> shadow (chain, reads)
    for c in chains_sorted:
        container = next((k for k in kept
                          if k not in annotation_only and _exact_subchain(c, k)), None)
        if container is None or c in annotated:
            # An annotated key is never a FOLDEE (it may be a container): the fold
            # is a read-level "these reads are 5' degradation of the container"
            # decision, and a genuine 5'-truncated annotated isoform's reads can
            # never refute that (they always start downstream of the extra
            # introns), while a zero-read annotated key has nothing to judge on.
            # Its existence is decided by EM + the final GTF gates, not here.
            # Measured (156-cell SIRV4 grid, 2026-09-03): folding annotated exact
            # sub-chains cost -3.0 F1 (Sn -3.8) in every guided condition.
            kept.append(c)
            reads[c] = set(vr.get(c, ()))
            continue
        cr = vr.get(c, set())
        if not span_guard or not cr:
            # no reads and no read evidence -> nothing can contradict the fold
            reads[container] |= cr
            folded_of.setdefault(container, []).append((c, set(cr)))
            continue
        # span guard: split c's reads into those that fold (5' degradation) and those whose
        # CIGAR reads exonically THROUGH one of the container's EXTRA introns.
        off = _subchain_offset(c, container)
        extra = container[:off] + container[off + len(c):]
        fold_r: Set[str] = set()
        keep_r: Set[str] = set()
        for rid in cr:
            aln = alns.get(rid)
            if aln is not None and any(_cigar_exonically_covers(aln, d, a) for d, a in extra):
                keep_r.add(rid)          # read through an extra intron -> distinct isoform
            else:
                fold_r.add(rid)
        if fold_r:
            reads[container] |= fold_r
            folded_of.setdefault(container, []).append((c, set(fold_r)))
        if keep_r:
            kept.append(c)               # sub-chain survives with only its contradicting reads
            reads[c] = keep_r

    def _mk(c: Chain) -> GenCandidate:
        return GenCandidate(
            IntronChain(introns=c), set(reads[c]),
            folded=[GenCandidate(IntronChain(introns=fc), set(fr))
                    for fc, fr in folded_of.get(c, [])])

    members = [_mk(c) for c in sorted(kept)]          # coord-sorted -> deterministic
    union: Set[str] = set()
    for c in kept:
        union |= reads[c]
    return ChainCluster(members=members, read_ids=union)
