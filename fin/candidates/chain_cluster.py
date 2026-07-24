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

    ``gtf_members`` are annotation transcripts attached to this family as source-tagged, ZERO-read
    structural hypotheses: each is ``(gtf_id, intron_chain)``. They carry NO reads (never in
    ``read_pool`` / ``variant_reads``) and get read mass ONLY if their exact model wins the downstream
    evidence assignment -- annotation is another competitor, never a snap target or a read-pool source.
    A family with empty ``variants`` but non-empty ``gtf_members`` is a GTF-only family (an annotated
    transcript that matched no read-derived structure; kept so it survives to the abundance filters).
    """
    variants: List[Chain]
    read_pool: Set[str] = field(default_factory=set)
    variant_reads: Dict[Chain, Set[str]] = field(default_factory=dict)
    gtf_members: List[Tuple[str, Chain]] = field(default_factory=list)

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


def _attach_gtf_variants(
    families: List[ChainFamily],
    gtf_variants: List[Tuple[str, Chain]],
    wobble_bp: int,
    cassette_max_exon_bp: int,
) -> None:
    """Attach each GTF hypothesis to the single most-related read family (in place).

    NON-BRIDGING: a GTF is added to at most ONE family and never merges two read families. A GTF
    related to no read family becomes its own GTF-only family (empty read pool). Mono (empty-chain)
    GTF entries are skipped (handled by the later mono finalizer). Deterministic: GTF entries are
    processed in sorted order and the best family is chosen by (most related variants, then smallest
    first-variant coords). GTF is a zero-read hypothesis: it is NEVER added to any ``read_pool`` or
    ``variant_reads``.
    """
    read_families = [f for f in families if f.variants]   # snapshot: a GTF-only family never attracts GTF
    gtf_only: List[ChainFamily] = []
    for gid, gchain in sorted(gtf_variants, key=lambda gv: (gv[1], gv[0])):
        if not gchain:                     # mono / empty-chain GTF -> deferred to the mono finalizer
            continue
        best = None                        # (n_related, first_variant, family)
        for fam in read_families:
            n = sum(1 for v in fam.variants
                    if _related(gchain, v, wobble_bp, cassette_max_exon_bp))
            if n == 0:
                continue
            if best is None or n > best[0] or (n == best[0] and fam.variants[0] < best[1]):
                best = (n, fam.variants[0], fam)
        if best is not None:
            best[2].gtf_members.append((gid, gchain))
        else:
            gtf_only.append(ChainFamily(
                variants=[], read_pool=set(), variant_reads={},
                gtf_members=[(gid, gchain)]))
    families.extend(gtf_only)


def cluster_families(
    read_chains: List[Tuple[Dict, IntronChain]],
    *,
    wobble_bp: int = 6,
    cassette_max_exon_bp: int = 70,
    gtf_variants: Optional[List[Tuple[str, Chain]]] = None,
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
      5. (optional) attaches each GTF hypothesis to the single best-related read family as a
         ZERO-read ``gtf_member`` -- NON-BRIDGING (a GTF never merges two read families) and
         never a read-pool source; a GTF matching no read family becomes a GTF-only family.

    Deterministic: families and their variants are coordinate-sorted; ``mono_reads`` is a
    plain set. Fold, path-completion, and mono resolution are SEPARATE later per-cluster
    steps -- this function makes no structural-collapse or read-ownership decision.

    Args:
        read_chains: per-read ``(read_dict, IntronChain)`` (read_dict needs "query_name";
            CIGAR-derived chains, used as-is).
        wobble_bp: per-junction tolerance for wobble / cassette / containment joins.
        cassette_max_exon_bp: max skipped-exon length for a cassette join (0 disables).
        gtf_variants: optional annotation hypotheses as ``(gtf_id, intron_chain)`` pairs. Each is
            attached to the single most-related read family (by the same ``_related`` predicate) as a
            zero-read ``gtf_member``; it NEVER bridges two read families and NEVER contributes to a
            ``read_pool``. Mono (empty-chain) GTF entries are ignored here (annotation of single-exon
            transcripts is handled by the later mono finalizer). The read families are formed BEFORE
            attachment, so annotation cannot change read-derived family structure.

    Returns:
        A :class:`ClusterResult` (multi-exon families + the interval mono bucket).
    """
    # 1. group reads by EXACT chain; single-exon (chain-less) reads -> the mono bucket.
    by_chain: Dict[Chain, Set[str]] = {}
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

    for i in range(len(multi)):
        for j in range(i + 1, len(multi)):
            a, b = multi[i], multi[j]
            if find(a) != find(b) and _related(a, b, wobble_bp, cassette_max_exon_bp):
                parent[find(a)] = find(b)

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
        ))

    # 4. (optional) attach GTF hypotheses to read families -- NON-BRIDGING, zero-read.
    if gtf_variants:
        _attach_gtf_variants(families, gtf_variants, wobble_bp, cassette_max_exon_bp)

    # Stable, total family order over every family's chains (read variants + gtf members),
    # so GTF-only families (empty ``variants``) also sort deterministically.
    families.sort(key=_family_sort_key)

    return ClusterResult(families=families, mono_reads=mono_reads)


def collapse(family: ChainFamily) -> ChainCluster:
    """Collapse a single (post-explore) family into a :class:`ChainCluster`, folding EXACT
    contiguous sub-chains into their longest exact container.

    The explicit, single-family FOLD step of the clustering redesign (replaces the fold that
    the old :func:`cluster_read_chains` did inline). Design (locked with the user + Codex
    gpt-5.6-sol xhigh):

    - **fold-all** (philosophy A): every variant that is an EXACT contiguous sub-chain of a
      longer variant folds into its LONGEST exact container, UNCONDITIONALLY (no read-support
      guard). dRNA 5' degradation ladders (same 3' end, successively shorter 5') are exactly
      these exact sub-chains; whether a folded sub-chain is a genuine short isoform is decided
      LATER by the post-EM 5'-TSS read-end recovery -- never guessed here.
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
    register before calling collapse). A GTF-only family (empty ``variants``) yields an empty cluster.
    Deterministic (coordinate-sorted); never mutates ``family``.
    """
    vr = family.variant_reads
    # longest-first: `next(exact-subchain of a kept chain)` then picks the LONGEST container.
    chains_sorted = sorted(family.variants, key=lambda c: (-len(c), -len(vr.get(c, ())), c))
    kept: List[Chain] = []
    reads: Dict[Chain, Set[str]] = {}
    folded_of: Dict[Chain, List[Tuple[Chain, Set[str]]]] = {}  # container -> shadow (chain, reads)
    for c in chains_sorted:
        container = next((k for k in kept if _exact_subchain(c, k)), None)
        if container is None:
            kept.append(c)
            reads[c] = set(vr.get(c, ()))
        else:
            cr = vr.get(c, set())
            reads[container] |= cr
            folded_of.setdefault(container, []).append((c, set(cr)))

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
