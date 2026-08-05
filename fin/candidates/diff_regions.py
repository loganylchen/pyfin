"""Cluster-level DIFF REGIONS: the loci that distinguish a cluster's candidates, computed
structurally from the candidates ALONE (no reads, no signal).

A cluster holds N multi-exon candidate transcripts that mostly agree but differ at a few places
(a wobbling junction, a short cassette exon, a truncation). This module locates those places once,
up front, as ``DiffRegion``s. Each diff region carries, per candidate, that candidate's local cDNA
sequence through the region -- its **allele** -- so downstream layers (mappy sequence match, then
krill signal) can ask "which allele does this read support?".

Two invariants, both learned the hard way:

1. **Always include shared context.** A diff region is NOT just the differing bases: it is the
   differing bases PLUS ``context`` bp of the flanking sequence that every candidate shares. Without
   context, a candidate that SKIPS the differing part (its transcript splices across it) would have
   "nothing covering" the region -- there would be no sequence to compare. With context, that skipper
   still carries the shared flank bases (its splice-junction sequence), so every candidate that reaches
   the locus has a concrete, comparable allele. ``context`` is one k-mer's worth of flank (3-5 bp).

2. **Exclude pure endpoints.** A locus is a diff region only where candidates disagree exon-vs-intron
   (``exon_any & intron_any``). A candidate merely being longer/shorter at a transcript END (no
   opposing splice) is a poly(A)/TSS boundary call, not an internal splice decision, and is left out --
   otherwise every ragged 5'/3' end would masquerade as a distinguishing feature.

Coincidental sequence matches at a junction (e.g. a skipped exon's first bases happening to equal the
spliced-in exon's first bases) are NOT resolved here; they are a property of the allele *sequences*
and are handled by whoever compares reads to alleles. This module's job is purely structural: where
the candidates differ, and what each one looks like there (with shared context).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from fin.candidates.dataclasses import TranscriptCandidate

# A diff region at or below this genomic span is a wobbling junction; a short exon is a cassette.
_WOBBLE_MAX_BP = 12
_SHORT_EXON_MAX_BP = 200


@dataclass
class DiffRegion:
    """One locus where a cluster's candidates structurally disagree.

    Attributes:
        g_lo, g_hi: genomic half-open span of the differing bases (context NOT included; the context
            is baked into each allele sequence instead).
        kind: ``"wobble"`` | ``"cassette"`` | ``"terminal"`` | ``"complex"`` -- informational.
        context: max bp of shared flank included on each side of every allele.
        alleles: candidate_id -> that candidate's cDNA over the diff plus the flank bases that ALL
            candidates here place at the same genomic position (context common to all, so allele fits
            are comparable), transcript 5'->3'. Only comparable candidates that form a non-empty allele
            appear; gross structural outliers and unreached candidates are omitted.
        groups: candidate_ids grouped by identical allele sequence (the partition a read could resolve).
    """

    g_lo: int
    g_hi: int
    kind: str
    context: int
    alleles: Dict[str, Optional[str]]
    groups: List[List[str]] = field(default_factory=list)


def _exons(cand: TranscriptCandidate) -> List[Tuple[int, int]]:
    """Genomic-order exon intervals (half-open) of a candidate."""
    introns = cand.intron_chain.introns
    if not introns:
        return [(cand.start, cand.end)]
    exons = [(cand.start, introns[0][0])]
    for k in range(len(introns) - 1):
        exons.append((introns[k][1], introns[k + 1][0]))
    exons.append((introns[-1][1], cand.end))
    return exons


def _genome_to_cdna_point(exons: Sequence[Tuple[int, int]], total_len: int, g: int) -> int:
    """cDNA offset (in GENOMIC-order frame, i.e. before any minus-strand flip) of genomic pos ``g``.

    A position inside an intron (or in the gap before an exon) snaps to the NEXT exon boundary -- the
    splice junction -- so a candidate that splices across ``g`` maps it to a single junction offset.
    Positions at/after the last exon map to ``total_len``.
    """
    off = 0
    for ex_s, ex_e in exons:
        if g < ex_s:
            return off                      # in the intron/gap before this exon -> junction
        if g < ex_e:
            return off + (g - ex_s)         # inside this exon
        off += ex_e - ex_s
    return off                              # at/after the last exon end


def _cdna_span(cand: TranscriptCandidate, g_lo: int, g_hi: int) -> Optional[Tuple[int, int]]:
    """(c_lo, c_hi) cDNA span (sequence 5'->3' frame) of the diff window ``[g_lo, g_hi)``, or None if
    the candidate does not reach it. A candidate that splices ACROSS the window (skips it) gets
    ``c_lo == c_hi`` -- the junction point."""
    if g_hi <= cand.start or g_lo >= cand.end:
        return None
    total = len(cand.sequence)
    if total == 0:
        return None
    exons = _exons(cand)
    c_lo = _genome_to_cdna_point(exons, total, g_lo)
    c_hi = _genome_to_cdna_point(exons, total, g_hi)
    if cand.strand == "-":                  # cDNA is reverse-complemented on the minus strand
        c_lo, c_hi = total - c_hi, total - c_lo
    if c_hi < c_lo:
        c_lo, c_hi = c_hi, c_lo
    return (c_lo, c_hi)


def _cdna_to_genome(exons: Sequence[Tuple[int, int]], strand: str, total: int, p: int) -> Optional[int]:
    """Genomic coordinate of cDNA (sequence 5'->3' frame) position ``p`` -- the inverse of the cDNA
    mapping, used to test whether two candidates' flank bases land on the SAME genomic position."""
    if p < 0 or p >= total:
        return None
    order = list(reversed(exons)) if strand == "-" else exons
    t = 0
    for s, e in order:
        L = e - s
        if t <= p < t + L:
            off = p - t
            return (e - 1 - off) if strand == "-" else (s + off)
        t += L
    return None


def _alleles_with_shared_context(grp, s_lo: int, s_hi: int, k: int) -> Dict[str, Optional[str]]:
    """Per-candidate allele = the candidate's cDNA over the diff ``[s_lo, s_hi)`` plus the SHARED
    flanking context: only those flank bases (up to ``k`` each side) that EVERY candidate in ``grp``
    places at the same genomic position. Walking out from the diff in cDNA, the moment the candidates'
    flank bases diverge genomically the context stops -- so the context is common to all (comparable),
    while the differing middle is what a read tells apart. ``None`` for a candidate that cannot form a
    non-empty allele (unreached, or a zero-width diff with no shared flank on either side)."""
    # info aligned with grp by index (TranscriptCandidate is unhashable -> no object keys).
    info: List[Optional[Tuple[Tuple[int, int], List[Tuple[int, int]], int]]] = []
    for c in grp:
        span = _cdna_span(c, s_lo, s_hi)
        info.append(None if span is None else (span, _exons(c), len(c.sequence)))
    present = [i for i in range(len(grp)) if info[i] is not None]
    if len(present) < 2:
        return {c.candidate_id: None for c in grp}

    def shared(side: str) -> int:
        m = 0
        for d in range(1, k + 1):
            gs = set()
            for i in present:
                (c_lo, c_hi), exons, total = info[i]
                p = (c_lo - d) if side == "left" else (c_hi + d - 1)
                g = _cdna_to_genome(exons, grp[i].strand, total, p)
                if g is None:
                    return m
                gs.add(g)
            if len(gs) == 1:                 # every candidate's d-th flank base is the same genomic pos
                m = d
            else:
                return m
        return m

    left, right = shared("left"), shared("right")
    out: Dict[str, Optional[str]] = {}
    for i, c in enumerate(grp):
        if info[i] is None:
            out[c.candidate_id] = None
            continue
        (c_lo, c_hi), _, total = info[i]
        lo = max(0, c_lo - left)
        hi = min(total, c_hi + right)
        out[c.candidate_id] = c.sequence[lo:hi] if hi > lo else None
    return out


def _classify(span: int, terminal: bool) -> str:
    if terminal:
        return "terminal"
    if span <= _WOBBLE_MAX_BP:
        return "wobble"
    if span <= _SHORT_EXON_MAX_BP:
        return "cassette"
    return "complex"


def _group(alleles: Dict[str, Optional[str]]) -> List[List[str]]:
    by_seq: Dict[str, List[str]] = {}
    for cid, seq in alleles.items():
        if seq is None:
            continue
        by_seq.setdefault(seq, []).append(cid)
    # deterministic: order groups by first (sorted) member
    return [sorted(m) for _, m in sorted(by_seq.items(), key=lambda kv: sorted(kv[1])[0])]


def _segment(cands: Sequence[TranscriptCandidate], w_lo: int, w_hi: int) -> List[Tuple[int, int]]:
    """Maximal genomic runs within ``[w_lo, w_hi)`` where ``cands`` disagree exon-vs-intron.

    Pure transcript-END differences (a candidate simply absent, no opposing splice) never set both an
    exon and an intron label at a base, so they do not appear -- endpoints are excluded by construction.
    """
    span = w_hi - w_lo
    if span <= 0 or len(cands) < 2:
        return []
    labels = np.zeros((len(cands), span), dtype=np.int8)
    for i, c in enumerate(cands):
        lo = max(c.start, w_lo) - w_lo
        hi = min(c.end, w_hi) - w_lo
        if hi > lo:
            labels[i, lo:hi] = 1
        for i_s, i_e in c.intron_chain.introns:
            il = max(i_s, w_lo) - w_lo
            ih = min(i_e, w_hi) - w_lo
            if ih > il:
                labels[i, il:ih] = 2
    diff = (labels == 1).any(axis=0) & (labels == 2).any(axis=0)
    loci: List[Tuple[int, int]] = []
    i = 0
    while i < span:
        if not diff[i]:
            i += 1
            continue
        j = i
        while j < span and diff[j]:
            j += 1
        loci.append((i + w_lo, j + w_lo))
        i = j
    return loci


def _flank_genomic(cand: TranscriptCandidate, g_lo: int, g_hi: int, k: int):
    """Genomic positions of the ``k`` EXONIC bases flanking ``[g_lo, g_hi)`` on each side, as
    ``(left, right)`` tuples -- or ``None`` if the candidate does not reach the locus.

    Two candidates are COMPARABLE at a locus iff both flank tuples are equal: they then share the same
    anchoring context on both sides, so a read's mapping/signal fit to one vs the other is measured
    against a common frame. A candidate whose splice reconverges at the same flanking exons (a wobble
    or a single small skip) is comparable; one whose intron runs far past the locus (a large multi-exon
    skip) lands on different flanking exons and is therefore NOT compared here -- such gross structural
    differences are trivially separable and handled outside this fine machinery.
    """
    if g_hi <= cand.start or g_lo >= cand.end:
        return None
    exons = _exons(cand)
    left: List[int] = []
    for s, e in exons:                              # highest k exonic bases strictly below g_lo
        hi = min(e, g_lo)
        if hi > s:
            left.extend(range(max(s, hi - k), hi))
    right: List[int] = []
    for s, e in exons:                              # lowest k exonic bases at/above g_hi
        lo = max(s, g_hi)
        if e > lo:
            right.extend(range(lo, min(e, lo + k)))
            if len(right) >= k:
                break
    return (tuple(left[-k:]), tuple(right[:k]))


def cluster_diff_regions(
    candidates: Sequence[TranscriptCandidate],
    *,
    context: int = 3,
) -> List[DiffRegion]:
    """Locate the FINE diff regions of a cluster of multi-exon candidates.

    Only comparable differences are emitted: a wobbling junction or a single short exon skip, among the
    subset of candidates that share the SAME flanking context on both sides (so a read's fit to one vs
    another is measured against a common anchor). Candidates separated by a gross structural difference
    (a large multi-exon skip, a fused/retained region, a truncation) land on different flanking exons
    and are dropped from the fine region -- they are trivially distinguishable by whole-block presence,
    not by these narrow context-anchored regions.

    Algorithm: (1) segment among all candidates into raw disagreement loci; (2) within each raw locus,
    group candidates by shared flank context; (3) re-segment TIGHTLY within each comparable group (a
    shared internal exon now breaks a lumped region apart) and emit a region per sub-locus that still
    carries >=2 distinct alleles.

    Args:
        candidates: the cluster's candidates. Mono-exon candidates (empty intron chain) contribute no
            structure; at least 2 multi-exon candidates are needed.
        context: bp of shared flanking sequence baked into each allele on each side (one k-mer worth).

    Returns:
        Diff regions sorted by genomic position. Empty if the candidates never comparably disagree.
    """
    cands = [c for c in candidates if c.intron_chain.introns]
    if len(cands) < 2:
        return []
    g_lo = min(c.start for c in cands)
    g_hi = max(c.end for c in cands)
    if g_hi <= g_lo:
        return []

    regions: List[DiffRegion] = []
    seen: set = set()
    for r_lo, r_hi in _segment(cands, g_lo, g_hi):
        groups: Dict[object, List[TranscriptCandidate]] = {}
        for c in cands:
            groups.setdefault(_flank_genomic(c, r_lo, r_hi, context), []).append(c)
        for fp, grp in groups.items():
            if fp is None or len(grp) < 2:            # unreachable, or a lone structural outlier
                continue
            for s_lo, s_hi in _segment(grp, r_lo, r_hi):
                raw = _alleles_with_shared_context(grp, s_lo, s_hi, context)
                alleles = {cid: s for cid, s in raw.items() if s is not None}
                if len(set(alleles.values())) < 2:
                    continue
                key = (s_lo, s_hi, tuple(sorted(alleles)))
                if key in seen:                       # a sub-locus shared by overlapping groups
                    continue
                seen.add(key)
                terminal = any(s_lo < c.start < s_hi or s_lo < c.end < s_hi
                               for c in grp if c.candidate_id in alleles)
                regions.append(DiffRegion(
                    g_lo=s_lo, g_hi=s_hi, kind=_classify(s_hi - s_lo, terminal),
                    context=context, alleles=alleles, groups=_group(alleles)))
    # tie-break on the allele-key so two flank groups emitting the same span order deterministically
    # regardless of input candidate order.
    regions.sort(key=lambda r: (r.g_lo, r.g_hi, tuple(sorted(r.alleles))))
    return regions
