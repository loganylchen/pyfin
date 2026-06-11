"""Stage F2 of fusion detection: per-arm splice-structure inference.

Given the chimeric reads resolved by Stage F1, this module behaves like a
mini-assembly run on each fusion arm:

  * Chimeric reads are clustered by partner pair (the two loci they connect)
    and breakpoint proximity.
  * For each arm of a cluster, candidate splice structures are inferred from the
    reads themselves (CIGAR intron chains + canonical "wobble" alternatives) and,
    in PARALLEL, from any overlapping annotation transcripts. Annotation chains
    are ADDED as extra candidate options — they do NOT snap or overwrite the
    read-derived structure (mirrors how assembly keeps GTF and novel candidates
    side by side).

Stage F3 takes the per-arm variant pools and stitches them across the breakpoint
into composite fusion-transcript candidates.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple

from fin.candidates.dataclasses import IntronChain
from fin.candidates.intron_chains import (
    extract_intron_chain,
    gtf_transcript_to_intron_chain,
)
from fin.fusion.chimeric import ArmAlignment, ChimericRead

logger = logging.getLogger(__name__)


@dataclass
class ArmVariant:
    """One candidate splice structure for a single fusion arm.

    Coordinates are GENOMIC (0-based, half-open). ``start``/``end`` are the arm's
    extent (used by Stage F3 to synthesize the spliced arm sequence). ``source``
    is "read" (inferred from read CIGAR, possibly wobble-shifted) or "annotation"
    (an overlapping GTF transcript offered as an additional option).
    """

    chrom: str
    strand: str
    start: int
    end: int
    intron_chain: IntronChain
    source: str  # "read" | "annotation"
    supporting_read_ids: Set[str] = field(default_factory=set)


@dataclass
class FusionPairCluster:
    """Chimeric reads connecting the same partner pair, with per-arm variants.

    ``signature`` = (chromA, strandA, chromB, strandB). ``breakpoint_a`` /
    ``breakpoint_b`` are the consensus (chrom, pos, strand) on each side.
    """

    signature: Tuple[str, str, str, str]
    breakpoint_a: Tuple[str, int, str]
    breakpoint_b: Tuple[str, int, str]
    reads: List[ChimericRead]
    arm_a_variants: List[ArmVariant] = field(default_factory=list)
    arm_b_variants: List[ArmVariant] = field(default_factory=list)


def cluster_chimeric_reads(
    reads: List[ChimericRead],
    max_dist: int = 500,
) -> List[FusionPairCluster]:
    """Cluster chimeric reads by partner-pair signature and breakpoint proximity.

    Two reads join the same cluster when they share the (chromA, strandA, chromB,
    strandB) signature AND both breakpoint positions are within ``max_dist`` bp.
    Single-linkage via union-find, with deterministic (signature-sorted) output.
    """
    if not reads:
        return []

    # Bucket by signature first.
    by_sig: Dict[Tuple[str, str, str, str], List[int]] = {}
    for i, r in enumerate(reads):
        sig = (r.breakpoint_a[0], r.breakpoint_a[2], r.breakpoint_b[0], r.breakpoint_b[2])
        by_sig.setdefault(sig, []).append(i)

    parent = list(range(len(reads)))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x: int, y: int) -> None:
        rx, ry = find(x), find(y)
        if rx == ry:
            return
        parent[max(rx, ry)] = min(rx, ry)

    for indices in by_sig.values():
        indices.sort(key=lambda i: (reads[i].breakpoint_a[1], reads[i].breakpoint_b[1], i))
        for a in range(len(indices)):
            ia = indices[a]
            for b in range(a + 1, len(indices)):
                ib = indices[b]
                if reads[ib].breakpoint_a[1] - reads[ia].breakpoint_a[1] > max_dist:
                    break  # sorted on breakpoint_a: no further read can be closer
                if abs(reads[ib].breakpoint_b[1] - reads[ia].breakpoint_b[1]) <= max_dist:
                    union(ia, ib)

    groups: Dict[int, List[ChimericRead]] = {}
    for i, r in enumerate(reads):
        groups.setdefault(find(i), []).append(r)

    clusters: List[FusionPairCluster] = []
    for members in groups.values():
        rep = members[0]
        sig = (rep.breakpoint_a[0], rep.breakpoint_a[2],
               rep.breakpoint_b[0], rep.breakpoint_b[2])
        pos_a = round(sum(m.breakpoint_a[1] for m in members) / len(members))
        pos_b = round(sum(m.breakpoint_b[1] for m in members) / len(members))
        clusters.append(FusionPairCluster(
            signature=sig,
            breakpoint_a=(rep.breakpoint_a[0], pos_a, rep.breakpoint_a[2]),
            breakpoint_b=(rep.breakpoint_b[0], pos_b, rep.breakpoint_b[2]),
            reads=members,
        ))

    clusters.sort(key=lambda c: (c.signature, c.breakpoint_a[1], c.breakpoint_b[1]))
    return clusters


def _read_derived_variants(
    arms: List[Tuple[str, ArmAlignment]],
    genome_seq: str,
    strand: str,
    motif_set,
    search_bp: int,
    max_chains_per_read: int,
    min_support: int,
) -> List[ArmVariant]:
    """Read-inferred chain variants for one arm (CIGAR chains + wobble)."""
    from fin.candidates.intron_chains import expand_canonical_chain_alternatives_v2

    if not arms:
        return []
    chrom = arms[0][1].chrom

    # Per-read original chain + minimal read dict for the wobble expander.
    read_chains: List[Tuple[Dict, IntronChain]] = []
    for qname, arm in arms:
        chain = extract_intron_chain(list(arm.cigartuples), arm.ref_start)
        read_chains.append(({"query_name": qname}, chain))

    alts: Dict[str, List[IntronChain]] = {}
    if search_bp > 0 and genome_seq and motif_set:
        alts = expand_canonical_chain_alternatives_v2(
            read_chains, genome_seq, strand, motif_set,
            search_bp=search_bp, max_chains_per_read=max_chains_per_read,
        )

    # A read supports every chain in its alternative set (its observed chain plus
    # canonical wobble shifts). Group reads by chain.
    chain_reads: Dict[IntronChain, Set[str]] = {}
    span: Dict[IntronChain, Tuple[int, int]] = {}
    for (rd, chain), (qname, arm) in zip(read_chains, arms):
        chains_for_read = alts.get(qname, [chain])
        for ch in chains_for_read:
            chain_reads.setdefault(ch, set()).add(qname)
            s, e = span.get(ch, (arm.ref_start, arm.ref_end))
            span[ch] = (min(s, arm.ref_start), max(e, arm.ref_end))

    variants: List[ArmVariant] = []
    for ch, qnames in chain_reads.items():
        if len(qnames) < min_support:
            continue
        s, e = span[ch]
        variants.append(ArmVariant(
            chrom=chrom, strand=strand, start=s, end=e,
            intron_chain=ch, source="read", supporting_read_ids=set(qnames),
        ))
    return variants


def _annotation_variants(
    chrom: str,
    strand: str,
    region_start: int,
    region_end: int,
    gtf_reader,
) -> List[ArmVariant]:
    """Annotation chains overlapping the arm region, as parallel candidate options."""
    if gtf_reader is None:
        return []
    try:
        txs = gtf_reader.get_transcripts_in_region(chrom, region_start, region_end)
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning("Fusion annotation lookup failed for %s: %s", chrom, exc)
        return []
    out: List[ArmVariant] = []
    for tx in txs:
        if tx.strand != strand:
            continue
        tx.sort_features()
        out.append(ArmVariant(
            chrom=chrom, strand=strand, start=tx.start, end=tx.end,
            intron_chain=gtf_transcript_to_intron_chain(tx),
            source="annotation", supporting_read_ids=set(),
        ))
    return out


def infer_arm_variants(
    cluster: FusionPairCluster,
    side: str,
    genome_by_chrom: Dict[str, str],
    gtf_reader=None,
    motif_set=None,
    search_bp: int = 4,
    max_chains_per_read: int = 16,
    min_support: int = 2,
) -> List[ArmVariant]:
    """Infer one arm's candidate splice variants (read-derived + annotation).

    Args:
        cluster: A :class:`FusionPairCluster` from :func:`cluster_chimeric_reads`.
        side: "a" (primary arm) or "b" (partner arm).
        genome_by_chrom: chrom -> full sequence string.
        gtf_reader: Optional opened GTF reader (annotation candidate source).
        motif_set: Canonical (donor, acceptor) motif set for wobble expansion.
        search_bp: Wobble window radius (0 disables).
        max_chains_per_read: Cap on wobble alternatives per read.
        min_support: Minimum reads for a READ-derived variant (annotation
            variants are exempt — they are offered as options regardless).

    Returns:
        Combined list of read-derived and annotation :class:`ArmVariant`.
    """
    arms: List[Tuple[str, ArmAlignment]] = []
    for r in cluster.reads:
        arm = r.arm_a if side == "a" else r.arm_b
        arms.append((r.query_name, arm))
    if not arms:
        return []
    chrom = arms[0][1].chrom
    strand = arms[0][1].strand
    genome_seq = genome_by_chrom.get(chrom, "") if genome_by_chrom else ""

    read_vars = _read_derived_variants(
        arms, genome_seq, strand, motif_set, search_bp,
        max_chains_per_read, min_support,
    )

    region_start = min(a.ref_start for _q, a in arms)
    region_end = max(a.ref_end for _q, a in arms)
    annot_vars = _annotation_variants(chrom, strand, region_start, region_end, gtf_reader)

    return read_vars + annot_vars


def assemble_fusion_arms(
    chimeric_reads: List[ChimericRead],
    genome_by_chrom: Dict[str, str],
    gtf_reader=None,
    motif_set=None,
    max_dist: int = 500,
    search_bp: int = 4,
    max_chains_per_read: int = 16,
    min_support: int = 2,
) -> List[FusionPairCluster]:
    """End-to-end F2: cluster chimeric reads and populate per-arm variant pools."""
    clusters = cluster_chimeric_reads(chimeric_reads, max_dist=max_dist)
    for cluster in clusters:
        cluster.arm_a_variants = infer_arm_variants(
            cluster, "a", genome_by_chrom, gtf_reader, motif_set,
            search_bp, max_chains_per_read, min_support,
        )
        cluster.arm_b_variants = infer_arm_variants(
            cluster, "b", genome_by_chrom, gtf_reader, motif_set,
            search_bp, max_chains_per_read, min_support,
        )
    return clusters
