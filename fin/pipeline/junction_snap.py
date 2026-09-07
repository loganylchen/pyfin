"""Read-supported junction consensus correction for finalized novel models."""
from __future__ import annotations

import logging
import math
from collections import Counter, defaultdict
from dataclasses import replace
from pathlib import Path
from typing import Callable, Dict, Optional, Sequence, Tuple

from fin.analysis.quantification import QuantResult
from fin.candidates.intron_chains import extract_intron_chain

logger = logging.getLogger(__name__)

JunctionSupport = Dict[Tuple[str, str], Counter]


def collect_junction_support(bam_path: str) -> Optional[JunctionSupport]:
    """Scan primary BAM records once and count exact junction observations."""
    if not bam_path or not Path(bam_path).exists():
        return None

    import pysam

    support: JunctionSupport = defaultdict(Counter)
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if (
                    read.is_unmapped
                    or read.is_secondary
                    or read.is_supplementary
                    or not read.cigartuples
                ):
                    continue
                strand = "-" if read.is_reverse else "+"
                chain = extract_intron_chain(
                    read.cigartuples, read.reference_start
                )
                for intron in chain.introns:
                    support[(read.reference_name, strand)][intron] += 1
    except Exception as exc:  # fail-open: correction is optional
        logger.warning(
            "Junction consensus scan failed (%s); snapping is disabled", exc
        )
        return None
    return support or None


def _snap_exons(
    exons: Tuple[Tuple[int, int], ...],
    support: Counter,
    *,
    tolerance: int,
    min_support: int,
    min_ratio: float,
    target_allowed: Optional[Callable[[Tuple[int, int]], bool]] = None,
) -> tuple[Tuple[Tuple[int, int], ...], int]:
    """Snap 0-based half-open exon junctions to stronger local read modes."""
    snapped = [list(exon) for exon in sorted(exons)]
    changed = 0
    for index in range(len(snapped) - 1):
        current = (snapped[index][1], snapped[index + 1][0])
        nearby = [
            (count, intron)
            for intron, count in support.items()
            if abs(intron[0] - current[0]) <= tolerance
            and abs(intron[1] - current[1]) <= tolerance
            and (target_allowed is None or target_allowed(intron))
        ]
        if not nearby:
            continue
        nearby.sort(
            key=lambda item: (
                -item[0],
                abs(item[1][0] - current[0])
                + abs(item[1][1] - current[1]),
                item[1],
            )
        )
        best_count, best = nearby[0]
        required = max(
            min_support,
            math.floor(support[current] * min_ratio) + 1,
        )
        if best == current or best_count < required:
            continue
        if not (snapped[index][0] < best[0] < best[1] < snapped[index + 1][1]):
            continue
        snapped[index][1] = best[0]
        snapped[index + 1][0] = best[1]
        changed += 1
    return tuple((start, end) for start, end in snapped), changed


def _introns(exons: Sequence[Tuple[int, int]]) -> Tuple[Tuple[int, int], ...]:
    ex = sorted(exons)
    return tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))


def _three_prime(qr: QuantResult) -> int:
    return qr.end if qr.strand == "+" else qr.start


def _five_prime(qr: QuantResult) -> int:
    return qr.start if qr.strand == "+" else qr.end


def _bounded_diameter_groups(
    items: Sequence[Tuple[int, QuantResult]], window: int
) -> list[list[QuantResult]]:
    """Partition by 5' coordinate so every group's DIAMETER is <= ``window``.

    Not connected components: chaining 0-20-40 through a 25 bp link would put
    endpoints 40 bp apart in one group. Sorted greedy from the smallest
    coordinate, closing a group as soon as the next item would exceed the
    diameter, is deterministic and gives the bound directly."""
    out: list[list[QuantResult]] = []
    cur: list[QuantResult] = []
    anchor = 0
    for pos, qr in sorted(items, key=lambda t: (t[0], t[1].candidate_id)):
        if cur and pos - anchor > window:
            out.append(cur)
            cur, anchor = [], pos
        if not cur:
            anchor = pos
        cur.append(qr)
    if cur:
        out.append(cur)
    return out


def _merge_group(group: list[QuantResult]) -> QuantResult:
    representative = max(
        group,
        key=lambda qr: (
            qr.num_assigned_reads,
            qr.abundance,
            qr.candidate_id,
        ),
    )
    assigned = tuple(sorted({rid for qr in group for rid in qr.assigned_read_ids}))
    hard_total = sum(qr.num_assigned_reads for qr in group)
    confidence = (
        sum(qr.confidence * qr.num_assigned_reads for qr in group) / hard_total
        if hard_total else 0.0
    )
    scored_fulllen = [
        qr for qr in group if qr.fulllen_frac >= 0.0 and qr.num_assigned_reads > 0
    ]
    fulllen_frac = (
        sum(qr.fulllen_frac * qr.num_assigned_reads for qr in scored_fulllen)
        / sum(qr.num_assigned_reads for qr in scored_fulllen)
        if scored_fulllen else -1.0
    )
    family_ids = sorted({qr.family_id for qr in group if qr.family_id is not None})
    return replace(
        representative,
        abundance=sum(qr.abundance for qr in group),
        confidence=confidence,
        num_assigned_reads=len(assigned) if assigned else hard_total,
        assigned_read_ids=assigned,
        family_id=family_ids[0] if family_ids else None,
        max_R=max((qr.max_R for qr in group), default=0.0),
        fulllen_frac=fulllen_frac,
    )


def snap_quant_results(
    results: Dict[str, QuantResult],
    observed: JunctionSupport,
    *,
    tolerance: int,
    min_support: int,
    min_ratio: float,
    genome_fasta: Optional[Dict[str, str]] = None,
    canonical_motifs: Sequence[str] = (),
    require_canonical: bool = False,
    return_redirects: bool = False,
    chain_merge: bool = False,
    chain_merge_window_bp: int = 25,
):
    """Correct novel junctions and merge models that become structurally equal.

    ``chain_merge`` additionally reconciles models that snapping gave the SAME
    CORRECTED INTRON CHAIN but not the same exon tuple. The default merge key is the full
    exon tuple, so two candidates that differed by a junction wobble -- and whose
    junctions snapping has just declared identical -- stay separate whenever
    their terminal coordinates differ at all. gffcompare then credits one of them
    and charges the other, which is a precision loss at zero recall cost.

    The extra merge is deliberately narrow, because "same corrected intron chain"
    is NOT "same transcript": genuine APA and alternative-TSS isoforms share an
    intron chain, and the chain-based scorer cannot see the difference, so a
    permissive rule would improve the number while destroying the distinction.
    An additional merge requires ALL of: every member novel and multi-exon on one
    chromosome and a known strand; the same post-snap intron chain; the collision
    was CREATED BY THE SNAP (the members' pre-snap chain sets are NON-EMPTY and
    pairwise disjoint, and at least one member was actually corrected); ONE non-null
    discovery family across every original member; identical strand-aware 3'
    coordinates; and a 5' spread no wider than ``chain_merge_window_bp``,
    measured as a group DIAMETER. These are operational equivalence criteria --
    they do not establish that two models are the same molecule species."""
    if chain_merge and chain_merge_window_bp < 0:
        raise ValueError(
            f"chain_merge_window_bp must be >= 0, got {chain_merge_window_bp}")
    motif_set = None
    if require_canonical:
        from fin.candidates.canonical import parse_motifs

        motif_set = parse_motifs(canonical_motifs)
    corrected: list[QuantResult] = []
    snapped_junctions = 0
    pre_chain: Dict[str, Tuple[Tuple[int, int], ...]] = {}
    snapped_ids: set = set()
    for qr in results.values():
        if qr.source != "novel" or len(qr.exons) < 2:
            corrected.append(qr)
            continue
        if chain_merge:
            pre_chain[qr.candidate_id] = _introns(qr.exons)
        target_allowed = None
        if require_canonical:
            from fin.candidates.canonical import chain_all_canonical

            genome = (genome_fasta or {}).get(qr.chrom, "")
            target_allowed = lambda intron: chain_all_canonical(
                (intron,), genome, qr.strand, motif_set
            )
        exons, changed = _snap_exons(
            qr.exons,
            observed.get((qr.chrom, qr.strand), Counter()),
            tolerance=tolerance,
            min_support=min_support,
            min_ratio=min_ratio,
            target_allowed=target_allowed,
        )
        corrected.append(replace(qr, exons=exons) if changed else qr)
        if changed and chain_merge:
            snapped_ids.add(qr.candidate_id)
        snapped_junctions += changed

    groups: dict[tuple, list[QuantResult]] = defaultdict(list)
    for qr in corrected:
        if qr.source == "novel" and len(qr.exons) >= 2:
            key = (qr.source, qr.chrom, qr.strand, qr.exons)
        else:
            # Exemption covers structural merging too: annotated/fusion/mono
            # models retain their identity even when coordinates coincide.
            key = (qr.source, qr.candidate_id)
        groups[key].append(qr)

    merged: Dict[str, QuantResult] = {}
    redirects: Dict[str, str] = {}
    members_of: Dict[str, list] = {}
    # Provenance of each first-pass survivor, kept per survivor and NOT pooled:
    # the pre-snap chains its members carried, and the discovery families they
    # came from (None included, so a missing family cannot be hidden by a
    # sibling that has one). The second pass reads these; the first pass's own
    # behaviour is untouched.
    origin_chains: Dict[str, frozenset] = {}
    origin_families: Dict[str, frozenset] = {}
    merged_models = 0
    for group in groups.values():
        result = _merge_group(group)
        merged[result.candidate_id] = result
        if chain_merge:
            members_of[result.candidate_id] = [m.candidate_id for m in group]
            origin_chains[result.candidate_id] = frozenset(
                pre_chain[m.candidate_id] for m in group
                if m.candidate_id in pre_chain)
            origin_families[result.candidate_id] = frozenset(
                m.family_id for m in group)
        for member in group:
            if member.candidate_id != result.candidate_id:
                redirects[member.candidate_id] = result.candidate_id
        merged_models += len(group) - 1
    if chain_merge:
        merged, extra, redirects = _merge_snap_collisions(
            merged, members_of, redirects, origin_chains, origin_families,
            snapped_ids, window=chain_merge_window_bp)
        merged_models += extra
    if return_redirects:
        return merged, snapped_junctions, merged_models, redirects
    return merged, snapped_junctions, merged_models


def _merge_snap_collisions(
    merged: Dict[str, QuantResult],
    members_of: Dict[str, list],
    redirects: Dict[str, str],
    origin_chains: Dict[str, frozenset],
    origin_families: Dict[str, frozenset],
    snapped_ids: set,
    *,
    window: int,
):
    """Second merge pass: models the SNAP gave one corrected intron chain.

    See :func:`snap_quant_results` for the rule and why it is this narrow. The
    survivor keeps the representative's own post-snap exon tuple and geometry, so
    no new endpoint estimate is invented here: the pass only decides which models
    are treated as one, and never where a transcript ends."""
    buckets: Dict[tuple, list] = defaultdict(list)
    passthrough: Dict[str, QuantResult] = {}
    for cid, qr in merged.items():
        fams = origin_families.get(cid, frozenset())
        # ONE non-null discovery family among the ORIGINAL members. The merged
        # survivor's own family_id cannot be used: `_merge_group` picks the
        # lexicographically first non-null family of its group, so a first-pass
        # merge across two families (or across a member with no family) would
        # otherwise present itself here as single-family.
        if (qr.source == "novel" and len(qr.exons) >= 2
                and qr.chrom and qr.strand in ("+", "-")
                and len(fams) == 1 and None not in fams):
            key = (qr.chrom, qr.strand, _introns(qr.exons), next(iter(fams)),
                   _three_prime(qr))
            buckets[key].append(qr)
        else:
            passthrough[cid] = qr

    out: Dict[str, QuantResult] = dict(passthrough)
    extra = 0
    for members in buckets.values():
        if len(members) == 1:
            out[members[0].candidate_id] = members[0]
            continue
        for group in _bounded_diameter_groups(
                [(_five_prime(qr), qr) for qr in members], window):
            if len(group) == 1:
                out[group[0].candidate_id] = group[0]
                continue
            # Provenance: the snap must have created the relationship between
            # EVERY pair being collapsed, not merely somewhere in the group. Two
            # models that already shared a chain before the snap are a
            # pre-existing endpoint distinction, and a third wobble model must
            # not supply the permission to collapse them -- so the survivors'
            # pre-snap chain sets have to be pairwise DISJOINT.
            sets = [origin_chains.get(qr.candidate_id, frozenset()) for qr in group]
            origin = [c for qr in group for c in members_of.get(qr.candidate_id, ())]
            disjoint = all(
                not (sets[i] & sets[j])
                for i in range(len(sets)) for j in range(i + 1, len(sets)))
            if (not disjoint or not all(sets)
                    or not any(c in snapped_ids for c in origin)):
                for qr in group:
                    out[qr.candidate_id] = qr
                continue
            result = _merge_group(group)
            # fulllen_frac is a fraction of reads full-length wrt THIS model's
            # ends; the abundance refit does not recompute it, and averaging
            # values measured against different ends over a grown read set would
            # ship a stale number. Mark it unavailable instead.
            result = replace(result, fulllen_frac=-1.0)
            out[result.candidate_id] = result
            for qr in group:
                if qr.candidate_id != result.candidate_id:
                    for c in members_of.get(qr.candidate_id, ()):
                        redirects[c] = result.candidate_id
                    redirects[qr.candidate_id] = result.candidate_id
            extra += len(group) - 1
    # absorbed survivors must not remain as redirect targets
    for src, dst in list(redirects.items()):
        seen = {src}
        while dst in redirects and dst not in seen:
            seen.add(dst)
            dst = redirects[dst]
        redirects[src] = dst
    return out, extra, redirects


def apply_junction_snap(
    config,
    results: Dict[str, QuantResult],
    genome_fasta: Optional[Dict[str, str]] = None,
    *,
    return_redirects: bool = False,
):
    """Apply the optional finalized-model correction, failing open on BAM errors."""
    if not getattr(config, "junction_snap", False):
        return (results, {}) if return_redirects else results
    observed = collect_junction_support(config.bam_path)
    if not observed:
        return (results, {}) if return_redirects else results
    snap_output = snap_quant_results(
        results,
        observed,
        tolerance=int(config.junction_snap_tolerance),
        min_support=int(config.junction_snap_min_support),
        min_ratio=float(config.junction_snap_min_ratio),
        genome_fasta=genome_fasta,
        canonical_motifs=getattr(config, "canonical_motifs", ()),
        require_canonical=bool(getattr(config, "canonical_gate", False)),
        return_redirects=return_redirects,
        chain_merge=bool(getattr(config, "snap_merge_chain", False)),
        chain_merge_window_bp=int(getattr(config, "snap_merge_chain_window_bp", 25)),
    )
    if return_redirects:
        corrected, snapped, merged, redirects = snap_output
    else:
        corrected, snapped, merged = snap_output
        redirects = {}
    logger.info(
        "Junction consensus snapped %d junctions and merged %d duplicate models",
        snapped,
        merged,
    )
    return (corrected, redirects) if return_redirects else corrected
