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
):
    """Correct novel junctions and merge models that become structurally equal."""
    motif_set = None
    if require_canonical:
        from fin.candidates.canonical import parse_motifs

        motif_set = parse_motifs(canonical_motifs)
    corrected: list[QuantResult] = []
    snapped_junctions = 0
    for qr in results.values():
        if qr.source != "novel" or len(qr.exons) < 2:
            corrected.append(qr)
            continue
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
    merged_models = 0
    for group in groups.values():
        result = _merge_group(group)
        merged[result.candidate_id] = result
        for member in group:
            if member.candidate_id != result.candidate_id:
                redirects[member.candidate_id] = result.candidate_id
        merged_models += len(group) - 1
    if return_redirects:
        return merged, snapped_junctions, merged_models, redirects
    return merged, snapped_junctions, merged_models


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
