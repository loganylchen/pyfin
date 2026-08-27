"""Run the normal CLI while auditing isoform-fraction cross-family suppression.

Set ``PYFIN_FAMILY_AUDIT`` to the JSON output path. This experiment-only wrapper
reconstructs discovery-like families on the exact pre-filter QuantResult set and
never changes the production drop set.
"""

from __future__ import annotations

from collections import defaultdict
import json
import os
from pathlib import Path
from typing import Iterable

from fin.analysis.quantification import _qr_overlap
from fin.candidates.chain_cluster import cluster_families
from fin.candidates.dataclasses import IntronChain
import fin.pipeline.selection as selection

_ORIGINAL = selection.isoform_fraction_drops


def _introns(qr) -> tuple[tuple[int, int], ...]:
    return tuple(
        (int(qr.exons[i][1]), int(qr.exons[i + 1][0]))
        for i in range(len(qr.exons) - 1)
    )


def _overlap_components(items: Iterable[object]) -> list[list[object]]:
    components: list[list[object]] = []
    current: list[object] = []
    current_end = -1
    for qr in sorted(items, key=lambda x: (x.start, x.end, x.candidate_id)):
        if current and qr.start >= current_end:
            components.append(current)
            current = []
            current_end = -1
        current.append(qr)
        current_end = max(current_end, qr.end)
    if current:
        components.append(current)
    return components


def _family_map(results) -> dict[str, str]:
    buckets = defaultdict(list)
    for qr in results.values():
        if len(qr.exons) >= 2 and qr.source != "fusion":
            buckets[(qr.chrom, qr.strand)].append(qr)

    mapping: dict[str, str] = {}
    for (chrom, strand), bucket in buckets.items():
        for component_index, component in enumerate(_overlap_components(bucket)):
            read_chains = [
                ({"query_name": qr.candidate_id}, IntronChain(_introns(qr)))
                for qr in component
                if qr.source == "novel"
            ]
            gtf_variants = [
                (qr.candidate_id, _introns(qr))
                for qr in component
                if qr.source == "gtf"
            ]
            clustered = cluster_families(
                read_chains,
                wobble_bp=6,
                cassette_max_exon_bp=70,
                gtf_variants=gtf_variants,
            )
            for family_index, family in enumerate(clustered.families):
                family_id = (
                    f"{chrom}:{strand}:{component_index}:{family_index}"
                )
                for candidate_id in family.read_pool:
                    mapping[candidate_id] = family_id
                for candidate_id, _chain in family.gtf_members:
                    mapping[candidate_id] = family_id
    return mapping


def _audit_isoform_fraction(results, min_fraction, locus="family"):
    drops = _ORIGINAL(results, min_fraction, locus=locus)
    output = os.environ.get("PYFIN_FAMILY_AUDIT")
    if not output or not drops:
        return drops

    family = _family_map(results)
    buckets = defaultdict(list)
    for qr in results.values():
        buckets[(qr.chrom, qr.strand)].append(qr)

    details = []
    for candidate_id in sorted(drops):
        qr = results[candidate_id]
        competitors = [
            other
            for other in buckets[(qr.chrom, qr.strand)]
            if other.candidate_id != candidate_id
            and _qr_overlap(qr, other)
            and other.abundance > qr.abundance
        ]
        dominant_abundance = max(other.abundance for other in competitors)
        dominants = [
            other
            for other in competitors
            if other.abundance == dominant_abundance
        ]
        candidate_family = family.get(candidate_id)
        same_family = any(
            candidate_family is not None
            and family.get(other.candidate_id) == candidate_family
            for other in dominants
        )
        details.append({
            "candidate_id": candidate_id,
            "candidate_abundance": qr.abundance,
            "candidate_family": candidate_family,
            "dominant_abundance": dominant_abundance,
            "dominant_ids": [other.candidate_id for other in dominants],
            "dominant_sources": sorted({other.source for other in dominants}),
            "dominant_families": sorted({
                family.get(other.candidate_id) or "family-less"
                for other in dominants
            }),
            "same_family_dominant": same_family,
        })

    payload = {
        "method": (
            "pre-isoform survivors reclustered with cluster_families; "
            "family IDs were not persisted by production"
        ),
        "min_fraction": min_fraction,
        "locus": locus,
        "total_drops": len(details),
        "same_family_dominant": sum(
            item["same_family_dominant"] for item in details
        ),
        "cross_family_only": sum(
            not item["same_family_dominant"] for item in details
        ),
        "details": details,
    }
    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return drops


selection.isoform_fraction_drops = _audit_isoform_fraction


if __name__ == "__main__":
    from fin.cli import main

    main()
