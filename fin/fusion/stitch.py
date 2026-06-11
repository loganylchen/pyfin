"""Stage F3 of fusion detection: cross-breakpoint stitching into candidates.

Takes the per-arm splice-variant pools from Stage F2 and stitches each
(arm-A variant x arm-B variant) combination across the breakpoint into a
composite fusion ``TranscriptCandidate(source="fusion")`` carrying:

  * a spliced cDNA sequence = spliced(arm A) ++ spliced(arm B) in transcribed
    orientation (each arm reverse-complemented per its own strand);
  * the consensus breakpoint on each side (``breakpoint_left``/``breakpoint_right``,
    ``fusion_junction``);
  * the union/intersection of supporting chimeric reads (see support model below).

Distinct splice combinations become distinct candidate variants; the signal +
EM machinery (unchanged) arbitrates between them downstream, exactly as it does
for competing assembly isoforms.

Support model (parallel to assembly's treatment of GTF vs novel):
  * read x read combo: supported by the chimeric reads that match BOTH arm
    structures (read-set intersection); kept only if >= ``min_support``.
  * any combo whose arm came from annotation: kept regardless of read support
    (annotation is an additional hypothesis for EM to test, like GTF passthrough).

This module replaces the old ``build_fusion_candidates``.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Set, Tuple

from fin.candidates.dataclasses import TranscriptCandidate
from fin.candidates.discovery import _build_spliced_sequence
from fin.fusion.arm_assembly import ArmVariant, FusionPairCluster

logger = logging.getLogger(__name__)


def _combo_support(
    var_a: ArmVariant,
    var_b: ArmVariant,
    cluster_read_ids: Set[str],
    min_support: int,
) -> Tuple[Set[str], bool]:
    """Return (supporting_read_ids, keep) for one (var_a, var_b) combination.

    read x read: intersection of arm read sets; keep iff >= min_support.
    annotation participating: exempt from min_support (kept), supported by the
    read-derived arm's reads, or the whole cluster when both arms are annotation.
    """
    a_read = var_a.source == "read"
    b_read = var_b.source == "read"
    if a_read and b_read:
        support = set(var_a.supporting_read_ids) & set(var_b.supporting_read_ids)
        return support, len(support) >= min_support
    if a_read:
        return set(var_a.supporting_read_ids), True
    if b_read:
        return set(var_b.supporting_read_ids), True
    # both annotation: circumstantial cluster-level support
    return set(cluster_read_ids), True


def stitch_cluster(
    cluster: FusionPairCluster,
    genome_by_chrom: Dict[str, str],
    min_support: int = 2,
) -> List[TranscriptCandidate]:
    """Stitch one cluster's arm-variant pools into fusion candidates."""
    if not cluster.arm_a_variants or not cluster.arm_b_variants:
        return []

    chrom_a, pos_a, strand_a = cluster.breakpoint_a
    chrom_b, pos_b, strand_b = cluster.breakpoint_b
    seq_a_chrom = genome_by_chrom.get(chrom_a, "") if genome_by_chrom else ""
    seq_b_chrom = genome_by_chrom.get(chrom_b, "") if genome_by_chrom else ""
    cluster_read_ids = {r.query_name for r in cluster.reads}

    # Dedup by (armA chain, armB chain); union supporting reads on collision.
    by_key: Dict[Tuple[Tuple, Tuple], TranscriptCandidate] = {}

    # Deterministic iteration: sort variants by (source, chain, start).
    def _vkey(v: ArmVariant):
        return (v.source, v.intron_chain.introns, v.start, v.end)

    a_sorted = sorted(cluster.arm_a_variants, key=_vkey)
    b_sorted = sorted(cluster.arm_b_variants, key=_vkey)

    for var_a in a_sorted:
        seq_a = _build_spliced_sequence(
            seq_a_chrom, var_a.start, var_a.end, var_a.intron_chain, var_a.strand
        )
        for var_b in b_sorted:
            support, keep = _combo_support(var_a, var_b, cluster_read_ids, min_support)
            if not keep:
                continue
            seq_b = _build_spliced_sequence(
                seq_b_chrom, var_b.start, var_b.end, var_b.intron_chain, var_b.strand
            )
            # A fusion requires BOTH arms; if either side fails to build (chrom
            # absent from genome_by_chrom, or unusable coordinates) the candidate
            # would be a half-fusion — drop it rather than emit a partial sequence.
            if not seq_a or not seq_b:
                continue
            sequence = seq_a + seq_b

            key = (var_a.intron_chain.introns, var_b.intron_chain.introns)
            existing = by_key.get(key)
            if existing is not None:
                existing.supporting_read_ids.update(support)
                continue

            cand = TranscriptCandidate(
                candidate_id="",  # assigned after dedup for stable numbering
                intron_chain=var_a.intron_chain,  # arm-A structure (decorative;
                # fusion is exempt from chain gates and length uses sequence)
                three_prime_pos=pos_b,
                sequence=sequence,
                source="fusion",
                supporting_read_ids=set(support),
                chrom=f"{chrom_a}::{chrom_b}",
                strand=".",
                start=pos_a,
                end=pos_b,
                fusion_junction=(pos_a, pos_b),
                breakpoint_left=(chrom_a, pos_a, strand_a),
                breakpoint_right=(chrom_b, pos_b, strand_b),
            )
            by_key[key] = cand

    # Stable candidate IDs: sort by (armA chain, armB chain) then index.
    out: List[TranscriptCandidate] = []
    for idx, key in enumerate(sorted(by_key.keys())):
        cand = by_key[key]
        cand.candidate_id = f"fusion_{chrom_a}_{pos_a}_{chrom_b}_{pos_b}_v{idx}"
        out.append(cand)
    return out


def build_fusion_candidates_v2(
    clusters: List[FusionPairCluster],
    genome_by_chrom: Dict[str, str],
    min_support: int = 2,
) -> List[TranscriptCandidate]:
    """Stitch every cluster's arm variants into fusion candidates (Stage F3 entry)."""
    out: List[TranscriptCandidate] = []
    for cluster in clusters:
        out.extend(stitch_cluster(cluster, genome_by_chrom, min_support=min_support))
    return out
