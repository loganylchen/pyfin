"""Fusion detection orchestrator: chains Stages F1 -> F2 -> F3.

This is the single entry point the pipeline runner calls. It composes the
read-driven fusion sub-pipeline:

  F1 (chimeric.collect_chimeric_reads): collect soft-clip chimeric reads,
     re-align the soft-clip arm against the genome, drop adapter-bridged chimeras.
  F2 (arm_assembly.assemble_fusion_arms): cluster by partner pair and infer each
     arm's splice variants (read-derived + annotation, in parallel).
  F3 (stitch.build_fusion_candidates_v2): stitch arm variants across the
     breakpoint into composite fusion TranscriptCandidates.
"""

from __future__ import annotations

from typing import Dict, List

from fin.candidates.dataclasses import TranscriptCandidate
from fin.fusion.arm_assembly import assemble_fusion_arms
from fin.fusion.chimeric import collect_chimeric_reads
from fin.fusion.stitch import build_fusion_candidates_v2


def detect_fusion_candidates(
    read_dicts: List[Dict],
    genome_aligner,
    genome_by_chrom: Dict[str, str],
    *,
    gtf_reader=None,
    motif_set=None,
    max_internal_gap_bp: int = 30,
    max_dist: int = 500,
    search_bp: int = 4,
    max_chains_per_read: int = 16,
    min_support: int = 2,
) -> List[TranscriptCandidate]:
    """Run the F1->F2->F3 fusion sub-pipeline and return fusion candidates.

    Args:
        read_dicts: Alignment dicts (``BamReader.alignment_to_dict``) for the
            region under consideration.
        genome_aligner: A genome-wide ``mappy.Aligner`` (or None) used to
            re-align soft-clip arms; build with ``build_genome_aligner``.
        genome_by_chrom: chrom -> full sequence string (for arm sequences).
        gtf_reader: Optional opened GTF reader (annotation candidate source).
        motif_set: Canonical (donor, acceptor) motif set for wobble expansion.
        max_internal_gap_bp: Adapter-chimera guard threshold (0 disables).
        max_dist: Breakpoint-proximity window for clustering.
        search_bp: Canonical wobble window radius (0 disables).
        max_chains_per_read: Cap on per-read wobble alternatives.
        min_support: Minimum reads for a read-derived fusion combo.

    Returns:
        List of fusion ``TranscriptCandidate`` (possibly empty).
    """
    chimeric = collect_chimeric_reads(
        read_dicts, genome_aligner, max_internal_gap_bp=max_internal_gap_bp
    )
    if not chimeric:
        return []
    clusters = assemble_fusion_arms(
        chimeric,
        genome_by_chrom,
        gtf_reader=gtf_reader,
        motif_set=motif_set,
        max_dist=max_dist,
        search_bp=search_bp,
        max_chains_per_read=max_chains_per_read,
        min_support=min_support,
    )
    return build_fusion_candidates_v2(clusters, genome_by_chrom, min_support=min_support)
