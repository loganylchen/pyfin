"""Read-driven fusion transcript detection.

The fusion sub-pipeline treats chimeric (soft-clipped) reads as a mini-assembly:

  - :func:`collect_chimeric_reads` (F1) collects soft-clip reads, re-aligns the
    soft-clip arm against the genome, and drops adapter-bridged chimera artifacts.
  - :func:`assemble_fusion_arms` (F2) clusters by partner pair and infers each
    arm's splice variants (read-derived + annotation, in parallel).
  - :func:`build_fusion_candidates_v2` (F3) stitches arm variants across the
    breakpoint into composite fusion ``TranscriptCandidate`` objects.
  - :func:`detect_fusion_candidates` chains F1->F2->F3 (the runner entry point).
  - :func:`build_genome_aligner` builds the genome-wide mappy aligner F1 needs.
"""

from __future__ import annotations

from fin.fusion.arm_assembly import (
    ArmVariant,
    FusionPairCluster,
    assemble_fusion_arms,
    cluster_chimeric_reads,
    infer_arm_variants,
)
from fin.fusion.chimeric import (
    ArmAlignment,
    ChimericRead,
    build_genome_aligner,
    collect_chimeric_reads,
)
from fin.fusion.detect import detect_fusion_candidates
from fin.fusion.stitch import build_fusion_candidates_v2, stitch_cluster

__all__ = [
    "ArmAlignment",
    "ChimericRead",
    "build_genome_aligner",
    "collect_chimeric_reads",
    "ArmVariant",
    "FusionPairCluster",
    "cluster_chimeric_reads",
    "infer_arm_variants",
    "assemble_fusion_arms",
    "stitch_cluster",
    "build_fusion_candidates_v2",
    "detect_fusion_candidates",
]
