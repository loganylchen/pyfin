"""Build fusion TranscriptCandidate objects from breakpoint clusters."""

from __future__ import annotations

import logging
from typing import Dict, List

logger = logging.getLogger(__name__)

try:
    from fin.fusion.breakpoints import Breakpoint
except ImportError:
    # Stub for when the parallel agent's breakpoints.py is not yet present.
    from dataclasses import dataclass, field
    from typing import Set

    @dataclass
    class Breakpoint:  # type: ignore[no-redef]
        chromA: str
        posA: int
        strandA: str
        chromB: str
        posB: int
        strandB: str
        support_count: int = 1
        supporting_read_ids: Set[str] = field(default_factory=set)

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate

# Complement table for reverse-complement helper
_COMPLEMENT: Dict[str, str] = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "a": "t", "t": "a", "c": "g", "g": "c",
    "N": "N", "n": "n",
}


def _revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return "".join(_COMPLEMENT.get(b, "N") for b in reversed(seq))


def build_fusion_candidates(
    clusters: List[Breakpoint],
    genome_fasta: Dict[str, str],
    gtf_reader=None,  # unused API stub for future gene-aware fusion naming
    flank_bp: int = 500,
) -> List[TranscriptCandidate]:
    """Build fusion TranscriptCandidate objects from a list of Breakpoint clusters.

    For each Breakpoint, extracts flanking genomic sequence from each side of
    the fusion junction, concatenates them, and returns a TranscriptCandidate
    with source='fusion'.

    Args:
        clusters: List of Breakpoint objects describing fusion junctions.
        genome_fasta: Mapping of chromosome name to sequence string (or any
            object supporting slice notation, e.g. a pyfaidx.Fasta handle).
        gtf_reader: Optional; currently unused. Reserved for future gene-aware
            fusion naming.
        flank_bp: Number of bases to extract on each side of the junction.

    Returns:
        List of TranscriptCandidate objects with source='fusion'.
    """
    candidates: List[TranscriptCandidate] = []

    for bp in clusters:
        # --- Left side: flank_bp bases upstream of posA ---
        left_start = max(0, bp.posA - flank_bp)
        seq_a = genome_fasta.get(bp.chromA)
        if seq_a is None:
            logger.warning("Fusion chrom %s not in genome_fasta; skipping breakpoint", bp.chromA)
            continue
        left_seq: str = str(seq_a[left_start : bp.posA])
        if bp.strandA == "-":
            left_seq = _revcomp(left_seq)

        # --- Right side: flank_bp bases downstream of posB ---
        right_end = bp.posB + flank_bp
        seq_b = genome_fasta.get(bp.chromB)
        if seq_b is None:
            logger.warning("Fusion chrom %s not in genome_fasta; skipping breakpoint", bp.chromB)
            continue
        right_seq: str = str(seq_b[bp.posB : right_end])
        if bp.strandB == "-":
            right_seq = _revcomp(right_seq)

        sequence = left_seq + right_seq

        candidate_id = f"fusion_{bp.chromA}_{bp.posA}_{bp.chromB}_{bp.posB}"

        tc = TranscriptCandidate(
            candidate_id=candidate_id,
            intron_chain=IntronChain(introns=()),
            three_prime_pos=bp.posB,
            sequence=sequence,
            source="fusion",
            supporting_read_ids=set(sorted(bp.supporting_read_ids)),
            chrom=f"{bp.chromA}::{bp.chromB}",
            strand=".",
            start=bp.posA,
            end=bp.posB,
            fusion_junction=(bp.posA, bp.posB),
            breakpoint_left=(bp.chromA, bp.posA, bp.strandA),
            breakpoint_right=(bp.chromB, bp.posB, bp.strandB),
        )
        candidates.append(tc)

    return candidates
