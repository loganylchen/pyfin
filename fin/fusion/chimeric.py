"""Stage F1 of fusion detection: chimeric-read collection and arm re-alignment.

A fusion read is a chimeric read: one contiguous stretch aligns to one locus
(the primary alignment) and a long soft-clipped stretch aligns elsewhere (the
partner arm). Rather than trusting the aligner's SA tag, we re-align the
soft-clip segment against the genome to locate the partner arm independently.

This module also guards against the dominant nanopore dRNA false-fusion mode:
**adapter-bridged chimeras**, where two unrelated RNA molecules are ligated
through an internal ONT adapter. The adapter sequence is DNA and the RNA
basecaller mis-calls it, so it maps to neither arm and shows up as an internal
unmapped gap (~adapter length) between the two aligned arms. We drop any read
whose internal unmapped gap is >= ``max_internal_gap_bp`` (see DeepChopper,
Nat. Commun. 2026).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from fin.io.interval_manager import is_fusion_read
from fin.scoring.mappy_preset import get_m1_preset

logger = logging.getLogger(__name__)

# pysam CIGAR op for soft-clip.
_CIGAR_SOFT_CLIP = 4


@dataclass
class ArmAlignment:
    """One arm of a chimeric read's alignment.

    Coordinates:
        chrom/ref_start/ref_end/strand are GENOMIC (0-based, half-open).
        q_start/q_end are READ coordinates (0-based, half-open) on the forward
        query sequence — the span of the read this arm covers.
        cigartuples is the arm's (op, length) CIGAR in pysam convention.
    """

    chrom: str
    ref_start: int
    ref_end: int
    strand: str
    q_start: int
    q_end: int
    cigartuples: Tuple[Tuple[int, int], ...]


@dataclass
class ChimericRead:
    """A surviving chimeric read with both arms resolved.

    ``arm_a`` is the read's primary alignment, ``arm_b`` the re-aligned soft-clip
    partner. ``breakpoint_a``/``breakpoint_b`` are (chrom, pos, strand) at the
    fusion junction on each side. ``internal_gap`` is the number of read bases
    between the two arms that mapped to neither (0 for a clean junction).
    """

    query_name: str
    arm_a: ArmAlignment
    arm_b: ArmAlignment
    breakpoint_a: Tuple[str, int, str]
    breakpoint_b: Tuple[str, int, str]
    internal_gap: int


def build_genome_aligner(genome_fasta_path: str, preset: Optional[str] = None):
    """Build a genome-wide ``mappy.Aligner`` for re-aligning soft-clip segments.

    Returns None if mappy is unavailable or the index cannot be built, so the
    caller degrades gracefully (no fusion calls rather than a crash).
    """
    try:
        import mappy
    except ImportError:
        logger.warning("mappy not available; fusion arm re-alignment disabled")
        return None
    try:
        return mappy.Aligner(fn_idx_in=genome_fasta_path, preset=preset or get_m1_preset())
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning("Failed to build genome aligner for fusion (%s)", exc)
        return None


def _clip_segments(read: Dict) -> List[Tuple[str, int, int]]:
    """Return soft-clip segments as (side, q_start, q_end) in read coordinates.

    ``side`` is "front" (5' of the primary aligned span) or "back" (3' of it).
    Only segments meeting the soft-clip threshold in ``is_fusion_read`` matter,
    but we return every present soft-clip and let the caller pick the longest.
    """
    cigar = read.get("cigartuples")
    seq = read.get("query_sequence")
    if not cigar or not seq:
        return []
    q_aln_start = read.get("query_alignment_start")
    q_aln_end = read.get("query_alignment_end")
    if q_aln_start is None or q_aln_end is None:
        return []

    segments: List[Tuple[str, int, int]] = []
    first_op, first_len = cigar[0]
    last_op, last_len = cigar[-1]
    if first_op == _CIGAR_SOFT_CLIP and first_len > 0:
        segments.append(("front", 0, q_aln_start))
    if last_op == _CIGAR_SOFT_CLIP and last_len > 0:
        segments.append(("back", q_aln_end, len(seq)))
    return segments


def _best_hit(aligner, segment_seq: str):
    """Return the best mappy hit for ``segment_seq`` (max matching length), or None."""
    best = None
    best_mlen = -1
    for hit in aligner.map(segment_seq):
        mlen = getattr(hit, "mlen", hit.q_en - hit.q_st)
        if mlen > best_mlen:
            best_mlen = mlen
            best = hit
    return best


def _arm_a_breakpoint_pos(side: str, ref_start: int, ref_end: int) -> int:
    """Genomic junction edge of the PRIMARY arm.

    pysam normalizes cigartuples and query_alignment_start/end to GENOME-FORWARD
    orientation for every read (``is_reverse`` only records the sequenced strand),
    so the read-coordinate-to-reference mapping is monotonic regardless of strand:
    the aligned span's read-low end is ``ref_start`` and its read-high end is
    ``ref_end``. For a back soft-clip the partner arm sits at higher read coords,
    so the junction is the primary's read-high end (``ref_end``); a front clip is
    the mirror image (``ref_start``).
    """
    return ref_end if side == "back" else ref_start


def _arm_b_breakpoint_pos(side: str, strand: str, r_st: int, r_en: int) -> int:
    """Genomic junction edge of the PARTNER arm (nearest the junction).

    The partner arm is found by re-aligning the soft-clip SEGMENT with mappy, so
    its orientation is the mappy hit ``strand`` (genuinely strand-aware, unlike the
    pysam-normalized primary). For a back clip the partner lies downstream in read
    coords; the junction is at the partner edge adjacent to the primary — its
    read-low end, i.e. ``r_st`` on '+' and ``r_en`` on '-'. A front clip mirrors.
    """
    if side == "back":
        return r_en if strand == "-" else r_st
    return r_st if strand == "-" else r_en


def _internal_gap(side: str, q_aln_start: int, q_aln_end: int,
                  seg_q_start: int, hit_q_st: int, hit_q_en: int) -> int:
    """Unmapped read bases sandwiched between the primary arm and the partner arm.

    The partner alignment covers, in read coords, ``[seg_q_start + hit_q_st,
    seg_q_start + hit_q_en)``. The internal gap is the stretch between the two
    arms that neither covers (adapter signature). Bases at the very read ends
    (outside both arms) are NOT counted — only the internal sandwich matters.
    """
    partner_read_start = seg_q_start + hit_q_st
    partner_read_end = seg_q_start + hit_q_en
    if side == "back":
        # primary ends at q_aln_end; partner begins at partner_read_start (>=).
        return max(0, partner_read_start - q_aln_end)
    # front: partner ends at partner_read_end; primary begins at q_aln_start.
    return max(0, q_aln_start - partner_read_end)


def collect_chimeric_reads(
    read_dicts: List[Dict],
    genome_aligner,
    max_internal_gap_bp: int = 30,
) -> List[ChimericRead]:
    """Resolve chimeric reads into two-armed alignments, dropping adapter chimeras.

    Args:
        read_dicts: Alignment dicts from ``BamReader.alignment_to_dict``.
        genome_aligner: A ``mappy.Aligner`` (or any object exposing ``.map`` that
            yields hits with ``ctg, r_st, r_en, strand, q_st, q_en, cigar``) used
            to re-align soft-clip segments. If None, returns [].
        max_internal_gap_bp: Drop a read whose internal unmapped gap (read coords,
            between the two arms) is >= this. 0 disables the guard.

    Returns:
        List of :class:`ChimericRead` for reads that pass the adapter-chimera guard.
    """
    if genome_aligner is None:
        return []

    out: List[ChimericRead] = []
    for read in read_dicts:
        if not is_fusion_read(read):
            continue
        chrom_a = read.get("reference_name")
        ref_start = read.get("reference_start")
        ref_end = read.get("reference_end")
        cigar = read.get("cigartuples")
        seq = read.get("query_sequence")
        qname = read.get("query_name")
        if None in (chrom_a, ref_start, ref_end, qname) or not cigar or not seq:
            continue
        q_aln_start = read.get("query_alignment_start")
        q_aln_end = read.get("query_alignment_end")
        if q_aln_start is None or q_aln_end is None:
            continue
        strand_a = "-" if read.get("is_reverse") else "+"

        # Pick the longest soft-clip segment as the partner-bearing arm.
        segments = _clip_segments(read)
        if not segments:
            continue
        side, seg_q_start, seg_q_end = max(segments, key=lambda s: s[2] - s[1])
        segment_seq = seq[seg_q_start:seg_q_end]
        if not segment_seq:
            continue

        hit = _best_hit(genome_aligner, segment_seq)
        if hit is None:
            continue

        gap = _internal_gap(side, q_aln_start, q_aln_end,
                            seg_q_start, hit.q_st, hit.q_en)
        if max_internal_gap_bp > 0 and gap >= max_internal_gap_bp:
            logger.debug(
                "Dropping adapter-bridged chimera %s: internal gap %d bp >= %d",
                qname, gap, max_internal_gap_bp,
            )
            continue

        strand_b = "+" if getattr(hit, "strand", 1) >= 0 else "-"
        arm_a = ArmAlignment(
            chrom=chrom_a, ref_start=ref_start, ref_end=ref_end, strand=strand_a,
            q_start=q_aln_start, q_end=q_aln_end, cigartuples=tuple(cigar),
        )
        arm_b = ArmAlignment(
            chrom=hit.ctg, ref_start=hit.r_st, ref_end=hit.r_en, strand=strand_b,
            q_start=seg_q_start + hit.q_st, q_end=seg_q_start + hit.q_en,
            cigartuples=tuple(getattr(hit, "cigar", ()) or ()),
        )
        # Breakpoint coordinate on each arm: the genomic edge nearest the
        # junction (strand-aware — read-forward coords flip vs genome on '-').
        pos_a = _arm_a_breakpoint_pos(side, ref_start, ref_end)
        pos_b = _arm_b_breakpoint_pos(side, strand_b, hit.r_st, hit.r_en)
        out.append(ChimericRead(
            query_name=qname,
            arm_a=arm_a,
            arm_b=arm_b,
            breakpoint_a=(chrom_a, pos_a, strand_a),
            breakpoint_b=(hit.ctg, pos_b, strand_b),
            internal_gap=gap,
        ))
    return out
