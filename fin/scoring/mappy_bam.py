"""Per-candidate mappy alignment → BAM.

Relocated from the (now-removed) f5c external-tools pipeline. The only remaining
consumer is the multi-sample quantify runner, which needs a per-candidate
``aligned.bam`` so the final assignment BAM can be assembled from EM hard
assignments. All signal scoring is in-memory krill; this module is alignment
(structure) only.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
from typing import List, Optional

from fin.candidates.dataclasses import TranscriptCandidate

logger = logging.getLogger(__name__)

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(s: str) -> str:
    return s.translate(_COMPLEMENT)[::-1]


def align_reads_with_mappy(
    candidate: TranscriptCandidate,
    fastq_path: Path,
    bam_path: Path,
    cached_reads: Optional[List[tuple]] = None,
) -> int:
    """Align all reads to a single candidate using mappy, write sorted+indexed BAM.

    Args:
        cached_reads: Optional pre-loaded list of (name, seq, qual) tuples.
            When provided, the FASTQ is not re-read from disk (avoids O(C * N)
            IO when aligning many candidates per interval).

    Returns:
        Number of primary alignments written.
    """
    import mappy
    import pysam

    # Build mappy index for this candidate's sequence
    aligner = mappy.Aligner(seq=candidate.sequence, preset="map-ont")
    if not aligner:
        raise RuntimeError(
            f"Failed to build mappy index for {candidate.candidate_id}"
        )

    # Parse FASTQ, align each read, collect alignments
    if cached_reads is not None:
        read_iter = cached_reads
    else:
        read_iter = mappy.fastx_read(str(fastq_path))

    alignments = []
    for name, seq, qual in read_iter:
        for hit in aligner.map(seq):
            if hit.is_primary:
                alignments.append((name, seq, qual, hit))

    logger.info(
        "mappy aligned %d reads to candidate %s",
        len(alignments),
        candidate.candidate_id,
    )

    # Write unsorted BAM first, then sort
    unsorted_bam = str(bam_path) + ".unsorted.bam"
    header = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "SQ": [
                {"SN": candidate.candidate_id, "LN": len(candidate.sequence)}
            ],
        }
    )

    with pysam.AlignmentFile(unsorted_bam, "wb", header=header) as outf:
        for name, seq, qual, hit in alignments:
            # mappy returns cigar as (length, op); pysam wants (op, length).
            # The CIGAR only covers the aligned region [hit.q_st, hit.q_en);
            # we must add soft-clip ops so CIGAR length matches SEQ length.
            cigar_ops = [(op, length) for length, op in hit.cigar]
            read_len = len(seq)

            if hit.strand == 1:
                bam_seq = seq
                bam_qual = qual
                left_clip = hit.q_st
                right_clip = read_len - hit.q_en
            else:
                # Reverse strand: BAM stores reverse complement; quals reversed.
                bam_seq = _revcomp(seq)
                bam_qual = qual[::-1] if qual else qual
                left_clip = read_len - hit.q_en
                right_clip = hit.q_st

            if left_clip > 0:
                cigar_ops.insert(0, (4, left_clip))
            if right_clip > 0:
                cigar_ops.append((4, right_clip))

            a = pysam.AlignedSegment(header)
            a.query_name = name
            a.query_sequence = bam_seq
            a.flag = 0 if hit.strand == 1 else 16
            a.reference_id = 0  # single contig
            a.reference_start = hit.r_st
            a.mapping_quality = hit.mapq
            a.cigar = cigar_ops
            if bam_qual:
                a.query_qualities = pysam.qualitystring_to_array(bam_qual)
            outf.write(a)

    # Sort and index
    pysam.sort("-o", str(bam_path), unsorted_bam)
    os.remove(unsorted_bam)
    pysam.index(str(bam_path))
    return len(alignments)
