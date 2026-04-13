"""Candidate transcript discovery per genomic interval."""

from __future__ import annotations

import hashlib
import logging
from typing import List, Optional, Set

from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.candidates.intron_chains import (
    extract_intron_chain,
    group_reads_by_three_prime_and_intron_chain,
    gtf_transcript_to_intron_chain,
    pick_representative_read,
)
from fin.io.interval_manager import (
    GenomicInterval,
    extract_strand_from_read,
    extract_three_prime_pos,
    is_fusion_read,
)

logger = logging.getLogger(__name__)


def _generate_novel_id() -> str:
    """Generate a short unique ID for a novel candidate."""
    import uuid

    return "novel_" + uuid.uuid4().hex[:8]


def _intron_chains_match(a: IntronChain, b: IntronChain) -> bool:
    """Check if two intron chains are identical."""
    return a.introns == b.introns


def _three_prime_within_threshold(pos_a: int, pos_b: int, threshold: int) -> bool:
    return abs(pos_a - pos_b) <= threshold


def discover_candidates(
    interval: GenomicInterval,
    bam_path: str,
    gtf_reader,
    genome_fasta: str,
    threshold: int = 24,
) -> CandidateSet:
    """Discover transcript candidates for a genomic interval.

    Combines GTF-annotated transcripts with novel isoforms discovered
    from read alignments using intron chain + 3' end clustering.

    Args:
        interval: Genomic interval to process.
        bam_path: Path to BAM file.
        gtf_reader: Opened GTFReader instance (already parsed).
        genome_fasta: Full chromosome sequence string for the interval's chromosome.
        threshold: 3' end clustering distance in bp.

    Returns:
        CandidateSet with discovered candidates.
    """
    from fin.io.io_bam import BamReader

    # 1. Fetch reads for interval, skip fusion reads
    all_read_ids: Set[str] = set()
    non_fusion_reads = []

    with BamReader(bam_path) as bam:
        region = interval.region_string
        for alignment in bam.fetch(
            reference=interval.chrom, start=interval.start, end=interval.end
        ):
            rd = bam.alignment_to_dict(alignment)
            if not rd or not rd.get("is_mapped"):
                continue
            read_id = rd.get("query_name")
            if read_id:
                all_read_ids.add(read_id)
            if is_fusion_read(rd):
                continue
            non_fusion_reads.append(rd)

    if not non_fusion_reads:
        return CandidateSet(interval=interval, candidates=[], read_ids=all_read_ids)

    # Determine strand from reads
    strand = interval.strand or "+"

    # 2. Group reads by 3' + intron chain → novel candidate groups
    groups = group_reads_by_three_prime_and_intron_chain(
        non_fusion_reads, strand, threshold=threshold
    )

    # 3. Get GTF transcripts in region
    gtf_transcripts = []
    if gtf_reader is not None:
        gtf_transcripts = gtf_reader.get_transcripts_in_region(
            interval.chrom, interval.start, interval.end
        )

    # 4. Build GTF candidates with intron chains
    gtf_candidates = []
    for tx in gtf_transcripts:
        chain = gtf_transcript_to_intron_chain(tx)
        tx.sort_features()
        # Get 3' position from GTF
        if tx.strand == "+":
            three_prime = tx.end
        else:
            three_prime = tx.start

        seq = tx.get_spliced_sequence(genome_fasta)

        gtf_candidates.append(
            TranscriptCandidate(
                candidate_id=tx.transcript_id,
                intron_chain=chain,
                three_prime_pos=three_prime,
                sequence=seq,
                source="gtf",
                supporting_read_ids=set(),
                chrom=interval.chrom,
                strand=tx.strand,
                start=tx.start,
                end=tx.end,
            )
        )

    # 5. Match novel groups vs GTF: same intron chain + 3' within threshold → merge
    candidates: List[TranscriptCandidate] = list(gtf_candidates)
    matched_group_keys: Set = set()

    for group_key, reads in groups.items():
        consensus_3prime, chain = group_key
        read_ids_in_group = {rd["query_name"] for rd in reads if "query_name" in rd}

        # Try to match with a GTF candidate
        matched = False
        for gtf_cand in gtf_candidates:
            if _intron_chains_match(chain, gtf_cand.intron_chain) and _three_prime_within_threshold(
                consensus_3prime, gtf_cand.three_prime_pos, threshold
            ):
                gtf_cand.supporting_read_ids.update(read_ids_in_group)
                matched = True
                matched_group_keys.add(group_key)
                break

        if not matched:
            # Novel candidate
            rep = pick_representative_read(reads)
            ref_seq = rep.get("reference_sequence", "")
            ref_start = rep.get("reference_start", interval.start)
            ref_end = rep.get("reference_end", interval.end)

            # 5' takes longest: use representative's span
            candidates.append(
                TranscriptCandidate(
                    candidate_id=_generate_novel_id(),
                    intron_chain=chain,
                    three_prime_pos=consensus_3prime,
                    sequence=ref_seq if ref_seq else "",
                    source="novel",
                    supporting_read_ids=read_ids_in_group,
                    chrom=interval.chrom,
                    strand=strand,
                    start=ref_start,
                    end=ref_end,
                )
            )

    # 6. Collapse: same 3' consensus + same intron chain → one candidate, 5' takes longest
    collapsed = _collapse_candidates(candidates, threshold)

    logger.info(
        "Interval %s: %d candidates (%d GTF, %d novel) from %d reads",
        interval.region_string,
        len(collapsed),
        sum(1 for c in collapsed if c.source == "gtf"),
        sum(1 for c in collapsed if c.source == "novel"),
        len(all_read_ids),
    )

    return CandidateSet(
        interval=interval, candidates=collapsed, read_ids=all_read_ids
    )


def _collapse_candidates(
    candidates: List[TranscriptCandidate], threshold: int
) -> List[TranscriptCandidate]:
    """Collapse candidates with same 3' consensus and intron chain.

    When duplicates are found, keep the one with the longest span (5' extension).
    Merge supporting read IDs.
    """
    if not candidates:
        return []

    # Group by (intron_chain, approximate 3' position)
    groups: dict = {}
    for cand in candidates:
        # Find existing group with matching chain and close 3'
        merged = False
        for key in list(groups.keys()):
            existing_3prime, existing_chain = key
            if (
                existing_chain == cand.intron_chain
                and abs(existing_3prime - cand.three_prime_pos) <= threshold
            ):
                existing = groups[key]
                # Keep the longer one (larger span)
                if cand.length > existing.length:
                    cand.supporting_read_ids.update(existing.supporting_read_ids)
                    groups[key] = cand
                else:
                    existing.supporting_read_ids.update(cand.supporting_read_ids)
                merged = True
                break

        if not merged:
            groups[(cand.three_prime_pos, cand.intron_chain)] = cand

    return list(groups.values())
