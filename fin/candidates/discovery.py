"""Candidate transcript discovery per genomic interval."""

from __future__ import annotations

import logging
from typing import List, Set

from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.candidates.intron_chains import (
    group_reads_by_three_prime_and_intron_chain,
    gtf_transcript_to_intron_chain,
    pick_representative_read,
)
from fin.io.interval_manager import (
    GenomicInterval,
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


_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _exons_from_chain(start: int, end: int, chain: IntronChain) -> List[tuple]:
    """Derive exon (start, end) blocks from an intron chain spanning [start, end)."""
    introns = chain.introns
    if not introns:
        return [(start, end)]
    exons = [(start, introns[0][0])]
    for i in range(len(introns) - 1):
        exons.append((introns[i][1], introns[i + 1][0]))
    exons.append((introns[-1][1], end))
    # Drop any zero-length blocks defensively
    return [(s, e) for s, e in exons if e > s]


def _build_spliced_sequence(
    genome_fasta: str,
    start: int,
    end: int,
    chain: IntronChain,
    strand: str,
) -> str:
    """Stitch a spliced transcript sequence from a genome chromosome string.

    Args:
        genome_fasta: Full chromosome sequence string.
        start: 0-based transcript start (left edge of first exon).
        end: 0-based exclusive transcript end (right edge of last exon).
        chain: Intron chain in genomic coordinates.
        strand: '+' or '-'. Reverse complement is applied for '-'.

    Returns:
        Spliced cDNA sequence, or empty string when bounds are unusable.
    """
    if not genome_fasta or end <= start:
        return ""
    chrom_len = len(genome_fasta)
    # Clip to chromosome bounds.
    s = max(0, min(start, chrom_len))
    e = max(0, min(end, chrom_len))
    if e <= s:
        return ""
    parts = []
    for ex_s, ex_e in _exons_from_chain(s, e, chain):
        ex_s_c = max(s, min(ex_s, chrom_len))
        ex_e_c = max(s, min(ex_e, chrom_len))
        if ex_e_c > ex_s_c:
            parts.append(genome_fasta[ex_s_c:ex_e_c])
    seq = "".join(parts)
    if strand == "-":
        seq = _reverse_complement(seq)
    return seq


def _three_prime_within_threshold(pos_a: int, pos_b: int, threshold: int) -> bool:
    return abs(pos_a - pos_b) <= threshold


def discover_gtf_only(
    interval: GenomicInterval,
    bam_path: str,
    gtf_reader,
    genome_fasta: str,
) -> CandidateSet:
    """Discover candidates using only GTF transcripts (no novel discovery).

    Args:
        interval: Genomic interval to process.
        bam_path: Path to BAM file.
        gtf_reader: Opened GTFReader instance (already parsed).
        genome_fasta: Full chromosome sequence string for the interval's chromosome.

    Returns:
        CandidateSet with GTF-only candidates and all read IDs from the interval.
    """
    from fin.io.io_bam import BamReader

    # Fetch all read IDs from BAM for the interval
    all_read_ids: Set[str] = set()
    with BamReader(bam_path) as bam:
        for alignment in bam.fetch(
            reference=interval.chrom, start=interval.start, end=interval.end
        ):
            rd = bam.alignment_to_dict(alignment)
            if not rd or not rd.get("is_mapped"):
                continue
            read_id = rd.get("query_name")
            if read_id:
                all_read_ids.add(read_id)

    # Get GTF transcripts in region
    gtf_transcripts = []
    if gtf_reader is not None:
        gtf_transcripts = gtf_reader.get_transcripts_in_region(
            interval.chrom, interval.start, interval.end
        )

    # Build TranscriptCandidate per GTF transcript
    gtf_candidates = []
    for tx in gtf_transcripts:
        chain = gtf_transcript_to_intron_chain(tx)
        tx.sort_features()
        three_prime = tx.end if tx.strand == "+" else tx.start
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

    logger.info(
        "Interval %s: %d GTF candidates, %d reads (gtf-only mode)",
        interval.region_string,
        len(gtf_candidates),
        len(all_read_ids),
    )

    return CandidateSet(
        interval=interval, candidates=gtf_candidates, read_ids=all_read_ids
    )


def discover_candidates(
    interval: GenomicInterval,
    bam_path: str,
    gtf_reader,
    genome_fasta: str,
    threshold: int = 24,
    min_novel_reads: int = 1,
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
            # Novel candidate: 5' takes longest read in the group; the spliced
            # cDNA sequence is stitched from genome FASTA + the intron chain.
            rep = pick_representative_read(reads)
            ref_start = rep.get("reference_start", interval.start)
            ref_end = rep.get("reference_end", interval.end)

            spliced_seq = _build_spliced_sequence(
                genome_fasta, ref_start, ref_end, chain, strand
            )
            if not spliced_seq:
                logger.warning(
                    "Skipping novel candidate at %s:%d-%d (%s): empty spliced "
                    "sequence from genome FASTA (chrom_len=%d).",
                    interval.chrom,
                    ref_start,
                    ref_end,
                    strand,
                    len(genome_fasta) if genome_fasta else 0,
                )
                continue

            candidates.append(
                TranscriptCandidate(
                    candidate_id=_generate_novel_id(),
                    intron_chain=chain,
                    three_prime_pos=consensus_3prime,
                    sequence=spliced_seq,
                    source="novel",
                    supporting_read_ids=read_ids_in_group,
                    chrom=interval.chrom,
                    strand=strand,
                    start=ref_start,
                    end=ref_end,
                )
            )

    # 6. Collapse: novel candidates with same intron chain + nearby 3' end
    collapsed = _collapse_candidates(candidates, threshold)

    # 7. (A2) Drop novel candidates below min_novel_reads threshold. GTF and
    # fusion candidates are not affected by this filter.
    if min_novel_reads > 1:
        before = len(collapsed)
        collapsed = [
            c
            for c in collapsed
            if c.source != "novel" or len(c.supporting_read_ids) >= min_novel_reads
        ]
        dropped = before - len(collapsed)
        if dropped:
            logger.info(
                "Interval %s: dropped %d novel candidates with <%d supporting reads",
                interval.region_string,
                dropped,
                min_novel_reads,
            )

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


def merge_fusion_candidates(
    candidate_set: CandidateSet,
    fusion_candidates: List[TranscriptCandidate],
) -> CandidateSet:
    """Append fusion candidates to an existing CandidateSet, deduplicating by ID.

    Args:
        candidate_set: Existing set of candidates (GTF/novel). Mutated in place.
        fusion_candidates: Fusion candidates to merge in.

    Returns:
        The same CandidateSet instance with fusion candidates appended.
    """
    existing_ids: Set[str] = {c.candidate_id for c in candidate_set.candidates}
    for fc in fusion_candidates:
        if fc.candidate_id not in existing_ids:
            candidate_set.candidates.append(fc)
            existing_ids.add(fc.candidate_id)

    all_ids = [c.candidate_id for c in candidate_set.candidates]
    if len(all_ids) != len(set(all_ids)):
        raise ValueError("Duplicate candidate_ids found after merge")

    return candidate_set


def _collapse_candidates(
    candidates: List[TranscriptCandidate], threshold: int
) -> List[TranscriptCandidate]:
    """Collapse novel candidates with same intron chain and nearby 3' end.

    Two-pass algorithm (A1):
      Pass 1: bucket NOVEL candidates by exact intron chain.
      Pass 2: within each chain-bucket, greedily merge candidates whose 3'
              ends fall within ``threshold`` of an existing representative.
    The representative is always the longest-span candidate (5' extension);
    supporting read IDs are unioned.

    GTF candidates are passed through unchanged (annotation-level duplicates
    are preserved so ground-truth recall metrics stay interpretable).
    Fusion candidates are also passed through unchanged.
    """
    if not candidates:
        return []

    gtf_passthrough: List[TranscriptCandidate] = []
    fusion_candidates: List[TranscriptCandidate] = []
    novel_by_chain: dict = {}

    for cand in candidates:
        if cand.source == "fusion":
            fusion_candidates.append(cand)
            continue
        if cand.source == "gtf":
            gtf_passthrough.append(cand)
            continue
        # novel: bucket by intron chain
        novel_by_chain.setdefault(cand.intron_chain, []).append(cand)

    collapsed_novel: List[TranscriptCandidate] = []
    for chain, bucket in novel_by_chain.items():
        # Within a chain bucket: greedy 3'-window merge. Sort by length desc
        # so the longest-span candidate naturally becomes the representative
        # of its 3' cluster (avoiding "first-wins" arbitrariness).
        bucket.sort(key=lambda c: c.length, reverse=True)
        reps: List[TranscriptCandidate] = []
        for cand in bucket:
            merged = False
            for rep in reps:
                if abs(rep.three_prime_pos - cand.three_prime_pos) <= threshold:
                    # rep is already at least as long (sort order). Only union
                    # the supporting reads; do NOT extend rep.start/rep.end to
                    # the union, because rep.sequence was built from rep's
                    # original genomic span and downstream scoring/output
                    # require coords and sequence to stay consistent (#P2-3).
                    rep.supporting_read_ids.update(cand.supporting_read_ids)
                    merged = True
                    break
            if not merged:
                reps.append(cand)
        collapsed_novel.extend(reps)

    return gtf_passthrough + collapsed_novel + fusion_candidates
