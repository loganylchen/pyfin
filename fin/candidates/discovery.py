"""Candidate transcript discovery per genomic interval."""

from __future__ import annotations

import logging
from typing import List, Set

from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.candidates.canonical import parse_motifs
from fin.candidates.intron_chains import (
    expand_canonical_chain_alternatives_v2,
    extract_intron_chain,
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


def _spatial_read_clusters(read_ids, spans):
    """Group chain-less (single-exon) reads into non-overlapping genomic loci via
    single-linkage over span overlap. Without this, every mono read in an interval
    pools into ONE candidate whose union span bridges distinct single-exon genes.
    Reads absent from ``spans`` are dropped. Returns a list of read-id lists."""
    items = sorted((spans[r][0], spans[r][1], r) for r in read_ids if r in spans)
    clusters = []
    cur = []
    cur_end = None
    for s, e, r in items:
        if cur and s > cur_end:                 # no overlap with the current locus
            clusters.append([x[2] for x in cur])
            cur = []
            cur_end = None
        cur.append((s, e, r))
        cur_end = e if cur_end is None else max(cur_end, e)
    if cur:
        clusters.append([x[2] for x in cur])
    return clusters


def _fold_mono_into_members(clusters, mono_ids, spans):
    """Fold each mono (single-exon) read wholly inside a multi-exon member's exon INTO that
    member, for the families path when no post-EM mono resolution runs (fold_monoexon semantics).

    Mirrors ``cluster_read_chains``'s generation-time mono fold: a mono read whose span lies
    wholly within one exon of a multi member (overlaps no intron) folds into the deterministically
    chosen best container (max by live read count, then #introns, then coords); reads in an intron
    or not exon-contained stay as mono. Multi footprints are computed ONCE from each member's own
    reads (pre-fold). Mutates member ``read_ids``; returns the remaining uncontained mono ids.
    """
    from fin.candidates.chain_cluster import _monoexon_in_exon
    containers = []                       # (member, lo, hi, chain_introns)
    for cl in clusters:
        for m in cl.members:
            ch = tuple(m.chain.introns)
            if not ch:
                continue
            pts = [spans[r] for r in m.read_ids if r in spans]
            if not pts:
                continue
            containers.append((m, min(p[0] for p in pts), max(p[1] for p in pts), ch))
    if not containers:
        return set(mono_ids)
    remaining = set()
    for rid in sorted(mono_ids):          # sorted -> deterministic under live-count tie-breaks
        sp = spans.get(rid)
        if sp is None:
            remaining.add(rid)
            continue
        s, e = sp
        hits = [(m, ch) for (m, lo, hi, ch) in containers
                if _monoexon_in_exon(s, e, ch, lo, hi)]
        if not hits:
            remaining.add(rid)
            continue
        best = max(hits, key=lambda t: (len(t[0].read_ids), len(t[1]), t[1]))
        best[0].read_ids.add(rid)
    return remaining


def _best_overlap_gtf(mono_gtf, lo: int, hi: int):
    """Pick the single-exon GTF candidate whose span most overlaps ``[lo, hi)``.

    Mono (empty-chain) candidates must NOT be matched to GTF by chain equality -- every
    single-exon transcript shares the empty chain, so a dict keyed by chain collapses
    them and routes all single-exon reads to whichever mono GTF was inserted first.
    Match by genomic overlap instead; returns None when nothing overlaps."""
    best = None
    best_ov = 0
    for gc, gs, ge in mono_gtf:
        ov = min(hi, ge) - max(lo, gs)
        if ov > best_ov:
            best_ov = ov
            best = gc
    return best


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


def _graph_groups(non_fusion_reads, strand, tol, min_edge_reads,
                  tss_brake=True, tss_tol=20, tss_min_reads=3, tss_frac=0.4):
    """Build novel-candidate read groups via intron-graph assembly: extend each
    read's chain to a maximal read-supported chain (:func:`assemble_chains`),
    group reads by that extended chain, and key each group by (TES, IntronChain)
    for the shared candidate-creation path. TES = the group's 3' terminus (max
    reference_end on '+', min reference_start on '-'). Truncated reads on the same
    maximal path land in ONE group -> one full-length candidate.

    ``tss_brake`` passes each read's 5'-end position to the assembler so genuine
    short isoforms sitting at a real TSS peak are not extended (merged) away."""
    from fin.candidates.junction_graph import assemble_chains

    reads_kept = []
    chains = []
    five_prime_pos = []
    for rd in non_fusion_reads:
        cigar = rd.get("cigartuples")
        ref_start = rd.get("reference_start")
        if cigar is None or ref_start is None:
            continue
        reads_kept.append(rd)
        chains.append(tuple(extract_intron_chain(cigar, ref_start).introns))
        # read 5'-end = the truncation-prone terminus (strand-aware)
        five_prime_pos.append(rd.get("reference_end") if strand == "-" else ref_start)
    if not chains:
        return {}
    _, extended = assemble_chains(
        chains, tol, min_edge_reads,
        five_prime_pos=(five_prime_pos if tss_brake else None),
        strand=strand, tss_tol=tss_tol,
        min_tss_reads=tss_min_reads, tss_frac=tss_frac,
    )

    by_chain: dict = {}
    for i, rd in enumerate(reads_kept):
        by_chain.setdefault(extended[i], []).append(rd)

    groups: dict = {}
    for ext_chain, reads in by_chain.items():
        if strand == "-":
            tes = min((rd.get("reference_start", 0) for rd in reads), default=0)
        else:
            tes = max((rd.get("reference_end", 0) for rd in reads), default=0)
        groups[(tes, IntronChain(introns=ext_chain))] = reads
    return groups


def _chain_cluster_candidates(
    interval, non_fusion_reads, all_read_ids, read_sequences, gtf_reader,
    genome_fasta, strand, min_novel_reads, *, wobble_bp, cassette_max_exon_bp,
    fold_monoexon_contained=False, fold_span_guard=False, clustering="read_chains",
) -> CandidateSet:
    """NEW generation: build candidates via intron-chain clustering.

    Group reads by intron chain (3' ignored), fold EXACT sub-chains, cluster
    wobble/cassette/containment (members kept); each surviving member chain is a NOVEL
    candidate whose span is the UNION of its reads' extents and whose sequence is
    stitched from the genome. A novel chain that EXACTLY matches a GTF transcript's
    chain merges its reads into that GTF candidate instead. No canonical search.

    ``clustering`` selects the clustering primitive (both yield the same downstream
    ``ChainCluster`` contract, so the emission below is shared):
      - ``"families"`` (production default): the clustering redesign -- :func:`cluster_families`
        (grouping only) then :func:`collapse` per family (exact-subchain fold, ``fold_span_guard``
        honoured). Mono reads go to a holding bucket: when the effective ``fold_monoexon_contained``
        flag is on they are generation-folded into containing multi members
        (:func:`_fold_mono_into_members`); otherwise the bucket is emitted standalone -- for post-EM
        ``mono_resolve_post_em`` resolution when that is enabled (the m2_em default), or as final
        standalone mono candidates when it is not. Both fold flags apply here.
      - ``"read_chains"``: the legacy :func:`cluster_read_chains` (fold + cluster inline).
    """
    read_chains = []
    spans: dict = {}
    for rd in non_fusion_reads:
        cigar = rd.get("cigartuples")
        ref_start = rd.get("reference_start")
        if cigar is None or ref_start is None:
            continue
        read_chains.append((rd, extract_intron_chain(cigar, ref_start)))
        q = rd.get("query_name")
        if q is not None:
            spans[q] = (ref_start, rd.get("reference_end", ref_start))

    if clustering == "families":
        from fin.candidates.chain_cluster import (
            ChainCluster, GenCandidate, cluster_families, collapse)
        from fin.candidates.dataclasses import IntronChain
        result = cluster_families(
            read_chains, wobble_bp=wobble_bp, cassette_max_exon_bp=cassette_max_exon_bp)
        clusters = [collapse(fam, span_guard=fold_span_guard, read_spans=spans)
                    for fam in result.families]
        mono_ids = set(result.mono_reads)
        # When the effective fold_monoexon_contained flag is on (the runner keeps it on unless a
        # post-EM mono_resolve will run), fold contained mono reads into a multi member at
        # generation to suppress single-exon truncation FPs; otherwise the bucket is emitted
        # standalone -- deferred to post-EM mono_resolve when enabled, else final mono candidates.
        if mono_ids and fold_monoexon_contained:
            mono_ids = _fold_mono_into_members(clusters, mono_ids, spans)
        if mono_ids:                     # remaining mono holding bucket -> one empty-chain cluster
            clusters.append(ChainCluster(
                members=[GenCandidate(IntronChain(introns=()), set(mono_ids))],
                read_ids=set(mono_ids)))
    else:
        from fin.candidates.chain_cluster import cluster_read_chains
        clusters = cluster_read_chains(
            read_chains, wobble_bp=wobble_bp, cassette_max_exon_bp=cassette_max_exon_bp,
            fold_monoexon_contained=fold_monoexon_contained,
            fold_span_guard=fold_span_guard)

    # GTF candidates (identical to the legacy path).
    gtf_candidates: List[TranscriptCandidate] = []
    gtf_by_chain: dict = {}
    if gtf_reader is not None:
        for tx in gtf_reader.get_transcripts_in_region(
                interval.chrom, interval.start, interval.end):
            # Interval is strand-separated (reads are strand-filtered into it), so an
            # opposite-strand annotation must NOT enter and compete for reads / absorb
            # same-chain support. (No effect when the interval is unstranded or there is
            # no GTF, e.g. p00.)
            if interval.strand is not None and tx.strand != interval.strand:
                continue
            chain = gtf_transcript_to_intron_chain(tx)
            tx.sort_features()
            three_prime = tx.end if tx.strand == "+" else tx.start
            gc = TranscriptCandidate(
                candidate_id=tx.transcript_id, intron_chain=chain,
                three_prime_pos=three_prime, sequence=tx.get_spliced_sequence(genome_fasta),
                source="gtf", supporting_read_ids=set(), chrom=interval.chrom,
                strand=tx.strand, start=tx.start, end=tx.end)
            gtf_candidates.append(gc)
            gtf_by_chain.setdefault(chain, gc)   # first GTF wins on an exact-chain tie

    # Single-exon (empty-chain) GTF transcripts, matched to mono members by genomic
    # overlap (NOT by empty-chain equality, which would collapse them all to one).
    mono_gtf = [(gc, gc.start, gc.end)
                for gc in gtf_candidates if not gc.intron_chain.introns]

    # Record each ChainCluster's members by candidate_id (NOT list index): downstream
    # gates (canonical_gate, junction-dominance, fusion) drop/reorder candidate_set
    # candidates after discovery, so an id survives those mutations where an index
    # would go stale. A member with no candidate (empty span/seq) contributes nothing.
    candidates: List[TranscriptCandidate] = list(gtf_candidates)
    cluster_id_lists: List[List[str]] = []
    # container candidate_id -> its folded SHADOW sub-chains (for 5'-TSS recovery).
    # Each shadow: (intron_chain_tuple, sorted_read_id_tuple). GTF ids can recur
    # across clusters (same annotation matched twice) -> merge the shadow lists.
    shadow_map: Dict[str, List[Tuple[Tuple[Tuple[int, int], ...], Tuple[str, ...]]]] = {}
    def _emit(chain, member_reads, start, end, gc):
        """Merge ``member_reads`` into GTF ``gc`` (if any) else emit a novel candidate.
        Returns the candidate_id, or None when the spliced sequence cannot be built."""
        if gc is not None:                        # GTF match -> merge
            gc.supporting_read_ids.update(member_reads)
            return gc.candidate_id
        seq = _build_spliced_sequence(genome_fasta, start, end, chain, strand)
        if not seq:
            return None
        cid = _generate_novel_id()
        candidates.append(TranscriptCandidate(
            candidate_id=cid, intron_chain=chain,
            three_prime_pos=(end if strand == "+" else start), sequence=seq,
            source="novel", supporting_read_ids=set(member_reads),
            chrom=interval.chrom, strand=strand, start=start, end=end))
        return cid

    for cl in clusters:
        ids: List[str] = []
        for m in cl.members:
            if m.chain.introns:                   # multi-exon: exact-chain GTF match
                starts = [spans[r][0] for r in m.read_ids if r in spans]
                ends = [spans[r][1] for r in m.read_ids if r in spans]
                if not starts:
                    continue
                cand_id = _emit(m.chain, m.read_ids, min(starts), max(ends),
                                gtf_by_chain.get(m.chain))
                if cand_id is None:
                    continue
                ids.append(cand_id)
                # folded exact-sub-chain shadows keyed by the candidate it became
                # (5'-TSS recovery). Mono members carry none.
                member_shadows = [
                    (tuple(f.chain.introns), tuple(sorted(f.read_ids)))
                    for f in m.folded if f.chain.introns
                ]
                if member_shadows:
                    shadow_map.setdefault(cand_id, []).extend(member_shadows)
            else:
                # mono (single-exon): split the pooled chain-less reads into per-locus
                # groups so distinct single-exon genes are not bridged, then overlap-match
                # each locus to a single-exon GTF (else emit a per-locus novel candidate).
                for group in _spatial_read_clusters(m.read_ids, spans):
                    gs = [spans[r][0] for r in group]
                    ge = [spans[r][1] for r in group]
                    lo, hi = min(gs), max(ge)
                    cand_id = _emit(m.chain, set(group), lo, hi,
                                    _best_overlap_gtf(mono_gtf, lo, hi))
                    if cand_id is not None:
                        ids.append(cand_id)
        cluster_id_lists.append(ids)

    # Drop low-support novels; prune cluster lists to surviving candidate ids.
    if min_novel_reads > 1:
        candidates = [c for c in candidates
                      if c.source != "novel"
                      or len(c.supporting_read_ids) >= min_novel_reads]
        keep_ids = {c.candidate_id for c in candidates}
        cluster_id_lists = [
            [cid for cid in ids if cid in keep_ids] for ids in cluster_id_lists
        ]
        # Prune shadow entries whose container candidate was dropped.
        shadow_map = {cid: sh for cid, sh in shadow_map.items() if cid in keep_ids}

    logger.info(
        "Interval %s: %d candidates (%d GTF, %d novel) from %d reads [chain-cluster]",
        interval.region_string, len(candidates),
        sum(1 for c in candidates if c.source == "gtf"),
        sum(1 for c in candidates if c.source == "novel"), len(all_read_ids))

    return CandidateSet(
        interval=interval, candidates=candidates, read_ids=all_read_ids,
        read_sequences=read_sequences, clusters=cluster_id_lists,
        shadows=shadow_map, read_spans=spans)


def discover_candidates(
    interval: GenomicInterval,
    bam_path: str,
    gtf_reader,
    genome_fasta: str,
    threshold: int = 24,
    min_novel_reads: int = 1,
    canonical_search_bp: int = 0,
    max_chains_per_read: int = 16,
    canonical_motifs: tuple = ("GT-AG", "GC-AG", "AT-AC"),
    denovo_wobble_tol: int = 0,
    denovo_wobble_shadow_ratio: float = 1.0,
    denovo_graph: bool = False,
    denovo_graph_tol: int = 6,
    denovo_graph_min_edge_reads: int = 2,
    denovo_graph_tss_brake: bool = True,
    denovo_graph_tss_tol: int = 20,
    denovo_graph_tss_min_reads: int = 3,
    denovo_graph_tss_frac: float = 0.4,
    chain_cluster: bool = False,
    chain_cluster_wobble_bp: int = 6,
    chain_cluster_cassette_max_exon_bp: int = 70,
    chain_cluster_fold_monoexon: bool = False,
    chain_cluster_fold_span_guard: bool = False,
    clustering: str = "read_chains",
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
        canonical_search_bp: If > 0, for every read scan ±N bp around each
            junction's CIGAR position for canonical GT-AG (donor/acceptor)
            motifs and emit the resulting chain alternatives in addition to
            the original chain. 0 disables.
        max_chains_per_read: Hard cap on the per-read alternative chain count
            during canonical search. Reads whose Cartesian-product expansion
            exceeds the cap fall back to their original chain.

    Returns:
        CandidateSet with discovered candidates.
    """
    from fin.io.io_bam import BamReader

    # 1. Fetch reads for interval, skip fusion reads
    all_read_ids: Set[str] = set()
    non_fusion_reads = []
    read_sequences: dict = {}

    with BamReader(bam_path) as bam:
        for alignment in bam.fetch(
            reference=interval.chrom, start=interval.start, end=interval.end
        ):
            rd = bam.alignment_to_dict(alignment)
            if not rd or not rd.get("is_mapped"):
                continue
            # Honor strand-separated intervals. bam.fetch returns reads on BOTH
            # strands overlapping the region, but each read belongs to its own
            # mapped strand's interval. Without this guard an antisense interval
            # that spatially contains a sense interval swallows the sense reads,
            # corrupting strand assignment (a candidate emitted on the wrong
            # strand can never match) and read-to-candidate quantification.
            if interval.strand is not None:
                read_strand = "+" if rd.get("is_forward") else "-"
                if read_strand != interval.strand:
                    continue
            read_id = rd.get("query_name")
            if read_id:
                all_read_ids.add(read_id)
                # Capture primary-alignment sequence only — scope is this
                # interval's fetch, satisfying the per-interval read-set
                # constraint. Used by R1 mappy-argmax bypass.
                if not rd.get("is_secondary") and not rd.get("is_supplementary"):
                    seq = rd.get("query_sequence")
                    if seq:
                        read_sequences.setdefault(read_id, seq)
            if is_fusion_read(rd):
                continue
            non_fusion_reads.append(rd)

    if not non_fusion_reads:
        return CandidateSet(
            interval=interval,
            candidates=[],
            read_ids=all_read_ids,
            read_sequences=read_sequences,
        )

    # Determine strand from reads
    strand = interval.strand or "+"

    # NEW generation-side clustering (production default). Group reads by intron chain
    # (3' ignored), collapse EXACT sub-chains, cluster wobble/cassette/containment
    # keeping members; NO canonical-search shadow generation. Bypasses the legacy
    # 3'+exact-chain grouping and canonical expansion entirely.
    if chain_cluster:
        return _chain_cluster_candidates(
            interval, non_fusion_reads, all_read_ids, read_sequences,
            gtf_reader, genome_fasta, strand, min_novel_reads,
            wobble_bp=chain_cluster_wobble_bp,
            cassette_max_exon_bp=chain_cluster_cassette_max_exon_bp,
            fold_monoexon_contained=chain_cluster_fold_monoexon,
            fold_span_guard=chain_cluster_fold_span_guard,
            clustering=clustering,
        )

    # 2. Group reads by 3' + intron chain → novel candidate groups.
    # Optional pre-processing:
    #   - canonical_search_bp>0: for each junction, scan the genome ±N bp for
    #     canonical GT-AG motifs and emit the resulting chain alternatives
    #     alongside the original (recovering ambiguous-GT-AG strict misses
    #     where multiple valid splice sites exist within a few bp).
    multi_chain_override = None

    if canonical_search_bp > 0:
        # Pre-extract CIGAR-derived chains once.
        read_chains = []
        for rd in non_fusion_reads:
            cigar = rd.get("cigartuples")
            ref_start = rd.get("reference_start")
            if cigar is None or ref_start is None:
                continue
            read_chains.append((rd, extract_intron_chain(cigar, ref_start)))

        multi_chain_override = expand_canonical_chain_alternatives_v2(
            read_chains, genome_fasta, strand,
            motif_set=parse_motifs(canonical_motifs),
            search_bp=canonical_search_bp,
            max_chains_per_read=max_chains_per_read,
        )

    if denovo_graph:
        # Intron-graph assembly: group reads by their EXTENDED (assembled) chain
        # so truncated partials become one full-length candidate. Takes precedence
        # over canonical-search expansion.
        groups = _graph_groups(
            non_fusion_reads, strand, denovo_graph_tol, denovo_graph_min_edge_reads,
            tss_brake=denovo_graph_tss_brake, tss_tol=denovo_graph_tss_tol,
            tss_min_reads=denovo_graph_tss_min_reads, tss_frac=denovo_graph_tss_frac,
        )
    else:
        groups = group_reads_by_three_prime_and_intron_chain(
            non_fusion_reads, strand, threshold=threshold,
            multi_chain_override=multi_chain_override,
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
        # Strand-separated interval: skip opposite-strand annotation (mirror of the
        # chain-cluster path). No effect when the interval is unstranded / no GTF.
        if interval.strand is not None and tx.strand != interval.strand:
            continue
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
            if denovo_graph:
                # Assembled candidate spans regions no single read covers -> use
                # the UNION of the group's reads' extents so the genome-stitched
                # spliced sequence covers the full extended intron chain.
                starts = [rd.get("reference_start") for rd in reads
                          if rd.get("reference_start") is not None]
                ends = [rd.get("reference_end") for rd in reads
                        if rd.get("reference_end") is not None]
                ref_start = min(starts) if starts else interval.start
                ref_end = max(ends) if ends else interval.end
            else:
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
    collapsed = _collapse_candidates(
        candidates, threshold, wobble_tol=denovo_wobble_tol,
        wobble_shadow_ratio=denovo_wobble_shadow_ratio)

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
        interval=interval,
        candidates=collapsed,
        read_ids=all_read_ids,
        read_sequences=read_sequences,
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


def _wobble_chain_match(a, b, tol: int) -> bool:
    """True iff two intron chains have the SAME number of introns and every
    corresponding (donor, acceptor) is within ``tol`` bp — same structure up to
    alignment wobble."""
    if a.num_introns != b.num_introns:
        return False
    for (s1, e1), (s2, e2) in zip(a.introns, b.introns):
        if abs(s1 - s2) > tol or abs(e1 - e2) > tol:
            return False
    return True


def _wobble_collapse_novel(
    novel: List[TranscriptCandidate], threshold: int, tol: int,
    shadow_ratio: float = 1.0,
) -> List[TranscriptCandidate]:
    """Greedy consensus merge of novel candidates that differ only by junction
    wobble (each junction within ``tol`` bp) AND share a 3' end within
    ``threshold``. The HIGHEST-read-support candidate is the consensus
    representative (the reads' own mode); wobbled shadows union their reads into
    it and vanish. Deterministic order: (-support, -length, start, 3', end, chain).

    ``shadow_ratio`` guards against merging genuine CLOSE isoforms (real
    alternative donor/acceptor sites a few bp apart, e.g. NAGNAG): a candidate is
    absorbed into a rep only if it is a true SHADOW — ``len(cand.reads) <=
    shadow_ratio * len(rep.reads)`` (rep always has >= support, given the sort).
    A minimap wobble shadow has few reads; a real close isoform has comparable
    support and stays separate. isoquant's bulge-removal uses the ½ analogue.
    ``shadow_ratio >= 1.0`` merges any wobble-match (no isoform protection).

    NOT annotation-snap: the consensus coordinate comes from the most-supported
    READS, never a reference. Merges shadows BEFORE they exist as separate
    candidates (reads never split-then-dropped — avoids the read-reassignment
    trap that made the post-hoc dominance filter backfire).
    """
    # Total, data-derived order (no reliance on input order or the id counter):
    # highest read support first (= consensus), then longest, then genomic
    # coords + chain. This makes both the representative choice AND the greedy
    # clustering deterministic regardless of the candidates' input order.
    # Total order (candidate_id last makes it strict, so representative choice
    # and clustering are fully input-order-independent).
    order = sorted(
        novel,
        key=lambda c: (-len(c.supporting_read_ids), -c.length, c.start,
                       c.three_prime_pos, c.end, tuple(c.intron_chain.introns),
                       c.candidate_id),
    )
    reps: List[TranscriptCandidate] = []
    for cand in order:
        merged = False
        for rep in reps:
            if (abs(rep.three_prime_pos - cand.three_prime_pos) > threshold
                    or not _wobble_chain_match(rep.intron_chain, cand.intron_chain, tol)):
                continue
            # EXACT-chain matches always merge (superset of the disabled exact
            # path). The shadow_ratio guard applies ONLY to inexact (wobbled)
            # matches — it protects genuine close isoforms, not exact duplicates.
            exact = rep.intron_chain == cand.intron_chain
            if exact or (len(cand.supporting_read_ids)
                         <= shadow_ratio * len(rep.supporting_read_ids)):
                rep.supporting_read_ids.update(cand.supporting_read_ids)
                merged = True
                break
        if not merged:
            reps.append(cand)
    return reps


def _exact_collapse_novel(
    novel: List[TranscriptCandidate], threshold: int
) -> List[TranscriptCandidate]:
    """Exact-chain collapse (THE default path — byte-identical when wobble is
    off). Bucket novel candidates by EXACT intron chain, then within each bucket
    greedily merge candidates whose 3' ends are within ``threshold`` of an
    existing (longest-span) representative, unioning supporting reads."""
    novel_by_chain: dict = {}
    for cand in novel:
        novel_by_chain.setdefault(cand.intron_chain, []).append(cand)
    collapsed: List[TranscriptCandidate] = []
    for chain, bucket in novel_by_chain.items():
        # Sort by length desc so the longest-span candidate is the representative
        # of its 3' cluster (avoiding "first-wins" arbitrariness).
        bucket.sort(key=lambda c: c.length, reverse=True)
        reps: List[TranscriptCandidate] = []
        for cand in bucket:
            merged = False
            for rep in reps:
                if abs(rep.three_prime_pos - cand.three_prime_pos) <= threshold:
                    # rep is already at least as long (sort order). Only union the
                    # supporting reads; do NOT extend rep.start/rep.end (rep.sequence
                    # was built from rep's own span; coords/sequence must stay
                    # consistent — #P2-3).
                    rep.supporting_read_ids.update(cand.supporting_read_ids)
                    merged = True
                    break
            if not merged:
                reps.append(cand)
        collapsed.extend(reps)
    return collapsed


def _collapse_candidates(
    candidates: List[TranscriptCandidate], threshold: int, wobble_tol: int = 0,
    wobble_shadow_ratio: float = 1.0,
) -> List[TranscriptCandidate]:
    """Collapse novel candidates with same intron chain and nearby 3' end.

    The exact-chain collapse (:func:`_exact_collapse_novel`) ALWAYS runs. When
    ``wobble_tol`` > 0, its reps are then wobble-collapsed
    (:func:`_wobble_collapse_novel`) — so wobble-on is a STRICT SUPERSET of the
    exact collapse (exact 3'-chaining fully preserved; wobble only adds inexact
    shadow absorption). ``wobble_tol == 0`` is byte-identical (exact only).

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
    novel: List[TranscriptCandidate] = []

    for cand in candidates:
        if cand.source == "fusion":
            fusion_candidates.append(cand)
            continue
        if cand.source == "gtf":
            gtf_passthrough.append(cand)
            continue
        novel.append(cand)

    # Exact collapse ALWAYS (byte-identical to the disabled path); wobble on top.
    collapsed_novel = _exact_collapse_novel(novel, threshold)
    if wobble_tol > 0:
        collapsed_novel = _wobble_collapse_novel(
            collapsed_novel, threshold, wobble_tol, wobble_shadow_ratio)

    return gtf_passthrough + collapsed_novel + fusion_candidates
