"""Candidate transcript discovery per genomic interval."""

from __future__ import annotations

import hashlib
import logging
from typing import Dict, List, Optional, Set, Tuple

from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.candidates.canonical import parse_motifs
from fin.candidates.chain_cluster import read_alignment_from_cigar
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


def _generate_novel_id(
    chrom: str,
    strand: str,
    start: int,
    end: int,
    chain: IntronChain,
) -> str:
    """Generate a stable structural ID for a novel candidate."""
    introns = ";".join(f"{donor}-{acceptor}" for donor, acceptor in chain.introns)
    payload = f"{chrom}|{strand}|{start}|{end}|{introns}".encode("ascii")
    return "novel_" + hashlib.blake2b(payload, digest_size=8).hexdigest()


def _generate_family_id(chrom: str, strand: str, chains) -> str:
    """Generate a stable ID from every distinct splice chain in one family."""
    normalized = sorted({
        tuple(chain.introns) if isinstance(chain, IntronChain) else tuple(chain)
        for chain in chains
    })
    chain_text = "/".join(
        ";".join(f"{donor}-{acceptor}" for donor, acceptor in chain)
        for chain in normalized
    )
    payload = f"{chrom}|{strand}|{chain_text}".encode("utf-8")
    return "fam_" + hashlib.blake2b(payload, digest_size=8).hexdigest()


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


COMPATIBLE_READ_END_TOL_BP = 50   # read 3' end may run past the annotated TES by this much
COMPATIBLE_READ_SLOP_BP = 10      # alignment slop allowed at the upstream model intron
COMPATIBLE_READ_INDEL_FLANK_BP = 0  # 0: no per-read indel policy (see below)


def chain_compatible_read_depths(
    gtf_models, reads, indels_by_read, *,
    end_tol_bp: int = COMPATIBLE_READ_END_TOL_BP,
    flank_bp: int = COMPATIBLE_READ_SLOP_BP,
    indel_flank_bp: int = COMPATIBLE_READ_INDEL_FLANK_BP,
) -> Dict[str, int]:
    """Distinct reads whose alignment is COMPATIBLE with an annotated chain.

    A read is compatible with an annotated multi-exon model when its own intron
    chain is EQUAL to the model's chain or a strand-aware contiguous 3' SUFFIX of
    it (the 5'-degraded dRNA fragment of a longer molecule: on '+' the chain's
    tail, on '-' its head, introns in ascending genomic order), its span does not
    run across the model intron immediately upstream of that suffix by more than
    ``flank_bp`` (an unspliced read over that intron contradicts the model, it is
    not degraded from it), and its 3' end lies within ``end_tol_bp`` of the
    annotated 3' end. The 3' condition is an OVERSHOOT limit only: a read may
    end anywhere inside the model's last exon (an upstream polyadenylation
    site is compatible with the chain) but not more than ``end_tol_bp`` past
    the annotated end. The 3' side is the read's END on '+' and its START on
    '-' -- the polyA side, the one dRNA sequences first and can be trusted. Reads are counted once per
    model by query name; a read compatible with several nested models counts
    for each of them. An exact chain match does NOT make a read a complete
    molecule of that model (a 5'-truncated read of a longer isoform can carry
    exactly the shorter annotated chain), so no exclusivity is applied: this is
    COMPATIBILITY support, not evidence that singles a model out, and callers
    must treat it as such.

    ``indel_flank_bp`` > 0 additionally requires every junction the read carries
    to be indel-clean within that flank. It is OFF by default: the read's
    junctions must already EQUAL the annotated ones, so an indel-shifted junction
    shows up as a mismatch, and on dRNA (~1 indel per 10-20 bases) requiring
    every junction of a read to be clean is a read-length filter, not an
    evidence filter -- measured on human H9 `full`, it cut the >= 3-read
    attested set from 1079 to 210 candidates at no precision gain (90% vs 95%
    truth, i.e. it removed 775 true transcripts to remove 61 false ones).

    ``gtf_models``: {intron tuple: [TranscriptCandidate,...]} for one strand-
    separated interval. ``reads``: iterable of (query_name, strand, intron tuple,
    ref_start, ref_end) for MAPQ>=1 primary alignments. Returns
    {candidate_id: depth}; models with no compatible read are absent.
    """
    by_terminal: Dict[Tuple[str, Tuple[int, int]], list] = {}
    for chain, models in gtf_models.items():
        if not chain:
            continue
        for g in models:
            strand = getattr(g, "strand", "+")
            term = chain[-1] if strand == "+" else chain[0]
            by_terminal.setdefault((strand, term), []).append((chain, g))
    if not by_terminal:
        return {}
    seen: Dict[str, set] = {}
    for q, strand, rchain, rs, re_ in reads:
        if not rchain:
            continue
        if indel_flank_bp > 0:
            ind = indels_by_read.get(q, ())
            if not all(
                not (lo - indel_flank_bp <= d <= hi + indel_flank_bp
                     or lo - indel_flank_bp <= a <= hi + indel_flank_bp)
                for d, a in rchain for lo, hi in ind
            ):
                continue
        term = rchain[-1] if strand == "+" else rchain[0]
        for chain, g in by_terminal.get((strand, term), ()):
            n = len(rchain)
            if n > len(chain):
                continue
            if strand == "+":
                if chain[len(chain) - n:] != rchain:
                    continue
                # the model intron just upstream of the suffix (5' side): the
                # read must start at or after its acceptor, allowing alignment slop
                if n < len(chain) and rs < chain[len(chain) - n - 1][1] - flank_bp:
                    continue
                if re_ > g.end + end_tol_bp:
                    continue
            else:
                if chain[:n] != rchain:
                    continue
                if n < len(chain) and re_ > chain[n][0] + flank_bp:
                    continue
                if rs < g.start - end_tol_bp:
                    continue
            seen.setdefault(g.candidate_id, set()).add(q)
    return {cid: len(qs) for cid, qs in seen.items()}


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


def junction_clean_depths(
    chain, member_read_ids, folded_read_ids, indels_by_read,
    flank_bp: int = 10,
) -> Tuple[int, ...]:
    """Per-junction indel-clean attestation depths of a complete-chain read
    set (5'->3' genomic junction order).

    Exact complete-chain reads = member reads minus folded sub-chain
    (degraded) reads: each such read's OWN raw alignment equals the
    candidate's full intron chain (traverses the 5'-most junction; the 5'
    terminus may sit anywhere inside the first exon - dRNA 5' degradation
    tolerated). A read indel-cleanly attests one junction when its alignment
    has no insertion/deletion within ``flank_bp`` of that junction's donor
    or acceptor (ESPRESSO-style PER-JUNCTION clean flank; cigar-level only -
    mismatches are not inspected). Reads absent from ``indels_by_read``
    (e.g. MAPQ-filtered at collection) never count.

    Why per-junction rather than per-read: a whole-read all-junction
    cleanliness AND was measured and rejected (genuine 16-17-junction SIRV
    chains with 7-12 exact reads had ZERO fully clean reads), while raw
    exact-count-only was also rejected (misalignment-born variant chains
    carry >=2 exact reads, but at their spurious junction EVERY exact read
    bears the indel signature: measured weakest depth 0-1 for all 17 revived
    junk chains vs >=2 for 12/13 genuine chains). Mono chains return ().
    """
    introns = tuple(chain.introns) if hasattr(chain, "introns") else tuple(chain)
    if not introns:
        return ()
    exact = [
        rid for rid in set(member_read_ids) - set(folded_read_ids)
        if rid in indels_by_read
    ]
    if not exact:
        return (0,) * len(introns)
    depths = []
    for donor, acceptor in introns:
        depth = 0
        for rid in exact:
            clean = True
            for lo, hi in indels_by_read[rid]:
                if (lo - flank_bp <= donor <= hi + flank_bp) or (
                    lo - flank_bp <= acceptor <= hi + flank_bp
                ):
                    clean = False
                    break
            if clean:
                depth += 1
        depths.append(depth)
    return tuple(depths)


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


def _read_extent(gc, contradicted, spans):
    """The read-extent sibling's span. ONLY a contradicted end is re-estimated,
    from the CONTRADICTING reads alone (``contradicted`` maps "tes"/"tss" to
    those read ids) and only outward (min/max with the annotated coordinate),
    so the span always contains the annotation's intron chain and a read that
    another annotated model explains never defines this model's end. The other
    end keeps the annotated value. On '+' the TES is ``end`` and the TSS
    ``start``; on '-' the reverse."""
    start, end = gc.start, gc.end
    tes_is_end = getattr(gc, "strand", "+") == "+"
    for which, rids in contradicted.items():
        sp = [spans[r] for r in rids if r in spans]
        if not sp:
            continue
        extend_end = (which == "tes") == tes_is_end
        if extend_end:
            end = max(end, max(e[1] for e in sp))
        else:
            start = min(start, min(e[0] for e in sp))
    return start, end


def _wobble_equivalent_reads(family_result, cluster_index, chain, wobble_bp, *, default):
    """Reads of every variant in the cluster's family that is a wobble of ``chain``
    (the chain itself included). ``default`` when the cluster is not a family
    (legacy clustering / the trailing mono bucket)."""
    if family_result is None or cluster_index >= len(family_result.families):
        return set(default)
    from fin.candidates.chain_cluster import _wobble

    fam = family_result.families[cluster_index]
    out: Set[str] = set()
    for v in fam.variants:
        if _wobble(v, chain, wobble_bp):
            out.update(fam.variant_reads.get(v, ()))
    return out or set(default)


def _chain_cluster_candidates(
    interval, non_fusion_reads, all_read_ids, read_sequences, gtf_reader,
    genome_fasta, strand, min_novel_reads, *, wobble_bp, cassette_max_exon_bp,
    fold_monoexon_contained=False, fold_span_guard=False, clustering="read_chains",
    collect_chain_evidence=True, junction_correction=None, junction_correction_apply=None,
    read_extent_contradicted=None, collect_compatible_depths=False,
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
    # Complete-chain / junction evidence is ALWAYS collected (one extra CIGAR
    # walk per MAPQ>=1 read, no alignment or signal work): per-read indel
    # positions, per-candidate indel-clean complete-chain depths, and the
    # interval-wide junction raw/clean registry. Collecting it is unconditional
    # so every run records it; whether any gate ACTS on it stays opt-in
    # (complete_chain_rescue_reads, containment/election/correction modes).
    # Nothing on the default path serializes these fields, so default outputs
    # are byte-identical.
    indels_by_read: dict = {}
    junction_counts: dict = {}     # (donor, acceptor) -> [raw, indel_clean]
    alignments: dict = {}          # query_name -> ReadAlignment (collapse span_guard)
    read_junctions: dict = {}      # query_name -> own CIGAR intron chain (evidence row)
    read_strands: dict = {}        # query_name -> mapped strand (evidence row)
    _FLANK = 10
    for rd in non_fusion_reads:
        cigar = rd.get("cigartuples")
        ref_start = rd.get("reference_start")
        if cigar is None or ref_start is None:
            continue
        chain = extract_intron_chain(cigar, ref_start)
        read_chains.append((rd, chain))
        q = rd.get("query_name")
        if q is not None:
            spans[q] = (ref_start, rd.get("reference_end", ref_start))
            alignments[q] = read_alignment_from_cigar(cigar, ref_start)
            read_junctions[q] = tuple(chain.introns)
            read_strands[q] = "+" if rd.get("is_forward", True) else "-"
            if int(rd.get("mapping_quality") or 0) >= 1:
                pos = ref_start
                ind = []
                jxs = []
                for op, ln in cigar:
                    if op in (0, 7, 8):    # M/=/X consume reference
                        pos += ln
                    elif op == 2:          # D: reference span
                        ind.append((pos, pos + ln))
                        pos += ln
                    elif op == 1:          # I: reference point
                        ind.append((pos, pos))
                    elif op == 3:          # N: junction
                        jxs.append((pos, pos + ln))
                        pos += ln
                indels_by_read[q] = tuple(ind)
                for (d, a) in jxs:
                    ent = junction_counts.setdefault((d, a), [0, 0])
                    ent[0] += 1
                    if all(
                        not (lo - _FLANK <= d <= hi + _FLANK
                             or lo - _FLANK <= a <= hi + _FLANK)
                        for lo, hi in ind
                    ):
                        ent[1] += 1

    # Step 5a ballots (additive junction correction). The callable returns
    # proposals only; the RAW read chains below are what clustering sees, so the
    # baseline candidate set is untouched. Proposals are instantiated as EXTRA
    # candidates after emission (see the end of this function).
    jc_result = junction_correction(read_chains, indels_by_read) if junction_correction else None

    # Annotation hypotheses for this (strand-separated) interval. Multi-exon
    # chains enter clustering as ordinary zero-read variant keys -- they join,
    # bridge, fold and are emitted exactly like read-derived chains; the only
    # difference is `source`. Single-exon annotation is matched to mono loci by
    # genomic overlap below. Several annotated transcripts may share one chain
    # (alternative TSS/TES): they all ride on that one key.
    gtf_models: Dict[Tuple[Tuple[int, int], ...], List[TranscriptCandidate]] = {}
    mono_gtf_candidates: List[TranscriptCandidate] = []
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
            if chain.introns:
                gtf_models.setdefault(tuple(chain.introns), []).append(gc)
            else:
                mono_gtf_candidates.append(gc)
    for models in gtf_models.values():
        models.sort(key=lambda g: (g.start, g.end, g.candidate_id))
    # Compatibility support for every annotated multi-exon model (exact or
    # strand-aware 3'-suffix reads, MAPQ>=1 only: those are the reads with an
    # indel record). Unlike the complete-chain depths (one extra CIGAR walk),
    # this is a per-read x per-model scan, so it is computed only when the
    # support-gate attestation will actually read it -- the default path pays
    # nothing.
    if gtf_models and collect_compatible_depths:
        depths = chain_compatible_read_depths(
            gtf_models,
            ((q, read_strands[q], read_junctions[q], spans[q][0], spans[q][1])
             for q in indels_by_read),
            indels_by_read)
        for models in gtf_models.values():
            for g in models:
                g.compatible_read_depth = depths.get(g.candidate_id, 0)

    family_result = None
    cluster_family_ids: List[Optional[str]] = []
    if clustering == "families":
        from fin.candidates.chain_cluster import (
            ChainCluster, GenCandidate, cluster_families, collapse)
        family_result = cluster_families(
            read_chains, wobble_bp=wobble_bp, cassette_max_exon_bp=cassette_max_exon_bp,
            gtf_variants=[(g.candidate_id, ch)
                          for ch in sorted(gtf_models) for g in gtf_models[ch]],
            # Contradiction evidence must come from EVERY read chain that
            # enters clustering. `junction_counts` is the MAPQ>=1 SUPPORT
            # registry (raw/clean depths for the election and the audits); using
            # it here would let a MAPQ=0 read shape the read-derived topology
            # while being invisible to the gate, so its junction would count as
            # "the reads are silent" when the reads in fact disagree.
            observed_junctions={j for _rd, ch in read_chains
                                for j in ch.introns})
        clusters = []
        for fam in family_result.families:
            cluster_family_ids.append(
                _generate_family_id(interval.chrom, strand, list(fam.variants)))
            clusters.append(collapse(fam, span_guard=fold_span_guard,
                                     read_alignments=alignments))
        mono_ids = set(family_result.mono_reads)
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

    # Trailing clusters (mono bucket / legacy path) carry no family id.
    cluster_family_ids.extend([None] * (len(clusters) - len(cluster_family_ids)))

    # Single-exon (empty-chain) GTF transcripts, matched to mono members by genomic
    # overlap (NOT by empty-chain equality, which would collapse them all to one).
    mono_gtf = [(gc, gc.start, gc.end) for gc in mono_gtf_candidates]

    # Record each ChainCluster's members by candidate_id (NOT list index): downstream
    # gates (canonical_gate, junction-dominance, fusion) drop/reorder candidate_set
    # candidates after discovery, so an id survives those mutations where an index
    # would go stale. A member with no candidate (empty span/seq) contributes nothing.
    # Mono annotation is always a candidate; multi-exon annotation is emitted only
    # through its surviving variant key (a folded annotated key is a dormant shadow,
    # exactly like a folded read chain).
    candidates: List[TranscriptCandidate] = list(mono_gtf_candidates)
    emitted_gtf: Set[str] = {g.candidate_id for g in mono_gtf_candidates}
    cluster_id_lists: List[List[str]] = []
    # container candidate_id -> its folded SHADOW sub-chains (for 5'-TSS recovery).
    # Each shadow: (intron_chain_tuple, sorted_read_id_tuple). GTF ids can recur
    # across clusters (same annotation matched twice) -> merge the shadow lists.
    shadow_map: Dict[str, List[Tuple[Tuple[Tuple[int, int], ...], Tuple[str, ...]]]] = {}
    extent_siblings: Dict[str, TranscriptCandidate] = {}
    annotation_ids = set(getattr(gtf_reader, "transcripts", None) or {}) if gtf_reader else set()
    def _emit(chain, member_reads, start, end, gc, family_id=None,
              junction_depths=(), extent_of=None):
        """Emit the annotation model(s) owning ``chain`` (annotated geometry, key
        reads as support) or a novel with its discovery family. Returns ids.
        ``extent_of``: emit instead the READ-EXTENT sibling of that annotated
        model -- its chain and junction trust (source "gtf"), the reads' span."""
        weakest = min(junction_depths) if junction_depths else 0
        if gc is not None:                        # annotated key -> its model(s)
            models = gc if isinstance(gc, list) else [gc]
            for g in models:
                g.supporting_read_ids.update(member_reads)
                g.raw_supporting_read_ids.update(member_reads)
                if family_id is not None:
                    g.family_id = family_id
                if weakest > g.complete_chain_weakest_indel_clean_depth or not (
                    g.complete_chain_junction_depths
                ):
                    g.complete_chain_junction_depths = tuple(junction_depths)
                    g.complete_chain_weakest_indel_clean_depth = int(weakest)
                if g.candidate_id not in emitted_gtf:
                    candidates.append(g)
                    emitted_gtf.add(g.candidate_id)
            return [g.candidate_id for g in models]
        seq = _build_spliced_sequence(genome_fasta, start, end, chain, strand)
        if not seq:
            return None
        if extent_of is not None:
            cid, source = f"{extent_of}.ext", "gtf"
            while cid in annotation_ids:          # never shadow a real annotation id
                cid += "_"
            prev = extent_siblings.get(cid)
            if prev is not None:
                # The annotation matched again (another read group / cluster):
                # ONE sibling, whose extent and support are the union of every
                # group that contradicted the annotated ends.
                prev.supporting_read_ids.update(member_reads)
                prev.raw_supporting_read_ids.update(member_reads)
                lo, hi = min(prev.start, start), max(prev.end, end)
                if (lo, hi) != (prev.start, prev.end):
                    seq2 = _build_spliced_sequence(genome_fasta, lo, hi, chain, strand)
                    if seq2:
                        prev.start, prev.end, prev.sequence = lo, hi, seq2
                        prev.three_prime_pos = hi if strand == "+" else lo
                return [cid]
        else:
            cid, source = _generate_novel_id(
                interval.chrom, strand, start, end, chain), "novel"
        candidates.append(TranscriptCandidate(
            candidate_id=cid, intron_chain=chain,
            three_prime_pos=(end if strand == "+" else start), sequence=seq,
            source=source, supporting_read_ids=set(member_reads),
            raw_supporting_read_ids=set(member_reads),
            chrom=interval.chrom, strand=strand, start=start, end=end,
            family_id=family_id, extent_of=extent_of,
            # a read-extent sibling inherits ITS annotation's chain-level
            # support; the count was span-conditioned on the annotated ends,
            # not the sibling's, so treat it as chain support only
            compatible_read_depth=(
                next((m.compatible_read_depth
                      for m in gtf_models.get(tuple(chain.introns), ())
                      if m.candidate_id == extent_of), 0)
                if extent_of is not None else 0),
            complete_chain_junction_depths=tuple(junction_depths),
            complete_chain_weakest_indel_clean_depth=(
                min(junction_depths) if junction_depths else 0
            )))
        if extent_of is not None:
            extent_siblings[cid] = candidates[-1]
            emitted_gtf.add(cid)
        return [cid]

    for cluster_index, cl in enumerate(clusters):
        ids: List[str] = []
        family_id = cluster_family_ids[cluster_index]
        for m in cl.members:
            if m.chain.introns:                   # multi-exon variant key
                models = gtf_models.get(tuple(m.chain.introns))
                starts = [spans[r][0] for r in m.read_ids if r in spans]
                ends = [spans[r][1] for r in m.read_ids if r in spans]
                if not starts and not models:
                    continue
                folded_reads: Set[str] = set()
                for f in m.folded:
                    folded_reads.update(f.read_ids)
                out_ids = _emit(
                    m.chain, m.read_ids,
                    min(starts) if starts else 0, max(ends) if ends else 0,
                    models, family_id,
                    junction_depths=junction_clean_depths(
                        m.chain, m.read_ids, folded_reads, indels_by_read),
                )
                if not out_ids:
                    continue
                ids.extend(out_ids)
                if models and read_extent_contradicted is not None:
                    # An annotated key is emitted with the ANNOTATED geometry; the
                    # reads that ride on it (its own and its folded wobble shadows)
                    # may contradict that extent -- a polyA read past the claimed
                    # TES, a read starting upstream of the claimed TSS. Then the
                    # reads' extent is a separate hypothesis: the same intron chain
                    # with the reads' span, competing in the assignment posterior.
                    # Without it a truncated annotation keeps the chain's reads
                    # only until alignment prefers a wobble sibling that happens
                    # to carry the full extent (the c_trunc failure mode). The
                    # evidence is every read of the family whose exact chain is a
                    # WOBBLE of the key (same introns within wobble_bp): those are
                    # the reads the annotated chain claims to explain, whether or
                    # not their own variant survives a later gate.
                    key_reads = _wobble_equivalent_reads(
                        family_result, cluster_index, tuple(m.chain.introns), wobble_bp,
                        default=set(m.read_ids) | folded_reads)
                    ext = [spans[r] for r in key_reads if r in spans]
                    # Every OTHER annotated model on this chain or a wobble of
                    # it can legitimately explain a read's ends (annotated
                    # APA / alternative TSS): such a read contradicts nothing.
                    from fin.candidates.chain_cluster import _wobble
                    alternatives = [
                        a for ch, ms in gtf_models.items()
                        if _wobble(ch, tuple(m.chain.introns), wobble_bp)
                        for a in ms
                    ]
                    for g in models:
                        which = ext and read_extent_contradicted(
                            g, key_reads, [a for a in alternatives if a is not g])
                        if not which:
                            continue
                        sib_start, sib_end = _read_extent(g, which, spans)
                        sib_ids = _emit(
                            m.chain, key_reads, sib_start, sib_end,
                            None, family_id,
                            junction_depths=junction_clean_depths(
                                m.chain, m.read_ids, folded_reads, indels_by_read),
                            extent_of=g.candidate_id,
                        )
                        if sib_ids:
                            ids.extend(sib_ids)
                # folded exact-sub-chain shadows keyed by the candidate(s) it became
                # (5'-TSS recovery). Mono members carry none.
                member_shadows = [
                    (tuple(f.chain.introns), tuple(sorted(f.read_ids)))
                    for f in m.folded if f.chain.introns
                ]
                if member_shadows:
                    for cid in out_ids:
                        shadow_map.setdefault(cid, []).extend(member_shadows)
            else:
                # mono (single-exon): split the pooled chain-less reads into per-locus
                # groups so distinct single-exon genes are not bridged, then overlap-match
                # each locus to a single-exon GTF (else emit a per-locus novel candidate).
                for group in _spatial_read_clusters(m.read_ids, spans):
                    gs = [spans[r][0] for r in group]
                    ge = [spans[r][1] for r in group]
                    lo, hi = min(gs), max(ge)
                    g_mono = _best_overlap_gtf(mono_gtf, lo, hi)
                    out_ids = _emit(m.chain, set(group), lo, hi, g_mono)
                    if out_ids:
                        ids.extend(out_ids)
                    if g_mono is None or read_extent_contradicted is None:
                        continue
                    # The reads' own extent is a SEPARATE hypothesis from an
                    # annotated single-exon model they overlap, emitted only
                    # when the reads CONTRADICT that model's ends under the
                    # sample's endpoint kernels (a polyA read past the claimed
                    # TES, or a read starting upstream of the claimed TSS) and
                    # no OTHER overlapping single-exon annotation explains them.
                    # EVERY overlapping annotation is examined -- not only the
                    # best-overlap one the group's reads were attached to --
                    # because the pooled span may pick the long model while the
                    # short one is the contradicted (and repairable) hypothesis.
                    # Otherwise the annotation IS the hypothesis: a duplicate
                    # would merely split its reads. The assignment's endpoint
                    # likelihood then decides between the two.
                    overlapping = [gc for gc, gs_, ge_ in mono_gtf if gs_ < hi and ge_ > lo]
                    for gc in overlapping:
                        which = read_extent_contradicted(
                            gc, group, [a for a in overlapping if a is not gc])
                        if not which:
                            continue
                        sib_start, sib_end = _read_extent(gc, which, spans)
                        sib_ids = _emit(m.chain, set(group), sib_start, sib_end, None,
                                        extent_of=gc.candidate_id)
                        if sib_ids:
                            ids.extend(sib_ids)
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

    # Step 5a instantiation: aggregate proposals by corrected chain. An existing
    # candidate (any source) gains correction support; a chain no candidate has
    # becomes a NEW novel candidate when >= new_chain_min_reads distinct reads
    # proposed exactly it. Baseline candidates are never removed or re-collapsed.
    if jc_result is not None and jc_result.proposals and junction_correction_apply:
        from fin.analysis.junction_correction import aggregate_proposals

        min_votes = int(getattr(junction_correction_apply, "min_reads", 2))
        by_chain_cand = {}
        for c in candidates:
            if c.intron_chain.introns:
                by_chain_cand.setdefault(tuple(c.intron_chain.introns), c)
        voter_span = {p.read_id: (p.reference_start, p.reference_end) for p in jc_result.proposals}
        n_new = n_merged = 0
        for chain_t, voters in aggregate_proposals(jc_result.proposals, min_reads=1).items():
            existing = by_chain_cand.get(chain_t)
            if existing is not None:
                existing.supporting_read_ids.update(voters)
                existing.correction_supporting_read_ids.update(voters)
                n_merged += 1
                continue
            if len(voters) < min_votes:
                continue
            starts = [voter_span[v][0] for v in voters]
            ends = [voter_span[v][1] for v in voters]
            st, en = min(starts), max(ends)
            seq = _build_spliced_sequence(genome_fasta, st, en, IntronChain(introns=chain_t), strand)
            if not seq:
                continue
            cid = _generate_novel_id(interval.chrom, strand, st, en, IntronChain(introns=chain_t))
            if cid in {c.candidate_id for c in candidates}:
                continue
            # family placement: the family of the voters' raw chains when unique,
            # else an independent correction family (never rewrites baseline ids)
            fam_ids = {by_chain_cand[p.original_chain].family_id
                       for p in jc_result.proposals
                       if p.read_id in voters and p.original_chain in by_chain_cand}
            fam_id = next(iter(fam_ids)) if len(fam_ids) == 1 else _generate_family_id(
                interval.chrom, strand, [IntronChain(introns=chain_t)])
            candidates.append(TranscriptCandidate(
                candidate_id=cid, intron_chain=IntronChain(introns=chain_t),
                three_prime_pos=(en if strand == "+" else st), sequence=seq,
                source="novel", supporting_read_ids=set(voters),
                correction_supporting_read_ids=set(voters),
                chrom=interval.chrom, strand=strand, start=st, end=en, family_id=fam_id))
            cluster_id_lists.append([cid])
            by_chain_cand[chain_t] = candidates[-1]
            n_new += 1
        if n_new or n_merged:
            logger.info("Interval %s: junction correction added %d candidate(s), "
                        "reinforced %d", interval.region_string, n_new, n_merged)

    logger.info(
        "Interval %s: %d candidates (%d GTF, %d novel) from %d reads [chain-cluster]",
        interval.region_string, len(candidates),
        sum(1 for c in candidates if c.source == "gtf"),
        sum(1 for c in candidates if c.source == "novel"), len(all_read_ids))

    return CandidateSet(
        interval=interval, candidates=candidates, read_ids=all_read_ids,
        read_sequences=read_sequences, clusters=cluster_id_lists,
        shadows=shadow_map, read_spans=spans, indels_by_read=indels_by_read,
        read_junctions=read_junctions, read_strands=read_strands,
        junction_support={j: (v[0], v[1]) for j, v in junction_counts.items()})


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
    collect_chain_evidence: bool = True,   # retained for API compatibility; evidence is always collected
    junction_correction=None,
    junction_correction_apply=None,
    read_extent_contradicted=None,
    collect_compatible_depths: bool = False,
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
            # Secondary mappings are alternative placements of the same query,
            # not independent structural support. Letting them enter generation
            # breaks the family invariant that one read belongs to one variant.
            if rd.get("is_secondary"):
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
            collect_chain_evidence=collect_chain_evidence,
            junction_correction=junction_correction,
            junction_correction_apply=junction_correction_apply,
            read_extent_contradicted=read_extent_contradicted,
            collect_compatible_depths=collect_compatible_depths,
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
                    candidate_id=_generate_novel_id(
                        interval.chrom, strand, ref_start, ref_end, chain
                    ),
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
