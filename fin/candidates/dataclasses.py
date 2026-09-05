"""Data classes for transcript candidate discovery."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple

from fin.io.interval_manager import GenomicInterval


@dataclass(frozen=True)
class IntronChain:
    """Hashable representation of a transcript's intron structure.

    Each intron is a (start, end) tuple in genomic coordinates.
    """

    introns: Tuple[Tuple[int, int], ...]

    @property
    def num_introns(self) -> int:
        return len(self.introns)

    @property
    def is_single_exon(self) -> bool:
        return len(self.introns) == 0


@dataclass
class TranscriptCandidate:
    """A candidate transcript (from GTF annotation, novel discovery, or fusion).

    Required fields (positional):
        candidate_id, intron_chain, three_prime_pos, sequence, source,
        supporting_read_ids, chrom, strand, start, end.

    Optional fusion fields (populated only when source == "fusion"):
        fusion_junction: Tuple of (left_genomic_pos, right_genomic_pos) marking
            where the two parental transcripts are joined.
        breakpoint_left: (chrom, pos, strand) for the left fusion partner.
        breakpoint_right: (chrom, pos, strand) for the right fusion partner.
    """

    candidate_id: str
    intron_chain: IntronChain
    three_prime_pos: int
    sequence: str
    source: str  # "gtf", "novel", or "fusion"
    supporting_read_ids: Set[str]
    chrom: str
    strand: str
    start: int
    end: int
    # Stable discovery splice-family identity. None for mono, fusion, and
    # legacy discovery paths that do not construct structural families.
    family_id: Optional[str] = None

    # Complete-splice-chain evidence (families discovery only; ALWAYS
    # collected - one CIGAR walk per read - and recorded on the candidate;
    # only the opt-in gates act on it).
    # An exact complete-chain read's OWN raw alignment (primary, MAPQ >= 1)
    # equals this candidate's full intron chain (every junction traversed at
    # the same coordinates; the 5' terminus may sit anywhere inside the
    # first exon - dRNA 5' degradation tolerated); folded sub-chain
    # (degraded) reads never count. A read INDEL-cleanly attests one
    # junction when it has no insertion/deletion within 10 bp of that
    # junction's boundaries (cigar-level only: mismatches/base quality are
    # NOT inspected - documented limitation).
    # complete_chain_junction_depths[i] = indel-clean attestation count of
    # junction i (5'->3' genomic order); the weakest field is min(depths).
    # () / 0 for mono, fusion, legacy paths, and disabled collection.
    complete_chain_junction_depths: Tuple[int, ...] = ()
    complete_chain_weakest_indel_clean_depth: int = 0
    # Provenance split of supporting_read_ids (which stays the union for every
    # existing consumer): reads whose RAW alignment chain produced this candidate
    # vs reads whose signal ballot proposed it (step-5a additive correction).
    # A correction vote is hypothesis support, not raw complete-chain attestation.
    raw_supporting_read_ids: Set[str] = field(default_factory=set)
    correction_supporting_read_ids: Set[str] = field(default_factory=set)

    # Fusion-specific (populated only when source == "fusion")
    fusion_junction: Optional[Tuple[int, int]] = None
    breakpoint_left: Optional[Tuple[str, int, str]] = None   # (chrom, pos, strand)
    breakpoint_right: Optional[Tuple[str, int, str]] = None

    @property
    def length(self) -> int:
        # Fusion candidates set start=posA on chromA and end=posB on chromB.
        # Subtracting them is meaningless (and can be negative) so we fall
        # back to the assembled sequence length, which is the only thing
        # downstream code actually uses for fusion candidates.
        if self.source == "fusion":
            return len(self.sequence) if self.sequence else 0
        return self.end - self.start

    @property
    def num_supporting_reads(self) -> int:
        return len(self.supporting_read_ids)


@dataclass
class CandidateSet:
    """Collection of transcript candidates for a genomic interval."""

    interval: GenomicInterval
    candidates: List[TranscriptCandidate]
    read_ids: Set[str]
    # read_id -> primary-alignment query_sequence, scoped to reads fetched
    # within this interval. Populated by discover_candidates(); used by R1
    # ablation (mappy argmax bypass) and any path that needs read FASTA.
    read_sequences: Dict[str, str] = field(default_factory=dict)
    # Generation-cluster grouping for quant_mode="cluster": each inner list holds
    # the candidate_ids that belong to one generation cluster. Ids (not indices)
    # survive post-discovery gate drops/reorders (canonical/dominance/fusion), which
    # would invalidate positional indices. Populated only by the chain-cluster
    # discovery path; None means "not chain-cluster mode" (the quant-cluster dispatch
    # then treats the whole candidate list as one unscoped cluster).
    clusters: Optional[List[List[str]]] = None
    # 5'-TSS short-isoform recovery provenance (chain-cluster path only). Maps a
    # container candidate_id to the SHADOW sub-chains folded into it during
    # generation: each shadow is (intron_chain_tuple, sorted_read_id_tuple). The
    # quant-cluster recovery step re-tests these folded truncations against the
    # container's read-5'-end pileup and re-activates real TSS peaks as short
    # isoforms. None means "not chain-cluster mode" (no recovery).
    shadows: Optional[
        Dict[str, List[Tuple[Tuple[Tuple[int, int], ...], Tuple[str, ...]]]]
    ] = None
    # Interval-level junction support (families path, always collected):
    # (donor, acceptor) -> (raw_reads, indel_clean_reads) over every MAPQ>=1
    # primary read junction observation. Recorded for every run; consumed by
    # the opt-in junction election / correction. None on non-families paths.
    junction_support: Optional[
        Dict[Tuple[int, int], Tuple[int, int]]
    ] = None
    # read_id -> (genomic reference_start, reference_end) for reads fetched in this
    # interval. Populated by the chain-cluster discovery path; used by the post-EM
    # mono-exon resolution (quant_mode="m2_em", mono_resolve_post_em) to test whether a
    # single-exon read's genomic span lies inside a surviving multi-exon candidate's
    # exon. None means "not chain-cluster mode" (mono-resolution no-ops).
    read_spans: Optional[Dict[str, Tuple[int, int]]] = None
    # read_id -> tuple of (lo, hi) reference indel spans (D as [lo,hi), I as a
    # point) for MAPQ>=1 primary reads of this interval; the raw CIGAR record
    # behind complete-chain depths. Lets post-EM consumers (folded-TSS
    # recovery) attest a re-emitted chain from CIGARs, never from EM weights.
    indels_by_read: Optional[Dict[str, tuple]] = None
    # read_id -> the read's OWN primary-CIGAR intron chain ((donor, acceptor), ...)
    # and read_id -> mapped strand. Together with read_spans/indels_by_read these
    # are the C_r / S_r / E_r observations of the per-molecule evidence row
    # (fin/analysis/molecule_evidence.py). Empty chain = single-exon read.
    read_junctions: Optional[Dict[str, Tuple[Tuple[int, int], ...]]] = None
    read_strands: Optional[Dict[str, str]] = None

    @property
    def num_candidates(self) -> int:
        return len(self.candidates)

    def get_candidate(self, candidate_id: str) -> Optional[TranscriptCandidate]:
        for c in self.candidates:
            if c.candidate_id == candidate_id:
                return c
        return None

    def candidate_ids(self) -> List[str]:
        return [c.candidate_id for c in self.candidates]
