"""Shared canonical splice-site primitives (motif parsing, search, filter).

These pure functions are the single source of truth for what counts as a
*canonical* splice junction. Sharing one implementation guarantees the
consistency constraint: a ±N bp canonical SEARCH and the canonical FILTER
always agree on the motif set ("search what you filter").

Coordinate frame: introns are ``(s, e)`` 0-based half-open genomic bounds. The
(donor, acceptor) dinucleotides are read in TRANSCRIPT orientation:
  + strand: donor = genome[s:s+2],          acceptor = genome[e-2:e].
  - strand: donor = revcomp(genome[e-2:e]), acceptor = revcomp(genome[s:s+2]).

Public API:
    MOTIFS_STRICT, MOTIFS_EXTENDED   -- canonical (donor, acceptor) pair sets.
    parse_motifs(tokens)             -- "GT-AG" tokens -> frozenset of pairs.
    canonical_intron_alts(...)       -- enumerate canonical (s', e') near (s, e).
    chain_all_canonical(introns,...) -- True if every intron is canonical.
"""
from __future__ import annotations

from typing import FrozenSet, List, Optional, Sequence, Tuple

# Allowed (donor, acceptor) dinucleotide pairs in TRANSCRIPT orientation.
#   strict   = major U2 spliceosome only.
#   extended = + minor U12 motifs (GC-AG, AT-AC).
MOTIFS_STRICT: FrozenSet[Tuple[str, str]] = frozenset({("GT", "AG")})
MOTIFS_EXTENDED: FrozenSet[Tuple[str, str]] = frozenset(
    {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
)

_RC = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def _revcomp2(s: str) -> str:
    """Reverse-complement a 2-base string (unknown bases -> 'N')."""
    return _RC.get(s[1], "N") + _RC.get(s[0], "N")


def parse_motifs(
    tokens: Optional[Sequence[str]],
) -> FrozenSet[Tuple[str, str]]:
    """Parse motif tokens (e.g. ["GT-AG", "GC-AG"]) into a frozenset of pairs.

    Each token is split on '-' into a (donor, acceptor) pair; each side must be
    exactly two characters and is upper-cased. Empty/None input -> {("GT", "AG")}.
    The returned set is shared by BOTH the canonical SEARCH (junction snap /
    expand) and the canonical FILTER (``chain_all_canonical``) so they always
    agree on what counts as canonical (the consistency constraint: search what
    you filter). Raises ValueError on a malformed token.
    """
    if not tokens:
        return frozenset({("GT", "AG")})
    pairs = set()
    for tok in tokens:
        halves = tok.split("-")
        if len(halves) != 2:
            raise ValueError(
                f"bad motif token {tok!r}: expected DONOR-ACCEPTOR (e.g. GT-AG)"
            )
        donor, acceptor = halves[0].strip().upper(), halves[1].strip().upper()
        if len(donor) != 2 or len(acceptor) != 2:
            raise ValueError(
                f"bad motif token {tok!r}: each side must be exactly 2 bases"
            )
        pairs.add((donor, acceptor))
    return frozenset(pairs)


def canonical_intron_alts(
    s: int,
    e: int,
    strand: str,
    genome: str,
    motif_set: FrozenSet[Tuple[str, str]],
    search_bp: int,
) -> List[Tuple[int, int]]:
    """Enumerate canonical (s', e') intron bounds within +-search_bp of (s, e).

    Junction coords are 0-based half-open intron bounds (same frame as
    ``chain_all_canonical``); the (donor, acceptor) dinucleotides are read in
    TRANSCRIPT orientation exactly as there. Only alternatives whose motif is in
    ``motif_set`` are kept, sorted by (|Δs| + |Δe|, s', e') so the nearest comes
    first. Does NOT include (s, e) unless it is itself canonical.
    """
    n = len(genome)
    alts = []
    for ns in range(s - search_bp, s + search_bp + 1):
        for ne in range(e - search_bp, e + search_bp + 1):
            if ns < 0 or ne > n or ne - 2 < 0 or ns + 2 > n or ne <= ns:
                continue
            up = genome[ns:ns + 2].upper()
            dn = genome[ne - 2:ne].upper()
            if strand == "+":
                donor, acceptor = up, dn
            else:
                donor, acceptor = _revcomp2(dn), _revcomp2(up)
            if (donor, acceptor) in motif_set:
                alts.append((abs(ns - s) + abs(ne - e), ns, ne))
    alts.sort()
    return [(ns, ne) for _, ns, ne in alts]


def chain_all_canonical(
    introns: Sequence[Tuple[int, int]],
    genome: str,
    strand: str,
    motif_set: FrozenSet[Tuple[str, str]] = MOTIFS_STRICT,
) -> bool:
    """True if EVERY junction in ``introns`` is a canonical splice.

    Pure function of (junction genomic coords, strand, genome). Does NOT modify
    any coordinate (no snap); the junctions were fixed at discovery. Intron
    (s, e), 0-based half-open; donor/acceptor read in TRANSCRIPT orientation (see
    module docstring). A junction is canonical iff (donor, acceptor) is in
    ``motif_set``. Empty chain (single-exon) trivially passes.
    """
    if not genome:
        return False
    n = len(genome)
    for s, e in introns:
        if s < 0 or e > n or e - 2 < 0 or s + 2 > n:
            return False
        up = genome[s:s + 2].upper()      # genomic 5' dinuc
        dn = genome[e - 2:e].upper()      # genomic 3' dinuc
        if strand == "+":
            donor, acceptor = up, dn
        else:
            donor, acceptor = _revcomp2(dn), _revcomp2(up)
        if (donor, acceptor) not in motif_set:
            return False
    return True
