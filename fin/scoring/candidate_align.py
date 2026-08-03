"""Family-scoped read -> candidate alignment goodness (the M1/assignment-redesign foundation).

Within one family, RE-MAP every read to every candidate transcript sequence with mappy (minimap2).
A read->transcript alignment is a much easier/cleaner problem than the splice-aware genome mapping
that produced the input BAM, so we align to the candidates' spliced sequences directly and judge
compatibility/identity from the QUALITY of that alignment.

Efficiency: build ONE multi-sequence mappy index per family (all candidate sequences in a single
temp FASTA) and map each read ONCE (`best_n` large). Near-identical (wobble) candidates are returned
as secondary hits, NOT suppressed (verified empirically); we ignore the primary/secondary flag and
key every hit by its candidate (`hit.ctg`). This is O(1 index build + reads map calls) per family,
vs O(reads x candidates) if a separate aligner were built per candidate.

Per-(read, candidate) goodness vector (all from mappy fields; no minimap2 CLI):
  AS        reconstructed map-ont alignment score (:func:`score_hit`). A single internal indel > 50bp
            (an exon-sized structural difference) makes the alignment structurally INCOMPATIBLE with
            that candidate -> score_hit returns None and the hit is dropped. Terminal containment (a
            missing 5'/3' exon) is a soft-clip, NOT an internal indel, so it is kept (lower coverage).
  q_cov     query coverage = aligned read fraction = (q_en - q_st) / read_len.
  event_id  event identity = mlen / (mlen + substitutions + gap_opens); a gap counts as ONE event
            regardless of length (minimap2's event-identity convention).
  aln_q     aligned query length (q_en - q_st).
  s5, s3    query-FORWARD soft-clip (unaligned read prefix/suffix, in the read's own orientation).
            The directional containment signal (a longer candidate's unique 5' exon becomes extra
            soft-clip on the shorter sibling). For a reverse hit (``strand == -1``) the transcript-5'/
            3' sense of s5/s3 is SWAPPED; the caller uses ``strand`` to orient 5'-specific logic.
  strand    hit strand (+1 forward, -1 reverse) relative to the candidate transcript sequence.
  r_st,r_en candidate-relative aligned span.
  cand3p    candidate 3' distance = candidate_len - r_en (dRNA is poly(A)-anchored: the alignment
            should reach the candidate 3' end).

Compatibility (permissive) vs identification (contrastive) decisions are made by the CALLER; this
module only produces the raw goodness. See :func:`is_compatible` for a NanoCount-style default gate.
"""
from __future__ import annotations

import os
import tempfile
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional

from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit


@dataclass
class Goodness:
    """Goodness of one read's best alignment to one candidate (see module docstring)."""
    candidate_id: str
    AS: float
    q_cov: float
    event_id: float
    aln_q: int
    s5: int
    s3: int
    strand: int
    r_st: int
    r_en: int
    cand3p: int


def _hit_goodness(cand_id: str, h, qlen: int, cand_len: int) -> Optional[Goodness]:
    """Goodness of a single mappy hit, or None if score_hit rejects it (>50bp internal indel =
    exon-sized structural incompatibility)."""
    s = score_hit(h)
    if s is None:
        return None
    ins = del_ = gap_opens = 0
    for ln, op in (h.cigar or ()):
        if op == 1:            # insertion (query has extra bases vs candidate)
            ins += ln
            gap_opens += 1
        elif op == 2:          # deletion (candidate has extra bases vs query)
            del_ += ln
            gap_opens += 1
    subs = max(0, h.NM - ins - del_)
    denom = h.mlen + subs + gap_opens
    aln_q = h.q_en - h.q_st
    return Goodness(
        candidate_id=cand_id, AS=float(s),
        q_cov=(aln_q / qlen if qlen else 0.0),
        event_id=(h.mlen / denom if denom else 0.0),
        aln_q=aln_q, s5=h.q_st, s3=qlen - h.q_en, strand=h.strand,
        r_st=h.r_st, r_en=h.r_en, cand3p=cand_len - h.r_en)


def align_reads_to_candidates(
    candidates: Iterable,
    read_seqs: Dict[str, str],
    *,
    preset: Optional[str] = None,
    best_n: int = 100,
) -> Dict[str, Dict[str, Goodness]]:
    """Map every read to every candidate ONCE via a single multi-sequence mappy index.

    Args:
        candidates: iterable of objects with ``candidate_id`` (str) and ``sequence`` (str). Candidates
            with an empty sequence are skipped.
        read_seqs: ``{read_id: sequence}``. Empty sequences are skipped.
        preset: mappy preset (default :func:`get_m1_preset`, i.e. map-ont).
        best_n: max hits per read to retain (>= number of candidates so every near-optimal candidate
            hit is returned; matches NanoCount's ``-N 100``).

    Returns:
        ``{read_id: {candidate_id: Goodness}}`` keeping the best-AS hit per candidate. A read with no
        acceptable hit to any candidate is absent from the result.
    """
    import mappy

    preset = preset or get_m1_preset()
    # Index each candidate under a SYNTHETIC safe contig name ("c<idx>") -- raw candidate_ids may
    # contain whitespace (mappy truncates FASTA headers at the first space) or duplicates; both would
    # corrupt the hit->candidate mapping. Dedup by candidate_id (first wins).
    name_cid: Dict[str, str] = {}
    name_len: Dict[str, int] = {}
    entries: List = []                       # (safe_name, sequence)
    seen: set = set()
    for c in candidates:
        seq = getattr(c, "sequence", None)
        cid = c.candidate_id
        if not seq or cid in seen:
            continue
        seen.add(cid)
        name = f"c{len(entries)}"
        name_cid[name] = cid
        name_len[name] = len(seq)
        entries.append((name, seq))
    if not entries or not read_seqs:
        return {}

    path = None
    try:
        with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as tmp:
            path = tmp.name
            for name, seq in entries:
                tmp.write(f">{name}\n{seq}\n")
        aligner = mappy.Aligner(path, preset=preset, best_n=best_n)
        out: Dict[str, Dict[str, Goodness]] = {}
        for rid, seq in read_seqs.items():
            if not seq:
                continue
            qlen = len(seq)
            best: Dict[str, Goodness] = {}
            for h in aligner.map(seq):
                cid = name_cid.get(h.ctg)
                if cid is None:              # defensive: hit to an unknown contig
                    continue
                g = _hit_goodness(cid, h, qlen, name_len[h.ctg])
                if g is None:
                    continue
                cur = best.get(cid)
                if cur is None or g.AS > cur.AS:
                    best[cid] = g
            if best:
                out[rid] = best
        return out
    finally:
        if path:
            try:
                os.unlink(path)
            except OSError:
                pass


def is_compatible(
    g: Goodness,
    best_event_id: float = 0.0,
    *,
    min_q_cov: float = 0.80,
    min_event_id: float = 0.80,
    event_id_margin: float = 0.05,
    max_cand3p: int = 50,
) -> bool:
    """Permissive NanoCount-style compatibility gate (existence is decided contrastively by the
    caller, NOT by this gate). A read is a plausible fit to the candidate iff it covers most of the
    read, fits with acceptable event identity (absolute floor OR within a margin of the read's best),
    and reaches the candidate 3' end (dRNA poly(A) anchoring). Candidate coverage is deliberately NOT
    gated -- a 5'-degraded read covering only a long transcript's 3' suffix is still compatible.
    """
    return (g.q_cov >= min_q_cov
            and g.event_id >= max(min_event_id, best_event_id - event_id_margin)
            and g.cand3p <= max_cand3p)
