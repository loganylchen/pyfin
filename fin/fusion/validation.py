"""Whole-molecule sequence validation against a complete stitched transcript."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

from fin.fusion.junctions import arm_sequence


@dataclass(frozen=True)
class TranscriptValidation:
    status: str
    reason: str
    score: Optional[float] = None
    query_coverage: float = 0.0
    identity: float = 0.0
    junction_anchors: Tuple[Tuple[int, int], ...] = ()
    # per template block: (aligned template positions, block length); an
    # unsupported block -- an intron carried into the template by a summary
    # CIGAR, an exon the read never aligned -- is what this exposes
    block_support: Tuple[Tuple[int, int], ...] = ()
    raw_query_coverage: float = 0.0
    # per boundary, per side: (aligned, matched, matches required, unaligned
    # run at the boundary) -- what a rejection for missing_junction_anchor
    # actually failed on
    anchor_detail: Tuple[Tuple[Tuple[int, int, int, int], Tuple[int, int, int, int]], ...] = ()


def stitched_template(path, genome):
    pieces = []
    boundaries = []
    offset = 0
    for arm in path.arms:
        sequence = arm_sequence(arm, genome)
        if not sequence:
            return "", ()
        for block in arm.blocks:
            if offset:
                boundaries.append(offset)
            offset += block.reference_end - block.reference_start
        pieces.append(sequence)
    return "".join(pieces), tuple(boundaries)


def template_blocks(path):
    """Template intervals of every observed block, in template order."""
    spans = []
    offset = 0
    for arm in path.arms:
        for block in arm.blocks:
            length = block.reference_end - block.reference_start
            spans.append((offset, offset + length))
            offset += length
    return tuple(spans)


def _flank_support(aligned, matched, positions):
    """(aligned, matched, unaligned run immediately adjoining the boundary)
    over an ordered tuple of template positions running AWAY from the
    boundary. The run stops at the first aligned position; deletions deeper
    in the flank are governed by the global indel and template-run limits."""
    a = m = 0
    run = 0
    counting = True
    for position in positions:
        if position in aligned:
            a += 1
            counting = False
            if position in matched:
                m += 1
        elif counting:
            run += 1
    return a, m, run


def validate_stitched_path(path, query, genome, *, min_anchor=12, min_identity=0.7,
                           min_coverage=0.9, max_indel=10, max_hits=32, domain=(0, 0),
                           anchor_window=None, anchor_max_edge_gap=4, min_block_fraction=0.7):
    """Require one forward alignment that represents the whole stitched template.

    Three things are checked on the SAME forward alignment, in template order:

    * every exon/fusion boundary is anchored on both sides: within the
      ``anchor_window`` (default ``2 * min_anchor``) flank on each side, at
      least ``min_anchor`` template positions are aligned, at least
      ``ceil(min_identity * aligned)`` of them match, and no more than
      ``anchor_max_edge_gap`` consecutive template positions right next to the
      boundary are unaligned. Deletions inside the flank are tolerated (they
      contribute nothing); a nanopore read with a 2-bp deletion near a real
      junction is not rejected for it;
    * every template block is REPRESENTED: at least ``min(min_anchor, L)``
      aligned positions and at least ``min_block_fraction`` of its length, no
      unaligned template run longer than ``max_indel`` inside the alignment,
      and no unrepresented template prefix/suffix longer than ``max_indel``.
      This closes the loophole where a template carrying 700 bp of intronic
      sequence "validated" because the read simply never aligned that part;
    * query coverage over the molecule-level ``domain`` (terminal bases no
      detailed alignment could place are excluded from the denominator).
    """
    import math
    import mappy

    if (min(min_anchor, max_hits) <= 0 or max_indel < 0 or anchor_max_edge_gap < 0
            or not 0 < min_identity <= 1 or not 0 < min_coverage <= 1 or not 0 < min_block_fraction <= 1):
        raise ValueError("invalid stitched-transcript validation settings")
    if anchor_window is None:
        anchor_window = 2 * min_anchor
    if anchor_window < min_anchor or min(domain) < 0:
        raise ValueError("invalid anchor window or coverage domain")
    blocks = template_blocks(path)
    target, boundaries = stitched_template(path, genome)
    query = query.upper().replace("U", "T") if query else ""
    if not target or not query:
        return TranscriptValidation("unresolved", "missing_sequence_or_template")
    if set(query + target) - set("ACGT"):
        return TranscriptValidation("unresolved", "ambiguous_sequence")
    aligner = mappy.Aligner(seq=target, preset="map-ont", best_n=max_hits)
    if not aligner:
        return TranscriptValidation("unresolved", "template_index_failure")
    results = []
    hits = []
    for hit in aligner.map(query):
        hits.append(hit)
        if len(hits) >= max_hits:
            return TranscriptValidation("unresolved", "template_hit_limit")
    for hit in hits:
        if hit.strand != 1 or not hit.cigar:
            continue
        query_position, target_position = hit.q_st, hit.r_st
        matches = mismatches = inserted = deleted = gap_cost = longest_gap = 0
        aligned_target = set()
        matched_target = set()
        invalid = False
        for length, operation in hit.cigar:
            if operation in (0, 7, 8):
                for offset in range(length):
                    position = target_position + offset
                    aligned_target.add(position)
                    if query[query_position + offset] == target[position]:
                        matches += 1
                        matched_target.add(position)
                    else:
                        mismatches += 1
                query_position += length
                target_position += length
            elif operation in (1, 2):
                longest_gap = max(longest_gap, length)
                gap_cost += 4 + 2 * length
                if operation == 1:
                    inserted += length
                    query_position += length
                else:
                    deleted += length
                    target_position += length
            else:
                invalid = True
                break
        if invalid or query_position != hit.q_en or target_position != hit.r_en:
            continue
        raw_coverage = (hit.q_en - hit.q_st) / len(query)
        lead, trail = domain
        lo, hi = lead, len(query) - trail
        covered = max(0, min(hit.q_en, hi) - max(hit.q_st, lo))
        coverage = covered / (hi - lo) if hi > lo else 0.0
        identity = matches / max(1, matches + mismatches + inserted + deleted)
        anchors = []
        detail = []
        anchored = True
        for boundary in boundaries:
            left = tuple(range(boundary - 1, max(-1, boundary - anchor_window - 1), -1))
            right = tuple(range(boundary, min(len(target), boundary + anchor_window)))
            la, lm, lgap = _flank_support(aligned_target, matched_target, left)
            ra, rm, rgap = _flank_support(aligned_target, matched_target, right)
            anchors.append((la, ra))
            detail.append(((la, lm, math.ceil(min_identity * la), lgap), (ra, rm, math.ceil(min_identity * ra), rgap)))
            for aligned, matched, gap in ((la, lm, lgap), (ra, rm, rgap)):
                if aligned < min_anchor or matched < math.ceil(min_identity * aligned) or gap > anchor_max_edge_gap:
                    anchored = False
        # template representation: every block, and the template ends
        support = []
        represented = True
        for start, end in blocks:
            aligned = sum(1 for position in range(start, end) if position in aligned_target)
            support.append((aligned, end - start))
            if aligned < min(min_anchor, end - start) or aligned < min_block_fraction * (end - start):
                represented = False
        ends_ok = hit.r_st <= max_indel and len(target) - hit.r_en <= max_indel
        longest_unaligned = 0
        run = 0
        for position in range(hit.r_st, hit.r_en):
            if position in aligned_target:
                run = 0
            else:
                run += 1
                longest_unaligned = max(longest_unaligned, run)
        clipped = len(query) - (hit.q_en - hit.q_st)
        score = 2 * matches - 4 * mismatches - gap_cost - 4 * clipped
        reason = ("insufficient_query_coverage" if coverage < min_coverage
                  else "low_sequence_identity" if identity < min_identity
                  else "large_indel" if longest_gap > max_indel or longest_unaligned > max_indel
                  else "unrepresented_template_end" if not ends_ok
                  else "insufficient_template_block_support" if not represented
                  else "missing_junction_anchor" if not anchored
                  else "whole_transcript_anchored")
        results.append(TranscriptValidation("validated" if reason == "whole_transcript_anchored" else "rejected",
                                             reason, score, coverage, identity, tuple(anchors),
                                             tuple(support), raw_coverage, tuple(detail)))
    if not results:
        return TranscriptValidation("rejected", "no_forward_template_alignment")
    return max(results, key=lambda result: (result.status == "validated", result.score,
                                            result.query_coverage, result.identity))
