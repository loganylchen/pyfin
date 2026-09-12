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


def validate_stitched_path(path, query, genome, *, min_anchor=12, min_identity=0.7,
                           min_coverage=0.9, max_indel=10, max_hits=32):
    """Require one forward alignment spanning every exon and fusion boundary."""
    import mappy

    if (min(min_anchor, max_hits) <= 0 or max_indel < 0
            or not 0 < min_identity <= 1 or not 0 < min_coverage <= 1):
        raise ValueError("invalid stitched-transcript validation settings")
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
        coverage = (hit.q_en - hit.q_st) / len(query)
        identity = matches / max(1, matches + mismatches + inserted + deleted)
        anchors = []
        anchored = True
        for boundary in boundaries:
            left = set(range(max(0, boundary - min_anchor), boundary))
            right = set(range(boundary, min(len(target), boundary + min_anchor)))
            counts = len(left & aligned_target), len(right & aligned_target)
            anchors.append(counts)
            if (min(counts) < min_anchor or len(left & matched_target) < min_identity * min_anchor
                    or len(right & matched_target) < min_identity * min_anchor):
                anchored = False
        clipped = len(query) - (hit.q_en - hit.q_st)
        score = 2 * matches - 4 * mismatches - gap_cost - 4 * clipped
        reason = ("insufficient_query_coverage" if coverage < min_coverage
                  else "low_sequence_identity" if identity < min_identity
                  else "large_indel" if longest_gap > max_indel
                  else "missing_junction_anchor" if not anchored
                  else "whole_transcript_anchored")
        results.append(TranscriptValidation("validated" if reason == "whole_transcript_anchored" else "rejected",
                                             reason, score, coverage, identity, tuple(anchors)))
    if not results:
        return TranscriptValidation("rejected", "no_forward_template_alignment")
    return max(results, key=lambda result: (result.status == "validated", result.score,
                                            result.query_coverage, result.identity))
