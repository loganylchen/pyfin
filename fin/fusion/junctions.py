"""Junction nomination and read-window contrasts, isolated from assembly gates."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Mapping, Optional, Tuple

from fin.candidates.canonical import MOTIFS_EXTENDED, canonical_intron_alts
from fin.fusion.alignment import JunctionAlignment, align_junction_window
from fin.fusion.annotation import FusionArm
from fin.utils.sequences import reverse_complement


@dataclass(frozen=True)
class JunctionContrast:
    positions: Tuple[int, int]
    alignment: JunctionAlignment
    template: str
    origin: Tuple[str, ...]


@dataclass(frozen=True)
class JunctionDecision:
    observed_positions: Tuple[int, int]
    selected_positions: Optional[Tuple[int, int]]
    status: str
    reason: str
    query_start: int
    query_end: int
    contrasts: Tuple[JunctionContrast, ...]
    margin: Optional[float]
    signal_status: str = "not_scored"


def _base_slice(genome, chrom, start, end):
    if chrom not in genome or start < 0 or end <= start:
        return ""
    sequence = genome[chrom]
    if end > len(sequence):
        return ""
    return sequence[start:end].upper()


def fusion_motif(genome, left: FusionArm, right: FusionArm, positions):
    """Read donor/acceptor dinucleotides in each arm's transcript orientation."""
    donor, acceptor = positions
    left_bases = (_base_slice(genome, left.chrom, donor, donor + 2) if left.strand == "+"
                  else reverse_complement(_base_slice(genome, left.chrom, donor - 2, donor)))
    right_bases = (_base_slice(genome, right.chrom, acceptor - 2, acceptor) if right.strand == "+"
                   else reverse_complement(_base_slice(genome, right.chrom, acceptor, acceptor + 2)))
    return left_bases, right_bases


def move_junction(left: FusionArm, right: FusionArm, positions):
    """Modify only the junction-facing exon boundaries, never read coordinates."""
    donor, acceptor = positions
    left_block = left.blocks[-1]
    right_block = right.blocks[0]
    if left.strand == "+":
        left_block = replace(left_block, reference_end=donor)
    else:
        left_block = replace(left_block, reference_start=donor)
    if right.strand == "+":
        right_block = replace(right_block, reference_start=acceptor)
    else:
        right_block = replace(right_block, reference_end=acceptor)
    moved_left = replace(left, blocks=left.blocks[:-1] + (left_block,))
    moved_right = replace(right, blocks=(right_block,) + right.blocks[1:])
    for arm in (moved_left, moved_right):
        for previous, following in zip(arm.blocks, arm.blocks[1:]):
            if ((arm.strand == "+" and previous.reference_end > following.reference_start)
                    or (arm.strand == "-" and previous.reference_start < following.reference_end)):
                raise ValueError("junction edit overlaps another exon")
    return moved_left, moved_right


def arm_sequence(arm: FusionArm, genome: Mapping[str, str]):
    pieces = []
    for block in arm.blocks:
        sequence = _base_slice(genome, arm.chrom, block.reference_start, block.reference_end)
        if not sequence:
            return ""
        pieces.append(sequence if arm.strand == "+" else reverse_complement(sequence))
    return "".join(pieces)


def nominate_junctions(
    left, right, genome, *, tolerance=3, supported_positions=(),
    annotation_positions=(), assembly_positions=(), internal_intron=False, max_candidates=128,
):
    """Keep the observed baseline and nominate sequence/evidence alternatives."""
    if tolerance < 0 or max_candidates <= 0:
        raise ValueError("invalid junction nomination limits")
    observed = left.exit_position, right.entry_position
    origins = {observed: {"observed"}}
    for source, positions_set in (("read_supported", supported_positions), ("annotation", annotation_positions),
                                 ("final_assembly_clean", assembly_positions)):
        for positions in positions_set:
            positions = tuple(positions)
            if len(positions) != 2:
                raise ValueError("junction hypothesis must contain two positions")
            if max(abs(value - baseline) for value, baseline in zip(positions, observed)) <= tolerance:
                origins.setdefault(positions, set()).add(source)
    if internal_intron:
        if left.chrom != right.chrom or left.strand != right.strand:
            raise ValueError("an internal intron must stay on the same chromosome and strand")
        if left.chrom in genome:
            intron = observed if left.strand == "+" else tuple(reversed(observed))
            alternatives = canonical_intron_alts(
                intron[0], intron[1], left.strand, genome[left.chrom], MOTIFS_EXTENDED, tolerance,
            )
            for alternative in alternatives:
                positions = alternative if left.strand == "+" else tuple(reversed(alternative))
                origins.setdefault(positions, set()).add("canonical")
    else:
        for donor in range(observed[0] - tolerance, observed[0] + tolerance + 1):
            for acceptor in range(observed[1] - tolerance, observed[1] + tolerance + 1):
                positions = donor, acceptor
                if fusion_motif(genome, left, right, positions) in MOTIFS_EXTENDED:
                    origins.setdefault(positions, set()).add("canonical")
    valid = {}
    for positions, sources in origins.items():
        try:
            moved_left, moved_right = move_junction(left, right, positions)
        except ValueError:
            continue
        if internal_intron:
            if ((left.strand == "+" and positions[0] >= positions[1])
                    or (left.strand == "-" and positions[0] <= positions[1])):
                continue
        if arm_sequence(moved_left, genome) and arm_sequence(moved_right, genome):
            valid[positions] = tuple(sorted(sources))
    if len(valid) > max_candidates:
        return {}, "candidate_limit"
    if observed not in valid:
        return {}, "invalid_observed_template"
    return dict(sorted(valid.items())), "complete_nomination"


def evaluate_junction(
    read_sequence, left, right, genome, *, tolerance=3, flank=35, min_anchor=12,
    min_identity=0.7, min_margin=6.0, supported_positions=(), annotation_positions=(), assembly_positions=(),
    internal_intron=False, max_candidates=128,
) -> JunctionDecision:
    """Compare all nominated coordinates on one fixed, fully consumed query window."""
    if flank < min_anchor or min_margin < 0:
        raise ValueError("invalid flank or decision margin")
    observed = left.exit_position, right.entry_position
    query_start = max(left.query_start, left.query_end - flank)
    query_end = min(right.query_end, right.query_start + flank)
    if (not read_sequence or query_start < 0 or query_end > len(read_sequence)
            or query_end <= query_start or left.query_end - right.query_start > min_anchor):
        return JunctionDecision(observed, None, "unresolved", "invalid_query_window",
                                query_start, query_end, (), None)
    candidates, reason = nominate_junctions(
        left, right, genome, tolerance=tolerance, supported_positions=supported_positions,
        annotation_positions=annotation_positions, assembly_positions=assembly_positions, internal_intron=internal_intron,
        max_candidates=max_candidates,
    )
    if not candidates:
        return JunctionDecision(observed, None, "unresolved", reason, query_start, query_end, (), None)
    query = read_sequence[query_start:query_end]
    contrasts = []
    padding = flank + tolerance
    for positions, origin in candidates.items():
        moved_left, moved_right = move_junction(left, right, positions)
        left_sequence = arm_sequence(moved_left, genome)[-padding:]
        right_sequence = arm_sequence(moved_right, genome)[:padding]
        template = left_sequence + right_sequence
        alignment = align_junction_window(
            query, template, len(left_sequence), min_anchor=min_anchor, min_identity=min_identity,
        )
        contrasts.append(JunctionContrast(positions, alignment, template, origin))
    contrasts = tuple(contrasts)
    scored = [contrast for contrast in contrasts if contrast.alignment.status == "scored"]
    if not scored:
        return JunctionDecision(observed, None, "unresolved", "no_anchored_hypothesis",
                                query_start, query_end, contrasts, None)
    ranked = sorted(scored, key=lambda contrast: (-contrast.alignment.score, contrast.positions))
    best = ranked[0]
    failed = [contrast for contrast in contrasts if contrast.alignment.status != "scored"]
    if any(contrast.alignment.reason in ("window_limit", "ambiguous_sequence", "invalid_window")
           for contrast in failed):
        return JunctionDecision(observed, None, "unresolved", "incomparable_hypotheses",
                                query_start, query_end, contrasts, None)
    if len(ranked) == 1:
        if best.positions != observed:
            return JunctionDecision(observed, None, "unresolved", "baseline_not_anchored",
                                    query_start, query_end, contrasts, None)
        return JunctionDecision(observed, observed, "retained", "observed_anchored",
                                query_start, query_end, contrasts, None)
    margin = best.alignment.score - ranked[1].alignment.score
    if margin <= 0 or margin < min_margin:
        return JunctionDecision(observed, None, "unresolved", "indistinguishable_coordinates",
                                query_start, query_end, contrasts, margin)
    if next(contrast for contrast in contrasts if contrast.positions == observed).alignment.status != "scored":
        return JunctionDecision(observed, None, "unresolved", "baseline_not_anchored",
                                query_start, query_end, contrasts, margin)
    status = "retained" if best.positions == observed else "correction_proposed"
    return JunctionDecision(observed, best.positions, status, "anchored_sequence_margin",
                            query_start, query_end, contrasts, margin)
