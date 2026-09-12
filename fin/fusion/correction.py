"""Transactional, sequence-supported junction proposals for complete read paths."""

from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import product
from typing import Optional, Tuple

from fin.fusion.annotation import AnnotatedReadPath
from fin.fusion.junctions import JunctionDecision, evaluate_junction, move_junction


@dataclass(frozen=True)
class PathJunction:
    kind: str
    arm_index: int
    block_index: int
    decision: JunctionDecision


@dataclass(frozen=True)
class PathCorrection:
    observed: AnnotatedReadPath
    proposed: Optional[AnnotatedReadPath]
    status: str
    reason: str
    junctions: Tuple[PathJunction, ...]
    validation: Tuple[PathJunction, ...] = ()


def correction_alternatives(result, *, min_margin=6.0, max_paths=64):
    """Retain near-tied coordinate combinations without declaring them resolved."""
    if max_paths <= 0 or min_margin < 0:
        raise ValueError("invalid correction alternative limits")
    if result.proposed is not None:
        return (result.proposed,), "complete"
    choices = []
    count = 1
    for item in result.junctions:
        decision = item.decision
        if decision.selected_positions is not None:
            positions = (decision.selected_positions,)
        elif decision.reason == "indistinguishable_coordinates":
            scored = [contrast for contrast in decision.contrasts if contrast.alignment.status == "scored"]
            best = max(contrast.alignment.score for contrast in scored)
            positions = tuple(contrast.positions for contrast in scored
                              if contrast.alignment.score == best
                              or best - contrast.alignment.score < min_margin)
        else:
            return (), decision.reason
        choices.append(positions)
        count *= len(positions)
        if count > max_paths:
            return (), "correction_path_limit"
    if not choices:
        return (), "no_junctions"
    alternatives = []
    for positions in product(*choices):
        arms = list(result.observed.arms)
        try:
            for item, selected in zip(result.junctions, positions):
                arm = arms[item.arm_index]
                if item.kind == "intron":
                    left, right = move_junction(
                        replace(arm, blocks=arm.blocks[:item.block_index + 1]),
                        replace(arm, blocks=arm.blocks[item.block_index + 1:]), selected,
                    )
                    arms[item.arm_index] = replace(arm, blocks=left.blocks + right.blocks)
                else:
                    arms[item.arm_index], arms[item.arm_index + 1] = move_junction(
                        arm, arms[item.arm_index + 1], selected,
                    )
        except ValueError:
            continue
        alternatives.append(replace(result.observed, arms=tuple(arms)))
    return tuple(alternatives), "complete" if alternatives else "incompatible_joint_edits"


def path_junctions(path):
    """Yield internal splice and fusion edges with spliced flanking context."""
    for arm_index, arm in enumerate(path.arms):
        for block_index in range(len(arm.blocks) - 1):
            yield ("intron", arm_index, block_index,
                   replace(arm, blocks=arm.blocks[:block_index + 1]),
                   replace(arm, blocks=arm.blocks[block_index + 1:]))
        if arm_index + 1 < len(path.arms):
            yield "fusion", arm_index, len(arm.blocks) - 1, arm, path.arms[arm_index + 1]


def junction_locus_key(kind, left, right):
    return (kind, left.chrom, left.strand, left.gene_ids,
            right.chrom, right.strand, right.gene_ids)


def observed_junction_registry(paths):
    """Observed coordinates nominate alternatives; they never supply score votes."""
    registry = {}
    for path in paths:
        for kind, _, _, left, right in path_junctions(path):
            key = junction_locus_key(kind, left, right)
            registry.setdefault(key, set()).add((left.exit_position, right.entry_position))
    return {key: tuple(sorted(positions)) for key, positions in sorted(registry.items())}


def annotation_junction_registry(paths, transcripts, *, tolerance=3):
    """Nominate nearby transcript introns and gene-specific fusion arm endpoints."""
    if tolerance < 0:
        raise ValueError("annotation nomination tolerance must be nonnegative")
    donors, acceptors, introns = {}, {}, {}
    for transcript in transcripts:
        if not transcript.gene_id or transcript.strand not in ("+", "-"):
            continue
        key = transcript.chrom, transcript.strand, transcript.gene_id
        exons = sorted(set(transcript.exons), reverse=transcript.strand == "-")
        for start, end in exons:
            donors.setdefault(key, set()).add(end if transcript.strand == "+" else start)
            acceptors.setdefault(key, set()).add(start if transcript.strand == "+" else end)
        for previous, following in zip(exons, exons[1:]):
            positions = ((previous[1], following[0]) if transcript.strand == "+"
                         else (previous[0], following[1]))
            introns.setdefault(key, set()).add(positions)
    registry = {}
    for path in paths:
        for kind, _, _, left, right in path_junctions(path):
            observed = left.exit_position, right.entry_position
            positions = set()
            if kind == "intron":
                for gene in set(left.gene_ids) & set(right.gene_ids):
                    positions.update(introns.get((left.chrom, left.strand, gene), ()))
            else:
                left_positions, right_positions = {observed[0]}, {observed[1]}
                for gene in left.gene_ids:
                    left_positions.update(position for position in donors.get((left.chrom, left.strand, gene), ())
                                          if abs(position - observed[0]) <= tolerance)
                for gene in right.gene_ids:
                    right_positions.update(position for position in acceptors.get((right.chrom, right.strand, gene), ())
                                           if abs(position - observed[1]) <= tolerance)
                positions.update(product(left_positions, right_positions))
                positions.discard(observed)
            nearby = {position for position in positions
                      if max(abs(value - baseline) for value, baseline in zip(position, observed)) <= tolerance}
            registry.setdefault(junction_locus_key(kind, left, right), set()).update(nearby)
    return {key: tuple(sorted(positions)) for key, positions in sorted(registry.items())}


def correct_read_path(
    path, read_sequence, genome, *, supported_registry=None, annotation_registry=None, assembly_registry=None,
    tolerance=3, flank=35, min_anchor=12, min_identity=0.7, min_margin=6.0,
    max_candidates=128,
) -> PathCorrection:
    """Propose a complete edit only after all edges and their joint result pass.

    These are per-molecule sequence proposals, not sample-level fusion calls.
    Original evidence remains immutable. Ambiguity at any edge withholds the
    transaction rather than silently publishing a partially corrected transcript.
    """
    supported_registry = supported_registry or {}
    annotation_registry = annotation_registry or {}
    assembly_registry = assembly_registry or {}
    decisions = []
    for kind, arm_index, block_index, left, right in path_junctions(path):
        key = junction_locus_key(kind, left, right)
        decision = evaluate_junction(
            read_sequence, left, right, genome, tolerance=tolerance, flank=flank,
            min_anchor=min_anchor, min_identity=min_identity, min_margin=min_margin,
            max_candidates=max_candidates, internal_intron=kind == "intron",
            supported_positions=supported_registry.get(key, ()),
            annotation_positions=annotation_registry.get(key, ()),
            assembly_positions=assembly_registry.get(key, ()),
        )
        decisions.append(PathJunction(kind, arm_index, block_index, decision))
    decisions = tuple(decisions)
    if not decisions:
        return PathCorrection(path, None, "not_applicable", "no_junctions", decisions)
    if any(item.decision.selected_positions is None for item in decisions):
        return PathCorrection(path, None, "unresolved", "unresolved_junction", decisions)
    arms = list(path.arms)
    try:
        for item in decisions:
            arm = arms[item.arm_index]
            if item.kind == "intron":
                left, right = move_junction(
                    replace(arm, blocks=arm.blocks[:item.block_index + 1]),
                    replace(arm, blocks=arm.blocks[item.block_index + 1:]),
                    item.decision.selected_positions,
                )
                arms[item.arm_index] = replace(arm, blocks=left.blocks + right.blocks)
            else:
                arms[item.arm_index], arms[item.arm_index + 1] = move_junction(
                    arm, arms[item.arm_index + 1], item.decision.selected_positions,
                )
    except ValueError:
        return PathCorrection(path, None, "unresolved", "incompatible_joint_edits", decisions)
    proposed = replace(path, arms=tuple(arms))
    validation = []
    for kind, arm_index, block_index, left, right in path_junctions(proposed):
        decision = evaluate_junction(
            read_sequence, left, right, genome, tolerance=0, flank=flank,
            min_anchor=min_anchor, min_identity=min_identity, min_margin=min_margin,
            internal_intron=kind == "intron", max_candidates=max_candidates,
        )
        validation.append(PathJunction(kind, arm_index, block_index, decision))
    validation = tuple(validation)
    if any(item.decision.status != "retained" for item in validation):
        return PathCorrection(path, None, "unresolved", "joint_template_not_anchored",
                              decisions, validation)
    status = "correction_proposed" if proposed != path else "retained"
    return PathCorrection(path, proposed, status, "all_junctions_sequence_validated",
                          decisions, validation)
