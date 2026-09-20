"""Which molecule, on which retained path, put THIS junction exactly here.

A reviewer could not previously get from a reported model to the molecules
behind one of its junctions. The raw material was all published -- model
membership in fusion_model_evidence.jsonl, chain support and alternative
placements in fusions.tsv, per-molecule selected paths in fusion_reads.jsonl,
and per-junction decisions in junction_corrections.jsonl -- but the JOIN was
missing, and the joins a reader would reach for first are unsound.

This is a provenance-preserving export, not a support adjudicator. It adds no
score, no ranking and no filter, and it decides nothing: every row states what
was established and, where it was not, why.

The three traps it exists to avoid, all named in Codex's design review:

* Chain membership is NOT precise-coordinate support. A molecule belongs to a
  model through chains and envelopes, and may well have put its own junction
  somewhere else. `edge_relation` is established from the molecule's own
  retained path or it is `unknown`; it is never inherited.
* A molecule may retain SEVERAL paths, supporting the displayed coordinates on
  one and different coordinates on another. Every path gets its own row and
  the ambiguity survives; molecules are counted distinctly, paths are not
  counted at all.
* "Reaches both ends" is not the intersection of the two endpoint read-id
  lists: one path may reach the 5' extent and a different path of the same
  molecule reach the 3'. A both-ends witness here is always ONE path.

Deliberately NOT in this version, also on review: no aggregate per-edge
ranking of alternatives. Alternatives stay inside the decision that produced
them, with their own scores, because a score is only comparable against
another from the same decision -- never across windows, derivations,
validation passes or molecules. The correction layer's `margin` is
best-minus-runner-up within ONE local decision, so it is exported under that
name rather than as "how much the alternative lost".
"""

from __future__ import annotations

import json

from fin.fusion.analysis import _geometry

SCHEMA_VERSION = "1"

# coordinates are 0-based half-open genomic positions, exactly as the model
# table reports them; an arm's exit/entry are strand-aware (the transcript
# leaves the upstream arm at its last base and enters the downstream arm at
# its first, which on a minus arm are the low and high ends respectively)
COORDINATE_CONVENTION = "0-based half-open, strand-aware exit/entry"


def _edge_of(arms, index):
    """The junction between arms[index] and arms[index + 1], as reported."""
    left, right = arms[index], arms[index + 1]
    return (left.chrom, left.exit_position, left.strand,
            right.chrom, right.entry_position, right.strand)


def _model_edge(model, index):
    return model.breakpoints[index]


def _relation(model_edge, path_edge, tolerance):
    """How a molecule's own placement stands to the one the model displays.

    Locus identity is part of the comparison, so two junctions that happen to
    share numbers on different contigs or strands are never related. Positions
    alone are compared only once the contigs and strands already agree.
    """
    if path_edge is None:
        return "unknown"
    m_left, m_exit, m_lstrand, m_right, m_entry, m_rstrand = model_edge
    p_left, p_exit, p_lstrand, p_right, p_entry, p_rstrand = path_edge
    if (m_left, m_lstrand, m_right, m_rstrand) != (p_left, p_lstrand, p_right, p_rstrand):
        return "different_locus"
    if (m_exit, m_entry) == (p_exit, p_entry):
        return "exact"
    if max(abs(m_exit - p_exit), abs(m_entry - p_entry)) <= tolerance:
        return "different_within_tolerance"
    return "different_beyond_tolerance"


def _membership(model, read_id):
    labels = []
    for name, ids in (("unique", model.unique_read_ids),
                      ("compatible", model.compatible_read_ids),
                      ("partial", model.partial_read_ids),
                      ("nearby", model.nearby_read_ids)):
        if read_id in ids:
            labels.append(name)
    return tuple(labels)


# The chain writer's own expression, imported rather than reimplemented: a
# copied recipe would have let fin/fusion/output.py drift while this file's
# test stayed green, which is exactly what Codex found when the id was
# duplicated here. The two artifacts now cannot disagree about what an id is.
from fin.fusion.output import correction_id as _correction_id


def _correction_records(audit, path, left_arm_index):
    """Provenance for this path's junction, joined WITHOUT using coordinates.

    The bridge is the selection code's own identity: an evaluated derivation
    whose `proposed` path has the same geometry as this retained path --
    ordered topology plus every block's reference AND query coordinates, not
    merely the fusion endpoints. Selection deduplicates geometries, so several
    evaluated derivations can map to one retained path; all of them are kept.
    Inside a derivation the junction is found by kind and arm index, never by
    matching positions.
    """
    wanted = _geometry(path)
    out = []
    for evaluated_index, item in enumerate(audit.evaluated):
        if item.proposed is None or _geometry(item.proposed) != wanted:
            continue
        record = {
            # identity of THIS derivation. Several evaluated derivations can
            # converge on one retained geometry, and without these two fields
            # they are indistinguishable in the export even when their
            # internal or query geometry and their alternatives differ.
            "evaluated_index": evaluated_index,
            "correction_id": None,
            "junction_arm_index": left_arm_index if left_arm_index is not None else None,
            "validation_status": item.validation.status,
            "validation_reason": item.validation.reason,
            "min_mapping_quality": item.min_mapping_quality,
            "junction_status": None,
            "junction_reason": None,
            # best-minus-runner-up INSIDE this one decision. Not the loss of
            # every alternative against the retained path, and not comparable
            # to a margin from any other decision.
            "decision_best_runner_up_margin": None,
            "correction_time_signal_status": None,
            "observed_arm_context": None,
        }
        correction = item.correction
        if correction is None:
            record["correction_status"] = "no_correction_record"
        else:
            record["correction_status"] = correction.status
            record["correction_reason"] = correction.reason
            record["correction_id"] = _correction_id(correction.observed)
            junction = None
            if left_arm_index is not None:
                junction = next((j for j in correction.junctions
                                 if j.kind == "fusion" and j.arm_index == left_arm_index), None)
            if left_arm_index is None:
                # [P2] no arm pair matched this model edge, so there is no
                # junction correspondence. Path-level derivation identity is
                # still useful; attaching SOME fusion junction's status and
                # margin would be attaching an unrelated decision.
                record["junction_status"] = "no_junction_correspondence"
            elif junction is None:
                record["junction_status"] = "junction_absent_from_correction"
            else:
                decision = junction.decision
                record["junction_status"] = decision.status
                record["junction_reason"] = decision.reason
                record["decision_best_runner_up_margin"] = decision.margin
                # the decision's own label at correction time -- NOT the
                # later signal audit, and not post-signal validation
                record["correction_time_signal_status"] = decision.signal_status
                observed = correction.observed
                if left_arm_index + 1 < len(observed.arms):
                    record["observed_arm_context"] = list(_edge_of(observed.arms, left_arm_index))
        out.append(record)
    return out


def _reach(arms):
    """(5' extent, 3' extent) of an arm sequence, in transcript orientation."""
    first, last = arms[0], arms[-1]
    five = first.exons[0][0] if first.strand == "+" else first.exons[-1][1]
    three = last.exons[-1][1] if last.strand == "+" else last.exons[0][0]
    return five, three


def _path_reach(path):
    first, last = path.arms[0], path.arms[-1]
    five = (first.blocks[0].reference_start if first.strand == "+"
            else first.blocks[0].reference_end)
    three = (last.blocks[-1].reference_end if last.strand == "+"
             else last.blocks[-1].reference_start)
    return (first.chrom, first.strand, five), (last.chrom, last.strand, three)


def build_junction_evidence(analysis, models, *, breakpoint_tolerance=6):
    """Edge witness rows and a same-path endpoint summary, both derived only.

    One row per (model, edge, molecule, retained path, matching arm offset).
    A molecule whose lineage cannot be established still gets a row: the row
    records the membership it does hold and says the placement is unknown,
    because silence there would read as absence of support rather than absence
    of evidence about support.
    """
    if breakpoint_tolerance < 0:
        raise ValueError("breakpoint tolerance must not be negative")
    audits = {audit.read_id: audit for audit in analysis.reads}
    rows, summaries = [], []
    for model in models:
        # every molecule the model associates with, in all four categories.
        # Enumerating only unique|compatible dropped partial-only molecules --
        # which can carry exact evidence for an individual edge, the whole
        # point of this export -- and hid nearby-only molecules' differing
        # placements entirely. Labels are reported; nothing is promoted into
        # full-model support by appearing here.
        molecules = sorted(set(model.unique_read_ids) | set(model.compatible_read_ids)
                           | set(model.partial_read_ids) | set(model.nearby_read_ids))
        five_reach, three_reach = _reach(model.arms)
        spanning = []
        for read_id in molecules:
            audit = audits.get(read_id)
            labels = _membership(model, read_id)
            if audit is None or not audit.selected:
                for edge_index in range(len(model.breakpoints)):
                    rows.append(_row(model, edge_index, read_id, labels, None, None, None,
                                     "no_retained_path" if audit is not None else "no_read_record",
                                     []))
                continue
            for path_index, path in enumerate(audit.selected):
                # every arm pair that matches, so a path fitting two offsets
                # stays ambiguous instead of being resolved arbitrarily
                for edge_index, model_edge in enumerate(model.breakpoints):
                    offsets = [i for i in range(len(path.arms) - 1)
                               if _relation(model_edge, _edge_of(path.arms, i),
                                            breakpoint_tolerance) != "different_locus"]
                    if not offsets:
                        rows.append(_row(model, edge_index, read_id, labels, path_index, None,
                                         None, "no_arm_pair_at_this_locus",
                                         _correction_records(audit, path, None)))
                        continue
                    for offset in offsets:
                        rows.append(_row(
                            model, edge_index, read_id, labels, path_index, offset,
                            _edge_of(path.arms, offset), None,
                            _correction_records(audit, path, offset),
                            ambiguous=len(offsets) > 1, path=path,
                            tolerance=breakpoint_tolerance))
                head, tail = _path_reach(path)
                reaches_five = head == (model.arms[0].chrom, model.arms[0].strand, five_reach)
                reaches_three = tail == (model.arms[-1].chrom, model.arms[-1].strand, three_reach)
                if reaches_five and reaches_three:
                    spanning.append({"read_id": read_id, "path_index": path_index})
        summaries.append({
            "schema_version": SCHEMA_VERSION,
            "record": "model_endpoint_summary",
            "model_id": model.model_id,
            "coordinate_convention": COORDINATE_CONVENTION,
            "five_prime_extent": five_reach,
            "three_prime_extent": three_reach,
            # ONE path whose outer coordinates coincide with the model's, on
            # both ends. Never the intersection of the two endpoint read-id
            # lists: a molecule can reach 5' on one path and 3' on another and
            # observe neither whole envelope.
            #
            # The contract is ENDPOINT COINCIDENCE ON ONE PATH, and no more.
            # It is deliberately not "observed the model envelope": a two-arm
            # A->C path with the same outer coordinates witnesses the ends of
            # an A->B->C model while both of that model's edge rows say
            # unknown. Read the edge rows for structure; this field answers
            # only "did one molecule's single path reach both extents".
            "same_path_both_ends_witnesses": spanning,
            "witness_contract": "endpoint_coincidence_on_one_path",
            "distinct_molecules_spanning": len({w["read_id"] for w in spanning}),
            "envelope_status": model.envelope_status,
        })
    return tuple(rows), tuple(summaries)


def _row(model, edge_index, read_id, labels, path_index, arm_offset, path_edge,
         failure, corrections, *, ambiguous=False, path=None, tolerance=6):
    model_edge = _model_edge(model, edge_index)
    query_gap = None
    if path is not None and arm_offset is not None:
        left, right = path.arms[arm_offset], path.arms[arm_offset + 1]
        # geometry of the SOURCE alignment across the junction. Not an
        # insertion, and not a corrected quantity.
        query_gap = right.blocks[0].query_start - left.blocks[-1].query_end
    return {
        "schema_version": SCHEMA_VERSION,
        "record": "edge_witness",
        "model_id": model.model_id,
        "edge_index": edge_index,
        "model_edge": list(model_edge),
        "coordinate_convention": COORDINATE_CONVENTION,
        "read_id": read_id,
        # membership in the model, which is NOT evidence about this placement
        "membership_labels": list(labels),
        "path_index": path_index,
        "arm_offset": arm_offset,
        "path_edge": list(path_edge) if path_edge is not None else None,
        "edge_relation": _relation(model_edge, path_edge, tolerance),
        "arm_offset_ambiguous": ambiguous,
        "query_gap_across_junction": query_gap,
        "lineage_status": "established" if failure is None else failure,
        "correction_provenance": corrections,
    }


def write_junction_evidence(analysis, models, output_dir, *, breakpoint_tolerance=6):
    """Sidecar only: no existing artifact is read back or rewritten."""
    from pathlib import Path

    rows, summaries = build_junction_evidence(analysis, models,
                                              breakpoint_tolerance=breakpoint_tolerance)
    path = Path(output_dir) / "fusion_junction_evidence.jsonl"
    with path.open("w") as handle:
        for record in list(summaries) + list(rows):
            handle.write(json.dumps(record, sort_keys=True, separators=(",", ":"),
                                    allow_nan=False) + "\n")
    return rows, summaries
