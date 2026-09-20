"""A molecule's membership in a model is not evidence about its junction.

The export exists so a reviewer can go from a reported model to the molecules
behind one of its junctions. Every test here is one of the adversarial rows
Codex set as release blockers for it, and the theme of all of them is the
same: the joins a reader would reach for first are unsound, and the exporter
must refuse them rather than produce a confident-looking row.

* chain membership is not precise-coordinate support;
* a molecule with several retained paths is not one opinion;
* equal numbers on different contigs or strands are not the same junction;
* "reaches both ends" needs ONE path, not two read-id lists intersected.
"""
import pytest

from fin.fusion.analysis import FusionAnalysis, FusionReadAudit, FusionSettings
from fin.fusion.annotation import AnnotatedReadPath, FusionArm
from fin.fusion.junction_evidence import (
    SCHEMA_VERSION, build_junction_evidence, write_junction_evidence,
)
from fin.fusion.records import AlignmentBlock
from fin.fusion.transcript_models import ArmResolution, FusedTranscriptModel, ModelArm

RESOLUTION = ArmResolution("first", "arm_novel", (), (), (), (), (), ())


def _arm(chrom, strand, exons, query_start=0):
    ordered = list(exons) if strand == "+" else list(reversed(exons))
    blocks, q = [], query_start
    for start, end in ordered:
        blocks.append(AlignmentBlock(start, end, q, q + (end - start), end - start))
        q += end - start
    return FusionArm(chrom, strand, tuple(blocks), ()), q


def _path(read_id, spec):
    arms, q = [], 0
    for chrom, strand, exons in spec:
        arm, q = _arm(chrom, strand, exons, q)
        arms.append(arm)
    return AnnotatedReadPath(read_id, tuple(arms), "interchromosomal_candidate", "resolved")


def _model(spec, unique=(), compatible=(), partial=(), nearby=(), model_id="FTM_t"):
    arms = tuple(ModelArm(chrom, strand, tuple(exons), RESOLUTION) for chrom, strand, exons in spec)
    breakpoints = tuple(
        (l.chrom, l.exons[-1][1] if l.strand == "+" else l.exons[0][0], l.strand,
         r.chrom, r.exons[0][0] if r.strand == "+" else r.exons[-1][1], r.strand)
        for l, r in zip(arms, arms[1:]))
    return FusedTranscriptModel(
        model_id, arms, breakpoints, "chain_x", ("chain_x",), tuple(unique), tuple(compatible),
        tuple(partial), (1,) * len(breakpoints), "exact_observed", "interchromosomal_candidate",
        "CAG_x", "annotated", "annotated", nearby_read_ids=tuple(nearby))


def _analysis(paths_by_read):
    reads = tuple(FusionReadAudit(read_id, "supported_path", "", path_count=len(paths),
                                  selected=tuple(paths))
                  for read_id, paths in paths_by_read.items())
    return FusionAnalysis((), reads, {}, FusionSettings(), "provided")


SPEC = [("chr1", "+", [(100, 200)]), ("chr2", "+", [(500, 600)])]


def _rows(model, paths_by_read, tolerance=6):
    rows, _ = build_junction_evidence(_analysis(paths_by_read), [model],
                                      breakpoint_tolerance=tolerance)
    return rows


# --------------------------------------------------------------------------
# membership is not support
# --------------------------------------------------------------------------


def test_a_molecule_at_the_reported_coordinates_is_exact():
    model = _model(SPEC, unique=("r0",))
    rows = _rows(model, {"r0": [_path("r0", SPEC)]})
    assert len(rows) == 1
    assert rows[0]["edge_relation"] == "exact"
    assert rows[0]["lineage_status"] == "established"
    assert rows[0]["membership_labels"] == ["unique"]


def test_a_member_whose_own_placement_is_shifted_is_never_exact_by_inheritance():
    """The first row of Codex's kill list. The molecule belongs to the model
    through chains and envelopes; its own junction is 3 bp away."""
    model = _model(SPEC, unique=("r0",))
    shifted = [("chr1", "+", [(100, 197)]), ("chr2", "+", [(503, 600)])]
    rows = _rows(model, {"r0": [_path("r0", shifted)]})
    assert rows[0]["edge_relation"] == "different_within_tolerance"
    assert rows[0]["path_edge"] == ["chr1", 197, "+", "chr2", 503, "+"]


def test_a_placement_beyond_tolerance_says_so_rather_than_vanishing():
    model = _model(SPEC, unique=("r0",))
    far = [("chr1", "+", [(100, 150)]), ("chr2", "+", [(560, 600)])]
    assert _rows(model, {"r0": [_path("r0", far)]})[0]["edge_relation"] \
        == "different_beyond_tolerance"


def test_a_member_with_no_retained_path_gets_a_row_that_says_unknown():
    """Silence would read as absence of support rather than absence of
    evidence about support."""
    model = _model(SPEC, unique=("r0",))
    rows, _ = build_junction_evidence(_analysis({}), [model])
    assert len(rows) == 1
    assert rows[0]["edge_relation"] == "unknown"
    assert rows[0]["lineage_status"] == "no_read_record"
    assert rows[0]["path_edge"] is None
    assert rows[0]["membership_labels"] == ["unique"]


def test_membership_labels_are_reported_but_never_drive_the_relation():
    model = _model(SPEC, compatible=("r0",))
    shifted = [("chr1", "+", [(100, 198)]), ("chr2", "+", [(502, 600)])]
    row = _rows(model, {"r0": [_path("r0", shifted)]})[0]
    assert row["membership_labels"] == ["compatible"]
    assert row["edge_relation"] != "exact"


# --------------------------------------------------------------------------
# several paths, and ambiguity
# --------------------------------------------------------------------------


def test_each_retained_path_of_one_molecule_gets_its_own_row():
    """A molecule may support the displayed coordinates on one path and
    different coordinates on another. Both survive; neither is chosen."""
    model = _model(SPEC, unique=("r0",))
    other = [("chr1", "+", [(100, 196)]), ("chr2", "+", [(504, 600)])]
    rows = _rows(model, {"r0": [_path("r0", SPEC), _path("r0", other)]})
    assert len(rows) == 2
    assert [r["path_index"] for r in rows] == [0, 1]
    assert {r["edge_relation"] for r in rows} == {"exact", "different_within_tolerance"}


def test_a_path_fitting_two_offsets_keeps_both_and_says_it_is_ambiguous():
    model = _model(SPEC, unique=("r0",))
    # a three-arm path whose first AND second junctions are chr1:+ -> chr2:+
    repeated = [("chr1", "+", [(100, 200)]), ("chr2", "+", [(500, 600)]),
                ("chr2", "+", [(700, 800)])]
    rows = _rows(model, {"r0": [_path("r0", [("chr1", "+", [(100, 200)]),
                                             ("chr2", "+", [(500, 600)]),
                                             ("chr1", "+", [(900, 1000)]),
                                             ("chr2", "+", [(1500, 1600)])])]})
    assert len(rows) == 2
    assert all(r["arm_offset_ambiguous"] for r in rows)
    assert sorted(r["arm_offset"] for r in rows) == [0, 2]
    assert [r["edge_relation"] for r in sorted(rows, key=lambda r: r["arm_offset"])] \
        == ["exact", "different_beyond_tolerance"]


# --------------------------------------------------------------------------
# a junction is a locus, not two numbers
# --------------------------------------------------------------------------


def test_the_same_numbers_on_another_contig_are_not_the_same_junction():
    model = _model(SPEC, unique=("r0",))
    elsewhere = [("chr9", "+", [(100, 200)]), ("chr22", "+", [(500, 600)])]
    row = _rows(model, {"r0": [_path("r0", elsewhere)]})[0]
    assert row["lineage_status"] == "no_arm_pair_at_this_locus"
    assert row["edge_relation"] == "unknown"


def test_the_same_numbers_on_the_other_strand_are_not_the_same_junction():
    model = _model(SPEC, unique=("r0",))
    flipped = [("chr1", "-", [(100, 200)]), ("chr2", "+", [(500, 600)])]
    assert _rows(model, {"r0": [_path("r0", flipped)]})[0]["lineage_status"] \
        == "no_arm_pair_at_this_locus"


def test_exit_and_entry_are_read_strand_aware():
    """On a minus arm the transcript leaves at the LOW end and enters at the
    HIGH end; comparing raw exon bounds would mismatch every minus junction."""
    spec = [("chr1", "-", [(100, 200)]), ("chr2", "-", [(500, 600)])]
    model = _model(spec, unique=("r0",))
    assert model.breakpoints == (("chr1", 100, "-", "chr2", 600, "-"),)
    assert _rows(model, {"r0": [_path("r0", spec)]})[0]["edge_relation"] == "exact"


# --------------------------------------------------------------------------
# both ends needs one path
# --------------------------------------------------------------------------


def test_one_path_reaching_both_ends_is_a_witness():
    model = _model(SPEC, unique=("r0",))
    _, summaries = build_junction_evidence(_analysis({"r0": [_path("r0", SPEC)]}), [model])
    assert summaries[0]["same_path_both_ends_witnesses"] == [{"read_id": "r0", "path_index": 0}]
    assert summaries[0]["distinct_molecules_spanning"] == 1


def test_one_path_reaching_five_and_another_reaching_three_is_not_a_witness():
    """Codex's row: intersecting the two endpoint read-id lists would call
    this molecule spanning. It observes neither whole envelope."""
    model = _model(SPEC, unique=("r0",))
    five_only = [("chr1", "+", [(100, 200)]), ("chr2", "+", [(500, 560)])]
    three_only = [("chr1", "+", [(140, 200)]), ("chr2", "+", [(500, 600)])]
    _, summaries = build_junction_evidence(
        _analysis({"r0": [_path("r0", five_only), _path("r0", three_only)]}), [model])
    assert summaries[0]["same_path_both_ends_witnesses"] == []
    assert summaries[0]["distinct_molecules_spanning"] == 0


def test_molecules_are_counted_distinctly_and_paths_are_not_counted():
    model = _model(SPEC, unique=("r0",))
    _, summaries = build_junction_evidence(
        _analysis({"r0": [_path("r0", SPEC), _path("r0", SPEC)]}), [model])
    assert len(summaries[0]["same_path_both_ends_witnesses"]) == 2      # both paths kept
    assert summaries[0]["distinct_molecules_spanning"] == 1             # one molecule


# --------------------------------------------------------------------------
# shape, determinism, and the absence of any new judgement
# --------------------------------------------------------------------------


def test_the_query_gap_is_labelled_as_source_geometry_not_an_insertion():
    model = _model(SPEC, unique=("r0",))
    row = _rows(model, {"r0": [_path("r0", SPEC)]})[0]
    assert "query_gap_across_junction" in row
    assert not any("insert" in key for key in row)


def test_no_score_ranking_or_filter_is_introduced():
    model = _model(SPEC, unique=("r0",))
    rows = _rows(model, {"r0": [_path("r0", SPEC)]})
    banned = ("confidence", "rank", "pass", "fail", "filter", "quality_score", "probability")
    assert not any(any(word in key for word in banned) for row in rows for key in row)
    # the one margin that is exported names the decision it belongs to
    assert all("margin" not in key or key == "decision_best_runner_up_margin"
               for row in rows for record in row["correction_provenance"] for key in record)


def test_every_record_carries_the_schema_version_and_coordinate_convention():
    model = _model(SPEC, unique=("r0",))
    rows, summaries = build_junction_evidence(_analysis({"r0": [_path("r0", SPEC)]}), [model])
    for record in list(rows) + list(summaries):
        assert record["schema_version"] == SCHEMA_VERSION
        assert "0-based half-open" in record["coordinate_convention"]


def test_the_result_does_not_depend_on_molecule_order():
    model = _model(SPEC, unique=("r0", "r1"))
    forward = _rows(model, {"r0": [_path("r0", SPEC)], "r1": [_path("r1", SPEC)]})
    reverse = _rows(model, {"r1": [_path("r1", SPEC)], "r0": [_path("r0", SPEC)]})
    assert forward == reverse


def test_a_negative_tolerance_is_rejected():
    with pytest.raises(ValueError):
        build_junction_evidence(_analysis({}), [_model(SPEC)], breakpoint_tolerance=-1)


def test_the_sidecar_is_written_and_reads_back(tmp_path):
    import json

    model = _model(SPEC, unique=("r0",))
    write_junction_evidence(_analysis({"r0": [_path("r0", SPEC)]}), [model], tmp_path)
    records = [json.loads(line) for line in
               (tmp_path / "fusion_junction_evidence.jsonl").open()]
    kinds = [r["record"] for r in records]
    assert kinds == ["model_endpoint_summary", "edge_witness"]
    assert {p.name for p in tmp_path.iterdir()} == {"fusion_junction_evidence.jsonl"}


# --------------------------------------------------------------------------
# Codex's round-2 findings: the provenance path was entirely unexercised --
# every fixture above leaves audit.evaluated empty, which made the margin-name
# assertion vacuous too.
# --------------------------------------------------------------------------


def _evaluated(path, *, arm_index=0, status="retained", reason="ok", margin=3.0,
               correction_status="retained", junctions=True, observed=None):
    from fin.fusion.correction import PathCorrection, PathJunction
    from fin.fusion.junctions import JunctionDecision
    from fin.fusion.analysis import EvaluatedPath
    from fin.fusion.validation import TranscriptValidation

    observed = observed if observed is not None else path
    decision = JunctionDecision((0, 0), (0, 0), status, reason, 0, 10, (), margin)
    entries = (PathJunction("fusion", arm_index, 0, decision),) if junctions else ()
    correction = PathCorrection(observed, path, correction_status, "joint", entries)
    return EvaluatedPath(observed, path, correction,
                         TranscriptValidation("validated", "ok", 100.0), 60)


def _audit(read_id, paths, evaluated=()):
    return FusionReadAudit(read_id, "supported_path", "", path_count=len(paths),
                           evaluated=tuple(evaluated), selected=tuple(paths))


def _rows_with(model, audits, tolerance=6):
    analysis = FusionAnalysis((), tuple(audits), {}, FusionSettings(), "provided")
    rows, summaries = build_junction_evidence(analysis, [model], breakpoint_tolerance=tolerance)
    return rows, summaries


def test_a_matching_derivation_is_joined_and_carries_its_own_identity():
    model = _model(SPEC, unique=("r0",))
    path = _path("r0", SPEC)
    rows, _ = _rows_with(model, [_audit("r0", [path], [_evaluated(path, margin=7.5)])])
    provenance = rows[0]["correction_provenance"]
    assert len(provenance) == 1
    assert provenance[0]["evaluated_index"] == 0
    assert provenance[0]["correction_id"] and len(provenance[0]["correction_id"]) == 20
    assert provenance[0]["junction_arm_index"] == 0
    assert provenance[0]["decision_best_runner_up_margin"] == 7.5


def test_several_derivations_on_one_retained_geometry_stay_distinguishable():
    """Selection deduplicates geometries, so provenance is many-to-one. Two
    derivations that agree on the retained path must not collapse into one
    record, or a reader cannot tell which decision they are looking at."""
    model = _model(SPEC, unique=("r0",))
    path = _path("r0", SPEC)
    other_observed = _path("r0", [("chr1", "+", [(100, 203)]), ("chr2", "+", [(497, 600)])])
    rows, _ = _rows_with(model, [_audit("r0", [path], [
        _evaluated(path, margin=7.5),
        _evaluated(path, margin=2.0, observed=other_observed),
    ])])
    provenance = rows[0]["correction_provenance"]
    assert [p["evaluated_index"] for p in provenance] == [0, 1]
    assert len({p["correction_id"] for p in provenance}) == 2      # different observed paths
    assert [p["decision_best_runner_up_margin"] for p in provenance] == [7.5, 2.0]


def test_equal_fusion_endpoints_with_different_internal_geometry_do_not_join():
    """The geometry key is every block's reference AND query coordinates, not
    just the fusion endpoints -- otherwise an unrelated derivation that happens
    to end at the same junction would supply this row's provenance."""
    model = _model(SPEC, unique=("r0",))
    retained = _path("r0", [("chr1", "+", [(100, 150), (160, 200)]), ("chr2", "+", [(500, 600)])])
    impostor = _path("r0", [("chr1", "+", [(100, 155), (165, 200)]), ("chr2", "+", [(500, 600)])])
    rows, _ = _rows_with(model, [_audit("r0", [retained], [_evaluated(impostor)])])
    assert rows[0]["edge_relation"] == "exact"
    assert rows[0]["correction_provenance"] == []


def test_an_unmatched_edge_never_borrows_another_junction_s_decision():
    """Codex reproduced this: a chr1->chr2 model and a chr8->chr9 path gave an
    `arm_offset` null row carrying an unrelated arm-zero decision's margin."""
    model = _model(SPEC, unique=("r0",))
    elsewhere = _path("r0", [("chr8", "+", [(100, 200)]), ("chr9", "+", [(500, 600)])])
    rows, _ = _rows_with(model, [_audit("r0", [elsewhere],
                                        [_evaluated(elsewhere, margin=12.0)])])
    row = rows[0]
    assert row["lineage_status"] == "no_arm_pair_at_this_locus"
    assert row["arm_offset"] is None
    for record in row["correction_provenance"]:
        assert record["junction_status"] == "no_junction_correspondence"
        assert record["decision_best_runner_up_margin"] is None
        assert record["junction_arm_index"] is None
        assert record["correction_id"]          # path-level identity is still useful


def test_a_correction_without_the_fusion_junction_says_so():
    model = _model(SPEC, unique=("r0",))
    path = _path("r0", SPEC)
    rows, _ = _rows_with(model, [_audit("r0", [path], [_evaluated(path, junctions=False)])])
    record = rows[0]["correction_provenance"][0]
    assert record["junction_status"] == "junction_absent_from_correction"
    assert record["decision_best_runner_up_margin"] is None


def test_the_signal_field_names_the_decision_it_came_from():
    """It is the correction decision's label, not the later signal audit, and
    must not be read as post-signal validation."""
    model = _model(SPEC, unique=("r0",))
    path = _path("r0", SPEC)
    rows, _ = _rows_with(model, [_audit("r0", [path], [_evaluated(path)])])
    record = rows[0]["correction_provenance"][0]
    assert "correction_time_signal_status" in record
    assert not any(key == "junction_signal_status" for key in record)


def test_an_audit_with_no_selected_path_is_distinct_from_no_read_record():
    """The earlier `no retained path` test actually exercised `no read
    record`; these are different failures and must not share a label."""
    model = _model(SPEC, unique=("r0",))
    rows, _ = _rows_with(model, [_audit("r0", [])])
    assert rows[0]["lineage_status"] == "no_retained_path"
    rows, _ = _rows_with(model, [])
    assert rows[0]["lineage_status"] == "no_read_record"


# --------------------------------------------------------------------------
# the molecule population
# --------------------------------------------------------------------------


def test_a_partial_only_molecule_is_exported_with_its_label():
    """Partial molecules can carry exact evidence for an individual edge --
    the whole point of this export. Enumerating only unique|compatible dropped
    them silently."""
    model = _model(SPEC, unique=("r0",), partial=("p0",))
    rows, _ = _rows_with(model, [_audit("r0", [_path("r0", SPEC)]),
                                 _audit("p0", [_path("p0", SPEC)])])
    by_read = {r["read_id"]: r for r in rows}
    assert set(by_read) == {"r0", "p0"}
    assert by_read["p0"]["membership_labels"] == ["partial"]
    assert by_read["p0"]["edge_relation"] == "exact"


def test_a_nearby_only_molecule_exposes_its_differing_placement():
    model = _model(SPEC, unique=("r0",), nearby=("n0",))
    shifted = [("chr1", "+", [(100, 196)]), ("chr2", "+", [(504, 600)])]
    rows, _ = _rows_with(model, [_audit("r0", [_path("r0", SPEC)]),
                                 _audit("n0", [_path("n0", shifted)])])
    by_read = {r["read_id"]: r for r in rows}
    assert by_read["n0"]["membership_labels"] == ["nearby"]
    assert by_read["n0"]["edge_relation"] == "different_within_tolerance"
    assert by_read["n0"]["path_edge"] == ["chr1", 196, "+", "chr2", 504, "+"]


# --------------------------------------------------------------------------
# the endpoint contract, stated and tested
# --------------------------------------------------------------------------


def test_endpoint_coincidence_is_the_whole_contract_and_says_so():
    """A two-arm A->C path with the model's outer coordinates witnesses the
    ENDS of a three-arm A->B->C model while both edge rows say unknown. That
    is what the field means; the contract is named in the record so it cannot
    be read as "observed the envelope"."""
    three = [("chr1", "+", [(100, 200)]), ("chr5", "+", [(300, 400)]),
             ("chr2", "+", [(500, 600)])]
    model = _model(three, unique=("r0",))
    shortcut = _path("r0", [("chr1", "+", [(100, 200)]), ("chr2", "+", [(500, 600)])])
    rows, summaries = _rows_with(model, [_audit("r0", [shortcut])])
    assert summaries[0]["witness_contract"] == "endpoint_coincidence_on_one_path"
    assert summaries[0]["distinct_molecules_spanning"] == 1
    assert {r["edge_relation"] for r in rows} == {"unknown"}


def test_a_minus_strand_multi_exon_reach_is_computed_in_transcript_order():
    spec = [("chr1", "-", [(100, 200), (300, 400)]), ("chr2", "-", [(500, 600), (700, 800)])]
    model = _model(spec, unique=("r0",))
    # transcript order on minus runs high->low: 5' extent is chr1 400, 3' is chr2 500
    assert (model.arms[0].exons[-1][1], model.arms[-1].exons[0][0]) == (400, 500)
    _, summaries = _rows_with(model, [_audit("r0", [_path("r0", spec)])])
    assert summaries[0]["five_prime_extent"] == 400
    assert summaries[0]["three_prime_extent"] == 500
    assert summaries[0]["distinct_molecules_spanning"] == 1


# --------------------------------------------------------------------------
# Codex's round-2 residue: three gaps he verified in memory but that had no
# persistent test. An implementation that happens to pass is not a guarantee.
# --------------------------------------------------------------------------


def test_signal_pruning_does_not_shift_provenance_onto_the_wrong_derivation():
    """The dangerous case: two retained paths, the signal validator removes
    the first, and the survivor becomes path_index 0. Its provenance must
    still be its OWN derivation -- evaluated_index 1, margin 7 -- and not the
    rejected one's margin 12, which a positional association would give it.
    """
    model = _model(SPEC, unique=("r0",))
    rejected = _path("r0", [("chr1", "+", [(100, 190)]), ("chr2", "+", [(510, 600)])])
    retained = _path("r0", SPEC)
    evaluated = [_evaluated(rejected, margin=12.0), _evaluated(retained, margin=7.0)]

    before, _ = _rows_with(model, [_audit("r0", [rejected, retained], evaluated)])
    assert [r["path_index"] for r in before] == [0, 1]

    # the validator drops the rejected path; `evaluated` keeps both entries
    after, _ = _rows_with(model, [_audit("r0", [retained], evaluated)])
    assert len(after) == 1
    assert after[0]["path_index"] == 0 and after[0]["edge_relation"] == "exact"
    provenance = after[0]["correction_provenance"]
    assert [p["evaluated_index"] for p in provenance] == [1]
    assert [p["decision_best_runner_up_margin"] for p in provenance] == [7.0]


def test_the_correction_id_equals_the_one_the_chain_writer_publishes():
    """Codex round 3: my first version compared the exporter against a COPIED
    recipe, so changing fin/fusion/output.py would have left it green -- it
    pinned today's algorithm, not the join. The exporter now imports the chain
    writer's own expression, so the two cannot diverge; this test keeps an
    INDEPENDENT recipe so that shared expression is still pinned, and asserts
    the exporter really does delegate to it."""
    import hashlib
    import json as _json_mod
    from dataclasses import asdict

    model = _model(SPEC, unique=("r0",))
    path = _path("r0", SPEC)
    item = _evaluated(path)
    rows, _ = _rows_with(model, [_audit("r0", [path], [item])])

    writer_id = hashlib.sha256(_json_mod.dumps(
        asdict(item.observed), sort_keys=True, separators=(",", ":"), allow_nan=False,
    ).encode()).hexdigest()[:20]
    assert rows[0]["correction_provenance"][0]["correction_id"] == writer_id

    # and it is the writer's function that produced it, not a second copy
    import fin.fusion.junction_evidence as module
    from fin.fusion.output import correction_id

    assert module._correction_id is correction_id
    assert correction_id(item.observed) == writer_id


def test_a_query_only_geometry_difference_blocks_the_provenance_join():
    """The geometry key includes QUERY coordinates. Two derivations identical
    in every reference coordinate but differing in where the alignment sits on
    the read are different derivations."""
    model = _model(SPEC, unique=("r0",))
    retained = _path("r0", SPEC)
    shifted_on_read = AnnotatedReadPath(
        "r0",
        tuple(FusionArm(arm.chrom, arm.strand,
                        tuple(AlignmentBlock(b.reference_start, b.reference_end,
                                             b.query_start + 5, b.query_end + 5, b.aligned_bases)
                              for b in arm.blocks), arm.gene_ids)
              for arm in retained.arms),
        retained.event_class, retained.identity_status)
    rows, _ = _rows_with(model, [_audit("r0", [retained], [_evaluated(shifted_on_read)])])
    assert rows[0]["edge_relation"] == "exact"
    assert rows[0]["correction_provenance"] == []


def test_provenance_is_looked_up_at_the_matching_arm_offset_not_at_zero():
    """A three-arm path whose SECOND junction is the model's. The decision
    that belongs to this row is the one at arm_index 1."""
    model = _model(SPEC, unique=("r0",))
    path = _path("r0", [("chr5", "+", [(10, 40)]), ("chr1", "+", [(100, 200)]),
                        ("chr2", "+", [(500, 600)])])
    rows, _ = _rows_with(model, [_audit("r0", [path], [_evaluated(path, arm_index=1, margin=4.0)])])
    assert rows[0]["arm_offset"] == 1 and rows[0]["edge_relation"] == "exact"
    record = rows[0]["correction_provenance"][0]
    assert record["junction_arm_index"] == 1
    assert record["decision_best_runner_up_margin"] == 4.0

    # and the same path against a correction whose fusion junction is at arm 0
    # must NOT borrow it
    rows, _ = _rows_with(model, [_audit("r0", [path], [_evaluated(path, arm_index=0, margin=9.0)])])
    record = rows[0]["correction_provenance"][0]
    assert record["junction_status"] == "junction_absent_from_correction"
    assert record["decision_best_runner_up_margin"] is None
