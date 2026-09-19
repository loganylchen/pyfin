"""Reported models that describe one junction geometry must be linked.

The whole K562 library published three BAG6>SLC44A4 rows -- 4 / 2 / 2
molecules, three different `alternative_group_id`s, `support_compatible` 0 on
all three -- two at byte-identical breakpoints and one 6 bp away at the
acceptor. Nothing linked them. The same shape appears at
chrX:49590731:+>chrX:49600723:+ (26 / 2) and chrX:135882535:->chrX:135864028:-
(7 / 1); three further duplicate-breakpoint pairs were already in one
alternative group because their molecules were mutually compatible.

The grouping here is geometric. It does not claim the members are isoforms of
one biological event, does not change which models are reported, and does not
move a molecule. The bound is on the breakpoints a model REPORTS: the two
BAG6 `...784` models retain an alternative at 31865787, so the family's
represented placements span 9 bp even though its representatives span 6.
"""
import pytest

from fin.fusion.event_families import (
    ALGORITHM_VERSION, assign_event_families, families_by_model,
)


class _Model:
    """Only the three attributes the partition reads."""

    def __init__(self, model_id, breakpoints, unique=(), compatible=()):
        self.model_id = model_id
        self.breakpoints = tuple(breakpoints)
        self.unique_read_ids = tuple(unique)
        self.compatible_read_ids = tuple(compatible)


def _bp(exit_position, entry_position, *, left="chr6", right="chr6", strand="-"):
    return (left, exit_position, strand, right, entry_position, strand)


# --------------------------------------------------------------------------
# the reported defect
# --------------------------------------------------------------------------


def test_the_three_bag6_rows_become_one_family_of_eight_molecules():
    models = [
        _Model("FTM_0f0c0ea97bb8ad0c", [_bp(31651655, 31865784)], unique=("a", "b", "c", "d")),
        _Model("FTM_3b8a7bfc9d14ce5d", [_bp(31651655, 31865784)], unique=("e", "f")),
        _Model("FTM_319b1b02a378c542", [_bp(31651655, 31865778)], unique=("g", "h")),
    ]
    families = assign_event_families(models, breakpoint_tolerance=6)
    assert len(families) == 1
    assert families[0].model_count == 3
    assert families[0].support == 8
    assert families[0].support_read_ids == ("a", "b", "c", "d", "e", "f", "g", "h")


def test_two_unlinked_siblings_at_identical_breakpoints_join():
    """chrX:49590731:+>chrX:49600723:+ -- 26 and 2 molecules, two groups."""
    models = [
        _Model("FTM_2463ed97f9f97dc8", [_bp(49590731, 49600723, left="chrX", right="chrX", strand="+")],
               unique=tuple(f"r{i}" for i in range(26))),
        _Model("FTM_7f7ee48a65ea8713", [_bp(49590731, 49600723, left="chrX", right="chrX", strand="+")],
               unique=("s0", "s1")),
    ]
    families = assign_event_families(models)
    assert len(families) == 1 and families[0].support == 28


def test_support_is_a_union_of_molecules_never_a_sum():
    shared = _Model("m1", [_bp(100, 900)], unique=("u1",), compatible=("c1", "c2"))
    other = _Model("m2", [_bp(100, 900)], unique=("u2",), compatible=("c1", "c2"))
    family = assign_event_families([shared, other])[0]
    assert family.support == 4                       # not 3 + 3
    assert family.support_read_ids == ("c1", "c2", "u1", "u2")


def test_partial_and_nearby_evidence_is_not_family_support():
    model = _Model("m1", [_bp(100, 900)], unique=("u1",))
    model.partial_read_ids = ("p1", "p2")            # present on the real model
    model.nearby_read_ids = ("n1",)
    assert assign_event_families([model])[0].support_read_ids == ("u1",)


# --------------------------------------------------------------------------
# what must NOT be grouped
# --------------------------------------------------------------------------


def test_a_different_contig_is_a_different_family():
    models = [_Model("m1", [_bp(100, 900, left="chr6", right="chr6")]),
              _Model("m2", [_bp(100, 900, left="chr6", right="chr7")])]
    assert len(assign_event_families(models)) == 2


def test_a_different_strand_is_a_different_family():
    models = [_Model("m1", [_bp(100, 900, strand="+")]),
              _Model("m2", [_bp(100, 900, strand="-")])]
    assert len(assign_event_families(models)) == 2


def test_a_two_arm_model_never_joins_a_three_arm_family():
    """An A->B model must not join A->B->C because one edge agrees."""
    two = _Model("m2arm", [_bp(100, 900)])
    three = _Model("m3arm", [_bp(100, 900), _bp(1500, 2500)])
    families = assign_event_families([two, three])
    assert len(families) == 2
    assert {f.model_count for f in families} == {1}


def test_exit_and_entry_are_compared_independently():
    """A donor 1 bp off and an acceptor 9 bp off must not average to 5."""
    models = [_Model("m1", [_bp(100, 900)]), _Model("m2", [_bp(101, 909)])]
    assert len(assign_event_families(models, breakpoint_tolerance=6)) == 2


def test_a_junction_beyond_tolerance_stays_separate():
    models = [_Model("m1", [_bp(100, 900)]), _Model("m2", [_bp(100, 907)])]
    assert len(assign_event_families(models, breakpoint_tolerance=6)) == 2


def test_gene_assignment_is_not_part_of_the_key():
    """Three of the six duplicate breakpoint strings in the reference run are
    at unannotated loci; a gene-based key would not group them at all."""
    a = _Model("m1", [_bp(100, 900)])
    b = _Model("m2", [_bp(100, 900)])
    a.gene_path, b.gene_path = (("BAG6",), ("SLC44A4",)), ((), ())
    assert len(assign_event_families([a, b])) == 1


# --------------------------------------------------------------------------
# the partition rule itself
# --------------------------------------------------------------------------


def test_bounded_diameter_not_single_linkage():
    """0, 6, 12 at tolerance 6 must NOT become one family: the extremes are
    12 apart. The specified rule gives {0,6} and {12}."""
    models = [_Model("m0", [_bp(100, 900)]),
              _Model("m6", [_bp(100, 906)]),
              _Model("m12", [_bp(100, 912)])]
    families = assign_event_families(models, breakpoint_tolerance=6)
    assert len(families) == 2
    sizes = {f.model_count for f in families}
    assert sizes == {2, 1}
    pair = next(f for f in families if f.model_count == 2)
    assert pair.model_ids == ("m0", "m6")            # lexicographic vector order


def test_every_pair_inside_a_family_is_within_tolerance():
    import random

    generator = random.Random(2026)
    for _ in range(200):
        models = [_Model(f"m{i}", [_bp(100, 900 + generator.randrange(0, 40))])
                  for i in range(generator.randrange(1, 9))]
        for family in assign_event_families(models, breakpoint_tolerance=6):
            vectors = family.breakpoint_vectors
            for first in vectors:
                for second in vectors:
                    assert max(abs(x - y) for x, y in zip(first, second)) <= 6


def test_adding_another_model_at_an_existing_geometry_cannot_repartition():
    """The partition runs over DISTINCT vectors, so publishing a second isoform
    at a geometry already present must not move anything."""
    base = [_Model("m0", [_bp(100, 900)]), _Model("m6", [_bp(100, 906)]),
            _Model("m12", [_bp(100, 912)])]
    before = [f.breakpoint_vectors for f in assign_event_families(base)]
    base.append(_Model("m6b", [_bp(100, 906)]))
    after = [f.breakpoint_vectors for f in assign_event_families(base)]
    assert before == after


def test_the_result_does_not_depend_on_input_order():
    import itertools

    models = [_Model("m0", [_bp(100, 900)]), _Model("m6", [_bp(100, 906)]),
              _Model("m12", [_bp(100, 912)]), _Model("mx", [_bp(100, 903)])]
    baseline = None
    for order in itertools.permutations(models):
        got = [(f.family_id, f.model_ids, f.breakpoint_vectors)
               for f in assign_event_families(list(order))]
        baseline = baseline if baseline is not None else got
        assert got == baseline


def test_zero_tolerance_groups_only_identical_geometries():
    models = [_Model("m1", [_bp(100, 900)]), _Model("m2", [_bp(100, 900)]),
              _Model("m3", [_bp(100, 901)])]
    families = assign_event_families(models, breakpoint_tolerance=0)
    assert sorted(f.model_count for f in families) == [1, 2]


def test_a_negative_tolerance_is_rejected():
    with pytest.raises(ValueError):
        assign_event_families([_Model("m1", [_bp(100, 900)])], breakpoint_tolerance=-1)


def test_a_model_without_a_junction_is_rejected_not_silently_grouped():
    with pytest.raises(ValueError, match="no junction"):
        assign_event_families([_Model("m1", [])])


# --------------------------------------------------------------------------
# identity
# --------------------------------------------------------------------------


def test_the_family_id_is_a_digest_of_geometry_only():
    """Support, gene labels and model ids must not redefine family identity --
    otherwise a family's id would change when a molecule is reassigned."""
    lean = [_Model("m1", [_bp(100, 900)], unique=("a",))]
    rich = [_Model("ZZZ_other_id", [_bp(100, 900)], unique=("a", "b", "c"))]
    assert assign_event_families(lean)[0].family_id == assign_event_families(rich)[0].family_id


def test_the_id_changes_with_the_tolerance_it_was_computed_under():
    one = assign_event_families([_Model("m1", [_bp(100, 900)])], breakpoint_tolerance=6)
    two = assign_event_families([_Model("m1", [_bp(100, 900)])], breakpoint_tolerance=4)
    assert one[0].family_id != two[0].family_id


def test_the_id_carries_the_algorithm_version():
    import fin.fusion.event_families as module

    before = assign_event_families([_Model("m1", [_bp(100, 900)])])[0].family_id
    original = module.ALGORITHM_VERSION
    try:
        module.ALGORITHM_VERSION = original + "x"
        after = assign_event_families([_Model("m1", [_bp(100, 900)])])[0].family_id
    finally:
        module.ALGORITHM_VERSION = original
    assert before != after and ALGORITHM_VERSION == original


def test_every_model_belongs_to_exactly_one_family():
    models = [_Model(f"m{i}", [_bp(100, 900 + 20 * i)]) for i in range(5)]
    families = assign_event_families(models)
    by_model = families_by_model(families)
    assert set(by_model) == {m.model_id for m in models}
    assert sum(f.model_count for f in families) == len(models)


# --------------------------------------------------------------------------
# the adversarial review's two findings
# --------------------------------------------------------------------------


def test_a_vector_eligible_for_two_families_joins_the_nearer_one():
    """Competition needs BOTH families to exist when the candidate is placed.
    Vectors are visited in sorted order, so varying only one coordinate puts
    the candidate last but leaves one family -- the first version of this test
    passed under "first eligible" and "last eligible" alike and pinned nothing.

    Two seeds 10 apart (so they cannot merge), then a candidate 6 from the
    first and 4 from the second. The nearer one must win, which rejects a
    "first eligible" rule.
    """
    models = [_Model("seed_a", [_bp(100, 910)]), _Model("seed_b", [_bp(101, 900)]),
              _Model("cand", [_bp(102, 904)])]
    families = assign_event_families(models, breakpoint_tolerance=6)
    assert len(families) == 2
    by_model = families_by_model(families)
    assert by_model["cand"].family_id == by_model["seed_b"].family_id
    assert by_model["cand"].family_id != by_model["seed_a"].family_id


def test_an_equidistant_vector_breaks_the_tie_on_the_smaller_vector():
    """Two seeds 12 apart, a candidate exactly 6 from each. The family holding
    the lexicographically smaller vector must win, which rejects a "last
    eligible" rule -- and makes the outcome independent of creation order."""
    models = [_Model("seed_low", [_bp(100, 912)]), _Model("seed_high", [_bp(101, 900)]),
              _Model("cand", [_bp(102, 906)])]
    families = assign_event_families(models, breakpoint_tolerance=6)
    assert len(families) == 2
    by_model = families_by_model(families)
    assert by_model["cand"].family_id == by_model["seed_low"].family_id
    assert by_model["cand"].family_id != by_model["seed_high"].family_id


def test_the_ranking_distance_is_the_worst_member_not_the_best():
    """Both earlier competition tests compared SINGLETON families, so ranking
    by max-over-members and ranking by min-over-members were indistinguishable.
    Codex found both mutants surviving the whole file.

    Family A ends up holding (100,910) and (102,906). The candidate (103,904)
    is 6 from A's far member but only 2 from its near one, and 4 from B. Under
    the worst-member rule B wins; under a best-member rule A would.
    """
    models = [_Model("a1", [_bp(100, 910)]), _Model("b1", [_bp(101, 900)]),
              _Model("a2", [_bp(102, 906)]), _Model("cand", [_bp(103, 904)])]
    by_model = families_by_model(assign_event_families(models, breakpoint_tolerance=6))
    assert by_model["a1"].family_id == by_model["a2"].family_id      # A formed first
    assert by_model["cand"].family_id == by_model["b1"].family_id
    assert by_model["cand"].family_id != by_model["a1"].family_id


def test_the_ranking_distance_examines_every_member_not_only_the_first():
    """Codex round 4: the worst-member fixture above happens to put A's
    farthest member first, so ranking on `members[0]` alone survived it.

    Here A holds (100,910) then (102,912), and the candidate (103,907) is 3
    from the first and 5 from the second, against 4 to B. Ranking on the whole
    membership sends it to B; ranking on the first member alone would keep it
    in A.
    """
    models = [_Model("a1", [_bp(100, 910)]), _Model("b1", [_bp(101, 903)]),
              _Model("a2", [_bp(102, 912)]), _Model("cand", [_bp(103, 907)])]
    by_model = families_by_model(assign_event_families(models, breakpoint_tolerance=6))
    assert by_model["a1"].family_id == by_model["a2"].family_id
    assert by_model["cand"].family_id == by_model["b1"].family_id
    assert by_model["cand"].family_id != by_model["a1"].family_id


def test_the_tie_key_is_the_family_s_smallest_vector_not_its_largest():
    """A holds (100,912) and (102,907); the candidate (103,906) is exactly 6
    from A's worst member and 6 from B, so the tie key decides. Smallest vector
    (100,912) < (101,900) sends it to A; a largest-vector key would compare
    (102,907) against (101,900) and send it to B."""
    models = [_Model("a1", [_bp(100, 912)]), _Model("b1", [_bp(101, 900)]),
              _Model("a2", [_bp(102, 907)]), _Model("cand", [_bp(103, 906)])]
    by_model = families_by_model(assign_event_families(models, breakpoint_tolerance=6))
    assert by_model["a1"].family_id == by_model["a2"].family_id
    assert by_model["cand"].family_id == by_model["a1"].family_id
    assert by_model["cand"].family_id != by_model["b1"].family_id


def test_the_writer_survives_a_generator_of_models(tmp_path):
    """assign_event_families consumes its argument. Handing the writer a
    generator used to publish a header and no rows, silently -- the family
    computation had eaten the models before the write loop saw them."""
    import csv

    from fin.fusion.transcript_output import write_transcript_outputs
    from tests.unit.test_fusion_transcript_models_unit import (
        GENOME, INDEX, _analysis, _path,
    )
    from fin.fusion.transcript_models import build_fused_transcript_models

    paths = [_path(f"r{i}", [("A", "+", [(150, 200), (300, 400)], "A"),
                             ("B", "+", [(500, 600)], "B")]) for i in range(2)]
    analysis = _analysis({p.read_id: [p] for p in paths}, [p.read_id for p in paths])
    models, collapsed = build_fused_transcript_models(analysis, INDEX, GENOME)
    assert len(models) == 1

    write_transcript_outputs(iter(models), collapsed, GENOME, tmp_path)
    rows = list(csv.DictReader((tmp_path / "fusion_models.tsv").open(), delimiter="\t"))
    assert len(rows) == 1
    assert (tmp_path / "fusion_models.fa").read_text().startswith(">")
