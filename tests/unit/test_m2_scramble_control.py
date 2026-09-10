"""The signal-free M2 control: what it must and must not change.

The control answers whether the tie-break's benefit is the SIGNAL CONTENT or
merely the act of concentrating a tied read onto one hypothesis. For that
comparison to mean anything, the permutation must leave every other input to
the decision untouched -- otherwise an arm difference could come from a changed
read population, a changed abstention, or a changed threshold crossing instead
of from the removed signal.

The decision is STRUCTURAL. `Assigner.assign` reduces the scored candidates to
per-structure minima and thresholds the margin between those minima, so the
quantity that must be preserved is the multiset of structure minima -- NOT the
multiset of candidate values. These tests are written against that, because a
candidate-level invariant can hold while the gate's decision flips.
"""
import hashlib

import pytest

from fin.pipeline.assignment import _scramble_contrast, structure_key


class FakeChain:
    def __init__(self, introns):
        self.introns = tuple(introns)


class FakeCand:
    def __init__(self, cid, introns=(), source="gtf"):
        self.candidate_id = cid
        self.intron_chain = FakeChain(introns)
        self.source = source


def gate_view(nlls, cands):
    """Reproduce what the diff-cover gate reads: per-structure minima, the
    winning structure, and the margin it thresholds."""
    by_group = {}
    for j, v in nlls.items():
        k = structure_key(cands[j])
        by_group[k] = v if k not in by_group else min(by_group[k], v)
    ordered = sorted(by_group.items(), key=lambda kv: kv[1])
    margin = ordered[1][1] - ordered[0][1] if len(ordered) > 1 else None
    return by_group, ordered[0][0], margin


# Two candidates sharing a chain (an annotated model and its read-extent
# sibling) plus a distinct structure -- the shape that broke the first version
# of this control.
DUP = [
    FakeCand("A", [(100, 200)]),
    FakeCand("A_extent", [(100, 200)]),
    FakeCand("B", [(100, 250)]),
]
SINGLE = [FakeCand(c, [(100, 100 * i + 200)]) for i, c in enumerate("ABCD")]


class TestOff:
    def test_zero_seed_returns_the_same_object(self):
        nlls = {0: 1.0, 1: 2.0, 2: 3.0}
        assert _scramble_contrast(nlls, 0, "r1", SINGLE) is nlls

    def test_negative_seed_is_off(self):
        nlls = {0: 1.0, 1: 2.0}
        assert _scramble_contrast(nlls, -7, "r1", SINGLE) is nlls

    def test_missing_candidates_is_off_rather_than_mis_grouped(self):
        # Without the candidates the control cannot know the structures, and
        # guessing would be exactly the confound it exists to avoid.
        nlls = {0: 1.0, 1: 2.0}
        assert _scramble_contrast(nlls, 11, "r1", None) is nlls

    def test_single_structure_is_untouched(self):
        # The gate abstains on a one-structure contrast either way, and there
        # is no structural choice to randomize.
        nlls = {0: 1.0, 1: 5.0}
        cands = [FakeCand("A", [(1, 2)]), FakeCand("A2", [(1, 2)])]
        assert _scramble_contrast(nlls, 11, "r1", cands) is nlls


class TestTheGateSeesAnIdenticalDecisionProblem:
    """The invariant that makes this a control rather than a second lever."""

    @pytest.mark.parametrize("seed", range(1, 40))
    def test_structure_minima_multiset_is_preserved_with_duplicates(self, seed):
        nlls = {0: 1.0, 1: 1.0, 2: 3.0}
        before, _, m_before = gate_view(nlls, DUP)
        out = _scramble_contrast(nlls, seed, "read-A", DUP)
        after, _, m_after = gate_view(out, DUP)
        assert sorted(after.values()) == sorted(before.values())
        assert m_after == m_before

    @pytest.mark.parametrize("seed", range(1, 40))
    def test_the_thresholded_margin_survives_every_seed(self, seed):
        # The exact failure the first version had: candidate multiset intact,
        # structure margin collapsing 2.0 -> 0.0, flipping hard-assign to
        # abstain and changing the covered-vote prior for OTHER reads.
        nlls = {0: 1.0, 1: 1.0, 2: 3.0}
        out = _scramble_contrast(nlls, seed, "read-A", DUP)
        assert gate_view(out, DUP)[2] == pytest.approx(2.0)

    @pytest.mark.parametrize("seed", range(1, 40))
    def test_hard_versus_abstain_classification_never_flips(self, seed):
        for values, thr in (
            ({0: 1.0, 1: 1.0, 2: 3.0}, 1.0),     # comfortably hard
            ({0: 1.0, 1: 1.0, 2: 1.4}, 1.0),     # comfortably abstaining
            ({0: 1.0, 1: 2.0, 2: 2.0}, 1.0),     # exactly at the threshold
        ):
            before = gate_view(values, DUP)[2] >= thr
            out = _scramble_contrast(values, seed, "read-A", DUP)
            assert (gate_view(out, DUP)[2] >= thr) is before

    @pytest.mark.parametrize("seed", range(1, 25))
    def test_the_scored_column_set_is_preserved(self, seed):
        nlls = {0: 1.0, 1: 2.0, 2: 3.0}
        assert set(_scramble_contrast(nlls, seed, "r", DUP)) == set(nlls)

    @pytest.mark.parametrize("seed", range(1, 25))
    def test_within_structure_order_survives(self, seed):
        # The extent sibling must not overtake its own annotated model just
        # because the structure moved; only the structure's rank changes.
        nlls = {0: 1.0, 1: 2.5, 2: 3.0}
        out = _scramble_contrast(nlls, seed, "r", DUP)
        assert out[0] < out[1]
        assert out[1] - out[0] == pytest.approx(1.5)

    @pytest.mark.parametrize("seed", range(1, 25))
    def test_singleton_structures_preserve_the_value_multiset_too(self, seed):
        # With no duplicates, structure minima ARE the candidate values, so the
        # older candidate-level invariant must still hold.
        nlls = {0: 1.0, 1: 2.0, 2: 3.0, 3: 9.0}
        out = _scramble_contrast(nlls, seed, "r", SINGLE)
        assert sorted(out.values()) == [1.0, 2.0, 3.0, 9.0]


class TestExactFloatingPointPreservation:
    """A margin preserved only to ~1e-16 is not preserved.

    The gate's test is ``margin >= margin_thr`` with no tolerance, so a
    one-ulp residue at the threshold flips a hard assignment to an abstention
    -- silently, and only for the reads sitting exactly on the boundary, which
    are the ones the control is least able to afford losing. Reconstructing a
    group's minimum as ``old_min + (new_min - old_min)`` does exactly that;
    these are the regressions for it, so the assertions are ``==``, never
    ``approx``.
    """

    # Values that are not exactly representable in binary floating point, at
    # magnitudes far enough apart that a shift loses low bits.
    TRIPLE = [
        FakeCand("a", [(100, 200)]),
        FakeCand("b", [(100, 201)]),
        FakeCand("c", [(100, 202)]),
    ]

    @pytest.mark.parametrize("seed", range(1, 40))
    def test_the_reported_cancellation_case(self, seed):
        # 0.1 / 1.1 / 100.0 with a structure margin of exactly 1.0: the shift
        # form produced 0.9999999999999942 here and dropped the read below the
        # production threshold.
        nlls = {0: 0.1, 1: 1.1, 2: 100.0}
        out = _scramble_contrast(nlls, seed, "r", self.TRIPLE)
        assert gate_view(out, self.TRIPLE)[2] == 1.0
        assert (gate_view(out, self.TRIPLE)[2] >= 1.0) is True

    @pytest.mark.parametrize("seed", range(1, 40))
    @pytest.mark.parametrize(
        "nlls",
        [
            {0: 0.1, 1: 0.2, 2: 0.3},
            {0: 1e-8, 1: 1.0 + 1e-8, 2: 1e6},
            {0: 0.30000000000000004, 1: 1.3000000000000003, 2: 7.7},
            {0: -0.1, 1: 0.9, 2: 1234.5678},
            {0: 3.3, 1: 3.3000000000000007, 2: 1e9},
        ],
    )
    def test_structure_minima_are_preserved_bit_for_bit(self, seed, nlls):
        before = sorted(gate_view(nlls, self.TRIPLE)[0].values())
        out = _scramble_contrast(nlls, seed, "r", self.TRIPLE)
        assert sorted(gate_view(out, self.TRIPLE)[0].values()) == before

    def test_the_known_limitation_is_pinned_not_claimed_away(self):
        # A tiny positive offset can round back onto a large-magnitude base, so
        # a member that was STRICTLY worse than its group minimum can become
        # equal to it. The gate reads minima only and is unaffected, but the
        # support gate collects every column at `v <= best`, so it can see one
        # extra winner. Pinned as a regression: if a future change makes this
        # exact or makes it worse, this test says so rather than the limitation
        # quietly changing shape.
        nlls = {0: -0.1, 1: -0.09999999999999999, 2: -100.0}
        assert nlls[1] > nlls[0]                      # strictly worse to start
        out = _scramble_contrast(nlls, 2, "r", DUP)
        assert out[0] == out[1] == -100.0             # collapsed to equal
        # ... while the quantity the gate actually reads is still exact.
        assert sorted(gate_view(out, DUP)[0].values()) == sorted(
            gate_view(nlls, DUP)[0].values())
        assert gate_view(out, DUP)[2] == gate_view(nlls, DUP)[2]

    @pytest.mark.parametrize("seed", range(1, 40))
    def test_a_minimum_member_lands_exactly_on_its_assigned_value(self, seed):
        # With duplicates the group minimum is shared by two candidates; both
        # must land on the assigned value exactly, or the `min` reduction picks
        # up the residue.
        nlls = {0: 0.1, 1: 0.1, 2: 1.1}
        out = _scramble_contrast(nlls, seed, "r", DUP)
        assert out[0] == out[1]
        assert sorted({out[0], out[2]}) == [0.1, 1.1]


class TestItActuallyRandomizesTheWinner:
    def test_the_winning_structure_does_move(self):
        nlls = {j: float(j) for j in range(4)}
        winners = {
            gate_view(_scramble_contrast(nlls, s, "read-A", SINGLE), SINGLE)[1]
            for s in range(1, 40)
        }
        assert len(winners) > 1

    def test_it_moves_with_duplicates_present(self):
        nlls = {0: 1.0, 1: 1.0, 2: 3.0}
        winners = {
            gate_view(_scramble_contrast(nlls, s, "read-A", DUP), DUP)[1]
            for s in range(1, 40)
        }
        assert len(winners) == 2

    def test_the_draw_depends_on_the_read(self):
        nlls = {j: float(j) for j in range(4)}
        winners = {
            gate_view(_scramble_contrast(nlls, 11, f"read-{i}", SINGLE),
                      SINGLE)[1]
            for i in range(40)
        }
        assert len(winners) > 1

    def test_members_of_one_structure_move_together(self):
        # Hard assignment gives mass to EVERY candidate of the winning
        # structure, so the two members must never end up in different
        # structures' slots.
        nlls = {0: 1.0, 1: 1.0, 2: 3.0}
        for s in range(1, 40):
            out = _scramble_contrast(nlls, s, "r", DUP)
            assert out[0] == out[1]


class TestDeterminism:
    def test_same_inputs_give_the_same_permutation(self):
        nlls = {0: 1.0, 1: 2.0, 2: 3.0}
        a = _scramble_contrast(nlls, 11, "read-A", SINGLE)
        b = _scramble_contrast(nlls, 11, "read-A", SINGLE)
        assert a == b

    def test_it_does_not_use_pythons_salted_hash(self):
        # PYTHONHASHSEED randomizes str.__hash__ per process, which would make
        # the six grid workers disagree and a rerun irreproducible.
        nlls = {j: float(j) for j in range(4)}
        winner = gate_view(
            _scramble_contrast(nlls, 11, "read-A", SINGLE), SINGLE)[1]
        keys = [structure_key(c) for c in SINGLE]
        expected = min(keys, key=lambda k: (hashlib.blake2b(
            f"11\tread-A\t{k}".encode(), digest_size=8).digest(), repr(k)))
        assert winner == expected

    def test_identity_is_the_structure_not_the_column(self):
        # Columns are per-interval packing positions. Keying on them would make
        # the draw depend on how an interval happened to be packed, so two
        # comparable runs could disagree.
        nlls = {0: 1.0, 1: 2.0}
        a = [FakeCand("X", [(1, 2)]), FakeCand("Y", [(3, 4)])]
        b = [FakeCand("Y", [(3, 4)]), FakeCand("X", [(1, 2)])]
        wa = gate_view(_scramble_contrast(nlls, 11, "r", a), a)[1]
        wb = gate_view(_scramble_contrast(nlls, 11, "r", b), b)[1]
        assert wa == wb


class TestStructureKeyIsShared:
    def test_mono_exon_candidates_are_their_own_structure(self):
        assert structure_key(FakeCand("m1")) != structure_key(FakeCand("m2"))

    def test_a_fusion_never_shares_a_structure_with_a_linear_transcript(self):
        # A fusion's intron_chain holds arm A only; sharing it would merge two
        # different hypotheses into one.
        chain = [(100, 200)]
        fus = FakeCand("f1", chain, source="fusion")
        lin = FakeCand("t1", chain)
        assert structure_key(fus) != structure_key(lin)

    def test_the_assigner_uses_this_one_definition(self):
        # If the gate re-derived its own grouping, the control could be
        # controlling for a partition the gate does not use.
        import inspect

        from fin.pipeline.assignment import Assigner

        src = inspect.getsource(Assigner.assign)
        assert "structure_key(c) for c in cand_list" in src
        assert "def _structure_key" not in src


class TestWiring:
    def test_the_production_default_is_off(self):
        from fin.pipeline.config import PipelineConfig

        assert PipelineConfig.m2_scramble_seed == 0

    def test_both_contrast_paths_apply_the_control(self):
        # A control covering only the batched path would make the arm a
        # mixture of scrambled and live decisions.
        import inspect

        from fin.pipeline import assignment

        src = inspect.getsource(assignment.tie_nll)
        assert "m2_scramble_seed" in src
        assert src.count("_scramble_contrast(nlls, scramble, rid, cand_list)") == 2

    def test_a_positive_seed_is_rejected_outside_m2_em(self):
        # argmax resolves ties through m2_resolve_tie, which this control does
        # not touch; a positive seed there would be a partial scramble.
        import inspect

        from fin import cli

        src = inspect.getsource(cli)
        assert 'm2_scramble_seed > 0 and quant_mode != "m2_em"' in src

    def test_a_positive_seed_is_rejected_when_the_metric_is_off(self):
        # `tie_nll` returns before scoring, so the seed would scramble nothing
        # and the arm would silently be the no-signal arm under another name.
        import inspect

        from fin import cli

        assert ('m2_scramble_seed > 0 and m2_metric == "off"'
                in inspect.getsource(cli))

    def test_a_positive_seed_is_rejected_without_the_diff_cover_gate(self):
        # The preserved invariant is the GATE's margin. With the gate off the
        # skeleton is rescaled by candidate-level extrema instead, which a
        # structure permutation does not preserve.
        import inspect

        from fin import cli

        assert ("m2_scramble_seed > 0 and not m2_diff_cover_gate"
                in inspect.getsource(cli))
