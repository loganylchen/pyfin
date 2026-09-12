"""Fixed-window sequence contrasts cannot win by dropping a junction arm."""

import pytest

from fin.fusion.alignment import align_junction_window


LEFT = "ACGTCAGTACCGATGCTAGC"
RIGHT = "TTCGAGCATGACCTAGGTCA"


def test_exact_junction_has_both_anchors():
    result = align_junction_window(LEFT + RIGHT, LEFT + RIGHT, len(LEFT))
    assert result.status == "scored"
    assert result.score == 80
    assert result.left_matches == result.left_aligned == 20
    assert result.right_matches == result.right_aligned == 20


def test_target_padding_is_not_counted_as_read_evidence():
    result = align_junction_window(LEFT + RIGHT, "GGGGG" + LEFT + RIGHT + "AAAAA", 25)
    assert result.status == "scored"
    assert result.target_start == 5
    assert result.target_end == 45
    assert result.score == 80


def test_wrong_junction_has_lower_score_than_read_supported_alternative():
    query = LEFT + RIGHT
    correct = align_junction_window(query, query, 20)
    shifted = align_junction_window(query, LEFT + "AA" + RIGHT, 22)
    assert shifted.status == "scored"
    assert shifted.score < correct.score
    assert shifted.deleted_bases == 2


def test_error_tolerant_anchors_do_not_require_perfect_matching():
    query = LEFT[:5] + "T" + LEFT[6:] + RIGHT[:5] + "T" + RIGHT[6:]
    result = align_junction_window(query, LEFT + RIGHT, 20)
    assert result.status == "scored"
    assert result.left_matches >= 18
    assert result.right_matches >= 18


def test_single_arm_mapping_does_not_validate_fusion():
    result = align_junction_window(LEFT, LEFT + RIGHT, 20)
    assert result.status == "unscorable"
    assert result.reason == "missing_junction_anchor"


def test_query_clipping_cannot_hide_an_internal_adapter():
    result = align_junction_window(LEFT + "A" * 30 + RIGHT, LEFT + RIGHT, 20)
    assert result.status == "unscorable"
    assert result.reason == "large_indel"


@pytest.mark.parametrize("query,target,junction,reason", [
    ("", LEFT + RIGHT, 20, "invalid_window"),
    (LEFT + RIGHT, LEFT + RIGHT, 0, "invalid_window"),
    (LEFT + "N" + RIGHT, LEFT + RIGHT, 20, "ambiguous_sequence"),
    ("A" * 300, "A" * 300, 150, "window_limit"),
])
def test_invalid_or_unbounded_windows_abstain(query, target, junction, reason):
    result = align_junction_window(query, target, junction)
    assert result.status == "unscorable"
    assert result.reason == reason


def test_insertion_penalty_matches_declared_affine_cost():
    exact = align_junction_window(LEFT + RIGHT, LEFT + RIGHT, 20)
    inserted = align_junction_window(LEFT + "AAA" + RIGHT, LEFT + RIGHT, 20)
    assert inserted.status == "scored"
    assert inserted.inserted_bases == 3
    assert exact.score - inserted.score == 4 + 3 * 2


def test_target_longer_than_query_cannot_obtain_extra_match_reward():
    result = align_junction_window(LEFT + RIGHT, "AAAAA" + LEFT + RIGHT + "GGGGG", 25)
    assert result.score <= 2 * len(LEFT + RIGHT)
