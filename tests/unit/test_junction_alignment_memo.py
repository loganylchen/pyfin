"""Memoising the junction DP must be an identity, not an approximation.

`align_junction_window` is the dominant cost of the fusion adjudication: 4,660
calls, 259 s of a 370 s phase, 79.5M `max()` calls inside a pure-Python
three-state affine recurrence. Codex's ranked advice was to compile the
recurrence rather than swap in another aligner -- mappy cannot even produce a
result for a window with no matching seeds, and, more sharply, two alignments
with the SAME score can carry different left/right anchor counts and therefore
different accept/reject verdicts, so score equality is not a sufficient gate.

Memoisation is the one option that needs no such gate: the core is a pure
function of the normalised window and the scoring parameters, returning a
frozen dataclass. These tests hold it to that.
"""
import random

import pytest

from fin.fusion.alignment import JunctionAlignment, _aligned_window, align_junction_window


def _fields(result):
    return {name: getattr(result, name) for name in result.__dataclass_fields__}


@pytest.fixture(autouse=True)
def clear_cache():
    _aligned_window.cache_clear()
    yield
    _aligned_window.cache_clear()


def _window(generator, query_length=40, target_length=60):
    target = "".join(generator.choices("ACGT", k=target_length))
    start = generator.randrange(0, max(1, target_length - query_length))
    query = list(target[start:start + query_length])
    for _ in range(generator.randrange(0, 5)):
        index = generator.randrange(len(query))
        query[index] = generator.choice("ACGT")
    return "".join(query), target, generator.randrange(1, target_length)


def test_a_repeated_window_returns_the_identical_object():
    generator = random.Random(7)
    query, target, junction = _window(generator)
    first = align_junction_window(query, target, junction)
    second = align_junction_window(query, target, junction)
    assert first is second
    assert _aligned_window.cache_info().hits == 1


def test_randomised_windows_agree_field_for_field_with_an_uncached_evaluation():
    generator = random.Random(11)
    for _ in range(300):
        query, target, junction = _window(generator)
        cached = align_junction_window(query, target, junction)
        _aligned_window.cache_clear()
        fresh = align_junction_window(query, target, junction)
        assert _fields(cached) == _fields(fresh)


def test_exhaustive_short_windows_agree():
    """The dense-tie regime Codex asked for: a two-letter alphabet over short
    strings, every junction offset, so repeats and equal-score alternatives are
    the common case rather than the exception."""
    alphabet = "AC"
    windows = 0
    for q in range(1, 4):
        for t in range(2, 5):
            for query in ("".join(p) for p in _products(alphabet, q)):
                for target in ("".join(p) for p in _products(alphabet, t)):
                    for junction in range(1, t):
                        cached = align_junction_window(query, target, junction, min_anchor=1)
                        _aligned_window.cache_clear()
                        fresh = align_junction_window(query, target, junction, min_anchor=1)
                        assert _fields(cached) == _fields(fresh)
                        windows += 1
    assert windows > 500


def _products(alphabet, length):
    if length == 0:
        yield ()
        return
    for head in alphabet:
        for tail in _products(alphabet, length - 1):
            yield (head,) + tail


@pytest.mark.parametrize("changed", [
    {"min_anchor": 1}, {"min_identity": 0.99}, {"max_indel": 0}, {"match_score": 3},
    {"mismatch_penalty": 1}, {"gap_open": 1}, {"gap_extend": 1},
])
def test_every_scoring_parameter_is_part_of_the_key(changed):
    """A parameter left out of the key would silently serve one caller's answer
    to another. Each of these must miss the cache, and the results must be
    computed independently rather than shared."""
    generator = random.Random(3)
    query, target, junction = _window(generator)
    baseline = align_junction_window(query, target, junction)
    before = _aligned_window.cache_info().hits
    other = align_junction_window(query, target, junction, **changed)
    assert _aligned_window.cache_info().hits == before
    assert other is not baseline


def test_max_window_is_decided_before_the_cache_is_reached():
    # max_window only gates admission, so it is not in the key -- which is only
    # safe because the rejection happens first. Prove the rejection still fires
    # for a window whose smaller sibling is already cached.
    generator = random.Random(5)
    query, target, junction = _window(generator, query_length=30, target_length=40)
    scored = align_junction_window(query, target, junction)
    assert scored.reason != "window_limit" and _aligned_window.cache_info().currsize == 1
    rejected = align_junction_window(query, target, junction, max_window=10)
    assert rejected == JunctionAlignment("unscorable", None, "window_limit")
    assert _aligned_window.cache_info().currsize == 1  # the rejection never reached it


def test_normalisation_happens_before_the_key_so_case_and_uracil_share_an_entry():
    generator = random.Random(13)
    query, target, junction = _window(generator)
    upper = align_junction_window(query, target, junction)
    lower = align_junction_window(query.lower(), target.lower(), junction)
    assert lower is upper
    uracil = align_junction_window(query.replace("T", "U"), target.replace("T", "U"), junction)
    assert uracil is upper


def test_the_unscorable_early_returns_are_not_cached_as_alignments():
    assert align_junction_window("", "ACGT", 2).reason == "invalid_window"
    assert align_junction_window("ACGT", "ACGT", 0).reason == "invalid_window"
    assert align_junction_window("ACGN", "ACGT", 2).reason == "ambiguous_sequence"
    assert _aligned_window.cache_info().currsize == 0
