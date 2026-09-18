"""Every coverage-qualified explanation of a molecule must reach adjudication.

`find_read_paths` enumerated from graph ROOTS only and emitted LEAVES only,
checking coverage afterwards. So an alignment that could be extended never
appeared as its own hypothesis: chr1 950M50S alone returned the chr1 path
(coverage 0.95); with a chr2 950S50M tail available, only chr1->chr2 came
back and the still-qualifying chr1-only explanation was gone. Downstream, the
molecule was then judged only against the fusion-shaped hypothesis -- the
normal-alignment competitor that `analysis` exists to weigh had been removed
before it could be weighed. Codex found it in the 2026-09-18 algorithm review;
reproduced on the live code before any change.
"""
import itertools

import pytest

from fin.fusion.paths import find_read_paths
from fin.fusion.records import segment_from_cigar


def _seg(chrom, start, cigar, source="bam_primary", mapq=60):
    return segment_from_cigar("r", chrom, start, "+", cigar, mapq, edit_distance=0, source=source)


def _window(chrom, q_start, q_end, read_length, source="fusion_genome_realign"):
    """A segment aligning query [q_start, q_end) of a read_length read."""
    parts = []
    if q_start:
        parts.append(f"{q_start}S")
    parts.append(f"{q_end - q_start}M")
    if read_length - q_end:
        parts.append(f"{read_length - q_end}S")
    return _seg(chrom, 1000, "".join(parts), source=source)


def _shapes(result):
    return sorted(tuple((s.chrom, s.query_start, s.query_end) for s in p.segments) for p in result.paths)


def test_a_normal_explanation_survives_the_arrival_of_a_placeable_tail():
    """The reproduction. Coverage 0.95 qualifies at 0.9 whether or not the
    tail can be placed somewhere else."""
    a = _seg("chr1", 100, "950M50S")
    alone = find_read_paths([a], min_mapping_quality=0, min_read_coverage=0.9)
    assert _shapes(alone) == [(("chr1", 0, 950),)]
    b = _seg("chr2", 5000, "950S50M", source="bam_supplementary")
    both = find_read_paths([a, b], min_mapping_quality=0, min_read_coverage=0.9)
    assert both.status == "ambiguous"
    assert _shapes(both) == [(("chr1", 0, 950),), (("chr1", 0, 950), ("chr2", 950, 1000))]


def test_a_qualifying_prefix_of_a_three_arm_walk_is_kept():
    a = _seg("chrA", 100, "600M400S")
    b = _seg("chrB", 100, "600S350M50S", source="bam_supplementary")
    c = _seg("chrC", 100, "950S50M", source="bam_supplementary")
    result = find_read_paths([a, b, c], min_mapping_quality=0, min_read_coverage=0.9)
    shapes = _shapes(result)
    assert (("chrA", 0, 600), ("chrB", 600, 950)) in shapes          # 0.95, the prefix
    assert (("chrA", 0, 600), ("chrB", 600, 950), ("chrC", 950, 1000)) in shapes
    assert (("chrA", 0, 600),) not in shapes                           # 0.60, fails


def test_a_qualifying_path_that_starts_at_a_non_root_node_is_kept():
    """A tiny 5' arm precedes the real alignment on the query. The walk that
    starts at the real alignment is a hypothesis of its own."""
    head = _seg("chr9", 100, "40M960S", source="bam_supplementary")
    body = _seg("chr1", 100, "40S960M")
    result = find_read_paths([head, body], min_mapping_quality=0, min_read_coverage=0.9)
    shapes = _shapes(result)
    assert (("chr1", 40, 1000),) in shapes                             # 0.96, non-root start
    assert (("chr9", 0, 40), ("chr1", 40, 1000)) in shapes


def test_the_result_does_not_depend_on_segment_input_order():
    a = _seg("chr1", 100, "950M50S")
    b = _seg("chr2", 5000, "950S50M", source="bam_supplementary")
    c = _seg("chr3", 700, "600M400S", source="fusion_genome_realign", mapq=30)
    baseline = None
    for order in itertools.permutations([a, b, c]):
        result = find_read_paths(list(order), min_mapping_quality=0, min_read_coverage=0.9)
        shapes = (result.status, _shapes(result))
        if baseline is None:
            baseline = shapes
        assert shapes == baseline


def test_no_path_is_emitted_twice():
    a = _seg("chr1", 100, "950M50S")
    b = _seg("chr2", 5000, "950S50M", source="bam_supplementary")
    result = find_read_paths([a, b], min_mapping_quality=0, min_read_coverage=0.9)
    shapes = _shapes(result)
    assert len(shapes) == len(set(shapes))


def test_the_path_limit_still_abstains_rather_than_truncating():
    # 12 mutually chainable segments: many qualifying prefixes and suffixes
    segs = [_window(f"chr{i + 1}", i * 80, (i + 1) * 80, 960) for i in range(12)]
    result = find_read_paths(segs, min_mapping_quality=0, min_read_coverage=0.5, max_paths=4)
    assert result.status == "unresolved" and result.reason == "path_limit"


def test_a_single_qualifying_path_is_still_unique():
    a = _seg("chr1", 100, "1000M")
    result = find_read_paths([a], min_mapping_quality=0, min_read_coverage=0.9)
    assert result.status == "unique" and _shapes(result) == [(("chr1", 0, 1000),)]


def test_prefixes_below_coverage_do_not_count_toward_the_limit():
    """A five-segment chain used to be ONE leaf path. It now has fifteen
    prefixes/suffixes, of which only the full walk covers 0.9. The limit must
    count what is kept, or a molecule that resolves today would become
    path_limit tomorrow."""
    segs = [_window(f"chr{i + 1}", i * 200, (i + 1) * 200, 1000) for i in range(5)]
    result = find_read_paths(segs, min_mapping_quality=0, min_read_coverage=0.9, max_paths=4)
    assert result.status == "unique"
    assert len(result.paths[0].segments) == 5


# ---------------------------------------------------------------------------
# Codex's gates: an exhaustive oracle, the exponential adversary, and the
# guarantee that nothing the old search returned is lost.
# ---------------------------------------------------------------------------


def _random_segments(generator, read_length, n):
    segs = []
    for i in range(n):
        q0 = generator.randrange(0, read_length - 60)
        q1 = generator.randrange(q0 + 40, min(read_length, q0 + 500) + 1)
        segs.append(_window(f"chr{generator.randrange(1, 6)}", q0, q1, read_length))
    return segs


def _oracle(segments, *, min_read_coverage, max_query_gap=29, max_query_overlap=10, domain=(0, 0)):
    """Every directed path under the successor rule, filtered by coverage --
    brute force, independent of the production traversal."""
    from fin.fusion.paths import ReadPath
    from fin.fusion.records import deduplicate_segments

    eligible = sorted(deduplicate_segments(segments),
                      key=lambda s: (s.query_start, s.query_end, s.geometry_key))

    def ok(prev, nxt):
        gap = nxt.query_start - prev.query_end
        return (nxt.query_start > prev.query_start and nxt.query_end > prev.query_end
                and -max_query_overlap <= gap <= max_query_gap)

    out = set()

    def walk(path):
        rp = ReadPath("r", tuple(path))
        if rp.effective_coverage(domain) >= min_read_coverage:
            out.add(tuple((s.chrom, s.query_start, s.query_end) for s in path))
        for s in eligible:
            if s.query_start > path[-1].query_start and ok(path[-1], s):
                walk(path + [s])

    for s in eligible:
        walk([s])
    return sorted(out)


def test_exhaustive_oracle_on_random_small_graphs():
    import random

    generator = random.Random(2026)
    checked = 0
    for _ in range(300):
        segs = _random_segments(generator, 1000, generator.randrange(1, 7))
        want = _oracle(segs, min_read_coverage=0.4)
        got = find_read_paths(segs, min_mapping_quality=0, min_read_coverage=0.4,
                              max_paths=10000, max_expansions=10**6)
        if got.status == "unresolved":
            assert got.reason == "insufficient_query_coverage" and want == [], (want, got.reason)
            continue
        assert _shapes(got) == want
        checked += 1
    assert checked > 50


def test_the_exponential_adversary_exhausts_the_work_budget_and_abstains():
    """20 consecutive 50-base layers, two placements each: about a million
    walks, every one of them failing 0.9 coverage, through a frontier that
    never grows large. A frontier bound would run for minutes and then say
    'insufficient coverage'; the work budget abstains at once. Codex's case."""
    segs = []
    for layer in range(20):
        for alt in range(2):
            segs.append(_window(f"chr{alt + 1}", layer * 50, (layer + 1) * 50 - 5, 2000))
    result = find_read_paths(segs, min_mapping_quality=0, min_read_coverage=0.9,
                             max_segments=128, max_paths=64, max_expansions=20000)
    assert result.status == "unresolved" and result.reason == "search_budget"
    assert result.paths == ()


def test_every_path_the_old_search_returned_is_still_returned():
    """The old traversal: roots to leaves only, coverage afterwards. Reproduced
    here as the oracle for the direction 'nothing is lost'."""
    import random

    from fin.fusion.paths import ReadPath
    from fin.fusion.records import deduplicate_segments

    def old_search(segments, min_read_coverage=0.6):
        eligible = sorted(deduplicate_segments(segments),
                          key=lambda s: (s.query_start, s.query_end, s.geometry_key))
        succ = [[] for _ in eligible]
        has_pred = set()
        for i, p in enumerate(eligible):
            for j in range(i + 1, len(eligible)):
                n = eligible[j]
                gap = n.query_start - p.query_end
                if n.query_start > p.query_start and n.query_end > p.query_end and -10 <= gap <= 29:
                    succ[i].append(j); has_pred.add(j)
        pending = [(i,) for i in reversed(range(len(eligible))) if i not in has_pred]
        done = []
        while pending:
            idx = pending.pop()
            if succ[idx[-1]]:
                pending.extend(idx + (k,) for k in reversed(succ[idx[-1]]))
            else:
                done.append(ReadPath("r", tuple(eligible[k] for k in idx)))
        return sorted(tuple((s.chrom, s.query_start, s.query_end) for s in p.segments)
                      for p in done if p.effective_coverage((0, 0)) >= min_read_coverage)

    generator = random.Random(7)
    for _ in range(300):
        segs = _random_segments(generator, 1000, generator.randrange(1, 7))
        old = old_search(segs)
        new = find_read_paths(segs, min_mapping_quality=0, min_read_coverage=0.6,
                              max_paths=10000, max_expansions=10**6)
        assert set(old) <= set(_shapes(new)), (old, _shapes(new))
