"""The negative control owed by the relaxed real-data rules.

`docs/FUSION_REAL_DATA_2026-09-15.md` relaxed three things so that real direct
RNA could be seen at all: junction anchors became windowed rather than
deletion-free, the read-coverage denominator lost its unplaceable terminal
bases, and a slightly-shorter path is no longer pruned on covered bases. Each
relaxation buys tolerance, and the thing tolerance can buy by accident is an
ADAPTER CHIMERA -- two unrelated molecules ligated in one read, which presents
as two arms on different loci exactly like a fusion.

The doc asserted these rules still exclude that case and called the negative
control owed. This is it. The claim under test is specific: an adapter sits
BETWEEN the arms, and neither relaxation reaches between the arms --
`coverage_domain` only ever forgives terminal bases, and a stitched template
has nothing to align the adapter against.
"""
import random

import pytest

from fin.fusion.annotation import AnnotatedReadPath, FusionArm
from fin.fusion.paths import coverage_domain, find_read_paths
from fin.fusion.records import AlignmentBlock, segment_from_cigar
from fin.fusion.validation import validate_stitched_path


def _genome(seed=17):
    generator = random.Random(seed)
    return {c: "".join(generator.choices("ACGT", k=4000)) for c in "AB"}


def _adapter(length, seed=101):
    return "".join(random.Random(seed).choices("ACGT", k=length))


def _path(arms):
    built, q = [], 0
    for chrom, exons in arms:
        blocks = []
        for start, end in exons:
            blocks.append(AlignmentBlock(start, end, q, q + end - start, end - start))
            q += end - start
        built.append(FusionArm(chrom, "+", tuple(blocks), (chrom,)))
    return AnnotatedReadPath("r", tuple(built), "interchromosomal_candidate", "resolved")


def _stitched(genome, arms):
    return "".join(genome[c][s:e] for c, exons in arms for s, e in exons)


ARMS = [("A", [(100, 400)]), ("B", [(1000, 1300)])]


def test_the_control_itself_validates_without_an_adapter():
    # the same two arms, ligated with nothing between them, must pass -- else
    # the rejections below would prove nothing
    genome = _genome()
    result = validate_stitched_path(_path(ARMS), _stitched(genome, ARMS), genome)
    assert result.status == "validated"


@pytest.mark.parametrize("length", [12, 20, 29])
def test_an_internal_adapter_shorter_than_the_query_gap_fails_validation(length):
    """A short adapter still lets the path search join the arms (the query gap
    is within `max_query_gap`), so the rejection has to come from sequence: the
    stitched template has no template for those bases."""
    genome = _genome()
    query = genome["A"][100:400] + _adapter(length) + genome["B"][1000:1300]
    result = validate_stitched_path(_path(ARMS), query, genome)
    assert result.status == "rejected"
    assert result.reason in ("large_indel", "missing_junction_anchor",
                             "insufficient_template_block_support", "insufficient_query_coverage")


def test_a_longer_adapter_is_not_even_a_path():
    # 40 bp between the arms exceeds max_query_gap=29, so no successor edge
    # joins them and no two-arm hypothesis is ever built
    a = segment_from_cigar("r", "A", 100, "+", "300M340S", 60, edit_distance=0, source="bam_primary")
    b = segment_from_cigar("r", "B", 1000, "+", "340S300M", 60, edit_distance=0, source="bam_supplementary")
    result = find_read_paths([a, b], min_mapping_quality=0, min_read_coverage=0.5)
    assert result.status == "unresolved"
    assert all(len(path.arms) < 2 for path in result.paths)


def test_the_coverage_domain_cannot_forgive_an_internal_adapter():
    """The domain is the relaxation that could plausibly launder this, so state
    it directly: it is computed from unplaceable TERMINAL bases only, and 20
    unplaceable bases in the middle of the read buy nothing."""
    a = segment_from_cigar("r", "A", 100, "+", "300M320S", 60, edit_distance=0, source="bam_primary")
    b = segment_from_cigar("r", "B", 1000, "+", "320S300M", 60, edit_distance=0, source="bam_supplementary")
    assert coverage_domain([a, b], 620) == (0, 0)


def test_an_adapter_at_the_end_is_forgiven_but_buys_no_second_arm():
    """The contrast that makes the rule legible: the SAME adapter sequence at
    the 3' end IS excluded from the denominator -- that is the direct-RNA tail
    the domain exists for -- and it still contributes no arm, because nothing
    placed it."""
    genome = _genome()
    arms = [("A", [(100, 1400)])]
    query = _stitched(genome, arms) + _adapter(80)
    strict = validate_stitched_path(_path(arms), query, genome, min_coverage=0.95)
    assert strict.reason == "insufficient_query_coverage"
    forgiven = validate_stitched_path(_path(arms), query, genome, min_coverage=0.95, domain=(0, 80))
    assert forgiven.status == "validated"
    primary = segment_from_cigar("r", "A", 100, "+", "1300M80S", 60, edit_distance=0, source="bam_primary")
    assert coverage_domain([primary], 1380) == (0, 80)


def test_the_tail_exemption_is_capped_as_a_fraction_of_the_molecule():
    """The same 80 bp on a 380-bp molecule is 21% and is NOT forgiven -- the
    bound is on the fraction, and an overhang is taken whole or not at all, so
    it is not quietly trimmed to the 20% that would fit. A short molecule
    mostly made of adapter therefore cannot buy its way past coverage."""
    short = segment_from_cigar("r", "A", 100, "+", "300M80S", 60, edit_distance=0, source="bam_primary")
    assert coverage_domain([short], 380) == (0, 0)


def test_a_scrambled_junction_that_still_anchors_is_accepted_and_this_is_the_floor():
    """What the windowed anchor DOES admit, recorded rather than assumed.

    Scrambling all 24 bp around the join does not reject: the aligner anchors
    19 and 22 template positions across the two flanks with 17 and 16 matching.
    Because the requirement is `ceil(min_identity x aligned)`, local damage that
    costs alignment lowers the bar with it -- so the protection against a flank
    dissolving entirely is the `min_anchor` FLOOR on aligned positions, not the
    identity ratio. Both flanks clear 12 here, and 19-aligned/17-matched is a
    real anchor, so accepting it is correct.

    Stated as a test because the asymmetry is easy to misread: the doc
    justified `x aligned` over `x min_anchor` on the grounds that the latter
    would accept a 37.5%-identity flank. This is the other end of that trade.
    """
    genome = _genome()
    query = list(_stitched(genome, ARMS))
    for index in range(288, 312):
        query[index] = "A" if query[index] != "A" else "C"
    result = validate_stitched_path(_path(ARMS), "".join(query), genome)
    assert result.status == "validated"
    for aligned, matched, required, edge_gap in result.anchor_detail[0]:
        assert aligned >= 12 and matched >= required and edge_gap == 0


def test_an_anchor_below_the_floor_is_rejected_however_well_it_matches():
    """The floor itself: a flank with fewer than `min_anchor` aligned template
    positions fails even at perfect identity, which is what stops the ratio
    above from being satisfied by a handful of bases."""
    genome = _genome()
    result = validate_stitched_path(_path(ARMS), _stitched(genome, ARMS), genome, min_anchor=400)
    assert result.status == "rejected"
    assert result.reason == "missing_junction_anchor"
