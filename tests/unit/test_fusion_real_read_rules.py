"""The rules that made the real BCR-ABL1 molecules visible.

Each was a specific failure on real direct RNA (see
docs/FUSION_REAL_DATA_2026-09-15.md): SA CIGARs are summaries that carry an
intron as a deletion; a 12/12 deletion-free anchor rejects nanopore reads
with a 2-bp deletion near a real junction; a template carrying 700 bp of
intron validated because the read simply never aligned that part; terminal
non-genomic tails failed the read-coverage gate; and a slightly-shorter
correct geometry was pruned on covered bases before any sequence judgement.
"""
import random

import pytest

from fin.fusion.annotation import AnnotatedReadPath, FusionArm
from fin.fusion.paths import coverage_domain, find_read_paths
from fin.fusion.records import AlignmentBlock, reconcile_sa_locators, segment_from_cigar
from fin.fusion.validation import validate_stitched_path


# ---------------------------------------------------------------------------
# SA locators
# ---------------------------------------------------------------------------


def _bam(cigar, chrom="chr22", start=1000, source="bam_supplementary", strand="+"):
    return segment_from_cigar("r", chrom, start, strand, cigar, 60, edit_distance=0, source=source)


def _sa(cigar, chrom="chr22", start=1000, strand="+"):
    return segment_from_cigar("r", chrom, start, strand, cigar, 60, edit_distance=0, source="sa")


def test_summarised_sa_cigar_is_superseded_by_the_spliced_record():
    # the real case: SA says 175M723D, the record says ...717N...
    record = _bam("11H22M1D49M2D28M1D2M717N26M1D9M3D22M2I15M1343H")
    locator = _sa("11S175M723D1343S")
    assert len(locator.blocks) == 1 and len(record.blocks) == 2
    kept, events = reconcile_sa_locators([record, locator])
    assert kept == (record,)
    assert len(events) == 1 and events[0]["event"] == "sa_geometry_superseded"
    assert events[0]["sa_blocks"] == 1 and events[0]["replacements"][0]["blocks"] == 2


def test_sa_without_corroboration_is_kept_as_it_is():
    locator = _sa("11S175M723D1343S")
    other = _bam("185S1262M82S", chrom="chr9", start=5000)
    kept, events = reconcile_sa_locators([other, locator])
    assert set(kept) == {other, locator} and events == ()


def test_query_overlap_alone_does_not_supersede_a_paralogous_placement():
    # same query interval, different reference position on the same chrom:
    # a paralog, not the same alignment
    record = _bam("11H175M1343H", start=1000)
    locator = _sa("11S175M1343S", start=50000)
    kept, events = reconcile_sa_locators([record, locator])
    assert set(kept) == {record, locator} and events == ()


def test_every_distinct_detailed_geometry_is_kept():
    a = _bam("11H100M75N75M1343H", source="bam_supplementary")
    b = _bam("11H100M77N75M1343H", source="fusion_genome_realign")
    locator = _sa("11S175M77D1343S")
    kept, events = reconcile_sa_locators([a, b, locator])
    assert set(kept) == {a, b}
    assert len(events[0]["replacements"]) == 2


def test_reconciliation_is_input_order_independent():
    record = _bam("11H22M1D49M2D28M1D2M717N26M1D9M3D22M2I15M1343H")
    locator = _sa("11S175M723D1343S")
    assert reconcile_sa_locators([record, locator])[0] == reconcile_sa_locators([locator, record])[0]


# ---------------------------------------------------------------------------
# coverage domain
# ---------------------------------------------------------------------------


def test_unplaced_terminal_tails_within_bounds_are_excluded():
    seg = _bam("11S1436M82S", source="bam_primary")
    assert coverage_domain([seg], 1529) == (11, 82)


def test_a_long_unplaced_tail_is_not_shaved():
    # the 400-bp tail exceeds the per-end bound and stays in the denominator
    # whole -- it is not trimmed to 100; the 11-bp head is its own eligible
    # overhang and is excluded on its own merits
    seg = _bam("11S1000M400S", source="bam_primary")
    assert coverage_domain([seg], 1411) == (11, 0)


def test_an_independently_placed_terminal_arm_is_never_excluded():
    primary = _bam("300S1000M", source="bam_primary")
    # a low-MAPQ alternative places the 300-bp head elsewhere: not a tail
    head = segment_from_cigar("r", "chr4", 100, "+", "300M1000S", 0, source="fusion_genome_realign", secondary=True)
    assert coverage_domain([primary, head], 1300) == (0, 0)


def test_sa_only_placements_do_not_define_the_domain():
    primary = _bam("100S1000M", source="bam_primary")
    locator = _sa("100M1000S", chrom="chr4", start=5)
    # the SA locator is not detailed evidence, so the 100-bp head is unplaced
    assert coverage_domain([primary, locator], 1100) == (100, 0)


def test_path_search_uses_the_domain_and_keeps_every_covered_path():
    a = _bam("100M300S", chrom="chrA", start=100, source="bam_primary")
    b = _bam("100S200M100S", chrom="chrB", start=100, source="bam_supplementary")
    b_short = _bam("100S190M110S", chrom="chrB", start=100, source="fusion_genome_realign")
    # 100 unplaced trailing bases: full-read coverage 0.75 fails at 0.9, the
    # domain (0, 100) makes it 1.0
    assert find_read_paths([a, b], min_mapping_quality=0, min_read_coverage=0.9).status == "unresolved"
    result = find_read_paths([a, b, b_short], min_mapping_quality=0, min_read_coverage=0.9, domain=(0, 100))
    assert result.status == "ambiguous"
    # the 190-bp alternative is NOT pruned for covering 10 bases fewer
    assert len(result.paths) == 2


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------


def _genome(seed=3):
    generator = random.Random(seed)
    return {c: "".join(generator.choices("ACGT", k=3000)) for c in "AB"}


def _path(arms):
    built = []
    q = 0
    for chrom, exons in arms:
        blocks = []
        for s, e in exons:
            blocks.append(AlignmentBlock(s, e, q, q + e - s, e - s))
            q += e - s
        built.append(FusionArm(chrom, "+", tuple(blocks), (chrom,)))
    return AnnotatedReadPath("r", tuple(built), "interchromosomal_candidate", "resolved")


def _query(genome, arms):
    return "".join(genome[c][s:e] for c, exons in arms for s, e in exons)


def test_error_free_read_validates_with_full_anchor_windows():
    genome = _genome()
    arms = [("A", [(100, 300), (500, 700)]), ("B", [(1000, 1300)])]
    result = validate_stitched_path(_path(arms), _query(genome, arms), genome)
    assert result.status == "validated"
    assert result.junction_anchors == ((24, 24), (24, 24))
    assert result.block_support == ((200, 200), (200, 200), (300, 300))


def test_a_two_base_deletion_beside_a_real_junction_is_tolerated():
    genome = _genome()
    arms = [("A", [(100, 300), (500, 700)]), ("B", [(1000, 1300)])]
    query = _query(genome, arms)
    # delete 2 read bases 5 bp before the A|B boundary (template position 400)
    query = query[:395] + query[397:]
    result = validate_stitched_path(_path(arms), query, genome)
    assert result.status == "validated", result
    assert min(min(pair) for pair in result.junction_anchors) >= 12


def test_low_identity_flank_fails_even_with_many_aligned_positions():
    genome = _genome()
    arms = [("A", [(100, 300), (500, 700)]), ("B", [(1000, 1300)])]
    query = list(_query(genome, arms))
    # scramble every other base in the 24 bp right of the fusion boundary
    for i in range(400, 424, 2):
        query[i] = "A" if query[i] != "A" else "C"
    result = validate_stitched_path(_path(arms), "".join(query), genome)
    assert result.status == "rejected"


def test_a_template_carrying_unobserved_intron_is_rejected_for_representation():
    genome = _genome()
    # the READ spliced out A[300:500); the TEMPLATE keeps it (one block
    # 100-700) -- the summarised-SA failure mode
    read_arms = [("A", [(100, 300), (500, 700)]), ("B", [(1000, 1300)])]
    template_arms = [("A", [(100, 700)]), ("B", [(1000, 1300)])]
    result = validate_stitched_path(_path(template_arms), _query(genome, read_arms), genome)
    assert result.status == "rejected"
    assert result.reason in ("insufficient_template_block_support", "large_indel", "unrepresented_template_end")


def test_an_unaligned_first_exon_is_rejected_as_unrepresented():
    genome = _genome()
    template_arms = [("A", [(100, 300), (500, 700)]), ("B", [(1000, 1300)])]
    read_arms = [("A", [(500, 700)]), ("B", [(1000, 1300)])]
    result = validate_stitched_path(_path(template_arms), _query(genome, read_arms), genome)
    assert result.status == "rejected"
    assert result.reason in ("insufficient_query_coverage", "unrepresented_template_end",
                             "insufficient_template_block_support")


def test_coverage_domain_excludes_a_non_genomic_tail_in_validation():
    genome = _genome()
    arms = [("A", [(100, 300), (500, 700)]), ("B", [(1000, 1300)])]
    query = _query(genome, arms) + "".join(random.Random(9).choices("ACGT", k=80))
    rejected = validate_stitched_path(_path(arms), query, genome, min_coverage=0.95)
    assert rejected.reason == "insufficient_query_coverage"
    accepted = validate_stitched_path(_path(arms), query, genome, min_coverage=0.95, domain=(0, 80))
    assert accepted.status == "validated"
    assert accepted.raw_query_coverage < accepted.query_coverage


# ---------------------------------------------------------------------------
# parallel realignment
# ---------------------------------------------------------------------------


def test_forked_realignment_returns_exactly_the_single_worker_result(tmp_path):
    """mappy releases no GIL for map(), so the only useful parallelism is
    forking after the index exists. The number of workers must not change a
    single byte of the analysis."""
    import random

    from fin.fusion.realignment import FusionGenomeAligner

    generator = random.Random(11)
    genome = {c: "".join(generator.choices("ACGT", k=4000)) for c in "AB"}
    reference = tmp_path / "genome.fa"
    reference.write_text("".join(f">{c}\n{s}\n" for c, s in genome.items()))
    molecules = []
    for index in range(6):
        start = 200 + 50 * index
        query = genome["A"][start:start + 400] + genome["B"][1000:1400]
        molecules.append((f"read{index}", query, ()))
    aligner = FusionGenomeAligner(str(reference), max_hits=16, max_windows=32, max_segments=32,
                                 min_anchor=25, min_identity=0.7, threads=1)
    with pytest.raises(ValueError, match="positive integer"):
        aligner.map_many(molecules, workers=0)
    serial = aligner.map_many(molecules, workers=1)
    forked = aligner.map_many(molecules, workers=3)
    assert set(serial) == set(forked) == {f"read{i}" for i in range(6)}
    for read_id in serial:
        assert serial[read_id] == forked[read_id]
