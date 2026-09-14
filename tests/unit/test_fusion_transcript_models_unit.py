"""Model-layer invariants on hand-built analyses (no aligner).

The integration fixture covers what real mappy paths produce. These pin the
cases Codex reproduced in memory that an aligner will not conveniently emit:
one molecule holding two selected paths under one splice key, same-gene
structures with different internal exons that spell the same sequence,
collapse followed by a microhomology merge, and same-sequence distinct loci.
"""
import json

import pytest

from fin.fusion.analysis import FusionAnalysis, FusionReadAudit, FusionSettings
from fin.fusion.annotation import AnnotatedReadPath, FusionArm
from fin.fusion.chains import FusionChain
from fin.fusion.records import AlignmentBlock
from fin.fusion.transcript_models import (
    TranscriptIndex, TranscriptModel, build_fused_transcript_models, model_sequence,
)


def _arm(chrom, strand, exons, gene, query_start=0):
    """FusionArm with blocks in TRANSCRIPT order (descending on minus)."""
    ordered = list(exons) if strand == "+" else list(reversed(exons))
    blocks = []
    q = query_start
    for s, e in ordered:
        blocks.append(AlignmentBlock(s, e, q, q + (e - s), e - s))
        q += e - s
    return FusionArm(chrom, strand, tuple(blocks), (gene,) if gene else ()), q


def _path(read_id, arms):
    built = []
    q = 0
    for chrom, strand, exons, gene in arms:
        arm, q = _arm(chrom, strand, exons, gene, q)
        built.append(arm)
    return AnnotatedReadPath(read_id, tuple(built), "intergenic_fusion_candidate", "resolved")


def _analysis(molecules, supported_ids, *, min_support=2):
    """molecules: {read_id: [path, ...]}; chains built one per splice topology
    by hand: every read id listed in `supported_ids` goes into a supported
    chain, the rest into insufficient-support chains, one chain per read."""
    reads = tuple(FusionReadAudit(rid, "supported_path", "", path_count=len(paths), selected=tuple(paths))
                  for rid, paths in molecules.items())
    chains = []
    if supported_ids:
        chains.append(FusionChain("fusion_chain_supported", molecules[supported_ids[0]][0], tuple(supported_ids), (),
                                  tuple(tuple(supported_ids) for _ in range(len(molecules[supported_ids[0]][0].arms) - 1)),
                                  (), "supported_rna_fusion", "sequence_resolved"))
    for rid in molecules:
        if rid not in supported_ids:
            chains.append(FusionChain(f"fusion_chain_{rid}", molecules[rid][0], (rid,), (),
                                      tuple((rid,) for _ in range(len(molecules[rid][0].arms) - 1)),
                                      (), "insufficient_full_chain_support", "sequence_resolved"))
    return FusionAnalysis(tuple(chains), reads, {}, FusionSettings(min_support=min_support), "provided")


GENOME = {"A": "ACGT" * 250, "B": "TTGCA" * 200, "D": "ACGT" * 250}
INDEX = TranscriptIndex([
    TranscriptModel("A1", "A", "gA", "A", "+", ((100, 200), (300, 400))),
    TranscriptModel("B1", "B", "gB", "B", "+", ((500, 600),)),
])


def test_one_molecule_with_two_paths_under_one_splice_key_is_judged_per_path():
    # long structure supported by two molecules; S1 holds TWO selected paths
    # with the same splice key: [(250,400)] (aligned across A1's intron, so it
    # cannot fit the long model) and [(350,400)] (fits). S2 = [(350,400)].
    long_a = _path("L1", [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    long_b = _path("L2", [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    s1_misfit = _path("S1", [("A", "+", [(250, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    s1_fit = _path("S1", [("A", "+", [(350, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    s2 = _path("S2", [("A", "+", [(350, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    analysis = _analysis({"L1": [long_a], "L2": [long_b], "S1": [s1_misfit, s1_fit], "S2": [s2]},
                         ["L1", "L2", "S1", "S2"])
    models, collapsed = build_fused_transcript_models(analysis, INDEX, GENOME)
    # the misfit observation keeps the short structure alive as a model
    assert len(models) == 2
    by_len = {len(model_sequence(m.arms, GENOME)[0]): m for m in models}
    assert set(by_len) == {300, 250}
    assert by_len[300].unique_read_ids == ("L1", "L2")
    # S1 fits the long model through one path and the short through the other:
    # compatible with both, unique to neither. S2 likewise.
    assert "S1" not in by_len[300].unique_read_ids + by_len[250].unique_read_ids
    assert "S1" in by_len[300].compatible_read_ids and "S1" in by_len[250].compatible_read_ids
    assert collapsed == ()


def test_same_gene_structures_with_different_internal_exons_are_not_merged_on_sequence():
    # chrA is a 4-mer repeat, so exons [(20,30)] and [(60,70)] spell the same
    # sequence; X and Y differ only there. Same breakpoint, same genes.
    x1 = _path("X1", [("A", "+", [(0, 10), (20, 30), (90, 100)], "A"), ("B", "+", [(500, 600)], "B")])
    x2 = _path("X2", [("A", "+", [(0, 10), (20, 30), (90, 100)], "A"), ("B", "+", [(500, 600)], "B")])
    y1 = _path("Y1", [("A", "+", [(0, 10), (60, 70), (90, 100)], "A"), ("B", "+", [(500, 600)], "B")])
    y2 = _path("Y2", [("A", "+", [(0, 10), (60, 70), (90, 100)], "A"), ("B", "+", [(500, 600)], "B")])
    analysis = _analysis({"X1": [x1], "X2": [x2], "Y1": [y1], "Y2": [y2]}, ["X1", "X2", "Y1", "Y2"])
    models, _ = build_fused_transcript_models(analysis, INDEX, GENOME)
    assert len(models) == 2
    sequences = {model_sequence(m.arms, GENOME)[0] for m in models}
    assert len(sequences) == 1               # same sequence...
    assert {m.arms[0].exons for m in models} == {((0, 10), (20, 30), (90, 100)), ((0, 10), (60, 70), (90, 100))}
    assert all(m.coordinate_status != "unresolved_microhomology_same_sequence" for m in models)
    assert sorted(m.unique_read_ids for m in models) == [("X1", "X2"), ("Y1", "Y2")]


def test_distinct_loci_with_identical_sequence_stay_distinct_models():
    a = [_path(f"A{i}", [("A", "+", [(100, 200)], "A"), ("B", "+", [(500, 600)], "B")]) for i in range(2)]
    d = [_path(f"D{i}", [("D", "+", [(100, 200)], "D"), ("B", "+", [(500, 600)], "B")]) for i in range(2)]
    analysis = _analysis({p.read_id: [p] for p in a + d}, [p.read_id for p in a + d])
    models, _ = build_fused_transcript_models(analysis, INDEX, GENOME)
    assert len(models) == 2
    assert {m.arms[0].chrom for m in models} == {"A", "D"}
    assert all(not m.coordinate_alternatives for m in models)


def test_collapsed_audit_targets_survive_a_microhomology_merge():
    # two same-sequence placements of one junction (breakpoint 400 vs 404 on a
    # 4-mer repeat genome) are merged; a degraded molecule compatible with the
    # placement that is NOT kept as the record must point at the SURVIVING id
    p1 = [_path(f"P{i}", [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")]) for i in range(2)]
    p2 = [_path(f"Q{i}", [("A", "+", [(100, 200), (300, 404)], "A"), ("B", "+", [(504, 600)], "B")]) for i in range(2)]
    deg = _path("DEG", [("A", "+", [(350, 404)], "A"), ("B", "+", [(504, 600)], "B")])
    analysis = _analysis({p.read_id: [p] for p in p1 + p2 + [deg]}, [p.read_id for p in p1 + p2])
    # a genuine microhomology: the 4 bases either side of the junction are the
    # same on both chromosomes, so both placements spell one sequence
    genome = {"A": "ACGT" * 250, "B": "ACGT" * 250}
    models, collapsed = build_fused_transcript_models(analysis, INDEX, genome)
    assert len(models) == 1
    model = models[0]
    assert model.coordinate_status == "unresolved_microhomology_same_sequence"
    assert len(model.coordinate_alternatives) == 1
    assert model.breakpoints == (("A", 400, "+", "B", 500, "+"),)     # annotation-agreeing placement is the record
    assert len(collapsed) == 1
    assert collapsed[0].read_ids == ("DEG",)
    assert collapsed[0].compatible_model_ids == (model.model_id,)
    assert "DEG" in model.unique_read_ids                               # fits the group and nothing else


def test_crossing_tss_envelope_is_resolved_after_widening():
    # observations [130,200) and [90,200): the envelope starts at 90, BEFORE
    # A1's annotated start (100); A1 must then be incompatible (the observed
    # exon is not inside the retained slice) -> no completion, never a
    # negative extension; the model still reports.
    a = _path("R1", [("A", "+", [(130, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    b = _path("R2", [("A", "+", [(90, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    for order in (["R1", "R2"], ["R2", "R1"]):
        analysis = _analysis({rid: [a if rid == "R1" else b] for rid in order}, order)
        models, _ = build_fused_transcript_models(analysis, INDEX, GENOME)
        assert len(models) == 1
        assert models[0].arms[0].exons == ((90, 200), (300, 400))
        assert models[0].arms[0].resolution.completion_exons == ()
        assert models[0].five_prime_status == "observed_boundary_5p_completeness_unknown"
        assert models[0].unique_read_ids == ("R1", "R2")


def test_minus_strand_arm_conversion_and_sequence():
    # a minus arm given in transcript order (descending) must come out with
    # ascending exons and a sequence = revcomp(exon 2) + revcomp(exon 1)
    genome = {"M": "A" * 100 + "C" * 10 + "A" * 90 + "G" * 10 + "A" * 790, "B": "T" * 1000}
    index = TranscriptIndex([TranscriptModel("M1", "M", "gM", "M", "-", ((100, 110), (200, 210))),
                             TranscriptModel("B1", "B", "gB", "B", "+", ((500, 600),))])
    paths = [_path(f"m{i}", [("M", "-", [(100, 110), (200, 210)], "M"), ("B", "+", [(500, 600)], "B")]) for i in range(2)]
    analysis = _analysis({p.read_id: [p] for p in paths}, [p.read_id for p in paths])
    models, _ = build_fused_transcript_models(analysis, index, genome)
    assert models[0].arms[0].exons == ((100, 110), (200, 210))
    seq, rows = model_sequence(models[0].arms, genome)
    assert seq[:20] == "C" * 10 + "G" * 10          # revcomp(G*10) then revcomp(C*10)
    assert [(r[4], r[5]) for r in rows[:2]] == [(200, 210), (100, 110)]
    assert models[0].arms[0].resolution.match_class == "arm_fsm"
    assert models[0].arms[0].resolution.breakpoint_exon_ranks == (("M1", 2),)
    assert models[0].breakpoints == (("M", 100, "-", "B", 500, "+"),)


def test_unavailable_annotated_completion_sequence_fails_the_writer(tmp_path):
    # the completion extends the first arm to A1's start [100,200); a genome
    # that stops short of it makes the OBSERVED sequence available but the
    # annotated completion unavailable -- that must raise, not vanish
    from fin.fusion.transcript_output import write_transcript_outputs

    paths = [_path(f"r{i}", [("A", "+", [(150, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")]) for i in range(2)]
    analysis = _analysis({p.read_id: [p] for p in paths}, [p.read_id for p in paths])
    models, collapsed = build_fused_transcript_models(analysis, INDEX, GENOME)
    assert models[0].five_prime_status == "annotated_extension_available"

    class Genome(dict):
        pass

    # observed slices resolve, the completion's [100,150) does not
    genome = Genome({"A": GENOME["A"], "B": GENOME["B"]})
    original_slice = __import__("fin.fusion.transcript_models", fromlist=["_slice"])._slice

    def failing(genome_, chrom, start, end):
        if chrom == "A" and start < 150:
            return ""
        return original_slice(genome_, chrom, start, end)

    import fin.fusion.transcript_models as module
    module_slice = module._slice
    module._slice = failing
    try:
        with pytest.raises(ValueError, match="annotated completion sequence is unavailable"):
            write_transcript_outputs(models, collapsed, genome, tmp_path)
    finally:
        module._slice = module_slice


# ---------------------------------------------------------------------------
# the second review's five reproduced defects
# ---------------------------------------------------------------------------


def test_two_decisive_placements_on_different_molecules_are_not_pooled():
    # each molecule settled on ONE placement, 3 bp apart; proximity is not
    # evidence, so neither becomes support for the other and nothing reports
    a = _path("A", [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    b = _path("B", [("A", "+", [(100, 200), (303, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    analysis = _analysis({"A": [a], "B": [b]}, ["A", "B"])
    models, _ = build_fused_transcript_models(analysis, INDEX, GENOME)
    assert models == ()


def test_placements_held_by_one_molecule_as_alternatives_do_merge_within_diameter():
    # both molecules carry BOTH placements (the corrector could not settle):
    # one model, both unique, the other placement listed
    a1 = _path("A", [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    a2 = _path("A", [("A", "+", [(100, 200), (303, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    b1 = _path("B", [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    b2 = _path("B", [("A", "+", [(100, 200), (303, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    analysis = _analysis({"A": [a1, a2], "B": [b1, b2]}, ["A", "B"])
    models, _ = build_fused_transcript_models(analysis, INDEX, GENOME)
    assert len(models) == 1
    assert models[0].unique_read_ids == ("A", "B")
    assert models[0].coordinate_status == "unresolved_nearby_placements"
    assert models[0].arms[0].exons == ((100, 200), (300, 400))       # annotation-agreeing record


def test_no_single_linkage_chaining_across_the_tolerance():
    # offsets 0, 6 and 12 with tolerance 6: a bridge molecule holding 0 and
    # 6, another holding 6 and 12. 0 and 12 are 12 apart: never one group.
    p0 = lambda r: _path(r, [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    p6 = lambda r: _path(r, [("A", "+", [(100, 200), (306, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    p12 = lambda r: _path(r, [("A", "+", [(100, 200), (312, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    analysis = _analysis({"X": [p0("X"), p6("X")], "Y": [p6("Y"), p12("Y")], "Z": [p0("Z")], "W": [p12("W")]},
                         ["X", "Y", "Z", "W"])
    models, _ = build_fused_transcript_models(analysis, INDEX, GENOME)
    assert len(models) >= 2                          # 0/6 and 6/12 are distinct populations
    for m in models:
        record = m.arms[0].exons[1][0]
        assert record in (300, 306, 312)
        # every listed alternative is within the tolerance of the record at
        # the INTERNAL acceptor -- the fusion breakpoint is identical across
        # all three, so checking it would be vacuous
        for alt in m.coordinate_alternatives:
            placement = json.loads(alt)
            assert abs(placement[0][2][1][0] - record) <= 6
    assert not any(len(m.unique_read_ids) == 4 for m in models)


def test_alias_chains_resolve_across_both_merge_passes():
    """A genuine two-stage retirement: Q -> P (same sequence) -> new id.

    P and Q are identical-sequence placements of one junction (400/500 vs
    404/504 on a 4-mer repeat genome), so the first pass retires Q behind P.
    R shares P's junction but shifts the internal acceptor 3 bp and reaches
    50 bp further at the 3' end; M holds both P and R, so the second pass
    merges them -- and because the aggregate envelope takes R's longer 3'
    extent, the surviving model gets a NEW id and P itself is retired. Q's
    alias must follow through both hops.
    """
    genome = {"A": "ACGT" * 250, "B": "ACGT" * 250}
    P = [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")]
    Q = [("A", "+", [(100, 200), (300, 404)], "A"), ("B", "+", [(504, 600)], "B")]
    R = [("A", "+", [(100, 200), (297, 400)], "A"), ("B", "+", [(500, 650)], "B")]
    molecules = {"P0": [_path("P0", P)], "P1": [_path("P1", P)],
                 "Q0": [_path("Q0", Q)], "Q1": [_path("Q1", Q)],
                 "M": [_path("M", P), _path("M", R)], "N": [_path("N", R)]}
    analysis = _analysis(molecules, list(molecules))
    models, collapsed = build_fused_transcript_models(analysis, INDEX, genome)   # must not raise
    assert len(models) == 1
    model = models[0]
    emitted = {m.model_id for m in models}
    for item in collapsed:
        assert set(item.compatible_model_ids) <= emitted
    # both retirements happened: the id is neither P's nor Q's own structure
    assert model.model_id not in {_model_id_of(P), _model_id_of(Q)}
    # the aggregate envelope kept R's 3' extent and P's internal placement
    assert model.arms[0].exons == ((100, 200), (300, 400))
    assert model.arms[-1].exons == ((500, 650),)
    # every retired placement is still listed, including the one inherited
    # from the first pass
    placements = {tuple(tuple(a[2][i]) for i in range(len(a[2]))) for alt in model.coordinate_alternatives
                  for a in json.loads(alt)}
    assert ((100, 200), (300, 404)) in placements          # Q, via the alias chain
    assert ((100, 200), (297, 400)) in placements          # R
    assert set(model.unique_read_ids) == set(molecules)
    # the extended 3' end is at 650; the annotated B exon ends at 600, so the
    # regenerated envelope must NOT inherit the members' "annotated boundary
    # observed" -- its own resolution does not support that claim
    assert model.arms[-1].exons == ((500, 650),)
    assert model.arms[-1].resolution.completion_exons == ()
    assert model.three_prime_status == "observed_end_completeness_unknown"


def _model_id_of(arms_spec):
    from fin.fusion.transcript_models import ModelArm, _model_id
    return _model_id(tuple(ModelArm(c, s, tuple(e), None) for c, s, e, _ in arms_spec))


def test_endpoint_reach_is_recomputed_against_the_final_envelope():
    # observation R1 starts at 100 (exact fit), R2 starts at 103 and holds
    # both placements of the acceptor; after the merge only R1 reaches the
    # surviving envelope's 5' extent
    r1 = _path("R1", [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    r2a = _path("R2", [("A", "+", [(103, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    r2b = _path("R2", [("A", "+", [(103, 200), (303, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    r3 = _path("R3", [("A", "+", [(103, 200), (303, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    r3b = _path("R3", [("A", "+", [(103, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    analysis = _analysis({"R1": [r1], "R2": [r2a, r2b], "R3": [r3, r3b]}, ["R1", "R2", "R3"])
    models, _ = build_fused_transcript_models(analysis, INDEX, GENOME)
    assert len(models) == 1
    assert models[0].arms[0].exons[0] == (100, 200)
    assert models[0].five_prime_observed_read_ids == ("R1",)


def test_mixed_sequence_group_does_not_keep_the_same_sequence_claim():
    genome = {"A": "ACGT" * 250, "B": "ACGT" * 250}
    p = [_path(f"P{i}", [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")]) for i in range(2)]
    q = [_path(f"Q{i}", [("A", "+", [(100, 200), (300, 404)], "A"), ("B", "+", [(504, 600)], "B")]) for i in range(2)]
    # a molecule holding the same-sequence survivor AND a 3-bp-shorter placement
    s1 = _path("S", [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    s2 = _path("S", [("A", "+", [(100, 200), (300, 397)], "A"), ("B", "+", [(500, 600)], "B")])
    t = _path("T", [("A", "+", [(100, 200), (300, 397)], "A"), ("B", "+", [(500, 600)], "B")])
    molecules = {x.read_id: [x] for x in p + q}
    molecules["S"] = [s1, s2]
    molecules["T"] = [t]
    analysis = _analysis(molecules, list(molecules))
    models, _ = build_fused_transcript_models(analysis, INDEX, genome)
    assert len(models) == 1
    assert models[0].coordinate_status == "unresolved_nearby_placements"


def test_nearby_decisive_molecule_is_reported_not_counted():
    a = [_path(f"A{i}", [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 600)], "B")]) for i in range(2)]
    n = _path("N", [("A", "+", [(100, 200), (302, 400)], "A"), ("B", "+", [(500, 600)], "B")])
    analysis = _analysis({**{x.read_id: [x] for x in a}, "N": [n]}, ["A0", "A1", "N"])
    models, _ = build_fused_transcript_models(analysis, INDEX, GENOME)
    assert len(models) == 1
    assert models[0].unique_read_ids == ("A0", "A1")
    assert models[0].nearby_read_ids == ("N",)
    assert models[0].coordinate_status != "unresolved_nearby_placements"


def test_a_regenerated_envelope_is_not_called_same_sequence_on_its_members():
    # an all-A reference: P and Q spell the same 250 bases, but the aggregate
    # record combines P's 5' extent with Q's 3' extent and is 297 -- a
    # sequence neither member has, so the stronger claim must not survive
    genome = {"A": "A" * 1000, "B": "A" * 1000}
    index = TranscriptIndex([TranscriptModel("A1", "A", "gA", "A", "+", ((100, 200), (300, 400))),
                             TranscriptModel("B1", "B", "gB", "B", "+", ((500, 600),))])
    P = [("A", "+", [(100, 200), (300, 400)], "A"), ("B", "+", [(500, 550)], "B")]
    Q = [("A", "+", [(150, 200), (297, 400)], "A"), ("B", "+", [(500, 597)], "B")]
    molecules = {"M": [_path("M", P), _path("M", Q)], "N": [_path("N", P), _path("N", Q)]}
    analysis = _analysis(molecules, list(molecules))
    models, _ = build_fused_transcript_models(analysis, index, genome)
    assert len(models) == 1
    model = models[0]
    assert len(model_sequence(model.arms, genome)[0]) == 297
    assert model.coordinate_status == "unresolved_nearby_placements"
