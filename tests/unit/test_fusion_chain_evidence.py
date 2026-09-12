"""Gene identity, complete-chain phasing and bounded breakpoint consensus."""

from dataclasses import replace
from itertools import permutations

from fin.fusion.annotation import GeneIndex, GeneModel, annotate_read_path
from fin.fusion.chains import build_fusion_chains
from fin.fusion.paths import ReadPath
from fin.fusion.records import segment_from_cigar
from fin.io.io_gtf import GTFTranscript


def _path(read_id, chromosomes, *, starts=None):
    total = len(chromosomes) * 100
    segments = []
    starts = starts or [100] * len(chromosomes)
    for index, (chrom, start) in enumerate(zip(chromosomes, starts)):
        front = index * 100
        back = total - front - 100
        cigar = (f"{front}S" if front else "") + "100M" + (f"{back}S" if back else "")
        segments.append(segment_from_cigar(read_id, chrom, start, "+", cigar, 60))
    return annotate_read_path(ReadPath(read_id, tuple(segments)), GeneIndex())


def test_gene_overlap_is_exonic_not_whole_gene_span():
    genes = GeneIndex([GeneModel("geneA", "chr1", "+", ((100, 200), (1000, 1100)))])
    intronic = segment_from_cigar("read", "chr1", 400, "+", "100M", 60)
    assert genes.labels("chr1", "+", intronic.blocks) == ()


def test_duplicate_terminal_extents_do_not_depend_on_record_order():
    original = _path("read", ("chrA", "chrB"))
    first_arm = original.arms[0]
    shorter = replace(original, arms=(replace(
        first_arm, blocks=(replace(first_arm.blocks[0], reference_start=110,
                                  query_start=10, aligned_bases=90),)
    ), original.arms[1]))
    forward = build_fusion_chains([original, shorter], min_full_support=1)
    backward = build_fusion_chains([shorter, original], min_full_support=1)
    assert forward == backward
    assert forward[0].full_support_read_ids == ("read",)


def test_overlapping_genes_preserve_identity_ambiguity():
    genes = GeneIndex([
        GeneModel("geneA", "chr1", "+", ((100, 200),)),
        GeneModel("geneB", "chr1", "+", ((100, 200),)),
    ])
    segment = segment_from_cigar("read", "chr1", 100, "+", "100M", 60)
    annotated = annotate_read_path(ReadPath("read", (segment,)), genes)
    assert annotated.arms[0].gene_ids == ("geneA", "geneB")
    assert annotated.identity_status == "ambiguous"


def test_antisense_overlapping_gene_is_not_a_matching_identity():
    genes = GeneIndex([GeneModel("antisense", "chr1", "-", ((100, 200),))])
    segment = segment_from_cigar("read", "chr1", 100, "+", "100M", 60)
    assert genes.labels("chr1", "+", segment.blocks) == ()


def test_isoform_exon_duplicates_do_not_inflate_gene_overlap():
    first = GTFTranscript("tx1", "geneA", "chr1", "+")
    first.exons = [(100, 125)]
    second = GTFTranscript("tx2", "geneA", "chr1", "+")
    second.exons = [(100, 125)]
    genes = GeneIndex.from_transcripts([first, second])
    segment = segment_from_cigar("read", "chr1", 100, "+", "100M", 60)
    assert genes.labels("chr1", "+", segment.blocks, min_fraction=0.4) == ()


def test_splice_cigar_can_cross_two_same_chromosome_genes():
    genes = GeneIndex([
        GeneModel("geneA", "chr1", "+", ((100, 200),)),
        GeneModel("geneB", "chr1", "+", ((1000, 1100),)),
    ])
    segment = segment_from_cigar("read", "chr1", 100, "+", "100M800N100M", 60)
    annotated = annotate_read_path(ReadPath("read", (segment,)), genes)
    assert annotated.gene_path == (("geneA",), ("geneB",))
    assert len(annotated.junctions) == 1
    assert annotated.event_class == "readthrough_or_splicing_candidate"


def test_same_gene_splicing_is_not_called_a_fusion():
    genes = GeneIndex([GeneModel("geneA", "chr1", "+", ((100, 200), (1000, 1100)))])
    segment = segment_from_cigar("read", "chr1", 100, "+", "100M800N100M", 60)
    annotated = annotate_read_path(ReadPath("read", (segment,)), genes)
    assert len(annotated.arms) == 1
    assert annotated.event_class == "ordinary_alignment"
    assert build_fusion_chains([annotated]) == ()


def test_reverse_spliced_alignment_reports_transcript_gene_order():
    genes = GeneIndex([
        GeneModel("geneA", "chr1", "-", ((100, 200),)),
        GeneModel("geneB", "chr1", "-", ((1000, 1100),)),
    ])
    segment = segment_from_cigar("read", "chr1", 100, "-", "100M800N100M", 60)
    annotated = annotate_read_path(ReadPath("read", (segment,)), genes)
    assert annotated.gene_path == (("geneB",), ("geneA",))


def test_three_genes_need_whole_molecule_support_not_pairwise_graph_connections():
    paths = [_path("ab1", ("chr1", "chr2")), _path("ab2", ("chr1", "chr2")),
             _path("bc1", ("chr2", "chr3")), _path("bc2", ("chr2", "chr3"))]
    chains = build_fusion_chains(paths)
    assert len(chains) == 2
    assert all(len(chain.representative.arms) == 2 for chain in chains)
    assert all(chain.status == "supported_candidate" for chain in chains)


def test_partial_reads_cannot_promote_one_complete_three_gene_read():
    paths = [_path("abc", ("chr1", "chr2", "chr3")),
             _path("ab1", ("chr1", "chr2")), _path("ab2", ("chr1", "chr2")),
             _path("bc1", ("chr2", "chr3")), _path("bc2", ("chr2", "chr3"))]
    chain = next(chain for chain in build_fusion_chains(paths) if len(chain.representative.arms) == 3)
    assert chain.status == "insufficient_full_chain_support"
    assert chain.full_support_read_ids == ("abc",)
    assert chain.partial_support_read_ids == ("ab1", "ab2", "bc1", "bc2")
    assert chain.edge_support_read_ids == (("ab1", "ab2", "abc"), ("abc", "bc1", "bc2"))


def test_complete_three_gene_molecules_support_all_edges():
    paths = [_path("abc1", ("chr1", "chr2", "chr3")),
             _path("abc2", ("chr1", "chr2", "chr3"))]
    chain, = build_fusion_chains(paths)
    assert chain.status == "supported_candidate"
    assert chain.edge_support_read_ids == (("abc1", "abc2"), ("abc1", "abc2"))


def test_multiple_fusion_chains_coexist_without_cartesian_products():
    paths = [_path("abc1", ("chr1", "chr2", "chr3")),
             _path("abc2", ("chr1", "chr2", "chr3")),
             _path("abd1", ("chr1", "chr2", "chr4")),
             _path("abd2", ("chr1", "chr2", "chr4"))]
    chains = build_fusion_chains(paths)
    assert len(chains) == 2
    assert all(chain.status == "supported_candidate" for chain in chains)
    assert {tuple(arm.chrom for arm in chain.representative.arms) for chain in chains} == {
        ("chr1", "chr2", "chr3"), ("chr1", "chr2", "chr4"),
    }


def test_ambiguous_molecule_does_not_support_both_paths():
    paths = [_path("ambiguous", ("chr1", "chr2")),
             _path("ambiguous", ("chr1", "chr3")),
             _path("one_good", ("chr1", "chr2"))]
    chain, = build_fusion_chains(paths)
    assert chain.full_support_read_ids == ("one_good",)
    assert chain.status == "insufficient_full_chain_support"


def test_duplicate_input_does_not_inflate_full_support():
    path = _path("one_molecule", ("chr1", "chr2"))
    chain, = build_fusion_chains([path] * 10)
    assert chain.full_support_read_ids == ("one_molecule",)
    assert chain.status == "insufficient_full_chain_support"


def test_single_linkage_bridge_cannot_merge_distant_breakpoints():
    paths = [_path("left", ("chr1", "chr2"), starts=[100, 100]),
             _path("bridge", ("chr1", "chr2"), starts=[105, 105]),
             _path("right", ("chr1", "chr2"), starts=[110, 110])]
    signatures = set()
    for order in permutations(paths):
        chains = build_fusion_chains(order, breakpoint_tolerance=6)
        assert len(chains) == 2
        assert not any("left" in chain.full_support_read_ids and "right" in chain.full_support_read_ids
                       for chain in chains)
        signatures.add(tuple((chain.chain_id, chain.full_support_read_ids) for chain in chains))
    assert len(signatures) == 1


def test_nearby_coordinate_modes_are_not_reported_as_exact_correction():
    paths = [_path("first", ("chr1", "chr2"), starts=[100, 100]),
             _path("second", ("chr1", "chr2"), starts=[102, 102])]
    chain, = build_fusion_chains(paths)
    assert chain.coordinate_status == "unresolved_nearby_modes"
    assert len(chain.raw_coordinate_modes) == 2


def test_changing_read_names_does_not_change_event_identity():
    first, = build_fusion_chains([_path("one", ("chr1", "chr2")), _path("two", ("chr1", "chr2"))])
    second, = build_fusion_chains([_path("other1", ("chr1", "chr2")), _path("other2", ("chr1", "chr2"))])
    assert first.chain_id == second.chain_id


def test_ambiguous_gene_labels_do_not_seed_a_confident_gene_chain():
    path = _path("read", ("chr1", "chr2"))
    ambiguous = replace(path, identity_status="ambiguous")
    assert build_fusion_chains([ambiguous]) == ()
