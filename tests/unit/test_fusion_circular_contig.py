"""A read across a circular reference's origin is not a rearrangement.

On the whole K562 library, 12 of the 69 fused-transcript models were on chrM,
and they held the two HIGHEST-support models of the entire run (390 and 91
molecules). Six of the twelve are a read continuing from the end of the
linearised mitochondrial reference to its start: chrM:16569:+ -> chrM:0:+, or
on the minus strand chrM:0:- -> chrM:16564:-. The classifier compares
coordinates, so exit > entry is 'not collinear' and the path lands as an
intrachromosomal rearrangement candidate. The reference is a circle; the
molecule is continuous.

The other six chrM models are INTERNAL jumps (chrM:12134 -> 13867, thousands
of bases from either end). Circularity explains none of them, so this must
not become a blanket chrM exclusion -- and the label says what it means:
the reference's topology explains the discontinuity, not that the RNA is
circular or artefactual. On the frozen run the origin-spanning junctions sit
within 7 bp of the ends and the internal ones no closer than 4,400 bp, so a
50 bp tolerance separates them with a wide margin.
"""
import pytest

from fin.fusion.annotation import AnnotatedReadPath, GeneIndex, GeneModel, annotate_read_path
from fin.fusion.chains import build_fusion_chains
from fin.fusion.paths import ReadPath
from fin.fusion.records import segment_from_cigar

CHRM = 16569
CIRCULAR = {"chrM": CHRM}


def _path(*segs):
    return ReadPath("r", tuple(segs))


def _seg(chrom, start, strand, cigar, source="bam_primary"):
    return segment_from_cigar("r", chrom, start, strand, cigar, 60, edit_distance=0, source=source)


def test_plus_strand_end_to_start_is_origin_spanning():
    # 300 bases ending exactly at the reference end, then 200 from position 0
    a = _seg("chrM", CHRM - 300, "+", "300M200S")
    b = _seg("chrM", 0, "+", "300S200M", source="bam_supplementary")
    path = annotate_read_path(_path(a, b), GeneIndex(), circular_contigs=CIRCULAR)
    assert len(path.arms) == 2
    assert path.event_class == "circular_reference_origin_spanning"


def test_minus_strand_start_to_end_is_origin_spanning():
    """The frozen run's 390-molecule model: chrM:0:- -> chrM:16564:-. On the
    minus strand the transcript runs high->low, so the arm that EXITS at the
    reference start continues from the reference end."""
    # CIGARs are reference-oriented: the read's first 300 bases map to
    # ref [0,300) on the minus strand, so in reference orientation they are
    # the reverse-complemented read's LAST 300 bases -> "200S300M"; the next
    # 200 read bases map to [CHRM-200, CHRM) -> "200M300S".
    a = _seg("chrM", 0, "-", "200S300M")
    b = _seg("chrM", CHRM - 200, "-", "200M300S", source="bam_supplementary")
    path = annotate_read_path(_path(a, b), GeneIndex(), circular_contigs=CIRCULAR)
    assert path.event_class == "circular_reference_origin_spanning"


def test_a_few_bases_short_of_the_ends_still_counts():
    # the frozen run: exit 7 bp before the end, entry at 0
    a = _seg("chrM", CHRM - 307, "+", "300M200S")
    b = _seg("chrM", 0, "+", "300S200M", source="bam_supplementary")
    assert annotate_read_path(_path(a, b), GeneIndex(), circular_contigs=CIRCULAR).event_class \
        == "circular_reference_origin_spanning"


def test_an_internal_mitochondrial_jump_is_not_origin_spanning():
    """chrM:12134 -> chrM:13867 is a deletion-shaped discontinuity thousands of
    bases from either end. Circularity has nothing to say about it, and it
    keeps whatever class the geometry gives it."""
    a = _seg("chrM", 11834, "+", "300M200S")
    b = _seg("chrM", 13867, "+", "300S200M", source="bam_supplementary")
    path = annotate_read_path(_path(a, b), GeneIndex(), circular_contigs=CIRCULAR)
    assert path.event_class != "circular_reference_origin_spanning"
    assert path.event_class == "readthrough_or_splicing_candidate"  # collinear, within distance


def test_a_contig_not_declared_circular_is_classified_as_before():
    a = _seg("chr1", 248956422 - 300, "+", "300M200S")
    b = _seg("chr1", 0, "+", "300S200M", source="bam_supplementary")
    with_decl = annotate_read_path(_path(a, b), GeneIndex(), circular_contigs=CIRCULAR)
    without = annotate_read_path(_path(a, b), GeneIndex())
    assert with_decl.event_class == without.event_class == "intrachromosomal_candidate"


def test_chrm_to_nuclear_is_not_exempted_by_containing_chrm():
    a = _seg("chrM", CHRM - 300, "+", "300M200S")
    b = _seg("chr1", 1000, "+", "300S200M", source="bam_supplementary")
    path = annotate_read_path(_path(a, b), GeneIndex(), circular_contigs=CIRCULAR)
    assert path.event_class == "interchromosomal_candidate"


def test_start_to_end_on_the_plus_strand_is_not_the_seam():
    """Plus strand, read leaves near the START and re-enters near the END.
    That is a forward jump of 16 kb -- collinear, an ordinary (if huge)
    intron-shaped junction -- not the seam, which on + runs end->start."""
    a = _seg("chrM", 0, "+", "300M200S")
    b = _seg("chrM", CHRM - 200, "+", "300S200M", source="bam_supplementary")
    path = annotate_read_path(_path(a, b), GeneIndex(), circular_contigs=CIRCULAR)
    assert path.event_class != "circular_reference_origin_spanning"
    assert path.event_class == "readthrough_or_splicing_candidate"


def test_an_origin_junction_plus_an_ordinary_splice_is_still_origin_spanning():
    """Three arms: an ordinary collinear splice inside chrM, then the origin.
    The topology explains the discontinuity; the splice is ordinary."""
    genes = GeneIndex([GeneModel("MT-A", "chrM", "+", ((15000, 15400),)),
                       GeneModel("MT-B", "chrM", "+", ((16000, CHRM),)),
                       GeneModel("MT-C", "chrM", "+", ((0, 600),))])
    a = _seg("chrM", 15000, "+", "400M600N569M200S")
    b = _seg("chrM", 0, "+", "969S200M", source="bam_supplementary")
    path = annotate_read_path(_path(a, b), genes, circular_contigs=CIRCULAR)
    assert path.event_class == "circular_reference_origin_spanning"


def test_an_origin_junction_plus_a_real_discontinuity_is_not_excused():
    """Origin crossing AND a reverse-order jump elsewhere: circularity explains
    one junction, not the other. The path keeps its rearrangement class."""
    a = _seg("chrM", CHRM - 300, "+", "300M400S")
    b = _seg("chrM", 0, "+", "300S200M200S", source="bam_supplementary")
    # after the seam the read jumps BACK to 100 while b exited at 200:
    # reverse order, 100 bp from the start is outside the 50 bp seam tolerance
    c = _seg("chrM", 100, "+", "500S200M", source="bam_supplementary")
    path = annotate_read_path(_path(a, b, c), GeneIndex(), circular_contigs=CIRCULAR)
    assert path.event_class == "intrachromosomal_candidate"


def test_chains_never_promote_an_origin_spanning_path_to_supported_fusion():
    paths = []
    for rid in ("m1", "m2", "m3"):
        a1 = segment_from_cigar(rid, "chrM", CHRM - 300, "+", "300M200S", 60, edit_distance=0, source="bam_primary")
        b1 = segment_from_cigar(rid, "chrM", 0, "+", "300S200M", 60, edit_distance=0, source="bam_supplementary")
        paths.append(annotate_read_path(ReadPath(rid, (a1, b1)), GeneIndex(), circular_contigs=CIRCULAR))
    chains = build_fusion_chains(paths, min_full_support=2)
    assert len(chains) == 1
    assert chains[0].status == "circular_reference_origin"
    assert chains[0].status != "supported_rna_fusion"
    assert chains[0].full_support_read_ids == ("m1", "m2", "m3")   # support is kept, not deleted


def test_the_default_declaration_covers_the_common_mitochondrial_names():
    from fin.fusion.analysis import FusionSettings

    assert set(FusionSettings().circular_contigs) >= {"chrM", "MT", "chrMT"}
    assert FusionSettings().circular_origin_tolerance_bp == 50
