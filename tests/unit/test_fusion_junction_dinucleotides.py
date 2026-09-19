"""The reference bases at a reported junction, published instead of a boolean.

`fusion_motif` already read the real donor/acceptor dinucleotides, used them
for one membership test against {GT-AG, GC-AG, AT-AC}, and threw them away.
The only thing published was the string "canonical" inside a contrast's
`origin` list, so a reader could not tell GT-AG from AT-AC and learned nothing
at all about a non-canonical junction.

What this annotation is not, and the tests below hold it to:

* Not evidence. A rearrangement can place GT-AG at a breakpoint by
  coincidence, and reading a motif at a coordinate this pipeline chose does
  not corroborate that coordinate.
* Not unbiased, and not independent evidence. `nominate_junctions` scans
  +/-tolerance and nominates a hypothesis wherever the bases spell a canonical
  motif, so canonical placements are GENERATED preferentially. Codex's
  correction: preferential generation is certain, enrichment among the
  placements that finally win is plausible but not established here, so the
  claim is motif-selection bias rather than a demonstrated over-representation.
* Not an adjudication of coordinate ambiguity. It describes the placement the
  row displays; other placements in `coordinate_alternatives` may spell
  something else.

Codex corrected two things I had assumed while designing this: an uncorrected
junction is NOT bias-free (retained means it survived comparison against
motif-enriched alternatives, not that it was never searched), and a model's
junctions are NOT simply its primary chain's corrected edges, because
envelopes aggregate across chains and merged placements can be chosen on
annotation agreement. The annotation is therefore computed from the model's
own reported breakpoints.
"""
import json

import pytest

from fin.fusion.junctions import junction_dinucleotides
from fin.utils.sequences import reverse_complement


# A genome whose bases are known by construction, so the expectations below are
# written down rather than computed by the function under test.
#   chr1 [0,60): positions 20..21 = "GT", 38..39 = "AG"
def _chr1():
    bases = list("C" * 60)
    bases[20:22] = list("GT")
    bases[38:40] = list("AG")
    return "".join(bases)


#   chr2 [0,60): positions 10..11 = "AC", 30..31 = "CT"  (the minus-strand
#   reverse complements are "GT" and "AG"). "CT"/"AG" is deliberately NOT
#   reverse-complement palindromic: Codex found the first version used "AT",
#   which is its own reverse complement and so could not detect an omitted
#   reverse-complement step on the acceptor side.
def _chr2():
    bases = list("T" * 60)
    bases[10:12] = list("AC")
    bases[30:32] = list("CT")
    return "".join(bases)


GENOME = {"chr1": _chr1(), "chr2": _chr2(), "chrN": "AC" + "N" * 8 + "GT", "tiny": "ACGT"}


# --------------------------------------------------------------------------
# the oracle: expectations written by hand, not produced by the helper
# --------------------------------------------------------------------------


def test_plus_to_plus_reads_forward_on_both_sides():
    """Donor = the two bases AFTER the upstream exon ends (chr1[20:22] = GT);
    acceptor = the two BEFORE the downstream exon begins (chr1[38:40] = AG)."""
    assert junction_dinucleotides(GENOME, "chr1", 20, "+", "chr1", 40, "+") == ("GT", "AG")


def test_minus_to_minus_reverse_complements_both_sides():
    """On minus the transcript runs high->low, so the donor window is the two
    bases BELOW the exit and is reverse-complemented: chr2[10:12] = "AC",
    revcomp = "GT". The acceptor window is the two ABOVE the entry:
    chr2[30:32] = "CT", revcomp = "AG" -- not a palindrome, so omitting the
    reverse complement here would be caught."""
    assert junction_dinucleotides(GENOME, "chr2", 12, "-", "chr2", 30, "-") == ("GT", "AG")


def test_plus_donor_with_minus_acceptor():
    assert junction_dinucleotides(GENOME, "chr1", 20, "+", "chr2", 30, "-") == ("GT", "AG")


def test_minus_donor_with_plus_acceptor():
    assert junction_dinucleotides(GENOME, "chr2", 12, "-", "chr1", 40, "+") == ("GT", "AG")


def test_an_interchromosomal_edge_reads_each_side_on_its_own_contig():
    donor, acceptor = junction_dinucleotides(GENOME, "chr1", 20, "+", "chr2", 12, "+")
    assert donor == "GT"                              # chr1[20:22]
    assert acceptor == GENOME["chr2"][10:12] == "AC"  # NOT a canonical pair


def test_a_noncanonical_pair_is_reported_as_it_is():
    """The point of the change: the bases are published whether or not they
    land in the motif set."""
    assert junction_dinucleotides(GENOME, "chr1", 0, "+", "chr1", 10, "+") == ("CC", "CC")


def test_a_non_n_ambiguity_code_is_complemented_not_passed_through():
    """Codex's [P2]: `reverse_complement` maps only A/C/G/T/N and leaves every
    other symbol alone, so R (A or G) on a minus arm came back as R where the
    transcript reads Y (C or T). The membership test that used to be the sole
    consumer rejected either spelling, so publishing the bases is what exposed
    it."""
    genome = {"x": "ARCC", "y": "CCWS"}
    # minus donor at 2 reads x[0:2] = "AR" -> revcomp = "YT"
    assert junction_dinucleotides(genome, "x", 2, "-", "chr1", 40, "+")[0] == "YT"
    # minus acceptor at 2 reads y[2:4] = "WS" -> revcomp = "SW"
    assert junction_dinucleotides(genome, "chr1", 20, "+", "y", 2, "-")[1] == "SW"


@pytest.mark.parametrize("code,complement", [
    ("R", "Y"), ("Y", "R"), ("S", "S"), ("W", "W"), ("K", "M"), ("M", "K"),
    ("B", "V"), ("V", "B"), ("D", "H"), ("H", "D"), ("N", "N"),
])
def test_every_iupac_code_complements_to_its_partner(code, complement):
    genome = {"z": "AA" + code + "A"}          # minus donor at 3 reads z[1:3]
    donor = junction_dinucleotides(genome, "z", 3, "-", "chr1", 40, "+")[0]
    assert donor == complement + "T"           # revcomp("A" + code) = comp(code) + "T"


def test_lower_case_reference_bases_are_upper_cased_as_before():
    """`_base_slice` upper-cases; soft-masked reference must not leak through
    as a differently-spelled dinucleotide."""
    genome = {"soft": "acgtacgt"}
    # donor soft[2:4] = "gt"; acceptor soft[4:6] = "ac"
    assert junction_dinucleotides(genome, "soft", 2, "+", "soft", 6, "+") == ("GT", "AC")


def test_ambiguous_reference_bases_survive_as_N():
    """`N` is what the reference says. Collapsing it into "missing" would lose
    the distinction between an unknown base and an unreadable window."""
    donor, acceptor = junction_dinucleotides(GENOME, "chrN", 4, "+", "chrN", 8, "+")
    assert donor == "NN" and acceptor == "NN"


# --------------------------------------------------------------------------
# missingness is null, never a guess
# --------------------------------------------------------------------------


@pytest.mark.parametrize("chrom,donor,strand", [
    ("absent", 20, "+"),          # no such contig
    ("tiny", 3, "+"),             # window runs off the end
    ("chr1", 1, "-"),             # minus-strand window runs off the start
])
def test_an_unreadable_donor_side_is_none(chrom, donor, strand):
    assert junction_dinucleotides(GENOME, chrom, donor, strand, "chr1", 40, "+")[0] is None


def test_one_unreadable_side_does_not_erase_the_other():
    donor, acceptor = junction_dinucleotides(GENOME, "chr1", 20, "+", "absent", 40, "+")
    assert donor == "GT" and acceptor is None


def test_the_membership_helper_keeps_its_old_empty_string_contract():
    """`fusion_motif` feeds a membership test; ("", "") is not a member. It
    must not start returning None and quietly change that test."""
    from fin.fusion.annotation import FusionArm
    from fin.fusion.junctions import fusion_motif
    from fin.fusion.records import AlignmentBlock

    arm = FusionArm("absent", "+", (AlignmentBlock(0, 10, 0, 10, 10),), ())
    assert fusion_motif(GENOME, arm, arm, (5, 5)) == ("", "")


# --------------------------------------------------------------------------
# what the published column carries
# --------------------------------------------------------------------------


def _writer_fixture(tmp_path, arms_spec):
    from fin.fusion.transcript_models import build_fused_transcript_models
    from fin.fusion.transcript_output import write_transcript_outputs
    from tests.unit.test_fusion_transcript_models_unit import INDEX, _analysis, _path

    genome = {"A": "ACGT" * 250, "B": "TTGCA" * 200, "D": "ACGT" * 250}
    paths = [_path(f"r{i}", arms_spec) for i in range(2)]
    analysis = _analysis({p.read_id: [p] for p in paths}, [p.read_id for p in paths])
    models, collapsed = build_fused_transcript_models(analysis, INDEX, genome)
    write_transcript_outputs(models, collapsed, genome, tmp_path)
    return models, genome


def test_the_column_holds_one_pair_per_inter_arm_edge_in_path_order(tmp_path):
    """Exactly arm_count - 1 entries. Internal splice junctions are junctions
    inside an arm, not between arms, and are not listed."""
    import csv

    models, genome = _writer_fixture(tmp_path, [
        ("A", "+", [(150, 200), (300, 400)], "A"),     # one internal splice
        ("B", "+", [(500, 600)], "B"),
    ])
    rows = list(csv.DictReader((tmp_path / "fusion_models.tsv").open(), delimiter="\t"))
    assert len(rows) == 1
    pairs = json.loads(rows[0]["junction_dinucleotides"])
    assert len(pairs) == int(rows[0]["arm_count"]) - 1 == 1
    assert all(len(pair) == 2 for pair in pairs)


def test_every_published_pair_recomputes_from_the_reference_and_the_row(tmp_path):
    """Codex's whole-output audit, in miniature: independently recompute each
    emitted pair from the genome and the row's own breakpoint coordinates."""
    import csv

    models, genome = _writer_fixture(tmp_path, [
        ("A", "+", [(150, 200), (300, 400)], "A"),
        ("B", "+", [(500, 600)], "B"),
    ])
    rows = {r["model_id"]: r for r in
            csv.DictReader((tmp_path / "fusion_models.tsv").open(), delimiter="\t")}
    for model in models:
        published = json.loads(rows[model.model_id]["junction_dinucleotides"])
        assert len(published) == len(model.breakpoints)
        for (left, donor, left_strand, right, acceptor, right_strand), pair in zip(
                model.breakpoints, published):
            # recomputed here from raw slices, not via junction_dinucleotides
            want_donor = (genome[left][donor:donor + 2] if left_strand == "+"
                          else reverse_complement(genome[left][donor - 2:donor]))
            want_acceptor = (genome[right][acceptor - 2:acceptor] if right_strand == "+"
                             else reverse_complement(genome[right][acceptor:acceptor + 2]))
            assert pair == [want_donor or None, want_acceptor or None]


def test_the_annotation_follows_the_displayed_breakpoints_field(tmp_path):
    """Codex's final-placement trap, rebuilt after he showed my first version
    was vacuous: it used two identical paths, so there was no competing
    placement and a regression reading somewhere else would still have passed.

    Note where the trap has to live. `write_transcript_outputs` never receives
    the chains at all, so "reads the parent chain's coordinates" is not
    reachable from it; what IS reachable is reading the arms with a different
    boundary convention instead of the placement the row displays. So the
    model here has `breakpoints` that DISAGREE with its arms.

    That state is SYNTHETIC and I should not dress it up: Codex checked all
    three construction sites and every one derives breakpoints and arms
    together -- initial construction and aggregation from the same arms, the
    same-sequence merge carrying both from the lead -- so no production merge
    produces this mismatch. An earlier version of this docstring claimed it
    could. The test is a field-authority contract on the writer: whichever
    field the `breakpoints` column prints is the field the annotation must
    read, enforced by making the alternative give a visibly different answer.
    """
    import csv

    from fin.fusion.transcript_models import ArmResolution, FusedTranscriptModel, ModelArm
    from fin.fusion.transcript_output import write_transcript_outputs

    from fin.fusion.transcript_models import _breakpoints_of

    genome = {"chr1": _chr1(), "chr2": _chr2()}
    resolution = ArmResolution("first", "arm_novel", (), (), (), (), (), ())
    # the arms' own boundaries are chr1 exit 42 and chr2 entry 30, which spell
    # CC / TT; the displayed placement is chr1 20 and chr2 12, which spell
    # GT / AC. Deriving the annotation from the arms instead of the displayed
    # field therefore produces a visibly different answer.
    arms = (ModelArm("chr1", "+", ((10, 42),), resolution),
            ModelArm("chr2", "+", ((30, 50),), resolution))
    displayed = (("chr1", 20, "+", "chr2", 12, "+"),)
    assert _breakpoints_of(arms) != displayed, "the trap needs the two to disagree"

    model = FusedTranscriptModel(
        "FTM_trap", arms, displayed, "fusion_chain_x", ("fusion_chain_x",),
        ("r0", "r1"), (), (), (1,), "exact_observed", "interchromosomal_candidate",
        "CAG_x", "annotated", "annotated",
    )
    write_transcript_outputs([model], [], genome, tmp_path)
    row = next(csv.DictReader((tmp_path / "fusion_models.tsv").open(), delimiter="\t"))

    assert row["breakpoints"] == "chr1:20:+>chr2:12:+"
    published = json.loads(row["junction_dinucleotides"])
    assert published == [["GT", "AC"]]                     # chr1[20:22], chr2[10:12]
    # and NOT what the arms' own boundaries would have spelled
    assert published != [["CC", "TT"]]


def test_the_evidence_record_carries_the_same_pairs(tmp_path):
    import csv

    models, genome = _writer_fixture(tmp_path, [
        ("A", "+", [(150, 200), (300, 400)], "A"),
        ("B", "+", [(500, 600)], "B"),
    ])
    rows = list(csv.DictReader((tmp_path / "fusion_models.tsv").open(), delimiter="\t"))
    evidence = [json.loads(line) for line in (tmp_path / "fusion_model_evidence.jsonl").open()]
    assert evidence[0]["junction_dinucleotides"] == json.loads(rows[0]["junction_dinucleotides"])


def test_no_mechanism_or_canonical_boolean_was_added(tmp_path):
    """Deliberate, on Codex's argument: the pair strictly contains more
    sequence information than set membership, and a boolean beside it would
    invite reading the annotation as a mechanism call."""
    import csv

    _writer_fixture(tmp_path, [("A", "+", [(150, 200), (300, 400)], "A"),
                               ("B", "+", [(500, 600)], "B")])
    header = (tmp_path / "fusion_models.tsv").open().readline().split("\t")
    assert not any("canonical" in name or "mechanism" in name for name in header)
