"""The support gate, asked of every event class, as one predicate.

`status` is a legacy composite whose meaning depends on which branch fired.
For a chain below the threshold it reports EVIDENCE; for a chain above it
reports GEOMETRY (readthrough, circular origin) or COORDINATE RESOLUTION
(supported_rna_fusion, supported_chain_unresolved_coordinates). Geometry is
terminal, so on the whole K562 library the largest chain in the run -- 537
molecules -- is a readthrough and can never carry a supported value, and only
3 of 595 chains carry `supported_rna_fusion`.

`support_status` answers the gate question and nothing else. Codex corrected
an overstatement of mine here: the answer IS recoverable from `status` today,
since `status != "insufficient_full_chain_support"` is the same predicate.
What the field removes is the need to know the legacy value set and its branch
ordering to recover it. The agreement between the two is itself pinned below.
The other two axes already had their own columns and are untouched.

Scope, from Codex's review of the design: this fixes the FILTERING INTERFACE,
not the inference. It admits no chain that was not already admitted -- the
model layer never filtered on chain status in the first place -- and it
publishes no mechanism claim.
"""
import pytest

from fin.fusion.annotation import AnnotatedReadPath, FusionArm
from fin.fusion.chains import build_fusion_chains
from fin.fusion.records import AlignmentBlock


def _arm(chrom, strand, exons, gene, query_start=0):
    ordered = list(exons) if strand == "+" else list(reversed(exons))
    blocks, q = [], query_start
    for start, end in ordered:
        blocks.append(AlignmentBlock(start, end, q, q + (end - start), end - start))
        q += end - start
    return FusionArm(chrom, strand, tuple(blocks), (gene,) if gene else ()), q


def _path(read_id, arms, event_class="interchromosomal_candidate"):
    built, q = [], 0
    for chrom, strand, exons, gene in arms:
        arm, q = _arm(chrom, strand, exons, gene, q)
        built.append(arm)
    return AnnotatedReadPath(read_id, tuple(built), event_class, "resolved")


def _chains(n, event_class, *, min_support=2, offsets=None):
    """n molecules on one geometry, all carrying `event_class`."""
    paths = []
    for index in range(n):
        shift = 0 if offsets is None else offsets[index]
        paths.append(_path(f"r{index}",
                           [("chrA", "+", [(100, 200), (300, 400 + shift)], "GA"),
                            ("chrB", "+", [(900, 1000)], "GB")],
                           event_class))
    return build_fusion_chains(paths, min_full_support=min_support)


EVENT_CLASSES = ["interchromosomal_candidate", "intrachromosomal_candidate",
                 "readthrough_or_splicing_candidate", "circular_reference_origin_spanning",
                 "intragenic_rearrangement_candidate"]


# --------------------------------------------------------------------------
# the invariant: one comparison, every event class
# --------------------------------------------------------------------------


@pytest.mark.parametrize("event_class", EVENT_CLASSES)
@pytest.mark.parametrize("min_support", [1, 2, 3])
@pytest.mark.parametrize("molecules", [1, 2, 3, 4])
def test_the_gate_field_is_exactly_the_count_comparison(event_class, min_support, molecules):
    """Codex's cross-product: below, at and above the threshold, for every
    event class, `support_status` must equal the count comparison and nothing
    else -- not geometry, not coordinates."""
    chains = _chains(molecules, event_class, min_support=min_support)
    assert chains, "fixture produced no chain"
    for chain in chains:
        expected = ("sufficient_full_chain_support"
                    if len(chain.full_support_read_ids) >= min_support
                    else "insufficient_full_chain_support")
        assert chain.support_status == expected


@pytest.mark.parametrize("event_class", EVENT_CLASSES)
def test_a_well_supported_chain_says_so_whatever_its_geometry(event_class):
    """The point of the change: readthrough and circular origin reach
    `sufficient_full_chain_support` exactly like a rearrangement candidate,
    while `status` still reports what it always did."""
    chain = _chains(4, event_class, min_support=2)[0]
    assert chain.support_status == "sufficient_full_chain_support"
    assert len(chain.full_support_read_ids) == 4


def test_the_legacy_status_values_are_unchanged():
    """Frozen: readthrough and circular origin still short-circuit `status`,
    and a supported rearrangement candidate is still `supported_candidate` at
    this layer (analysis.finalize_chains refines it afterwards)."""
    assert _chains(4, "readthrough_or_splicing_candidate")[0].status == "readthrough_or_splicing"
    assert _chains(4, "circular_reference_origin_spanning")[0].status == "circular_reference_origin"
    assert _chains(4, "interchromosomal_candidate")[0].status == "supported_candidate"
    assert _chains(1, "readthrough_or_splicing_candidate")[0].status == "insufficient_full_chain_support"


def test_partial_support_cannot_satisfy_the_gate():
    """A molecule covering only part of the chain is partial evidence. It is
    reported, and it must not move the gate."""
    full = [_path(f"f{i}", [("chrA", "+", [(100, 200), (300, 400)], "GA"),
                            ("chrB", "+", [(900, 1000)], "GB"),
                            ("chrC", "+", [(1500, 1600)], "GC")]) for i in range(1)]
    partial = [_path(f"p{i}", [("chrA", "+", [(100, 200), (300, 400)], "GA"),
                               ("chrB", "+", [(900, 1000)], "GB")]) for i in range(3)]
    chains = build_fusion_chains(full + partial, min_full_support=2)
    for chain in chains:
        assert chain.support_status == (
            "sufficient_full_chain_support" if len(chain.full_support_read_ids) >= 2
            else "insufficient_full_chain_support")
        if chain.partial_support_read_ids:
            assert chain.support_status == "insufficient_full_chain_support"


def test_coordinate_resolution_does_not_touch_the_gate():
    """Two molecules on nearby-but-different placements still clear a
    threshold of 2; their coordinate ambiguity is a different axis."""
    chains = _chains(2, "interchromosomal_candidate", min_support=2, offsets=[0, 3])
    chain = chains[0]
    assert chain.support_status == "sufficient_full_chain_support"
    assert chain.coordinate_status == "unresolved_nearby_modes"


def test_the_field_is_required_not_defaulted():
    """A default would read as a quiet finding about the evidence."""
    import dataclasses

    from fin.fusion.chains import FusionChain

    field = {f.name: f for f in dataclasses.fields(FusionChain)}["support_status"]
    assert field.default is dataclasses.MISSING
    assert field.default_factory is dataclasses.MISSING


def test_no_mechanism_field_is_published():
    """Deliberate: mapping geometry into causal words adds no evidence, cannot
    say which other explanations stay compatible, and this pipeline cannot
    separate trans-splicing from a genomic rearrangement."""
    import dataclasses

    from fin.fusion.chains import FusionChain

    names = {f.name for f in dataclasses.fields(FusionChain)}
    assert not any("mechanism" in name for name in names)


def test_the_new_field_agrees_with_what_the_legacy_status_already_implied():
    """Codex corrected me here: the gate IS recoverable from `status` today --
    `status != "insufficient_full_chain_support"` is the same predicate. The
    field exists so a consumer need not know the legacy value set, not because
    the information was unavailable. Pinning the agreement makes a future
    divergence a test failure rather than a silent contradiction between two
    published columns.
    """
    for event_class in EVENT_CLASSES:
        for molecules in (1, 2, 3):
            for min_support in (1, 2, 3):
                for chain in _chains(molecules, event_class, min_support=min_support):
                    legacy_says_supported = chain.status != "insufficient_full_chain_support"
                    assert legacy_says_supported == (chain.support_status
                                                     == "sufficient_full_chain_support")


def _finalized(molecules, event_class, *, min_support=2, offsets=None):
    """The same fixtures, taken through finalize_chains -- which is what the
    published columns actually carry."""
    from fin.fusion.analysis import FusionReadAudit, FusionSettings, finalize_chains

    paths = []
    for index in range(molecules):
        shift = 0 if offsets is None else offsets[index]
        paths.append(_path(f"r{index}",
                           [("chrA", "+", [(100, 200), (300, 400 + shift)], "GA"),
                            ("chrB", "+", [(900, 1000)], "GB")],
                           event_class))
    audits = {p.read_id: FusionReadAudit(p.read_id, "supported_path", "", path_count=1,
                                         selected=(p,)) for p in paths}
    return finalize_chains(paths, audits, FusionSettings(min_support=min_support,
                                                         breakpoint_tolerance=6))


def test_the_agreement_survives_finalization_not_only_the_builder():
    """Codex's round-2 note: the builder-level cross-product above does not by
    itself pin agreement between the two PUBLISHED columns, because
    finalize_chains rewrites `status` afterwards. Run the same cross-product
    through finalization -- including the unresolved-coordinate branch, which
    the builder-level test never reaches."""
    seen = set()
    for event_class in EVENT_CLASSES:
        for molecules in (1, 2, 3):
            for min_support in (1, 2, 3):
                for offsets in (None, [0, 3, 0][:molecules]):
                    finalized = _finalized(molecules, event_class,
                                           min_support=min_support, offsets=offsets)
                    assert finalized, (event_class, molecules, min_support, offsets)
                    for chain in finalized:
                        seen.add(chain.status)
                        assert chain.status != "supported_candidate"
                        assert ((chain.status != "insufficient_full_chain_support")
                                == (chain.support_status == "sufficient_full_chain_support"))
                        assert chain.support_status == (
                            "sufficient_full_chain_support"
                            if len(chain.full_support_read_ids) >= min_support
                            else "insufficient_full_chain_support")
    # the fixtures really did reach both finalized supported branches, so the
    # assertions above are not vacuous on the unresolved side
    assert {"supported_rna_fusion", "supported_chain_unresolved_coordinates"} <= seen


def test_the_builder_s_intermediate_value_is_replaced_by_the_exact_expected_status():
    """`supported_candidate` is internal. Exact coordinates finalize to
    `supported_rna_fusion`; nearby-mode ambiguity to
    `supported_chain_unresolved_coordinates`. Asserting the replacement, not
    merely that something changed."""
    exact = _finalized(3, "interchromosomal_candidate", min_support=2)
    assert [c.status for c in exact] == ["supported_rna_fusion"]
    assert exact[0].support_status == "sufficient_full_chain_support"

    ambiguous = _finalized(2, "interchromosomal_candidate", min_support=2, offsets=[0, 3])
    assert [c.status for c in ambiguous] == ["supported_chain_unresolved_coordinates"]
    assert ambiguous[0].support_status == "sufficient_full_chain_support"


def test_the_two_published_columns_carry_their_own_values_through_the_writer(tmp_path):
    """Codex round 3: the finalization tests cannot catch a header/value
    misalignment in the serializer -- a transposed pair would still satisfy
    every assertion above. Write the TSV and read it back BY COLUMN NAME.
    """
    import csv

    from fin.fusion.analysis import FusionAnalysis, FusionReadAudit, FusionSettings, finalize_chains
    from fin.fusion.output import write_fusion_outputs

    genome = {"chrA": "ACGT" * 500, "chrB": "TTGCA" * 400}
    paths = ([_path(f"s{i}", [("chrA", "+", [(100, 200), (300, 400)], "GA"),
                              ("chrB", "+", [(900, 1000)], "GB")]) for i in range(3)]
             + [_path("lone", [("chrA", "+", [(100, 200), (300, 700)], "GA"),
                               ("chrB", "+", [(1200, 1300)], "GB")])])
    audits = tuple(FusionReadAudit(p.read_id, "supported_path", "", path_count=1, selected=(p,))
                   for p in paths)
    settings = FusionSettings(min_support=2)
    chains = finalize_chains(paths, {a.read_id: a for a in audits}, settings)
    write_fusion_outputs(FusionAnalysis(chains, audits, {}, settings, "provided"), genome, tmp_path)

    rows = list(csv.DictReader((tmp_path / "fusions.tsv").open(), delimiter="\t"))
    assert len(rows) == len(chains)
    by_id = {c.chain_id: c for c in chains}
    for row in rows:
        chain = by_id[row["chain_id"]]
        assert row["status"] == chain.status
        assert row["chain_support_status"] == chain.support_status
        assert row["status"] != row["chain_support_status"] or \
            chain.status == "insufficient_full_chain_support"
        assert int(row["full_read_support"]) == len(chain.full_support_read_ids)
        assert (row["chain_support_status"] == "sufficient_full_chain_support") == \
            (int(row["full_read_support"]) >= settings.min_support)
    # both sides of the gate really are present in this fixture
    assert {r["chain_support_status"] for r in rows} == {
        "sufficient_full_chain_support", "insufficient_full_chain_support"}
