"""TSS evidence: hazard geometry, degradation null, three-way verdicts."""
import random

import pytest

from fin.analysis.tss_evidence import (
    VERDICT_SUPPORTED,
    VERDICT_UNIDENTIFIABLE,
    VERDICT_UNSUPPORTED,
    build_hazard_profile,
    classify_identifiability,
    genomic_to_offset,
    pooled_background_hazard,
    read_five_prime_offset,
    spliced_length,
    evaluate_internal_tss,
)

EXONS = [(1000, 1200), (1500, 1700), (2000, 2300)]   # spliced length 700


def test_spliced_length_and_offsets_plus_strand():
    assert spliced_length(EXONS) == 700
    assert genomic_to_offset(1000, EXONS, "+") == 0
    assert genomic_to_offset(1199, EXONS, "+") == 199
    assert genomic_to_offset(1500, EXONS, "+") == 200      # intron skipped
    assert genomic_to_offset(2000, EXONS, "+") == 400
    assert genomic_to_offset(1300, EXONS, "+") is None     # intronic -> abstain


def test_offsets_are_mirrored_on_minus_strand():
    # '-' strand: transcript 5' end is the genomic right edge.
    assert genomic_to_offset(2300, EXONS, "-") == 0
    assert genomic_to_offset(1000, EXONS, "-") == 700
    assert genomic_to_offset(2000, EXONS, "-") == 300


def test_read_five_prime_offset_is_strand_aware():
    assert read_five_prime_offset((1000, 2300), EXONS, "+") == 0
    assert read_five_prime_offset((1000, 2300), EXONS, "-") == 0
    # a '+'-strand read starting mid-transcript is 5'-truncated
    assert read_five_prime_offset((2000, 2300), EXONS, "+") == 400
    # overhang beyond the model's 5' terminus clamps to 0, not None
    assert read_five_prime_offset((900, 2300), EXONS, "+") == 0


def test_hazard_uses_reached_reads_as_denominator():
    # 5' offsets: three reads stop at 0, two at 100, five at 400
    offsets = [0, 0, 0, 100, 100] + [400] * 5
    prof = build_hazard_profile(offsets, tx_len=700, bin_bp=25)
    # bin 16 covers 400..424: 5 terminations, all 10 reads reached it
    assert prof.ends[16] == 5
    assert prof.at_risk[16] == 10
    assert prof.hazard(16) == pytest.approx(0.5)
    # bin 0 covers 0..24: 3 terminations, and only those 3 reached it
    assert prof.ends[0] == 3
    assert prof.at_risk[0] == 3
    assert prof.hazard(0) == pytest.approx(1.0)
    assert prof.hazard(99) is None      # never reached -> no hazard


def _degradation_only(n, tx_len, hazard, rng, tss_offset=0):
    """Reads from ONE transcript: geometric termination walking 3'->5'."""
    out = []
    for _ in range(n):
        d = tx_len
        while d > tss_offset:
            if rng.random() < hazard:
                break
            d -= 1
        out.append(max(d, tss_offset))
    return out


def test_degradation_only_locus_is_not_called_supported():
    """The critical false-positive control: smooth decay must not fire."""
    rng = random.Random(1)
    offsets = _degradation_only(400, 700, 0.004, rng)
    bg = pooled_background_hazard(
        [build_hazard_profile(offsets, 700, bin_bp=25)])
    fired = 0
    for d0 in range(100, 600, 25):
        ev = evaluate_internal_tss(
            candidate_id=f"s{d0}", parent_id="L", offsets=offsets,
            tss_offset=d0, background_hazard=bg, n_bootstrap=500,
        )
        if ev.verdict == VERDICT_SUPPORTED:
            fired += 1
    # a handful of bins out of 20 may pass at alpha=0.05 by construction;
    # the requirement is that degradation alone does not systematically fire
    assert fired <= 3, f"degradation-only locus fired {fired} times"


def test_genuine_internal_tss_is_supported():
    rng = random.Random(2)
    long_reads = _degradation_only(200, 700, 0.004, rng)
    # a real short isoform: its own population starts exactly at offset 300
    short_reads = _degradation_only(120, 700, 0.004, rng, tss_offset=300)
    offsets = long_reads + short_reads
    bg = pooled_background_hazard(
        [build_hazard_profile(long_reads, 700, bin_bp=25)])
    ev = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=offsets, tss_offset=300,
        background_hazard=bg, n_bootstrap=1000,
    )
    assert ev.verdict == VERDICT_SUPPORTED, ev.reason
    assert ev.n_peak >= 3
    assert ev.effect_size > 0
    assert ev.mixture_pi > 0


def test_low_depth_abstains_rather_than_guessing():
    ev = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=[300, 300, 300],
        tss_offset=300, background_hazard=0.01, min_at_risk=10,
    )
    assert ev.verdict == VERDICT_UNIDENTIFIABLE
    assert "insufficient_depth" in ev.reason


def test_dominant_start_without_upstream_reads_is_supported():
    """All reads starting at d0 with nothing upstream supports the short model.

    Degradation is smooth and requires the parent to be present; it cannot
    deposit every read into one bin while leaving the parent's 5' region
    completely unobserved. Measured on SIRV this pattern (e.g. 49/49 reads)
    marks real nested transcripts, so calling it unidentifiable was wrong.
    """
    offsets = [300] * 40
    ev = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=offsets, tss_offset=300,
        background_hazard=0.01,
    )
    assert ev.verdict == VERDICT_SUPPORTED
    assert ev.reason == "dominant_start_no_upstream_parent_unobserved"
    assert ev.n_upstream_of_tss == 0


def test_diffuse_start_without_upstream_reads_abstains():
    """A weak/diffuse start with no upstream evidence stays unidentifiable."""
    offsets = [300, 340, 380, 420, 460, 500, 540, 580, 620, 660, 700, 740]
    ev = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=offsets, tss_offset=300,
        background_hazard=0.01, min_at_risk=1,
    )
    assert ev.verdict == VERDICT_UNIDENTIFIABLE
    assert "no_upstream_reads" in ev.reason


def test_missing_background_model_abstains():
    """With no pooled AND no estimable local null the test must abstain."""
    ev = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=[0] * 20 + [300] * 20,
        tss_offset=300, background_hazard=0.0, use_local_background=False,
    )
    assert ev.verdict == VERDICT_UNIDENTIFIABLE
    assert ev.reason == "no_background_model"


def test_local_background_absorbs_a_degradation_hotspot():
    """A bin that is high but sits in an equally-high neighbourhood is not a TSS."""
    from fin.analysis.tss_evidence import local_background_hazard
    rng = random.Random(9)
    # uniformly elevated termination across a broad region -> local median high
    offsets = [0] * 40 + [rng.randrange(400, 900) for _ in range(300)]
    local = local_background_hazard(offsets, 600, bin_bp=25,
                                    neighbourhood_bins=12)
    assert local is not None and local > 0.0
    ev = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=offsets, tss_offset=600,
        background_hazard=1e-4, n_bootstrap=400,
    )
    # the local null must be far above the tiny pooled value it was given
    assert ev.background_hazard >= local


def test_identifiability_ladder():
    long_ex = [(1000, 1200), (1500, 1700), (2000, 2300)]
    # same chain, same TES, internal TSS only -> hardest rung
    tss_only = [(1600, 1700), (2000, 2300)]
    assert classify_identifiability(tss_only, long_ex, "+") == "tss_only"
    # own TES: 3' end differs (reliable end in dRNA)
    own_tes = [(1000, 1200), (1500, 1700), (2000, 2100)]
    assert classify_identifiability(own_tes, long_ex, "+") == "own_tes"
    # own junction: a junction the parent does not have
    own_j = [(1000, 1200), (1600, 1700)]
    assert classify_identifiability(own_j, long_ex, "+") == "own_junction"


def test_evidence_row_schema_is_stable():
    ev = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=[0] * 20 + [300] * 20,
        tss_offset=300, background_hazard=0.01, n_bootstrap=200,
    )
    row = ev.as_row()
    for key in ("candidate_id", "parent_id", "tss_offset", "identifiability",
                "verdict", "reason", "n_reads_locus", "n_at_risk", "n_peak",
                "n_upstream_of_tss", "peak_fraction", "background_hazard",
                "expected_peak", "effect_size", "mixture_pi", "llr",
                "p_value"):
        assert key in row
    assert row["verdict"] in (VERDICT_SUPPORTED, VERDICT_UNSUPPORTED,
                              VERDICT_UNIDENTIFIABLE)


def test_config_validates_tss_mode_and_requires_endpoint_refine(tmp_path):
    from fin.pipeline.config import PipelineConfig
    bam = tmp_path / "x.bam"
    bam.touch()
    bad = PipelineConfig(bam_path=str(bam), tss_evidence_mode="bogus",
                         endpoint_refine=True, post_selection_refit=True)
    with pytest.raises(ValueError, match="tss_evidence_mode"):
        bad.validate()
    orphan = PipelineConfig(bam_path=str(bam), tss_evidence_mode="audit")
    with pytest.raises(ValueError, match="requires --endpoint-refine"):
        orphan.validate()
    ok = PipelineConfig(bam_path=str(bam), tss_evidence_mode="audit",
                        endpoint_refine=True, post_selection_refit=True)
    ok.validate()
    assert ok.tss_evidence_mode == "audit"


def test_cli_exposes_tss_evidence_mode():
    from click.testing import CliRunner
    from fin.cli import main
    out = CliRunner().invoke(main, ["--help"]).output
    assert "--tss-evidence-mode" in out
    bad = CliRunner().invoke(main, ["--tss-evidence-mode", "bogus"])
    assert bad.exit_code != 0


def _split_fixture(n_alt: int):
    """A parent with reads at its own TSS plus a block at an internal start."""
    from fin.analysis.quantification import QuantResult

    # NOTE: plan_endpoint_splits only re-cuts the FIRST exon, so an
    # alternative TSS must fall inside it; a start in a downstream exon would
    # delete whole exons and is a chain change, not an endpoint change.
    exons = ((1000, 1300), (1500, 3000))
    ends, ids = {}, []
    for i in range(30):                       # parent-start reads
        rid = f"p{i}"
        ends[rid] = (1000, 3000)
        ids.append(rid)
    for i in range(n_alt):                    # alternative-start reads
        rid = f"a{i}"
        ends[rid] = (1150, 3000)
        ids.append(rid)
    qr = QuantResult(
        candidate_id="L", abundance=float(len(ids)), confidence=1.0,
        num_assigned_reads=len(ids), source="novel", chrom="chr1",
        strand="+", start=1000, end=3000, exons=exons,
        assigned_read_ids=tuple(ids),
    )
    return {"L": qr}, ends


def test_audit_mode_records_verdicts_without_changing_the_split():
    from fin.analysis.endpoint_refine import plan_endpoint_splits

    results, ends = _split_fixture(20)
    base = plan_endpoint_splits(results, ends)
    audited = plan_endpoint_splits(results, ends, tss_evidence_mode="audit",
                                   background_hazard=0.01)
    assert [s.candidate_id for s in base.replacements.get("L", [])] == \
           [s.candidate_id for s in audited.replacements.get("L", [])]
    assert audited.tss_verdicts.get("L")
    v = audited.tss_verdicts["L"][0]
    assert v["verdict"] in (VERDICT_SUPPORTED, VERDICT_UNSUPPORTED,
                            VERDICT_UNIDENTIFIABLE)


def test_require_mode_keeps_every_abstained_state():
    """Invariant: an `unidentifiable` verdict never deletes an endpoint state.

    Insufficient evidence is not evidence of absence -- only an explicitly
    `unsupported` start may be removed by `require`.
    """
    from fin.analysis.endpoint_refine import plan_endpoint_splits

    for n_alt in (8, 12, 20, 40):
        results, ends = _split_fixture(n_alt)
        plan = plan_endpoint_splits(results, ends,
                                    tss_evidence_mode="require",
                                    background_hazard=0.01)
        kept = {s.candidate_id for s in plan.replacements.get("L", [])}
        for v in plan.tss_verdicts.get("L", []):
            if v["verdict"] == VERDICT_UNIDENTIFIABLE:
                assert v["candidate_id"] in kept, (
                    f"abstained state {v['candidate_id']} was dropped")


def test_require_mode_can_remove_an_unsupported_start():
    """A start that is refuted by the degradation null may be removed."""
    from fin.analysis.endpoint_refine import plan_endpoint_splits

    results, ends = _split_fixture(12)
    # a very permissive null (high background) refutes a modest peak
    plan = plan_endpoint_splits(results, ends, tss_evidence_mode="require",
                                background_hazard=0.95)
    verdicts = plan.tss_verdicts.get("L", [])
    assert verdicts, "expected the TSS test to run"
    unsupported = [v for v in verdicts if v["verdict"] == VERDICT_UNSUPPORTED]
    kept = {s.candidate_id for s in plan.replacements.get("L", [])}
    for v in unsupported:
        assert v["candidate_id"] not in kept


def test_broad_peak_is_rejected_as_too_wide():
    """A wide degradation hotspot must not pass as a sharp TSS."""
    rng = random.Random(11)
    # 60 reads spread across the whole 25bp bin -> large MAD
    offsets = [0] * 60 + [300 + rng.randrange(0, 25) for _ in range(60)]
    ev = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=offsets, tss_offset=300,
        background_hazard=0.005, n_bootstrap=400, max_peak_mad=2.0,
    )
    assert ev.verdict == VERDICT_UNSUPPORTED
    assert "peak_too_broad" in ev.reason
    assert ev.peak_mad > 2.0


def test_parent_unobserved_flag_is_recorded():
    ev = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=[300] * 40, tss_offset=300,
        background_hazard=0.01,
    )
    assert ev.parent_unobserved is True
    assert ev.as_row()["parent_unobserved"] is True


def test_determinism_same_seed_same_pvalue():
    kw = dict(candidate_id="S", parent_id="L",
              offsets=[0] * 30 + [300] * 15, tss_offset=300,
              background_hazard=0.01, n_bootstrap=500)
    assert evaluate_internal_tss(**kw).p_value == evaluate_internal_tss(**kw).p_value


def test_evidence_rows_are_strict_json_serialisable():
    """Diagnostics land in endpoint_refine.json; NaN is not valid JSON."""
    import json

    rows = []
    # structural 3'-end call (no calibrated p)
    from fin.analysis.tss_evidence import evaluate_tes_support
    rows.append(evaluate_tes_support(
        candidate_id="s", parent_id="L",
        read_three_prime_offsets=[100] * 20 + [900] * 20,
        candidate_tes_offset=900, parent_tes_offset=100).as_row())
    # abstention
    rows.append(evaluate_internal_tss(
        candidate_id="s", parent_id="L", offsets=[1, 2],
        tss_offset=1, background_hazard=0.01).as_row())
    # hazard-tested call
    rows.append(evaluate_internal_tss(
        candidate_id="s", parent_id="L",
        offsets=list(range(0, 400, 4)) + [300] * 40,
        tss_offset=300, background_hazard=0.01).as_row())
    blob = json.dumps(rows, allow_nan=False)   # must not raise
    back = json.loads(blob)
    assert len(back) == 3
    for r in back:
        assert r["p_value"] is None or isinstance(r["p_value"], float)


def test_broad_peak_rejected_even_with_no_upstream_reads():
    """Regression: the n_up == 0 shortcut must not bypass the width guard."""
    rng = random.Random(5)
    # every read starts inside one bin but SPREAD across it, nothing upstream
    offsets = [300 + rng.randrange(0, 25) for _ in range(60)]
    ev = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=offsets, tss_offset=300,
        background_hazard=0.005, max_peak_mad=2.0,
    )
    assert ev.n_upstream_of_tss == 0
    assert ev.parent_unobserved is True
    assert ev.verdict == VERDICT_UNSUPPORTED
    assert "peak_too_broad" in ev.reason


def test_bonferroni_selection_correction_scales_with_search_space():
    """A peak found by scanning many bins must pay for the search."""
    offsets = list(range(0, 500, 5)) + [300] * 12
    one = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=offsets, tss_offset=300,
        background_hazard=0.02, n_eligible_bins=1, n_bootstrap=2000)
    many = evaluate_internal_tss(
        candidate_id="S", parent_id="L", offsets=offsets, tss_offset=300,
        background_hazard=0.02, n_eligible_bins=40, n_bootstrap=2000)
    assert many.p_value >= one.p_value
    assert many.p_value <= 1.0
    # exactly Bonferroni, not a simulated max-statistic
    assert many.p_value == pytest.approx(min(1.0, one.p_value * 40), rel=1e-9)
    # a single-bin search space is still a corrected (n=1) test
    assert one.selection_corrected is True
    assert many.selection_corrected is True


def test_locus_fdr_combines_alternatives_then_corrects_across_loci():
    """Within a locus: Bonferroni. Across loci: BY. One vote per locus."""
    from fin.analysis.tss_evidence import TssEvidence, apply_grouped_fdr

    def _ev(cid, parent, p):
        e = TssEvidence(candidate_id=cid, parent_id=parent, tss_offset=1,
                        identifiability="tss_only",
                        verdict=VERDICT_SUPPORTED, reason="x")
        e.p_value, e.selection_corrected = p, True
        return e

    # locus A has 3 alternatives; without locus grouping it would get 3 votes
    evs = [_ev("a1", "A", 0.01), _ev("a2", "A", 0.02), _ev("a3", "A", 0.03),
           _ev("b1", "B", 0.9)]
    apply_grouped_fdr(evs, alpha=0.05)
    # every row of locus A shares ONE locus q-value
    qs = {e.q_value for e in evs if e.parent_id == "A"}
    assert len(qs) == 1
    # BY over 2 loci with C(2)=1.5 is strictly harsher than raw p
    assert evs[0].q_value > 0.01


def test_structural_calls_never_enter_the_locus_fdr_family():
    from fin.analysis.tss_evidence import (
        TssEvidence, apply_grouped_fdr, evaluate_tes_support,
    )

    structural = evaluate_tes_support(
        candidate_id="s", parent_id="L",
        read_three_prime_offsets=[100] * 20 + [900] * 20,
        candidate_tes_offset=900, parent_tes_offset=100)
    assert structural.p_value is None
    hazard = TssEvidence(candidate_id="h", parent_id="L2", tss_offset=1,
                         identifiability="tss_only",
                         verdict=VERDICT_SUPPORTED, reason="x")
    hazard.p_value, hazard.selection_corrected = 1e-9, True
    apply_grouped_fdr([structural, hazard], alpha=0.05)
    assert structural.q_value is None
    assert structural.verdict == VERDICT_SUPPORTED   # untouched by FDR
    assert hazard.q_value is not None


def test_strong_alternative_does_not_certify_its_weak_sibling():
    """A parent q-value must not carry weak alternatives through.

    Regression: propagating only the parent-level q meant one strong start
    certified every other start on the same parent.
    """
    from fin.analysis.tss_evidence import TssEvidence, apply_grouped_fdr

    def _ev(cid, p):
        e = TssEvidence(candidate_id=cid, parent_id="SAME_PARENT",
                        tss_offset=1, identifiability="tss_only",
                        verdict=VERDICT_SUPPORTED, reason="x")
        e.p_value, e.selection_corrected = p, True
        return e

    strong, weak = _ev("strong", 1e-9), _ev("weak", 0.40)
    apply_grouped_fdr([strong, weak], alpha=0.05)
    assert strong.verdict == VERDICT_SUPPORTED
    assert weak.verdict == VERDICT_UNSUPPORTED
    assert "p_within" in weak.reason
    # both see the SAME parent q, but only the strong one clears p_within
    assert strong.q_value == weak.q_value
    assert strong.p_within <= 0.05 < weak.p_within
    assert "p_within" in weak.as_row()
