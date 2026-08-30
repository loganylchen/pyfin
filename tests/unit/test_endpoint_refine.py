"""EndpointRefine: strand-aware splits, guards, routing, requantification."""
import pytest

from fin.analysis.abundance_refit import (
    ResponsibilityLedger,
    refit_survivor_abundance,
)
from fin.analysis.endpoint_refine import (
    apply_endpoint_splits,
    plan_endpoint_splits,
)
from fin.analysis.quantification import QuantResult
from fin.pipeline.config import PipelineConfig


def _qr(cid, exons, *, strand="+", reads=(), source="novel"):
    exons = tuple(exons)
    return QuantResult(
        candidate_id=cid, abundance=float(len(reads)), confidence=1.0,
        num_assigned_reads=len(reads), source=source, chrom="chr1",
        strand=strand, start=exons[0][0], end=exons[-1][1], exons=exons,
        assigned_read_ids=tuple(reads),
    )


def _ends(groups):
    """groups: [(prefix, n, (start, end))] -> read_ends map + read ids."""
    ends, ids = {}, []
    for prefix, n, span in groups:
        for i in range(n):
            rid = f"{prefix}{i}"
            ends[rid] = span
            ids.append(rid)
    return ends, ids


BASE_EXONS = [(1000, 1300), (1500, 2000)]


def test_plus_strand_split_on_distinct_tes_modes():
    ends, ids = _ends([
        ("a", 5, (1000, 2000)),   # primary state
        ("b", 5, (1000, 2600)),   # longer TES mode (exterior, genuine)
    ])
    qr = _qr("x", [(1000, 1300), (1500, 2600)], reads=ids)
    plan = plan_endpoint_splits({"x": qr}, ends)
    assert "x" in plan.replacements
    subs = plan.replacements["x"]
    assert len(subs) == 2
    spans = sorted((s.start, s.end) for s in subs)
    assert spans == [(1000, 2000), (1000, 2600)]
    # every routed read points at one of the new states
    assert set(plan.read_routes["x"].values()) == {s.candidate_id for s in subs}


def test_minus_strand_orientation():
    # '-' strand: TSS = span end, TES = span start.
    ends, ids = _ends([
        ("a", 5, (1000, 2000)),
        ("b", 5, (400, 2000)),    # alternative TES (exterior on '-')
    ])
    qr = _qr("m", [(1000, 1300), (1500, 2000)], strand="-", reads=ids)
    plan = plan_endpoint_splits({"m": qr}, ends)
    assert "m" in plan.replacements
    spans = sorted((s.start, s.end) for s in plan.replacements["m"])
    assert spans == [(400, 2000), (1000, 2000)]


def test_unsupported_pair_is_rejected():
    ends, ids = _ends([
        ("a", 8, (1000, 2000)),
        ("b", 2, (1000, 2600)),   # below min_end_reads=3
    ])
    qr = _qr("x", BASE_EXONS, reads=ids)
    plan = plan_endpoint_splits({"x": qr}, ends)
    assert not plan.replacements


def test_degradation_interior_tss_needs_double_support():
    # '+': interior (downstream) 5' modes are the truncation direction.
    ends_weak, ids_weak = _ends([
        ("a", 6, (1000, 2000)),
        ("b", 4, (1210, 2000)),   # interior TSS, support 4 < 2*3
    ])
    qr = _qr("x", BASE_EXONS, reads=ids_weak)
    assert not plan_endpoint_splits({"x": qr}, ends_weak).replacements

    ends_strong, ids_strong = _ends([
        ("a", 6, (1000, 2000)),
        ("b", 6, (1210, 2000)),   # interior TSS with 2x support -> allowed
    ])
    qr2 = _qr("x", BASE_EXONS, reads=ids_strong)
    assert "x" in plan_endpoint_splits({"x": qr2}, ends_strong).replacements


def test_split_cap_and_determinism():
    ends, ids = _ends([
        ("a", 6, (1000, 2000)),
        ("b", 5, (1000, 2600)),
        ("c", 4, (700, 2000)),    # third state - beyond max_splits=2
    ])
    qr = _qr("x", [(1000, 1300), (1500, 2600)], reads=ids)
    p1 = plan_endpoint_splits({"x": qr}, ends)
    p2 = plan_endpoint_splits({"x": qr}, ends)
    assert len(p1.replacements["x"]) == 2
    assert [s.candidate_id for s in p1.replacements["x"]] == \
           [s.candidate_id for s in p2.replacements["x"]]
    assert all(s.candidate_id.startswith("novel_") for s in p1.replacements["x"])


def test_gtf_and_mono_are_exempt():
    ends, ids = _ends([("a", 5, (1000, 2000)), ("b", 5, (1000, 2600))])
    gtf = _qr("g", BASE_EXONS, reads=ids, source="gtf")
    mono = _qr("m", [(1000, 2000)], reads=ids)
    plan = plan_endpoint_splits({"g": gtf, "m": mono}, ends)
    assert not plan.replacements


def test_refit_routes_conserve_mass_and_requantify():
    ends, ids = _ends([
        ("a", 5, (1000, 2000)),
        ("b", 5, (1000, 2600)),
    ])
    qr = _qr("x", [(1000, 1300), (1500, 2600)], reads=ids)
    plan = plan_endpoint_splits({"x": qr}, ends)
    split = apply_endpoint_splits({"x": qr}, plan)
    assert set(split) == {s.candidate_id for s in plan.replacements["x"]}

    ledger = ResponsibilityLedger(
        weights={rid: {"x": 1.0} for rid in ids},
        input_read_ids=tuple(ids),
    )
    refitted, diag = refit_survivor_abundance(
        split, [ledger],
        split_routes=plan.read_routes, split_primary=plan.primary,
    )
    assert diag["selection_orphaned_reads"] == 0
    assert diag["mass_balance_error"] < 1e-9
    total = sum(r.abundance for r in refitted.values())
    assert total == pytest.approx(len(ids))
    by_span = {(r.start, r.end): r.abundance for r in refitted.values()}
    assert by_span[(1000, 2000)] == pytest.approx(5.0)
    assert by_span[(1000, 2600)] == pytest.approx(5.0)


def test_idempotence_after_split():
    ends, ids = _ends([
        ("a", 5, (1000, 2000)),
        ("b", 5, (1000, 2600)),
    ])
    qr = _qr("x", [(1000, 1300), (1500, 2600)], reads=ids)
    plan = plan_endpoint_splits({"x": qr}, ends)
    split = apply_endpoint_splits({"x": qr}, plan)
    # re-planning on the split states with per-state reads yields no new split
    for sub in split.values():
        sub_plan = plan_endpoint_splits({sub.candidate_id: sub}, ends)
        assert not sub_plan.replacements


def test_diagnostic_details_expose_route_audit():
    """Details must reconcile the split against what the refit will route."""
    ends, ids = _ends([
        ("a", 5, (1000, 2000)),
        ("b", 5, (1000, 2600)),
        ("c", 2, (1000, 3300)),   # unsupported pair -> unrouted
    ])
    qr = _qr("x", [(1000, 1300), (1500, 3300)], reads=ids)
    plan = plan_endpoint_splits({"x": qr}, ends)
    assert len(plan.details) == 1
    d = plan.details[0]
    for key in ("candidate_id", "n_end_reads", "n_assigned_reads", "tss_modes",
                "tes_modes", "polya_available", "states", "routed_reads",
                "route_counts_by_state", "unrouted_to_primary",
                "primary_state"):
        assert key in d, key
    assert sum(d["route_counts_by_state"].values()) == d["routed_reads"]
    assert d["unrouted_to_primary"] == d["n_assigned_reads"] - d["routed_reads"]
    assert d["unrouted_to_primary"] >= 0
    assert set(d["route_counts_by_state"]) <= {
        s.candidate_id for s in plan.replacements["x"]}
    assert d["primary_state"] == plan.primary["x"]


def test_unrouted_reads_follow_primary_in_refit():
    """unrouted_to_primary is not a leak: refit sends those reads to primary."""
    from fin.analysis.abundance_refit import (
        ResponsibilityLedger,
        refit_survivor_abundance,
    )
    ends, ids = _ends([
        ("a", 5, (1000, 2000)),
        ("b", 5, (1000, 2600)),
        ("c", 2, (1000, 3300)),
    ])
    qr = _qr("x", [(1000, 1300), (1500, 3300)], reads=ids)
    plan = plan_endpoint_splits({"x": qr}, ends)
    split = apply_endpoint_splits({"x": qr}, plan)
    ledger = ResponsibilityLedger(
        weights={rid: {"x": 1.0} for rid in ids},
        input_read_ids=tuple(ids),
    )
    refitted, diag = refit_survivor_abundance(
        split, [ledger], split_routes=plan.read_routes,
        split_primary=plan.primary,
    )
    assert diag["selection_orphaned_reads"] == 0
    assert diag["mass_balance_error"] < 1e-9
    assert sum(r.abundance for r in refitted.values()) == pytest.approx(len(ids))


def test_config_requires_refit_for_endpoint_refine(tmp_path):
    bam = tmp_path / "x.bam"
    bam.touch()
    bad = PipelineConfig(bam_path=str(bam), endpoint_refine=True)
    with pytest.raises(ValueError):
        bad.validate()
    ok = PipelineConfig(bam_path=str(bam), endpoint_refine=True,
                        post_selection_refit=True)
    ok.validate()
    assert ok.post_selection_refit_effective is True


# --- TSS-evidence rung routing (production must match the study) ----------

def _verdict_rows(plan):
    return [r for rows in plan.tss_verdicts.values() for r in rows]


def test_rung_routing_shared_tss_own_tes_uses_three_prime_rule():
    """Shared start + own 3' end -> decided on the RELIABLE end, not a hazard.

    There is no internal TSS to detect here, so the verdict must come from
    `evaluate_tes_support` and must NOT carry a hazard p-value.
    """
    ends, ids = _ends([
        ("a", 30, (1000, 2000)),   # primary: start 1000, TES 2000
        ("b", 30, (1000, 2600)),   # same start, different TES
    ])
    qr = _qr("x", [(1000, 1300), (1500, 2600)], reads=ids)
    plan = plan_endpoint_splits({"x": qr}, ends, tss_evidence_mode="audit",
                                background_hazard=0.01)
    rows = _verdict_rows(plan)
    assert rows, "expected a verdict for the alternative state"
    r = rows[0]
    # EXACT route: the 3'-end rule, not the 5' hazard test.
    assert r["identifiability"] == "own_tes"
    assert r["reason"] == "distinct_3prime_cluster_on_the_reliable_end"
    assert r["verdict"] == "supported"
    # structural rule -> NO calibrated probability, excluded from BH
    assert r["p_value"] is None
    assert r["q_value"] is None
    assert r["selection_corrected"] is False
    # the 3'-end rule counts ALL reads at the locus, none discarded
    assert r["n_at_risk"] == 60
    assert r["n_peak"] == 30      # the 30 reads ending at the candidate TES


def test_rung_routing_shared_tes_distinct_tss_uses_hazard_test():
    """Same 3' end, different start -> the hard `tss_only` rung."""
    ends, ids = _ends([
        ("a", 40, (1000, 2000)),   # start at the model's 5' end
        ("b", 40, (1200, 2000)),   # internal start, SAME TES
    ])
    qr = _qr("x", [(1000, 1300), (1500, 2000)], reads=ids)
    plan = plan_endpoint_splits({"x": qr}, ends, tss_evidence_mode="audit",
                                background_hazard=0.005)
    rows = _verdict_rows(plan)
    assert rows
    tss = [r for r in rows if r["identifiability"] == "tss_only"]
    assert tss, "shared TES + distinct start must take the hazard route"
    r = tss[0]
    # hazard route -> a calibrated, selection-corrected probability EXISTS
    assert isinstance(r["p_value"], float)
    assert 0.0 < r["p_value"] <= 1.0
    assert r["selection_corrected"] is True
    # same TES -> no read is filtered out, all 80 remain at risk
    assert r["n_at_risk"] == 80


def test_rung_routing_own_tes_and_distinct_tss_partitions_reads_first():
    """Different start AND different 3' end -> reads partitioned by 3' end."""
    # NOTE: the alternative start must lie in the FIRST exon (1000-1300);
    # plan_endpoint_splits only re-cuts exon 1 by design.
    ends, ids = _ends([
        ("a", 30, (1000, 2000)),
        ("b", 30, (1200, 2600)),   # own start AND own TES
    ])
    qr = _qr("x", [(1000, 1300), (1500, 2600)], reads=ids)
    plan = plan_endpoint_splits({"x": qr}, ends, tss_evidence_mode="audit",
                                background_hazard=0.005)
    rows = _verdict_rows(plan)
    assert rows
    r = rows[0]
    assert r["identifiability"] == "own_tes"
    # Distinct start AND distinct TES -> reads are PARTITIONED by 3' end
    # first: only the 30 reads ending at 2600 may vote on this start.
    assert r["n_at_risk"] == 30, f"expected the 30 TES-matched reads, got {r}"
    assert r["n_peak"] == 30
    assert r["reason"] == "dominant_start_no_upstream_parent_unobserved"
    assert r["verdict"] == "supported"
    assert r["parent_unobserved"] is True
    # parent-unobserved is a STRUCTURAL rule -> no calibrated probability
    assert r["p_value"] is None
    assert r["q_value"] is None


def test_audit_mode_never_changes_the_split_but_require_can():
    ends, ids = _ends([
        ("a", 30, (1000, 2000)),
        ("b", 30, (1000, 2600)),
    ])
    qr = _qr("x", [(1000, 1300), (1500, 2600)], reads=ids)
    off = plan_endpoint_splits({"x": qr}, ends)
    aud = plan_endpoint_splits({"x": qr}, ends, tss_evidence_mode="audit",
                               background_hazard=0.01)
    assert {k: [s.candidate_id for s in v] for k, v in off.replacements.items()} \
        == {k: [s.candidate_id for s in v] for k, v in aud.replacements.items()}
    assert aud.tss_verdicts, "audit must still record verdicts"


def test_require_reverts_model_when_every_alternative_is_refuted():
    """If all alternative starts are unsupported the model stays unsplit."""
    # Strong enough to pass EndpointRefine's own guards (so a split exists to
    # revert), but refutable once tested against an overwhelming background.
    ends, ids = _ends([
        ("a", 30, (1000, 2000)),
        ("b", 30, (1200, 2000)),   # internal start, SAME TES -> tss_only
    ])
    qr = _qr("x", [(1000, 1300), (1500, 2000)], reads=ids)
    off = plan_endpoint_splits({"x": qr}, ends)
    assert "x" in off.replacements, "precondition: the model splits when off"
    plan = plan_endpoint_splits({"x": qr}, ends, tss_evidence_mode="require",
                                background_hazard=0.5)  # huge background
    verdicts = [r["verdict"] for r in _verdict_rows(plan)]
    assert verdicts, "require must still record verdicts"
    # Not vacuous: EVERY alternative must actually be refuted here...
    assert all(v == "unsupported" for v in verdicts), verdicts
    # ...and the model must therefore have REVERTED to unsplit.
    assert "x" not in plan.replacements
    assert "x" not in plan.read_routes
    assert "x" not in plan.primary
