"""P2 containment decision layer (containment_evidence_mode).

Contract under test:
1. Pair discovery/routing: contiguous-subchain shorts route to own_tes
   (distinct 3' end) or tss_only (same TES, internal TSS); mono shorts nest
   by span; non-contained pairs are ignored.
2. A short may be contradicted ONLY when the parent's unique junctions carry
   >= 2 indel-clean exact complete-chain reads; otherwise the pair abstains.
3. Tri-state actions: require mode drops contradicted shorts and exempts
   supported ones; audit mode never acts; unidentifiable never acts.
4. Shard merge is deterministic and attaches grouped q-values to
   hazard-tested rows only.
5. Config validation: mode choices; m2_em requirement; rescue >= 0.
"""

from __future__ import annotations

import pytest

from fin.analysis.containment_evidence import (
    ContainmentRow,
    _classify_pair,
    _parent_unique_min_depth,
    containment_interval_actions,
    evaluate_containment,
    finalize_containment,
    write_containment_tsv,
)
from fin.analysis.quantification import QuantResult
from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.pipeline.config import PipelineConfig


def _cand(cid, chain, start, end, source="novel", strand="+", depths=None):
    c = TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=tuple(chain)),
        three_prime_pos=end if strand == "+" else start,
        sequence="A",
        source=source,
        supporting_read_ids=set(),
        chrom="chr1",
        strand=strand,
        start=start,
        end=end,
    )
    if depths is not None:
        c.complete_chain_junction_depths = tuple(depths)
        c.complete_chain_weakest_indel_clean_depth = min(depths) if depths else 0
    return c


def _qr(cid, assigned):
    return QuantResult(
        candidate_id=cid, abundance=float(len(assigned)), confidence=1.0,
        num_assigned_reads=len(assigned), source="novel",
        assigned_read_ids=tuple(assigned),
    )


# L: exons (1000-1200)(1500-1700)(2000-2200); chain j1=(1200,1500) j2=(1700,2000)
L_CHAIN = ((1200, 1500), (1700, 2000))


class TestClassifyPair:
    def test_suffix_same_tes_is_tss_only(self):
        L = _cand("L", L_CHAIN, 1000, 2200)
        S = _cand("S", (L_CHAIN[1],), 1550, 2200)
        assert _classify_pair(S, L) == "tss_only"

    def test_distinct_tes_is_own_tes(self):
        L = _cand("L", L_CHAIN, 1000, 2200)
        S = _cand("S", (L_CHAIN[0],), 1000, 1690)
        assert _classify_pair(S, L) == "own_tes"

    def test_mono_span_route(self):
        L = _cand("L", L_CHAIN, 1000, 2200)
        S = _cand("S", (), 2000, 2100)
        assert _classify_pair(S, L) == "mono_span"

    def test_non_subchain_not_contained(self):
        L = _cand("L", L_CHAIN, 1000, 2200)
        S = _cand("S", ((1200, 1400),), 1000, 2200)
        assert _classify_pair(S, L) is None

    def test_equal_chain_same_endpoints_not_contained(self):
        L = _cand("L", L_CHAIN, 1000, 2200)
        S = _cand("S", L_CHAIN, 1000, 2200)
        assert _classify_pair(S, L) is None

    def test_equal_chain_interior_start_is_contained(self):
        """The central internal-TSS case: identical chain, same TES, 5'
        start strictly inside the parent."""
        L = _cand("L", L_CHAIN, 1000, 2200)
        S = _cand("S", L_CHAIN, 1150, 2200)
        assert _classify_pair(S, L) == "equal_chain"


class TestParentEvidence:
    def test_unique_junction_depth(self):
        L = _cand("L", L_CHAIN, 1000, 2200, depths=(5, 3))
        S = _cand("S", (L_CHAIN[1],), 1550, 2200)
        assert _parent_unique_min_depth(S, L) == 5  # j1 is unique to L

    def test_no_collected_evidence(self):
        L = _cand("L", L_CHAIN, 1000, 2200)
        S = _cand("S", (L_CHAIN[1],), 1550, 2200)
        assert _parent_unique_min_depth(S, L) == -1


def _run_pair(read_specs, s_start=1550, s_end=2200, depths=(5, 5)):
    """read_specs: list of (read_id, span_start, span_end, owner)."""
    L = _cand("L", L_CHAIN, 1000, 2200, depths=depths)
    S = _cand("S", (L_CHAIN[1],), s_start, s_end)
    spans = {r[0]: (r[1], r[2]) for r in read_specs}
    ql = _qr("L", [r[0] for r in read_specs if r[3] == "L"])
    qs = _qr("S", [r[0] for r in read_specs if r[3] == "S"])
    rows = evaluate_containment(
        [L, S], [ql, qs], spans, interval_name="chr1:1000-2200")
    return rows


class TestEvaluate:
    def test_parent_unsupported_abstains(self):
        reads = [(f"r{i}", 1000, 2200, "L") for i in range(20)]
        rows = _run_pair(reads, depths=(1, 1))
        assert len(rows) == 1
        assert rows[0].verdict == "unidentifiable"
        assert rows[0].reason.startswith("parent_unsupported")

    def test_degradation_only_contradicts(self):
        # 5' ends spread smoothly along the parent -> no spike at S's start
        reads = [("r%d" % i, 1000 + i * 55, 2200, "L") for i in range(20)]
        rows = _run_pair(reads)
        assert len(rows) == 1
        assert rows[0].rung == "tss_only"
        # smooth decay must never yield 'supported'
        assert rows[0].verdict in ("contradicted", "unidentifiable")

    def test_sharp_spike_supports(self):
        reads = [("u%d" % i, 1000 + (i % 3), 2200, "L") for i in range(12)]
        reads += [("s%d" % i, 1550 + (i % 2), 2200, "S") for i in range(10)]
        rows = _run_pair(reads)
        assert len(rows) == 1
        assert rows[0].verdict == "supported", rows[0].reason

    def test_distinct_tes_cluster_supports(self):
        # S ends at 1690 (own TES); 6 reads end there, 12 read through
        reads = [("u%d" % i, 1000, 2200, "L") for i in range(12)]
        reads += [("s%d" % i, 1000, 1690 + (i % 2), "S") for i in range(6)]
        L = _cand("L", L_CHAIN, 1000, 2200, depths=(5, 5))
        S = _cand("S", (L_CHAIN[0],), 1000, 1690)
        spans = {r[0]: (r[1], r[2]) for r in reads}
        ql = _qr("L", [r[0] for r in reads if r[3] == "L"])
        qs = _qr("S", [r[0] for r in reads if r[3] == "S"])
        rows = evaluate_containment(
            [L, S], [ql, qs], spans, interval_name="i")
        assert rows and rows[0].rung == "own_tes"
        assert rows[0].verdict == "supported", rows[0].reason

    def test_absent_tes_cluster_contradicts(self):
        reads = [("u%d" % i, 1000, 2195 + (i % 5), "L") for i in range(18)]
        L = _cand("L", L_CHAIN, 1000, 2200, depths=(5, 5))
        S = _cand("S", (L_CHAIN[0],), 1000, 1690)
        spans = {r[0]: (r[1], r[2]) for r in reads}
        rows = evaluate_containment(
            [L, S], [_qr("L", [r[0] for r in reads]), _qr("S", [])],
            spans, interval_name="i")
        assert rows and rows[0].rung == "own_tes"
        assert rows[0].verdict == "contradicted", rows[0].reason


class TestIntervalActions:
    def _rows(self):
        r1 = ContainmentRow(
            interval="i", short_id="a", parent_id="L", short_source="novel",
            parent_source="novel", rung="tss_only", verdict="contradicted",
            reason="x", test_family="hazard")
        r2 = ContainmentRow(
            interval="i", short_id="b", parent_id="L", short_source="gtf",
            parent_source="novel", rung="own_tes", verdict="supported",
            reason="y", test_family="structural_tes")
        r3 = ContainmentRow(
            interval="i", short_id="c", parent_id="L", short_source="novel",
            parent_source="novel", rung="tss_only",
            verdict="unidentifiable", reason="z", test_family="hazard")
        r4 = ContainmentRow(
            interval="i", short_id="d", parent_id="L", short_source="gtf",
            parent_source="novel", rung="own_tes", verdict="contradicted",
            reason="w", test_family="structural_tes")
        r5 = ContainmentRow(
            interval="i", short_id="e", parent_id="L", short_source="novel",
            parent_source="novel", rung="equal_chain", verdict="supported",
            reason="v", test_family="hazard")
        return [r1, r2, r3, r4, r5]

    def test_audit_never_acts(self):
        drops, exempts = containment_interval_actions(self._rows(), "audit")
        assert drops == set() and exempts == set()

    def test_require_structural_acts_hazard_defers(self):
        rows = self._rows()
        drops, exempts = containment_interval_actions(rows, "require")
        # structural TES contradiction drops now; hazard contradiction
        # NEVER drops at interval level (needs the global FDR pass)
        assert drops == {"d"}
        # structural supported + provisionally supported hazard both exempt
        assert exempts == {"b", "e"}
        assert rows[0].action == "none"      # hazard contradicted: deferred
        assert rows[3].action == "drop"
        assert rows[2].action == "none"


def _qr_simple(cid, verdict=""):
    q = QuantResult(
        candidate_id=cid, abundance=1.0, confidence=1.0,
        num_assigned_reads=1, source="novel")
    q.containment_verdict = verdict
    return q


class TestFinalize:
    def _mk_shards(self, work, rows):
        shard_dir = work / "_containment_shards"
        shard_dir.mkdir()
        for i, r in enumerate(rows):
            write_containment_tsv(shard_dir / f"{i}.tsv", [r])

    def test_merge_deterministic_q_attached_and_hazard_actions(self, tmp_path):
        r1 = ContainmentRow(
            interval="i2", short_id="s2", parent_id="p2",
            short_source="novel", parent_source="novel", rung="tss_only",
            verdict="supported", reason="ok", test_family="hazard",
            p_value=0.0001)
        r2 = ContainmentRow(
            interval="i1", short_id="s1", parent_id="p1",
            short_source="gtf", parent_source="novel", rung="own_tes",
            verdict="contradicted", reason="tes",
            test_family="structural_tes")
        self._mk_shards(tmp_path, [r1, r2])
        agg = {"s1": _qr_simple("s1"), "s2": _qr_simple("s2")}
        agg, dropped, out = finalize_containment(tmp_path, agg, "require")
        lines = (tmp_path / "containment_evidence.tsv").read_text().splitlines()
        assert len(lines) == 3
        assert lines[1].startswith("i1\ts1")      # globally sorted
        assert lines[2].startswith("i2\ts2")
        assert not (tmp_path / "_containment_shards").exists()
        assert lines[2].split("\t")[-1] != ""     # hazard q attached
        assert lines[1].split("\t")[-1] == ""     # structural has no q
        assert dropped == []                       # supported hazard kept
        assert agg["s2"].containment_verdict == "supported"

    def test_fdr_demotion_changes_require_action(self, tmp_path):
        # one strong parent + one weak parent: BY across 2 parents demotes
        # the weak supported row (q > alpha) -> its short is dropped
        rows = [
            ContainmentRow(
                interval="i", short_id="weak", parent_id="pw",
                short_source="novel", parent_source="novel",
                rung="tss_only", verdict="supported", reason="ok",
                test_family="hazard", p_value=0.04),  # BY q=0.06 > 0.05
            ContainmentRow(
                interval="i", short_id="strong", parent_id="ps",
                short_source="novel", parent_source="novel",
                rung="tss_only", verdict="supported", reason="ok",
                test_family="hazard", p_value=1e-8),
        ]
        self._mk_shards(tmp_path, rows)
        agg = {"weak": _qr_simple("weak"), "strong": _qr_simple("strong")}
        agg, dropped, _ = finalize_containment(tmp_path, agg, "require")
        assert dropped == ["weak"]                # demoted by grouped BY
        assert "weak" not in agg
        assert agg["strong"].containment_verdict == "supported"

    def test_audit_mode_never_mutates(self, tmp_path):
        rows = [ContainmentRow(
            interval="i", short_id="x", parent_id="p",
            short_source="novel", parent_source="novel", rung="tss_only",
            verdict="contradicted", reason="r", test_family="hazard",
            p_value=0.9)]
        self._mk_shards(tmp_path, rows)
        agg = {"x": _qr_simple("x")}
        agg, dropped, _ = finalize_containment(tmp_path, agg, "audit")
        assert dropped == [] and "x" in agg
        assert agg["x"].containment_verdict == ""  # audit sets nothing

    def test_header_only_artifact_when_no_pairs(self, tmp_path):
        agg = {}
        agg, dropped, out = finalize_containment(tmp_path, agg, "audit")
        lines = (tmp_path / "containment_evidence.tsv").read_text().splitlines()
        assert len(lines) == 1 and lines[0].startswith("interval\t")


class TestEqualChainEvaluation:
    def _run(self, read_specs, s_start=1600):
        L = _cand("L", L_CHAIN, 1000, 2200, depths=(5, 5))
        S = _cand("S", L_CHAIN, s_start, 2200)
        spans = {r[0]: (r[1], r[2]) for r in read_specs}
        ql = _qr("L", [r[0] for r in read_specs if r[3] == "L"])
        qs = _qr("S", [r[0] for r in read_specs if r[3] == "S"])
        return evaluate_containment(
            [L, S], [ql, qs], spans, interval_name="i")

    def test_spike_supports_equal_chain_short(self):
        reads = [("u%d" % i, 1000 + (i % 3), 2200, "L") for i in range(12)]
        reads += [("s%d" % i, 1600 + (i % 2), 2200, "S") for i in range(10)]
        rows = self._run(reads)
        assert rows and rows[0].rung == "equal_chain"
        assert rows[0].test_family == "hazard"
        assert rows[0].verdict == "supported", rows[0].reason

    def test_smooth_decay_with_upstream_contradicts(self):
        reads = [("r%d" % i, 1000 + i * 55, 2200, "L") for i in range(20)]
        rows = self._run(reads)
        assert rows and rows[0].rung == "equal_chain"
        assert rows[0].verdict in ("contradicted", "unidentifiable")

    def test_no_upstream_cannot_contradict(self):
        # every read starts at/after the short start: the parent's 5' region
        # is unobserved, so contradiction is unfounded -> abstain (or the
        # dominant-start shortcut supports the short model; never
        # contradicted).
        reads = [("r%d" % i, 1650 + i * 9, 2200, "L") for i in range(20)]
        rows = self._run(reads)
        assert rows and rows[0].verdict != "contradicted"

    def test_novel_start_pays_search_correction(self):
        """Same data: a GTF-provided start uses one eligible bin, a novel
        (discovered) start is Bonferroni-corrected over the parent frame ->
        its p-value can only be larger."""
        reads = [("u%d" % i, 1000 + (i % 3), 2200, "L") for i in range(12)]
        reads += [("s%d" % i, 1600 + (i % 2), 2200, "S") for i in range(8)]
        L = _cand("L", L_CHAIN, 1000, 2200, depths=(5, 5))
        spans = {r[0]: (r[1], r[2]) for r in reads}
        ql = _qr("L", [r[0] for r in reads if r[3] == "L"])
        out = {}
        for src in ("gtf", "novel"):
            S = _cand("S", L_CHAIN, 1600, 2200, source=src)
            qs = _qr("S", [r[0] for r in reads if r[3] == "S"])
            rows = evaluate_containment(
                [L, S], [ql, qs], spans, interval_name="i")
            assert rows and rows[0].p_value is not None
            out[src] = rows[0].p_value
        assert out["novel"] >= out["gtf"]


class TestAggregationPreservesVerdict:
    def test_verdict_survives_cross_interval_merge(self):
        from fin.analysis.quantification import aggregate_across_intervals

        q = _qr_simple("c1", verdict="supported")
        merged = aggregate_across_intervals([[q]])
        assert merged["c1"].containment_verdict == "supported"


class TestConfigValidation:
    def test_bad_mode_rejected(self, tmp_path):
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        cfg = PipelineConfig(bam_path=str(bam),
                             containment_evidence_mode="bogus")
        with pytest.raises(ValueError, match="containment_evidence_mode"):
            cfg.validate()

    def test_requires_m2_em(self, tmp_path):
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        cfg = PipelineConfig(bam_path=str(bam), quant_mode="argmax",
                             containment_evidence_mode="require")
        with pytest.raises(ValueError, match="requires --quant-mode m2_em"):
            cfg.validate()

    def test_negative_rescue_rejected(self, tmp_path):
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        cfg = PipelineConfig(bam_path=str(bam),
                             complete_chain_rescue_reads=-1)
        with pytest.raises(ValueError, match="complete_chain_rescue_reads"):
            cfg.validate()

    def test_cli_default_off(self):
        from fin.cli import main

        opt = next(
            p for p in main.params if p.name == "containment_evidence_mode"
        )
        assert opt.default == "off"


class TestMultiParentVerdictPrecedence:
    """A short model is tested once per PARENT, so it owns several rows.

    Resolving them by row order (a dict comprehension keeping whichever row
    sorted last) silently discarded `supported` exemptions -- the P2 finding
    from the 2026-09-02 Codex review. Resolution is by explicit precedence:
    supported > contradicted > unsupported > unidentifiable.
    """

    def _mk(self, work, rows):
        shard_dir = work / "_containment_shards"
        shard_dir.mkdir(exist_ok=True)
        for i, r in enumerate(rows):
            write_containment_tsv(shard_dir / f"{i}.tsv", [r])

    def _row(self, short_id, parent_id, verdict, p=0.0001):
        return ContainmentRow(
            interval="i", short_id=short_id, parent_id=parent_id,
            short_source="novel", parent_source="novel", rung="tss_only",
            verdict=verdict, reason="r", test_family="hazard", p_value=p)

    def test_supported_survives_a_later_sorting_row(self, tmp_path):
        # "pa" sorts before "pz": the last row is the unsupported one, which the
        # old last-row-wins merge would have kept.
        self._mk(tmp_path, [self._row("s", "pa", "supported"),
                            self._row("s", "pz", "unsupported", p=0.9)])
        agg = {"s": _qr_simple("s")}
        agg, dropped, _ = finalize_containment(tmp_path, agg, "require")
        assert dropped == []
        assert agg["s"].containment_verdict == "supported"

    def test_supported_outranks_contradicted_and_is_not_dropped(self, tmp_path):
        # positive evidence against one parent is not erased by another
        # comparison failing against a different parent.
        self._mk(tmp_path, [self._row("s", "p1", "contradicted", p=0.9),
                            self._row("s", "p2", "supported")])
        agg = {"s": _qr_simple("s")}
        agg, dropped, _ = finalize_containment(tmp_path, agg, "require")
        assert dropped == [] and "s" in agg
        assert agg["s"].containment_verdict == "supported"

    def test_contradicted_still_drops_when_nothing_supports(self, tmp_path):
        self._mk(tmp_path, [self._row("s", "p1", "contradicted", p=0.9),
                            self._row("s", "p2", "unidentifiable", p=0.9)])
        agg = {"s": _qr_simple("s")}
        agg, dropped, _ = finalize_containment(tmp_path, agg, "require")
        assert dropped == ["s"] and "s" not in agg

    def test_resolution_is_independent_of_shard_order(self, tmp_path):
        rows = [self._row("s", "p1", "contradicted", p=0.9),
                self._row("s", "p2", "supported")]
        self._mk(tmp_path, rows)
        agg1, d1, _ = finalize_containment(tmp_path, {"s": _qr_simple("s")}, "require")
        self._mk(tmp_path, list(reversed(rows)))
        agg2, d2, _ = finalize_containment(tmp_path, {"s": _qr_simple("s")}, "require")
        assert d1 == d2
        assert agg1["s"].containment_verdict == agg2["s"].containment_verdict
