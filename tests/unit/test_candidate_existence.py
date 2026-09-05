"""Unit tests for the audit-only candidate existence test.

The properties pinned here are the ones the audit must prove BEFORE anything is
allowed to consume the verdict (CANDIDATE_EXISTENCE_DESIGN.md section 6):
determinism/order-invariance, selection honesty, evidence-based components (not
posterior-based), correct held-out accounting, and — the red line — that a
candidate is never dropped and "not significant" is never read as absence.
"""
from __future__ import annotations

import numpy as np
import pytest

from fin.analysis.candidate_existence import (
    CandidateExistenceRow,
    component_folds,
    STATUS_SUPPORTED,
    STATUS_UNRESOLVED,
    apply_sample_fdr,
    build_conflict_components,
    evaluate_candidate_existence,
    merge_existence_shards,
    read_fold,
    write_existence_shard,
)
from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval

MISSING = 1e6


def _cand(cid, source="gtf", introns=((100, 200), (300, 400))):
    return TranscriptCandidate(
        candidate_id=cid, intron_chain=IntronChain(introns=introns),
        three_prime_pos=500, sequence="ACGT" * 30, source=source,
        supporting_read_ids=set(), chrom="chr1", strand="+", start=0, end=500,
        family_id="fam_1")


class _QO:
    """Minimal stand-in for QuantOutput carrying only what the test consumes."""

    def __init__(self, cand_list, kept_read_ids, d_tx, ties, nlls, cover):
        self.cand_list = cand_list
        self.kept_read_ids = kept_read_ids
        self.structural_d_tx_raw = d_tx
        self.d_tx = d_tx
        self.ties_by_read = ties
        self.nlls_by_read = nlls
        self.cover_by_read = cover
        self.R = np.zeros((len(kept_read_ids), len(cand_list)))


def _interval():
    return GenomicInterval(chrom="chr1", start=0, end=1000, strand="+")


def _scenario(n_reads=40, favour=0.0, source_a="gtf", source_b="gtf"):
    """Two competing multi-exon candidates. ``favour`` is the per-read distance
    advantage given to candidate 0 (0.0 = the two are indistinguishable)."""
    cands = [_cand("ENST_A", source_a, ((100, 200), (300, 400))),
             _cand("ENST_B", source_b, ((100, 200), (320, 400)))]
    reads = [f"r{i}" for i in range(n_reads)]
    d = np.full((n_reads, 2), MISSING)
    for i in range(n_reads):
        d[i, 0] = 0.0
        d[i, 1] = favour
    ties = {i: [0, 1] for i in range(n_reads)}
    nlls = {i: {0: 0.0, 1: favour} for i in range(n_reads)}
    cover = {i: True for i in range(n_reads)}
    return _QO(cands, reads, d, ties, nlls, cover)


class TestFolds:
    def test_component_folds_are_balanced_and_deterministic(self):
        ids = [f"r{i}" for i in range(11)]
        f, k = component_folds("cmp_x", ids)
        assert k == 5
        counts = [list(f.values()).count(x) for x in range(k)]
        assert max(counts) - min(counts) <= 1        # balanced
        assert set(f) == set(ids)                    # every read placed once
        f2, k2 = component_folds("cmp_x", list(reversed(ids)))
        assert (f, k) == (f2, k2)                    # order-invariant
        # a different component gets a different (still deterministic) split
        f3, _ = component_folds("cmp_y", ids)
        assert f3 == component_folds("cmp_y", ids)[0]

    def test_component_folds_degrade_to_loo_then_abstain(self):
        for n in (2, 3, 4, 5):
            f, k = component_folds("c", [f"r{i}" for i in range(n)])
            assert k == n and sorted(f.values()) == list(range(n))
        assert component_folds("c", ["only"])[1] == 0
        assert component_folds("c", [])[1] == 0

    def test_fold_is_deterministic_and_spread(self):
        ids = [f"read{i}" for i in range(500)]
        a = [read_fold(r) for r in ids]
        b = [read_fold(r) for r in reversed(ids)][::-1]
        assert a == b
        assert set(a) == {0, 1, 2, 3, 4}          # all folds used
        counts = [a.count(k) for k in range(5)]
        assert min(counts) > 50                    # roughly balanced


class TestComponents:
    def test_edges_need_a_covered_scored_contrast(self):
        qo = _scenario()
        comps = build_conflict_components(qo, None)
        assert len(comps) == 1 and len(comps[0].columns) == 2

    def test_uncovered_reads_still_build_the_component(self):
        # coverage is a SCORING condition, not an edge condition: two candidates a
        # read supports compete for its mass whether or not it can tell them apart.
        qo = _scenario()
        qo.cover_by_read = {i: False for i in qo.cover_by_read}
        comps = build_conflict_components(qo, None)
        assert len(comps) == 1 and len(comps[0].columns) == 2
        # ...but such a read contributes no evidence
        rows = {r.candidate_id: r for r in evaluate_candidate_existence(
            qo, None, comps, _interval(), sample_id="s1")}
        assert rows["ENST_A"].identifiability_reason == "no_held_out_discriminating_read"
        assert rows["ENST_A"].n_discriminating_reads == 0

    def test_component_is_not_built_from_posteriors(self):
        # a dense R must not create a component when the evidence maps are empty
        qo = _scenario()
        qo.ties_by_read, qo.nlls_by_read = {}, {}
        qo.R = np.full((len(qo.kept_read_ids), 2), 0.5)
        assert build_conflict_components(qo, None) == []

    def test_component_id_is_order_invariant(self):
        qo = _scenario()
        a = build_conflict_components(qo, None)[0]
        qo2 = _scenario()
        qo2.cand_list = list(reversed(qo2.cand_list))
        qo2.d_tx = qo2.d_tx[:, ::-1]
        qo2.ties_by_read = {i: [1, 0] for i in qo2.ties_by_read}
        qo2.nlls_by_read = {i: {1: 0.0, 0: 0.0} for i in qo2.nlls_by_read}
        b = build_conflict_components(qo2, None)[0]
        assert a.component_id == b.component_id


class TestEvaluation:
    def _rows(self, qo):
        comps = build_conflict_components(qo, None)
        return {r.candidate_id: r for r in evaluate_candidate_existence(
            qo, None, comps, _interval(), sample_id="s1")}

    def test_indistinguishable_candidates_contribute_no_evidence(self):
        # identical likelihood columns: adding either explains nothing new.
        rows = self._rows(_scenario(favour=0.0))
        for r in rows.values():
            assert r.identifiable
            assert abs(r.delta_predictive_ll) < 1e-6

    def test_a_distinguishing_candidate_earns_positive_evidence(self):
        # candidate 0 fits every read strictly better -> removing it must hurt.
        rows = self._rows(_scenario(favour=3.0))
        assert rows["ENST_A"].delta_predictive_ll > 0
        assert rows["ENST_A"].n_discriminating_reads > 0
        assert rows["ENST_A"].n_tested_folds == 5

    def test_novel_candidates_are_not_selection_honest_in_v1(self):
        rows = self._rows(_scenario(favour=3.0, source_a="novel"))
        a = rows["ENST_A"]
        assert a.selection_honest is False and a.identifiable is False
        assert a.identifiability_reason == "not_selection_honest_v1"
        assert a.delta_predictive_ll is None

    def test_mono_exon_annotation_is_not_tested(self):
        qo = _scenario(favour=3.0)
        qo.cand_list[0] = _cand("ENST_A", "gtf", introns=())
        rows = self._rows(qo)
        assert rows["ENST_A"].identifiability_reason == "not_multi_exon"

    def test_every_read_is_held_out_exactly_once(self):
        qo = _scenario(favour=3.0, n_reads=60)
        rows = self._rows(qo)
        # each read contributes to exactly one fold's delta, and all reads here
        # discriminate, so the count equals the read total.
        assert rows["ENST_A"].n_discriminating_reads == 60
        assert rows["ENST_A"].n_candidate_folds == 5

    def test_small_components_degrade_to_leave_one_out(self):
        # 3 reads: a global hash-mod-5 would leave most folds empty; the adaptive
        # scheme uses K=3 so every read is still scored exactly once.
        qo = _scenario(favour=3.0, n_reads=3)
        rows = self._rows(qo)
        r = rows["ENST_A"]
        assert r.n_candidate_folds == 3 and r.n_tested_folds == 3
        assert r.n_discriminating_reads == 3 and r.identifiable

    def test_single_read_component_has_no_honest_split(self):
        qo = _scenario(favour=3.0, n_reads=1)
        rows = self._rows(qo)
        assert rows["ENST_A"].identifiability_reason == "no_honest_split"
        assert rows["ENST_A"].identifiable is False

    def test_gate_style_one_hot_evidence_is_read_from_the_raw_kernel(self):
        # Production writes a DECISION matrix into d_tx: a distinguishing read
        # becomes one-hot. Scoring that would discard the strongest evidence, so
        # the module must read structural_d_tx_raw instead.
        qo = _scenario(favour=3.0, n_reads=20)
        qo.d_tx = np.full_like(qo.structural_d_tx_raw, MISSING)
        qo.d_tx[:, 0] = 0.0                       # one-hot on the winner
        rows = self._rows(qo)
        assert rows["ENST_A"].n_discriminating_reads == 20
        assert rows["ENST_A"].delta_predictive_ll > 0

    def test_no_bootstrap_means_no_p_value(self):
        # without a null nothing may claim a p-value; with one, it is valid and
        # the status is still decided only at sample-level finalize.
        qo = _scenario(favour=3.0)
        comps = build_conflict_components(qo, None)
        rows = {r.candidate_id: r for r in evaluate_candidate_existence(
            qo, None, comps, _interval(), sample_id="s1", bootstrap_replicates=0)}
        for r in rows.values():
            assert r.p_value is None and r.p_valid is False
            assert r.status == STATUS_UNRESOLVED and r.decision == "no_action"
        rows = {r.candidate_id: r for r in evaluate_candidate_existence(
            qo, None, comps, _interval(), sample_id="s1", bootstrap_replicates=99)}
        a = rows["ENST_A"]
        assert a.p_valid and a.p_value is not None and a.p_minus is not None
        assert a.null_method == "label_permutation_molecule_bootstrap"
        assert a.status == STATUS_UNRESOLVED and a.decision == "no_action"

    def test_tss_and_polya_channels_are_empty_not_zero(self):
        rows = self._rows(_scenario(favour=3.0))
        for r in rows.values():
            assert r.tss_llr is None and r.tes_polya_llr is None

    def test_row_order_does_not_change_statistics(self):
        qo = _scenario(favour=3.0, n_reads=50)
        a = self._rows(qo)["ENST_A"].delta_predictive_ll
        qo2 = _scenario(favour=3.0, n_reads=50)
        order = list(reversed(range(len(qo2.kept_read_ids))))
        qo2.kept_read_ids = [qo2.kept_read_ids[i] for i in order]
        qo2.d_tx = qo2.d_tx[order]
        qo2.ties_by_read = {n: qo2.ties_by_read[i] for n, i in enumerate(order)}
        qo2.nlls_by_read = {n: qo2.nlls_by_read[i] for n, i in enumerate(order)}
        qo2.cover_by_read = {n: True for n in range(len(order))}
        b = self._rows(qo2)["ENST_A"].delta_predictive_ll
        assert a == pytest.approx(b, rel=1e-9)


class TestFdr:
    """Directional hierarchical control: component BY, within-component Holm."""

    def _row(self, cid, p, comp="c1", *, p_minus=0.9, cp=None, cm=None,
             identifiable=True, honest=True, delta=1.0, neg=""):
        return CandidateExistenceRow(
            candidate_id=cid, component_id=comp, p_value=p, p_minus=p_minus,
            component_p_plus=(p if cp is None else cp),
            component_p_minus=(p_minus if cm is None else cm),
            p_valid=p is not None, identifiable=identifiable,
            selection_honest=honest, delta_predictive_ll=delta, negative_channel=neg)

    def test_family_excludes_invalid_rows(self):
        rows = [self._row("a", 0.001), self._row("b", None, comp="c2"),
                self._row("c", 0.001, comp="c3", identifiable=False),
                self._row("d", 0.001, comp="c4", honest=False)]
        out = apply_sample_fdr(rows, alpha=0.05)
        assert all(r.fdr_family_size == 1 for r in out)     # one component corrected
        assert out[1].status == STATUS_UNRESOLVED and out[2].status == STATUS_UNRESOLVED

    def test_component_by_is_monotone_and_conservative(self):
        rows = [self._row(f"c{i}", p, comp=f"comp{i}") for i, p in enumerate(
            [0.0001, 0.002, 0.02, 0.2, 0.6])]
        out = apply_sample_fdr(rows, alpha=0.05)
        qs = [r.component_q_plus for r in sorted(out, key=lambda r: r.p_value)]
        assert all(b >= a - 1e-12 for a, b in zip(qs, qs[1:]))
        assert all(q >= p for q, p in zip(qs, sorted(r.p_value for r in rows)))
        assert out[0].status == STATUS_SUPPORTED and out[-1].status == STATUS_UNRESOLVED

    def test_holm_within_component(self):
        # two candidates in ONE component: omnibus passes, Holm decides each
        rows = [self._row("a", 0.001, cp=0.001), self._row("b", 0.03, cp=0.001)]
        out = {r.candidate_id: r for r in apply_sample_fdr(rows, alpha=0.05)}
        assert out["a"].status == STATUS_SUPPORTED
        assert out["a"].holm_p_plus == pytest.approx(0.002)          # 2 * 0.001
        assert out["b"].holm_p_plus == pytest.approx(0.03) and out["b"].status == STATUS_UNRESOLVED

    def test_nonsignificant_is_unresolved_never_dropped(self):
        out = apply_sample_fdr([self._row("a", 0.9)], alpha=0.05)
        assert out[0].status == STATUS_UNRESOLVED and out[0].decision == "no_action"

    def test_supported_needs_positive_evidence(self):
        rows = [self._row("a", 1e-4, delta=-5.0), self._row("b", 1e-4, comp="c2", delta=5.0)]
        out = {r.candidate_id: r for r in apply_sample_fdr(rows, alpha=0.05)}
        assert out["a"].status == STATUS_UNRESOLVED and out["b"].status == STATUS_SUPPORTED

    def test_contradicted_needs_negative_channel(self):
        # a negative effect with a tiny p- but NO physical negative channel stays unresolved
        rows = [self._row("a", 0.9, p_minus=1e-4, cm=1e-4, delta=-6.0),
                self._row("b", 0.9, comp="c2", p_minus=1e-4, cm=1e-4, delta=-6.0,
                          neg="tes_readthrough_polya")]
        out = {r.candidate_id: r for r in apply_sample_fdr(rows, alpha=0.05)}
        assert out["a"].status == STATUS_UNRESOLVED
        assert out["b"].status == "contradicted"

    def test_directional_alpha_is_half(self):
        out = apply_sample_fdr([self._row("a", 0.02)], alpha=0.05)
        assert out[0].directional_alpha == pytest.approx(0.025)
        assert out[0].status == STATUS_SUPPORTED          # 0.02 <= 0.025 (single component)


class TestArtifact:
    def test_shard_roundtrip_applies_sample_fdr(self, tmp_path):
        rows = [CandidateExistenceRow(candidate_id="a", component_id="c1", p_value=1e-6,
                                      p_minus=0.9, component_p_plus=1e-6, component_p_minus=0.9,
                                      p_valid=True, identifiable=True, selection_honest=True,
                                      delta_predictive_ll=4.0, chrom="chr1"),
                CandidateExistenceRow(candidate_id="b", component_id="c2", p_value=0.5,
                                      p_minus=0.5, component_p_plus=0.5, component_p_minus=0.5,
                                      p_valid=True, identifiable=True, selection_honest=True,
                                      delta_predictive_ll=0.1, chrom="chr1")]
        write_existence_shard(tmp_path, _interval(), rows)
        out = merge_existence_shards(tmp_path)
        assert out is not None
        lines = [l.split("\t") for l in open(out).read().splitlines()]
        header, body = lines[0], lines[1:]
        idx = {k: i for i, k in enumerate(header)}
        by = {r[idx["candidate_id"]]: r for r in body}
        assert by["a"][idx["status"]] == STATUS_SUPPORTED
        assert by["b"][idx["status"]] == STATUS_UNRESOLVED
        assert by["a"][idx["fdr_family_size"]] == "2"
        assert all(r[idx["decision"]] == "no_action" for r in body)
        assert not (tmp_path / "_existence_shards").exists()

    def test_no_rows_writes_nothing(self, tmp_path):
        write_existence_shard(tmp_path, _interval(), [])
        assert merge_existence_shards(tmp_path) is None
