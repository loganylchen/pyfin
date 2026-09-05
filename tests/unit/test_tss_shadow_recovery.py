"""Wiring tests for post-EM folded-TSS shadow recovery on the production
``_quant_m2_em`` seam (fin/analysis/tss_shadow_recovery.py).

Real mappy + real EM on a synthetic genome; signal-free (``m2_metric="off"``
so no krill is needed). A parent 3-exon chain carries one folded shadow, the
exact 5'-suffix chain. In the SHARP fixture the shadow reads pile their 5'
ends at one position inside exon 2 (a genuine TSS) on top of a degradation
ramp of parent reads; in the BROAD fixture the same number of shadow reads
start all over exon 2 (degradation), which must NOT be recovered.
"""
from __future__ import annotations

import random
from pathlib import Path

import numpy as np
import pytest

from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner

mappy = pytest.importorskip("mappy")

_RNG = random.Random(20260903)
GENOME = "".join(_RNG.choice("ACGT") for _ in range(6000))
EXONS = [(500, 900), (1500, 1900), (2500, 3100)]
CHAIN = ((900, 1500), (1900, 2500))
SHADOW_CHAIN = ((1900, 2500),)          # exact 5' suffix of CHAIN (+ strand)
TSS = 1700                               # sharp short-isoform start inside exon 2


def _spliced(start, end):
    parts = []
    for s, e in EXONS:
        lo, hi = max(s, start), min(e, end)
        if hi > lo:
            parts.append(GENOME[lo:hi])
    return "".join(parts)


def _parent():
    return TranscriptCandidate(
        candidate_id="parent", intron_chain=IntronChain(introns=CHAIN),
        three_prime_pos=3100, sequence=_spliced(500, 3100), source="novel",
        supporting_read_ids=set(), chrom="chr1", strand="+", start=500, end=3100,
        family_id="fam_x",
    )


def _build(shadow_starts):
    """CandidateSet with the parent + folded shadow; read spans/sequences for
    6 full-length parent reads, an 8-read 5' degradation ramp, and the shadow
    reads starting at ``shadow_starts``."""
    spans, seqs = {}, {}
    ids = []
    for i in range(6):
        rid = f"full{i}"; spans[rid] = (500, 3100); seqs[rid] = _spliced(500, 3100); ids.append(rid)
    for k, st in enumerate((560, 620, 690, 760, 830, 1530, 1580, 1620)):
        rid = f"deg{k}"; spans[rid] = (st, 3100); seqs[rid] = _spliced(st, 3100); ids.append(rid)
    shadow = []
    for k, st in enumerate(shadow_starts):
        rid = f"sh{k}"; spans[rid] = (st, 3100); seqs[rid] = _spliced(st, 3100); ids.append(rid); shadow.append(rid)
    parent = _parent()
    parent.supporting_read_ids.update(ids)
    cs = CandidateSet(
        interval=_interval(), candidates=[parent], read_ids=set(ids),
        read_sequences=seqs, clusters=[["parent"]],
        shadows={"parent": [(SHADOW_CHAIN, tuple(sorted(shadow)))]},
        read_spans=spans, indels_by_read={rid: () for rid in ids},
    )
    return cs, ids


def _interval():
    return GenomicInterval(chrom="chr1", start=0, end=6000, strand="+")


def _runner(tmp_path, mode):
    cfg = PipelineConfig(
        bam_path="/tmp/x.bam", quant_mode="m2_em", use_gpu=False, m2_metric="off",
        work_dir=str(tmp_path), folded_tss_recovery_mode=mode,
    )
    r = PipelineRunner(cfg)
    r._genome_fasta = {"chr1": GENOME}
    r._observed_junctions = lambda interval: None      # support gates fail-open
    calls = []
    orig = r._assigner.assign

    def _spy(*a, **k):
        qo = orig(*a, **k)
        calls.append(qo)
        return qo

    r._assigner.assign = _spy
    return r, calls


def _rows(tmp_path):
    from fin.analysis.tss_shadow_recovery import SHARD_DIR, TSV_FIELDS
    rows = []
    for f in sorted((Path(tmp_path) / SHARD_DIR).glob("*.tsv")):
        for line in f.read_text().splitlines():
            rows.append(dict(zip(TSV_FIELDS, line.split("\t"))))
    return rows


SHARP = [TSS] * 14
BROAD = [1520 + 26 * k for k in range(14)]       # 1520..1858, spread over exon 2


def test_sharp_shadow_is_recovered_and_em_rerun_once(tmp_path):
    cs, ids = _build(SHARP)
    runner, calls = _runner(tmp_path, "require")
    results, _ = runner._quant_m2_em(cs, ids, _interval())
    assert len(calls) == 2                               # exactly one rerun
    qo = calls[-1]
    rec = [c for c in cs.candidates if c.candidate_id != "parent"]
    assert len(rec) == 1
    r = rec[0]
    assert tuple(r.intron_chain.introns) == SHADOW_CHAIN
    assert (r.start, r.end) == (TSS, 3100)
    assert r.source == "novel" and r.family_id == "fam_x"
    assert r.supporting_read_ids == {f"sh{k}" for k in range(14)}
    assert r.candidate_id in [c.candidate_id for c in qo.cand_list]
    assert cs.clusters == [["parent", r.candidate_id]]
    # second-pass EM owns the reads: every row still sums to one, and the
    # total mass over parent+recovered equals the kept reads (no manual transfer)
    np.testing.assert_allclose(qo.R.sum(axis=1), 1.0, atol=1e-9)
    assert abs(qo.R.sum() - len(qo.kept_read_ids)) < 1e-9
    # the frozen 5' positional likelihood lets the 14 in-bin shadow reads WIN the
    # recovered model by argmax (hard 14 / 14); the 14 upstream reads are excluded
    # from it (MISSING). Soft mass ~ 14 * 1/(1+h) with h ~ 0.019 -> ~13.7-13.9.
    js = [c.candidate_id for c in qo.cand_list].index(r.candidate_id)
    jp = [c.candidate_id for c in qo.cand_list].index("parent")
    hard = qo.hard_assignments
    assert int((hard == js).sum()) == 14 and int((hard == jp).sum()) == 14
    assert 13.5 < float(qo.R[:, js].sum()) < 14.0
    for i, rid in enumerate(qo.kept_read_ids):
        if rid.startswith("full") or rid.startswith("deg"):
            assert qo.R[i, js] == 0.0                     # upstream: not from S
    # ...and it SURVIVES the ordinary selection cascade (no exemption: it passes
    # max_soft_mass_ratio on its own hard reads).
    out = {q.candidate_id: q for q in results}
    assert set(out) == {"parent", r.candidate_id}
    assert out[r.candidate_id].num_assigned_reads == 14
    rows = _rows(tmp_path)
    assert len(rows) == 1 and rows[0]["verdict"] == "supported"
    assert rows[0]["reemitted_candidate_id"] == r.candidate_id
    assert int(rows[0]["tss_offset"]) == 400 + (TSS - 1500)


def test_audit_records_but_does_not_rerun(tmp_path):
    cs, ids = _build(SHARP)
    runner, calls = _runner(tmp_path, "audit")
    runner._quant_m2_em(cs, ids, _interval())
    assert len(calls) == 1
    assert [c.candidate_id for c in cs.candidates] == ["parent"]
    rows = _rows(tmp_path)
    assert len(rows) == 1 and rows[0]["verdict"] == "supported"


def test_broad_shadow_is_not_recovered(tmp_path):
    cs, ids = _build(BROAD)
    runner, calls = _runner(tmp_path, "require")
    runner._quant_m2_em(cs, ids, _interval())
    assert len(calls) == 1
    assert [c.candidate_id for c in cs.candidates] == ["parent"]
    rows = _rows(tmp_path)
    assert len(rows) == 1 and rows[0]["verdict"] != "supported"


def test_off_mode_is_inert(tmp_path):
    cs, ids = _build(SHARP)
    runner, calls = _runner(tmp_path, "off")
    runner._quant_m2_em(cs, ids, _interval())
    assert len(calls) == 1
    assert not (Path(tmp_path) / "_tss_recovery_shards").exists()


def test_eligibility_rules_abstain():
    from fin.analysis.tss_shadow_recovery import recover_folded_tss_candidates

    class _QO:
        def __init__(self, cand_list, kept, R):
            self.cand_list, self.kept_read_ids, self.R = cand_list, kept, R

    cs, ids = _build(SHARP)
    parent = cs.candidates[0]
    R = np.ones((len(ids), 1))
    qo = _QO([parent], ids, R)
    # not a pure 5' truncation (prefix on + strand) -> abstain
    cs.shadows = {"parent": [(((900, 1500),), tuple(f"sh{k}" for k in range(14)))]}
    out = recover_folded_tss_candidates(cs, qo, _interval(), GENOME)
    assert out.candidates == [] and out.rows[0].reason == "not_pure_5p_truncation"
    # own-TES shadow (3' ends far from the parent's) -> abstain
    cs.shadows = {"parent": [(SHADOW_CHAIN, tuple(f"sh{k}" for k in range(14)))]}
    for k in range(14):
        cs.read_spans[f"sh{k}"] = (TSS, 2900)
    out = recover_folded_tss_candidates(cs, qo, _interval(), GENOME)
    assert out.candidates == [] and out.rows[0].reason.startswith("tes_disagrees")
    # parent with zero first-pass mass -> abstain
    for k in range(14):
        cs.read_spans[f"sh{k}"] = (TSS, 3100)
    out = recover_folded_tss_candidates(cs, _QO([parent], ids, np.zeros((len(ids), 1))),
                                        _interval(), GENOME)
    assert out.candidates == [] and out.rows[0].reason == "parent_zero_abundance"
    # fewer than 3 shadow reads -> abstain
    cs.shadows = {"parent": [(SHADOW_CHAIN, ("sh0", "sh1"))]}
    out = recover_folded_tss_candidates(cs, qo, _interval(), GENOME)
    assert out.candidates == [] and out.rows[0].reason.startswith("too_few_shadow_reads")


def test_first_pass_is_byte_identical_with_recovery_on(tmp_path):
    # the first assignment never sees the positional term: R/hard identical to off.
    cs_off, ids = _build(SHARP)
    r_off, calls_off = _runner(tmp_path / "off", "off")
    r_off._quant_m2_em(cs_off, ids, _interval())
    cs_req, _ = _build(SHARP)
    r_req, calls_req = _runner(tmp_path / "req", "require")
    r_req._quant_m2_em(cs_req, ids, _interval())
    assert np.array_equal(calls_off[0].R, calls_req[0].R)
    assert np.array_equal(calls_off[0].hard_assignments, calls_req[0].hard_assignments)


class TestApplyStartLikelihood:
    """Unit rules of the positional term (fin/pipeline/assignment.py)."""

    def _setup(self):
        from fin.pipeline.assignment import MISSING
        parent = _parent()
        short = TranscriptCandidate(
            candidate_id="short", intron_chain=IntronChain(introns=SHADOW_CHAIN),
            three_prime_pos=3100, sequence=_spliced(TSS, 3100), source="novel",
            supporting_read_ids=set(), chrom="chr1", strand="+", start=TSS, end=3100)
        other = TranscriptCandidate(
            candidate_id="other", intron_chain=IntronChain(introns=()),
            three_prime_pos=3100, sequence="ACGT" * 50, source="novel",
            supporting_read_ids=set(), chrom="chr1", strand="+", start=0, end=200)
        cands = [parent, short, other]
        spans = {"peak": (TSS, 3100), "up": (600, 3100), "down": (1800, 3100),
                 "three": (TSS, 3100), "nospan": (TSS, 3100)}
        model = dict(parent_id="parent", recovered_id="short", tss_bin=(400 + TSS - 1500) // 25,
                     bin_bp=25, hazard=0.019, exons=EXONS, strand="+",
                     read_spans={k: v for k, v in spans.items() if k != "nospan"})
        reads = ["peak", "up", "down", "three", "nospan"]
        ties = {0: [0, 1], 1: [0, 1], 2: [0, 1], 3: [0, 1, 2], 4: [0, 1]}
        d = np.full((5, 3), MISSING)
        for i, t in ties.items():
            for j in t:
                d[i, j] = 0.0
        return d, ties, reads, cands, [model], MISSING

    def test_rules(self):
        from fin.pipeline.assignment import apply_start_likelihood
        d, ties, reads, cands, models, MISSING = self._setup()
        before = d.copy()
        n = apply_start_likelihood(d, ties, reads, cands, models)
        assert n == 2
        assert d[0, 0] == pytest.approx(-np.log(0.019)) and d[0, 1] == 0.0   # peak bin
        assert d[1, 1] == MISSING and d[1, 0] == 0.0                         # upstream
        assert np.array_equal(d[2], before[2])                               # downstream
        assert np.array_equal(d[3], before[3])                               # 3-way tie
        assert np.array_equal(d[4], before[4])                               # no span
        assert np.array_equal(d[:, 2], before[:, 2])                         # other col

    def test_support_gate_evidence_mirrors_the_contrast(self):
        from fin.pipeline.assignment import apply_start_likelihood
        from fin.scoring.m2_junction_nll import support_gate_drops
        d, ties, reads, cands, models, MISSING = self._setup()
        nlls = {2: {0: 1.0, 1: 3.0}}                     # downstream read: untouched
        apply_start_likelihood(d, ties, reads, cands, models, nlls_by_read=nlls)
        assert nlls[0] == {0: pytest.approx(-np.log(0.019)), 1: 0.0}   # peak: S best
        assert nlls[1] == {0: 0.0}                                       # upstream: P
        assert nlls[2] == {0: 1.0, 1: 3.0}
        assert 3 not in nlls and 4 not in nlls
        # the recovered model is now some read's best -> the support gate keeps it
        assert support_gate_drops(cands[:2], {k: v for k, v in ties.items() if k < 3},
                                  nlls) == set()
        # and a scored junction contrast keeps its values plus the penalty on P
        d2, ties2, reads2, cands2, models2, _ = self._setup()
        nlls2 = {0: {0: 0.5, 1: 0.7}}
        apply_start_likelihood(d2, ties2, reads2, cands2, models2, nlls_by_read=nlls2)
        assert nlls2[0] == {0: pytest.approx(0.5 - np.log(0.019)), 1: 0.7}

    def test_positional_term_is_not_divided_by_sigma2(self):
        # the term must be applied after junction sigma2 normalisation: an assign()
        # with m2_metric=off (no junction distances) yields exactly -log(h) on the
        # parent cell, not 1.0 (which is what median-normalising it would give).
        from fin.pipeline import assignment as asg
        captured = {}
        orig = asg.em_with_coherence

        def spy(**kw):
            captured["d"] = np.array(kw["dist_read_to_tx"])
            return orig(**kw)

        cs, ids = _build(SHARP)
        d, ties, reads, cands, models, MISSING = self._setup()
        cs.candidates.append(cands[1])
        models[0]["read_spans"] = cs.read_spans
        runner, _ = _runner(Path("/tmp"), "off")
        from unittest.mock import patch
        with patch.object(asg, "em_with_coherence", side_effect=spy):
            qo = runner._assigner.assign(
                cs, ids, tie_nll_fn=runner._tie_nll, eff_lengths_fn=runner._eff_lengths,
                start_likelihood_context=models)
        jp = [c.candidate_id for c in qo.cand_list].index("parent")
        i = qo.kept_read_ids.index("sh0")
        assert captured["d"][i, jp] == pytest.approx(-np.log(0.019), rel=1e-5)


class TestHardAssignTieBreak:
    """The explicit hard-assignment rule (fin/pipeline/assignment.py::hard_assign).

    Pins that assignment policy is NOT carried by candidate list order.
    """

    class _C:
        def __init__(self, cid, source):
            self.candidate_id, self.source = cid, source

    def _cands(self):
        return [self._C("novel_z", "novel"), self._C("ENST1", "gtf"),
                self._C("novel_a", "novel"), self._C("ENST0", "gtf")]

    def test_exact_tie_prefers_annotation_then_candidate_id(self):
        from fin.pipeline.assignment import hard_assign
        R = np.array([
            [0.5, 0.5, 0.0, 0.0],    # novel vs annotation -> annotation
            [0.9, 0.05, 0.05, 0.0],  # clear winner, no tie-break involved
            [0.0, 0.5, 0.0, 0.5],    # two annotations -> lower candidate_id
            [0.5, 0.0, 0.5, 0.0],    # two novels -> lower candidate_id
        ])
        hard, n_ties, n_anno = hard_assign(R, self._cands())
        assert hard.tolist() == [1, 0, 3, 2]
        assert n_ties == 3 and n_anno == 2

    def test_result_is_invariant_to_candidate_order(self):
        from fin.pipeline.assignment import hard_assign
        cands = self._cands()
        R = np.array([[0.5, 0.5, 0.0, 0.0], [0.0, 0.5, 0.0, 0.5],
                      [0.5, 0.0, 0.5, 0.0]])
        hard, _, _ = hard_assign(R, cands)
        for perm in ([2, 0, 3, 1], [3, 2, 1, 0]):
            h2, _, _ = hard_assign(R[:, perm], [cands[j] for j in perm])
            assert ([cands[perm[j]].candidate_id for j in h2]
                    == [cands[j].candidate_id for j in hard])

    def test_positive_evidence_is_never_overridden_by_annotation(self):
        from fin.pipeline.assignment import hard_assign
        R = np.array([[0.6, 0.4, 0.0, 0.0]])          # novel strictly better
        hard, n_ties, _ = hard_assign(R, self._cands())
        assert hard.tolist() == [0] and n_ties == 0

    def test_empty_matrix(self):
        from fin.pipeline.assignment import hard_assign
        hard, n_ties, n_anno = hard_assign(np.zeros((0, 3)), self._cands()[:3])
        assert hard.shape == (0,) and n_ties == 0 and n_anno == 0


def test_audit_never_mutates_the_candidate_set(tmp_path):
    """The evaluator is pure. An already-existing recovered id must be reported
    as a merge, not applied -- support counts feed gates that are on by default,
    so mutating in audit mode would break byte-inertness."""
    from fin.analysis.tss_shadow_recovery import recover_folded_tss_candidates

    class _QO:
        def __init__(self, cand_list, kept, R):
            self.cand_list, self.kept_read_ids, self.R = cand_list, kept, R

    cs, ids = _build(SHARP)
    parent = cs.candidates[0]
    before = {c.candidate_id: set(c.supporting_read_ids) for c in cs.candidates}
    qo = _QO([parent], ids, np.ones((len(ids), 1)))
    out = recover_folded_tss_candidates(cs, qo, _interval(), GENOME)
    # the evaluator returned candidates/merges but changed nothing in place
    assert {c.candidate_id: set(c.supporting_read_ids) for c in cs.candidates} == before
    assert len(cs.candidates) == 1


def test_existing_recovered_id_is_returned_as_a_merge(tmp_path):
    from fin.analysis.tss_shadow_recovery import recover_folded_tss_candidates

    class _QO:
        def __init__(self, cand_list, kept, R):
            self.cand_list, self.kept_read_ids, self.R = cand_list, kept, R

    cs, ids = _build(SHARP)
    parent = cs.candidates[0]
    qo = _QO([parent], ids, np.ones((len(ids), 1)))
    first = recover_folded_tss_candidates(cs, qo, _interval(), GENOME)
    assert first.candidates, "expected a recovered candidate to exist"
    # now pretend it is already in the set: the second pass must MERGE, not mutate
    cs.candidates.append(first.candidates[0])
    before = set(first.candidates[0].supporting_read_ids)
    second = recover_folded_tss_candidates(cs, qo, _interval(), GENOME)
    assert second.candidates == []
    assert second.support_merges and set(first.candidates[0].supporting_read_ids) == before
