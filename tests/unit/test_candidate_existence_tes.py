"""The 3'/polyA channel of the existence test (spec A7 / E3).

Pinned here: strand-aware distances; PASS / FAIL / MISSING semantics; a
polyA read ending at a candidate's TES is evidence FOR it; a polyA read that
reads THROUGH a candidate's claimed TES is evidence AGAINST it; genuine APA
isoforms coexist; and the synthetic c_merge component -- where the structural
channel alone cannot refute a merged annotation -- is separated by this channel.
"""
from __future__ import annotations

import numpy as np
import pytest

from fin.analysis.candidate_existence import (
    build_conflict_components, evaluate_candidate_existence,
)
from fin.analysis.molecule_evidence import MoleculeEvidence
from fin.analysis.sample_nuisance import (
    EndpointModel, SampleNuisanceModel, endpoint_distances, fit_endpoint_model,
    fit_sample_nuisance, tes_distance, tes_polya_loglik,
)
from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval

MISSING = 1e6


def _cand(cid, start, end, introns, strand="+"):
    return TranscriptCandidate(
        candidate_id=cid, intron_chain=IntronChain(introns=introns),
        three_prime_pos=end if strand == "+" else start, sequence="ACGT" * 50,
        source="gtf", supporting_read_ids=set(), chrom="chr1", strand=strand,
        start=start, end=end, family_id="f")


def _mol(rid, three, strand="+", qc="PASS", tie=("A", "B")):
    rs, re = (100, three) if strand == "+" else (three, 5000)
    return MoleculeEvidence(
        read_id=rid, strand=strand, alignment_scores={c: 10.0 for c in tie},
        raw_tie_ids=tuple(tie), cigar_junctions=(), indel_spans=(),
        reference_start=rs, reference_end=re,
        five_prime_pos=rs if strand == "+" else re,
        three_prime_pos=three, m2_nll={}, covered_diff_regions={},
        polya_qc=qc, polya_length=40.0 if qc == "PASS" else None)


def _model():
    return EndpointModel(b_narrow=1.0, b_broad=80.0, eta=0.04, n_fit=1000)


class TestEndpointModel:
    def test_fit_recovers_mixture(self):
        import random
        rng = random.Random(3)
        d = []
        for _ in range(4000):
            b = 60 if rng.random() < 0.06 else 1.0
            d.append(int(round(rng.expovariate(1 / b))) * rng.choice([-1, 1]))
        m = fit_endpoint_model(d)
        assert 0.6 < m.b_narrow < 1.6
        assert m.b_broad > 20
        assert 0.02 < m.eta < 0.12

    def test_split_half_distances_measure_against_the_other_half(self):
        # 8 reads at ~100, split alternately: each half defines a mode at 100 and
        # the other half is measured against it (small distances); a lone far end
        # is measured against the other half's mode (large distance).
        ends = {("c", "+"): [100, 100, 101, 100, 101, 100, 100, 101, 500]}
        d = endpoint_distances(ends)
        assert d and max(abs(x) for x in d) >= 399
        assert sum(1 for x in d if abs(x) <= 1) >= 6

    def test_minus_strand_flips_sign(self):
        # on -, transcription runs toward lower coordinates: an end BELOW the
        # mode is "past" it -> positive distance. Split-half: the mode must be
        # defined by >= 2 reads of the other half.
        d_plus = endpoint_distances({("c", "+"): [100, 100, 100, 100, 130]})
        d_minus = endpoint_distances({("c", "-"): [100, 100, 100, 100, 70]})
        assert 30 in d_plus and 30 in d_minus

    def test_claimed_tes_kernel_uses_the_annotation(self):
        pm = {f"r{i}": (30.0, "PASS") for i in range(30)}
        ends = {f"r{i}": ("c", "+", 1003 + (i % 3)) for i in range(30)}   # cluster ~1004
        m = fit_sample_nuisance(pm, ends, claimed_tes={("c", "+"): [1000]})
        # distances 3..5 to the CLAIMED TES are what the kernel must absorb
        assert m.endpoint.n_fit == 30 and m.endpoint.b_narrow >= 2.0

    def test_readthrough_is_far_less_likely_than_landing(self):
        m = _model()
        assert m.loglik(0) > m.loglik(6) > m.loglik(60)
        assert m.narrow_posterior(0) > 0.99 and m.narrow_posterior(60) < 0.01


class TestPerReadTerm:
    def test_strand_aware_distance(self):
        a_plus = _cand("A", 100, 2000, ((500, 900),))
        a_minus = _cand("A", 100, 2000, ((500, 900),), strand="-")
        assert tes_distance(_mol("r", 2030), a_plus) == 30            # past on +
        assert tes_distance(_mol("r", 70, strand="-"), a_minus) == 30  # past on -
        assert tes_distance(_mol("r", 1980), a_plus) == -20

    def test_fail_and_missing_abstain(self):
        m = _model(); a = _cand("A", 100, 2000, ((500, 900),))
        assert tes_polya_loglik(_mol("r", 2000, qc="FAIL"), a, m) is None
        assert tes_polya_loglik(_mol("r", 2000, qc="MISSING"), a, m) is None
        assert tes_polya_loglik(_mol("r", 2000), a, m) is not None

    def test_unavailable_model_abstains(self):
        m = EndpointModel(b_narrow=1, b_broad=80, eta=0.04, n_fit=1)
        assert tes_polya_loglik(_mol("r", 2000), _cand("A", 100, 2000, ((500, 900),)), m) is None


class _QO:
    def __init__(self, cands, kept, kernel, ties, nlls, cover):
        self.cand_list, self.kept_read_ids = cands, kept
        self.structural_d_tx_raw = kernel
        self.ties_by_read, self.nlls_by_read, self.cover_by_read = ties, nlls, cover


def _component(cands, molecules, kernel_rows):
    kept = [m.read_id for m in molecules]
    n, k = len(kept), len(cands)
    kernel = np.full((n, k), MISSING)
    for i, row in enumerate(kernel_rows):
        for j, v in row.items():
            kernel[i, j] = v
    ties = {i: [j for j in range(k) if kernel[i, j] < MISSING] for i in range(n)}
    nlls = {i: {j: kernel[i, j] for j in ties[i]} for i in range(n)}
    cover = {i: True for i in range(n)}
    return _QO(cands, kept, kernel, ties, nlls, cover)


def _interval():
    return GenomicInterval(chrom="chr1", start=0, end=6000, strand="+")


def _run(qo, molecules, nuisance):
    comps = build_conflict_components(qo, None)
    rows = evaluate_candidate_existence(
        qo, None, comps, _interval(), sample_id="s",
        molecules={m.read_id: m for m in molecules}, nuisance=nuisance)
    return {r.candidate_id: r for r in rows}


class TestMergeComponent:
    """c_merge: annotation M spans two real transcripts A (TES 2000) and B
    (TSS after A's TES, TES 4000). Structurally M explains both A's and B's
    junction evidence, so the structural channel cannot refute it. But every
    polyA read of A ends at 2000 -- reading THROUGH M's claimed TES (4000) is
    exactly what M predicts and never happens; ending at 2000 is a read-through
    of nothing for A and a 2000 bp shortfall for M."""

    def _setup(self, polya=True):
        A = _cand("A", 100, 2000, ((500, 900),))
        M = _cand("M", 100, 4000, ((500, 900), (2100, 2500)))
        mols, rows = [], []
        for i in range(12):                               # A's reads: end at ~2000
            mols.append(_mol(f"a{i}", 2000 + (i % 3) - 1, qc="PASS" if polya else "FAIL",
                             tie=("A", "M")))
            rows.append({0: 0.0, 1: 0.0})                 # structurally indistinguishable
        nuisance = SampleNuisanceModel(endpoint=_model(), polya_call_rate=0.7, n_scorable=100)
        return _component([A, M], mols, rows), mols, nuisance

    def test_structural_channel_cannot_separate_merge(self):
        qo, mols, nuisance = self._setup(polya=False)     # all abstain on 3'
        rows = _run(qo, mols, nuisance)
        assert abs(rows["A"].delta_predictive_ll) < 1e-6
        assert abs(rows["M"].delta_predictive_ll) < 1e-6
        assert rows["M"].n_polya_pass == 0 and rows["M"].tes_polya_llr == pytest.approx(0.0)

    def test_polya_channel_refutes_the_merge(self):
        qo, mols, nuisance = self._setup(polya=True)
        rows = _run(qo, mols, nuisance)
        assert rows["A"].delta_predictive_ll > 0          # removing A hurts: nothing else ends at 2000
        assert rows["M"].delta_predictive_ll <= 1e-9      # M adds nothing the data wants
        assert rows["A"].tes_polya_llr > 0
        assert rows["M"].n_polya_pass == 12 and rows["M"].n_tes_near == 0
        assert rows["A"].n_tes_near == 12 and rows["A"].n_tes_readthrough == 0
        assert rows["M"].negative_channel == "" or rows["M"].n_tes_readthrough == 0

    def test_readthrough_is_recorded_as_negative_channel(self):
        # a short claimed TES with polyA reads running past it
        S = _cand("S", 100, 1500, ((500, 900),))
        L = _cand("L", 100, 2000, ((500, 900),))
        mols = [_mol(f"r{i}", 2000 + (i % 2), tie=("S", "L")) for i in range(10)]
        qo = _component([S, L], mols, [{0: 0.0, 1: 0.0}] * 10)
        rows = _run(qo, mols, SampleNuisanceModel(endpoint=_model(), polya_call_rate=0.7, n_scorable=50))
        assert rows["S"].n_tes_readthrough == 10
        assert rows["S"].negative_channel == "tes_readthrough_polya"
        assert rows["S"].delta_predictive_ll <= 1e-9 and rows["L"].delta_predictive_ll > 0


class TestApaCoexistence:
    def test_two_real_polya_sites_both_supported(self):
        # genuine alternative polyadenylation: half the reads end at 2000, half at 3000
        P = _cand("P", 100, 2000, ((500, 900),))
        Q = _cand("Q", 100, 3000, ((500, 900),))
        mols = ([_mol(f"p{i}", 2000 + (i % 2), tie=("P", "Q")) for i in range(10)]
                + [_mol(f"q{i}", 3000 - (i % 2), tie=("P", "Q")) for i in range(10)])
        qo = _component([P, Q], mols, [{0: 0.0, 1: 0.0}] * 20)
        rows = _run(qo, mols, SampleNuisanceModel(endpoint=_model(), polya_call_rate=0.7, n_scorable=50))
        assert rows["P"].delta_predictive_ll > 0 and rows["Q"].delta_predictive_ll > 0
        assert rows["P"].n_tes_near == 10 and rows["Q"].n_tes_near == 10


class TestSampleNuisance:
    def test_fit_from_calls_and_ends_with_interval_leave_out(self):
        pm = {f"r{i}": (30.0, "PASS") for i in range(120)}
        pm.update({f"f{i}": (5.0, "FAIL") for i in range(30)})
        ends = {f"r{i}": ("c", "+", 1000 + (i % 3) - 1) for i in range(60)}
        ends.update({f"r{i}": ("c", "+", 5000 + (i % 3) - 1) for i in range(60, 120)})
        m = fit_sample_nuisance(pm, ends)
        assert m.polya_available and m.polya_call_rate == pytest.approx(0.8)
        assert m.endpoint.n_fit == 120
        m2 = fit_sample_nuisance(pm, ends, exclude_reads=[f"r{i}" for i in range(60)])
        assert m2.endpoint.n_fit == 60                 # the left-out interval is absent

    def test_no_calls_means_channel_unavailable(self):
        m = fit_sample_nuisance({}, {})
        assert not m.polya_available
