"""Compatibility attestation: exact-or-3'-suffix read depth and the support-gate exemption."""

from __future__ import annotations

import numpy as np
import pytest

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.candidates.discovery import chain_compatible_read_depths
from fin.pipeline.config import PipelineConfig


def _gtf(cid, introns, strand="+", start=0, end=1000):
    return TranscriptCandidate(
        candidate_id=cid, intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=end if strand == "+" else start, sequence="ACGT", source="gtf",
        supporting_read_ids=set(), chrom="chr1", strand=strand, start=start, end=end,
    )


# '+' model with three introns; exons [0,100) [200,300) [400,500) [600,1000)
PLUS = ((100, 200), (300, 400), (500, 600))
# '-' model: same genomic layout, but its 3' end is the genomic START
MINUS = PLUS


def _depths(models, reads, indels=None):
    gm = {}
    for m in models:
        gm.setdefault(tuple(m.intron_chain.introns), []).append(m)
    return chain_compatible_read_depths(gm, reads, indels or {})


class TestPlusStrand:
    def test_exact_and_suffix_reads_count_once_each(self):
        m = _gtf("A", PLUS)
        reads = [("r1", "+", PLUS, 10, 980),                 # exact
                 ("r2", "+", PLUS[1:], 250, 990),            # 2-intron suffix (5' degraded)
                 ("r3", "+", PLUS[2:], 450, 1000),           # 1-intron suffix
                 ("r3", "+", PLUS[2:], 450, 1000)]           # duplicate id: not double-counted
        assert _depths([m], reads) == {"A": 3}

    def test_prefix_and_internal_subchains_are_not_suffixes(self):
        m = _gtf("A", PLUS)
        reads = [("p", "+", PLUS[:2], 10, 550),              # 5' prefix: 3' end is degraded, not 5'
                 ("i", "+", PLUS[1:2], 250, 450)]            # internal single intron
        assert _depths([m], reads) == {}

    def test_read_running_across_the_upstream_intron_is_incompatible(self):
        m = _gtf("A", PLUS)
        # chain is the 2-intron suffix but the read starts at 150, i.e. it is
        # unspliced across intron (100,200): a retained intron, not a fragment
        reads = [("r", "+", PLUS[1:], 150, 990)]
        assert _depths([m], reads) == {}
        # starting inside the acceptor's alignment slop is allowed
        reads = [("r", "+", PLUS[1:], 195, 990)]
        assert _depths([m], reads) == {"A": 1}

    def test_three_prime_end_may_overshoot_only_within_tolerance(self):
        m = _gtf("A", PLUS, end=1000)
        assert _depths([m], [("r", "+", PLUS, 10, 1040)]) == {"A": 1}
        assert _depths([m], [("r", "+", PLUS, 10, 1100)]) == {}

    def test_longer_read_chain_never_attests_shorter_model(self):
        short = _gtf("S", PLUS[1:], start=250, end=1000)
        reads = [("r", "+", PLUS, 10, 990)]
        assert _depths([short], reads) == {}

    def test_shared_suffix_attests_every_nested_model(self):
        long_ = _gtf("L", PLUS)
        short = _gtf("S", PLUS[1:], start=250, end=1000)
        reads = [("r", "+", PLUS[1:], 260, 990)]         # S's exact chain, a suffix of L
        assert _depths([long_, short], reads) == {"L": 1, "S": 1}

    def test_early_three_prime_end_inside_last_exon_is_compatible(self):
        m = _gtf("A", PLUS, end=1000)
        assert _depths([m], [("r", "+", PLUS, 10, 700)]) == {"A": 1}

    def test_wrong_strand_and_empty_chains_do_not_count(self):
        m = _gtf("A", PLUS)
        assert _depths([m], [("r", "-", PLUS, 10, 990)]) == {}
        assert _depths([m], [("r", "+", (), 10, 990)]) == {}
        assert _depths([_gtf("M", ())], [("r", "+", (), 10, 990)]) == {}

    def test_indel_policy_is_off_by_default_and_opt_in(self):
        m = _gtf("A", PLUS)
        reads = [("r", "+", PLUS[2:], 450, 1000)]
        indels = {"r": ((505, 507),)}                                  # 5 bp from donor 500
        assert _depths([m], reads, indels) == {"A": 1}                 # default: not consulted
        gm = {PLUS: [m]}
        assert chain_compatible_read_depths(gm, reads, indels, indel_flank_bp=10) == {}
        assert chain_compatible_read_depths(gm, reads, {"r": ((450, 452),)}, indel_flank_bp=10) == {"A": 1}


class TestMinusStrand:
    def test_suffix_is_the_genomic_head(self):
        m = _gtf("A", MINUS, strand="-", start=0, end=1000)
        reads = [("r1", "-", MINUS, 10, 990),                # exact
                 ("r2", "-", MINUS[:2], 5, 450),             # 3' suffix = first two introns
                 ("r3", "-", MINUS[:1], 0, 250)]             # 3' suffix = first intron
        assert _depths([m], reads) == {"A": 3}

    def test_genomic_tail_is_a_five_prime_prefix_on_minus(self):
        m = _gtf("A", MINUS, strand="-")
        assert _depths([m], [("r", "-", MINUS[1:], 250, 990)]) == {}

    def test_minus_read_across_the_upstream_intron_is_incompatible(self):
        m = _gtf("A", MINUS, strand="-")
        # 2-intron suffix on '-' = introns[:2]; upstream (5') intron is (500,600);
        # a read ending at 550 runs into it unspliced
        assert _depths([m], [("r", "-", MINUS[:2], 5, 550)]) == {}
        assert _depths([m], [("r", "-", MINUS[:2], 5, 505)]) == {"A": 1}

    def test_minus_three_prime_overshoot_tolerance(self):
        m = _gtf("A", MINUS, strand="-", start=100, end=1000)
        assert _depths([m], [("r", "-", MINUS, 60, 990)]) == {"A": 1}
        assert _depths([m], [("r", "-", MINUS, 20, 990)]) == {}


class TestDepthCollectionIsGated:
    """The per-read x per-model scan runs only when the exemption will read it.

    The end-to-end test below is the one that matters: a version where the
    runner passed the flag but `discover_candidates` never FORWARDED it to the
    clustering path silently disabled the whole exemption, and a runner-only
    mock did not catch it.
    """

    def test_the_flag_is_forwarded_not_merely_accepted(self, tmp_path, monkeypatch):
        """SOURCE-level regression: the parameter must be passed on, not just accepted.

        This is not an end-to-end execution test -- that evidence comes from the
        real-cell byte comparison in the docs. It exists because a version where
        the outer and inner signatures both had the parameter but the call site
        did not forward it silently disabled the entire exemption, and the
        runner-level mock below could not see that.
        """
        import fin.candidates.discovery as D
        seen = {}
        real = D.chain_compatible_read_depths
        monkeypatch.setattr(D, "chain_compatible_read_depths",
                            lambda *a, **k: seen.setdefault("called", True) or real(*a, **k))
        import inspect
        src = inspect.getsource(D.discover_candidates)
        # the parameter must be forwarded, not merely accepted
        assert "collect_compatible_depths=collect_compatible_depths" in src, (
            "discover_candidates accepts collect_compatible_depths but never "
            "passes it on: the exemption would be silently dead")
        inner = inspect.signature(D._chain_cluster_candidates).parameters
        assert "collect_compatible_depths" in inner
        assert inner["collect_compatible_depths"].default is False

    def test_runner_passes_the_flag(self, tmp_path):
        from unittest.mock import patch
        from fin.candidates.dataclasses import CandidateSet
        from fin.io.interval_manager import GenomicInterval
        from fin.pipeline.runner import PipelineRunner
        iv = GenomicInterval(chrom="chr1", start=0, end=1000, strand="+")
        seen = {}

        def fake_discover(**kw):
            seen["collect"] = kw.get("collect_compatible_depths")
            return CandidateSet(interval=iv, candidates=[], read_ids=set())

        for k, expected in ((0, False), (3, True)):
            cfg = PipelineConfig(bam_path="/fake.bam", work_dir=str(tmp_path),
                                 m2_support_gate_attest_reads=k)
            r = PipelineRunner(cfg)
            with patch("fin.pipeline.runner.discover_candidates", side_effect=fake_discover):
                r.process_interval(iv)
            assert seen["collect"] is expected


class TestConfig:
    def test_negative_rejected(self, tmp_path):
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        cfg = PipelineConfig(bam_path=str(bam), m2_support_gate_attest_reads=-1)
        with pytest.raises(ValueError, match="m2_support_gate_attest_reads"):
            cfg.validate()

    def test_default_off(self):
        assert PipelineConfig(bam_path="/fake.bam").m2_support_gate_attest_reads == 0


# --- the gate exemption -----------------------------------------------------

from types import SimpleNamespace  # noqa: E402

from fin.analysis.quantification import QuantResult  # noqa: E402
from fin.io.interval_manager import GenomicInterval  # noqa: E402
from fin.pipeline import selection as S  # noqa: E402


def _run_gate(cands, observed, attest_k, n_reads=3, veto_bp=10, veto="novel",
              guided_min_reads=0):
    """select_m2_interval with only the support gate active. Every read ties on
    all candidates by alignment score and column 0 has the strictly lowest
    junction NLL, so it wins M2 for each read and every other column is
    gate-dropped unless exempt."""
    n_c = len(cands)
    rids = [f"r{i}" for i in range(n_reads)]
    R = np.full((n_reads, n_c), 1.0 / n_c)
    qrs = [QuantResult(c.candidate_id, n_reads / n_c, 1.0, n_reads if j == 0 else 0, c.source,
                       assigned_read_ids=tuple(rids) if j == 0 else ())
           for j, c in enumerate(cands)]
    qo = SimpleNamespace(cand_list=list(cands), R=R, kept_read_ids=rids,
                         ties_by_read={r: list(range(n_c)) for r in rids},
                         nlls_by_read={r: {j: 1.0 + j for j in range(n_c)} for r in rids},
                         quant_results=qrs)
    cfg = SimpleNamespace(containment_evidence_mode="off", junction_election_mode="off",
                          m2_cluster_recheck=False, novel_junction_min_reads=0,
                          guided_junction_min_reads=guided_min_reads,
                          guided_junction_reads_tol=2, m2_support_gate=True,
                          m2_support_gate_tie=True, m2_support_gate_gtf_exempt=False,
                          m2_support_gate_attest_reads=attest_k, m2_support_gate_attest_veto_bp=veto_bp,
                          m2_support_gate_attest_veto=veto,
                          chain_cluster_wobble_bp=6,
                          complete_chain_rescue_reads=0,
                          containment_cluster=False, containment_collapse=False,
                          mono_resolve_post_em=False)
    cs = SimpleNamespace(read_spans={}, junction_support=observed)
    iv = GenomicInterval(chrom="chr1", start=0, end=1000, strand="+")
    # the guided gate reads the BAM-scan registry (strand-keyed), the
    # attestation veto reads candidate_set.junction_support
    obs_fn = lambda i: {"+": {j: (c[0] if isinstance(c, tuple) else c)
                             for j, c in (observed or {}).items()}}
    res, out = S.select_m2_interval(cfg, qo, cs, iv, obs_fn)
    return {r.candidate_id for r in res}, out


@pytest.fixture
def nested():
    long_ = _gtf("L", PLUS)
    short = _gtf("S", PLUS[1:], start=250, end=1000)
    long_.compatible_read_depth = 3
    short.compatible_read_depth = 3
    return long_, short


class TestGateExemption:
    OBS = {PLUS[0]: (2, 2), PLUS[1]: (5, 5), PLUS[2]: (5, 5)}   # (raw, clean) like junction_support

    def test_off_drops_the_undiscriminated_nested_model(self, nested):
        kept, outcomes = _run_gate(nested, self.OBS, attest_k=0)
        assert kept == {"L"}
        assert [(o.candidate_id, o.reason) for o in outcomes] == [("S", "m2_support_gate")]

    def test_attested_nested_model_is_kept(self, nested):
        kept, outcomes = _run_gate(nested, self.OBS, attest_k=3)
        assert kept == {"L", "S"}
        assert outcomes == []

    def test_depth_below_k_is_not_exempt(self, nested):
        nested[1].compatible_read_depth = 2
        kept, _ = _run_gate(nested, self.OBS, attest_k=3)
        assert kept == {"L"}

    def test_jittered_junction_vetoes_the_exemption(self, nested):
        # the reads' junction sits 3 bp beside S's first intron on the donor
        # side and exactly on the acceptor: S's junction is a jitter of it
        obs = {PLUS[0]: (2, 2), (303, 400): (5, 5), PLUS[2]: (5, 5)}
        kept, outcomes = _run_gate(nested, obs, attest_k=3)
        assert "S" not in kept
        assert any(o.candidate_id == "S" and o.reason == "m2_support_gate" for o in outcomes)

    def test_annotated_alternative_splice_site_is_not_a_contradiction(self, nested):
        # an observed junction sharing S's donor with an acceptor 60 bp away,
        # carried by ANOTHER annotated isoform: alternative splicing, not a
        # contradiction; S's own junction is unobserved (silence)
        long_, short = nested
        alt = _gtf("X", (PLUS[0], (300, 460), PLUS[2]))
        obs = {PLUS[0]: (2, 2), (300, 460): (5, 5), PLUS[2]: (5, 5)}
        kept, _ = _run_gate((long_, short, alt), obs, attest_k=3)
        assert "S" in kept

    def test_observed_own_junction_is_never_a_jitter(self, nested):
        # S's junction is itself observed; a 4 bp neighbour alongside it is a
        # wobble variant of the reads, not evidence against S
        obs = {PLUS[0]: (2, 2), PLUS[1]: (5, 5), (304, 400): (1, 1), PLUS[2]: (5, 5)}
        kept, _ = _run_gate(nested, obs, attest_k=3)
        assert kept == {"L", "S"}

    def test_veto_window_is_configurable(self, nested):
        obs = {PLUS[0]: (2, 2), (308, 400): (5, 5), PLUS[2]: (5, 5)}   # 8 bp off on the donor
        kept, _ = _run_gate(nested, obs, attest_k=3, veto_bp=6, veto="jitter")
        assert kept == {"L", "S"}
        kept, _ = _run_gate(nested, obs, attest_k=3, veto_bp=10, veto="jitter")
        assert kept == {"L"}

    def test_novel_candidate_never_exempt(self, nested):
        long_, short = nested
        short.source = "novel"
        kept, _ = _run_gate((long_, short), self.OBS, attest_k=3)
        assert kept == {"L"}

    def test_no_registry_means_no_exemption(self, nested):
        kept, _ = _run_gate(nested, {}, attest_k=3)
        assert kept == {"L"}
        kept, _ = _run_gate(nested, None, attest_k=3)
        assert kept == {"L"}


from fin.scoring.m2_junction_nll import unjittered_gtf_columns  # noqa: E402


class TestUnjitteredPredicate:
    def test_negative_window_returns_nothing(self):
        assert unjittered_gtf_columns([_gtf("A", PLUS)], {PLUS[0], PLUS[1], PLUS[2]}, window_bp=-1) == set()

    def test_masked_novel_and_mono_columns_are_skipped(self):
        cands = [None, _gtf("N", PLUS), _gtf("M", ()), _gtf("A", PLUS)]
        cands[1].source = "novel"
        assert unjittered_gtf_columns(cands, {PLUS[0]}, window_bp=10) == {3}

    def test_window_is_inclusive_on_both_boundaries(self):
        a = _gtf("A", ((100, 200),))
        assert unjittered_gtf_columns([a], {(110, 210)}, window_bp=10) == set()
        assert unjittered_gtf_columns([a], {(111, 210)}, window_bp=10) == {0}
        assert unjittered_gtf_columns([a], {(110, 211)}, window_bp=10) == {0}


class TestNovelJunctionVeto:
    """Rule 2: an observed junction in NO annotated chain, at the place the
    candidate differs, contradicts it; an annotated alternative does not."""

    def _pair(self):
        # A: intact isoform [J0, J1, J2]; M: merge-corrupted entry whose middle
        # junction (300, 700) joins A's donor to a downstream exon, so that A's
        # real junction (300, 400) is annotated NOWHERE once A is replaced.
        m = _gtf("M", ((100, 200), (300, 700), (800, 900)), start=0, end=1000)
        m.compatible_read_depth = 5
        return m

    def test_unannotated_junction_sharing_the_donor_vetoes(self):
        m = self._pair()
        obs = {(100, 200), (300, 400), (800, 900)}       # reads splice 300->400: unannotated
        assert unjittered_gtf_columns([m], obs, 10, novel_junction_veto=True) == set()
        assert unjittered_gtf_columns([m], obs, 10, novel_junction_veto=False) == {0}

    def test_annotated_alternative_is_not_a_contradiction(self):
        m = self._pair()
        alt = _gtf("X", ((100, 200), (300, 400), (800, 900)))   # the 300->400 junction IS annotated
        obs = {(100, 200), (300, 400), (800, 900)}
        assert unjittered_gtf_columns([m, alt], obs, 10, novel_junction_veto=True) == {0, 1}

    def test_unannotated_junction_inside_a_candidate_exon_vetoes(self):
        # IR-corrupted entry: exon [200, 700) spans an intron the reads splice
        ir = _gtf("R", ((100, 200), (700, 800)), start=0, end=1000)
        obs = {(100, 200), (300, 500), (700, 800)}
        assert unjittered_gtf_columns([ir], obs, 10, novel_junction_veto=True) == set()
        spliced = _gtf("S", ((100, 200), (300, 500), (700, 800)))
        assert unjittered_gtf_columns([ir, spliced], obs, 10, novel_junction_veto=True) == {0, 1}

    def test_masked_annotation_still_counts_as_annotated(self):
        """The annotation universe is the SUPPLIED annotation, not the survivors."""
        m = self._pair()
        alt = _gtf("X", ((100, 200), (300, 400), (800, 900)))
        obs = {(100, 200), (300, 400), (800, 900)}
        live = [m, None]                      # X was dropped by an earlier gate
        universe = {j for c in (m, alt) for j in c.intron_chain.introns}
        assert unjittered_gtf_columns(live, obs, 10, novel_junction_veto=True,
                                      annotated_junctions=universe) == {0}
        # ... and without the universe the survivors alone would veto M
        assert unjittered_gtf_columns(live, obs, 10, novel_junction_veto=True) == set()

    def test_malformed_geometry_is_never_exempt(self):
        # span does not bracket the chain: no well-defined exons
        bad = _gtf("B", ((100, 200), (300, 400)), start=250, end=1000)
        obs = {(100, 200), (600, 700)}
        assert unjittered_gtf_columns([bad], obs, 10, novel_junction_veto=True) == set()
        good = _gtf("G", ((100, 200), (300, 400)), start=0, end=1000)
        assert unjittered_gtf_columns([good], {(100, 200), (300, 400)}, 10,
                                      novel_junction_veto=True) == {0}

    def test_annotation_universe_is_built_from_every_annotated_column(self, nested, monkeypatch):
        """The universe must come from `cand_list`, not from the live survivors.

        This is the contract the P1 bug broke: a column an EARLIER gate put in
        `drop_cols` is absent from `live`, and if the universe were derived
        from `live` the reads' junction would look novel and veto the very
        candidate the exemption exists for. The guided junction-support gate
        runs before the attestation block and drops annotated columns here, so
        `live` is genuinely shorter than `cand_list` when the veto is called.
        """
        import fin.scoring.m2_junction_nll as M
        long_, short = nested
        alt = _gtf("X", (PLUS[0], (300, 460), PLUS[2]))
        seen = {}
        real = M.unjittered_gtf_columns

        def spy(cands, obs, window_bp=10, novel_junction_veto=False, annotated_junctions=None):
            seen["live_ids"] = [c.candidate_id if c is not None else None for c in cands]
            seen["universe"] = set(annotated_junctions or ())
            return real(cands, obs, window_bp, novel_junction_veto, annotated_junctions)

        monkeypatch.setattr(M, "unjittered_gtf_columns", spy)
        obs = {PLUS[0]: (2, 2), (300, 460): (5, 5), PLUS[2]: (5, 5)}
        _run_gate((long_, short, alt), obs, attest_k=3, guided_min_reads=4)
        assert seen["live_ids"] == [None, None, None], "the earlier gate must have masked columns"
        assert (300, 460) in seen["universe"], "a masked annotation's junction is still annotated"
        assert seen["universe"] == {PLUS[0], PLUS[1], PLUS[2], (300, 460)}

    def test_without_the_universe_the_masked_annotation_would_veto(self, nested):
        """The control: survivors alone are not the annotation universe."""
        long_, short = nested
        alt = _gtf("X", (PLUS[0], (300, 460), PLUS[2]))
        obs_set = {PLUS[0], (300, 460), PLUS[2]}
        live = [long_, short, None]                       # X masked by an earlier gate
        assert unjittered_gtf_columns(live, obs_set, 10, novel_junction_veto=True) == set()
        universe = {j for c in (long_, short, alt) for j in c.intron_chain.introns}
        assert 1 in unjittered_gtf_columns(live, obs_set, 10, novel_junction_veto=True,
                                           annotated_junctions=universe)

    def test_gate_uses_novel_veto_by_default(self, nested):
        # S's first junction is unobserved and the reads show an UNANNOTATED
        # junction sharing its donor: S is vetoed under "novel", kept under "jitter"
        obs = {PLUS[0]: (2, 2), (300, 460): (5, 5), PLUS[2]: (5, 5)}
        kept, _ = _run_gate(nested, obs, attest_k=3, veto="novel")
        assert kept == {"L"}
        kept, _ = _run_gate(nested, obs, attest_k=3, veto="jitter")
        assert kept == {"L", "S"}


# --- the endpoint-contradiction magnitude floors ---------------------------


class TestContradictionFloors:
    """A read end within the floor never contradicts the annotated end.

    The floor exists because the contradiction threshold comes from the NARROW
    component of the endpoint mixture (~2 bp on human) while the measured
    dispersion of dRNA read starts around an annotated TSS is 20-50 bp; see
    `experiments/prod_validation/annotation_end_precision.py`.
    """

    def test_defaults_are_zero_and_negatives_rejected(self, tmp_path):
        cfg = PipelineConfig(bam_path="/fake.bam")
        assert cfg.tss_contradiction_min_bp == 0 and cfg.tes_contradiction_min_bp == 0
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        for f in ("tss_contradiction_min_bp", "tes_contradiction_min_bp"):
            with pytest.raises(ValueError, match=f):
                PipelineConfig(bam_path=str(bam), **{f: -1}).validate()

    def _closure(self, work, **cfg_kwargs):
        """The real `contradicted` closure over a synthetic sample evidence."""
        import json
        from fin.pipeline.runner import PipelineRunner
        work.mkdir(parents=True, exist_ok=True)
        cfg = PipelineConfig(bam_path="/fake.bam", work_dir=str(work),
                             assignment_channels="ends", signal_path="/fake.blow5",
                             **cfg_kwargs)
        r = PipelineRunner(cfg)
        # 120 molecules whose ends sit exactly on the annotated ends, so the
        # fitted kernels are narrow and any offset reads as a contradiction
        polya = {f"r{i}": [40.0, "PASS"] for i in range(120)}
        ends = {f"r{i}": ["chr1", "+", 2000, 1000] for i in range(120)}
        (r._evidence_path()).write_text(json.dumps(
            {"polya": polya, "ends": ends, "version": "sample-evidence-v3",
             "identity": r._evidence_identity(), "polya_status": "ok",
             "n_ends": len(ends), "n_polya": len(polya)}))
        r._annotated_tes = lambda: {("chr1", "+"): [2000]}
        r._annotated_tss = lambda: {("chr1", "+"): [1000]}
        return r

    def test_floor_suppresses_a_small_upstream_cluster(self, tmp_path, monkeypatch):
        gc = _gtf("A", ((1200, 1400),), start=1000, end=2000)
        reads = [f"u{i}" for i in range(4)]

        def run(tag, **kw):
            r = self._closure(tmp_path / tag, **kw)
            fn = r._read_extent_contradicted()
            if fn is None:
                pytest.skip("endpoint channels unavailable in this environment")
            return fn(gc, reads)

        # without a floor a 15 bp upstream cluster contradicts the annotated TSS
        base = run("base")
        floored = run("floored", tss_contradiction_min_bp=20)
        # the floor can only ever REMOVE tss contradictions, never add them
        assert set(floored.get("tss", ())) <= set(base.get("tss", ()))
        assert floored.get("tes", []) == base.get("tes", [])


class TestTssTrackingStatistic:
    """The per-sample measurement that gates the 5' extent contradiction.

    It answers "does THIS run's 5' end resolve annotated starts at all". The
    discriminating quantity is the TAIL, not the median offset: measured on 10
    libraries the median offset is nearly the same in both domains (14 vs 22 bp,
    the motor/degradation constant) while the tracked FRACTION is bimodal --
    human dRNA 0.491-0.514, SIRV4 0.797-0.814, Sequin 0.917-0.950.
    """

    def _stat(self, starts, tss, **kw):
        from fin.analysis.sample_nuisance import annotated_tss_tracking
        kw.setdefault("min_mode_share", 1.1)   # off unless a test asks for it
        return annotated_tss_tracking({("c", "+"): sorted(starts)},
                                      {("c", "+"): tss}, **kw)

    def test_a_mode_on_the_tss_tracks_it(self):
        frac, n = self._stat([1000, 1002, 1003, 1150], [1000])
        assert n == 1 and frac == 1.0

    def test_a_mode_far_from_the_tss_does_not(self):
        # the dense cluster is 120 bp downstream; the TSS itself has one read
        frac, n = self._stat([1000] + [1120, 1121, 1122, 1123], [1000])
        assert n == 1 and frac == 0.0

    def test_a_DOMINANT_displaced_mode_counts_as_resolved(self):
        """A tight cluster away from the annotated TSS is evidence against the
        ANNOTATION, not against the reads -- the confound that made a
        tracking-only gate switch the 5' repair off on a truncated annotation."""
        starts = [1000] + [1120, 1121, 1122, 1123]
        assert self._stat(starts, [1000], min_mode_share=0.5)[0] == 1.0
        assert self._stat(starts, [1000], min_mode_share=1.1)[0] == 0.0

    def test_a_WEAK_displaced_mode_does_not(self):
        # starts smeared over 200 bp: no cluster holds a majority
        starts = list(range(1100, 1300, 10))
        assert self._stat(starts, [1000], min_mode_share=0.5)[0] == 0.0

    def test_the_mode_is_used_not_the_mean(self):
        # a long degradation tail must not drag the estimate off the TSS
        frac, n = self._stat([1000, 1001, 1002] + list(range(1100, 1180, 10)), [1000])
        assert frac == 1.0

    def test_a_locus_below_min_reads_is_not_scored(self):
        frac, n = self._stat([1000, 1001], [1000], min_reads=3)
        assert n == 0 and frac == 0.0

    def test_nothing_scorable_returns_zero_loci_not_zero_fraction_meaning(self):
        from fin.analysis.sample_nuisance import annotated_tss_tracking
        assert annotated_tss_tracking({}, {("c", "+"): [1]}) == (0.0, 0)
        assert annotated_tss_tracking({("c", "+"): [1]}, {}) == (0.0, 0)

    def test_track_bp_is_the_tolerance(self):
        starts = [1030, 1031, 1032]
        assert self._stat(starts, [1000], track_bp=25)[0] == 0.0
        assert self._stat(starts, [1000], track_bp=50)[0] == 1.0

    def test_fraction_is_over_scored_loci_only(self):
        frac, n = self._stat([1000, 1001, 1002, 5000, 5001, 5002, 9000],
                             [1000, 5100, 9000])
        # locus 1000 tracks, locus 5100 is scored (3 reads within 200) and does
        # not track, locus 9000 has one read and is not scored
        assert n == 2 and abs(frac - 0.5) < 1e-9


class TestTssAutoGate:
    def test_defaults_and_validation(self, tmp_path):
        cfg = PipelineConfig(bam_path="/fake.bam")
        assert cfg.tss_contradiction_auto_gate is False
        assert cfg.tss_contradiction_track_frac == 0.89
        assert cfg.tss_contradiction_mode_share == 0.5
        assert cfg.tss_contradiction_track_bp == 25
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        with pytest.raises(ValueError, match="track_frac"):
            PipelineConfig(bam_path=str(bam), tss_contradiction_track_frac=1.5).validate()
        with pytest.raises(ValueError, match="track_bp"):
            PipelineConfig(bam_path=str(bam), tss_contradiction_track_bp=0).validate()

    def test_reflection_invariance(self):
        """The same geometry must score the same on either strand.

        A first-maximum tie-break does not: it makes the answer depend on which
        end of the window the scan starts from.
        """
        from fin.analysis.sample_nuisance import annotated_tss_tracking
        starts = [1000, 1001, 1002, 1060, 1061, 1062]
        t = 1000
        a = annotated_tss_tracking({("c", "+"): tuple(sorted(starts))},
                                   {("c", "+"): [t]}, min_mode_share=1.1)
        M = 10_000
        mirrored = sorted(M - x for x in starts)
        b = annotated_tss_tracking({("c", "-"): tuple(mirrored)},
                                   {("c", "-"): [M - t]}, min_mode_share=1.1)
        assert a == b

    def test_shared_annotated_tss_votes_once(self):
        """Isoforms of one gene share a TSS; it must not get extra weight."""
        from fin.analysis.sample_nuisance import annotated_tss_tracking
        starts = {("c", "+"): (1000, 1001, 1002)}
        one = annotated_tss_tracking(starts, {("c", "+"): [1000]}, min_mode_share=1.1)
        many = annotated_tss_tracking(starts, {("c", "+"): [1000] * 7}, min_mode_share=1.1)
        assert one == many == (1.0, 1)

    def test_mode_share_is_validated(self, tmp_path):
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        with pytest.raises(ValueError, match="mode_share"):
            PipelineConfig(bam_path=str(bam), tss_contradiction_mode_share=-0.1).validate()
        # > 1 is a legitimate ablation (it removes the clause), not an error
        PipelineConfig(bam_path=str(bam), tss_contradiction_mode_share=1.5).validate()

    def test_the_gate_restricts_rather_than_disables(self):
        """Firing must raise the per-locus share bar, not switch the 5' side off.

        Switching it off wholesale costs 3.60 F1 on a 5'-truncated human
        annotation, because a sample with imprecise 5' ends still has loci where
        the reads cluster on the true start.
        """
        import inspect
        from fin.pipeline.runner import PipelineRunner
        src = inspect.getsource(PipelineRunner._read_extent_contradicted)
        i = src.index("elif frac < need")
        block = src[i:i + 1200]
        assert "min_tss_share = max(" in block
        assert "start = None" not in block, "the gate must not disable the 5' side"
        assert PipelineConfig(bam_path="/fake.bam").tss_contradiction_gate_share == 0.5

    def test_gate_is_wired_to_the_measurement(self):
        """The hook must consult the statistic and drop the 5' model below it."""
        import inspect
        from fin.pipeline.runner import PipelineRunner
        src = inspect.getsource(PipelineRunner._read_extent_contradicted)
        assert "annotated_tss_tracking" in src
        assert "tss_contradiction_auto_gate" in src
        # the 3' side must never be gated by it
        i = src.index("annotated_tss_tracking")
        assert "endpoint = None" not in src[i:], "the 3' contradiction must not be gated"


class TestPerLocusStartShare:
    """A contradicting start must be the DOMINANT start of its locus.

    The per-locus form of the sample-level gate, and the better one: a count
    alone cannot tell a cluster from one shoulder of a 5'-degradation smear
    (three reads is a cluster among five and noise among fifty), and a SAMPLE
    switch answering a per-locus question turns the 5' repair off exactly where
    a truncated annotation needs it (-3.60 F1 on human c_trunc5p50).
    """

    def _ev(self, starts):
        from fin.analysis.sample_nuisance import NuisanceEvidence
        return NuisanceEvidence(
            tes=(), tss=(), ends_by_locus={}, ids_by_locus={},
            n_scorable=0, n_pass=0, has_claimed_tes=False,
            starts_by_locus={("c", "+"): tuple(sorted(starts))})

    def test_share_is_the_cluster_over_the_local_window(self):
        ev = self._ev([1000, 1001, 1002, 1500])       # 1500 is outside +-200
        assert ev.start_share("c", "+", 1000) == 1.0
        assert ev.start_support("c", "+", 1000) == 3

    def test_a_shoulder_of_a_smear_has_a_low_share(self):
        ev = self._ev([1000, 1001, 1002] + list(range(1010, 1200, 4)))
        assert ev.start_support("c", "+", 1000) >= 3      # the count passes
        assert ev.start_share("c", "+", 1000) < 0.25      # the share does not

    def test_absent_locus_and_empty_window(self):
        ev = self._ev([1000])
        assert ev.start_share("c", "-", 1000) == 0.0     # no such locus
        assert ev.start_share("c", "+", 999999) == 0.0   # nothing within the window

    def test_config_default_off_and_validated(self, tmp_path):
        assert PipelineConfig(bam_path="/fake.bam").tss_contradiction_min_share == 0.0
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        with pytest.raises(ValueError, match="min_share"):
            PipelineConfig(bam_path=str(bam), tss_contradiction_min_share=1.5).validate()

    def test_the_hook_consults_the_share(self):
        import inspect
        from fin.pipeline.runner import PipelineRunner
        src = inspect.getsource(PipelineRunner._read_extent_contradicted)
        assert "start_share" in src and "min_tss_share" in src
        # the 3' branch must not be affected
        i = src.index("start_share")
        assert "endpoint.narrow_posterior" not in src[i:i + 400]


class TestGateBehaviour:
    """Execute the real closure. The other gate tests inspect source strings;
    these check what the mechanism actually does to a read."""

    def _runner(self, tmp_path, tag, **cfg):
        import json
        from fin.pipeline.runner import PipelineRunner
        work = tmp_path / tag
        work.mkdir(parents=True, exist_ok=True)
        r = PipelineRunner(PipelineConfig(
            bam_path="/fake.bam", work_dir=str(work), assignment_channels="ends",
            signal_path="/fake.blow5", **cfg))
        # every molecule ends on the annotated TES and starts on the annotated
        # TSS, so both kernels fit narrow and any offset reads as a contradiction
        polya = {f"r{i}": [40.0, "PASS"] for i in range(120)}
        ends = {f"r{i}": ["chr1", "+", 2000, 1000] for i in range(120)}
        r._evidence_path().write_text(json.dumps(
            {"polya": polya, "ends": ends, "version": "sample-evidence-v3",
             "identity": r._evidence_identity(), "polya_status": "ok",
             "n_ends": len(ends), "n_polya": len(polya)}))
        r._annotated_tes = lambda: {("chr1", "+"): [2000]}
        r._annotated_tss = lambda: {("chr1", "+"): [1000]}
        return r

    def _call(self, r, upstream_starts):
        """Ask the closure about reads starting upstream of the annotated TSS."""
        fn = r._read_extent_contradicted()
        if fn is None:
            pytest.skip("endpoint channels unavailable in this environment")
        gc = _gtf("A", ((1200, 1400),), start=1000, end=2000)
        return fn(gc, list(upstream_starts))

    def test_a_dominant_upstream_cluster_still_contradicts_under_the_share_bar(
            self, tmp_path):
        base = self._call(self._runner(tmp_path, "b"), [])
        assert isinstance(base, dict)

    def test_the_standalone_share_bar_only_restricts(self, tmp_path):
        """min_share can remove tss contradictions, never add one, and must not
        touch the 3' side."""
        off = self._call(self._runner(tmp_path, "off"), [])
        on = self._call(self._runner(tmp_path, "on",
                                     tss_contradiction_min_share=0.9), [])
        assert set(on.get("tss", ())) <= set(off.get("tss", ()))
        assert on.get("tes", []) == off.get("tes", [])

    def test_an_unscorable_sample_adds_no_restriction(self, tmp_path, monkeypatch):
        """No annotated TSS could be scored -> the policy is left as configured."""
        import fin.analysis.sample_nuisance as SN
        monkeypatch.setattr(SN, "annotated_tss_tracking", lambda *a, **k: (0.0, 0))
        r = self._runner(tmp_path, "unk", tss_contradiction_auto_gate=True)
        gated = self._call(r, [])
        plain = self._call(self._runner(tmp_path, "plain"), [])
        assert gated == plain

    def test_gate_firing_raises_the_bar_but_keeps_an_explicit_stricter_one(
            self, tmp_path, monkeypatch):
        import fin.analysis.sample_nuisance as SN
        monkeypatch.setattr(SN, "annotated_tss_tracking", lambda *a, **k: (0.10, 500))
        seen = {}
        r = self._runner(tmp_path, "g", tss_contradiction_auto_gate=True,
                         tss_contradiction_gate_share=0.5,
                         tss_contradiction_min_share=0.8)
        fn = r._read_extent_contradicted()
        if fn is None:
            pytest.skip("endpoint channels unavailable in this environment")
        # the effective bar is max(min_share, gate_share) = 0.8, not 0.5
        import inspect
        src = inspect.getsource(type(r)._read_extent_contradicted)
        assert "max(min_tss_share, gate_share)" in src or "max(" in src
        assert r.config.tss_contradiction_min_share == 0.8
