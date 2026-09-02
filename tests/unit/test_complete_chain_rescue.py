"""P2 complete splice-chain evidence rescue (complete_chain_rescue_reads).

Contract:
1. Default 0 disables everything: no field changes any gate decision.
2. weakest_clean_junction_depth: the candidate metric is the MINIMUM, over
   its junctions, of the number of exact complete-chain reads whose
   alignment has no indel within 10 bp of that junction; folded (degraded
   sub-chain) reads never count; mono chains score 0.
3. apply_complete_chain_rescue removes sufficiently supported candidates
   (novel or gtf: a chain-correct endpoint-truncated annotation entry
   qualifies; junction-corrupt chains never carry a nonzero metric).
4. Soft-mass global filter honors the rescue for QuantResults.
"""

from __future__ import annotations

from fin.analysis.quantification import QuantResult
from fin.candidates.discovery import junction_clean_depths
from fin.pipeline.config import PIPELINE_PROFILES, PipelineConfig
from fin.pipeline.runner import PipelineRunner
from fin.pipeline.selection import (
    COMPLETE_CHAIN_RESCUE_NOVEL_MONO_JUNC_MIN,
    apply_complete_chain_rescue,
    complete_chain_rescue_ok,
)


CHAIN = ((100, 200), (300, 400))  # two junctions: 100-200 and 300-400


class TestJunctionCleanDepths:
    def test_all_clean(self):
        ind = {"a": (), "b": ()}
        assert junction_clean_depths(CHAIN, {"a", "b"}, set(), ind) == (2, 2)

    def test_folded_reads_excluded(self):
        ind = {"a": (), "b": ()}
        assert junction_clean_depths(CHAIN, {"a", "b"}, {"b"}, ind) == (1, 1)

    def test_weakest_junction_governs(self):
        # read a: clean everywhere; read b: deletion at 393-396 (within 10
        # of acceptor 400) -> junction 2 depth 1, junction 1 depth 2
        ind = {"a": (), "b": ((393, 396),)}
        assert junction_clean_depths(CHAIN, {"a", "b"}, set(), ind) == (2, 1)

    def test_spurious_junction_signature_zero(self):
        # every exact read carries an indel at the first junction
        ind = {"a": ((95, 96),), "b": ((104, 105),)}
        assert junction_clean_depths(CHAIN, {"a", "b"}, set(), ind) == (0, 2)

    def test_far_indels_ignored(self):
        ind = {"a": ((250, 252),)}  # 50 bp from both junctions
        assert junction_clean_depths(CHAIN, {"a"}, set(), ind) == (1, 1)

    def test_mono_chain_empty(self):
        assert junction_clean_depths((), {"a"}, set(), {"a": ()}) == ()

    def test_no_exact_reads_zeros(self):
        assert junction_clean_depths(CHAIN, {"a"}, {"a"}, {"a": ()}) == (0, 0)

    def test_mapq_filtered_read_never_counts(self):
        # a read absent from indels_by_read (collection-filtered) is ignored
        assert junction_clean_depths(CHAIN, {"x"}, set(), {}) == (0, 0)


class _Chain:
    def __init__(self, n):
        self.introns = tuple((i * 100, i * 100 + 50) for i in range(n))


class _Cand:
    def __init__(self, n, source="novel", n_junc=3):
        self.complete_chain_weakest_indel_clean_depth = n
        self.source = source
        self.intron_chain = _Chain(n_junc)


class TestRescuePolicy:
    def test_gtf_bar_is_k(self):
        assert complete_chain_rescue_ok("gtf", 1, 2, 2) is True
        assert complete_chain_rescue_ok("gtf", 5, 2, 2) is True
        assert complete_chain_rescue_ok("gtf", 5, 1, 2) is False

    def test_novel_multi_bar_is_k_plus_one(self):
        assert complete_chain_rescue_ok("novel", 2, 3, 2) is True
        assert complete_chain_rescue_ok("novel", 2, 2, 2) is False

    def test_novel_single_junction_artifact_bar(self):
        m = COMPLETE_CHAIN_RESCUE_NOVEL_MONO_JUNC_MIN
        assert complete_chain_rescue_ok("novel", 1, m, 2) is True
        assert complete_chain_rescue_ok("novel", 1, m - 1, 2) is False

    def test_disabled_and_exempt_sources(self):
        assert complete_chain_rescue_ok("gtf", 3, 50, 0) is False
        assert complete_chain_rescue_ok("fusion", 3, 50, 2) is False
        assert complete_chain_rescue_ok("novel", 0, 50, 2) is False


class TestApplyRescue:
    def test_disabled_returns_same(self):
        d = {0, 1}
        out, n = apply_complete_chain_rescue(d, [_Cand(5), _Cand(5)], 0)
        assert out is d and n == 0

    def test_rescues_only_supported(self):
        cands = [_Cand(0), _Cand(2), _Cand(3), _Cand(7)]
        out, n = apply_complete_chain_rescue({0, 1, 2, 3}, cands, 2)
        assert out == {0, 1} and n == 2  # novel multi bar is K+1=3

    def test_gtf_rescued_when_chain_attested(self):
        """A chain-correct (e.g. endpoint-truncated) GTF entry with every
        junction cleanly attested must be rescued; an unattested one not."""
        cands = [_Cand(2, source="gtf"), _Cand(0, source="gtf")]
        out, n = apply_complete_chain_rescue({0, 1}, cands, 2)
        assert out == {1} and n == 1

    def test_novel_single_junction_needs_artifact_bar(self):
        cands = [
            _Cand(8, n_junc=1),   # observed artifact depth -> not rescued
            _Cand(COMPLETE_CHAIN_RESCUE_NOVEL_MONO_JUNC_MIN, n_junc=1),
        ]
        out, n = apply_complete_chain_rescue({0, 1}, cands, 2)
        assert out == {0} and n == 1

    def test_empty_dropset(self):
        out, n = apply_complete_chain_rescue(set(), [_Cand(9)], 2)
        assert out == set() and n == 0


def _qr(cid, source, abundance, hard, nccr=0, nex=3):
    return QuantResult(
        candidate_id=cid, abundance=abundance, confidence=0.0,
        num_assigned_reads=hard, source=source,
        exons=tuple((i * 20, i * 20 + 10) for i in range(nex)),
        complete_chain_weakest_indel_clean_depth=nccr,
    )


def _runner(**overrides) -> PipelineRunner:
    cfg = PipelineConfig(
        bam_path="/tmp/x.bam",
        min_abundance=0.0, min_gtf_abundance=0.0,
        min_isoform_fraction=0.0, min_fulllen_fraction=0.0,
        min_polya5p_reads=0, max_soft_mass_ratio=2.0,
        **overrides,
    )
    runner = PipelineRunner.__new__(PipelineRunner)
    runner.config = cfg
    runner._gtf_reader = None
    return runner


class TestSoftMassRescue:
    def test_rescue_keeps_supported_novel(self):
        runner = _runner(complete_chain_rescue_reads=2)
        agg = {
            # soft 10 vs hard 2 -> ratio 5 >= 2.0 would drop; rescued
            # (novel multi-junction bar K+1=3, nccr 3)
            "n_true": _qr("n_true", "novel", 10.0, 2, nccr=3),
            # same ratio, below the bar -> still dropped
            "n_junk": _qr("n_junk", "novel", 10.0, 2, nccr=2),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert "n_true" in out
        assert "n_junk" not in out

    def test_gtf_source_not_rescued_by_flag(self):
        """The soft-mass filter itself exempts gtf, so nothing changes for
        gtf rows; the rescue must not alter that behavior either."""
        runner = _runner(complete_chain_rescue_reads=2)
        agg = {"g": _qr("g", "gtf", 10.0, 2, nccr=9)}
        out = runner._finalize_and_write(agg, None, None)
        assert "g" in out  # gtf exempt from soft-mass regardless

    def test_default_off_drops_both(self):
        runner = _runner()
        agg = {
            "n_true": _qr("n_true", "novel", 10.0, 2, nccr=3),
            "n_junk": _qr("n_junk", "novel", 10.0, 2, nccr=2),
        }
        out = runner._finalize_and_write(agg, None, None)
        assert "n_true" not in out and "n_junk" not in out


class TestDefaults:
    def test_config_default_zero(self):
        assert PipelineConfig(bam_path="/t.bam").complete_chain_rescue_reads == 0

    def test_cli_default_zero(self):
        from fin.cli import main

        opt = next(
            p for p in main.params if p.name == "complete_chain_rescue_reads"
        )
        assert opt.default == 0

    def test_no_profile_enables_it(self):
        for name, values in PIPELINE_PROFILES.items():
            assert "complete_chain_rescue_reads" not in values, name

    def test_quantresult_field_default_zero(self):
        qr = QuantResult(
            candidate_id="x", abundance=0.0, confidence=0.0,
            num_assigned_reads=0, source="novel",
        )
        assert qr.complete_chain_weakest_indel_clean_depth == 0
