"""Tests for the M1-keep SPLIT production quant (``_quant_argmax_keep``).

argmax_keep assigns each read to ALL of its simultaneously-best-AS candidates,
splitting the per-read abundance 1/K across the K tied winners (read mass
conserved). This mirrors the ablation harness ``_quant_tie(mode="split")`` ==
config "M1-split".
"""

from __future__ import annotations

import mappy
import pytest

from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval
from fin.pipeline.config import PipelineConfig
from fin.pipeline.runner import PipelineRunner

# A unique read region and two distractor candidates. The read sequence below
# aligns perfectly to candidates carrying it verbatim.
READ = (
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "TTGGCCAATTGGCCAATTGGCCAATTGGCCAATTGGCCAA"
)
# A clearly different sequence so a read never aligns to it.
OTHER = (
    "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"
    "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
)


def _cand(cid, seq, source="novel"):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=()),
        three_prime_pos=0,
        sequence=seq,
        source=source,
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=0,
        end=len(seq),
    )


def _runner():
    cfg = PipelineConfig(bam_path="x", r1_variant="argmax_keep")
    return PipelineRunner(cfg)


def _interval():
    return GenomicInterval(chrom="chr1", start=0, end=100, strand="+")


def _quant(cands, read_seqs):
    cs = CandidateSet(
        interval=_interval(),
        candidates=cands,
        read_ids=set(read_seqs),
        read_sequences=read_seqs,
    )
    runner = _runner()
    results = runner._quant_argmax_keep(cs, list(read_seqs), _interval())
    return {r.candidate_id: r for r in results}


@pytest.mark.skipif(
    mappy is None, reason="mappy required for argmax_keep quant"
)
class TestArgmaxKeepSplit:
    def test_unique_winner_gets_full_mass(self):
        cands = [_cand("c_read", READ), _cand("c_other", OTHER)]
        out = _quant(cands, {"r1": READ})
        assert out["c_read"].abundance == pytest.approx(1.0)
        assert out["c_other"].abundance == pytest.approx(0.0)
        assert out["c_read"].max_R == pytest.approx(1.0)

    def test_two_way_tie_splits_half_half(self):
        # Two identical candidate sequences -> a read ties on both -> 0.5 each.
        cands = [_cand("c_a", READ), _cand("c_b", READ)]
        out = _quant(cands, {"r1": READ})
        assert out["c_a"].abundance == pytest.approx(0.5)
        assert out["c_b"].abundance == pytest.approx(0.5)
        # Read mass is conserved across the tie.
        assert out["c_a"].abundance + out["c_b"].abundance == pytest.approx(1.0)
        assert out["c_a"].max_R == pytest.approx(0.5)
        assert "r1" in out["c_a"].assigned_read_ids
        assert "r1" in out["c_b"].assigned_read_ids

    def test_mass_conserved_over_many_reads(self):
        cands = [_cand("c_a", READ), _cand("c_b", READ), _cand("c_other", OTHER)]
        reads = {f"r{i}": READ for i in range(4)}
        out = _quant(cands, reads)
        total = sum(out[c].abundance for c in out)
        assert total == pytest.approx(4.0)
        assert out["c_other"].abundance == pytest.approx(0.0)

    def test_no_hit_read_contributes_nothing(self):
        cands = [_cand("c_read", READ)]
        out = _quant(cands, {"r1": OTHER})
        assert out["c_read"].abundance == pytest.approx(0.0)
        assert out["c_read"].max_R == pytest.approx(0.0)
