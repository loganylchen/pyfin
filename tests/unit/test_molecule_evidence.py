"""The per-molecule evidence row: one row per molecule, absent channels abstain,
resampling keeps a molecule whole."""
from __future__ import annotations

import numpy as np

from fin.analysis.molecule_evidence import (
    POLYA_FAIL, POLYA_MISSING, POLYA_PASS, MoleculeEvidence,
    build_molecule_evidence, polya_state, resample_molecules,
)
from fin.candidates.dataclasses import CandidateSet, IntronChain, TranscriptCandidate
from fin.io.interval_manager import GenomicInterval


def _cand(cid):
    return TranscriptCandidate(
        candidate_id=cid, intron_chain=IntronChain(introns=((100, 200),)),
        three_prime_pos=500, sequence="ACGT" * 20, source="novel",
        supporting_read_ids=set(), chrom="chr1", strand="+", start=0, end=500)


class _QO:
    def __init__(self, cands, kept, raw, ties, nlls, cover):
        self.cand_list, self.kept_read_ids, self.raw = cands, kept, raw
        self.ties_by_read, self.nlls_by_read, self.cover_by_read = ties, nlls, cover


def _cs(spans, strands=None, junctions=None, indels=None):
    return CandidateSet(
        interval=GenomicInterval(chrom="chr1", start=0, end=1000, strand="+"),
        candidates=[_cand("a"), _cand("b")], read_ids=set(spans),
        read_spans=spans, read_strands=strands or {}, read_junctions=junctions or {},
        indels_by_read=indels or {})


class TestPolyaState:
    def test_pass_needs_positive_length(self):
        assert polya_state((30.0, "PASS")) == (POLYA_PASS, 30.0)
        assert polya_state((0.0, "PASS"))[0] == POLYA_FAIL
        assert polya_state((None, "PASS"))[0] == POLYA_FAIL

    def test_fail_and_missing_are_distinct_but_both_abstain(self):
        assert polya_state((12.0, "FAIL"))[0] == POLYA_FAIL
        assert polya_state(None)[0] == POLYA_MISSING
        assert polya_state((None, None))[0] == POLYA_MISSING


class TestBuild:
    def _setup(self):
        cands = [_cand("a"), _cand("b")]
        kept = ["r1", "r2"]
        raw = np.array([[10.0, 10.0], [9.0, -np.inf]])
        ties = {0: [0, 1], 1: [0]}
        nlls = {0: {0: 1.5, 1: 2.5}}
        cover = {0: True}
        cs = _cs({"r1": (10, 400), "r2": (50, 420)},
                 strands={"r1": "+", "r2": "-"},
                 junctions={"r1": ((100, 200),)}, indels={"r1": ((150, 152),)})
        return _QO(cands, kept, raw, ties, nlls, cover), cs

    def test_one_row_per_molecule_with_all_channels(self):
        qo, cs = self._setup()
        rows = build_molecule_evidence(qo, cs, polya_map={"r1": (40.0, "PASS")})
        assert set(rows) == {"r1", "r2"}
        r1 = rows["r1"]
        assert r1.raw_tie_ids == ("a", "b") and r1.alignment_scores == {"a": 10.0, "b": 10.0}
        assert r1.m2_nll == {"a": 1.5, "b": 2.5}
        assert r1.covered_diff_regions == {"a": ("diff",), "b": ("diff",)}
        assert r1.cigar_junctions == ((100, 200),) and r1.indel_spans == ((150, 152),)
        assert r1.polya_qc == POLYA_PASS and r1.polya_available
        assert (r1.five_prime_pos, r1.three_prime_pos) == (10, 400)   # + strand

    def test_minus_strand_swaps_endpoints(self):
        qo, cs = self._setup()
        rows = build_molecule_evidence(qo, cs)
        r2 = rows["r2"]
        assert r2.strand == "-"
        assert (r2.five_prime_pos, r2.three_prime_pos) == (420, 50)

    def test_absent_channels_abstain(self):
        qo, cs = self._setup()
        rows = build_molecule_evidence(qo, cs)          # no polya_map at all
        r2 = rows["r2"]
        assert r2.raw_tie_ids == ("a",) and not r2.structural_available("b")
        assert r2.m2_nll == {} and not r2.signal_available("a")
        assert r2.covered_diff_regions == {}
        assert r2.polya_qc == POLYA_MISSING and not r2.polya_available
        assert "b" not in r2.alignment_scores                  # -inf never stored

    def test_fail_call_abstains_like_missing(self):
        qo, cs = self._setup()
        rows = build_molecule_evidence(qo, cs, polya_map={"r1": (5.0, "FAIL")})
        assert rows["r1"].polya_qc == POLYA_FAIL and not rows["r1"].polya_available

    def test_reads_without_span_are_not_rows(self):
        qo, cs = self._setup()
        del cs.read_spans["r2"]
        assert set(build_molecule_evidence(qo, cs)) == {"r1"}

    def test_rows_are_immutable(self):
        qo, cs = self._setup()
        r = build_molecule_evidence(qo, cs)["r1"]
        try:
            r.polya_qc = "PASS"
        except Exception:
            return
        raise AssertionError("row was mutable")


class TestResample:
    def test_whole_row_copies_with_stable_ids(self):
        qo, cs = TestBuild()._setup()
        rows = list(build_molecule_evidence(qo, cs, polya_map={"r1": (40.0, "PASS")}).values())
        boot = resample_molecules(rows, [0, 0, 1], replicate=3)
        assert [b.read_id for b in boot] == ["r1#boot3#0", "r1#boot3#1", "r2#boot3#0"]
        # every channel travels with the molecule
        assert boot[0].m2_nll == rows[0].m2_nll and boot[0].polya_qc == rows[0].polya_qc
        assert boot[1].three_prime_pos == rows[0].three_prime_pos
        # deterministic
        assert [b.read_id for b in resample_molecules(rows, [0, 0, 1], replicate=3)] == \
               [b.read_id for b in boot]
