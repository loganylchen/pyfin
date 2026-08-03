"""Unit tests for family-scoped read->candidate alignment goodness.

Synthetic transcripts + reads exercise the two structural signatures the goodness must capture:
  - 5' containment: a read reaching the long candidate's unique 5' exon maps better to it (higher AS,
    extra 5' soft-clip on the shorter sibling) -> identifies the long one; a read starting in the
    shared suffix maps equally to both -> ambiguous.
  - internal exon difference: a read that would need a >50bp internal indel to fit a candidate is
    structurally INCOMPATIBLE with it (score_hit rejects -> no hit).
"""
import random

from fin.scoring.candidate_align import align_reads_to_candidates, is_compatible


class _Cand:
    def __init__(self, cid, seq):
        self.candidate_id = cid
        self.sequence = seq


def _rand_seq(n, seed):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


BASE = _rand_seq(1500, 42)          # "candA" full transcript
SHORT = BASE[200:]                   # "candB": 5'-shorter (missing the first 200 bp exon)


class TestContainment:
    def _run(self, reads):
        cands = [_Cand("A", BASE), _Cand("B", SHORT)]
        return align_reads_to_candidates(cands, reads, best_n=50)

    def test_full_read_identifies_long(self):
        # a read spanning the full long transcript reaches A's unique 5' region.
        out = self._run({"full": BASE})
        g = out["full"]
        assert "A" in g and "B" in g
        assert g["A"].AS > g["B"].AS + 100          # clearly better on the long one
        assert g["B"].s5 >= 180                       # ~200 bp of unique 5' clips on the short one
        assert g["A"].s5 <= 5
        assert g["A"].q_cov > g["B"].q_cov            # more of the read explained by A

    def test_truncated_read_is_ambiguous(self):
        # a 5'-degraded read starting in the shared suffix cannot tell A from B.
        out = self._run({"trunc": BASE[200:]})
        g = out["trunc"]
        assert "A" in g and "B" in g
        assert abs(g["A"].AS - g["B"].AS) < 60        # ~equal -> ambiguous
        assert g["A"].q_cov > 0.95 and g["B"].q_cov > 0.95
        assert is_compatible(g["A"]) and is_compatible(g["B"])


class TestInternalExonDifference:
    def test_read_with_extra_internal_exon_incompatible_with_skip(self):
        # candSkip lacks a 60bp internal block that the read includes -> a >50bp internal INSERTION
        # against candSkip -> score_hit rejects -> the read is not compatible with candSkip.
        skip = BASE[:700] + BASE[760:]                # removed a 60bp internal exon
        cands = [_Cand("full", BASE), _Cand("skip", skip)]
        out = align_reads_to_candidates(cands, {"r": BASE}, best_n=50)
        g = out["r"]
        assert "full" in g                            # clean match to the including candidate
        assert "skip" not in g                        # rejected (exon-sized internal indel)


class TestGoodnessFields:
    def test_perfect_match_fields(self):
        out = align_reads_to_candidates([_Cand("A", BASE)], {"r": BASE}, best_n=5)
        g = out["r"]["A"]
        assert g.q_cov > 0.99 and g.event_id > 0.99 and g.s5 == 0 and g.cand3p <= 5
        assert is_compatible(g)

    def test_empty_inputs(self):
        assert align_reads_to_candidates([], {"r": "ACGT"}) == {}
        assert align_reads_to_candidates([_Cand("A", BASE)], {}) == {}
        assert align_reads_to_candidates([_Cand("A", "")], {"r": BASE}) == {}

    def test_candidate_id_with_whitespace_is_safe(self):
        # mappy truncates FASTA headers at whitespace; synthetic contig names must keep the mapping
        # intact so the id round-trips and cand3p stays non-negative.
        cid = "gene A|tx 1"
        out = align_reads_to_candidates([_Cand(cid, BASE)], {"r": BASE}, best_n=5)
        assert cid in out["r"]
        assert out["r"][cid].cand3p >= 0 and is_compatible(out["r"][cid])

    def test_duplicate_candidate_ids_deduped(self):
        # two candidates sharing an id (different sequences) must not collapse into a corrupt entry.
        out = align_reads_to_candidates(
            [_Cand("dup", BASE), _Cand("dup", BASE[300:])], {"r": BASE}, best_n=5)
        g = out["r"]["dup"]
        assert g.cand3p >= 0            # keyed to the first (full) sequence, not corrupted
        assert g.strand in (1, -1)
