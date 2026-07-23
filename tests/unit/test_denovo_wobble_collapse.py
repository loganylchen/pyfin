"""Unit tests for the de-novo wobble-tolerant collapse (EXPERIMENT, default-off).

When ``wobble_tol > 0``, ``_collapse_candidates`` merges NOVEL candidates whose
intron chains match within ``tol`` bp per junction (and 3' within threshold) into
the highest-read-support consensus, BEFORE they become separate candidates. This
attacks wobble shadows (a junction and its ±few-bp minimap variant surviving as
separate novel transcripts). ``wobble_tol == 0`` is byte-identical to the exact
path.
"""
from __future__ import annotations

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.candidates.discovery import (
    _collapse_candidates,
    _wobble_chain_match,
    _wobble_collapse_novel,
)


def _cand(introns, *, tp=None, reads=None, source="novel", start=1, cid="c"):
    end = (introns[-1][1] + 100) if introns else 200
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=tuple(introns)),
        three_prime_pos=tp if tp is not None else end,
        sequence="ACGT" * 50,
        source=source,
        supporting_read_ids=set(reads or set()),
        chrom="chr1",
        strand="+",
        start=start,
        end=end,
    )


class TestWobbleMatch:
    def test_identical(self):
        a = IntronChain(((100, 200),)); b = IntronChain(((100, 200),))
        assert _wobble_chain_match(a, b, 6)

    def test_within_tol(self):
        a = IntronChain(((100, 200),)); b = IntronChain(((103, 197),))
        assert _wobble_chain_match(a, b, 6)      # both boundaries within 6
        assert not _wobble_chain_match(a, b, 2)  # 3 bp > tol=2

    def test_beyond_tol(self):
        a = IntronChain(((100, 200),)); b = IntronChain(((110, 190),))
        assert not _wobble_chain_match(a, b, 6)  # 10 bp shift

    def test_different_intron_count(self):
        a = IntronChain(((100, 200),))
        b = IntronChain(((100, 200), (400, 500)))
        assert not _wobble_chain_match(a, b, 6)

    def test_mono_exon_match(self):
        a = IntronChain(()); b = IntronChain(())
        assert _wobble_chain_match(a, b, 6)  # no junctions -> trivially match


class TestWobbleCollapse:
    def test_wobbled_pair_merges_to_consensus(self):
        # c1 = 3 reads (consensus), c2 = 1 read, junctions 3bp apart, 3' 5bp apart
        c1 = _cand([(100, 200)], tp=300, reads={"r1", "r2", "r3"}, cid="c1")
        c2 = _cand([(103, 197)], tp=305, reads={"r4"}, cid="c2")
        out = _collapse_candidates([c1, c2], threshold=24, wobble_tol=6)
        assert len(out) == 1
        rep = out[0]
        assert rep.candidate_id == "c1"                 # higher support = consensus rep
        assert rep.intron_chain.introns == ((100, 200),)  # consensus keeps its coords
        assert rep.supporting_read_ids == {"r1", "r2", "r3", "r4"}  # reads unioned

    def test_off_keeps_wobbles_separate(self):
        # wobble_tol=0 -> exact-chain path -> the two wobbled chains stay separate
        c1 = _cand([(100, 200)], tp=300, reads={"r1"}, cid="c1")
        c2 = _cand([(103, 197)], tp=300, reads={"r2"}, cid="c2")
        out = _collapse_candidates([c1, c2], threshold=24, wobble_tol=0)
        assert len(out) == 2

    def test_off_exact_chain_merges_to_longest(self):
        # tol=0 exact path: SAME chain + 3' within threshold -> merge; rep is the
        # longest-span candidate (existing policy); reads unioned. Different chain
        # -> separate. Guards the byte-identical disabled path.
        long_c = _cand([(100, 200)], tp=300, reads={"r1"}, start=1, cid="long")
        short_c = _cand([(100, 200)], tp=305, reads={"r2"}, start=150, cid="short")
        diff = _cand([(100, 300)], tp=400, reads={"r3"}, cid="diff")  # different chain
        out = _collapse_candidates([short_c, long_c, diff], threshold=24, wobble_tol=0)
        ids = {c.candidate_id for c in out}
        assert ids == {"long", "diff"}                          # short merged into long
        rep = next(c for c in out if c.candidate_id == "long")
        assert rep.supporting_read_ids == {"r1", "r2"}

    def test_shadow_ratio_protects_real_close_isoforms(self):
        # Two wobble-matching candidates with COMPARABLE support are genuine close
        # isoforms, not shadows: at ratio=0.5 they must stay SEPARATE.
        c1 = _cand([(100, 200)], tp=300, reads={"a", "b", "c"}, cid="c1")
        c2 = _cand([(103, 197)], tp=300, reads={"d", "e", "f"}, cid="c2")  # equal support
        out = _wobble_collapse_novel([c1, c2], 24, 6, shadow_ratio=0.5)
        assert len(out) == 2  # 3 <= 0.5*3 is False -> not merged

    def test_shadow_ratio_merges_weak_shadow(self):
        # A weak shadow (1 read) vs a strong rep (4 reads): 1 <= 0.5*4=2 -> merge.
        rep = _cand([(100, 200)], tp=300, reads={"a", "b", "c", "d"}, cid="rep")
        shadow = _cand([(103, 197)], tp=300, reads={"x"}, cid="shadow")
        out = _wobble_collapse_novel([rep, shadow], 24, 6, shadow_ratio=0.5)
        assert len(out) == 1
        assert out[0].candidate_id == "rep"
        assert out[0].supporting_read_ids == {"a", "b", "c", "d", "x"}

    def test_beyond_tol_stays_separate(self):
        c1 = _cand([(100, 200)], tp=300, reads={"r1", "r2"}, cid="c1")
        c2 = _cand([(112, 188)], tp=300, reads={"r3"}, cid="c2")  # 12bp shift
        out = _collapse_candidates([c1, c2], threshold=24, wobble_tol=6)
        assert len(out) == 2

    def test_different_three_prime_stays_separate(self):
        # wobble-match but 3' ends far apart (> threshold) -> different transcripts
        c1 = _cand([(100, 200)], tp=300, reads={"r1"}, cid="c1")
        c2 = _cand([(101, 199)], tp=1000, reads={"r2"}, cid="c2")
        out = _collapse_candidates([c1, c2], threshold=24, wobble_tol=6)
        assert len(out) == 2

    def test_wobble_on_is_superset_of_exact(self):
        # Only exact-chain candidates (no inexact wobbles): wobble-on must yield
        # the SAME candidate set as wobble-off — the exact collapse runs first, so
        # its 3'-chaining is preserved; wobble adds nothing. Guards Codex's
        # 0/20/40 non-transitivity counterexample.
        def build():
            return [_cand([(100, 200)], tp=t, reads={f"r{t}"}, cid=f"c{t}")
                    for t in (0, 20, 40)]
        off = _collapse_candidates(build(), 24, wobble_tol=0)
        on = _collapse_candidates(build(), 24, wobble_tol=6, wobble_shadow_ratio=0.5)
        assert {c.candidate_id for c in off} == {c.candidate_id for c in on}

    def test_gtf_and_fusion_passthrough(self):
        g = _cand([(100, 200)], source="gtf", cid="ENST1")
        f = _cand([(100, 200)], source="fusion", cid="fus")
        n = _cand([(103, 197)], reads={"r1"}, cid="c1")
        out = _collapse_candidates([g, f, n], threshold=24, wobble_tol=6)
        ids = {c.candidate_id for c in out}
        assert "ENST1" in ids and "fus" in ids  # not collapsed

    def test_exact_chain_merges_under_wobble_even_comparable_support(self):
        # Codex correctness gap: EXACT same chain + 3' within threshold must ALWAYS
        # merge with wobble on (superset of the disabled exact path) — the
        # shadow_ratio guard applies only to INEXACT (wobbled) matches.
        c1 = _cand([(100, 200)], tp=300, reads={"a", "b", "c"}, cid="c1")
        c2 = _cand([(100, 200)], tp=305, reads={"d", "e", "f"}, cid="c2")  # same chain
        out = _wobble_collapse_novel([c1, c2], 24, 6, shadow_ratio=0.5)
        assert len(out) == 1
        assert out[0].supporting_read_ids == {"a", "b", "c", "d", "e", "f"}

    def test_determinism(self):
        # Fresh objects per call (the merge mutates supporting_read_ids in place).
        def build():
            return [_cand([(100 + i, 200 - i)], tp=300, reads={f"r{i}"}, cid=f"c{i}")
                    for i in range(5)]  # all within 6bp of c0
        a = _wobble_collapse_novel(build(), 24, 6, shadow_ratio=1.0)
        b = _wobble_collapse_novel(list(reversed(build())), 24, 6, shadow_ratio=1.0)
        assert [c.candidate_id for c in a] == [c.candidate_id for c in b]

    def test_determinism_on_full_key_ties(self):
        # Two candidates equal on support/length/start/chain but different 3' + id:
        # total order (candidate_id last) makes the rep deterministic regardless of
        # input order. Fresh objects each call (merge mutates in place).
        def build():
            return [_cand([(100, 200)], tp=302, reads={"r1"}, start=1, cid="x"),
                    _cand([(100, 200)], tp=300, reads={"r2"}, start=1, cid="y")]
        r1 = _wobble_collapse_novel(build(), 24, 6)
        r2 = _wobble_collapse_novel(list(reversed(build())), 24, 6)
        assert len(r1) == 1 and len(r2) == 1
        assert r1[0].candidate_id == r2[0].candidate_id == "y"  # lower 3' wins
