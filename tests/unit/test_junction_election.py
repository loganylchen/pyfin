"""P3 K-way signal-ballot junction election (junction_election_mode).

Contract:
1. Nomination: accused junction (indel-clean <= 1) with >= 1 competitor
   within +-15 bp holding >= max(3, 3x) clean support; accused variants
   sharing a competitor merge into ONE K-way contest; carriers of a
   champion junction are never accused; every nominee needs an emitted
   carrier; contests beyond the cap are reported as truncated.
2. Ballots (injected scorer): a valid ballot needs >= 2 scored hypotheses,
   >= JE_BALLOT_MIN_EVENTS events each, and a best-vs-runner-up margin >=
   m2_tiebreak_margin; the vote goes to the argmin hypothesis.
3. Verdicts per accused variant: quorum + supermajority winner defeats
   accused variants whose own share is <= 1 - majority; upheld when the
   accused itself wins; split/no_quorum/no_ballots/window_impure retain.
   Only 'defeated' acts, and only in require mode.
4. Window purity: representatives disagreeing on non-own exon boundaries
   inside the ballot windows abstain the whole contest.
5. Ballot ablation modes: screen defeats on nomination alone; shuffle
   permutes scores deterministically per read.
6. Artifact merge deterministic; header-only when nothing nominated.
7. Config validation: mode/ballot choices; m2_em + signal + families
   discovery required; CLI defaults.
"""

from __future__ import annotations

import pytest

from fin.analysis.junction_election import (
    ElectionRow,
    merge_election_shards,
    nominate_contests,
    run_junction_election,
    write_election_tsv,
)
from fin.analysis.quantification import QuantResult
from fin.candidates.dataclasses import (
    CandidateSet,
    GenomicInterval,
    IntronChain,
    TranscriptCandidate,
)
from fin.pipeline.config import PipelineConfig


def _cand(cid, chain, start=1000, end=3000, source="novel"):
    return TranscriptCandidate(
        candidate_id=cid,
        intron_chain=IntronChain(introns=tuple(chain)),
        three_prime_pos=end,
        sequence="ACGT" * 200,
        source=source,
        supporting_read_ids=set(),
        chrom="chr1",
        strand="+",
        start=start,
        end=end,
    )


def _qr(cid, abundance, assigned=()):
    return QuantResult(
        candidate_id=cid, abundance=float(abundance), confidence=1.0,
        num_assigned_reads=len(assigned), source="novel",
        assigned_read_ids=tuple(assigned),
    )


TRUE_J = (1500, 2000)
FAKE_J = (1503, 1998)           # within 15 bp of TRUE_J
FAKE_J2 = (1497, 2003)          # second accused variant, same neighbourhood
FAR_J = (2500, 2800)


class TestNomination:
    def test_basic_contest(self):
        strong = _cand("strong", (TRUE_J,))
        weak = _cand("weak", (FAKE_J,))
        support = {TRUE_J: (30, 25), FAKE_J: (3, 0)}
        quant = {"strong": _qr("strong", 30), "weak": _qr("weak", 3)}
        contests, trunc = nominate_contests([strong, weak], quant, support)
        assert trunc == 0 and len(contests) == 1
        c = contests[0]
        assert c["nominees"] == sorted([FAKE_J, TRUE_J])
        assert list(c["accused"]) == [FAKE_J]
        assert [x.candidate_id for x in c["accused"][FAKE_J]] == ["weak"]
        assert c["reps"][TRUE_J].candidate_id == "strong"

    def test_two_accused_merge_into_one_kway_contest(self):
        strong = _cand("strong", (TRUE_J,))
        w1 = _cand("w1", (FAKE_J,))
        w2 = _cand("w2", (FAKE_J2,))
        support = {TRUE_J: (30, 25), FAKE_J: (3, 0), FAKE_J2: (2, 1)}
        quant = {"strong": _qr("strong", 30), "w1": _qr("w1", 3),
                 "w2": _qr("w2", 2)}
        contests, _ = nominate_contests([strong, w1, w2], quant, support)
        assert len(contests) == 1
        c = contests[0]
        assert set(c["accused"]) == {FAKE_J, FAKE_J2}
        assert c["nominees"] == sorted([FAKE_J2, FAKE_J, TRUE_J])

    def test_no_contest_when_accused_clean_enough(self):
        support = {TRUE_J: (30, 25), FAKE_J: (5, 2)}
        cands = [_cand("s", (TRUE_J,)), _cand("w", (FAKE_J,))]
        quant = {"s": _qr("s", 30), "w": _qr("w", 5)}
        assert nominate_contests(cands, quant, support) == ([], 0)

    def test_no_contest_when_competitor_weak(self):
        support = {TRUE_J: (3, 2), FAKE_J: (3, 0)}
        cands = [_cand("s", (TRUE_J,)), _cand("w", (FAKE_J,))]
        quant = {"s": _qr("s", 3), "w": _qr("w", 3)}
        assert nominate_contests(cands, quant, support) == ([], 0)

    def test_far_junction_not_a_competitor(self):
        support = {FAR_J: (30, 25), FAKE_J: (3, 0)}
        cands = [_cand("s", (FAR_J,)), _cand("w", (FAKE_J,))]
        quant = {"s": _qr("s", 30), "w": _qr("w", 3)}
        assert nominate_contests(cands, quant, support) == ([], 0)

    def test_carrier_of_champion_never_accused(self):
        both = _cand("both", (FAKE_J, TRUE_J))
        strong = _cand("strong", (TRUE_J,))
        support = {TRUE_J: (30, 25), FAKE_J: (3, 0)}
        quant = {"both": _qr("both", 5), "strong": _qr("strong", 30)}
        assert nominate_contests([both, strong], quant, support) == ([], 0)

    def test_champion_junction_must_have_carrier(self):
        weak = _cand("weak", (FAKE_J,))
        support = {TRUE_J: (30, 25), FAKE_J: (3, 0)}
        quant = {"weak": _qr("weak", 3)}
        assert nominate_contests([weak], quant, support) == ([], 0)


def _setup(mode, votes, ballot="signal", extra_chain=None):
    """One-contest interval. votes[read_id] = junction the read's signal
    prefers (TRUE_J/FAKE_J) or None (unusable)."""
    strong = _cand("strong", (TRUE_J,) if extra_chain is None
                   else (TRUE_J, extra_chain))
    weak = _cand("weak", (FAKE_J,))
    reads = sorted(votes)
    quant = [
        _qr("strong", 30, assigned=[r for r in reads if not r.startswith("w")]),
        _qr("weak", 3, assigned=[r for r in reads if r.startswith("w")]),
    ]
    cs = CandidateSet(
        interval=GenomicInterval(chrom="chr1", start=0, end=5000, strand="+"),
        candidates=[strong, weak],
        read_ids=set(reads),
        read_sequences={r: "ACGT" * 100 for r in reads},
        junction_support={TRUE_J: (30, 25), FAKE_J: (3, 0)},
    )
    cs.read_spans = {r: (900, 3100) for r in reads}

    def scorer(rid, seq, cand, windows):
        v = votes.get(rid)
        if v is None:
            return float("nan"), 0
        prefers_this = (v == TRUE_J) == (cand.candidate_id == "strong")
        return (1.0 if prefers_this else 2.0), 5

    cfg = PipelineConfig(
        bam_path="/tmp/x.bam", junction_election_mode=mode,
        junction_election_ballot=ballot,
    )
    return cfg, [strong, weak], quant, cs, scorer


class TestElection:
    def test_defeated_drops_in_require(self):
        votes = {f"s{i}": TRUE_J for i in range(8)}
        votes["w0"] = FAKE_J
        cfg, cands, quant, cs, scorer = _setup("require", votes)
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == {"weak"}
        assert len(rows) == 1
        r = rows[0]
        assert r.verdict == "defeated" and r.action == "drop"
        assert r.votes_winner == 8 and r.votes_accused == 1
        assert r.winner_junction == f"{TRUE_J[0]}-{TRUE_J[1]}"

    def test_audit_never_drops(self):
        votes = {f"s{i}": TRUE_J for i in range(8)}
        cfg, cands, quant, cs, scorer = _setup("audit", votes)
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == set()
        assert rows[0].verdict == "defeated" and rows[0].action == "none"

    def test_upheld_keeps(self):
        votes = {f"w{i}": FAKE_J for i in range(8)}
        cfg, cands, quant, cs, scorer = _setup("require", votes)
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == set() and rows[0].verdict == "upheld"

    def test_split_abstains(self):
        votes = {f"s{i}": TRUE_J for i in range(5)}
        votes.update({f"w{i}": FAKE_J for i in range(4)})
        cfg, cands, quant, cs, scorer = _setup("require", votes)
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == set() and rows[0].verdict == "split"

    def test_no_quorum_abstains(self):
        votes = {f"s{i}": TRUE_J for i in range(4)}
        cfg, cands, quant, cs, scorer = _setup("require", votes)
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == set() and rows[0].verdict == "no_quorum"

    def test_unusable_signal_abstains(self):
        votes = {f"s{i}": TRUE_J for i in range(6)}
        votes.update({f"x{i}": None for i in range(10)})
        cfg, cands, quant, cs, scorer = _setup("require", votes)
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert rows[0].n_ballots == 6 and drops == {"weak"}

    def test_margin_too_small_abstains(self):
        votes = {f"s{i}": TRUE_J for i in range(8)}
        cfg, cands, quant, cs, _ = _setup("require", votes)
        cfg.m2_tiebreak_margin = 10.0    # bigger than the 1.0 score gap

        def close_scorer(rid, seq, cand, windows):
            return (1.0 if cand.candidate_id == "strong" else 1.5), 5

        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=close_scorer)
        assert drops == set() and rows[0].verdict == "no_quorum"

    def test_min_events_enforced(self):
        votes = {f"s{i}": TRUE_J for i in range(8)}
        cfg, cands, quant, cs, _ = _setup("require", votes)

        def thin_scorer(rid, seq, cand, windows):
            return (1.0 if cand.candidate_id == "strong" else 2.0), 1

        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=thin_scorer)
        assert drops == set() and rows[0].verdict == "no_quorum"

    def test_window_impure_abstains(self):
        # strong rep carries an extra junction whose acceptor falls inside
        # the ballot windows -> fingerprints differ -> abstain
        votes = {f"s{i}": TRUE_J for i in range(8)}
        cfg, cands, quant, cs, scorer = _setup(
            "require", votes, extra_chain=(1400, 1495))
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == set() and rows[0].verdict == "window_impure"

    def test_no_junction_support_no_election(self):
        cfg, cands, quant, cs, scorer = _setup("require", {"s0": TRUE_J})
        cs.junction_support = None
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == set() and rows == []


class TestDefenseBurden:
    def _cfg_defense(self, votes):
        cfg, cands, quant, cs, scorer = _setup("require", votes)
        cfg.junction_election_burden = "defense"
        return cfg, cands, quant, cs, scorer

    def test_exonerated_accused_survives(self):
        votes = {f"w{i}": FAKE_J for i in range(8)}
        cfg, cands, quant, cs, scorer = self._cfg_defense(votes)
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == set() and rows[0].verdict == "upheld"

    def test_abstention_follows_screen_and_drops(self):
        votes = {f"s{i}": TRUE_J for i in range(4)}   # below quorum
        cfg, cands, quant, cs, scorer = self._cfg_defense(votes)
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == {"weak"}
        assert rows[0].verdict == "no_quorum" and rows[0].action == "drop"

    def test_convicted_drops(self):
        votes = {f"s{i}": TRUE_J for i in range(8)}
        cfg, cands, quant, cs, scorer = self._cfg_defense(votes)
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == {"weak"} and rows[0].verdict == "defeated"

    def test_prosecution_default_retains_on_abstention(self):
        votes = {f"s{i}": TRUE_J for i in range(4)}
        cfg, cands, quant, cs, scorer = _setup("require", votes)
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == set() and rows[0].verdict == "no_quorum"

    def test_audit_defense_never_drops(self):
        votes = {f"s{i}": TRUE_J for i in range(4)}
        cfg, cands, quant, cs, scorer = _setup("audit", votes)
        cfg.junction_election_burden = "defense"
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == set() and rows[0].action == "none"


class TestBallotModes:
    def test_screen_defeats_without_ballots(self):
        cfg, cands, quant, cs, scorer = _setup(
            "require", {"s0": TRUE_J}, ballot="screen")
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        assert drops == {"weak"}
        assert rows[0].verdict == "defeated" and rows[0].n_ballots == 0

    def test_shuffle_is_deterministic(self):
        votes = {f"s{i}": TRUE_J for i in range(12)}
        out = []
        for _ in range(2):
            cfg, cands, quant, cs, scorer = _setup(
                "require", votes, ballot="shuffle")
            drops, rows = run_junction_election(
                cfg, cands, quant, cs, interval_name="i", scorer=scorer)
            out.append((tuple(sorted(drops)), rows[0].votes_winner,
                        rows[0].votes_accused, rows[0].verdict))
        assert out[0] == out[1]


class TestArtifact:
    def test_merge_deterministic_and_header_only(self, tmp_path):
        shard_dir = tmp_path / "_jelection_shards"
        shard_dir.mkdir()
        r1 = ElectionRow(
            interval="i2", contest_id="i2#0", accused_junction="10-20",
            nominees="10-20,12-22", accused_clean=0,
            best_competitor_clean=9, accused_ids="a",
            verdict="defeated", action="drop")
        r2 = ElectionRow(
            interval="i1", contest_id="i1#0", accused_junction="5-9",
            nominees="5-9,6-9", accused_clean=1, best_competitor_clean=8,
            accused_ids="b")
        write_election_tsv(shard_dir / "z.tsv", [r1])
        write_election_tsv(shard_dir / "a.tsv", [r2])
        out = merge_election_shards(tmp_path)
        lines = (tmp_path / "junction_election.tsv").read_text().splitlines()
        assert out and len(lines) == 3
        assert lines[1].startswith("i1\t") and lines[2].startswith("i2\t")
        assert not shard_dir.exists()
        assert merge_election_shards(tmp_path / "missing") is None
        (tmp_path / "wd").mkdir()
        merge_election_shards(tmp_path / "wd")
        lines2 = (tmp_path / "wd" / "junction_election.tsv").read_text().splitlines()
        assert len(lines2) == 1 and lines2[0].startswith("interval\t")


class TestConfig:
    def test_bad_mode(self, tmp_path):
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        cfg = PipelineConfig(bam_path=str(bam), junction_election_mode="zzz")
        with pytest.raises(ValueError, match="junction_election_mode"):
            cfg.validate()

    def test_bad_ballot(self, tmp_path):
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        cfg = PipelineConfig(
            bam_path=str(bam), junction_election_ballot="zzz")
        with pytest.raises(ValueError, match="junction_election_ballot"):
            cfg.validate()

    def test_requires_m2_em(self, tmp_path):
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        sig = tmp_path / "x.blow5"; sig.write_bytes(b"")
        cfg = PipelineConfig(
            bam_path=str(bam), signal_path=str(sig), quant_mode="argmax",
            junction_election_mode="audit")
        with pytest.raises(ValueError, match="requires --quant-mode m2_em"):
            cfg.validate()

    def test_requires_signal(self, tmp_path):
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        cfg = PipelineConfig(
            bam_path=str(bam), junction_election_mode="require")
        with pytest.raises(ValueError, match="requires --signal"):
            cfg.validate()

    def test_requires_families_discovery(self, tmp_path):
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        sig = tmp_path / "x.blow5"; sig.write_bytes(b"")
        cfg = PipelineConfig(
            bam_path=str(bam), signal_path=str(sig),
            junction_election_mode="audit", chain_cluster_discovery=False)
        with pytest.raises(ValueError, match="families"):
            cfg.validate()

    def test_bad_burden(self, tmp_path):
        bam = tmp_path / "x.bam"; bam.write_bytes(b"")
        cfg = PipelineConfig(
            bam_path=str(bam), junction_election_burden="zzz")
        with pytest.raises(ValueError, match="junction_election_burden"):
            cfg.validate()

    def test_cli_defaults(self):
        from fin.cli import main

        opt = next(
            p for p in main.params if p.name == "junction_election_mode")
        assert opt.default == "off"
        opt2 = next(
            p for p in main.params if p.name == "junction_election_ballot")
        assert opt2.default == "signal"
        opt3 = next(
            p for p in main.params if p.name == "junction_election_burden")
        assert opt3.default == "prosecution"


class TestElectorateCoversAllCarriers:
    """The electorate must be drawn from EVERY candidate carrying a nominated
    junction, not just its representative.

    A junction's support is spread over its carriers; taking only the
    representative's reads under-counts a well-supported junction into
    `no_quorum`, which the defense burden then converts into a drop. This is the
    P1 finding from the 2026-09-02 Codex review.
    """

    def _two_carrier_setup(self, mode="require"):
        # TRUE_J is carried by two candidates that split its reads; FAKE_J by one.
        strong_a = _cand("strong_a", (TRUE_J,))
        strong_b = _cand("strong_b", (TRUE_J,))
        weak = _cand("weak", (FAKE_J,))
        a_reads = [f"a{i}" for i in range(4)]
        b_reads = [f"b{i}" for i in range(4)]
        w_reads = ["w0"]
        reads = a_reads + b_reads + w_reads
        quant = [_qr("strong_a", 30, assigned=a_reads),
                 _qr("strong_b", 20, assigned=b_reads),
                 _qr("weak", 3, assigned=w_reads)]
        cs = CandidateSet(
            interval=GenomicInterval(chrom="chr1", start=0, end=5000, strand="+"),
            candidates=[strong_a, strong_b, weak],
            read_ids=set(reads),
            read_sequences={r: "ACGT" * 100 for r in reads},
            junction_support={TRUE_J: (30, 25), FAKE_J: (3, 0)},
        )
        cs.read_spans = {r: (900, 3100) for r in reads}
        seen = []

        def scorer(rid, seq, cand, windows):
            seen.append(rid)
            prefers = (not rid.startswith("w")) == (
                cand.candidate_id != "weak")
            return (1.0 if prefers else 2.0), 5

        cfg = PipelineConfig(bam_path="/tmp/x.bam", junction_election_mode=mode,
                             junction_election_ballot="signal")
        return cfg, [strong_a, strong_b, weak], quant, cs, scorer, seen

    def test_all_carriers_reads_enter_the_electorate(self):
        cfg, cands, quant, cs, scorer, seen = self._two_carrier_setup()
        drops, rows = run_junction_election(
            cfg, cands, quant, cs, interval_name="i", scorer=scorer)
        # the representative alone owns 4 reads; both carriers together own 8,
        # and only the full electorate reaches the quorum of 5.
        assert rows[0].n_electorate == 9
        voters = {r for r in seen if not r.startswith("w")}
        assert voters == {f"a{i}" for i in range(4)} | {f"b{i}" for i in range(4)}
        assert rows[0].verdict == "defeated" and drops == {"weak"}

    def test_contest_records_every_carrier(self):
        from fin.analysis.junction_election import nominate_contests

        cfg, cands, quant, cs, scorer, _ = self._two_carrier_setup()
        contests, _ = nominate_contests(
            cands, {q.candidate_id: q for q in quant}, cs.junction_support)
        assert len(contests) == 1
        ids = {c.candidate_id for c in contests[0]["carriers"][TRUE_J]}
        assert ids == {"strong_a", "strong_b"}


class TestScorerIsProcessCached:
    def test_same_config_returns_the_same_scorer(self):
        from fin.analysis import junction_election as je

        je._SCORER_CACHE.clear()
        calls = []

        def fake_build(config, num_thread=0):
            calls.append(config)
            return lambda *a: (1.0, 5)

        orig = je._build_signal_scorer
        je._build_signal_scorer = fake_build
        try:
            cfg = PipelineConfig(bam_path="/tmp/x.bam", krill_pore="rna002")
            a = je.default_signal_scorer(cfg)
            b = je.default_signal_scorer(cfg)
            assert a is b and len(calls) == 1     # built once per process
            cfg2 = PipelineConfig(bam_path="/tmp/x.bam", krill_pore="rna004")
            je.default_signal_scorer(cfg2)
            assert len(calls) == 2                # a different pore rebuilds
        finally:
            je._build_signal_scorer = orig
            je._SCORER_CACHE.clear()
