"""Unit tests for step-5a per-read junction correction.

The signal ballot is stubbed with a deterministic fake aligner so the structural
contract (registry, nomination screen, hypothesis construction, chain validity,
apply-all-at-once, audit inertness, determinism) is tested without krill.
"""
from __future__ import annotations

import numpy as np
import pytest

from fin.analysis.junction_correction import (
    JunctionRecord,
    aggregate_proposals,
    build_registry,
    build_transcript,
    chain_is_valid,
    correct_read_chains,
    merge_correction_shards,
    nominate,
    substitute,
    write_correction_tsv,
)
from fin.candidates.dataclasses import IntronChain


class _Cfg:
    junction_correction_mode = "require"
    junction_correction_tol_bp = 15
    junction_correction_flank_bp = 4
    junction_correction_margin = 10.0
    junction_correction_runner_up_margin = 10.0


class _Interval:
    chrom, start, end, strand = "chr1", 0, 5000, "+"
    region_string = "chr1:0-5000"


def _read(name, start=100, end=1200, reverse=False):
    return {"query_name": name, "reference_start": start, "reference_end": end,
            "is_reverse": reverse}


def _rc(name, chain, **kw):
    return (_read(name, **kw), IntronChain(introns=tuple(chain)))


# --------------------------------------------------------------------- registry
def test_registry_counts_raw_and_indel_clean_separately():
    rcs = [_rc("a", [(200, 300)]), _rc("b", [(200, 300)]), _rc("c", [(200, 300)])]
    indels = {"a": (), "b": ((295, 297),), "c": ()}      # b has an indel 3 bp from 300
    reg = build_registry(rcs, indels)
    assert reg[(200, 300)].raw == 3
    assert reg[(200, 300)].clean == 2
    assert reg[(200, 300)].read_ids == {"a", "b", "c"}


def test_registry_read_without_indel_record_counts_raw_only():
    reg = build_registry([_rc("a", [(200, 300)])], {})
    assert reg[(200, 300)].raw == 1 and reg[(200, 300)].clean == 0


# --------------------------------------------------------------------- nomination
def _reg(**kw):
    return {j: JunctionRecord(raw=r, clean=c) for j, (r, c) in kw.items()}


def test_nomination_has_no_support_threshold():
    """Every registry neighbour within tolerance is an alternative; the registry only
    records, the signal decides."""
    reg = {(200, 300): JunctionRecord(raw=10, clean=1),
           (200, 303): JunctionRecord(raw=1, clean=0)}
    assert nominate((200, 300), reg, tol_bp=15) == [(200, 303)]


def test_nomination_respects_tolerance_on_both_coordinates():
    reg = {(200, 300): JunctionRecord(raw=2, clean=0),
           (200, 316): JunctionRecord(raw=30, clean=9),   # acceptor 16 bp away
           (184, 300): JunctionRecord(raw=30, clean=9)}   # donor 16 bp away
    assert nominate((200, 300), reg, tol_bp=15) == []
    reg[(200, 314)] = JunctionRecord(raw=30, clean=9)
    assert nominate((200, 300), reg, tol_bp=15) == [(200, 314)]


def test_nomination_is_deterministically_sorted_and_kway():
    reg = {(200, 300): JunctionRecord(raw=1, clean=0),
           (200, 305): JunctionRecord(raw=30, clean=9),
           (198, 300): JunctionRecord(raw=30, clean=9),
           (203, 297): JunctionRecord(raw=2, clean=2)}
    assert nominate((200, 300), reg, tol_bp=15) == [(198, 300), (200, 305), (203, 297)]


# --------------------------------------------------------------------- hypotheses
def test_substitute_and_chain_validity():
    ch = ((200, 300), (400, 500))
    assert substitute(ch, (200, 300), (200, 303)) == ((200, 303), (400, 500))
    assert chain_is_valid(100, 900, ((200, 303), (400, 500)))
    assert not chain_is_valid(100, 900, ((400, 500), (200, 300)))     # out of order
    assert not chain_is_valid(100, 900, ((200, 450), (400, 500)))     # overlapping
    assert not chain_is_valid(100, 250, ((200, 300),))                # past the read end


def test_build_transcript_plus_strand_maps_positions():
    g = "ACGT" * 200
    seq, gmap = build_transcript(g, 10, 40, ((20, 30),), "+")
    assert seq == g[10:20] + g[30:40]
    assert list(gmap[:10]) == list(range(10, 20))
    assert list(gmap[10:]) == list(range(30, 40))


def test_build_transcript_minus_strand_is_reverse_complemented():
    g = "AAAACCCCGGGGTTTT" * 20
    seq, gmap = build_transcript(g, 10, 40, ((20, 30),), "-")
    plus, pmap = build_transcript(g, 10, 40, ((20, 30),), "+")
    comp = str.maketrans("ACGT", "TGCA")
    assert seq == plus.translate(comp)[::-1]
    assert list(gmap) == list(pmap)[::-1]


# --------------------------------------------------------------------- driver
class _FakeAligner:
    """Returns a fixed NLL per hypothesis sequence via a caller-supplied rule."""

    def __init__(self, prefer):
        self.prefer = prefer          # callable(sequence) -> mean current
        self.calls = 0

    def align(self, items):
        out = []
        for it in items:
            self.calls += 1
            n = 400
            mu = self.prefer(it["sequence"])
            # one event per transcript base over 60..459: covers every contested
            # junction region in these fixtures plus >= 4 bp shared flank each side
            out.append(dict(
                status=0, init_scale=1.0, init_shift=0.0, var=1.0,
                position=np.arange(n, dtype=np.int64) + 60,
                start_idx=np.arange(n, dtype=np.int64) * 10,
                end_idx=(np.arange(n, dtype=np.int64) + 1) * 10,
                model_mean=np.full(n, mu), model_stdv=np.full(n, 2.0),
                event_level_mean=np.full(n, 100.0), event_stdv=np.full(n, 1.0),
                model_kmer=np.array(["AAAAA"] * n)))
        return out


def _signal(_read_id):
    return dict(signal=np.zeros(500, dtype="float32"), digitisation=8192.0,
                offset=0.0, range=1500.0, sample_rate=4000.0)


def _setup(margin=10.0, mode="require"):
    """3 reads carry the wrong junction (200,300); 12 carry (200,303)."""
    rcs = [_rc(f"bad{i}", [(200, 300)]) for i in range(3)]
    rcs += [_rc(f"good{i}", [(200, 303)]) for i in range(12)]
    indels = {r["query_name"]: () for r, _ in rcs}
    cfg = _Cfg()
    cfg.junction_correction_margin = margin
    cfg.junction_correction_mode = mode
    return rcs, indels, cfg


def _genome():
    return "ACGTTGCA" * 500


def _run(rcs, indels, cfg, aligner, reader=_signal):
    return correct_read_chains(
        rcs, indels, config=cfg, interval=_Interval(), genome_seq=_genome(),
        signal_reader=reader, aligner=aligner, krill_mod=None)


def test_correction_is_additive_it_never_touches_the_input_chains():
    rcs, indels, cfg = _setup()
    aligner = _FakeAligner(lambda seq: 100.0 if len(seq) == 997 else 130.0)   # alt (200-303) fits
    snapshot = [tuple(c.introns) for _, c in rcs]
    res = _run(rcs, indels, cfg, aligner)
    assert [tuple(c.introns) for _, c in rcs] == snapshot
    # the three "bad" reads propose the (200,303) chain; the 12 good reads see the
    # alternative (200,300) fit WORSE and propose nothing
    props = {p.read_id: p.corrected_chain for p in res.proposals}
    assert props == {f"bad{i}": ((200, 303),) for i in range(3)}
    assert all(p.own_margin >= 10 for p in res.proposals)
    assert aggregate_proposals(res.proposals, min_reads=2) == {((200, 303),): ["bad0", "bad1", "bad2"]}
    assert aggregate_proposals(res.proposals, min_reads=4) == {}


def test_audit_and_require_produce_the_same_proposals():
    rcs, indels, cfg = _setup(mode="audit")
    aligner = _FakeAligner(lambda seq: 100.0 if len(seq) == 997 else 130.0)
    res = _run(rcs, indels, cfg, aligner)
    assert len(res.proposals) == 3 and res.rows


def test_off_mode_does_no_work_at_all():
    rcs, indels, cfg = _setup(mode="off")
    aligner = _FakeAligner(lambda seq: 100.0)
    res = _run(rcs, indels, cfg, aligner)
    assert res.proposals == [] and res.rows == [] and aligner.calls == 0


def test_no_proposal_when_the_read_own_hypothesis_fits_better():
    rcs, indels, cfg = _setup()
    aligner = _FakeAligner(lambda seq: 100.0 if len(seq) == 1000 else 130.0)   # own (200-300) fits
    res = _run(rcs, indels, cfg, aligner)
    assert all(p.read_id.startswith("good") for p in res.proposals)          # bad reads keep their own
    assert all(r.decision == "abstain" for r in res.rows if r.read_id.startswith("bad"))


def test_own_margin_below_threshold_abstains():
    rcs, indels, cfg = _setup(margin=1e9)
    aligner = _FakeAligner(lambda seq: 100.0 if len(seq) == 997 else 130.0)
    res = _run(rcs, indels, cfg, aligner)
    assert res.proposals == []
    assert any(r.reason == "own_margin_below_threshold" for r in res.rows)


def test_runner_up_margin_abstains_when_two_alternatives_tie():
    """Three-way contest: own (200,300) vs (200,303) vs (200,306); the two alternatives
    score identically -> best cannot beat the runner-up -> abstain."""
    rcs = [_rc("r0", [(200, 300)]), _rc("a", [(200, 303)]), _rc("b", [(200, 306)])]
    indels = {r["query_name"]: () for r, _ in rcs}
    cfg = _Cfg()
    aligner = _FakeAligner(lambda seq: 130.0 if len(seq) == 1000 else 100.0)  # both alts fit equally
    res = _run(rcs, indels, cfg, aligner)
    assert all(p.read_id != "r0" for p in res.proposals)
    assert any(r.read_id == "r0" and r.reason == "runner_up_margin_below_threshold" for r in res.rows)


def test_one_proposal_per_read_even_with_two_correctable_junctions():
    rcs = [_rc("r0", [(200, 300), (500, 600)]), _rc("x", [(200, 303), (500, 600)]),
           _rc("y", [(200, 300), (500, 604)])]
    indels = {r["query_name"]: () for r, _ in rcs}
    cfg = _Cfg()
    # every alternative sequence (shorter than own 900-nt-ish chain) fits better
    own_len = len(_genome()[100:200] + _genome()[300:500] + _genome()[600:1200])
    aligner = _FakeAligner(lambda seq: 130.0 if len(seq) == own_len else 100.0)
    res = _run(rcs, indels, cfg, aligner)
    mine = [p for p in res.proposals if p.read_id == "r0"]
    assert len(mine) == 1
    assert sum(1 for p in mine[0].corrected_chain if p not in ((200, 300), (500, 600))) == 1  # single edit


def test_reads_with_no_signal_abstain_and_are_recorded():
    rcs, indels, cfg = _setup()
    aligner = _FakeAligner(lambda seq: 100.0)
    res = _run(rcs, indels, cfg, aligner, reader=lambda rid: None)
    assert res.proposals == []
    assert res.rows and all(r.reason == "no_signal" for r in res.rows)


def test_decisions_do_not_depend_on_read_order():
    rcs, indels, cfg = _setup()
    aligner = _FakeAligner(lambda seq: 100.0 if len(seq) == 997 else 130.0)
    a = _run(list(rcs), indels, cfg, aligner)
    b = _run(list(reversed(rcs)), indels, cfg, aligner)
    key = lambda p: (p.read_id, p.corrected_chain)
    assert sorted(map(key, a.proposals)) == sorted(map(key, b.proposals))
    assert sorted(r.as_row() for r in a.rows) == sorted(r.as_row() for r in b.rows)


def test_no_nominations_means_no_krill_call():
    """No registry neighbour within tolerance -> not even the capability probe runs."""
    rcs = [_rc(f"a{i}", [(200, 300)]) for i in range(5)]
    indels = {r["query_name"]: () for r, _ in rcs}
    aligner = _FakeAligner(lambda seq: 100.0)
    res = _run(rcs, indels, _Cfg(), aligner)
    assert res.proposals == [] and res.rows == [] and aligner.calls == 0


# --------------------------------------------------------------------- artifact
def test_shard_merge_is_deterministic_and_removes_the_shard_dir(tmp_path):
    from fin.analysis.junction_correction import CorrectionRow

    d = tmp_path / "_jcorrection_shards"
    for k, iv in enumerate(("chr1:0-10", "chr1:10-20")):
        write_correction_tsv(d / f"s{k}.tsv", [CorrectionRow(
            interval=iv, read_id=f"r{k}", strand="+", own_junction="1-2",
            alt_junction="1-3", shift_donor=0, shift_acceptor=1, own_raw=1,
            own_clean=0, alt_raw=9, alt_clean=5)])
    path = merge_correction_shards(tmp_path)
    lines = open(path).read().strip().split("\n")
    assert lines[0].split("\t")[0] == "interval" and len(lines) == 3
    assert lines[1:] == sorted(lines[1:])
    assert not d.exists()


# --------------------------------------------------------------------- discovery append
def test_discovery_appends_a_new_candidate_and_never_loses_baseline(monkeypatch):
    """Instantiation inside discovery: baseline candidates untouched; a chain proposed
    by >= min_reads distinct reads becomes an extra novel candidate carrying its voters
    as correction support only."""
    from fin.analysis.junction_correction import CorrectionProposal, JunctionCorrectionResult
    from fin.candidates import discovery as disc

    genome = "ACGTTGCA" * 500
    reads = []
    for i in range(3):
        reads.append({"query_name": f"bad{i}", "reference_start": 100, "reference_end": 1200,
                      "reference_end_": 1200, "is_reverse": False, "mapping_quality": 60,
                      "cigartuples": [(0, 100), (3, 100), (0, 900)]})           # chain (200,300)
    for i in range(4):
        reads.append({"query_name": f"other{i}", "reference_start": 100, "reference_end": 1200,
                      "is_reverse": False, "mapping_quality": 60,
                      "cigartuples": [(0, 100), (3, 100), (0, 200), (3, 100), (0, 600)]})  # (200,300),(500,600)

    class _Iv:
        chrom, start, end, strand = "chr1", 0, 5000, "+"
        region_string = "chr1:0-5000"

    def ballot(read_chains, indels_by_read):
        res = JunctionCorrectionResult()
        for rid in ("bad0", "bad1"):
            res.proposals.append(CorrectionProposal(
                read_id=rid, original_chain=((200, 300),), corrected_chain=((200, 303),),
                changed_junction=(200, 300), replacement=(200, 303), own_margin=12.0,
                runner_up_margin=float("inf"), confidence=12.0, reference_start=100,
                reference_end=1200))
        res.proposals.append(CorrectionProposal(                       # lone voter -> no candidate
            read_id="bad2", original_chain=((200, 300),), corrected_chain=((200, 297),),
            changed_junction=(200, 300), replacement=(200, 297), own_margin=12.0,
            runner_up_margin=float("inf"), confidence=12.0, reference_start=100, reference_end=1200))
        return res

    class _Apply:
        min_reads = 2

    cs0 = disc._chain_cluster_candidates(_Iv(), reads, {r["query_name"] for r in reads}, {}, None,
                                         genome, "+", 1, wobble_bp=6, cassette_max_exon_bp=70,
                                         clustering="families")
    cs1 = disc._chain_cluster_candidates(_Iv(), reads, {r["query_name"] for r in reads}, {}, None,
                                         genome, "+", 1, wobble_bp=6, cassette_max_exon_bp=70,
                                         clustering="families", junction_correction=ballot,
                                         junction_correction_apply=_Apply())
    chains0 = {tuple(c.intron_chain.introns): c for c in cs0.candidates}
    chains1 = {tuple(c.intron_chain.introns): c for c in cs1.candidates}
    assert set(chains0) <= set(chains1)                                   # nothing lost
    assert set(chains1) - set(chains0) == {((200, 303),)}                 # exactly one new chain
    new = chains1[((200, 303),)]
    assert new.source == "novel"
    assert new.correction_supporting_read_ids == {"bad0", "bad1"}
    assert new.raw_supporting_read_ids == set()
    assert new.complete_chain_weakest_indel_clean_depth == 0
    # every baseline candidate keeps its raw reads untouched (here (200,300) was an
    # exact sub-chain of (200,300),(500,600) and is folded into it by collapse)
    assert set(chains0) == {((200, 300), (500, 600))}
    for ch, c0 in chains0.items():
        assert chains1[ch].raw_supporting_read_ids == c0.raw_supporting_read_ids
        assert chains1[ch].supporting_read_ids == c0.supporting_read_ids
    # audit (no apply object) appends nothing
    cs2 = disc._chain_cluster_candidates(_Iv(), reads, {r["query_name"] for r in reads}, {}, None,
                                         genome, "+", 1, wobble_bp=6, cassette_max_exon_bp=70,
                                         clustering="families", junction_correction=ballot)
    assert {tuple(c.intron_chain.introns) for c in cs2.candidates} == set(chains0)
