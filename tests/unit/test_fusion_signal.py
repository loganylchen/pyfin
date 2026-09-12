"""Shared signal arithmetic, strand projection and conservative fusion gates."""

from dataclasses import replace
import random
from types import SimpleNamespace
import sys

import numpy as np
import pytest

from fin.fusion.analysis import FusionAnalysis, FusionReadAudit, FusionSettings, finalize_chains
from fin.fusion.annotation import AnnotatedReadPath, FusionArm
from fin.fusion.correction import PathCorrection, PathJunction, path_junctions
from fin.fusion.junctions import JunctionContrast, JunctionDecision
from fin.fusion.records import AlignmentBlock
from fin.fusion.signal import (
    FusionSignalSettings, KrillFusionBackend, _move_path, score_signal_contrast,
    validate_fusion_signal, virtual_candidates,
)
from fin.fusion.validation import stitched_template
from fin.scoring.m2_junction_nll import _mean_nll_in_window, diff_junction_windows, tx2genome_array


def _fixture(strands=("+", "+", "+")):
    generator = random.Random(273)
    genome = {chrom: "".join(generator.choices("ACGT", k=1200)) for chrom in "ABC"}
    arms = []
    for index, (chrom, strand) in enumerate(zip("ABC", strands)):
        intervals = [(200, 350), (500, 650)] if index == 1 else [(200, 500)]
        if strand == "-":
            intervals.reverse()
        blocks = []
        query = index * 300
        for start, end in intervals:
            blocks.append(AlignmentBlock(start, end, query, query + end - start, end - start))
            query += end - start
        arms.append(FusionArm(chrom, strand, tuple(blocks), (chrom,)))
    path = AnnotatedReadPath("read", tuple(arms), "interchromosomal_candidate", "resolved")
    return genome, path


def _record(candidate, deviation=0.0):
    size = len(candidate.sequence)
    return dict(variant_label=candidate.candidate_id, status=0, position=list(range(size)),
                reference_kmer=["AAAAA"] * size, model_kmer=["AAAAA"] * size,
                event_level_mean=[deviation] * size, model_mean=[0.0] * size, model_stdv=[1.0] * size)


@pytest.mark.parametrize("strands", [("+", "+", "+"), ("+", "-", "+"), ("-", "-", "-")])
def test_virtual_axis_preserves_full_sequence_and_exon_lengths(strands):
    genome, path = _fixture(strands)
    junction = tuple(path_junctions(path))[1]
    positions = junction[3].exit_position + 1, junction[4].entry_position
    candidates = virtual_candidates([path, _move_path(path, junction, positions)], genome)
    for candidate in candidates:
        projected = tx2genome_array(candidate, np.arange(len(candidate.sequence)))
        assert np.all(projected >= 0)
        assert np.all(np.diff(projected) > 0)
        assert np.count_nonzero(np.diff(projected) > 1) == 3
    assert candidates[0].sequence == stitched_template(path, genome)[0]
    assert len(diff_junction_windows(list(candidates), flank=8)) == 2


def _contrast():
    genome, path = _fixture()
    junction = next(path_junctions(path))
    alternate = _move_path(path, junction, (501, 200))
    return virtual_candidates([path, alternate], genome)


@pytest.mark.parametrize("deviations,status", [((0, 2), "supported"), ((2, 0), "conflict"), ((1, 1), "abstain")])
def test_signal_contrast_reuses_exact_assembly_sum(deviations, status):
    candidates = _contrast()
    records = [_record(candidate, deviation) for candidate, deviation in zip(candidates, deviations)]
    result = score_signal_contrast(candidates, records, FusionSignalSettings())
    assert result["status"] == status
    for candidate, record, score in zip(candidates, records, result["scores"]):
        expected, count = _mean_nll_in_window(record, candidate, result["windows"], reduce="sum")
        assert (score["summed_nll"], score["events"]) == (expected, count)


@pytest.mark.parametrize("failure", ["missing", "duplicate", "failed", "one_side", "invalid_events"])
def test_incomplete_signal_cannot_win(failure):
    candidates = _contrast()
    records = [_record(candidate) for candidate in candidates]
    if failure == "missing":
        records.pop()
    elif failure == "duplicate":
        records.append(records[0])
    elif failure == "failed":
        records[1]["status"] = 1
    elif failure == "one_side":
        records[1]["position"] = [0] * len(records[1]["position"])
    else:
        records[1]["model_stdv"] = [float("nan")] * len(records[1]["position"])
    assert score_signal_contrast(candidates, records, FusionSignalSettings())["status"] == "abstain"


def _analysis():
    genome, path = _fixture()
    junction = next(path_junctions(path))
    selected = junction[3].exit_position, junction[4].entry_position
    alternative = selected[0] + 1, selected[1]
    contrasts = tuple(JunctionContrast(positions, SimpleNamespace(status="scored"), "", ())
                      for positions in (selected, alternative))
    decision = JunctionDecision(selected, selected, "retained", "test", 0, 0, contrasts, 10)
    correction = PathCorrection(path, path, "retained", "test", (PathJunction("fusion", 0, 0, decision),))
    read = FusionReadAudit("read", "supported_path", "whole_chain_validated", selected=(path,), corrections=(correction,))
    settings = FusionSettings(min_support=1, correction_tolerance=0)
    chains = finalize_chains([path], {"read": read}, settings)
    analysis = FusionAnalysis(chains, (read,), {}, settings, "provided")
    evidence = SimpleNamespace(molecules=(SimpleNamespace(read_id="read", sequence=stitched_template(path, genome)[0]),))
    return analysis, evidence, genome


def test_signal_conflict_removes_full_chain_vote_without_mutating_original():
    analysis, evidence, genome = _analysis()
    result = validate_fusion_signal(analysis, evidence, genome, "unused", backend=lambda read, query, candidates: [
        _record(candidate, 2 if index == 0 else 0) for index, candidate in enumerate(candidates)])
    assert not result.chains
    assert result.reads[0].reason == "signal_junction_conflict"
    assert analysis.reads[0].status == "supported_path"
    assert len(analysis.chains) == 1
    assert len(result.reads[0].signal) == 3


def test_signal_resolves_existing_alternatives_without_discarding_supported_read():
    analysis, evidence, genome = _analysis()
    path = analysis.reads[0].selected[0]
    alternate = _move_path(path, next(path_junctions(path)), (501, 200))
    analysis = replace(analysis, reads=(replace(analysis.reads[0], selected=(path, alternate),
                                               reason="coordinate_ambiguity"),))
    truth = stitched_template(path, genome)[0]
    result = validate_fusion_signal(analysis, evidence, genome, "unused", backend=lambda read, query, candidates: [
        _record(candidate, 0 if candidate.sequence == truth else 2) for candidate in candidates])
    assert result.reads[0].selected == (path,)
    assert result.reads[0].reason == "signal_resolved_coordinates"
    assert result.chains[0].coordinate_status == "signal_resolved"
    assert result.chains[0].full_support_read_ids == ("read",)
    assert any(not audit["retained_alternative"] for audit in result.reads[0].signal)


def test_backend_missing_abstains_and_keeps_sequence_evidence():
    analysis, evidence, genome = _analysis()

    def fail(*args):
        raise ImportError("krill unavailable")

    result = validate_fusion_signal(analysis, evidence, genome, "unused", backend=fail)
    assert result.chains == analysis.chains
    assert result.signal_status == "abstained"
    assert result.reads[0].signal[0]["reason"] == "signal_backend_or_record_error"


def test_hypothesis_limit_abstains_entire_read():
    analysis, evidence, genome = _analysis()
    result = validate_fusion_signal(analysis, evidence, genome, "unused", FusionSignalSettings(max_hypotheses=1),
                                    backend=lambda *args: pytest.fail("must not score truncated hypotheses"))
    assert result.chains == analysis.chains
    assert all(audit["reason"] == "signal_hypothesis_limit" for audit in result.reads[0].signal)


def test_identical_templates_never_resolve_microhomology():
    analysis, evidence, genome = _analysis()
    path = analysis.reads[0].selected[0]
    decision = analysis.reads[0].corrections[0].junctions[0].decision
    contrast = replace(decision.contrasts[1], positions=(501, 201))
    decision = replace(decision, contrasts=(decision.contrasts[0], contrast))
    correction = replace(analysis.reads[0].corrections[0], junctions=(PathJunction("fusion", 0, 0, decision),))
    analysis = replace(analysis, reads=(replace(analysis.reads[0], corrections=(correction,)),))
    genome = dict(genome, A=genome["A"][:500] + genome["B"][200] + genome["A"][501:])
    result = validate_fusion_signal(analysis, evidence, genome, "unused",
                                    backend=lambda *args: pytest.fail("identical templates must abstain"))
    assert result.reads[0].selected == (path,)
    assert result.reads[0].signal[0]["reason"] == "identical_transcript_alternatives"


@pytest.mark.parametrize("settings", [FusionSignalSettings(flank=0), FusionSignalSettings(min_margin=float("nan")),
                                     FusionSignalSettings(min_events=0), FusionSignalSettings(max_hypotheses=0)])
def test_invalid_signal_settings_rejected(settings):
    with pytest.raises(ValueError):
        settings.validate()


def test_backend_uses_real_mappy_offsets_and_caches_reordered_templates(monkeypatch):
    genome, path = _fixture()
    query = stitched_template(path, genome)[0]
    candidates = virtual_candidates([path, _move_path(path, next(path_junctions(path)), (501, 200))], genome)
    calls = []

    def align(signal_path, variants, **kwargs):
        calls.append((signal_path, variants, kwargs))
        return {"read": [dict(variant_label=label, status=0) for label in variants["read"]]}

    monkeypatch.setitem(sys.modules, "krill", SimpleNamespace(align_reads_variants=align))
    monkeypatch.setattr("fin.scoring.krill_aligner.make_krill_aligner", lambda *args, **kwargs: (object(), False))
    backend = KrillFusionBackend("signal.blow5", FusionSignalSettings())
    first = backend("read", query, candidates)
    second = backend("read", query, tuple(reversed(candidates)))
    assert len(calls) == 1
    assert calls[0][0] == "signal.blow5"
    assert calls[0][2]["pore"] == "rna002"
    assert calls[0][2]["start"] == {"read": {"hypothesis_0": 0, "hypothesis_1": 0}}
    assert first == list(reversed(second))


def test_incomplete_sequence_nominations_abstain_before_signal(monkeypatch):
    analysis, evidence, genome = _analysis()
    monkeypatch.setattr("fin.fusion.signal.evaluate_junction", lambda *args, **kwargs:
                        SimpleNamespace(reason="candidate_limit"))
    result = validate_fusion_signal(analysis, evidence, genome, "unused",
                                    backend=lambda *args: pytest.fail("must not score partial nominations"))
    assert result.chains == analysis.chains
    assert all(audit["reason"] == "sequence_nomination_candidate_limit" for audit in result.reads[0].signal)
