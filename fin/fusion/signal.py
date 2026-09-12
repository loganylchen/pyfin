"""Fusion-local eventalign validation using assembly's summed Gaussian NLL.

An artificial coordinate axis allows the existing exon/event projection and
differing-junction windows to represent arbitrary chromosome/strand arms.
These transient candidates never enter assembly discovery, scoring or EM.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math

from fin.candidates.dataclasses import IntronChain, TranscriptCandidate
from fin.fusion.analysis import finalize_chains
from fin.fusion.chains import coordinate_key
from fin.fusion.correction import path_junctions
from fin.fusion.junctions import evaluate_junction, move_junction
from fin.fusion.validation import stitched_template


@dataclass(frozen=True)
class FusionSignalSettings:
    pore: str = "rna002"
    use_gpu: bool = False
    flank: int = 8
    min_margin: float = 1.0
    min_events: int = 3
    max_hypotheses: int = 128

    def validate(self):
        if (not self.pore or min(self.flank, self.min_events, self.max_hypotheses) < 1
                or not math.isfinite(self.min_margin) or self.min_margin <= 0):
            raise ValueError("invalid fusion signal settings")


def virtual_candidates(paths, genome):
    """Map each ordered arm to a fixed, disjoint forward-coordinate interval."""
    if not paths or not paths[0].arms:
        raise ValueError("signal contrast needs nonempty paths")
    shape = tuple((arm.chrom, arm.strand, len(arm.blocks)) for arm in paths[0].arms)
    if any(tuple((arm.chrom, arm.strand, len(arm.blocks)) for arm in path.arms) != shape for path in paths):
        raise ValueError("signal contrast requires identical arm topology")
    axes = []
    offset = 0
    for arm_index in range(len(shape)):
        blocks = [block for path in paths for block in path.arms[arm_index].blocks]
        low = min(block.reference_start for block in blocks)
        high = max(block.reference_end for block in blocks)
        axes.append((low, high, offset))
        offset += high - low + 128
    candidates = []
    for index, path in enumerate(paths):
        exons = []
        for arm, (low, high, offset) in zip(path.arms, axes):
            for block in arm.blocks:
                exons.append((offset + block.reference_start - low, offset + block.reference_end - low)
                             if arm.strand == "+" else
                             (offset + high - block.reference_end, offset + high - block.reference_start))
        if any(previous[1] >= following[0] for previous, following in zip(exons, exons[1:])):
            raise ValueError("invalid virtual exon order")
        sequence, _ = stitched_template(path, genome)
        if not sequence or len(sequence) != sum(end - start for start, end in exons):
            raise ValueError("invalid signal transcript sequence")
        candidates.append(TranscriptCandidate(
            f"hypothesis_{index}", IntronChain(tuple((previous[1], following[0])
                                                     for previous, following in zip(exons, exons[1:]))),
            exons[-1][1], sequence, "fusion_signal", {path.read_id}, "fusion_virtual", "+",
            exons[0][0], exons[-1][1],
        ))
    return tuple(candidates)


class KrillFusionBackend:
    """Own backend state and per-label offsets, with no assembly cache writes."""

    def __init__(self, signal_path, settings):
        self.signal_path = str(signal_path)
        self.settings = settings
        self.aligner = None
        self._read_id = None
        self._records = {}

    def __call__(self, read_id, query, candidates):
        import krill
        import mappy
        from fin.scoring.krill_aligner import make_krill_aligner, krill_thread_count
        from fin.scoring.mappy_score import score_hit

        if self._read_id != read_id:
            self._read_id, self._records = read_id, {}
        if all(candidate.sequence in self._records for candidate in candidates):
            return [dict(self._records[candidate.sequence], variant_label=candidate.candidate_id)
                    for candidate in candidates]
        if self.aligner is None:
            self.aligner, self.gpu = make_krill_aligner(
                krill, self.settings.pore, self.settings.use_gpu,
                hmm_confidence=False, num_thread=krill_thread_count(),
            )
        if self.aligner is None:
            raise RuntimeError("krill aligner unavailable")
        variants, starts = {}, {}
        for candidate in candidates:
            aligner = mappy.Aligner(seq=candidate.sequence, preset="map-ont", best_n=32)
            if aligner is None:
                raise RuntimeError("signal template index unavailable")
            hits = []
            for hit in aligner.map(query):
                hits.append(hit)
                if len(hits) >= 32:
                    raise RuntimeError("signal template mapping limit")
            scored = [(score_hit(hit), hit) for hit in hits if hit.strand == 1 and hit.cigar
                      and (hit.q_en - hit.q_st) >= 0.9 * len(query)]
            scored = [(score, hit) for score, hit in scored if score is not None]
            if not scored:
                raise RuntimeError("signal template lacks complete query alignment")
            best = max(score for score, _ in scored)
            offsets = {(hit.r_st, hit.r_en) for score, hit in scored if score == best}
            if len(offsets) != 1:
                raise RuntimeError("ambiguous signal template alignment")
            start, end = offsets.pop()
            variants[candidate.candidate_id] = candidate.sequence[start:end]
            starts[candidate.candidate_id] = start
        output = krill.align_reads_variants(
            self.signal_path, {read_id: variants}, pore=self.settings.pore, use_gpu=self.gpu,
            num_thread=krill_thread_count(), aligner=self.aligner, start={read_id: starts},
        )
        records = output.get(read_id, ())
        by_label = {record.get("variant_label"): record for record in records}
        if len(by_label) == len(records) and set(by_label) == set(variants):
            for candidate in candidates:
                self._records[candidate.sequence] = by_label[candidate.candidate_id]
        return records


def score_signal_contrast(candidates, records, settings):
    """All hypotheses must cover every differing junction with finite events."""
    import numpy as np
    from fin.scoring.m2_junction_nll import (
        _mean_nll_in_window, _zrecords, diff_junction_windows, tx2genome_array,
    )

    windows = diff_junction_windows(list(candidates), flank=settings.flank)
    audit = dict(status="abstain", reason="no_differing_junction", scores=[], windows=windows,
                 metric="summed_llr", window_coordinate_system="virtual_forward_arm_axis")
    if not windows:
        return audit
    labels = [record.get("variant_label") for record in records]
    if len(labels) != len(set(labels)) or set(labels) != {candidate.candidate_id for candidate in candidates}:
        return dict(audit, reason="incomplete_or_duplicate_signal_hypotheses")
    by_label = dict(zip(labels, records))
    variable = [index for index in range(len(candidates[0].intron_chain.introns))
                if len({candidate.intron_chain.introns[index] for candidate in candidates}) > 1]
    for candidate in candidates:
        record = by_label[candidate.candidate_id]
        if record.get("status", -1) != 0:
            return dict(audit, reason="eventalign_failed_hypothesis")
        zscores, deviations = _zrecords(record)
        positions = tx2genome_array(candidate, np.asarray(record["position"], dtype=np.int64))
        finite = np.isfinite(zscores) & np.isfinite(deviations) & (positions >= 0)
        covered = positions[finite]
        for index in variable:
            donors, acceptors = zip(*(other.intron_chain.introns[index] for other in candidates))
            left = (min(donors) - settings.flank, min(donors))
            right = (max(acceptors), max(acceptors) + settings.flank)
            if not np.any(covered <= left[0]) or not np.any(covered >= right[1]):
                return dict(audit, reason="insufficient_shared_flank_events")
            for start, end in (left, right):
                if np.count_nonzero((covered >= start) & (covered < end)) < settings.min_events:
                    return dict(audit, reason="insufficient_shared_flank_events")
        score, count = _mean_nll_in_window(record, candidate, windows, reduce="sum")
        if count < settings.min_events or not math.isfinite(score):
            return dict(audit, reason="insufficient_finite_events")
        audit["scores"].append(dict(label=candidate.candidate_id, summed_nll=score, events=count))
    ranked = sorted(audit["scores"], key=lambda item: (item["summed_nll"], item["label"]))
    margin = ranked[1]["summed_nll"] - ranked[0]["summed_nll"]
    winner = ranked[0]["label"]
    if margin < settings.min_margin:
        return dict(audit, reason="signal_tie", margin=margin)
    status = "supported" if winner == candidates[0].candidate_id else "conflict"
    return dict(audit, status=status, reason="summed_llr_contrast", winner=winner, margin=margin)


def _move_path(path, junction, positions):
    kind, arm_index, _, left, right = junction
    left, right = move_junction(left, right, positions)
    arms = list(path.arms)
    if kind == "intron":
        arms[arm_index] = replace(arms[arm_index], blocks=left.blocks + right.blocks)
    else:
        arms[arm_index], arms[arm_index + 1] = left, right
    return replace(path, arms=tuple(arms))


def validate_fusion_signal(analysis, evidence, genome, signal_path, settings=None, *, backend=None):
    """Validate sequence-selected paths; conflicts veto, abstentions never vote.

    Unlike assembly's AS-tie gate, this tests junction coordinate alternatives
    around a sequence-supported fusion. It reuses the scorer, not that gate.
    Signal never promotes an unvalidated sequence path or breaks identical-
    transcript microhomology ties.
    """
    settings = settings or FusionSignalSettings()
    settings.validate()
    backend = backend if backend is not None else KrillFusionBackend(signal_path, settings)
    molecules = {molecule.read_id: molecule for molecule in evidence.molecules}
    reads = []
    for read in analysis.reads:
        audits = []
        used = 0
        if read.status != "supported_path":
            reads.append(read)
            continue
        for alternative_index, path in enumerate(read.selected):
            for junction in path_junctions(path):
                kind, arm_index, block_index, left, right = junction
                audit = dict(alternative_index=alternative_index, kind=kind, arm_index=arm_index,
                             block_index=block_index, path_coordinates=coordinate_key(path),
                             status="abstain", reason="no_alternative")
                positions = {(left.exit_position, right.entry_position)}
                for correction in read.corrections:
                    if tuple((arm.chrom, arm.strand, len(arm.blocks)) for arm in correction.observed.arms) != tuple(
                            (arm.chrom, arm.strand, len(arm.blocks)) for arm in path.arms):
                        continue
                    for item in correction.junctions:
                        if (item.kind, item.arm_index, item.block_index) == (kind, arm_index, block_index):
                            positions.update(contrast.positions for contrast in item.decision.contrasts
                                             if contrast.alignment.status == "scored")
                decision = evaluate_junction(
                    molecules[read.read_id].sequence, left, right, genome,
                    tolerance=analysis.settings.correction_tolerance, flank=analysis.settings.junction_flank,
                    min_anchor=analysis.settings.junction_anchor, min_identity=analysis.settings.min_identity,
                    internal_intron=kind == "intron", max_candidates=analysis.settings.max_junction_candidates,
                )
                if decision.reason in {"candidate_limit", "incomparable_hypotheses", "invalid_query_window",
                                       "invalid_observed_template", "no_anchored_hypothesis", "baseline_not_anchored"}:
                    audits.append(dict(audit, reason="sequence_nomination_" + decision.reason))
                    continue
                positions.update(contrast.positions for contrast in decision.contrasts
                                 if contrast.alignment.status == "scored")
                baseline = left.exit_position, right.entry_position
                ordered = [baseline] + sorted(positions - {baseline})
                audit["positions"] = ordered
                used += len(ordered)
                if used > settings.max_hypotheses:
                    audits.append(dict(audit, reason="signal_hypothesis_limit"))
                    continue
                if len(ordered) < 2:
                    audits.append(audit)
                    continue
                try:
                    candidates = virtual_candidates([_move_path(path, junction, position) for position in ordered], genome)
                    if len({candidate.sequence for candidate in candidates}) != len(candidates):
                        audits.append(dict(audit, reason="identical_transcript_alternatives"))
                        continue
                    records = backend(read.read_id, molecules[read.read_id].sequence, candidates)
                    audit.update(score_signal_contrast(candidates, records, settings))
                    audit["effective_gpu"] = getattr(backend, "gpu", None)
                except Exception as error:
                    audit.update(reason="signal_backend_or_record_error", error=f"{type(error).__name__}: {error}")
                audits.append(audit)
        if used > settings.max_hypotheses:
            audits = [dict(audit, status="abstain", reason="signal_hypothesis_limit") for audit in audits]
        rejected = {audit["alternative_index"] for audit in audits if audit["status"] == "conflict"}
        selected = tuple(path for index, path in enumerate(read.selected) if index not in rejected)
        audits = [dict(audit, retained_alternative=audit["alternative_index"] not in rejected) for audit in audits]
        reason = read.reason
        if not selected:
            reason = "signal_junction_conflict"
        elif rejected:
            reason = "signal_resolved_coordinates" if len({coordinate_key(path) for path in selected}) == 1 else "coordinate_ambiguity"
        reads.append(replace(read, signal=tuple(audits), selected=selected,
                             status=read.status if selected else "unresolved", reason=reason,
                             assembly_junctions=tuple(dict(audit, retained_alternative=audit["alternative_index"] not in rejected)
                                                      for audit in read.assembly_junctions)))
    by_read = {read.read_id: read for read in reads}
    chains = finalize_chains([read.selected[0] for read in reads if read.status == "supported_path"],
                             by_read, analysis.settings)
    statuses = [audit["status"] for read in reads for audit in read.signal]
    status = ("no_applicable_reads" if not statuses else "abstained" if set(statuses) == {"abstain"}
              else "scored" if "abstain" not in statuses else "partially_scored")
    return replace(analysis, reads=tuple(reads), chains=chains, signal_status=status)
