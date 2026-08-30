"""Assignment layer: read-to-candidate assignment for quant_mode="m2_em".

Extracted from the runner's ``_quant_m2_em`` so the M1/M2/EM assignment orchestration lives in
its own module (intervals, clustering, EM, quantification are already standalone modules; this is
the last piece that was glued into the runner).

- Step A: ``tie_nll`` — the per-read best-AS tie junction-NLL scoring (moved verbatim from
  ``PipelineRunner._tie_nll``; ``self.config`` -> ``config``).
- Step B: ``Assigner.assign`` — the pure assignment pipeline (mappy AS matrix + kept-read set ->
  ``tie_nll`` -> d_tx skeleton -> zero-coherence assignment -> optional krill tie-break ->
  ``quantify_transcripts`` + ``max_R`` stamp) returning an immutable ``QuantOutput``. NO candidate
  dropping happens here; every structural gate/fold stays in the selection layer.
  ``quantify_transcripts`` is stamped here (right after assignment) rather than
  after the drop cascade: byte-identical, because the intervening gates only accumulate a
  ``drop_cols`` column-index set and never mutate ``R`` / ``hard_assignments`` / ``cand_list``.

Byte-identical extraction. The runner keeps thin ``_tie_nll`` / ``_eff_lengths`` methods (passed into
``assign`` at call time) so the existing ``patch.object(PipelineRunner, ...)`` test seams survive; the
``em_with_coherence`` / ``quantify_transcripts`` seams move to this module (tests re-point here).
"""
from __future__ import annotations

import logging
import math
from collections import defaultdict
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from fin.analysis.assignments import em_with_coherence
from fin.analysis.quantification import quantify_transcripts
from fin.scoring.krill_tiebreak import krill_tiebreak
from fin.scoring.m2_junction_nll import MISSING

logger = logging.getLogger(__name__)

# Junction-window flank (bp) for the M2 tie discrimination window (matches runner).
_M2_EM_FLANK = 2
# d_tx MISSING-cell pad added to (hi-lo) so an unscored tie cell sits just beyond the
# worst scored margin (matches the runner's original constant).
_M2_EM_PAD = 1.0
_EVENTALIGN_BASES = frozenset("ACGTUacgtu")


def _eventalignable_variants(variants: Dict[str, str]) -> bool:
    """True when every candidate slice is non-empty and krill-compatible."""
    return bool(variants) and all(
        seq and all(base in _EVENTALIGN_BASES for base in seq)
        for seq in variants.values()
    )


def _raw_tie_sets(raw: np.ndarray, n_c: int) -> Tuple[Dict[int, List[int]], int]:
    """Return exact best-AS tie sets without invoking the signal backend."""
    ties: Dict[int, List[int]] = {}
    n_ties = 0
    for i in range(raw.shape[0]):
        row = raw[i]
        finite = np.isfinite(row)
        if not finite.any():
            continue
        best_as = float(row[finite].max())
        tie = [j for j in range(n_c) if finite[j] and row[j] >= best_as - 1e-9]
        ties[i] = tie
        if len(tie) >= 2:
            n_ties += 1
    return ties, n_ties


def _m2_score_region(config, tied_cands) -> Tuple[Optional[set], List[Tuple[int, int]]]:
    """Select the production scoring footprint for one M1 tie."""
    from fin.scoring.m2_junction_nll import (
        class_junction_window_set,
        diff_junction_windows,
    )

    metric = getattr(config, "m2_metric", "mean")
    if metric in ("summed_llr", "sqrt_count_mean_llr"):
        # Containment/cassette contrasts have asymmetric candidate-private event
        # populations. The validated summed contrast is the same-intron-count
        # wobble niche; other structural classes abstain. The experimental
        # sqrt-count variant shares the footprint but rescales by effective
        # event count, so unequal per-candidate coverage cannot masquerade as
        # signal strength.
        intron_counts = {len(c.intron_chain.introns) for c in tied_cands}
        if len(intron_counts) != 1:
            return None, []
        flank = int(getattr(config, "m2_summed_llr_flank", 6))
        return None, diff_junction_windows(tied_cands, flank=flank)
    if metric == "mean":
        gset = class_junction_window_set(
            tied_cands,
            flank=_M2_EM_FLANK,
            k=config.m2_tiebreak_junction_k,
        )
        return gset, []
    return None, []


def _record_abstention(stats, rid, tie, cand_list, reason):
    """Machine-readable abstention record (capture mode only)."""
    stats[f"abstain_{reason}"] += 1
    records = stats.get("_records")
    if isinstance(records, list):
        records.append({
            "kind": "abstention",
            "read_id": str(rid),
            "reason": reason,
            "tie_candidates": [
                getattr(cand_list[j], "candidate_id", str(j)) for j in tie
            ],
        })


def _finalize_contrast(metric, nlls, evs, stats, *, read_id=None,
                       cand_ids=None, coverage=None):
    """Record contrast observability and apply the sqrt-count rescale.

    For the summed/sqrt-count footprint the two best hypotheses may rest on
    unequal effective event counts; the stats make that visible per run. The
    experimental ``sqrt_count_mean_llr`` metric multiplies each candidate's MEAN
    window NLL by sqrt(min event count over scored hypotheses), so the margin
    behaves like a z-scale difference of means instead of an undivided sum
    whose magnitude tracks event count. Ordering within a read is preserved
    (positive common scale); ``summed_llr`` values pass through unchanged.
    """
    if metric not in ("summed_llr", "sqrt_count_mean_llr"):
        return nlls
    ordered = sorted(nlls.items(), key=lambda kv: kv[1])
    (j1, v1), (j2, v2) = ordered[0], ordered[1]
    e1, e2 = evs.get(j1, 0), evs.get(j2, 0)
    stats["decided"] += 1
    stats["same_event_count"] += 1 if e1 == e2 else 0
    stats["diff_event_count"] += 0 if e1 == e2 else 1
    stats["min_ev_sum"] += min(e1, e2)
    if metric == "sqrt_count_mean_llr":
        scale = math.sqrt(max(min(evs[j] for j in nlls), 1))
        nlls = {j: v * scale for j, v in nlls.items()}
        ordered = sorted(nlls.items(), key=lambda kv: kv[1])
    margin = ordered[1][1] - ordered[0][1]
    stats["margin_sum"] += margin
    records = stats.get("_records")
    if isinstance(records, list):
        # Per-comparison experimental record. v1/v2 are the metric's raw
        # per-hypothesis values BEFORE any sqrt-count rescale; mean and sum
        # NLL are both materialized (same event set, so mean = sum / n_ev).
        # The same-intron-count guard already held or this contrast would
        # not exist.
        if metric == "summed_llr":
            sum1, sum2 = float(v1), float(v2)
            mean1 = sum1 / e1 if e1 else float("nan")
            mean2 = sum2 / e2 if e2 else float("nan")
        else:  # sqrt_count_mean_llr: raw values are means
            mean1, mean2 = float(v1), float(v2)
            sum1, sum2 = mean1 * e1, mean2 * e2
        j_best, j_runner = int(ordered[0][0]), int(ordered[1][0])
        records.append({
            "kind": "comparison",
            "read_id": None if read_id is None else str(read_id),
            "winner_col": j_best,
            "winner_id": (cand_ids or {}).get(j_best),
            "runner_id": (cand_ids or {}).get(j_runner),
            "margin": round(float(margin), 6),
            "nll_mean_best": round(mean1, 6),
            "nll_mean_runner": round(mean2, 6),
            "nll_sum_best": round(sum1, 6),
            "nll_sum_runner": round(sum2, 6),
            "nll_mean_delta": round(mean2 - mean1, 6),
            "nll_sum_delta": round(sum2 - sum1, 6),
            "ev_best": int(e1),
            "ev_runner": int(e2),
            "n_scored": len(nlls),
            "coverage": coverage,
            "same_intron_count": True,
        })
    return nlls


def _log_contrast_stats(metric, stats, config=None):
    """Observability for the tight-window metrics.

    Always logs one compact line per interval. When
    ``config.m2_contrast_stats_jsonl`` is enabled (experimental), also
    appends a machine-readable JSON record under
    ``<work_dir>/m2_contrasts/`` (one file per interval batch, pid+counter
    named, aggregate afterwards) so calibration/reliability analysis does not
    have to parse log text.
    """
    if metric not in ("summed_llr", "sqrt_count_mean_llr"):
        return
    decided = int(stats.get("decided", 0))
    # ANY abstention reason must produce output: a locus that abstained
    # entirely (no krill backend, no mappy alignment, missing batch rows) is
    # exactly the case calibration analysis must be able to see.
    if not decided and not any(
        k.startswith("abstain_") and v for k, v in stats.items()
    ):
        return
    logger.info(
        "M2 %s contrasts: decided=%d abstained_lt2=%d abstain_window=%d "
        "abstain_count=%d abstain_payload=%d same_events=%d diff_events=%d "
        "mean_min_ev=%.1f mean_margin=%.2f",
        metric, decided, int(stats.get("abstain_lt2_scored", 0)),
        int(stats.get("abstain_no_window", 0)),
        int(stats.get("abstain_unequal_intron_count", 0)),
        int(stats.get("abstain_invalid_payload", 0)),
        int(stats.get("same_event_count", 0)),
        int(stats.get("diff_event_count", 0)),
        stats.get("min_ev_sum", 0.0) / decided if decided else 0.0,
        stats.get("margin_sum", 0.0) / decided if decided else 0.0,
    )
    if config is not None and getattr(config, "m2_contrast_stats_jsonl", False):
        try:
            import json
            import os
            import time
            from pathlib import Path

            out_dir = Path(config.work_dir) / "m2_contrasts"
            out_dir.mkdir(parents=True, exist_ok=True)
            lines = []
            for rec in stats.get("_records") or []:
                rec = dict(rec)
                rec.setdefault("kind", "comparison")
                rec["metric"] = metric
                lines.append(json.dumps(rec, sort_keys=True))
            aggregate = {
                "kind": "aggregate", "metric": metric,
                **{k: float(v) for k, v in sorted(stats.items())
                   if not k.startswith("_")},
            }
            lines.append(json.dumps(aggregate, sort_keys=True))
            path = out_dir / f"{os.getpid()}_{time.monotonic_ns()}.jsonl"
            path.write_text("\n".join(lines) + "\n")
        except Exception:  # experimental diagnostics never break a run
            logger.exception("m2 contrast JSONL emission failed")


def tie_nll(config, kept_read_ids, read_seqs, cand_list, aligners, raw):
    """Per-read best-AS tie junction-NLL (the m2_resolve_tie niche).

    For each read whose raw mappy AS is simultaneously best across >=2 candidates, build the
    configured class junction-discrimination window and score each tied candidate with either
    legacy wide-window mean NLL or tight-window summed LLR. Only tied, window-scorable cells are
    touched -- never the dense read×candidate matrix. A valid refinement requires at least two
    scored hypotheses; technical missingness therefore causes abstention, not a hard winner.

    Returns (nlls_by_read {i: {j: nll}}, ties_by_read {i: [j,...]}, n_ties, n_refined,
    cover_by_read {i: bool}). ``cover_by_read`` is populated only when the diff-region coverage
    gate (config.m2_diff_cover_gate) is ON; otherwise it is empty.
    """
    metric = getattr(config, "m2_metric", "mean")
    gate = bool(getattr(config, "m2_diff_cover_gate", False))
    n_c = len(cand_list)
    if metric == "off":
        ties, n_ties = _raw_tie_sets(raw, n_c)
        return {}, ties, n_ties, 0, {}

    import krill

    from fin.scoring.krill_aligner import (
        krill_thread_count,
        make_krill_aligner,
    )
    from fin.scoring.m2_junction_nll import (
        _mean_nll_in_gset,
        _mean_nll_in_window,
        event_genomic_positions,
        read_cand_mean_nll,
        read_straddles,
        wobble_diff_spans,
    )
    from fin.scoring.mappy_score import score_hit

    num_threads = krill_thread_count()
    nlls_by_read: Dict[int, Dict[int, float]] = {}
    ties_by_read: Dict[int, List[int]] = {}
    cover_by_read: Dict[int, bool] = {}
    stats: Dict[str, object] = defaultdict(float)
    if getattr(config, "m2_contrast_stats_jsonl", False):
        stats["_records"] = []
    krill_aligner, eff_gpu = make_krill_aligner(
        krill, config.krill_pore, config.use_gpu,
        hmm_confidence=False, num_thread=num_threads,
    )
    if krill_aligner is None:
        # No krill signal backend. Gate OFF: preserve the legacy behavior
        # (empty ties -> _quant_m2_em leaves rows MISSING, byte-identical).
        # Gate ON: still expose the AS tie sets (derivable from ``raw`` alone,
        # no signal needed) so the gate's unique-best / drop logic is
        # well-defined; without them every row would be all-MISSING and the
        # EM would normalize each read to a spurious uniform vote.
        if gate:
            for i in range(len(kept_read_ids)):
                row = raw[i]
                finite = np.isfinite(row)
                if not finite.any():
                    continue
                best_as = float(row[finite].max())
                ties_by_read[i] = [
                    j for j in range(n_c) if finite[j] and row[j] >= best_as - 1e-9
                ]
        for i, tie in ties_by_read.items():
            if len(tie) >= 2:
                _record_abstention(
                    stats, kept_read_ids[i], tie, cand_list, "no_backend",
                )
        _log_contrast_stats(metric, stats, config)
        return nlls_by_read, ties_by_read, 0, 0, cover_by_read

    pore = config.krill_pore
    sig_path = config.signal_path

    # --- Pass 1: tie detection + per-(read, candidate) mappy best-hit slices.
    # Each tied cell is reduced to the SAME (sliced sequence, start offset)
    # the singular read_cand_mean_nll would have eventaligned, so the batched
    # eventalign below is byte-identical to the per-call path -- only the
    # Python/GIL/per-call overhead is removed.
    gset_by_read: Dict[int, set] = {}
    windows_by_read: Dict[int, List[Tuple[int, int]]] = {}
    spans_by_read: Dict[int, List[Tuple[int, int]]] = {}  # diff-cover gate spans
    reads_variants: Dict[str, Dict[str, str]] = {}  # rid -> {cand_id: seq}
    starts: Dict[str, Dict[str, int]] = {}          # rid -> {cand_id: r_st}
    cidj_by_read: Dict[str, Dict[str, int]] = {}     # rid -> {cand_id: col j}
    n_ties = 0
    n_invalid_payloads = 0
    for i, rid in enumerate(kept_read_ids):
        row = raw[i]
        finite = np.isfinite(row)
        if not finite.any():
            continue
        best_as = float(row[finite].max())
        tie = [j for j in range(n_c) if finite[j] and row[j] >= best_as - 1e-9]
        ties_by_read[i] = tie
        if len(tie) < 2:
            continue
        n_ties += 1
        tied_cands = [cand_list[j] for j in tie]
        gset, windows = _m2_score_region(config, tied_cands)
        if not gset and not windows:
            _record_abstention(
                stats, rid, tie, cand_list,
                "unequal_intron_count"
                if len({len(c.intron_chain.introns) for c in tied_cands}) != 1
                else "no_window",
            )
            continue
        if gate:
            spans_by_read[i] = wobble_diff_spans(
                tied_cands, flank=_M2_EM_FLANK
            )
        seq = read_seqs[rid]
        per_seqs: Dict[str, str] = {}
        per_starts: Dict[str, int] = {}
        per_cidj: Dict[str, int] = {}
        for j in tie:
            aln = aligners[j]
            if aln is None:
                continue
            best_hit = None
            best_hit_as = None
            for h in aln.map(seq):
                v = score_hit(h)
                if v is None:
                    continue
                if best_hit_as is None or v > best_hit_as:
                    best_hit_as, best_hit = v, h
            if best_hit is None:
                continue
            cand = cand_list[j]
            per_seqs[cand.candidate_id] = cand.sequence[best_hit.r_st:best_hit.r_en]
            per_starts[cand.candidate_id] = best_hit.r_st
            per_cidj[cand.candidate_id] = j
        if not per_seqs:
            _record_abstention(stats, rid, tie, cand_list, "no_alignment")
            continue
        if not _eventalignable_variants(per_seqs):
            # One invalid hypothesis must not poison eventalign for the entire
            # locus. Preserve this read's complete M1 tie and abstain from M2;
            # replacing ambiguous bases or scoring only a subset would turn
            # technical missingness into biological evidence.
            n_invalid_payloads += 1
            _record_abstention(stats, rid, tie, cand_list, "invalid_payload")
            continue
        if gset:
            gset_by_read[i] = gset
        if windows:
            windows_by_read[i] = windows
        reads_variants[rid] = per_seqs
        starts[rid] = per_starts
        cidj_by_read[rid] = per_cidj

    if n_invalid_payloads:
        logger.warning(
            "M2 abstained on %d tied reads with unsupported bases in candidate "
            "eventalign sequences; remaining reads stay batched",
            n_invalid_payloads,
        )
    if not reads_variants:
        # Every tied read abstained in pass 1; still emit the record set so a
        # fully-abstaining locus is visible in the calibration data.
        _log_contrast_stats(metric, stats, config)
        return nlls_by_read, ties_by_read, n_ties, 0, cover_by_read

    # --- Pass 2: ONE batched eventalign over every tied pair (per-pair start
    # via the {label: offset} form); per-read singular fallback on absence or
    # batch failure. Both branches use the configured mean/summed reducer.
    rid_to_i = {rid: i for i, rid in enumerate(kept_read_ids)}
    n_refined = 0
    use_batch = hasattr(krill, "align_reads_variants")
    batch_out = None
    if use_batch:
        try:
            batch_out = krill.align_reads_variants(
                sig_path, reads_variants, pore=pore,
                use_gpu=eff_gpu, num_thread=num_threads,
                aligner=krill_aligner, start=starts,
            )
        except Exception as exc:  # noqa: BLE001 - krill raises broad errors
            logger.warning("M2 batch eventalign failed (%s); per-read fallback", exc)
            use_batch = False

    if use_batch and batch_out is not None:
        for rid in cidj_by_read:
            if rid not in batch_out:
                i_missing = rid_to_i.get(rid)
                _record_abstention(
                    stats, rid,
                    ties_by_read.get(i_missing, []) if i_missing is not None else [],
                    cand_list, "batch_output_missing",
                )
        for rid, res_list in batch_out.items():
            i = rid_to_i.get(rid)
            if i is None:
                continue
            gset = gset_by_read.get(i)
            windows = windows_by_read.get(i, [])
            if metric == "mean" and gset is None:
                continue
            if metric in ("summed_llr", "sqrt_count_mean_llr") and not windows:
                _record_abstention(
                    stats, rid, ties_by_read.get(i, []), cand_list, "no_window",
                )
                continue
            by_label = {x.get("variant_label"): x for x in res_list}
            nlls: Dict[int, float] = {}
            evs: Dict[int, int] = {}
            spans = spans_by_read.get(i) if gate else None
            span_cov = [False] * len(spans) if spans else None
            for cid, j in cidj_by_read.get(rid, {}).items():
                res = by_label.get(cid)
                if res is None or res.get("status", -1) != 0:
                    continue
                if metric == "summed_llr":
                    nll, n_ev = _mean_nll_in_window(
                        res, cand_list[j], windows, reduce="sum"
                    )
                elif metric == "sqrt_count_mean_llr":
                    nll, n_ev = _mean_nll_in_window(
                        res, cand_list[j], windows, reduce="mean"
                    )
                else:
                    nll, n_ev = _mean_nll_in_gset(res, cand_list[j], gset)
                if n_ev > 0 and np.isfinite(nll):
                    nlls[j] = float(nll)
                    evs[j] = int(n_ev)
                if span_cov is not None:
                    egp = event_genomic_positions(res, cand_list[j])
                    if egp:
                        for s_idx, (lo, hi) in enumerate(spans):
                            if not span_cov[s_idx] and read_straddles(egp, lo, hi):
                                span_cov[s_idx] = True
            if len(nlls) >= 2:
                cover_val = (all(span_cov) if span_cov else True) if gate else None
                nlls = _finalize_contrast(
                    metric, nlls, evs, stats, read_id=rid,
                    cand_ids={j: cand_list[j].candidate_id for j in nlls},
                    coverage=cover_val,
                )
                n_refined += 1
                nlls_by_read[i] = nlls
                if gate and cover_val is not None:
                    # covered iff every wobbling diff span is straddled by the
                    # read on >=1 candidate (empty spans -> vacuously covered).
                    cover_by_read[i] = cover_val
            elif evs or cidj_by_read.get(rid):
                _record_abstention(
                    stats, rid, ties_by_read.get(i, []), cand_list,
                    "lt2_scored",
                )
        _log_contrast_stats(metric, stats, config)
        return nlls_by_read, ties_by_read, n_ties, n_refined, cover_by_read

    for rid, per_cidj in cidj_by_read.items():
        i = rid_to_i.get(rid)
        if i is None:
            continue
        gset = gset_by_read.get(i)
        windows = windows_by_read.get(i, [])
        seq = read_seqs[rid]
        nlls = {}
        evs = {}
        for j in per_cidj.values():
            nll, n_ev = read_cand_mean_nll(
                rid, seq, cand_list[j], windows, krill_aligner, aligners[j],
                sig_path, pore, gset=gset, use_gpu=eff_gpu,
                num_thread=num_threads,
                reduce="sum" if metric == "summed_llr" else "mean",
            )
            if n_ev > 0 and np.isfinite(nll):
                nlls[j] = float(nll)
                evs[j] = int(n_ev)
        if len(nlls) >= 2:
            nlls = _finalize_contrast(
                metric, nlls, evs, stats, read_id=rid,
                cand_ids={j: cand_list[j].candidate_id for j in nlls},
                coverage=None,
            )
            n_refined += 1
            nlls_by_read[i] = nlls
        elif evs or per_cidj:
            _record_abstention(
                stats, rid, ties_by_read.get(i, []), cand_list, "lt2_scored",
            )
    # Per-read fallback path does not expose eventalign positions, so the
    # coverage map is left empty here, so fallback contrasts never seed the
    # covered-read redistribution prior.
    _log_contrast_stats(metric, stats, config)
    return nlls_by_read, ties_by_read, n_ties, n_refined, cover_by_read


@dataclass
class QuantOutput:
    """Immutable result of the m2_em assignment pipeline (pre-selection).

    Carries everything the runner's selection tail needs: the read x candidate
    responsibility matrix ``R`` and ``hard_assignments``, the column-aligned
    pre-selection ``quant_results`` (``max_R`` already stamped), the candidate
    and kept-read axes, M2 tie evidence, and log counters. Consumers treat it as
    read-only. ``quant_results`` is column-aligned to ``R`` (enumerate index ==
    R column) at construction; selection must not reorder it in place.
    """

    R: np.ndarray
    hard_assignments: np.ndarray
    quant_results: List  # List[QuantResult]
    cand_list: List
    kept_read_ids: List[str]
    ties_by_read: Dict[int, List[int]]
    nlls_by_read: Dict[int, Dict[int, float]]
    cover_by_read: Dict[int, bool]
    read_seqs: Dict[str, str]
    n_ties: int
    n_refined: int


class Assigner:
    """m2_em read-to-candidate assignment (no dropping). See module docstring.

    Holds the pipeline config. ``assign`` receives the runner's ``_tie_nll`` / ``_eff_lengths`` bound
    methods at call time (``tie_nll_fn`` / ``eff_lengths_fn``) so the ``patch.object(PipelineRunner,
    ...)`` test seams keep intercepting them.
    """

    def __init__(self, config):
        self.config = config

    def assign(
        self,
        candidate_set,
        read_ids: List[str],
        *,
        tie_nll_fn: Callable,
        eff_lengths_fn: Callable,
    ) -> Optional[QuantOutput]:
        """Run the m2_em assignment pipeline; return a QuantOutput or None (empty interval)."""
        import mappy

        from fin.scoring.mappy_preset import get_m1_preset
        from fin.scoring.mappy_score import score_hit

        cand_list = list(candidate_set.candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        reads_iter = [(rid, read_sequences.get(rid, "")) for rid in read_ids]
        reads_iter = [(rid, seq) for rid, seq in reads_iter if seq]
        if not reads_iter:
            return None

        read_seqs = {rid: seq for rid, seq in reads_iter}
        n_c = len(cand_list)
        max_iter_em = (
            self.config.em_max_iter_override
            if self.config.em_max_iter_override is not None
            else self.config.em_max_iter
        )

        # --- Per-candidate mappy aligners + raw AS matrix (tie detection). This
        #     single AS pass also yields the kept-read set, so the previously
        #     separate mappy_multimap_responsibilities pass (which recomputed the
        #     identical AS only to derive kept reads, never feeding the EM) is
        #     removed: ~40% of the per-interval cost on dense loci. ---
        preset = get_m1_preset()
        aligners = [
            mappy.Aligner(seq=c.sequence, preset=preset) if c.sequence else None
            for c in cand_list
        ]
        all_ids = [rid for rid, _ in reads_iter]
        raw_all = np.full((len(all_ids), n_c), -np.inf, dtype=np.float64)
        for i, rid in enumerate(all_ids):
            seq = read_seqs[rid]
            for j, aln in enumerate(aligners):
                if aln is None:
                    continue
                best = None
                for h in aln.map(seq):
                    v = score_hit(h)
                    if v is None:
                        continue
                    if best is None or v > best:
                        best = v
                if best is not None:
                    raw_all[i, j] = best

        # Keep reads with >=1 candidate at AS>0 (matches the removed
        # mappy_multimap_responsibilities keep rule: row.sum()>0 over best>=0 AS).
        # Mirror its MAPPY_R1_MIN_AS env knob so the kept set stays identical when
        # that R1 tuning threshold is set.
        import os as _os
        _min_as = float(_os.environ.get("MAPPY_R1_MIN_AS", "0") or "0")
        keep_rows = ((raw_all > 0.0) & (raw_all >= _min_as)).any(axis=1)
        if not keep_rows.any():
            return None
        kept_read_ids = [rid for rid, k in zip(all_ids, keep_rows) if k]
        raw = raw_all[keep_rows]
        n_r = len(kept_read_ids)

        # --- Junction NLL on the per-read AS-tie cells only. ---
        nlls_by_read, ties_by_read, n_ties, n_refined, cover_by_read = tie_nll_fn(
            kept_read_ids, read_seqs, cand_list, aligners, raw
        )

        # --- d_tx skeleton. Gate OFF: the pure soft NLL-graded skeleton (M1/AS
        #     picks the tie set; junction-NLL is the only graded distance; cells
        #     outside the tie stay MISSING). Gate ON: the diff-region coverage
        #     decision matrix with proportional redistribution of ambiguous reads
        #     (no read is ever dropped -> recall-safe). ---
        gate = bool(getattr(self.config, "m2_diff_cover_gate", False))
        d_tx = np.full((n_r, n_c), MISSING, dtype=np.float64)
        if not gate:
            for i in range(n_r):
                tie = ties_by_read.get(i)
                if not tie:
                    continue
                if len(tie) < 2:
                    d_tx[i, tie[0]] = 0.0
                    continue
                nlls = nlls_by_read.get(i)
                if not nlls or len(nlls) < 2:
                    # A contrast needs two scored hypotheses. Technical
                    # missingness abstains to the flat M1 tie.
                    for j in tie:
                        d_tx[i, j] = 0.0
                    continue
                lo = min(nlls.values())
                hi = max(nlls.values())
                miss = (hi - lo) + _M2_EM_PAD
                for j in tie:
                    d_tx[i, j] = (nlls[j] - lo) if j in nlls else miss
            small = (d_tx > 0) & (d_tx < 1e5)
            diffs = d_tx[small]
            sigma2 = max(float(np.median(diffs)) if diffs.size else 1.0, 1e-3)
            d_tx[small] = d_tx[small] / sigma2
        else:
            # Diff-region coverage gate. margin = runner-up NLL - best NLL.
            #
            # Pass A: classify each tied read and accumulate a per-candidate prior
            # (covered_vote) from ONLY the covered + distinguishing reads -- the
            # gold-standard evidence of the locus isoform ratio. Ambiguous reads
            # (no junction signal, or covered-but-indistinguishable that did not
            # straddle every diff span) are deferred to Pass B.
            metric = getattr(self.config, "m2_metric", "mean")
            if metric in ("summed_llr", "sqrt_count_mean_llr"):
                margin_thr = float(
                    getattr(self.config, "m2_summed_llr_margin", 2.0)
                )
            else:
                margin_thr = float(
                    getattr(self.config, "m2_diff_cover_margin", 0.0)
                )
            covered_vote = np.zeros(n_c, dtype=np.float64)
            fuzzy: List[int] = []
            for i in range(n_r):
                tie = ties_by_read.get(i)
                if not tie:
                    continue
                if len(tie) < 2:
                    d_tx[i, tie[0]] = 0.0
                    continue
                nlls = nlls_by_read.get(i)
                if not nlls or len(nlls) < 2:
                    fuzzy.append(i)  # no valid contrast -> abstain to prior/flat
                    continue
                ordered = sorted(nlls.items(), key=lambda kv: kv[1])
                best_j = ordered[0][0]
                margin = ordered[1][1] - ordered[0][1]
                # The redistribution PRIOR (covered_vote) must be fed only by reads
                # with EXPLICIT coverage: the per-read _tie_nll fallback leaves
                # cover_by_read empty (coverage uncomputable), so a missing entry
                # must NOT seed covered_vote (else the fallback path biases fuzzy
                # reads off unverified coverage).
                covered_for_vote = cover_by_read.get(i) is True
                if margin >= margin_thr:
                    # M2 distinguishes -> hard assign to the lowest-NLL candidate.
                    d_tx[i, best_j] = 0.0  # every other tie cell stays MISSING
                    if covered_for_vote:
                        # Only covered+distinguishing reads define the locus ratio.
                        covered_vote[best_j] += 1.0
                else:
                    # Indistinguishable (margin < thr), whether or not covered:
                    # defer to Pass B and redistribute by the covered-read ratio
                    # (flat 1/K only when the tie has no covered prior). The margin
                    # is thus the SOLE decider of hard-assign vs ratio-follow.
                    fuzzy.append(i)

            # Pass B: redistribute each ambiguous read across its tie in proportion
            # to the covered_vote ratio of those candidates (one-shot prior, not an
            # EM feedback loop). d_tx = -log(p) so softmax(-d_tx) reproduces p at
            # sigma=1; a zero-prior candidate stays MISSING. If NONE of the tie
            # candidates earned covered votes, fall back to the flat 1/K split
            # (== legacy behavior -> signal-dead loci are never starved).
            for i in fuzzy:
                tie = ties_by_read[i]
                masses = np.array([covered_vote[j] for j in tie], dtype=np.float64)
                s = masses.sum()
                if s > 0:
                    p = masses / s
                    for idx, j in enumerate(tie):
                        if p[idx] > 0:
                            d_tx[i, j] = -float(np.log(p[idx]))
                        # p==0 -> leave MISSING (no covered support for this cand)
                else:
                    for j in tie:
                        d_tx[i, j] = 0.0  # flat 1/K (recall-safe fallback)
            # ON-path d_tx values are already calibrated (0 / -log(p) / MISSING) for
            # sigma=1; no sigma2 normalization (that is the OFF skeleton's step).
        d_tx = d_tx.astype(np.float32)

        # Production M1/M2 assignment has no read-by-read coherence term.
        # Keep the generic EM engine at beta=0 for abundance-feedback support.
        dist_read_to_read = np.zeros((n_r, n_r), dtype=np.float32)

        R, hard_assignments, _ = em_with_coherence(
            dist_read_to_tx=d_tx,
            dist_read_to_read=dist_read_to_read,
            sigma=1.0,
            beta=0.0,
            max_iter=max_iter_em,
            tol=self.config.em_tol,
            verbose=False,
            use_gpu=self.config.use_gpu,
            abundance_feedback=self.config.abundance_feedback,
            abundance_length_norm=self.config.abundance_length_norm,
            eff_lengths=eff_lengths_fn(cand_list),
        )

        if self.config.krill_tiebreak:
            R = krill_tiebreak(
                R=R, read_ids=kept_read_ids, read_seqs=read_seqs,
                candidates=cand_list, signal_path=self.config.signal_path,
                pore=self.config.krill_pore,
                ambig_threshold=self.config.tiebreak_ambig_threshold,
                use_gpu=self.config.use_gpu,
            )
            hard_assignments = R.argmax(axis=1)

        quant_results = quantify_transcripts(
            R, hard_assignments, cand_list, kept_read_ids
        )
        # Stamp max_R BEFORE any filtering: quant_results is column-aligned to R
        # here (enumerate index == R column). Filtering would break that alignment.
        for j, qr in enumerate(quant_results):
            qr.max_R = float(R[:, j].max()) if R.shape[0] > 0 else 0.0

        return QuantOutput(
            R=R,
            hard_assignments=hard_assignments,
            quant_results=quant_results,
            cand_list=cand_list,
            kept_read_ids=kept_read_ids,
            ties_by_read=ties_by_read,
            nlls_by_read=nlls_by_read,
            cover_by_read=cover_by_read,
            read_seqs=read_seqs,
            n_ties=n_ties,
            n_refined=n_refined,
        )
