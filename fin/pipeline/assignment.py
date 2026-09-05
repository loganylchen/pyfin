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


def hard_assign(R, cand_list, *, posterior_tie_atol: float = 1e-12,
                allowed_by_read=None):
    """Hard read->candidate assignment with an EXPLICIT, source-aware tie-break.

    ``R.argmax(axis=1)`` breaks an exact posterior tie by lowest COLUMN INDEX, so
    for years the decision was silently carried by the order discovery happened to
    append candidates in (annotation was appended first). That is assignment policy
    hidden in a list order: reordering emission changed which model won ties, hence
    ``num_assigned_reads``, hence which models survived the soft/hard-ratio gate
    (measured: Pr -0.93 on the 156-cell SIRV4 grid from emission order alone).

    The rule is made explicit here (Codex gpt-5.6-sol, 2026-09-04):

      1. take each read's maximum posterior;
      2. the tie set is every column within ``atol`` of it;
      3. one winner -> that column;
      4. otherwise prefer an ANNOTATED candidate: when the molecular evidence is
         numerically identical, the hypothesis the user supplied is the MAP choice
         under a weak prior. This is a decision rule for EXACT ties only -- it never
         overrides positive read evidence, and it cannot fire de novo (no GTF);
      5. remaining ties break on ``candidate_id``, never on column index, so the
         result is invariant to emission order.

    The tolerance is ABSOLUTE on purpose: ties are judged inside one read's
    normalised posterior row (rows sum to 1), not across transcript abundances,
    and a relative tolerance would swallow small but real posterior differences.
    At float32 precision 1e-12 captures exact equality only -- which is the
    stated policy.

    ``allowed_by_read`` (optional, the per-read raw-AS tie sets) restricts the tie
    set to columns the read could compete on at all. On the current M2 path a
    posterior outside that set cannot arise (cells outside the tie stay MISSING),
    so this is a defended invariant rather than a behaviour change: an empty
    intersection is a bug, and raises rather than silently falling back.

    Returns (hard_assignments, n_ties, n_annotation_ties).
    """
    import numpy as _np

    if posterior_tie_atol < 0:
        raise ValueError(f"posterior_tie_atol must be >= 0, got {posterior_tie_atol}")
    n_r, n_c = R.shape
    if n_r == 0 or n_c == 0:
        return _np.zeros(n_r, dtype=_np.int64), 0, 0
    if not _np.isfinite(R).all():
        raise ValueError("hard_assign: R contains non-finite values")
    # rank key per column: annotated first, then candidate_id -- both independent
    # of position in cand_list.
    order = sorted(range(n_c), key=lambda j: (
        getattr(cand_list[j], "source", "") != "gtf",
        getattr(cand_list[j], "candidate_id", str(j)),
    ))
    rank = _np.empty(n_c, dtype=_np.int64)
    for pos, j in enumerate(order):
        rank[j] = pos
    best = R.max(axis=1)
    tied = (best[:, None] - R) <= posterior_tie_atol
    if allowed_by_read:
        allowed = _np.zeros_like(tied)
        for i in range(n_r):
            cols = allowed_by_read.get(i)
            if not cols:
                allowed[i, :] = True          # no recorded tie set -> unrestricted
            else:
                allowed[i, list(cols)] = True
        restricted = tied & allowed
        if not restricted.any(axis=1).all():
            raise ValueError(
                "hard_assign: a read's posterior winner lies outside its raw-AS "
                "tie set (assignment invariant violated)")
        tied = restricted
    n_tied = tied.sum(axis=1)
    # pick the tied column with the smallest explicit rank
    # argmin over the RANK of tied columns, then map that rank back to a column.
    masked = _np.where(tied, rank[None, :], n_c + 1)
    hard = _np.argmin(masked, axis=1).astype(_np.int64)
    n_ties = int((n_tied > 1).sum())
    n_anno = int(sum(
        1 for i in range(n_r)
        if n_tied[i] > 1 and getattr(cand_list[hard[i]], "source", "") == "gtf"
    ))
    return hard, n_ties, n_anno


def apply_endpoint_likelihood(d_tx, ties_by_read, kept_read_ids, cand_list, ctx) -> int:
    """Add the endpoint channels (3'/polyA, and optionally the one-sided 5' start
    term) to every finite cell of a read's raw tie set.

    For read r with a PASS call and candidate c in its raw tie set,
    ``d[r,c] += -log[(1-eta) g(d_rc) + eta b(d_rc)]`` with ``d_rc`` the
    strand-aware signed distance from c's claimed TES to r's 3' end (positive =
    the read runs past the end). Then the row is re-anchored so its best cell
    stays at 0 (the EM is shift-invariant per row; this keeps the MISSING
    sentinel meaningful). Reads without a PASS call are untouched: absence of
    a call is abstention, never evidence. Returns the number of rows touched.
    """
    from fin.analysis.molecule_evidence import POLYA_PASS, polya_state

    endpoint = ctx.get("endpoint")
    start = ctx.get("start")
    polya = ctx.get("polya") or {}
    ends = ctx.get("ends") or {}
    starts = ctx.get("starts") or {}
    tes_on = endpoint is not None and getattr(endpoint, "available", False)
    tss_on = start is not None and getattr(start, "available", False)
    if not tes_on and not tss_on:
        return 0
    n_touched = 0
    for i, rid in enumerate(kept_read_ids):
        tie = ties_by_read.get(i)
        if not tie or len(tie) < 2:
            continue
        vals = {j: 0.0 for j in tie if d_tx[i, j] < 1e5}
        if len(vals) < 2:
            continue
        touched = False
        # 3' / polyA: only a PASS read carries endpoint evidence
        if tes_on:
            qc, _ = polya_state(polya.get(rid))
            loc = ends.get(rid)
            if qc == POLYA_PASS and loc is not None:
                _c, _s, e_r = loc
                for j in vals:
                    c = cand_list[j]
                    strand = getattr(c, "strand", "+")
                    t_c = c.end if strand == "+" else c.start
                    d = (e_r - t_c) if strand == "+" else (t_c - e_r)
                    vals[j] += -endpoint.loglik(d)
                touched = True
        # 5' start: a read that begins UPSTREAM of a model's claimed TSS is
        # explained by that model only through the error kernel; reads at or
        # downstream of the TSS abstain (degradation is a separate channel).
        if tss_on:
            loc5 = starts.get(rid)
            if loc5 is not None:
                _c, _s, s_r = loc5
                for j in vals:
                    c = cand_list[j]
                    strand = getattr(c, "strand", "+")
                    t_c = c.start if strand == "+" else c.end
                    u = (t_c - s_r) if strand == "+" else (s_r - t_c)   # > 0 = upstream
                    if u > 0:
                        vals[j] += -start.loglik_upstream(u) - (-start.loglik_upstream(0))
                touched = True
        if not touched:
            continue
        lo = min(vals.values())
        for j, v in vals.items():
            d_tx[i, j] += (v - lo)
        n_touched += 1
    return n_touched


def apply_start_likelihood(d_tx, ties_by_read, kept_read_ids, cand_list, models,
                           nlls_by_read=None) -> int:
    """Add the frozen read-5'-end positional NLL of folded-TSS recovery to ``d_tx``.

    Design (Codex gpt-5.6-sol, 2026-09-03): the recovery verdict only decides that a
    second model (the recovered short isoform S of parent P) exists; which model each
    read resembles is a per-read LIKELIHOOD, and EM then decides mass. For a read whose
    raw best-AS tie set is EXACTLY ``{P, S}``, with 5' offset bin ``b_i`` in P's
    transcript frame and the supported TSS bin ``b0``:

      - ``b_i == b0``: L(P) = h_P(b0) (P's frozen degradation hazard of terminating in
        that bin), L(S) = 1  ->  d[P] += -log h_P, d[S] += 0.
      - ``b_i <  b0``: a read starting 5' of S's TSS cannot come from S -> d[S] = MISSING.
      - ``b_i >  b0``: abstain (no validated short-degradation kernel yet) -> no change.

    Ties with a third candidate, reads without an offset, and every other cell are left
    untouched, so the term can never leak outside the (P, S) contrast. Mutates ``d_tx``
    in place; returns the number of rows touched.

    ``nlls_by_read`` (the per-read graded contrast the M2/M1 support gate reads) is
    updated with the SAME contrast so the gate judges "best model for this read" on
    the full evidence: a read with no scored junction contrast gets {P: pen, S: 0}
    (peak bin) or {P: 0} (upstream); a read with a scored contrast has the penalty
    added to its parent cell. The positional evidence is per-read evidence, not an
    exemption: a recovered model still has to be some read's best.
    """
    from fin.analysis.tss_evidence import read_five_prime_offset

    col = {c.candidate_id: j for j, c in enumerate(cand_list)}
    n_touched = 0
    for m in models:
        jp, js = col.get(m["parent_id"]), col.get(m["recovered_id"])
        if jp is None or js is None:
            continue
        pair = {jp, js}
        pen = -math.log(max(float(m["hazard"]), 1e-12))
        for i, rid in enumerate(kept_read_ids):
            tie = ties_by_read.get(i)
            if not tie or set(tie) != pair:
                continue
            span = m["read_spans"].get(rid)
            if span is None:
                continue
            off = read_five_prime_offset(span, m["exons"], m["strand"])
            if off is None:
                continue
            b_i = int(off) // int(m["bin_bp"])
            nl = None if nlls_by_read is None else nlls_by_read.get(i)
            if b_i == int(m["tss_bin"]):
                d_tx[i, jp] += pen
                if nlls_by_read is not None:
                    if not nl:
                        nlls_by_read[i] = {jp: pen, js: 0.0}
                    elif jp in nl:
                        nl[jp] = nl[jp] + pen
            elif b_i < int(m["tss_bin"]):
                d_tx[i, js] = MISSING
                if nlls_by_read is not None:
                    if not nl:
                        nlls_by_read[i] = {jp: 0.0}
                    else:
                        nl.pop(js, None)
            else:
                continue
            n_touched += 1
    return n_touched


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
    # The UNNORMALISED per-cell distances that produced ``R`` (0 / graded / MISSING),
    # kept because ``R`` is row-normalised and cannot be inverted. NOTE this is the
    # DECISION matrix the EM consumes: with the diff-coverage gate on (the default)
    # a distinguishing read is collapsed to one-hot and a fuzzy read carries the
    # component-wide covered_vote ratio. Read-only for consumers.
    d_tx: object = None
    # Per-read junction-NLL contrasts over the raw best-AS tie set, built BEFORE
    # the gate turned them into decisions. This is the EVIDENCE matrix the
    # candidate-existence test reads; see fin/analysis/candidate_existence.py.
    structural_d_tx_raw: object = None
    # Raw mappy AS matrix over kept reads x candidates (-inf = no hit): the A_r
    # channel of the molecule evidence row. Read-only.
    raw: object = None


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
        start_likelihood_context=None,
        endpoint_context=None,
    ) -> Optional[QuantOutput]:
        """Run the m2_em assignment pipeline; return a QuantOutput or None (empty interval).

        ``start_likelihood_context`` (default None = no new branch, byte-identical) is the
        frozen read-5'-end model list produced by folded-TSS shadow recovery
        (:func:`fin.analysis.tss_shadow_recovery.recover_folded_tss_candidates`) and is
        passed ONLY on the second assignment of an interval that re-emitted a shadow. It
        adds a positional NLL term to ``d_tx`` (see :func:`apply_start_likelihood`).
        """
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
        # RAW structural kernel, built BEFORE any gate transform: for each read,
        # the per-candidate junction NLL over its raw best-AS tie set (MISSING
        # elsewhere). The gated `d_tx` below is a DECISION matrix -- a
        # distinguishing read becomes one-hot and a fuzzy read carries the
        # component-wide covered_vote ratio rather than its own likelihood -- so
        # it cannot be used as evidence about whether a candidate exists. The
        # existence test consumes this matrix instead. Read-only for consumers.
        structural_d_tx_raw = np.full((n_r, n_c), MISSING, dtype=np.float64)
        for i in range(n_r):
            tie = ties_by_read.get(i)
            if not tie:
                continue
            scored = nlls_by_read.get(i) or {}
            if len(scored) >= 2:
                lo = min(scored.values())
                for j, v in scored.items():
                    structural_d_tx_raw[i, j] = v - lo
            elif len(tie) == 1:
                structural_d_tx_raw[i, tie[0]] = 0.0
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
            # The junction channel compares STRUCTURES. Candidates that share an
            # intron chain (an annotated model and its read-extent sibling) are
            # one structural hypothesis here: they score identically, so the
            # margin, the hard assignment and the covered-vote prior are all
            # taken over unique chains, and every member of the winning chain
            # receives the decision (the endpoint channel and the mixture
            # weights then split the read within the group). Without this the
            # twin would always be the "runner-up" at margin 0, the read would
            # be declared indistinguishable, and -- the twins never earning a
            # vote -- its mass would flow to wobble variants that did.
            # Single-exon candidates have no junctions to compare and keep
            # their own identity. Singleton groups == legacy behaviour.
            def _structure_key(c):
                if getattr(c, "source", "") == "fusion":
                    # a fusion's intron_chain holds arm A only; its structure is
                    # its own (both arms), never shared with a linear transcript
                    return ("fusion", c.candidate_id)
                introns = tuple(getattr(c.intron_chain, "introns", ()) or ())
                return introns if introns else ("mono", c.candidate_id)
            group_of = [_structure_key(c) for c in cand_list]
            covered_vote: Dict[object, float] = {}
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
                by_group: Dict[object, float] = {}
                for j, v in nlls.items():
                    k = group_of[j]
                    by_group[k] = v if k not in by_group else min(by_group[k], v)
                if len(by_group) < 2:
                    fuzzy.append(i)  # one structure scored -> no contrast
                    continue
                ordered = sorted(by_group.items(), key=lambda kv: kv[1])
                best_k = ordered[0][0]
                margin = ordered[1][1] - ordered[0][1]
                # The redistribution PRIOR (covered_vote) must be fed only by reads
                # with EXPLICIT coverage: the per-read _tie_nll fallback leaves
                # cover_by_read empty (coverage uncomputable), so a missing entry
                # must NOT seed covered_vote (else the fallback path biases fuzzy
                # reads off unverified coverage).
                covered_for_vote = cover_by_read.get(i) is True
                if margin >= margin_thr:
                    # M2 distinguishes -> hard assign to the lowest-NLL structure.
                    for j in tie:
                        if group_of[j] == best_k:
                            d_tx[i, j] = 0.0  # every other tie cell stays MISSING
                    if covered_for_vote:
                        # Only covered+distinguishing reads define the locus ratio.
                        covered_vote[best_k] = covered_vote.get(best_k, 0.0) + 1.0
                else:
                    # Indistinguishable (margin < thr), whether or not covered:
                    # defer to Pass B and redistribute by the covered-read ratio
                    # (flat 1/K only when the tie has no covered prior). The margin
                    # is thus the SOLE decider of hard-assign vs ratio-follow.
                    fuzzy.append(i)

            # Pass B: redistribute each ambiguous read across the STRUCTURES in
            # its tie in proportion to their covered_vote ratio (one-shot prior,
            # not an EM feedback loop). d_tx = -log(p) so softmax(-d_tx)
            # reproduces p at sigma=1; a zero-prior structure stays MISSING. If
            # NONE of the tie's structures earned covered votes, fall back to the
            # flat 1/K split (== legacy behavior -> signal-dead loci are never
            # starved).
            for i in fuzzy:
                tie = ties_by_read[i]
                groups = []
                for j in tie:
                    if group_of[j] not in groups:
                        groups.append(group_of[j])
                masses = np.array([covered_vote.get(k, 0.0) for k in groups], dtype=np.float64)
                s = masses.sum()
                if s > 0:
                    p = {k: m / s for k, m in zip(groups, masses)}
                    size = {k: sum(1 for j in tie if group_of[j] == k) for k in groups}
                    for j in tie:
                        k = group_of[j]
                        if p[k] > 0:
                            # the STRUCTURE's share p_k is split evenly over its
                            # members so softmax(-d) conserves p_k per structure
                            d_tx[i, j] = -float(np.log(p[k] / size[k]))
                        # p==0 -> leave MISSING (no covered support for this structure)
                else:
                    for j in tie:
                        d_tx[i, j] = 0.0  # flat 1/K (recall-safe fallback)
            # ON-path d_tx values are already calibrated (0 / -log(p) / MISSING) for
            # sigma=1; no sigma2 normalization (that is the OFF skeleton's step).
        if endpoint_context:
            # 3'/polyA channel (Part A8): calibrated -log p, added AFTER the junction
            # sigma2 normalisation. A read without a PASS call abstains.
            apply_endpoint_likelihood(d_tx, ties_by_read, kept_read_ids, cand_list,
                                      endpoint_context)
        if start_likelihood_context:
            # Positional (read 5' end) NLL is calibrated -log p already: it is added
            # AFTER the junction sigma2 normalisation, never divided by it.
            apply_start_likelihood(
                d_tx, ties_by_read, kept_read_ids, cand_list, start_likelihood_context,
                nlls_by_read=nlls_by_read)
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

        # Hard assignment is recomputed here UNCONDITIONALLY, after every
        # transform of R, by the explicit rule (never by column order).
        hard_assignments, n_hard_ties, n_anno_ties = hard_assign(
            R, cand_list, allowed_by_read=ties_by_read)
        if n_hard_ties:
            logger.debug("hard assignment: %d exact ties (%d resolved to annotation)",
                         n_hard_ties, n_anno_ties)

        quant_results = quantify_transcripts(
            R, hard_assignments, cand_list, kept_read_ids
        )
        # Stamp max_R BEFORE any filtering: quant_results is column-aligned to R
        # here (enumerate index == R column). Filtering would break that alignment.
        for j, qr in enumerate(quant_results):
            qr.max_R = float(R[:, j].max()) if R.shape[0] > 0 else 0.0

        return QuantOutput(
            d_tx=d_tx,
            structural_d_tx_raw=structural_d_tx_raw,
            raw=raw,
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
