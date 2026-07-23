"""Assignment layer: read-to-candidate assignment for quant_mode="m2_em".

Extracted from the runner's ``_quant_m2_em`` so the M1/M2/EM assignment orchestration lives in
its own module (intervals, clustering, EM, quantification are already standalone modules; this is
the last piece that was glued into the runner).

Step A (this file): ``tie_nll`` — the per-read best-AS tie junction-NLL scoring (moved verbatim
from ``PipelineRunner._tie_nll``; ``self.config`` -> ``config``). Byte-identical.
"""
from __future__ import annotations

import logging
from typing import Dict, List, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Junction-window flank (bp) for the M2 tie discrimination window (matches runner).
_M2_EM_FLANK = 2


def tie_nll(config, kept_read_ids, read_seqs, cand_list, aligners, raw):
    """Per-read best-AS tie junction-NLL (the m2_resolve_tie niche).

    For each read whose raw mappy AS is simultaneously best across >=2 candidates, build the
    class junction-discrimination window over that tie set and score each tied candidate with the
    per-event junction mean-NLL (non-HMM krill eventalign over the window). Only tied,
    window-scorable cells are touched -- never the dense read×candidate matrix.

    Returns (nlls_by_read {i: {j: nll}}, ties_by_read {i: [j,...]}, n_ties, n_refined,
    cover_by_read {i: bool}). ``cover_by_read`` is populated only when the diff-region coverage
    gate (config.m2_diff_cover_gate) is ON; otherwise it is empty.
    """
    import krill

    from fin.scoring.krill_aligner import (
        krill_thread_count,
        make_krill_aligner,
    )
    from fin.scoring.m2_junction_nll import (
        _mean_nll_in_gset,
        class_junction_window_set,
        event_genomic_positions,
        read_cand_mean_nll,
        read_straddles,
        wobble_diff_spans,
    )
    from fin.scoring.mappy_score import score_hit

    gate = bool(getattr(config, "m2_diff_cover_gate", False))
    num_threads = krill_thread_count()
    n_c = len(cand_list)
    nlls_by_read: Dict[int, Dict[int, float]] = {}
    ties_by_read: Dict[int, List[int]] = {}
    cover_by_read: Dict[int, bool] = {}
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
        return nlls_by_read, ties_by_read, 0, 0, cover_by_read

    k = config.m2_tiebreak_junction_k
    pore = config.krill_pore
    sig_path = config.signal_path

    # --- Pass 1: tie detection + per-(read, candidate) mappy best-hit slices.
    # Each tied cell is reduced to the SAME (sliced sequence, start offset)
    # the singular read_cand_mean_nll would have eventaligned, so the batched
    # eventalign below is byte-identical to the per-call path -- only the
    # Python/GIL/per-call overhead is removed.
    gset_by_read: Dict[int, set] = {}
    spans_by_read: Dict[int, List[Tuple[int, int]]] = {}  # diff-cover gate spans
    reads_variants: Dict[str, Dict[str, str]] = {}  # rid -> {cand_id: seq}
    starts: Dict[str, Dict[str, int]] = {}          # rid -> {cand_id: r_st}
    cidj_by_read: Dict[str, Dict[str, int]] = {}     # rid -> {cand_id: col j}
    n_ties = 0
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
        gset = class_junction_window_set(
            [cand_list[j] for j in tie], flank=_M2_EM_FLANK, k=k
        )
        if not gset:
            continue
        if gate:
            spans_by_read[i] = wobble_diff_spans(
                [cand_list[j] for j in tie], flank=_M2_EM_FLANK
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
            continue
        gset_by_read[i] = gset
        reads_variants[rid] = per_seqs
        starts[rid] = per_starts
        cidj_by_read[rid] = per_cidj

    if not reads_variants:
        return nlls_by_read, ties_by_read, n_ties, 0, cover_by_read

    # --- Pass 2: ONE batched eventalign over every tied pair (per-pair start
    # via the {label: offset} form); per-read singular fallback on absence or
    # batch failure. Both branches score with the same _mean_nll_in_gset.
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
        for rid, res_list in batch_out.items():
            i = rid_to_i.get(rid)
            if i is None:
                continue
            gset = gset_by_read.get(i)
            if gset is None:
                continue
            by_label = {x.get("variant_label"): x for x in res_list}
            nlls: Dict[int, float] = {}
            spans = spans_by_read.get(i) if gate else None
            span_cov = [False] * len(spans) if spans else None
            for cid, j in cidj_by_read.get(rid, {}).items():
                res = by_label.get(cid)
                if res is None or res.get("status", -1) != 0:
                    continue
                nll, n_ev = _mean_nll_in_gset(res, cand_list[j], gset)
                if n_ev > 0 and np.isfinite(nll):
                    nlls[j] = float(nll)
                if span_cov is not None:
                    egp = event_genomic_positions(res, cand_list[j])
                    if egp:
                        for s_idx, (lo, hi) in enumerate(spans):
                            if not span_cov[s_idx] and read_straddles(egp, lo, hi):
                                span_cov[s_idx] = True
            if nlls:
                n_refined += 1
                nlls_by_read[i] = nlls
                if gate:
                    # covered iff every wobbling diff span is straddled by the
                    # read on >=1 candidate (empty spans -> vacuously covered).
                    cover_by_read[i] = all(span_cov) if span_cov else True
        return nlls_by_read, ties_by_read, n_ties, n_refined, cover_by_read

    for rid, per_cidj in cidj_by_read.items():
        i = rid_to_i.get(rid)
        if i is None:
            continue
        gset = gset_by_read.get(i)
        seq = read_seqs[rid]
        nlls = {}
        for j in per_cidj.values():
            nll, n_ev = read_cand_mean_nll(
                rid, seq, cand_list[j], [], krill_aligner, aligners[j],
                sig_path, pore, gset=gset, use_gpu=eff_gpu,
                num_thread=num_threads,
            )
            if n_ev > 0 and np.isfinite(nll):
                nlls[j] = float(nll)
        if nlls:
            n_refined += 1
            nlls_by_read[i] = nlls
    # Per-read fallback path does not expose eventalign positions, so the
    # coverage map is left empty here; the d_tx builder defaults missing reads
    # to covered=True (recall-safe: never drops on the rare fallback path).
    return nlls_by_read, ties_by_read, n_ties, n_refined, cover_by_read
