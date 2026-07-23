"""Assignment layer: read-to-candidate assignment for quant_mode="m2_em".

Extracted from the runner's ``_quant_m2_em`` so the M1/M2/EM assignment orchestration lives in
its own module (intervals, clustering, EM, quantification are already standalone modules; this is
the last piece that was glued into the runner).

- Step A: ``tie_nll`` — the per-read best-AS tie junction-NLL scoring (moved verbatim from
  ``PipelineRunner._tie_nll``; ``self.config`` -> ``config``).
- Step B: ``Assigner.assign`` — the pure assignment pipeline (mappy AS matrix + kept-read set ->
  ``tie_nll`` -> d_tx skeleton -> optional M3 coherence -> EM -> optional krill tie-break ->
  ``quantify_transcripts`` + ``max_R`` stamp) returning an immutable ``QuantOutput``. NO candidate
  dropping happens here — every structural gate/fold stays in the runner's selection tail (moved to
  the SelectionEngine in M3). ``quantify_transcripts`` is stamped here (right after EM) rather than
  after the drop cascade: byte-identical, because the intervening gates only accumulate a
  ``drop_cols`` column-index set and never mutate ``R`` / ``hard_assignments`` / ``cand_list``.

Byte-identical extraction. The runner keeps thin ``_tie_nll`` / ``_eff_lengths`` methods (passed into
``assign`` at call time) so the existing ``patch.object(PipelineRunner, ...)`` test seams survive; the
``em_with_coherence`` / ``quantify_transcripts`` seams move to this module (tests re-point here).
"""
from __future__ import annotations

import logging
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


@dataclass
class QuantOutput:
    """Immutable result of the m2_em assignment pipeline (pre-selection).

    Carries everything the runner's selection tail (and, in M3, the SelectionEngine) needs:
    the read x candidate responsibility matrix ``R`` and ``hard_assignments``, the column-aligned
    pre-selection ``quant_results`` (``max_R`` already stamped), the candidate axis ``cand_list``,
    the kept-read axis, the M2 tie evidence (``ties_by_read`` / ``nlls_by_read`` / ``cover_by_read``),
    and the log counters. Consumers treat it as READ-ONLY: selection materializes a NEW result rather
    than mutating this object (M3). ``quant_results`` is column-aligned to ``R`` (enumerate index ==
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
    beta_use: float


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
                if not nlls:
                    # nothing scorable -> flat tie (EM 1/K split among ties)
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
            margin_thr = float(getattr(self.config, "m2_diff_cover_margin", 0.0))
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
                if not nlls:
                    fuzzy.append(i)  # no junction signal -> defer (was: drop)
                    continue
                ordered = sorted(nlls.items(), key=lambda kv: kv[1])
                best_j = ordered[0][0]
                margin = (ordered[1][1] - ordered[0][1]) if len(ordered) > 1 else float("inf")
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

        # --- M3 read×read junction-window DTW coherence (opt-in; DTW is costly).
        #     Each read anchored to its d_tx-argmin (best) candidate; σ3 (median
        #     non-zero DTW) calibrates it to O(1) so β mixes cleanly. ---
        if self.config.m3_coherence:
            from fin.scoring.m3_junction_coherence import build_m3_coherence

            winner_col = np.asarray(d_tx).argmin(axis=1).astype(np.int64)
            no_data = np.asarray(d_tx).min(axis=1) >= 1e5  # all-MISSING rows
            winner_col[no_data] = -1
            dist_read_to_read = build_m3_coherence(
                kept_read_ids, read_seqs, cand_list, winner_col,
                self.config.signal_path, pore=self.config.krill_pore,
                junction_k=self.config.m2_tiebreak_junction_k, flank=_M2_EM_FLANK,
                use_gpu=self.config.use_gpu,
            )
            nz = dist_read_to_read[dist_read_to_read > 0]
            sigma3 = max(float(np.median(nz)) if nz.size else 1.0, 1e-3)
            dist_read_to_read = (dist_read_to_read / sigma3).astype(np.float32)
            beta_use = self.config.em_beta
        else:
            dist_read_to_read = np.zeros((n_r, n_r), dtype=np.float32)
            beta_use = 0.0

        R, hard_assignments, _ = em_with_coherence(
            dist_read_to_tx=d_tx,
            dist_read_to_read=dist_read_to_read,
            sigma=1.0,
            beta=beta_use,
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
            beta_use=beta_use,
        )
