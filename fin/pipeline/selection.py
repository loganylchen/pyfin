"""Selection layer: decide which candidates survive, in one place.

Historically "which candidates are real" had no single home — it was smeared across three
incompatible mechanisms: pre-EM gates mutating the candidate list, a post-EM ``drop_cols`` column
set buried inside ``_quant_m2_em``, and finalize-time dict deletion. This module is the single home
for the post-EM INTERVAL scope (M3 Step 1): the cascade that runs after assignment, within one
interval, over the ``QuantOutput`` of the m2_em Assigner.

Design (behavior-preserving; see the pipeline reorganization plan + Codex review 2):
- ``SelectionOutcome`` is the uniform record of every keep/drop/fold decision (free FN/FP-attribution
  provenance, replacing the ad-hoc ``logger.info`` drop counters).
- ``select_m2_interval`` reproduces the CURRENT sequential cascade EXACTLY — each gate evaluates
  against the survivor set left by the earlier gates (``exclude=drop_cols``), the two folds
  (containment-collapse, mono-resolve) move reads/mass into their parent before the final apply, and
  ordering is load-bearing (it is the byte-identity guarantee). The cleaner "evaluate all on one base
  snapshot + declared precedence" model is a separate, behavior-CHANGING follow-up, not this step.

Scope: the DEFAULT ``m2_em`` path only. ``argmax`` / ``m1_em`` / ``cluster`` keep their own survival
handling. The PRE_ASSIGN (canonical / junction-dominance) and GLOBAL (finalize filters) scopes move
here in later M3 steps.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from fin.analysis.quantification import (
    fulllen_fraction_drops,
    isoform_fraction_drops,
    mono_exon_drops,
    soft_mass_ratio_drops,
)

logger = logging.getLogger(__name__)


@dataclass
class SelectionOutcome:
    """One keep/drop/fold decision for a candidate.

    ``action`` is ``"drop"`` or ``"fold"`` (surviving candidates are simply not recorded).
    ``fold_into`` is the parent candidate_id for a fold. ``scope`` is the selection phase
    (``"interval"`` here). ``reason`` names the gate that made the decision.
    """

    candidate_id: str
    action: str
    reason: str
    scope: str = "interval"
    fold_into: Optional[str] = None


def select_m2_interval(
    config,
    qo,
    candidate_set,
    interval,
    observed_junctions_fn: Callable,
) -> Tuple[List, List[SelectionOutcome]]:
    """Run the post-EM INTERVAL selection cascade over an m2_em ``QuantOutput``.

    Byte-identical move of the runner's former ``_quant_m2_em`` selection tail. Returns
    ``(surviving_quant_results, outcomes)``. ``observed_junctions_fn`` is the runner's
    ``_observed_junctions`` bound method (passed at call time so the test seam survives); it is
    invoked lazily and only when a gate actually needs directly-observed read junctions.
    """
    cand_list = qo.cand_list
    R = qo.R
    ties_by_read = qo.ties_by_read
    nlls_by_read = qo.nlls_by_read
    quant_results = qo.quant_results

    outcomes: List[SelectionOutcome] = []
    drop_cols: set = set()

    def _record_drops(cols, reason: str) -> None:
        for j in cols:
            outcomes.append(SelectionOutcome(
                cand_list[j].candidate_id, "drop", reason,
            ))

    # --- Cluster-internal wobble recheck (post-EM, data-layer shadow kill).
    #     Cluster ALL multi-exon candidates by structure (same intron count +
    #     every junction within +-bp); within a cluster the highest-EM-abundance
    #     candidate anchors (often the true isoform, frequently a GTF
    #     passthrough), and a NOVEL sibling whose abundance/anchor < fraction is
    #     a wobble shadow -> drop. GTF/fusion never dropped (they only anchor).
    #     Pure abundance evidence: GTF participates in the abundance race but is
    #     NOT used as a correctness oracle, so this stays robust when the
    #     annotation is wrong (and self-disables with no GTF: a real novel
    #     anchors its own cluster). OFF -> no drops (byte-identical). ---
    # Strand-keyed {strand: Counter{(donor,acceptor): n_reads}} of directly-
    # observed read junctions; built lazily (reopens BAM) and shared by the
    # cluster-recheck GTF guard and the Lever-2 junction-support gate.
    observed_jct = None
    if getattr(config, "m2_cluster_recheck", False):
        from fin.scoring.m2_junction_nll import structural_wobble_clusters

        frac = config.m2_cluster_recheck_fraction
        if frac <= 0.0:
            frac = config.min_isoform_fraction
        if frac > 0.0:
            em_ab = np.asarray(R).sum(axis=0)
            bp = config.m2_cluster_recheck_bp
            cassette_bp = getattr(
                config, "m2_cluster_recheck_cassette_max_exon_bp", 0
            )
            clusters = structural_wobble_clusters(
                cand_list, bp, cassette_max_exon_bp=cassette_bp
            )
            from fin.scoring.m2_junction_nll import (
                gtf_guard_needed,
                wobble_shadow_drops,
            )

            gtf_drop = bool(getattr(
                config, "m2_cluster_recheck_novel_displaces_gtf", True
            ))
            # Build the directly-observed read junctions for the GTF read-support
            # guard ONLY when a clustered low-abundance GTF sibling could actually
            # be dropped — this reopens the BAM, so skip it for no-GTF intervals
            # (and avoids touching the BAM for mocked/placeholder-path callers).
            if gtf_drop and gtf_guard_needed(cand_list, em_ab, clusters, frac):
                observed_jct = observed_junctions_fn(interval)
            _d = wobble_shadow_drops(
                cand_list, em_ab, clusters, frac,
                gtf_drop_enabled=gtf_drop,
                observed_jct=observed_jct,
                jct_tol=int(getattr(config, "m2_cluster_recheck_jct_tol", 2)),
                gtf_min_jct_reads=int(getattr(
                    config, "m2_cluster_recheck_gtf_min_jct_reads", 1
                )),
            )
            drop_cols |= _d
            _record_drops(_d, "cluster_recheck")

    # --- Lever 2: per-junction read-support gate. Drop a NOVEL multi-exon
    #     candidate if ANY of its junctions is spliced by fewer than
    #     novel_junction_min_reads directly-observed reads (CIGAR introns,
    #     strand-keyed, matched within novel_junction_reads_tol bp) — a novel
    #     junction must be carried by >= N reads, not just 1. gtf/fusion/mono
    #     exempt. <=1 disables (byte-identical). Reuses observed_jct if the
    #     cluster-recheck already built it; otherwise builds it here. ---
    _njmr = int(getattr(config, "novel_junction_min_reads", 0))
    if _njmr > 1:
        from fin.scoring.m2_junction_nll import junction_support_drops
        if observed_jct is None:
            observed_jct = observed_junctions_fn(interval)
        _d = junction_support_drops(
            cand_list, observed_jct, min_reads=_njmr,
            tol=int(getattr(config, "novel_junction_reads_tol", 2)),
        )
        drop_cols |= _d
        _record_drops(_d, "novel_junction_support")

    # --- GUIDED junction-support gate (EXPERIMENT, default-off). Mirror of
    #     Lever 2 for GTF-passthrough candidates: drop a guided multi-exon
    #     candidate if ANY of its junctions lacks >= guided_junction_min_reads
    #     directly-observed reads (CIGAR introns within guided_junction_reads
    #     _tol bp). Extends the coordinate-EXACT read-support check — which
    #     Lever 2 applies only to novel candidates — to GTF junctions, so a
    #     jitter-corrupted annotation junction (coords > tol from the true
    #     site, thus zero exact support) is dropped instead of surviving the
    #     coordinate-inexact M2/M1 gate. Reuses observed_jct. 0 disables
    #     (byte-identical). NOT recall-safe (also drops low-coverage annotated
    #     junctions); default-off pending real-data validation. ---
    _gjmr = int(getattr(config, "guided_junction_min_reads", 0))
    if _gjmr > 0:
        from fin.scoring.m2_junction_nll import guided_junction_support_drops
        if observed_jct is None:
            observed_jct = observed_junctions_fn(interval)
        _d = guided_junction_support_drops(
            cand_list, observed_jct, min_reads=_gjmr,
            tol=int(getattr(config, "guided_junction_reads_tol", 2)),
        )
        drop_cols |= _d
        _record_drops(_d, "guided_junction_support")

    # --- M2/M1 read-support gate. A multi-exon candidate (GTF or novel;
    #     fusion/mono exempt) must earn >=1 read's support: either it is some
    #     read's M1 SOLE best-AS (tie set == just this candidate), or it is
    #     some read's M2-best (lowest junction-NLL in that read's tie). A
    #     jitter-corrupted GTF junction wins neither (loses the M2 contest to
    #     the true-junction candidate and never holds a read's sole AS), so it
    #     is dropped; a genuine isoform always earns one or the other. ---
    if getattr(config, "m2_support_gate", False):
        from fin.scoring.m2_junction_nll import support_gate_drops

        _d = support_gate_drops(
            cand_list, ties_by_read, nlls_by_read,
            tie_ok=bool(getattr(config, "m2_support_gate_tie", True)),
        )
        drop_cols |= _d
        _record_drops(_d, "m2_support_gate")

    # --- Containment-cluster drop: drop a NOVEL candidate whose intron chain is a
    #     contiguous SUB-CHAIN of a longer candidate (a truncation / exon-skip shadow
    #     the same-intron-count wobble cluster never groups) when it is a low-support
    #     shadow by BOTH EM abundance AND supporting-read count. Runs AFTER all the
    #     structural/support gates so ``exclude=drop_cols`` covers every already-doomed
    #     parent — a shadow is never folded into a parent that is itself dropped. The
    #     read-support guard keeps most genuine low-abundance short/alt-TSS isoforms;
    #     gtf/fusion never dropped. DEFAULT-ON (--no-containment-cluster disables). ---
    if getattr(config, "containment_cluster", False):
        from fin.scoring.m2_junction_nll import containment_cluster_drops
        em_ab_c = np.asarray(R).sum(axis=0)
        read_counts = [len(getattr(c, "supporting_read_ids", ()) or ())
                       for c in cand_list]
        _d = containment_cluster_drops(
            cand_list, em_ab_c, read_counts,
            wobble_bp=int(getattr(config, "containment_cluster_wobble_bp", 6)),
            min_ab_ratio=float(getattr(
                config, "containment_cluster_min_ab_ratio", 0.3)),
            min_read_ratio=float(getattr(
                config, "containment_cluster_min_read_ratio", 0.3)),
            max_shadow_reads=int(getattr(
                config, "containment_cluster_max_shadow_reads", 0)),
            exclude=set(drop_cols),
        )
        drop_cols |= _d
        _record_drops(_d, "containment_cluster")

    # --- Lever 1: containment / 5'-truncation collapse (post-EM mass-fold).
    #     Fold a NOVEL candidate whose intron chain is a pure 3' suffix of a
    #     longer candidate (a 5'-truncation shadow) INTO that parent: the
    #     shadow's hard reads (read-id union) + soft EM mass move to the parent,
    #     then the shadow joins drop_cols. Parents already in drop_cols
    #     (wobble/support) are EXCLUDED, so a fold never targets an
    #     ALREADY-dropped interval candidate; a folded parent CAN still be
    #     removed by a later _finalize_and_write filter (isoform-fraction /
    #     soft-mass / fulllen / polyA), which would lose the absorbed reads --
    #     an accepted limitation given this lever is default-off. Chained
    #     containment resolves to the terminal longest parent. OFF (default)
    #     -> no folds (byte-identical). NOT recall-safe by
    #     construction (a 3'-suffix also fits a real alt-TSS isoform); see
    #     containment_shadow_drops -> kept default-off pending real-data tuning.
    if getattr(config, "containment_collapse", False):
        from fin.scoring.m2_junction_nll import containment_shadow_drops

        em_ab_c = np.asarray(R).sum(axis=0)
        fold = containment_shadow_drops(
            cand_list, em_ab_c,
            tol_bp=int(config.containment_3p_tol_bp),
            min_ratio=float(config.containment_min_abundance_ratio),
            exclude=drop_cols,
        )
        if fold:
            qr_by_id = {q.candidate_id: q for q in quant_results}
            n_folded = 0
            for shadow_col, parent_col in fold.items():
                sq = qr_by_id.get(cand_list[shadow_col].candidate_id)
                pq = qr_by_id.get(cand_list[parent_col].candidate_id)
                if sq is None or pq is None:
                    continue
                # Hard reads are disjoint by construction (EM argmax mask), so
                # the union count == sum, but union is robust either way. Soft
                # mass is added so aggregate_across_intervals (single interval:
                # unique/weight ratio == 1) reports parent + shadow abundance.
                union = set(pq.assigned_read_ids) | set(sq.assigned_read_ids)
                pq.assigned_read_ids = tuple(sorted(union))
                pq.num_assigned_reads = len(union)
                pq.abundance += sq.abundance
                # Keep the parent's max EM responsibility consistent with the
                # absorbed reads (max_R is a reporting/FP-analysis metric, not a
                # filter input; confidence is left as the parent's own mean).
                pq.max_R = max(pq.max_R, sq.max_R)
                drop_cols.add(shadow_col)
                outcomes.append(SelectionOutcome(
                    cand_list[shadow_col].candidate_id, "fold",
                    "containment_collapse",
                    fold_into=cand_list[parent_col].candidate_id,
                ))
                n_folded += 1
            if n_folded:
                logger.info(
                    "m2_em interval %s: containment folded %d "
                    "5'-truncation shadows", interval.region_string, n_folded,
                )

    # --- Post-EM mono-exon resolution (mono_resolve_post_em). Re-resolve each
    #     SURVIVING single-exon candidate's reads against the SURVIVING multi-exon
    #     candidates by strict strand-aware exonic containment: a read wholly inside
    #     one exon of exactly one surviving multi folds into it (a degradation
    #     fragment); inside several -> the highest-EM-abundance one; uncovered reads
    #     stay on the mono candidate, which is dropped if it retains fewer than
    #     mono_resolve_min_reads. Runs AFTER every structural gate, so it can only
    #     target survivors (guard 1: never fold into a doomed candidate). Containment
    #     is strict + strand-aware with terminal slop (guard 2). Multi-cover assigns
    #     the whole hard read to the top-abundance candidate, not a fractional split
    #     (guard 3); the uncovered read-support floor gates genuine mono calls. ---
    if getattr(config, "mono_resolve_post_em", False):
        from fin.scoring.m2_junction_nll import mono_resolve_drops
        mono_drops = mono_resolve_drops(
            cand_list, quant_results, np.asarray(R).sum(axis=0),
            getattr(candidate_set, "read_spans", None) or {},
            exclude=drop_cols,
            slop_bp=int(getattr(config, "mono_resolve_slop_bp", 10)),
            min_reads=int(getattr(config, "mono_resolve_min_reads", 2)),
        )
        if mono_drops:
            drop_cols |= mono_drops
            _record_drops(mono_drops, "mono_resolve")
            logger.info(
                "m2_em interval %s: mono-resolve dropped %d under-supported "
                "single-exon candidates", interval.region_string, len(mono_drops),
            )

    if drop_cols:
        drop_ids = {cand_list[j].candidate_id for j in drop_cols}
        quant_results = [q for q in quant_results if q.candidate_id not in drop_ids]
        logger.info(
            "m2_em interval %s: post-EM gates dropped %d candidates "
            "(cluster-recheck / support / containment)",
            interval.region_string, len(drop_ids),
        )

    return quant_results, outcomes


def select_global(
    config,
    aggregated: Dict[str, "object"],
    apply_polya5p_fn: Callable,
) -> Tuple[Dict[str, "object"], List[SelectionOutcome]]:
    """Run the GLOBAL post-aggregate selection cascade over the aggregated transcript dict.

    Byte-identical move of the filter section of the runner's ``_finalize_and_write``: the A4
    per-source abundance floor, then isoform-fraction / soft-mass / mono-exon / full-length /
    polyA drops, each gated (as before) by ``enable_score_filter`` and its own threshold, over the
    survivors of the earlier filters (sequential, locus-relative). ``apply_polya5p_fn`` is the
    runner's ``_apply_polya5p_filter`` bound method (needs the signal reader), passed at call time.
    Returns ``(filtered_aggregated, outcomes)``. Diagnostic dumps, gene-id resolution, and the
    GTF/TSV/BEDPE writers stay in the runner (M4 splits the writers).
    """
    outcomes: List[SelectionOutcome] = []

    def _record(cids, reason: str) -> None:
        for cid in cids:
            outcomes.append(SelectionOutcome(cid, "drop", reason, scope="global"))

    # R5 ablation: enable_score_filter=False disables all post-EM filters so
    # the ablation row sees unfiltered EM output (AC7).
    _score_filter_on = getattr(config, "enable_score_filter", True)

    # A4: post-EM abundance floor with per-source thresholds (on soft EM
    # abundance). NOVEL transcripts face min_abundance (the `fin` CLI defaults
    # it to 3). GTF transcripts face their own lighter floor min_gtf_abundance
    # (CLI default 1) so a GTF candidate whose EM soft-mass is below 1 read is
    # dropped — pyfin is not a copy-annotation tool — while genuine annotated
    # transcripts with real EM support are kept. When floor_gtf_abundance is
    # set (--floor-gtf-abundance) GTF is raised to the NOVEL floor — precisely
    # max(min_gtf_abundance, min_abundance) so the toggle never LOWERS an
    # explicit GTF floor — reproducing the SIRV gffcompare T>=3 with-GTF
    # operating point. Fusion is ALWAYS exempt so opt-in fusion calls are
    # never silently dropped. SIRV-tuned: these floors are overfit and gut
    # genuine low-abundance isoform recall on real data; disable via
    # --min-abundance 0 / --min-gtf-abundance 0.
    # NOTE: the separate ablation harness (fin/ablation/runner.py) keeps GTF
    # exempt regardless — that divergence is intentional (no-GTF de novo).
    floor_gtf = config.floor_gtf_abundance
    gtf_floor = (
        max(config.min_gtf_abundance, config.min_abundance)
        if floor_gtf
        else config.min_gtf_abundance
    )
    if _score_filter_on and (config.min_abundance > 0.0 or gtf_floor > 0.0):
        before = len(aggregated)
        kept = {
            cid: qr
            for cid, qr in aggregated.items()
            if qr.source == "fusion"
            or (
                qr.abundance >= gtf_floor
                if qr.source == "gtf"
                else qr.abundance >= config.min_abundance
            )
        }
        dropped = before - len(kept)
        if dropped:
            _record(set(aggregated) - set(kept), "min_abundance")
            logger.info(
                "Dropped %d transcripts with abundance < floor "
                "(novel < %.3f, gtf < %.3f)",
                dropped,
                config.min_abundance,
                gtf_floor,
            )
        aggregated = kept

    # Minimum isoform fraction (locus-relative abundance) filter. Drops
    # NOVEL multi-exon transcripts whose abundance is below
    # `min_isoform_fraction` of the dominant overlapping novel isoform at
    # the same locus — the standard Cufflinks --min-isoform-fraction /
    # StringTie -f heuristic (the low-fraction tail is enriched for
    # incompletely-spliced precursors and assembly artifacts).
    # gtf/fusion/mono candidates are EXEMPT. SIRV WARNING: the conservative
    # literature default (0.01) is used; SIRV's F1-optimal ~0.4 is overfit.
    if _score_filter_on and config.min_isoform_fraction > 0.0:
        drop_ids = isoform_fraction_drops(
            aggregated, config.min_isoform_fraction
        )
        if drop_ids:
            aggregated = {
                cid: qr
                for cid, qr in aggregated.items()
                if cid not in drop_ids
            }
            _record(drop_ids, "isoform_fraction")
            logger.info(
                "Dropped %d novel transcripts with isoform fraction < %.3f",
                len(drop_ids),
                config.min_isoform_fraction,
            )

    # Soft-mass / hard-read ratio filter. Drops a NOVEL multi-exon candidate
    # whose EM soft abundance is inflated far above its honest hard-read
    # (argmax) count — the signature of a wobble shadow that steals little
    # real read support but accumulates fractional soft crumbs from a
    # high-abundance structural near-copy (anchor). Catches the HIGH-relative
    # -abundance shadows that the isoform_fraction / cluster-recheck levers
    # cannot (a shadow at 32% of its anchor passes a fraction gate). Pure EM
    # evidence; gtf/fusion/mono exempt. max_soft_mass_ratio<=0 disables.
    if _score_filter_on and config.max_soft_mass_ratio > 0.0:
        drop_ids = soft_mass_ratio_drops(
            aggregated, config.max_soft_mass_ratio
        )
        if drop_ids:
            aggregated = {
                cid: qr
                for cid, qr in aggregated.items()
                if cid not in drop_ids
            }
            _record(drop_ids, "soft_mass_ratio")
            logger.info(
                "Dropped %d novel transcripts with soft/hard read ratio >= %.2f",
                len(drop_ids),
                config.max_soft_mass_ratio,
            )

    # Lever 3 — mono-exon (single-exon) read-support gate. Drop a NOVEL
    # single-exon candidate whose hard read count < min_mono_exon_reads OR
    # genomic length < min_mono_exon_length. Suppresses single-exon de novo
    # noise (IsoQuant drops novel unspliced by default for ONT) WITHOUT a
    # blanket drop — a high-support/long real intronless gene survives.
    # gtf/fusion/multi-exon exempt. Needs the master switch AND >=1 threshold
    # > 0; otherwise no drops (byte-identical).
    if (
        _score_filter_on
        and getattr(config, "drop_mono_exon_novel", False)
        and (
            config.min_mono_exon_reads > 0
            or config.min_mono_exon_length > 0
        )
    ):
        drop_ids = mono_exon_drops(
            aggregated,
            min_reads=config.min_mono_exon_reads,
            min_len=config.min_mono_exon_length,
        )
        if drop_ids:
            aggregated = {
                cid: qr
                for cid, qr in aggregated.items()
                if cid not in drop_ids
            }
            _record(drop_ids, "mono_exon")
            logger.info(
                "Dropped %d novel single-exon transcripts (mono gate: "
                "min_reads=%d min_len=%d)",
                len(drop_ids),
                config.min_mono_exon_reads,
                config.min_mono_exon_length,
            )

    # Full-length end-coherence filter (FLAIR/TALON-style full-length read
    # support). Drops NOVEL multi-exon (>=2 intron) transcripts whose
    # fraction of full-length assigned reads (read genomic 5' AND 3' both
    # within fulllen_window_bp of the candidate's ends) is below
    # min_fulllen_fraction. The fulllen_frac METRIC itself is signal-free
    # (BAM primary-alignment spans only), but it is computed per-interval
    # over the quant assignment population (the m2_em/m1_em EM hard
    # assignments or, under argmax, the M2-tiebreak picks); candidates with
    # fulllen_frac < 0 (unreachable or never scored) are EXEMPT and never
    # dropped. gtf/fusion/mono exempt. ORTHOGONAL to the isoform-fraction
    # filter above (the two stack on the competitor metric). SIRV WARNING:
    # the default 0.1 is SIRV-tuned; re-tune or disable on real data.
    if _score_filter_on and config.min_fulllen_fraction > 0.0:
        drop_ids = fulllen_fraction_drops(
            aggregated, config.min_fulllen_fraction
        )
        if drop_ids:
            aggregated = {
                cid: qr
                for cid, qr in aggregated.items()
                if cid not in drop_ids
            }
            _record(drop_ids, "fulllen_fraction")
            logger.info(
                "Dropped %d novel transcripts with full-length fraction < %.3f",
                len(drop_ids),
                config.min_fulllen_fraction,
            )

    # polyA + 5'-proximity candidate-retention filter. Drop a candidate
    # unless >= min_polya5p_reads of its assigned reads BOTH have a krill
    # whole-read polyA tail (qc PASS & length > min_polya_length) AND map
    # with their genomic 5' end within polya5p_window_bp of the candidate's
    # 5' end. GTF candidates are exempt by default (polya5p_exempt_gtf=True,
    # like fusion); only novel candidates are gated. Requires signal; no-ops
    # when absent.
    if (
        _score_filter_on
        and getattr(config, "min_polya5p_reads", 0) > 0
        and config.signal_path
        and Path(config.signal_path).exists()
    ):
        before_keys = set(aggregated)
        aggregated = apply_polya5p_fn(aggregated)
        _record(before_keys - set(aggregated), "polya5p")

    return aggregated, outcomes
