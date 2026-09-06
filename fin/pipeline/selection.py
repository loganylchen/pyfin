"""Selection layer: decide which candidates survive, in one place.

Historically "which candidates are real" had no single home — it was smeared across three
incompatible mechanisms: pre-EM gates mutating the candidate list, a post-EM ``drop_cols`` column
set buried inside ``_quant_m2_em``, and finalize-time dict deletion. This module is the single home
for the post-assignment INTERVAL scope: the cascade that runs within one
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
handling. PRE_ASSIGN (canonical/junction-dominance) and GLOBAL filters are also
owned here.
"""
from __future__ import annotations

import logging
from collections import defaultdict
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

    recorded: set = set()

    def _record_drops(cols, reason: str) -> None:
        # One decisive outcome per candidate: a column already folded/dropped
        # by an earlier step is not re-recorded by a later gate.
        for j in sorted(cols):
            if j in recorded:
                continue
            recorded.add(j)
            outcomes.append(SelectionOutcome(
                cand_list[j].candidate_id, "drop", reason,
            ))

    # --- Extent-pair resolution. An annotated model whose ends the reads
    #     contradicted was emitted twice: with the annotated span and with the
    #     reads' span (TranscriptCandidate.extent_of). Those are two estimates
    #     of ONE model's extent parameter, never two transcripts, so exactly one
    #     survives: the one the posterior gave more mass, with the other's reads
    #     and mass pooled into it. Resolved FIRST so every later gate sees the
    #     model once, with all its evidence, as it would have without the
    #     hypothesis. No siblings (the default path) -> nothing happens;
    #     otherwise R / ties / nlls are REPLACED by pooled copies so every
    #     mass- or support-based gate below reads the merged column. ---
    folded, R, ties_by_read, nlls_by_read = resolve_extent_pairs(
        cand_list, R, ties_by_read, nlls_by_read, quant_results, drop_cols, outcomes)
    if folded:
        recorded.update(drop_cols)
        logger.info(
            "m2_em interval %s: resolved %d annotated extent pair(s)",
            interval.region_string, folded,
        )

    # --- P2 containment decision layer (opt-in). Evaluate every contained
    #     short-in-long candidate pair with the degradation-aware TES/TSS
    #     tests BEFORE any drop, while both hypotheses and their EM read
    #     assignments still exist. audit: record only (byte-inert).
    #     require: 'contradicted' shorts (independently supported parent +
    #     rejecting test) are dropped here; 'supported' shorts are exempt
    #     from the competition gates below; abstention never acts. ---
    containment_exempt_ids: set = set()
    _cmode = getattr(config, "containment_evidence_mode", "off")
    # Hypotheses already folded away (an extent alternative) are not evidence
    # for the opt-in decision layers below: they compare live models only.
    live_cands = [c for j, c in enumerate(cand_list) if j not in drop_cols]
    live_qrs = [q for j, q in enumerate(quant_results) if j not in drop_cols]
    if _cmode != "off":
        from fin.analysis.containment_evidence import (
            containment_interval_actions,
            evaluate_containment,
            write_containment_tsv,
        )

        c_rows = evaluate_containment(
            live_cands, live_qrs,
            getattr(candidate_set, "read_spans", None) or {},
            interval_name=interval.region_string,
        )
        c_drop_ids, containment_exempt_ids = containment_interval_actions(
            c_rows, _cmode)
        verdict_by_id = {r.short_id: r.verdict for r in c_rows}
        for q in quant_results:
            v = verdict_by_id.get(q.candidate_id)
            if v:
                q.containment_verdict = v
        if c_rows:
            shard_dir = Path(config.work_dir) / "_containment_shards"
            shard_dir.mkdir(parents=True, exist_ok=True)
            name = (f"{interval.chrom}_{interval.start}_{interval.end}_"
                    f"{interval.strand or 'u'}.tsv")
            write_containment_tsv(shard_dir / name, c_rows)
        if c_drop_ids:
            cols = {
                j for j, c in enumerate(cand_list)
                if c.candidate_id in c_drop_ids
            }
            drop_cols |= cols
            _record_drops(cols, "containment_contradicted_tes")
            logger.info(
                "m2_em interval %s: containment TES evidence contradicted "
                "%d contained model(s)", interval.region_string, len(cols),
            )

    # --- P3 signal-ballot junction election (opt-in). Nearby junction
    #     variants with a cigar-level screen imbalance are put to a raw-
    #     signal per-read vote (the validated M2 wobble NLL over contested
    #     windows). Supermajority + quorum elect; 'defeated' chains drop in
    #     require mode; every insufficient outcome retains. audit only
    #     writes the artifact (byte-inert). ---
    _jmode = getattr(config, "junction_election_mode", "off")
    if _jmode != "off":
        from fin.analysis.junction_election import (
            run_junction_election,
            write_election_tsv,
        )

        je_drop_ids, je_rows = run_junction_election(
            config, live_cands, live_qrs, candidate_set,
            interval_name=interval.region_string,
        )
        if je_rows:
            shard_dir = Path(config.work_dir) / "_jelection_shards"
            shard_dir.mkdir(parents=True, exist_ok=True)
            name = (f"{interval.chrom}_{interval.start}_{interval.end}_"
                    f"{interval.strand or 'u'}.tsv")
            write_election_tsv(shard_dir / name, je_rows)
        if je_drop_ids:
            cols = {
                j for j, c in enumerate(cand_list)
                if c.candidate_id in je_drop_ids
            }
            drop_cols |= cols
            _record_drops(cols, "junction_election_defeated")
            logger.info(
                "m2_em interval %s: junction election defeated %d "
                "candidate(s)", interval.region_string, len(cols),
            )

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
            if containment_exempt_ids:
                _d = {
                    j for j in _d
                    if cand_list[j].candidate_id
                    not in containment_exempt_ids
                }
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
        # Opt-in complete-chain rescue: direct structural evidence beats the
        # per-read M1/M2 contest. A candidate (novel OR gtf) whose EVERY
        # junction is cleanly attested by >= K exact complete-chain reads
        # cannot be starved out of the per-read contests - e.g. an
        # endpoint-truncated annotation entry whose intron chain is correct
        # loses every alignment-score tie to full-length competitors, and a
        # true novel chain loses them to a corrupt merged annotation. The
        # per-junction weakest-depth metric preserves the gate's
        # anti-corruption purpose: junction-corrupt (jitter/flip) chains have
        # NO exact complete-chain reads, and misalignment-born variant chains
        # carry the indel signature at their spurious junction (measured
        # weakest depth 0-1 for all 17 junk chains vs >=2 for 12/13 genuine
        # chains).
        if getattr(config, "m2_support_gate_gtf_exempt", False) and _d:
            from fin.scoring.m2_junction_nll import uncontradicted_gtf_columns
            if observed_jct is None:
                observed_jct = observed_junctions_fn(interval)
            exempt = uncontradicted_gtf_columns(
                cand_list, observed_jct,
                wobble_bp=int(getattr(config, "chain_cluster_wobble_bp", 6)))
            if exempt & _d:
                logger.info(
                    "m2_em interval %s: %d annotated candidate(s) exempt from the "
                    "support gate (no junction contradicted by the reads)",
                    interval.region_string, len(exempt & _d))
                _d = _d - exempt
        _d, n_rescued = apply_complete_chain_rescue(
            _d, cand_list,
            int(getattr(config, "complete_chain_rescue_reads", 0)),
        )
        if n_rescued:
            logger.info(
                "m2_em interval %s: complete-chain rescue kept %d "
                "support-gate candidates", interval.region_string, n_rescued,
            )
        if containment_exempt_ids:
            before_exempt = len(_d)
            _d = {
                j for j in _d
                if cand_list[j].candidate_id not in containment_exempt_ids
            }
            if len(_d) != before_exempt:
                logger.info(
                    "m2_em interval %s: containment evidence kept %d "
                    "supported contained model(s) at the support gate",
                    interval.region_string, before_exempt - len(_d),
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
        if containment_exempt_ids:
            _d = {
                j for j in _d
                if cand_list[j].candidate_id not in containment_exempt_ids
            }
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


# Single-junction NOVEL chains are the classic RT/template-switching artifact
# class: on the SIRV4 grid a 1-junction artifact reached weakest depth 8
# while genuine merge-freed 1-junction chains span 4-34, so no low threshold
# separates them. This higher bar (dev-fitted on the six H9 samples; kills
# the observed artifact, keeps 7/12 genuine single-junction recoveries)
# applies ONLY to novel single-junction chains. Independent validation is
# required before any default promotion.
COMPLETE_CHAIN_RESCUE_NOVEL_MONO_JUNC_MIN = 10


def resolve_extent_pairs(cand_list, R, ties_by_read, nlls_by_read, quant_results,
                         drop_cols: set, outcomes: List[SelectionOutcome]):
    """Keep one of {annotated span, read span} per annotated model; fold the other.

    The winner is the column with more EM mass (ties -> the annotated span, so
    an uninformative posterior leaves the annotation exactly as it was). The
    loser's hard reads and soft mass move to the winner, its column mass, tie
    membership and junction NLLs are re-pointed at the winner (so every later
    gate that reads ``R`` / ``ties_by_read`` / ``nlls_by_read`` sees ONE model
    with all its evidence), and the loser joins ``drop_cols`` with a ``fold``
    outcome. Returns ``(n_pairs, R, ties_by_read, nlls_by_read)``; the inputs
    are returned untouched when there is no pair.
    """
    col = {c.candidate_id: j for j, c in enumerate(cand_list)}
    pairs = []
    for j, c in enumerate(cand_list):
        parent = getattr(c, "extent_of", None)
        if parent and parent in col and col[parent] not in drop_cols and j not in drop_cols:
            pairs.append((col[parent], j))
    if not pairs:
        return 0, R, ties_by_read, nlls_by_read
    R = np.array(R, dtype=float, copy=True)
    mass = R.sum(axis=0)
    qr_by_id = {q.candidate_id: q for q in quant_results}
    remap: Dict[int, int] = {}
    n = 0
    for i, j in pairs:
        win, lose = (j, i) if mass[j] > mass[i] else (i, j)
        wq = qr_by_id.get(cand_list[win].candidate_id)
        lq = qr_by_id.get(cand_list[lose].candidate_id)
        if wq is None or lq is None:
            continue
        union = set(wq.assigned_read_ids) | set(lq.assigned_read_ids)
        n_w, n_l = wq.num_assigned_reads, lq.num_assigned_reads
        if n_w + n_l > 0:                       # read-weighted mean assignment probability
            wq.confidence = (wq.confidence * n_w + lq.confidence * n_l) / (n_w + n_l)
        wq.assigned_read_ids = tuple(sorted(union))
        wq.num_assigned_reads = len(union)
        wq.abundance += lq.abundance
        wq.max_R = max(wq.max_R, lq.max_R)
        R[:, win] += R[:, lose]
        R[:, lose] = 0.0
        remap[lose] = win
        drop_cols.add(lose)
        outcomes.append(SelectionOutcome(
            cand_list[lose].candidate_id, "fold", "extent_alternative",
            fold_into=cand_list[win].candidate_id,
        ))
        n += 1
    if remap:
        new_ties: Dict[int, List[int]] = {}
        for r, cols in (ties_by_read or {}).items():
            seen: List[int] = []
            for c in cols:
                c = remap.get(c, c)
                if c not in seen:
                    seen.append(c)
            new_ties[r] = seen
        new_nlls: Dict[int, Dict[int, float]] = {}
        for r, d in (nlls_by_read or {}).items():
            m: Dict[int, float] = {}
            for c, v in d.items():
                c = remap.get(c, c)
                m[c] = v if c not in m else min(m[c], v)
            new_nlls[r] = m
        ties_by_read, nlls_by_read = new_ties, new_nlls
    return n, R, ties_by_read, nlls_by_read


def complete_chain_rescue_ok(
    source: str, n_junctions: int, weakest_depth: int, min_reads: int
) -> bool:
    """Shared rescue policy for the support-gate and soft-mass exemptions.

    - gtf: weakest depth >= min_reads (a chain-correct endpoint-truncated
      annotation entry qualifies; junction-corrupt annotation cannot - no
      exact complete-chain read attests a corrupt junction);
    - novel, >=2 junctions: weakest depth >= min_reads + 1 (one extra read
      over the annotation bar for structures with no annotation prior);
    - novel, single junction: weakest depth >=
      COMPLETE_CHAIN_RESCUE_NOVEL_MONO_JUNC_MIN (artifact class, above);
    - anything else (fusion, mono): never rescued.
    """
    if min_reads <= 0 or n_junctions <= 0 or weakest_depth <= 0:
        return False
    if source == "gtf":
        return weakest_depth >= min_reads
    if source == "novel":
        if n_junctions >= 2:
            return weakest_depth >= min_reads + 1
        return weakest_depth >= COMPLETE_CHAIN_RESCUE_NOVEL_MONO_JUNC_MIN
    return False


def apply_complete_chain_rescue(
    drop_cols: set, cand_list, min_reads: int
) -> Tuple[set, int]:
    """Remove from ``drop_cols`` every candidate accepted by
    :func:`complete_chain_rescue_ok` (weakest-junction clean attestation by
    exact complete-splice-chain reads; see
    TranscriptCandidate.complete_chain_weakest_indel_clean_depth). ``min_reads <= 0``
    disables the rescue and returns ``drop_cols`` unchanged (byte-identical
    default path). Returns (filtered_drop_cols, n_rescued)."""
    if min_reads <= 0 or not drop_cols:
        return drop_cols, 0
    kept = {
        col
        for col in drop_cols
        if not complete_chain_rescue_ok(
            getattr(cand_list[col], "source", ""),
            len(getattr(
                getattr(cand_list[col], "intron_chain", None),
                "introns", (),
            ) or ()),
            getattr(cand_list[col], "complete_chain_weakest_indel_clean_depth", 0),
            min_reads,
        )
    }
    return kept, len(drop_cols) - len(kept)


def unify_extent_pairs(
    interval_results: List[List["object"]],
) -> Tuple[Dict[str, str], List[SelectionOutcome]]:
    """The "one model, one extent" invariant ACROSS intervals, decided before
    aggregation. An annotation that spans two intervals can keep its annotated
    span in one and elect its read-extent sibling in the other. The pair is
    resolved once on total abundance (ties -> the annotation) and the loser's
    per-interval results are RELABELLED to the winner's id in place, so
    :func:`aggregate_across_intervals` deduplicates their reads and soft mass
    under one identity with its ordinary per-read policy (no mass is invented
    or approximated here). Returns ``({loser_id: winner_id}, fold outcomes)``;
    the redirects also route the interval ledgers in the survivor refit.
    """
    mass: Dict[str, float] = defaultdict(float)
    parent_of: Dict[str, str] = {}
    geometry: Dict[str, tuple] = {}
    for results in interval_results:
        for qr in results:
            mass[qr.candidate_id] += float(qr.abundance)
            ext = getattr(qr, "extent_of", None)
            if ext:
                parent_of[qr.candidate_id] = ext
            elif qr.candidate_id not in geometry:
                geometry[qr.candidate_id] = (qr.start, qr.end, tuple(qr.exons))
    redirects: Dict[str, str] = {}
    outcomes: List[SelectionOutcome] = []
    for sib, parent in sorted(parent_of.items()):
        if parent not in mass:
            continue                                   # the pair never split across intervals
        sib_wins = mass[sib] > mass[parent]
        loser, winner = (parent, sib) if sib_wins else (sib, parent)
        if not sib_wins and parent not in geometry:
            continue
        for results in interval_results:
            for qr in results:
                if qr.candidate_id != loser:
                    continue
                qr.candidate_id = winner
                if sib_wins:
                    qr.extent_of = parent            # aggregate takes the union extent
                else:
                    qr.extent_of = None
                    qr.start, qr.end, qr.exons = geometry[parent]
        redirects[loser] = winner
        outcomes.append(SelectionOutcome(
            loser, "fold", "extent_alternative", scope="global", fold_into=winner))
    return redirects, outcomes


def strip_zero_mass_gtf_survivors(
    aggregated: Dict[str, "object"], *, enabled: bool
) -> Tuple[Dict[str, "object"], List[str]]:
    """P0 strict-positive GTF gate, post-refit enforcement.

    The A4 gate tests pre-refit EM soft mass, but the post-selection survivor
    renormalization can reassign every read away from a GTF model and leave an
    exact-zero survivor. Zero abundance here means zero retained read mass
    (the refit conserves mass), so removing the row is mass- and read-neutral
    and never changes any other survivor's abundance. Returns the (possibly
    same) dict and the sorted dropped candidate ids (reason: gtf_zero_mass).
    """
    if not enabled:
        return aggregated, []
    zero_ids = sorted(
        cid
        for cid, qr in aggregated.items()
        if qr.source == "gtf" and qr.abundance <= 0.0
    )
    if not zero_ids:
        return aggregated, []
    zero_set = set(zero_ids)
    kept = {
        cid: qr for cid, qr in aggregated.items() if cid not in zero_set
    }
    return kept, zero_ids


def enforce_gtf_positive_after_refit(
    aggregated: Dict[str, "object"],
    diagnostics: Optional[dict],
    *,
    enabled: bool,
) -> Tuple[Dict[str, "object"], List[str]]:
    """Apply the post-refit strict-positive GTF sweep and keep diagnostics
    consistent with the emitted survivor set.

    ``refit_survivor_abundance`` reports ``surviving_candidates`` for its
    input set; when the sweep then removes zero-mass GTF survivors, the
    written diagnostics must record the drop explicitly so the JSON always
    matches the emitted output. Fields are added only when the gate is
    enabled so default runs stay byte-identical.
    """
    aggregated, zero_ids = strip_zero_mass_gtf_survivors(
        aggregated, enabled=enabled
    )
    if enabled and diagnostics is not None:
        diagnostics["post_refit_zero_gtf_dropped"] = len(zero_ids)
        diagnostics["post_refit_zero_gtf_ids"] = list(zero_ids)
        diagnostics["output_surviving_candidates"] = len(aggregated)
    return aggregated, zero_ids


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
    # P0 strict-positive gate (require_gtf_positive_mass): a GTF-source model
    # must additionally carry EM soft-mass > 0, even when every floor is 0 —
    # zero-mass catalog entries are unexpressed pass-through, not evidence.
    floor_gtf = config.floor_gtf_abundance
    gtf_floor = (
        max(config.min_gtf_abundance, config.min_abundance)
        if floor_gtf
        else config.min_gtf_abundance
    )
    require_gtf_positive = bool(
        getattr(config, "require_gtf_positive_mass", False)
    )
    if _score_filter_on and (
        config.min_abundance > 0.0
        or gtf_floor > 0.0
        or require_gtf_positive
    ):
        before = len(aggregated)
        strict_novel = bool(
            getattr(config, "strict_novel_abundance_floor", False)
        )
        kept = {
            cid: qr
            for cid, qr in aggregated.items()
            if qr.source == "fusion"
            or (
                (
                    qr.abundance >= gtf_floor
                    and (not require_gtf_positive or qr.abundance > 0.0)
                )
                if qr.source == "gtf"
                else (
                    qr.abundance > config.min_abundance
                    if strict_novel
                    else qr.abundance >= config.min_abundance
                )
            )
        }
        dropped = before - len(kept)
        if dropped:
            dropped_ids = set(aggregated) - set(kept)
            zero_mass_ids = {
                cid
                for cid in dropped_ids
                if require_gtf_positive
                and aggregated[cid].source == "gtf"
                and aggregated[cid].abundance <= 0.0
                and aggregated[cid].abundance >= gtf_floor
            }
            if zero_mass_ids:
                _record(zero_mass_ids, "gtf_zero_mass")
            _record(dropped_ids - zero_mass_ids, "min_abundance")
            logger.info(
                "Dropped %d transcripts failing abundance floor "
                "(keep novel %s %.3f, gtf >= %.3f%s)",
                dropped,
                ">" if strict_novel else ">=",
                config.min_abundance,
                gtf_floor,
                " and > 0" if require_gtf_positive else "",
            )
        aggregated = kept


    # Minimum isoform fraction (locus-relative abundance) filter. Drops
    # NOVEL multi-exon transcripts whose abundance is below
    # `min_isoform_fraction` of the dominant transcript in its persisted
    # splice family (or overlap fallback) — the standard Cufflinks heuristic /
    # StringTie -f heuristic (the low-fraction tail is enriched for
    # incompletely-spliced precursors and assembly artifacts).
    # gtf/fusion/mono candidates are EXEMPT. SIRV WARNING: the conservative
    # literature default (0.01) is used; SIRV's F1-optimal ~0.4 is overfit.
    if _score_filter_on and config.min_isoform_fraction > 0.0:
        drop_ids = isoform_fraction_drops(
            aggregated,
            config.min_isoform_fraction,
            locus=getattr(config, "isoform_fraction_locus", "family"),
        )
        if drop_ids:
            aggregated = {
                cid: qr
                for cid, qr in aggregated.items()
                if cid not in drop_ids
            }
            _record(drop_ids, "isoform_fraction")
            logger.info(
                "Dropped %d novel transcripts with isoform fraction < %.3f "
                "(%s locus)",
                len(drop_ids),
                config.min_isoform_fraction,
                getattr(config, "isoform_fraction_locus", "family"),
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
        # Opt-in complete-chain rescue (mirrors the interval support-gate
        # rescue): a candidate whose every junction is cleanly attested by
        # >= K exact complete-chain reads is direct structural evidence, not
        # EM-borrowed soft mass - a competing corrupt annotation model
        # inflates the true chain's soft/hard ratio, and dropping the true
        # chain for that competition artifact is wrong. (The soft-mass
        # filter itself only ever drops novel candidates.)
        _ccr = int(getattr(config, "complete_chain_rescue_reads", 0))
        _cmode_g = getattr(config, "containment_evidence_mode", "off")
        if (_ccr > 0 or _cmode_g == "require") and drop_ids:
            rescued = {
                cid
                for cid in drop_ids
                if complete_chain_rescue_ok(
                    getattr(aggregated[cid], "source", ""),
                    max(len(getattr(aggregated[cid], "exons", ()) or ()) - 1, 0),
                    getattr(
                        aggregated[cid], "complete_chain_weakest_indel_clean_depth", 0
                    ),
                    _ccr,
                )
                or (
                    _cmode_g == "require"
                    and getattr(aggregated[cid], "containment_verdict", "")
                    == "supported"
                )
            }
            if rescued:
                drop_ids = drop_ids - rescued
                logger.info(
                    "Complete-chain rescue kept %d soft-mass candidates",
                    len(rescued),
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


def canonical_gate_select(config, candidate_set, chrom_seq: str) -> List[SelectionOutcome]:
    """PRE_ASSIGN Stage B: drop NOVEL multi-exon candidates with non-canonical junctions.

    Byte-identical move of the runner's former ``_apply_canonical_gate`` body. Mutates
    ``candidate_set.candidates`` in place, keeping a candidate iff it is GTF/fusion (EXEMPT),
    single-exon (trivially canonical), or every internal junction motif is canonical. ``read_ids``
    is left untouched (reads off a dropped candidate re-compete in quant). Returns the drop outcomes.
    """
    from fin.candidates.canonical import chain_all_canonical, parse_motifs

    outcomes: List[SelectionOutcome] = []
    if not chrom_seq:
        # No genome sequence → cannot evaluate motifs; skip the gate rather
        # than drop everything.
        return outcomes
    motif_set = parse_motifs(config.canonical_motifs)
    kept: List = []
    dropped = 0
    for c in candidate_set.candidates:
        if c.source in ("gtf", "fusion"):
            kept.append(c)
            continue
        introns = c.intron_chain.introns
        if not introns:
            kept.append(c)  # mono: no internal junction
            continue
        if chain_all_canonical(introns, chrom_seq, c.strand, motif_set):
            kept.append(c)
        else:
            dropped += 1
            outcomes.append(SelectionOutcome(
                c.candidate_id, "drop", "canonical_gate", scope="pre_assign",
            ))
    if dropped:
        logger.info(
            "Canonical gate: dropped %d/%d novel candidates (interval %s)",
            dropped,
            len(candidate_set.candidates),
            candidate_set.interval.region_string,
        )
    candidate_set.candidates = kept
    return outcomes


def junction_dominance_select(
    config, candidate_set, interval, observed_junctions_fn: Callable
) -> List[SelectionOutcome]:
    """PRE_ASSIGN junction-first gate: drop NOVEL candidates with a weak or non-dominant junction.

    Byte-identical move of the runner's former ``_apply_junction_dominance_gate`` body. Builds the
    directly-observed read junctions (via ``observed_junctions_fn``, passed at call time so the test
    seam survives) and drops a novel multi-exon candidate whose junction has < min_reads reads or is
    dominated by a stronger DIFFERENT junction within the window. Mutates ``candidate_set.candidates``
    in place; ``read_ids`` untouched. Returns the drop outcomes.
    """
    from fin.scoring.m2_junction_nll import dominant_junction_drops

    outcomes: List[SelectionOutcome] = []
    observed = observed_junctions_fn(interval)
    drops = dominant_junction_drops(
        candidate_set.candidates, observed,
        min_reads=int(config.junction_dominance_min_reads),
        window=int(config.junction_dominance_window_bp),
        tol=int(getattr(config, "junction_dominance_tol_bp", 2)),
    )
    if drops:
        kept = [c for i, c in enumerate(candidate_set.candidates)
                if i not in drops]
        for i, c in enumerate(candidate_set.candidates):
            if i in drops:
                outcomes.append(SelectionOutcome(
                    c.candidate_id, "drop", "junction_dominance", scope="pre_assign",
                ))
        logger.info(
            "Junction-dominance gate: dropped %d/%d novel candidates "
            "(interval %s)", len(drops), len(candidate_set.candidates),
            interval.region_string,
        )
        candidate_set.candidates = kept
    return outcomes
