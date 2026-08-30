"""EndpointRefine (EXPERIMENTAL, off by default): split endpoint states.

A surviving novel multi-exon model represents exactly one splice chain with
one TSS/TES pair, but dRNA read ends can show several genuine endpoint states
for the same chain (measured upper bound: 62/4,477 = 1.38% of missed T3
multi-exon truths were endpoint-collapsed chains). This pass splits such a
model into at most ``max_splits`` endpoint states when the read-end evidence
supports distinct modes, with:

* strand-aware 5'/3' semantics ('+': TSS = span start, TES = span end;
  '-' mirrored);
* degradation protection: an alternative TSS mode INTERIOR to the primary
  (downstream on '+', upstream on '-') is the 5'-truncation direction and
  must clear ``degradation_guard`` times the normal support;
* supported endpoint PAIRS: a split exists only when >= ``min_end_reads``
  assigned reads jointly hit its (TSS mode, TES mode) pair within
  ``window_bp``, and carries >= ``min_pair_frac`` of end-mapped reads;
* a hard ``max_splits`` cap (primary state included);
* deterministic stable IDs (BLAKE2b over chain+endpoints, same format as
  discovery IDs);
* MANDATORY post-split requantification: the split plan emits per-read
  routes that the final-survivor abundance refit consumes, so every read's
  mass is re-dealt over the split states and mass conservation is enforced
  by the existing refit invariants. Splitting without requantification is
  not possible through this API.

Poly(A)-supported TES is a v2 hook: ``polya_read_ids`` (reads with detected
poly(A) tails) strengthens TES modes when provided; the current pipeline
wiring passes None, so v1 is signal-free and relies on end-mode sharpness
plus the degradation guard. This limitation is deliberate and documented.
"""
from __future__ import annotations

import hashlib
import logging
from dataclasses import dataclass, field, replace
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

logger = logging.getLogger(__name__)


def _endpoint_id(chrom: str, strand: str, introns: Tuple[Tuple[int, int], ...],
                 start: int, end: int) -> str:
    payload = f"{chrom}:{strand}:{introns}:{start}-{end}".encode()
    return "novel_" + hashlib.blake2b(payload, digest_size=8).hexdigest()


def _introns(qr) -> Tuple[Tuple[int, int], ...]:
    exons = sorted(qr.exons)
    return tuple((exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1))


def _cluster_modes(positions: List[int], window_bp: int) -> List[Tuple[int, int]]:
    """Greedy 1-D clustering -> [(mode_position, support)], deterministic."""
    if not positions:
        return []
    modes: List[List[int]] = []
    for pos in sorted(positions):
        if modes and pos - modes[-1][-1] <= window_bp:
            modes[-1].append(pos)
        else:
            modes.append([pos])
    out = []
    for cluster in modes:
        mid = cluster[len(cluster) // 2] if len(cluster) % 2 else cluster[len(cluster) // 2 - 1]
        out.append((int(mid), len(cluster)))
    return out


@dataclass
class EndpointSplitPlan:
    """Replacements plus the read routes the refit must honor."""

    replacements: Dict[str, List] = field(default_factory=dict)
    read_routes: Dict[str, Dict[str, str]] = field(default_factory=dict)
    primary: Dict[str, str] = field(default_factory=dict)
    details: List[dict] = field(default_factory=list)
    tss_verdicts: Dict[str, List[dict]] = field(default_factory=dict)
    # (parent_id, state_id, evidence|None, skip_reason|None) collected before
    # the run-wide Benjamini-Hochberg correction is applied.
    pending_tss: List[tuple] = field(default_factory=list)

    @property
    def n_splits(self) -> int:
        return sum(len(v) for v in self.replacements.values())


def plan_endpoint_splits(
    results: Mapping[str, object],
    read_ends: Mapping[str, Tuple[int, int]],
    *,
    window_bp: int = 25,
    min_end_reads: int = 3,
    min_pair_frac: float = 0.15,
    max_splits: int = 2,
    degradation_guard: float = 2.0,
    polya_read_ids: Optional[set] = None,
    tss_evidence_mode: str = "off",
    background_hazard: float = 0.0,
) -> EndpointSplitPlan:
    """Compute endpoint splits for novel multi-exon survivors.

    Pure planning - ``results`` is not mutated. Only ``source == "novel"``
    multi-exon candidates with assigned reads are considered.
    """
    plan = EndpointSplitPlan()
    for cid in sorted(results):
        qr = results[cid]
        if getattr(qr, "source", "") != "novel":
            continue
        introns = _introns(qr)
        if not introns:
            continue
        assigned = tuple(getattr(qr, "assigned_read_ids", ()) or ())
        spans = [(rid, read_ends.get(rid)) for rid in assigned]
        spans = [(rid, s) for rid, s in spans if s is not None]
        if len(spans) < 2 * min_end_reads:
            continue

        plus = qr.strand == "+"
        five = [(rid, s[0] if plus else s[1]) for rid, s in spans]
        three = [(rid, s[1] if plus else s[0]) for rid, s in spans]
        n_span = len(spans)

        tss_modes = _cluster_modes([p for _, p in five], window_bp)
        tes_modes = _cluster_modes([p for _, p in three], window_bp)
        if not tss_modes or not tes_modes:
            continue

        # Primary endpoint state = the candidate's own ends.
        primary_tss = qr.start if plus else qr.end
        primary_tes = qr.end if plus else qr.start

        def keep_mode(mode, support, primary_pos, is_tss):
            if support < min_end_reads:
                return False
            if abs(mode - primary_pos) <= window_bp:
                return True  # the primary state itself
            if is_tss:
                interior = mode > primary_pos if plus else mode < primary_pos
                if interior:
                    # 5'-truncation (degradation) direction: demand more.
                    needed = max(min_end_reads,
                                 int(degradation_guard * min_end_reads))
                    return support >= needed
            else:
                if polya_read_ids is not None:
                    return support >= min_end_reads  # polyA handled at pair level
            return True

        tss_keep = [(m, s) for m, s in tss_modes
                    if keep_mode(m, s, primary_tss, True)]
        tes_keep = [(m, s) for m, s in tes_modes
                    if keep_mode(m, s, primary_tes, False)]
        if not tss_keep or not tes_keep:
            continue

        # Route reads to their nearest kept (tss, tes) pair within window.
        pair_support: Dict[Tuple[int, int], List[str]] = {}
        for (rid, p5), (_rid2, p3) in zip(five, three):
            best5 = min(tss_keep, key=lambda ms: (abs(ms[0] - p5), ms[0]))
            best3 = min(tes_keep, key=lambda ms: (abs(ms[0] - p3), ms[0]))
            if abs(best5[0] - p5) > window_bp or abs(best3[0] - p3) > window_bp:
                continue
            pair = (best5[0], best3[0])
            if polya_read_ids is not None and rid not in polya_read_ids:
                # v2 hook: unsupported-TES reads do not vote for non-primary TES.
                if abs(pair[1] - primary_tes) > window_bp:
                    continue
            pair_support.setdefault(pair, []).append(rid)

        qualified = sorted(
            (
                (pair, rids) for pair, rids in pair_support.items()
                if len(rids) >= min_end_reads
                and len(rids) >= min_pair_frac * n_span
            ),
            key=lambda item: (-len(item[1]), item[0]),
        )[:max_splits]
        if len(qualified) < 2:
            continue  # nothing to split - a single state is the status quo

        subs = []
        routes: Dict[str, str] = {}
        exons = sorted(qr.exons)
        seen_ids = set()
        detail_states = []
        for pair, rids in qualified:
            tss, tes = pair
            if plus:
                new_start, new_end = tss, tes
            else:
                new_start, new_end = tes, tss
            # Endpoint sanity: ends must stay outside the intron chain.
            if new_start >= exons[0][1] or new_end <= exons[-1][0]:
                continue
            new_exons = ((new_start, exons[0][1]),) + tuple(exons[1:-1]) + (
                (exons[-1][0], new_end),)
            sub_id = _endpoint_id(qr.chrom, qr.strand, introns,
                                  new_start, new_end)
            if sub_id in seen_ids or sub_id in results:
                continue
            seen_ids.add(sub_id)
            sub = replace(
                qr,
                candidate_id=sub_id,
                start=new_start,
                end=new_end,
                exons=new_exons,
                abundance=float(len(rids)),
                num_assigned_reads=len(rids),
                assigned_read_ids=tuple(sorted(rids)),
            )
            subs.append(sub)
            for rid in rids:
                routes[rid] = sub_id
            detail_states.append({
                "state_id": sub_id,
                "tss": pair[0], "tes": pair[1],
                "pair_reads": len(rids),
                "polya_supported_reads": (
                    sum(1 for r in rids if r in polya_read_ids)
                    if polya_read_ids is not None else None
                ),
            })
        # --- TSS evidence gate (experimental) -----------------------------
        # `audit` records the degradation-null verdict without changing the
        # split; `require` keeps only endpoint states whose alternative TSS
        # beats the local degradation background. An `unidentifiable` verdict
        # NEVER drops a state: insufficient evidence is not evidence of
        # absence, so the model is left as it was.
        if tss_evidence_mode in ("audit", "require") and len(subs) >= 2:
            from fin.analysis.tss_evidence import (
                evaluate_internal_tss,
                evaluate_tes_support,
                genomic_to_offset,
            )

            primary_id = subs[0].candidate_id
            primary_tes_g = subs[0].end if plus else subs[0].start
            for sub in subs:
                if sub.candidate_id == primary_id:
                    continue
                s5 = sub.start if plus else sub.end
                d0 = genomic_to_offset(s5, qr.exons, qr.strand)
                if d0 is None:
                    plan.pending_tss.append(
                        (cid, sub.candidate_id, None, "offset_unmappable"))
                    continue
                # Anchor on the RELIABLE 3' end first: only reads sharing this
                # state's TES may vote on its TSS. Mixing 3' populations
                # inflates the apparent 5' peak.
                sub_tes_g = sub.end if plus else sub.start
                same_tes = abs(sub_tes_g - primary_tes_g) <= window_bp
                primary_5p = subs[0].start if plus else subs[0].end
                d_primary = genomic_to_offset(primary_5p, qr.exons, qr.strand)
                shared_tss = (d_primary is not None
                              and abs(d0 - d_primary) < window_bp)

                if shared_tss and not same_tes:
                    # Rung 2: the state shares the primary's START and differs
                    # only at the 3' end. There is NO internal TSS to detect;
                    # the decision belongs to dRNA's reliable end.
                    three = []
                    for rid in assigned:
                        sp = read_ends.get(rid)
                        if sp is None:
                            continue
                        three.append(sp[1] if plus else sp[0])
                    ev = evaluate_tes_support(
                        candidate_id=sub.candidate_id, parent_id=cid,
                        read_three_prime_offsets=three,
                        candidate_tes_offset=sub_tes_g,
                        parent_tes_offset=primary_tes_g,
                        window_bp=window_bp,
                    )
                else:
                    # Restrict to reads sharing this state's 3' end before
                    # letting them vote on its 5' end.
                    offsets = []
                    for rid in assigned:
                        sp = read_ends.get(rid)
                        if sp is None:
                            continue
                        if not same_tes:
                            r3 = sp[1] if plus else sp[0]
                            if abs(r3 - sub_tes_g) > 2 * window_bp:
                                continue
                        r5 = sp[0] if plus else sp[1]
                        o = genomic_to_offset(r5, qr.exons, qr.strand)
                        if o is not None:
                            offsets.append(o)
                    # Deterministic search-space size: EndpointRefine can only
                    # propose a start inside the FIRST exon, so the number of
                    # bins a peak could have been selected from is that exon's
                    # length in bins -- not the observed read span, which
                    # varies with depth.
                    first_exon = qr.exons[0] if plus else qr.exons[-1]
                    exon_len = max(int(first_exon[1] - first_exon[0]), 0)
                    n_bins = max(exon_len // window_bp + 1, 1)
                    ev = evaluate_internal_tss(
                        candidate_id=sub.candidate_id, parent_id=cid,
                        offsets=offsets, tss_offset=d0,
                        background_hazard=background_hazard,
                        identifiability="tss_only" if same_tes else "own_tes",
                        bin_bp=window_bp,
                        n_eligible_bins=n_bins,
                    )
                plan.pending_tss.append((cid, sub.candidate_id, ev, None))

        if len(subs) < 2:
            continue
        plan.replacements[cid] = subs
        plan.read_routes[cid] = routes
        plan.primary[cid] = subs[0].candidate_id  # largest-support state
        # Route audit: what the refit will actually see. Reads assigned to
        # the original candidate but not explicitly routed fall back to the
        # primary state inside refit_survivor_abundance, so both numbers are
        # needed to reconcile the split against the refit accounting.
        route_counts: Dict[str, int] = {}
        for target in routes.values():
            route_counts[target] = route_counts.get(target, 0) + 1
        plan.details.append({
            "candidate_id": cid,
            "n_end_reads": n_span,
            "_placeholder_bh": None,
            "n_assigned_reads": len(assigned),
            "tss_modes": [list(m) for m in tss_keep],
            "tes_modes": [list(m) for m in tes_keep],
            "polya_available": polya_read_ids is not None,
            "states": detail_states,
            "routed_reads": len(routes),
            "route_counts_by_state": dict(sorted(route_counts.items())),
            "unrouted_to_primary": len(assigned) - len(routes),
            "primary_state": subs[0].candidate_id,
        })

    # --- run-wide locus-level FDR, then enforce `require` -----------------
    # Correction must span every LOCUS in the run: a per-test p<=0.05 over
    # ~150 candidate starts would manufacture ~7 false positives. Alternatives
    # within one locus are combined first (Bonferroni), then BY is applied
    # across loci, which is valid under arbitrary dependence.
    if plan.pending_tss:
        from fin.analysis.tss_evidence import apply_grouped_fdr

        evs = [e for _p, _s, e, _r in plan.pending_tss if e is not None]
        apply_grouped_fdr(evs)
        verdict_by_state: Dict[str, str] = {}
        for parent, state, ev, skip in plan.pending_tss:
            if ev is None:
                row = {"candidate_id": state, "parent_id": parent,
                       "verdict": "unidentifiable", "reason": skip or "skipped"}
            else:
                row = ev.as_row()
            plan.tss_verdicts.setdefault(parent, []).append(row)
            verdict_by_state[state] = row["verdict"]

        if tss_evidence_mode == "require":
            for parent in list(plan.replacements):
                subs = plan.replacements[parent]
                primary = plan.primary.get(parent)
                kept = [s for s in subs
                        if s.candidate_id == primary
                        or verdict_by_state.get(s.candidate_id) != "unsupported"]
                if len(kept) < 2:
                    # every alternative start was refuted: keep the model
                    # unsplit rather than emitting a lone primary state
                    del plan.replacements[parent]
                    plan.read_routes.pop(parent, None)
                    plan.primary.pop(parent, None)
                    continue
                plan.replacements[parent] = kept
                keep_ids = {s.candidate_id for s in kept}
                plan.read_routes[parent] = {
                    r: t for r, t in plan.read_routes.get(parent, {}).items()
                    if t in keep_ids
                }
    for d in plan.details:
        d.pop("_placeholder_bh", None)
    return plan


def apply_endpoint_splits(results: Dict[str, object], plan: EndpointSplitPlan):
    """Replace each split candidate with its endpoint states (new dict)."""
    out = {}
    for cid, qr in results.items():
        if cid in plan.replacements:
            for sub in plan.replacements[cid]:
                out[sub.candidate_id] = sub
        else:
            out[cid] = qr
    return out
