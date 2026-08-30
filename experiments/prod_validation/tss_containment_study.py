#!/usr/bin/env python3
"""SIRV ground-truth study: can we tell a real nested short transcript from a
5'-degradation artifact of the long one?

DESIGN (leakage-guarded; the reference is used ONLY for labels):

POSITIVES  Real SIRV nested pairs (short chain is a contiguous sub-chain of the
           long one, short span inside the long span) where BOTH transcripts
           are expressed in this sample. Scored once per (short transcript,
           locus) so one short model matched to several long parents cannot
           inflate the count.

NEGATIVES  Two independent kinds:
  (a) pseudo-TSS: an offset sampled INSIDE a long transcript at a position
      where truth has no transcript start, depth-matched to the positives;
  (b) degradation-only loci: single-isoform SIRV/ERCC transcripts with no
      nested partner at all -- every internal offset there is a negative.

The discriminator under test is the conditional-termination-hazard test in
`fin.analysis.tss_evidence`. The incumbent baseline is the EndpointRefine-style
heuristic (raw 5'-end mode count with a fixed read floor), so the comparison
answers "is the model actually better than counting modes?".

Reported: PR-AUC / ROC-AUC, FPR at fixed recall, verdict confusion including
the abstention class, and a read-depth stress ladder (down-sampling) that
checks the method ABSTAINS at low depth instead of inventing transcripts.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import random
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).parents[2]))

from fin.analysis.tss_evidence import (  # noqa: E402
    TssEvidence,
    VERDICT_SUPPORTED,
    build_hazard_profile,
    classify_identifiability,
    evaluate_internal_tss,
    evaluate_tes_support,
    pooled_background_hazard,
    read_five_prime_offset,
    spliced_length,
)


# ---------------------------------------------------------------- annotation
def load_transcripts(gtf: Path) -> Dict[str, dict]:
    ex: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    meta: Dict[str, Tuple[str, str, str]] = {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        ex[m.group(1)].append((int(f[3]) - 1, int(f[4])))
        g = re.search(r'gene_id "([^"]+)"', f[8])
        meta[m.group(1)] = (f[0], f[6], g.group(1) if g else m.group(1))
    out = {}
    for tid, exons in ex.items():
        exons = sorted(exons)
        chrom, strand, gene = meta[tid]
        out[tid] = {
            "id": tid, "chrom": chrom, "strand": strand, "gene": gene,
            "exons": exons,
            "start": exons[0][0], "end": exons[-1][1],
            "len": spliced_length(exons),
            "chain": tuple((exons[i][1], exons[i + 1][0])
                           for i in range(len(exons) - 1)),
        }
    return out


def is_subchain(a: Sequence, b: Sequence) -> bool:
    if len(a) == 0:
        return True
    if len(a) > len(b):
        return False
    return any(tuple(b[i:i + len(a)]) == tuple(a)
               for i in range(len(b) - len(a) + 1))


def nested_pairs(tx: Dict[str, dict]) -> List[Tuple[str, str]]:
    pairs = []
    by_key = defaultdict(list)
    for t in tx.values():
        by_key[(t["chrom"], t["strand"])].append(t)
    for key, group in by_key.items():
        for s in group:
            for l in group:
                if s["id"] == l["id"]:
                    continue
                if not (l["start"] <= s["start"] and s["end"] <= l["end"]):
                    continue
                if (s["end"] - s["start"]) >= (l["end"] - l["start"]):
                    continue
                if not is_subchain(s["chain"], l["chain"]):
                    continue
                pairs.append((s["id"], l["id"]))
    return pairs


# ---------------------------------------------------------------------- reads
def read_spans(bam: Path) -> Dict[Tuple[str, str], List[Tuple[int, int]]]:
    """(chrom, strand) -> primary read spans."""
    import pysam

    out = defaultdict(list)
    with pysam.AlignmentFile(str(bam), "rb") as fh:
        for r in fh.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            if r.reference_start is None or r.reference_end is None:
                continue
            strand = "-" if r.is_reverse else "+"
            out[(r.reference_name, strand)].append(
                (int(r.reference_start), int(r.reference_end)))
    return out


def offsets_for_model_tes_grouped(
    model: dict, spans: Sequence[Tuple[int, int]], tes_genomic: int,
    *, tol: int = 50,
) -> List[int]:
    """5' offsets of reads whose 3' end matches a SPECIFIC TES.

    For an `own_tes` candidate the reliable dRNA end must partition the reads
    first: testing a short model's TSS against reads that belong to a
    different 3' end mixes two populations and inflates the apparent peak.
    """
    plus = model["strand"] == "+"
    out = []
    for s, e in spans:
        if e < model["start"] or s > model["end"]:
            continue
        if s < model["start"] - 25 or e > model["end"] + 25:
            continue
        r3 = e if plus else s
        if abs(r3 - tes_genomic) > tol:
            continue
        off = read_five_prime_offset((s, e), model["exons"], model["strand"])
        if off is not None:
            out.append(off)
    return out


def offsets_for_model(model: dict, spans: Sequence[Tuple[int, int]]) -> List[int]:
    """5' offsets of reads compatible with (contained in) the model's span."""
    out = []
    for s, e in spans:
        if e < model["start"] or s > model["end"]:
            continue
        # require the read to sit inside the model's genomic footprint so an
        # unrelated neighbouring gene cannot pollute the hazard
        if s < model["start"] - 25 or e > model["end"] + 25:
            continue
        off = read_five_prime_offset((s, e), model["exons"], model["strand"])
        if off is not None:
            out.append(off)
    return out


def expressed(nanocount: Path, threshold: float) -> set:
    keep = set()
    for row in csv.DictReader(open(nanocount), delimiter="\t"):
        try:
            if float(row["est_count"]) >= threshold:
                keep.add(row["transcript_name"].split(".")[0])
        except (KeyError, TypeError, ValueError):
            continue
    return keep


# ------------------------------------------------------------------- baseline
def heuristic_mode_score(offsets: Sequence[int], tss_offset: int,
                         window: int = 25) -> int:
    """Incumbent: raw count of 5' ends within a window of the candidate TSS."""
    return sum(1 for o in offsets if abs(o - tss_offset) <= window)


# ----------------------------------------------------------------- evaluation
def pr_auc(scored: List[Tuple[float, int]]) -> Tuple[float, float]:
    """(PR-AUC, ROC-AUC), TIE-AWARE.

    Ties must be resolved as a group: advancing positives before negatives at
    an identical score manufactures a perfect AUC. Points are therefore
    grouped by score, the curve advances once per distinct threshold, and the
    ROC contribution of a tied block uses the 0.5 tie credit (trapezoid).
    """
    if not scored:
        return 0.0, 0.0
    P = sum(1 for _, y in scored if y == 1)
    N = len(scored) - P
    if P == 0 or N == 0:
        return 0.0, 0.0
    groups: Dict[float, List[int]] = defaultdict(list)
    for s, y in scored:
        groups[s].append(y)
    tp = fp = 0
    prev_recall = prev_fpr = 0.0
    pr_area = roc_area = 0.0
    for s in sorted(groups, reverse=True):
        blk = groups[s]
        tp += sum(1 for y in blk if y == 1)
        fp += sum(1 for y in blk if y == 0)
        recall, fpr = tp / P, fp / N
        precision = tp / (tp + fp) if (tp + fp) else 1.0
        pr_area += precision * (recall - prev_recall)
        # trapezoid == 0.5 tie credit inside the tied block
        roc_area += 0.5 * (recall + prev_recall) * (fpr - prev_fpr)
        prev_recall, prev_fpr = recall, fpr
    return pr_area, roc_area


def fpr_at_recall(scored: List[Tuple[float, int]], target: float) -> Optional[float]:
    """FPR at the first THRESHOLD reaching the target recall (tie-aware)."""
    P = sum(1 for _, y in scored if y == 1)
    N = len(scored) - P
    if P == 0 or N == 0:
        return None
    groups: Dict[float, List[int]] = defaultdict(list)
    for s, y in scored:
        groups[s].append(y)
    tp = fp = 0
    for s in sorted(groups, reverse=True):
        blk = groups[s]
        tp += sum(1 for y in blk if y == 1)
        fp += sum(1 for y in blk if y == 0)
        if tp / P >= target:
            return fp / N
    return 1.0


# ------------------------------------------------------- simulation (tss_only)
def _simulate_degradation(n: int, tx_len: int, hazard: float,
                          rng: random.Random, tss_offset: int = 0) -> List[int]:
    """Reads from ONE transcript whose 5' end is at ``tss_offset``.

    Walks 3'->5' and terminates with per-base probability ``hazard``; a read
    that never terminates reaches the transcript's own start.
    """
    out = []
    for _ in range(n):
        d = tx_len
        while d > tss_offset:
            if rng.random() < hazard:
                break
            d -= 1
        out.append(max(d, tss_offset))
    return out


def simulation_study(background: float, *, bin_bp: int, seed: int) -> dict:
    """Ground-truth simulation of the HARDEST rung: same chain, same TES,
    internal TSS only.

    SIRV contains no expressed `tss_only` nested pair, so the only way to
    characterise that case with known truth is to generate it. Per-base
    degradation hazard is derived from the measured per-bin background so the
    simulated reads decay like the real sample.
    """
    rng = random.Random(seed)
    per_base = 1.0 - (1.0 - min(max(background, 1e-9), 0.5)) ** (1.0 / bin_bp)
    tx_len, d0 = 3000, 1200
    out = {"per_base_hazard": per_base, "tx_len": tx_len, "tss_offset": d0,
           "cells": []}
    for depth in (10, 20, 50, 100, 200):
        for short_frac in (0.0, 0.1, 0.25, 0.5):
            tp = fp = tn = fn = abst_pos = abst_neg = 0
            for rep in range(40):
                r = random.Random(seed * 1000 + depth * 17 + int(short_frac * 100) * 7 + rep)
                n_short = int(round(depth * short_frac))
                n_long = depth - n_short
                offs = _simulate_degradation(n_long, tx_len, per_base, r)
                offs += _simulate_degradation(n_short, tx_len, per_base, r,
                                              tss_offset=d0)
                ev = evaluate_internal_tss(
                    candidate_id="sim", parent_id="simL", offsets=offs,
                    tss_offset=d0, background_hazard=background,
                    identifiability="tss_only", bin_bp=bin_bp,
                    n_bootstrap=400, seed=r.randrange(10**6),
                )
                truth = short_frac > 0.0
                if ev.verdict == "unidentifiable":
                    if truth:
                        abst_pos += 1
                    else:
                        abst_neg += 1
                elif ev.verdict == VERDICT_SUPPORTED:
                    if truth:
                        tp += 1
                    else:
                        fp += 1
                else:
                    if truth:
                        fn += 1
                    else:
                        tn += 1
            decided_pos = tp + fn
            decided_neg = fp + tn
            out["cells"].append({
                "depth": depth, "short_fraction": short_frac,
                "tp": tp, "fp": fp, "tn": tn, "fn": fn,
                "abstain_pos": abst_pos, "abstain_neg": abst_neg,
                "recall": (tp / decided_pos) if decided_pos else None,
                "fpr": (fp / decided_neg) if decided_neg else None,
            })
    return out


# ----------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True, type=Path)
    ap.add_argument("--bam", required=True, type=Path)
    ap.add_argument("--nanocount", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--expr-threshold", type=float, default=3.0)
    ap.add_argument("--bin-bp", type=int, default=25)
    ap.add_argument("--negatives-per-positive", type=int, default=4)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    tx = load_transcripts(args.gtf)
    spans = read_spans(args.bam)
    expr = expressed(args.nanocount, args.expr_threshold) if args.nanocount else None
    pairs = nested_pairs(tx)
    contained = {s for s, _ in pairs}

    print(f"transcripts={len(tx)} nested_pairs={len(pairs)} "
          f"distinct_contained={len(contained)}")

    # --- degradation background from loci with NO nested partner -----------
    clean_profiles = []
    clean_ids = []
    for tid, t in tx.items():
        if tid in contained:
            continue
        if any(l == tid for _, l in pairs):     # is a parent of something
            continue
        offs = offsets_for_model(t, spans.get((t["chrom"], t["strand"]), []))
        if len(offs) >= 20:
            clean_profiles.append(
                build_hazard_profile(offs, t["len"], bin_bp=args.bin_bp))
            clean_ids.append(tid)
    background = pooled_background_hazard(clean_profiles)
    print(f"degradation background from {len(clean_profiles)} single-isoform "
          f"loci: hazard={background:.6f} per {args.bin_bp}bp bin")

    rows = []
    scored_model: List[Tuple[float, int]] = []
    scored_heur: List[Tuple[float, int]] = []

    # --- POSITIVES: real nested shorts, ONE ROW PER SHORT TRANSCRIPT -------
    # A short model contained in several long parents is ONE biological
    # observation, not several; counting each parent separately would inflate
    # the positive set and the apparent performance.
    from fin.analysis.tss_evidence import genomic_to_offset
    best_by_short: Dict[str, dict] = {}
    for s_id, l_id in pairs:
        if expr is not None:
            if s_id.split(".")[0] not in expr or l_id.split(".")[0] not in expr:
                continue
        s, l = tx[s_id], tx[l_id]
        s5 = s["start"] if s["strand"] == "+" else s["end"]
        d0 = genomic_to_offset(s5, l["exons"], l["strand"])
        if d0 is None:
            continue
        ident = classify_identifiability(s["exons"], l["exons"], l["strand"])
        all_spans = spans.get((l["chrom"], l["strand"]), [])
        plus = l["strand"] == "+"
        # Route to the evidence that can actually decide this rung. When the
        # contained model starts where the parent starts (d0 inside the first
        # bin) there is NO internal TSS to detect, and the decision belongs to
        # dRNA's reliable 3' end.
        if ident == "own_tes" and d0 < args.bin_bp:
            s_tes_g = s["end"] if plus else s["start"]
            l_tes_g = l["end"] if plus else l["start"]
            three = []
            for a, b in all_spans:
                if b < l["start"] or a > l["end"]:
                    continue
                if a < l["start"] - 25 or b > l["end"] + 25:
                    continue
                three.append(b if plus else a)
            ev = evaluate_tes_support(
                candidate_id=s_id, parent_id=l_id,
                read_three_prime_offsets=three,
                candidate_tes_offset=s_tes_g, parent_tes_offset=l_tes_g,
                window_bp=args.bin_bp,
            )
            offs = three
            heur_at = s_tes_g
        else:
            if ident == "own_tes":
                s_tes = s["end"] if plus else s["start"]
                offs = offsets_for_model_tes_grouped(l, all_spans, s_tes)
            else:
                offs = offsets_for_model(l, all_spans)
            ev = evaluate_internal_tss(
                candidate_id=s_id, parent_id=l_id, offsets=offs,
                tss_offset=d0, background_hazard=background,
                identifiability=ident, bin_bp=args.bin_bp,
            )
            heur_at = d0
        row = ev.as_row()
        row.update(label=1, kind="real_nested", dedup_key=s_id,
                   locus_id=l.get("gene") or l_id,
                   heuristic=heuristic_mode_score(offs, heur_at))
        # keep the parent giving the most reads: the best-powered test
        prev = best_by_short.get(s_id)
        if prev is None or row["n_at_risk"] > prev["n_at_risk"]:
            best_by_short[s_id] = row

    for row in best_by_short.values():
        rows.append(row)
        scored_model.append((row["effect_size"], 1))
        scored_heur.append((float(row["heuristic"]), 1))

    n_pos = len(best_by_short)
    print(f"positives (expressed real nested, DEDUPED per short transcript): "
          f"{n_pos}")

    # --- NEGATIVE (a): depth-matched pseudo-TSS inside long transcripts ----
    # A negative position must not coincide with the TSS of ANY real
    # transcript sharing the locus -- not merely the nested partners. The
    # first pass of this study used the narrow definition and 13 of its 15
    # apparent false positives turned out to be the method correctly finding
    # a different real transcript's start, i.e. mislabelled negatives.
    from fin.analysis.tss_evidence import genomic_to_offset
    real_tss_offsets = defaultdict(set)
    for pid, p in tx.items():
        for tid, t in tx.items():
            if tid == pid or t["chrom"] != p["chrom"] or t["strand"] != p["strand"]:
                continue
            t5 = t["start"] if t["strand"] == "+" else t["end"]
            d = genomic_to_offset(t5, p["exons"], p["strand"])
            if d is not None:
                real_tss_offsets[pid].add(d)

    parents = sorted({l for _, l in pairs})
    n_neg_a = 0
    for l_id in parents:
        l = tx[l_id]
        offs = offsets_for_model(l, spans.get((l["chrom"], l["strand"]), []))
        if len(offs) < 20:
            continue
        forbidden = real_tss_offsets.get(l_id, set())
        for _ in range(args.negatives_per_positive):
            for _try in range(40):
                d0 = rng.randrange(int(0.1 * l["len"]), int(0.9 * l["len"]) or 1)
                if all(abs(d0 - f) > 3 * args.bin_bp for f in forbidden):
                    break
            else:
                continue
            ev = evaluate_internal_tss(
                candidate_id=f"pseudo_{l_id}_{d0}", parent_id=l_id,
                offsets=offs, tss_offset=d0, background_hazard=background,
                identifiability="tss_only", bin_bp=args.bin_bp,
            )
            row = ev.as_row()
            row.update(label=0, kind="pseudo_tss",
                       dedup_key=f"pseudo_{l_id}_{d0}",
                       locus_id=l.get("gene") or l_id,
                       heuristic=heuristic_mode_score(offs, d0))
            rows.append(row)
            scored_model.append((ev.effect_size, 0))
            scored_heur.append((float(row["heuristic"]), 0))
            n_neg_a += 1

    # --- NEGATIVE (b): HARD negatives -- the strongest fake TSS available.
    # For each single-isoform transcript (truth says only ONE start exists),
    # take the INTERNAL offset with the largest 5'-end pileup, excluding the
    # transcript's own TSS region. This is precisely the position a counting
    # heuristic is most likely to call a novel short transcript, so it is the
    # fair adversarial negative; random offsets are trivially easy.
    n_neg_b = 0
    for tid in clean_ids:
        t = tx[tid]
        offs = offsets_for_model(t, spans.get((t["chrom"], t["strand"]), []))
        if len(offs) < 20:
            continue
        forbidden_here = real_tss_offsets.get(tid, set())
        counts: Dict[int, int] = defaultdict(int)
        for o in offs:
            if o < 2 * args.bin_bp:       # the transcript's own TSS
                continue
            # never place a "negative" where a real transcript actually starts
            if any(abs(o - f) <= 2 * args.bin_bp for f in forbidden_here):
                continue
            counts[(o // args.bin_bp) * args.bin_bp] += 1
        if not counts:
            continue
        d0 = max(sorted(counts), key=lambda b: counts[b])
        # This offset was DISCOVERED by scanning every eligible bin, so the
        # selection correction must pay for that search. `counts` holds
        # exactly the bins that were scanned.
        ev = evaluate_internal_tss(
            candidate_id=f"hardneg_{tid}_{d0}", parent_id=tid,
            offsets=offs, tss_offset=d0, background_hazard=background,
            identifiability="tss_only", bin_bp=args.bin_bp,
            n_eligible_bins=len(counts),
        )
        row = ev.as_row()
        row.update(label=0, kind="hard_negative_max_pileup",
                   dedup_key=f"hardneg_{tid}_{d0}",
                   locus_id=t.get("gene") or tid,
                   heuristic=heuristic_mode_score(offs, d0))
        rows.append(row)
        scored_model.append((ev.effect_size, 0))
        scored_heur.append((float(row["heuristic"]), 0))
        n_neg_b += 1

    print(f"negatives: pseudo_tss={n_neg_a} hard_max_pileup={n_neg_b}")

    # --- metrics ----------------------------------------------------------
    m_pr, m_roc = pr_auc(scored_model)
    h_pr, h_roc = pr_auc(scored_heur)
    # --- rung-stratified metrics + BH-adjusted verdicts ------------------
    # A single AUC over all rows would mix TES-rule scores (positives) with
    # hazard scores (negatives); those are NOT calibrated to a common scale,
    # so the pooled number is meaningless. Report per rung instead.
    from fin.analysis.tss_evidence import apply_grouped_fdr

    def _locus(row) -> str:
        # The locus is the GTF gene_id, attached explicitly to every row.
        # There is NO fallback: silently using candidate_id would make
        # correlated rows from one gene look like independent evidence, which
        # is exactly the bias this grouping exists to remove.
        lid = row.get("locus_id")
        if not lid:
            raise ValueError(
                f"row {row.get('candidate_id')!r} has no locus_id; every row "
                "must carry its GTF gene_id for locus-level aggregation")
        return str(lid)

    by_rung: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for r in rows:
        rung = r.get("identifiability") or "unknown"
        by_rung[rung][f"label{r['label']}_{r['verdict']}"] += 1
    rung_report = {k: dict(v) for k, v in sorted(by_rung.items())}

    # Locus-grouped: one vote per (locus, label, rung); a locus counts as
    # `supported` if ANY of its candidate rows was supported. Row-level counts
    # treat correlated candidates from one locus as independent evidence.
    loci: Dict[tuple, set] = defaultdict(set)
    for r in rows:
        key = (r.get("identifiability") or "unknown", r["label"], _locus(r))
        loci[key].add(r["verdict"])
    locus_report: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for (rung, label, _lid), verdicts in loci.items():
        if VERDICT_SUPPORTED in verdicts:
            v = VERDICT_SUPPORTED
        elif "unsupported" in verdicts:
            v = "unsupported"
        else:
            v = "unidentifiable"
        locus_report[rung][f"label{label}_{v}"] += 1
    locus_report = {k: dict(v) for k, v in sorted(locus_report.items())}

    # Raw (per-test) verdicts are what the table above reports. Recompute the
    # same rows under run-wide BH so both are on the record.
    adjusted_objs = []
    for r in rows:
        e = TssEvidence(
            candidate_id=str(r["candidate_id"]), parent_id=_locus(r),
            tss_offset=int(r["tss_offset"]),
            identifiability=str(r["identifiability"]),
            verdict=str(r["verdict"]), reason=str(r["reason"]),
        )
        pv = r.get("p_value")
        e.p_value = None if pv is None else float(pv)
        e.selection_corrected = bool(r.get("selection_corrected"))
        adjusted_objs.append(e)
    apply_grouped_fdr(adjusted_objs)
    # AGGREGATION RULE (documented, applied consistently):
    #   * transcript-row view  = one entry per candidate row;
    #   * locus view           = one entry per (rung, label, gene_id), where a
    #                            locus counts as `supported` if ANY of its rows
    #                            is supported, else `unsupported` if any row is
    #                            unsupported, else `unidentifiable`.
    adjusted_counts: Dict[str, int] = defaultdict(int)
    for r, e in zip(rows, adjusted_objs):
        r["verdict_adjusted"] = e.verdict
        r["q_value"] = "" if e.q_value is None else e.q_value
        r["p_within"] = "" if e.p_within is None else e.p_within
        adjusted_counts[f"label{r['label']}_{e.verdict}"] += 1

    adj_loci: Dict[tuple, set] = defaultdict(set)
    for r in rows:
        adj_loci[(r.get("identifiability") or "unknown", r["label"],
                  _locus(r))].add(r["verdict_adjusted"])
    adj_locus_counts: Dict[str, int] = defaultdict(int)
    for (_rung, label, _lid), vs in adj_loci.items():
        v = (VERDICT_SUPPORTED if VERDICT_SUPPORTED in vs
             else "unsupported" if "unsupported" in vs else "unidentifiable")
        adj_locus_counts[f"label{label}_{v}"] += 1

    summary = {
        "metrics_by_rung_raw": rung_report,
        "metrics_by_rung_locus_grouped": locus_report,
        "verdicts_after_locus_fdr_rows": dict(sorted(adjusted_counts.items())),
        "verdicts_after_locus_fdr_loci": dict(sorted(adj_locus_counts.items())),
        "aggregation_rule": (
            "rows = one entry per candidate; loci = one entry per "
            "(rung,label,gene_id), supported if ANY row supported, "
            "else unsupported if any row unsupported, else "
            "unidentifiable"),
        "auc_removed_because": (
            "A pooled PR/ROC would mix evidence rungs (3'-end rule scores for "
            "own_tes positives vs 5' hazard scores for tss_only negatives), "
            "which are not calibrated to a common scale. SIRV has no "
            "expressed tss_only positive, so no valid model-vs-heuristic AUC "
            "exists on this data; the verdict tables are the result."
        ),
        "schema_version": 1,
        "background_hazard": background,
        "bin_bp": args.bin_bp,
        "n_transcripts": len(tx),
        "n_nested_pairs": len(pairs),
        "n_positives": n_pos,
        "n_neg_pseudo_tss": n_neg_a,
        "n_neg_degradation_only": n_neg_b,
    }

    # verdict confusion including abstention
    conf = defaultdict(int)
    for r in rows:
        conf[f"label{r['label']}_{r['verdict']}"] += 1
    summary["verdict_confusion"] = dict(sorted(conf.items()))
    # per-kind verdicts so the hard negatives are visible separately
    bykind = defaultdict(int)
    for r in rows:
        bykind[f"{r['kind']}::{r['verdict']}"] += 1
    summary["verdict_by_kind"] = dict(sorted(bykind.items()))
    summary["simulation_tss_only"] = simulation_study(
        background, bin_bp=args.bin_bp, seed=args.seed)
    ident = defaultdict(int)
    for r in rows:
        if r["label"] == 1:
            ident[r["identifiability"]] += 1
    summary["positive_identifiability"] = dict(sorted(ident.items()))

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "tss_containment_rows.tsv").open("w") as fh:
        cols = list(rows[0]) if rows else []
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")
    (args.out / "tss_containment_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, default=float) + "\n")

    print(json.dumps(summary, indent=2, sort_keys=True, default=float))


if __name__ == "__main__":
    main()
