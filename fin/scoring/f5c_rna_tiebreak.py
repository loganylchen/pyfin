"""f5c_rna-based tiebreak: resolve M1 ties using in-memory eventalign match score.

For ambiguous reads (R_best < threshold):
  1. For each tied candidate, run f5c_rna align with candidate's mapped region
  2. Compute match score = matched_events / max_events_across_candidates
     - matched = events with |z| < 2 AND model_kmer != NNNNN
     - FAIL (status!=0) → matched = 0
  3. Redistribute R proportionally by match scores, preserving original row sum
"""
from __future__ import annotations

import logging
from typing import Dict, List, Optional

import mappy
import numpy as np

from fin.candidates.dataclasses import TranscriptCandidate

logger = logging.getLogger(__name__)


def _build_m2_f5c_rna(
    read_ids: List[str],
    read_seqs: Dict[str, str],
    candidates: List,
    signal_path: str,
    pore: str = "rna002",
    as_distance: bool = True,
) -> np.ndarray:
    """Build (n_reads, n_cands) matrix from f5c_rna match scores.

    If as_distance=True: distance = 1 - (matched/max_ev). Lower = better.
    If as_distance=False: returns row-normalized match scores (sum=1 per row).
        Can be used directly as R matrix.
    FAIL → 0.
    """
    import krill as f5c_rna

    n_reads = len(read_ids)
    n_cands = len(candidates)
    dist = np.ones((n_reads, n_cands), dtype=np.float32)

    # Build mappy aligners
    cand_aligners = {}
    for j, c in enumerate(candidates):
        seq = c.sequence if hasattr(c, 'sequence') else None
        if seq and len(seq) > 20:
            a = mappy.Aligner(seq=seq, preset="map-ont")
            if a:
                cand_aligners[j] = a

    # Build per-read mapped sequences for all reads at once
    reads_variants = {}   # {rid: {cand_id: seq}}
    rid_cand_idx = {}     # {rid: {cand_id: j}}

    for i, rid in enumerate(read_ids):
        rseq = read_seqs.get(rid, "")
        if not rseq:
            continue
        per_read = {}
        per_read_j = {}
        for j, a in cand_aligners.items():
            best = None
            for hit in a.map(rseq):
                if best is None or hit.mlen > best.mlen:
                    best = hit
            if best:
                cid = candidates[j].candidate_id
                per_read[cid] = candidates[j].sequence[best.r_st:best.r_en]
                per_read_j[cid] = j
        if per_read:
            reads_variants[rid] = per_read
            rid_cand_idx[rid] = per_read_j

    if not reads_variants:
        return dist

    # Try batch API first (align_reads_variants, plural), fall back to per-read
    import os
    num_threads = int(os.environ.get("F5C_RNA_THREADS", "0"))
    aligner = f5c_rna.Aligner(pore=pore, use_gpu=False,
                               hmm_confidence=True, num_thread=num_threads)
    use_batch = hasattr(f5c_rna, 'align_reads_variants')

    logger.info("f5c_rna scoring: %d reads, %d total pairs, batch=%s",
                len(reads_variants),
                sum(len(v) for v in reads_variants.values()),
                use_batch)

    if use_batch:
        try:
            all_results = f5c_rna.align_reads_variants(
                signal_path, reads_variants,
                pore=pore, num_thread=num_threads,
            )
        except Exception as e:
            logger.warning("f5c_rna batch failed, falling back to per-read: %s", e)
            use_batch = False

    if not use_batch:
        # Per-read fallback with align_read_variants (singular)
        all_results = {}
        for rid, seq_dict in reads_variants.items():
            try:
                res = f5c_rna.align_read_variants(
                    blow5_path=signal_path, read_id=rid,
                    sequences=seq_dict, pore=pore,
                    use_gpu=False, aligner=aligner,
                )
                all_results[rid] = res
            except Exception:
                continue

    # Parse results and fill dist matrix
    rid_to_i = {rid: i for i, rid in enumerate(read_ids)}
    for rid, res_list in all_results.items():
        i = rid_to_i.get(rid)
        if i is None:
            continue
        cand_map = rid_cand_idx.get(rid, {})

        cand_matched = {}
        cand_events = {}
        for res in res_list:
            label = res.get("variant_label", "?")
            status = res.get("status", -1)
            ev = res.get("event_level_mean", np.array([]))
            n = len(ev)
            if status != 0 or n == 0:
                cand_matched[label] = 0
                cand_events[label] = 0
            else:
                z = res["standardized_level"]
                mk = res.get("model_kmer", [])
                nnnnn = np.array([k == "NNNNN" for k in mk])
                matched = int(np.sum((np.abs(z) < 2.0) & ~nnnnn))
                cand_matched[label] = matched
                cand_events[label] = n

        max_ev = max(cand_events.values()) if cand_events else 1
        if max_ev == 0:
            continue

        for cid, j in cand_map.items():
            matched = cand_matched.get(cid, 0)
            # mismatch_ratio: lower = better match. FAIL → 1.0
            dist[i, j] = 1.0 - (matched / max_ev) if max_ev > 0 else 1.0

    # Center each row: best (lowest mismatch) → 0, like M1's max(AS)-AS
    for i in range(n_reads):
        row = dist[i]
        valid = row[row < 0.999]  # exclude default 1.0 (no data)
        if len(valid) > 0:
            row_min = valid.min()
            dist[i] = dist[i] - row_min

    if as_distance:
        return dist
    else:
        # Convert to support score: 1/(1+dist), then normalize row-sum=1
        scores = 1.0 / (1.0 + dist.astype(np.float64))
        row_sums = scores.sum(axis=1, keepdims=True)
        row_sums = np.where(row_sums > 0, row_sums, 1.0)
        return (scores / row_sums).astype(np.float32)


def f5c_rna_tiebreak(
    R: np.ndarray,
    read_ids: List[str],
    read_seqs: Dict[str, str],
    candidates: List[TranscriptCandidate],
    signal_path: str,
    pore: str = "rna002",
    ambig_threshold: float = 0.90,
) -> np.ndarray:
    """Resolve M1 ties using f5c_rna in-memory eventalign.

    Returns updated R matrix (copy).
    """
    import krill as f5c_rna

    R_new = R.copy()
    n_reads, n_cands = R.shape
    if n_cands < 2:
        return R_new

    # Build mappy aligners
    aligners = {}
    for j, c in enumerate(candidates):
        if c.sequence and len(c.sequence) > 20:
            a = mappy.Aligner(seq=c.sequence, preset="map-ont")
            if a:
                aligners[j] = a

    # Create f5c_rna aligner (reuse for all reads)
    aligner = f5c_rna.Aligner(pore=pore, use_gpu=False, hmm_confidence=True)

    n_ambiguous = 0
    n_tiebroken = 0

    for i in range(n_reads):
        best_R = float(R[i].max())
        if best_R >= ambig_threshold:
            continue

        n_ambiguous += 1
        rid = read_ids[i]
        rseq = read_seqs.get(rid, "")
        if not rseq:
            continue

        # Find tied candidates
        tied_indices = np.where(R[i] > 0.1 * best_R)[0]
        if len(tied_indices) < 2:
            continue

        # Get mappy hits → mapped region sequences
        mapped_seqs = {}
        for j in tied_indices:
            if j not in aligners:
                continue
            best_hit = None
            for hit in aligners[j].map(rseq):
                if best_hit is None or hit.mlen > best_hit.mlen:
                    best_hit = hit
            if best_hit:
                seq = candidates[j].sequence[best_hit.r_st:best_hit.r_en]
                mapped_seqs[candidates[j].candidate_id] = (j, seq)

        if len(mapped_seqs) < 2:
            continue

        # Run f5c_rna for all tied candidates at once
        seq_dict = {cid: seq for cid, (j, seq) in mapped_seqs.items()}
        try:
            results = f5c_rna.align_read_variants(
                blow5_path=signal_path,
                read_id=rid,
                sequences=seq_dict,
                pore=pore,
                use_gpu=False,
                aligner=aligner,
            )
        except Exception as e:
            logger.debug("f5c_rna failed for %s: %s", rid[:20], e)
            continue

        # Compute matched events per candidate
        cand_matched = {}
        cand_events = {}
        for res in results:
            label = res.get("variant_label", "?")
            status = res.get("status", -1)
            ev = res.get("event_level_mean", np.array([]))
            n = len(ev)

            if status != 0 or n == 0:
                cand_matched[label] = 0
                cand_events[label] = 0
            else:
                z = res["standardized_level"]
                model_kmers = res.get("model_kmer", [])
                nnnnn_mask = np.array([k == "NNNNN" for k in model_kmers])
                good_z = np.abs(z) < 2.0
                matched = int(np.sum(good_z & ~nnnnn_mask))
                cand_matched[label] = matched
                cand_events[label] = n

        # Score = matched / max_events
        max_ev = max(cand_events.values()) if cand_events else 1
        if max_ev == 0:
            continue

        scores = {}
        for cid, (j, seq) in mapped_seqs.items():
            scores[j] = cand_matched.get(cid, 0) / max_ev

        # Check if scores differ enough
        score_vals = list(scores.values())
        if max(score_vals) - min(score_vals) < 0.02:
            continue

        # Soft redistribution preserving original row sum
        tied_set = set(tied_indices)
        original_tied_sum = sum(float(R[i, j]) for j in tied_set)
        total_s = sum(scores.values())
        if total_s <= 0:
            continue

        for j in tied_indices:
            if j in scores:
                R_new[i, j] = original_tied_sum * (scores[j] / total_s)
            else:
                R_new[i, j] = 0.0

        n_tiebroken += 1

    logger.info(
        "f5c_rna tiebreak: %d ambiguous, %d tiebroken",
        n_ambiguous, n_tiebroken,
    )
    return R_new
