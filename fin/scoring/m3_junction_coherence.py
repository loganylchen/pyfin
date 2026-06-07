"""M3 junction-window read x read coherence (winner-anchored).

M3 is a read-to-read coherence matrix fed to :func:`em_with_coherence` (the
``dist_read_to_read`` slot, weight ``beta``): reads that look alike over the
discriminative junction window should be pulled onto the same transcript.

Unlike M2 (a per-(read, candidate) signal *distance*), a read x read matrix
needs ONE signal extraction per read -- but the junction window is defined per
candidate. We resolve this by anchoring each read to its OWN winner candidate:

  * M3-1: winner = the read's M1 argmax candidate (``argmin`` of the M1 distance).
  * M3-2: winner = the read's M2 best candidate (``argmin`` of the M2 distance).

The window itself is the *class-level* genomic union
(:func:`class_junction_window_set`), so every read in a wobble class is scored
over the same set of genomic positions; only the projection frame (each read's
winner ``tx2genome``) is per-read. For each read we run krill eventalign once
against its winner, keep the ``event_level_mean`` of events whose genomic
projection lands in the class window (ordered by genomic position), robustly
normalize, and take pairwise krill global DTW *within each class*. Cross-class
and no-window pairs are 0 (no coherence signal).

``event_level_mean`` is the raw per-event current (reference-independent), so
two reads of the same true transcript, both winner-anchored correctly, yield
homologous window-signal sequences -> low DTW; a mis-assigned read extracts
through the wrong frame -> a different event set -> higher DTW.
"""
from __future__ import annotations

import logging
import math
from typing import Dict, List, Optional

import numpy as np

from fin.candidates.dataclasses import TranscriptCandidate
from fin.scoring.diff_region_dtw import cluster_candidates_by_chain
from fin.scoring.m2_junction_nll import class_junction_window_set, _tx2genome
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit

logger = logging.getLogger(__name__)


def _window_signal(
    res: dict,
    cand: TranscriptCandidate,
    gset: set,
) -> Optional[np.ndarray]:
    """Robust-normalized ``event_level_mean`` sequence over in-window events.

    Events are ordered by genomic position (ascending) so two reads traverse
    the same genomic direction. Returns ``None`` when no event lands in
    ``gset``.
    """
    pos = res.get("position", [])
    elm = res.get("event_level_mean", [])
    if len(pos) == 0 or len(elm) == 0:
        return None
    pos = np.asarray(pos, dtype=np.int64)
    elm = np.asarray(elm, dtype=np.float64)
    pts: List[tuple] = []
    for i in range(len(pos)):
        if not math.isfinite(elm[i]):
            continue
        gen = _tx2genome(cand, int(pos[i]))
        if gen is None or gen not in gset:
            continue
        pts.append((gen, elm[i]))
    if not pts:
        return None
    pts.sort(key=lambda t: t[0])
    seg = np.asarray([v for _, v in pts], dtype=np.float32)
    med = float(np.median(seg))
    mad = float(np.median(np.abs(seg - med)))
    scale = mad if mad > 0 else 1.0
    return ((seg - med) / scale).astype(np.float32, copy=False)


def build_m3_coherence(
    read_ids: List[str],
    read_seqs: Dict[str, str],
    candidates: List[TranscriptCandidate],
    winner_col: np.ndarray,
    sig_path: str,
    pore: str = "rna002",
    flank: int = 2,
    chain_wobble: int = 20,
    junction_k: int = 10,
    use_gpu: bool = True,
) -> np.ndarray:
    """Read x read junction-window DTW coherence, anchored to each read's winner.

    Args:
        read_ids: Row/col order.
        read_seqs: read_id -> query sequence.
        candidates: Column order for ``winner_col``.
        winner_col: (n_reads,) int; candidate index each read is anchored to
            (its M1 or M2 argmin). ``< 0`` -> read is left uncoupled (zeros).
        sig_path: blow5/slow5 signal path for krill.
        pore: krill pore model.
        flank: diff-window padding (bp), forwarded to the window builder.
        chain_wobble: junction wobble tolerance (bp) for class clustering.
        junction_k: tx-bp on each side of the junction (window half-width).

    Returns:
        (n_reads, n_reads) float64; symmetric, block-diagonal within each
        wobble class, 0 on the diagonal / cross-class / no-window pairs.
    """
    n = len(read_ids)
    m3 = np.zeros((n, n), dtype=np.float64)
    if n == 0 or not candidates:
        return m3

    winner_col = np.asarray(winner_col, dtype=np.int64)

    # Cluster candidates; only classes with a junction window can couple reads.
    classes = cluster_candidates_by_chain(candidates, chain_wobble=chain_wobble)
    cand_to_class: Dict[int, int] = {}
    for ci, members in enumerate(classes):
        for j in members:
            cand_to_class[j] = ci
    class_gset: Dict[int, set] = {}
    for ci, members in enumerate(classes):
        if len(members) < 2:
            continue
        cands = [candidates[j] for j in members]
        gset = class_junction_window_set(cands, flank=flank, k=junction_k)
        if gset:
            class_gset[ci] = gset
    if not class_gset:
        return m3  # no wobble competition anywhere

    # Signal-path deps imported lazily so the pure no-window paths above stay
    # importable without krill (e.g. unit tests).
    try:
        import krill
        import mappy

        from fin.scoring.krill_aligner import (
            krill_thread_count,
            make_krill_aligner,
        )
    except Exception as exc:  # krill/mappy not importable -> no coherence
        logger.warning("build_m3_coherence: signal deps import failed (%s); skipping", exc)
        return m3

    num_thread = krill_thread_count()
    krill_aligner, use_gpu = make_krill_aligner(
        krill, pore, use_gpu, hmm_confidence=False, num_thread=num_thread
    )
    if krill_aligner is None:
        return m3
    preset = get_m1_preset()
    mappy_by_col: Dict[int, "mappy.Aligner"] = {}

    # Per-read winner-anchored window signal; group reads by class.
    seg_by_read: Dict[int, np.ndarray] = {}
    rows_by_class: Dict[int, List[int]] = {}
    for i, rid in enumerate(read_ids):
        j = int(winner_col[i])
        if j < 0:
            continue
        ci = cand_to_class.get(j, -1)
        gset = class_gset.get(ci)
        if gset is None:
            continue  # winner's class has no discrimination window
        seq = read_seqs.get(rid)
        cand = candidates[j]
        if not seq or not cand.sequence:
            continue
        aln = mappy_by_col.get(j)
        if aln is None:
            aln = mappy.Aligner(seq=cand.sequence, preset=preset)
            mappy_by_col[j] = aln

        best_hit = None
        best_as = None
        for h in aln.map(seq):
            v = score_hit(h)
            if v is None:
                continue
            if best_as is None or v > best_as:
                best_as, best_hit = v, h
        if best_hit is None:
            continue
        try:
            recs = krill.align_read_variants(
                sig_path, rid,
                {cand.candidate_id: cand.sequence[best_hit.r_st:best_hit.r_en]},
                pore=pore, use_gpu=use_gpu, num_thread=num_thread,
                aligner=krill_aligner, start=best_hit.r_st,
            )
        except Exception as e:  # noqa: BLE001 - krill raises broad errors
            logger.debug("align_read_variants failed %s/%s: %s", rid, cand.candidate_id, e)
            continue
        res = {x.get("variant_label"): x for x in recs}.get(cand.candidate_id)
        if res is None or res.get("status", -1) != 0:
            continue
        seg = _window_signal(res, cand, gset)
        if seg is None or seg.size == 0:
            continue
        seg_by_read[i] = seg
        rows_by_class.setdefault(ci, []).append(i)

    # Within-class pairwise krill global DTW; cross-class stays 0.
    n_pairs = 0
    for ci, rows in rows_by_class.items():
        if len(rows) < 2:
            continue
        segs = [np.ascontiguousarray(seg_by_read[i], dtype=np.float32) for i in rows]
        # krill global DTW on the effective device (use_gpu is post-fallback:
        # True only if GPU eventalign succeeded above, so GPU DTW is safe).
        dmat = krill.dtw_pairwise_varlen(segs, use_gpu=use_gpu, num_threads=num_thread)
        for a in range(len(rows)):
            ia = rows[a]
            for b in range(a + 1, len(rows)):
                ib = rows[b]
                d = float(dmat[a, b])
                if math.isfinite(d):
                    m3[ia, ib] = d
                    m3[ib, ia] = d
                    n_pairs += 1

    logger.info(
        "build_m3_coherence: %d windowed classes, %d coupled reads, %d read-pairs",
        len(class_gset), len(seg_by_read), n_pairs,
    )
    return m3
