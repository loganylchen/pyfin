"""Diff-region DTW for read-to-read coherence (m4 source for R4/R5 ablation).

This module computes a read-to-read distance matrix by restricting DTW to
genomic regions where candidates disagree (one calls exon, another calls
intron). Shared exons contribute no information to coherence, so excluding
them denoises the m4 matrix used in `em_with_coherence`.

Public functions:
    extract_diff_regions(candidates) -> List[(g_start, g_end)] in genomic coords
    genomic_region_to_cdna(candidate, region_genome) -> Optional[(c_start, c_end)]
    cdna_region_to_signal_range(score, region_cdna) -> Optional[(s_start, s_end)]
    compute_diff_region_m4(...) -> np.ndarray (N x N)
"""

from __future__ import annotations

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np

from fin.candidates.dataclasses import TranscriptCandidate
from fin.scoring.eventalign_parser import ReadCandidateScore

logger = logging.getLogger(__name__)


def extract_diff_regions(
    candidates: List[TranscriptCandidate],
    k_max: int = 15,
    flank_w: int = 1,
) -> List[Tuple[int, int]]:
    """Return junction-anchored diff regions with shared flanking context.

    A diff region is a maximal contiguous run of bases where at least one
    candidate calls exon AND at least one calls intron, subject to:

    - **Width filter** (``k_max``): region length must be ≤ ``k_max`` bp.
      This focuses on nanopore-noise-scale junction wobble and excludes
      large structural differences (cassette exons, alt 5'/3' UTRs).
    - **Shared-flank filter** (``flank_w``): for the ``flank_w`` bases
      immediately on each side of the region, every candidate that covers
      those flanking bases must (a) be covered at every flank base and
      (b) agree on the same exon/intron label. At least two such "spanning"
      candidates must exist, and the disagreement inside the region must
      come from those spanning candidates. This guarantees the candidates
      being compared share identical context on both sides of the wobble.

    Per-base labels:
        0 = not covered by the candidate
        1 = exon
        2 = intron

    Args:
        candidates: List of TranscriptCandidate (must have intron_chain).
        k_max: Maximum allowed width of a diff region (bp). Default 15.
        flank_w: Width of the shared-context flank required on each side
            (bp). Default 1.

    Returns:
        Sorted list of (g_start, g_end) tuples (half-open) in genomic
        coordinates. Empty if fewer than 2 candidates or no diff regions
        survive filtering.
    """
    if len(candidates) < 2:
        return []

    g_lo = min(c.start for c in candidates)
    g_hi = max(c.end for c in candidates)
    if g_hi <= g_lo:
        return []

    span = g_hi - g_lo
    n = len(candidates)

    # Per-candidate per-base label array.
    labels = np.zeros((n, span), dtype=np.int8)
    for i, cand in enumerate(candidates):
        c_lo = max(cand.start, g_lo) - g_lo
        c_hi = min(cand.end, g_hi) - g_lo
        if c_hi <= c_lo:
            continue
        labels[i, c_lo:c_hi] = 1  # exon
        for i_start, i_end in cand.intron_chain.introns:
            i_lo = max(i_start, g_lo) - g_lo
            i_hi = min(i_end, g_hi) - g_lo
            if i_hi <= i_lo:
                continue
            labels[i, i_lo:i_hi] = 2  # intron

    exon_any = (labels == 1).any(axis=0)
    intron_any = (labels == 2).any(axis=0)
    diff_mask = exon_any & intron_any
    if not diff_mask.any():
        return []

    # Raw diff runs (maximal contiguous diff bases).
    raw: List[Tuple[int, int]] = []
    in_region = False
    seg_start = 0
    for i in range(span):
        if diff_mask[i] and not in_region:
            seg_start = i
            in_region = True
        elif not diff_mask[i] and in_region:
            raw.append((seg_start, i))
            in_region = False
    if in_region:
        raw.append((seg_start, span))

    # Apply k_max width and shared-flank filters.
    filtered: List[Tuple[int, int]] = []
    for lo, hi in raw:
        if hi - lo > k_max:
            continue
        l_lo, l_hi = lo - flank_w, lo
        r_lo, r_hi = hi, hi + flank_w
        if l_lo < 0 or r_hi > span:
            continue

        spanning_indices: List[int] = []
        l_label_ref: Optional[int] = None
        r_label_ref: Optional[int] = None
        ok = True
        for ci in range(n):
            l_seg = labels[ci, l_lo:l_hi]
            r_seg = labels[ci, r_lo:r_hi]
            # Skip candidates that don't fully cover both flanks.
            if (l_seg == 0).any() or (r_seg == 0).any():
                continue
            # Within this candidate, each flank must be a single label
            # (no exon/intron transition inside the flank itself).
            l_set = set(l_seg.tolist())
            r_set = set(r_seg.tolist())
            if len(l_set) != 1 or len(r_set) != 1:
                ok = False
                break
            l_v = int(l_seg[0])
            r_v = int(r_seg[0])
            if l_label_ref is None:
                l_label_ref, r_label_ref = l_v, r_v
            elif l_v != l_label_ref or r_v != r_label_ref:
                ok = False
                break
            spanning_indices.append(ci)
        if not ok or len(spanning_indices) < 2:
            continue

        # The exon/intron disagreement inside the region must come from
        # the spanning (shared-flank) candidates.
        inner_exon = False
        inner_intron = False
        for ci in spanning_indices:
            inner = labels[ci, lo:hi]
            if not inner_exon and (inner == 1).any():
                inner_exon = True
            if not inner_intron and (inner == 2).any():
                inner_intron = True
            if inner_exon and inner_intron:
                break
        if not (inner_exon and inner_intron):
            continue

        filtered.append((lo + g_lo, hi + g_lo))

    return filtered


def genomic_region_to_cdna(
    candidate: TranscriptCandidate,
    region_genome: Tuple[int, int],
) -> Optional[Tuple[int, int]]:
    """Map a genomic region to cDNA-relative coordinates for `candidate`.

    Walks ``candidate.intron_chain``, clips the region to the candidate's
    exonic portion (intronic portion is excluded), and returns cDNA-relative
    half-open coordinates. Accounts for strand: on '-' strand, cDNA positions
    are flipped relative to the spliced transcript length.

    Args:
        candidate: TranscriptCandidate with intron_chain and start/end.
        region_genome: (g_start, g_end) half-open in genomic coordinates.

    Returns:
        (c_start, c_end) half-open cDNA coordinates, or None if no exonic
        overlap with this candidate.
    """
    g_start, g_end = region_genome
    if g_end <= g_start:
        return None

    introns = candidate.intron_chain.introns
    # Build exon list from candidate.start/end and intron_chain.
    if not introns:
        exons = [(candidate.start, candidate.end)]
    else:
        exons = [(candidate.start, introns[0][0])]
        for k in range(len(introns) - 1):
            exons.append((introns[k][1], introns[k + 1][0]))
        exons.append((introns[-1][1], candidate.end))

    # Walk exons left-to-right (genomic order), accumulate cDNA offset.
    cdna_offset = 0
    matched: List[Tuple[int, int]] = []  # (cdna_start, cdna_end) in left-to-right cDNA frame
    for ex_s, ex_e in exons:
        ex_len = ex_e - ex_s
        if ex_len <= 0:
            continue
        ov_s = max(ex_s, g_start)
        ov_e = min(ex_e, g_end)
        if ov_e > ov_s:
            c_s = cdna_offset + (ov_s - ex_s)
            c_e = cdna_offset + (ov_e - ex_s)
            matched.append((c_s, c_e))
        cdna_offset += ex_len

    if not matched:
        return None

    spliced_len = cdna_offset  # total exonic bases = full cDNA length
    c_start = min(s for s, _ in matched)
    c_end = max(e for _, e in matched)

    # On the minus strand, cDNA is reverse-complemented: position i in the
    # genomic-order frame corresponds to (spliced_len - 1 - i) in cDNA order.
    if candidate.strand == "-":
        c_start_rc = spliced_len - c_end
        c_end_rc = spliced_len - c_start
        c_start, c_end = c_start_rc, c_end_rc

    if c_end <= c_start:
        return None
    return (c_start, c_end)


def cdna_region_to_signal_range(
    score: ReadCandidateScore,
    region_cdna: Tuple[int, int],
) -> Optional[Tuple[int, int]]:
    """Map a cDNA region to a (signal_start_idx, signal_end_idx) range for a read.

    Iterates the opt-in ``score.events`` list (per-event records) and returns
    the min/max signal index over events whose reference_position falls within
    ``region_cdna`` (half-open).

    Args:
        score: ReadCandidateScore with events populated (collect_events=True).
        region_cdna: (c_start, c_end) half-open cDNA range.

    Returns:
        (sig_lo, sig_hi) or None if no events fall in the region.
    """
    if not score.events:
        return None
    c_start, c_end = region_cdna
    sig_lo: Optional[int] = None
    sig_hi: Optional[int] = None
    for pos, s_lo, s_hi in score.events:
        if c_start <= pos < c_end:
            if sig_lo is None or s_lo < sig_lo:
                sig_lo = s_lo
            if sig_hi is None or s_hi > sig_hi:
                sig_hi = s_hi
    if sig_lo is None or sig_hi is None or sig_hi <= sig_lo:
        return None
    return (sig_lo, sig_hi)


def _dtw_distance(seq1: np.ndarray, seq2: np.ndarray, use_gpu: bool) -> float:
    """Single-pair DTW distance with GPU fast path and CPU fallback."""
    if len(seq1) == 0 or len(seq2) == 0:
        return float("nan")
    if use_gpu:
        try:
            from fin._dtw import dtw_pairwise_varlen, is_available

            if is_available() and len(seq1) >= 16 and len(seq2) >= 16:
                pw = dtw_pairwise_varlen([np.asarray(seq1, dtype=np.float32),
                                          np.asarray(seq2, dtype=np.float32)],
                                         use_open_end=True)
                return float(pw[0, 1])
        except (ImportError, RuntimeError) as e:
            logger.debug("GPU DTW unavailable, using CPU: %s", e)

    from fin.scoring.signal_dtw import _cpu_dtw

    return float(_cpu_dtw(
        np.asarray(seq1, dtype=np.float32),
        np.asarray(seq2, dtype=np.float32),
    ))


def _best_candidate_per_read(
    read_id: str,
    candidates: List[TranscriptCandidate],
    scores: Dict[Tuple[str, str], ReadCandidateScore],
) -> Optional[Tuple[TranscriptCandidate, ReadCandidateScore]]:
    """Return the (candidate, score) with the highest total_log_likelihood for read."""
    best: Optional[Tuple[TranscriptCandidate, ReadCandidateScore]] = None
    for cand in candidates:
        s = scores.get((read_id, cand.candidate_id))
        if s is None:
            continue
        if best is None or s.total_log_likelihood > best[1].total_log_likelihood:
            best = (cand, s)
    return best


def _read_signal(signal_reader, read_id: str, signal_format: str) -> Optional[np.ndarray]:
    """Fetch the raw signal for a read; returns None on miss."""
    if signal_reader is None:
        return None
    try:
        if signal_format == "pod5":
            result = signal_reader.get_calibrated_signal(read_id)
        else:
            result = signal_reader.get_picoamp_signal(read_id)
    except Exception as e:
        logger.debug("Signal fetch failed for %s: %s", read_id, e)
        return None
    if result is None:
        return None
    sig, _meta = result
    return np.asarray(sig, dtype=np.float32)


def compute_diff_region_m4(
    read_ids: List[str],
    candidates: List[TranscriptCandidate],
    scores_by_pair: Dict[Tuple[str, str], ReadCandidateScore],
    signal_reader,
    interval_start: int,
    interval_end: int,
    signal_format: str = "slow5",
    use_gpu: bool = True,
    normalize: bool = True,
) -> np.ndarray:
    """Compute read-to-read distance matrix using diff-region DTW.

    For each read pair (i, j):
      1. Find each read's best-scoring candidate (host transcript).
      2. For each genomic diff region, map to cDNA (per each read's host)
         and then to signal range via events.
      3. Run DTW(sig_i_region, sig_j_region); aggregate across regions by
         simple mean. Skip a region if either read has no eventalign coverage.

    Cells with no covered region remain NaN; the caller's NaN adapter
    (runner) is responsible for sanitizing before EM (AC8).

    Args:
        read_ids: Ordered list of read IDs (N rows/cols).
        candidates: List of TranscriptCandidate spanning this interval.
        scores_by_pair: Mapping (read_id, candidate_id) -> ReadCandidateScore
            with ``events`` populated (collect_events=True upstream).
        signal_reader: Open Pod5Reader or Slow5Reader.
        interval_start, interval_end: Genomic span of the interval (used only
            for sanity logging; diff regions are computed from candidates).
        signal_format: "pod5" or "slow5"/"blow5".
        use_gpu: Whether to attempt GPU DTW.
        normalize: Per-read robust z-score (median / MAD) before DTW.

    Returns:
        np.ndarray shape (N, N) dtype float32 with pairwise distances; NaN
        where no diff region was scored.
    """
    n = len(read_ids)
    m4 = np.full((n, n), np.nan, dtype=np.float32)
    if n == 0:
        return m4

    diff_regions = extract_diff_regions(candidates)
    if not diff_regions:
        logger.info(
            "compute_diff_region_m4: no diff regions in interval %d-%d",
            interval_start, interval_end,
        )
        # No diff regions -> return zero matrix (no coherence signal). The
        # runner's NaN adapter will treat all-NaN rows as zero anyway, but
        # returning zeros here is a cheaper communication.
        return np.zeros((n, n), dtype=np.float32)

    # Precompute per-read: (host candidate, host score) and the cached
    # full-read signal (avoid re-fetching for every pair).
    host_by_read: Dict[str, Optional[Tuple[TranscriptCandidate, ReadCandidateScore]]] = {}
    signal_by_read: Dict[str, Optional[np.ndarray]] = {}
    for rid in read_ids:
        host = _best_candidate_per_read(rid, candidates, scores_by_pair)
        host_by_read[rid] = host
        signal_by_read[rid] = _read_signal(signal_reader, rid, signal_format)

    # For each read, precompute its per-region signal segments.
    # segments_by_read[rid] is a dict: region_idx -> np.ndarray segment.
    segments_by_read: Dict[str, Dict[int, np.ndarray]] = {rid: {} for rid in read_ids}
    for rid in read_ids:
        host = host_by_read[rid]
        sig = signal_by_read[rid]
        if host is None or sig is None:
            continue
        cand, score = host
        for r_idx, region_g in enumerate(diff_regions):
            region_c = genomic_region_to_cdna(cand, region_g)
            if region_c is None:
                continue
            sig_range = cdna_region_to_signal_range(score, region_c)
            if sig_range is None:
                continue
            s_lo, s_hi = sig_range
            s_lo = max(0, s_lo)
            s_hi = min(len(sig), s_hi)
            if s_hi <= s_lo:
                continue
            seg = sig[s_lo:s_hi]
            if normalize:
                med = float(np.median(seg))
                mad = float(np.median(np.abs(seg - med)))
                scale = mad if mad > 0 else 1.0
                seg = ((seg - med) / scale).astype(np.float32, copy=False)
            segments_by_read[rid][r_idx] = seg

    # Pairwise: per-region DTW averaged.
    for i in range(n):
        m4[i, i] = 0.0
        segs_i = segments_by_read[read_ids[i]]
        if not segs_i:
            continue
        for j in range(i + 1, n):
            segs_j = segments_by_read[read_ids[j]]
            if not segs_j:
                continue
            common = set(segs_i.keys()) & set(segs_j.keys())
            if not common:
                continue
            dists: List[float] = []
            for r_idx in common:
                d = _dtw_distance(segs_i[r_idx], segs_j[r_idx], use_gpu=use_gpu)
                if np.isfinite(d):
                    dists.append(d)
            if dists:
                mean_d = float(np.mean(dists))
                m4[i, j] = mean_d
                m4[j, i] = mean_d

    return m4
