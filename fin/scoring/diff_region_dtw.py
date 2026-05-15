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
    cluster_candidates_by_chain(...) -> List[List[int]]  (class -> cand indices)
    compute_class_partitioned_m4(...) -> np.ndarray (N x N)
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
    pad_bp: int = 6,
) -> List[Tuple[int, int]]:
    """Return padded diff regions for SAME-CLASS candidates.

    Same-class candidates differ only by per-junction wobble (≤ chain_wobble
    bp per donor/acceptor; see :func:`cluster_candidates_by_chain`). This
    function locates each junction-wobble disagreement and returns a genomic
    window that brackets the full affected intron with ``pad_bp`` flanking
    exonic context on both sides.

    Algorithm:
        1. Build per-base exon/intron labels per candidate.
        2. Find ``diff_mask = exon_any & intron_any`` and segment into raw
           contiguous clusters.
        3. For each raw cluster, take the **union** of every overlapping
           intron interval from any candidate. This grows the window from
           "just the disagreeing bases" to "the entire affected intron".
        4. Pad the union by ``pad_bp`` on each side (clipped to the
           candidate span).
        5. Merge overlapping/adjacent padded windows.

    **The caller is responsible for passing candidates from the same class.**
    Cross-class structural differences (cassette exons, alt UTRs) are out of
    scope: they are resolved by the chain partition itself, not by signal
    DTW.

    Per-base labels (internal):
        0 = not covered by the candidate
        1 = exon
        2 = intron

    Args:
        candidates: List of same-class TranscriptCandidate.
        pad_bp: Flanking exonic context (bp) added to each side of the
            unioned intron window. Default 6.

    Returns:
        Sorted list of (g_start, g_end) tuples (half-open) in genomic
        coordinates. Empty if fewer than 2 candidates or no disagreement.
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

    # Raw diff runs (maximal contiguous bases of exon/intron disagreement).
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

    # For each raw cluster, union all overlapping intron intervals across
    # candidates, then pad ±pad_bp.
    windows: List[Tuple[int, int]] = []
    for lo, hi in raw:
        g_cluster_lo = lo + g_lo
        g_cluster_hi = hi + g_lo
        u_lo = g_cluster_lo
        u_hi = g_cluster_hi
        for cand in candidates:
            for i_s, i_e in cand.intron_chain.introns:
                # Half-open overlap with the raw cluster.
                if i_e > g_cluster_lo and i_s < g_cluster_hi:
                    if i_s < u_lo:
                        u_lo = i_s
                    if i_e > u_hi:
                        u_hi = i_e
        pad_lo = max(g_lo, u_lo - pad_bp)
        pad_hi = min(g_hi, u_hi + pad_bp)
        if pad_hi > pad_lo:
            windows.append((pad_lo, pad_hi))

    if not windows:
        return []

    # Merge overlapping/adjacent windows (sorted by start).
    windows.sort()
    merged: List[Tuple[int, int]] = [windows[0]]
    for s, e in windows[1:]:
        prev_s, prev_e = merged[-1]
        if s <= prev_e:
            merged[-1] = (prev_s, max(prev_e, e))
        else:
            merged.append((s, e))
    return merged


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


def _any_cand_signal_range(
    read_id: str,
    region_g: Tuple[int, int],
    candidates: List[TranscriptCandidate],
    scores: Dict[Tuple[str, str], ReadCandidateScore],
    sig_len: int,
) -> Optional[Tuple[int, int]]:
    """Find a signal range for (read, diff_region) using ANY candidate.

    The single-host projection used by the original implementation fails when
    the read's max-LL host places the diff region inside an intron — then
    ``genomic_region_to_cdna`` returns None and the region is skipped, which
    on some intervals cascades to off-diag NaN frac = 100%.

    This helper iterates every candidate that has an eventalign score for
    ``read_id``, projects the genomic region through that candidate's exon
    structure, and picks the candidate whose cDNA→signal projection captures
    the most eventalign events in the region. Returns the resulting signal
    range clipped to ``[0, sig_len)``, or None if no candidate qualifies.

    Args:
        read_id: The read whose signal we want to slice.
        region_g: Genomic diff region (g_start, g_end), half-open.
        candidates: All candidates in this interval (any may carry the region
            in its exonic body).
        scores: Mapping (read_id, candidate_id) -> ReadCandidateScore with
            ``events`` populated.
        sig_len: Length of the read's raw signal (for clipping).

    Returns:
        (s_lo, s_hi) signal-sample range or None.
    """
    best_n = 0
    best_sig: Optional[Tuple[int, int]] = None
    for cand in candidates:
        score = scores.get((read_id, cand.candidate_id))
        if score is None or not score.events:
            continue
        region_c = genomic_region_to_cdna(cand, region_g)
        if region_c is None:
            continue
        c_lo, c_hi = region_c
        sig_lo: Optional[int] = None
        sig_hi: Optional[int] = None
        cnt = 0
        for pos, s_lo, s_hi in score.events:
            if c_lo <= pos < c_hi:
                cnt += 1
                if sig_lo is None or s_lo < sig_lo:
                    sig_lo = s_lo
                if sig_hi is None or s_hi > sig_hi:
                    sig_hi = s_hi
        if cnt > best_n and sig_lo is not None and sig_hi is not None:
            s_lo_c = max(0, int(sig_lo))
            s_hi_c = min(int(sig_len), int(sig_hi))
            if s_hi_c > s_lo_c:
                best_n = cnt
                best_sig = (s_lo_c, s_hi_c)
    return best_sig


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
    chain_wobble: int = 6,
    pad_bp: int = 6,
) -> np.ndarray:
    """Compute read-to-read distance matrix using class-partitioned diff-region DTW.

    Pipeline:
      1. Cluster candidates by intron-chain wobble (``chain_wobble`` bp per
         junction). This produces transcript-form classes; cross-class
         differences are structural (cassette exons, alt UTRs) and resolved
         by the partition itself.
      2. Assign each read to a class via its max-LL host candidate
         (highest ``total_log_likelihood`` over candidates in
         ``scores_by_pair``).
      3. For each class:
         a. Extract padded diff regions across the class's candidates
            (per-junction wobble + ``pad_bp`` flanking exon).
         b. For each read in the class and each diff region, find ANY
            candidate (within the class) whose cDNA→signal projection
            captures the most eventalign events. Skip the region for that
            read if no candidate qualifies.
         c. Pairwise DTW between reads in the same class, averaged across
            common diff regions.

    Output structure:
      - ``M[i, i] = 0`` (diagonal).
      - ``M[i, j] = 0`` when reads i, j are in different classes (no
        information; class separation is handled by chain partition).
      - ``M[i, j] = 0`` when reads are in a singleton class (no intra-class
        disagreement → coherent by construction).
      - ``M[i, j] = mean DTW distance`` over common diff regions otherwise.
      - ``M[i, j] = 0`` if no diff region is commonly covered for the pair.

    Args:
        read_ids: Ordered list of read IDs (N rows/cols).
        candidates: All candidates in this interval.
        scores_by_pair: Mapping (read_id, candidate_id) -> ReadCandidateScore
            with ``events`` populated (collect_events=True upstream). Also
            provides ``total_log_likelihood`` for class assignment.
        signal_reader: Open Pod5Reader or Slow5Reader.
        interval_start, interval_end: Genomic span (used for logging only).
        signal_format: "pod5" or "slow5"/"blow5".
        use_gpu: Whether to attempt GPU DTW.
        normalize: Per-segment robust z-score (median / MAD) before DTW.
        chain_wobble: ±bp tolerance per junction for class clustering.
        pad_bp: Flanking exon pad (bp) added on each side of unioned intron.

    Returns:
        np.ndarray shape (N, N) dtype float32, block-diagonal with zeros
        in cross-class entries.
    """
    n = len(read_ids)
    m3 = np.zeros((n, n), dtype=np.float32)
    if n == 0 or not candidates:
        return m3

    # 1. Cluster candidates by chain.
    classes = cluster_candidates_by_chain(candidates, chain_wobble=chain_wobble)
    if not classes:
        return m3

    cand_idx_to_class: Dict[int, int] = {}
    for ci, members in enumerate(classes):
        for j in members:
            cand_idx_to_class[j] = ci

    # 2. Assign each read to a class via its max-LL host candidate.
    read_to_class: Dict[str, int] = {}
    for rid in read_ids:
        best_ll = -np.inf
        best_cand_idx: Optional[int] = None
        for cand_idx, cand in enumerate(candidates):
            s = scores_by_pair.get((rid, cand.candidate_id))
            if s is None:
                continue
            if s.total_log_likelihood > best_ll:
                best_ll = s.total_log_likelihood
                best_cand_idx = cand_idx
        if best_cand_idx is not None:
            read_to_class[rid] = cand_idx_to_class.get(best_cand_idx, -1)
        else:
            read_to_class[rid] = -1

    # Group read row indices by class.
    reads_by_class: Dict[int, List[int]] = {}
    for i, rid in enumerate(read_ids):
        ci = read_to_class[rid]
        if ci >= 0:
            reads_by_class.setdefault(ci, []).append(i)

    # Cache full-read signal once per read (any class).
    signal_by_read: Dict[str, Optional[np.ndarray]] = {}
    for rid in read_ids:
        signal_by_read[rid] = _read_signal(signal_reader, rid, signal_format)

    # 3. Process each class independently. Singleton classes and classes with
    # no intra-class disagreement contribute zeros (already initialized).
    any_region = False
    for ci, class_cand_idxs in enumerate(classes):
        row_idxs = reads_by_class.get(ci, [])
        if len(row_idxs) < 2:
            continue  # singleton or empty: zeros suffice
        class_cands = [candidates[j] for j in class_cand_idxs]
        diff_regions = extract_diff_regions(class_cands, pad_bp=pad_bp)
        if not diff_regions:
            continue  # class is fully coherent; zeros suffice
        any_region = True

        # Precompute per-read segments for this class.
        segments_by_row: Dict[int, Dict[int, np.ndarray]] = {}
        for i in row_idxs:
            rid = read_ids[i]
            sig = signal_by_read[rid]
            if sig is None:
                continue
            segs: Dict[int, np.ndarray] = {}
            for r_idx, region_g in enumerate(diff_regions):
                seg_range = _any_cand_signal_range(
                    rid, region_g, class_cands, scores_by_pair, len(sig),
                )
                if seg_range is None:
                    continue
                s_lo, s_hi = seg_range
                seg = sig[s_lo:s_hi]
                if seg.size == 0:
                    continue
                if normalize:
                    med = float(np.median(seg))
                    mad = float(np.median(np.abs(seg - med)))
                    scale = mad if mad > 0 else 1.0
                    seg = ((seg - med) / scale).astype(np.float32, copy=False)
                segs[r_idx] = seg
            segments_by_row[i] = segs

        # Pairwise DTW within this class.
        for a, i in enumerate(row_idxs):
            segs_i = segments_by_row.get(i, {})
            if not segs_i:
                continue
            for j in row_idxs[a + 1:]:
                segs_j = segments_by_row.get(j, {})
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
                    m3[i, j] = mean_d
                    m3[j, i] = mean_d

    if not any_region:
        logger.info(
            "compute_diff_region_m4: no diff regions across %d classes in interval %d-%d",
            len(classes), interval_start, interval_end,
        )
    return m3


# ---------------------------------------------------------------------------
# Class-partitioned m4 (intron-chain clustering + within-class full-read DTW)
# ---------------------------------------------------------------------------


def _chains_match_wobble(
    a: Tuple[Tuple[int, int], ...],
    b: Tuple[Tuple[int, int], ...],
    wobble: int,
) -> bool:
    """Wobble-tolerant equality between two intron chains.

    Same number of introns, and every donor/acceptor pair differs by at most
    ``wobble`` bp.
    """
    if len(a) != len(b):
        return False
    for (sa, ea), (sb, eb) in zip(a, b):
        if abs(sa - sb) > wobble or abs(ea - eb) > wobble:
            return False
    return True


def cluster_candidates_by_chain(
    candidates: List[TranscriptCandidate],
    chain_wobble: int = 6,
) -> List[List[int]]:
    """Group candidates into classes of wobble-tolerant intron-chain equality.

    Greedy single-linkage by representative: each new candidate joins the first
    existing class whose representative chain matches within ``chain_wobble``
    bp per donor and per acceptor (and same intron count). Otherwise starts a
    new class. Single-exon candidates form their own class iff every other
    candidate is multi-exon.

    Args:
        candidates: List of TranscriptCandidate.
        chain_wobble: ± bp tolerance per junction. Default 6.

    Returns:
        List of classes; each class is a list of candidate indices (into
        ``candidates``). Indices in a class are sorted ascending.
    """
    classes: List[List[int]] = []
    reps: List[Tuple[Tuple[int, int], ...]] = []
    for i, cand in enumerate(candidates):
        chain = cand.intron_chain.introns
        placed = False
        for ci, rep in enumerate(reps):
            if _chains_match_wobble(chain, rep, chain_wobble):
                classes[ci].append(i)
                placed = True
                break
        if not placed:
            reps.append(chain)
            classes.append([i])
    return classes


def _assign_reads_to_classes(
    read_ids: List[str],
    cand_ids: List[str],
    classes: List[List[int]],
    dist_read_cand: np.ndarray,
) -> List[int]:
    """Map each read to a class index via argmin distance to any candidate.

    A read is assigned to the class containing its single best-scoring
    (lowest-distance) candidate. Reads whose row is entirely NaN/inf or
    sentinel-large get class = -1 (unassigned).

    Args:
        read_ids: Ordered read IDs (rows of dist_read_cand).
        cand_ids: Ordered candidate IDs (columns of dist_read_cand).
        classes: Output of ``cluster_candidates_by_chain`` over the SAME
            ordered candidate list as ``cand_ids``.
        dist_read_cand: (n_reads, n_cands) distance matrix (lower=better).
            Missing entries should be NaN or ≥1e5.

    Returns:
        List of class indices, length ``len(read_ids)``. ``-1`` for reads
        with no usable distance.
    """
    if dist_read_cand.size == 0 or len(classes) == 0:
        return [-1] * len(read_ids)
    n_cands = len(cand_ids)
    cand_to_class = [-1] * n_cands
    for ci, members in enumerate(classes):
        for j in members:
            if 0 <= j < n_cands:
                cand_to_class[j] = ci

    assigned: List[int] = []
    M = np.asarray(dist_read_cand, dtype=np.float64)
    # Treat sentinel-large and non-finite as "missing".
    M = np.where(np.isfinite(M) & (M < 1e5), M, np.inf)
    for i in range(len(read_ids)):
        row = M[i]
        if not np.any(np.isfinite(row)):
            assigned.append(-1)
            continue
        j = int(np.argmin(row))
        ci = cand_to_class[j] if 0 <= j < n_cands else -1
        assigned.append(ci)
    return assigned


def compute_class_partitioned_m4(
    read_ids: List[str],
    candidates: List[TranscriptCandidate],
    cand_ids: List[str],
    dist_read_cand: np.ndarray,
    signal_reader,
    scores_by_pair: Optional[Dict[Tuple[str, str], ReadCandidateScore]] = None,
    signal_format: str = "slow5",
    use_gpu: bool = True,
    normalize: bool = True,
    chain_wobble: int = 6,
    inter_class_fill: object = "max",
) -> np.ndarray:
    """Class-partitioned read-to-read DTW (alternative m4 source).

    Pipeline:
      1. Cluster ``candidates`` by intron-chain wobble (``chain_wobble`` bp).
      2. Assign each read to a class via ``argmin(dist_read_cand)`` row-wise
         (M1 mappy or M2 eventalign — caller chooses).
      3. Extract per-read signal (eventalign signal-range if a host
         ReadCandidateScore is available, else full read) and z-normalize.
      4. Compute pairwise DTW only between reads in the same class.
      5. Inter-class entries are filled with either ``max(intra_class_dists)``
         (when ``inter_class_fill == "max"``) or the float constant provided.
      6. Self-distance = 0 on the diagonal.

    Args:
        read_ids: Ordered list of read IDs (N rows/cols).
        candidates: Candidate list (defines cluster membership).
        cand_ids: Column order of ``dist_read_cand``; must align by index to
            the candidate list (or a subset; indices that fall outside
            ``candidates`` are ignored).
        dist_read_cand: (N_reads, N_cands) distance matrix; lower = better
            match. Used solely for class assignment.
        signal_reader: Pod5/Slow5 reader.
        scores_by_pair: Optional eventalign scores; when present and the
            assigned host pair exists, signal is sliced to the eventalign
            range for that pair (mirrors the dense whole-read m4 path).
        signal_format: ``"pod5"`` or ``"slow5"``/``"blow5"``.
        use_gpu: Attempt GPU DTW.
        normalize: Robust z-score per segment.
        chain_wobble: Junction wobble tolerance (bp) for clustering.
        inter_class_fill: ``"max"`` (default; use max of computed intra-class
            distances; falls back to 0 when no intra-class pair was computed)
            or a numeric constant.

    Returns:
        np.ndarray shape (N, N) dtype float32.
    """
    n = len(read_ids)
    m4 = np.zeros((n, n), dtype=np.float32)
    if n == 0:
        return m4
    if not candidates:
        return m4

    classes = cluster_candidates_by_chain(candidates, chain_wobble=chain_wobble)
    if not classes:
        return m4

    # Align cand_ids to the candidates list order. Caller is expected to pass
    # cand_ids = [c.candidate_id for c in candidates], but we don't enforce it.
    id_to_idx = {c.candidate_id: i for i, c in enumerate(candidates)}
    cand_order_idx = [id_to_idx.get(cid, -1) for cid in cand_ids]

    # Build "classes-over-cand_ids-columns" so _assign_reads_to_classes works
    # on the dist matrix column order.
    classes_in_col_order: List[List[int]] = [[] for _ in classes]
    cand_idx_to_class: Dict[int, int] = {}
    for ci, members in enumerate(classes):
        for j in members:
            cand_idx_to_class[j] = ci
    for col, cand_i in enumerate(cand_order_idx):
        ci = cand_idx_to_class.get(cand_i)
        if ci is not None and ci >= 0:
            classes_in_col_order[ci].append(col)

    read_class = _assign_reads_to_classes(
        read_ids, cand_ids, classes_in_col_order, dist_read_cand,
    )

    # Fetch per-read signal segments (best-host-based slice if scores known).
    segments: Dict[str, np.ndarray] = {}
    for i, rid in enumerate(read_ids):
        ci = read_class[i]
        if ci < 0:
            continue
        sig = _read_signal(signal_reader, rid, signal_format)
        if sig is None or len(sig) == 0:
            continue
        s_lo, s_hi = 0, len(sig)
        # If we have a ReadCandidateScore for the best (read, candidate) in
        # this class, use its eventalign signal range to trim.
        if scores_by_pair:
            host_score: Optional[ReadCandidateScore] = None
            best_ll = -float("inf")
            for j in classes[ci]:
                cid = candidates[j].candidate_id
                s = scores_by_pair.get((rid, cid))
                if s is None:
                    continue
                if s.total_log_likelihood > best_ll:
                    best_ll = s.total_log_likelihood
                    host_score = s
            if host_score is not None:
                lo = max(0, int(host_score.signal_start_idx))
                hi = min(len(sig), int(host_score.signal_end_idx))
                if hi > lo:
                    s_lo, s_hi = lo, hi
        seg = sig[s_lo:s_hi]
        if seg.size == 0:
            continue
        if normalize:
            med = float(np.median(seg))
            mad = float(np.median(np.abs(seg - med)))
            scale = mad if mad > 0 else 1.0
            seg = ((seg - med) / scale).astype(np.float32, copy=False)
        else:
            seg = np.asarray(seg, dtype=np.float32)
        segments[rid] = seg

    # Intra-class pairwise DTW.
    intra_dists: List[float] = []
    # Group read indices by class.
    by_class: Dict[int, List[int]] = {}
    for i, ci in enumerate(read_class):
        if ci < 0:
            continue
        by_class.setdefault(ci, []).append(i)

    for ci, members in by_class.items():
        if len(members) < 2:
            continue
        for a_pos in range(len(members)):
            ia = members[a_pos]
            seg_a = segments.get(read_ids[ia])
            if seg_a is None:
                continue
            for b_pos in range(a_pos + 1, len(members)):
                ib = members[b_pos]
                seg_b = segments.get(read_ids[ib])
                if seg_b is None:
                    continue
                d = _dtw_distance(seg_a, seg_b, use_gpu=use_gpu)
                if not np.isfinite(d):
                    continue
                m4[ia, ib] = d
                m4[ib, ia] = d
                intra_dists.append(float(d))

    # Inter-class fill.
    if isinstance(inter_class_fill, (int, float)) and not isinstance(inter_class_fill, bool):
        fill_value = float(inter_class_fill)
    else:
        # "max" (default) or anything else: max intra-class, else 0.
        fill_value = float(max(intra_dists)) if intra_dists else 0.0

    for i in range(n):
        ci = read_class[i]
        for j in range(i + 1, n):
            cj = read_class[j]
            same_class = (ci >= 0 and cj >= 0 and ci == cj)
            if same_class:
                continue
            m4[i, j] = fill_value
            m4[j, i] = fill_value

    # Diagonal is 0 (already from np.zeros).
    return m4
