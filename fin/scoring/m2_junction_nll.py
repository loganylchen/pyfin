"""M2 junction NLL: per-event signal distance over internal-junction diff windows.

This is the production-shaped form of the validated "version B" wobble
discriminator (see ``experiments/_gt_zval_multi.py``; SIRV interval 10582-11643,
6 wobble pairs, pooled 86.2%). The signal engine is krill ``align_read_variants``
(in-memory eventalign, ``hmm_confidence=False``).

The three essential ingredients of the metric:
  1. NON-HMM eventalign (clean per-kmer ``standardized_level`` z for every event;
     ``model_kmer == reference_kmer`` is forced, so the only signal quantity is |z|).
  2. A SHARED genomic discrimination window restricted to *internal* intron
     boundaries (excludes 5'/3'/UTR structural differences), built from the
     exon symmetric difference between candidates of the same wobble class.
  3. per-event NLL = ``0.5·z² + log(model_stdv)``, **mean over the window**
     (÷n; the sum is read-length confounded).

M2 is an M1-structural-mask refinement: only (read, candidate) cells that M1
accepts (mappy aligns with no single indel > ``M1_MAX_INDEL_BP``) are scored;
M1-rejected cells get the sentinel ``MISSING`` distance. Cells whose candidate
class has no discrimination window (singletons / no internal diff) fall back to
the M1 distance, so M2 degrades to M1 where there is no wobble competition.
"""
from __future__ import annotations

import logging
import math
from functools import lru_cache
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from fin.candidates.dataclasses import TranscriptCandidate
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit

logger = logging.getLogger(__name__)

MISSING = 1e6  # sentinel distance for M1-rejected / no-coverage cells


# ---------------------------------------------------------------------------
# Exon / coordinate helpers (ported from experiments/_gt_zval_multi.py)
# ---------------------------------------------------------------------------


@lru_cache(maxsize=16384)
def _build_exons_key(
    introns: Tuple[Tuple[int, int], ...], start: int, end: int
) -> Tuple[Tuple[int, int], ...]:
    """Cached exon geometry keyed on the hashable candidate signature.

    The exon layout is a pure function of (intron_chain, start, end), but the
    per-event NLL loop calls ``_tx2genome`` (hence ``_exon_tx_map`` ->
    ``_build_exons``) ~15M times against the SAME handful of candidates. Caching
    on the geometry tuple collapses those rebuilds to one per unique candidate.
    """
    exons: List[Tuple[int, int]] = []
    pos = start
    for s, e in sorted(introns):
        if s > pos:
            exons.append((pos, s))
        pos = max(pos, e)
    if end > pos:
        exons.append((pos, end))
    return tuple(exons)


@lru_cache(maxsize=16384)
def _exon_tx_map_key(
    introns: Tuple[Tuple[int, int], ...], start: int, end: int, strand: str
) -> Tuple[Tuple[int, int, int], ...]:
    """Cached (g_start, g_end, tx_offset) map; see :func:`_build_exons_key`."""
    exons = _build_exons_key(introns, start, end)
    order = tuple(reversed(exons)) if strand == "-" else exons
    out: List[Tuple[int, int, int]] = []
    t = 0
    for s, e in order:
        out.append((s, e, t))
        t += e - s
    return tuple(out)


def _build_exons(c: TranscriptCandidate) -> List[Tuple[int, int]]:
    """Genomic exon intervals (sorted ascending) from intron_chain + start/end."""
    return list(_build_exons_key(c.intron_chain.introns, c.start, c.end))


def _exon_tx_map(c: TranscriptCandidate) -> List[Tuple[int, int, int]]:
    """Map each genomic exon to its transcript (cDNA) offset.

    Exons are walked in transcript 5'->3' order (reversed for '-' strand), so
    the returned offsets index into ``c.sequence`` (the spliced transcript).
    Returns list of (g_start, g_end, tx_offset).
    """
    return list(_exon_tx_map_key(c.intron_chain.introns, c.start, c.end, c.strand))


def _tx2genome(c: TranscriptCandidate, p: int) -> Optional[int]:
    """Map a transcript-frame position ``p`` (into c.sequence) to genomic."""
    minus = c.strand == "-"
    for s, e, t in _exon_tx_map_key(c.intron_chain.introns, c.start, c.end, c.strand):
        if t <= p < t + (e - s):
            off = p - t
            return (e - 1 - off) if minus else (s + off)
    return None


@lru_cache(maxsize=16384)
def _tx2genome_table(
    introns: Tuple[Tuple[int, int], ...], start: int, end: int, strand: str
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, bool, int]:
    """Cached arrays for a vectorized transcript->genome projection.

    Returns ``(tx_starts, s_arr, e_arr, minus, total_len)`` where exon k spans
    tx offsets ``[tx_starts[k], tx_starts[k] + (e-s))``. Mirrors :func:`_tx2genome`
    exactly; built once per unique candidate geometry.
    """
    em = _exon_tx_map_key(introns, start, end, strand)
    if not em:
        return (np.empty(0, np.int64), np.empty(0, np.int64),
                np.empty(0, np.int64), strand == "-", 0)
    s_arr = np.array([s for s, _, _ in em], dtype=np.int64)
    e_arr = np.array([e for _, e, _ in em], dtype=np.int64)
    tx_starts = np.array([t for _, _, t in em], dtype=np.int64)
    total = int(tx_starts[-1] + (e_arr[-1] - s_arr[-1]))
    return tx_starts, s_arr, e_arr, strand == "-", total


def tx2genome_array(c: TranscriptCandidate, positions: np.ndarray) -> np.ndarray:
    """Vectorized :func:`_tx2genome` over many tx-frame positions.

    Returns an int64 array of genomic positions; entries for out-of-range tx
    positions (``_tx2genome`` would return None) are set to ``-1``. Exons tile the
    transcript contiguously in tx order, so a single ``searchsorted`` on the exon
    tx-start boundaries selects each position's exon.
    """
    tx_starts, s_arr, e_arr, minus, total = _tx2genome_table(
        c.intron_chain.introns, c.start, c.end, c.strand
    )
    p = np.asarray(positions, dtype=np.int64)
    out = np.full(p.shape, -1, dtype=np.int64)
    if tx_starts.size == 0:
        return out
    valid = (p >= 0) & (p < total)
    if not valid.any():
        return out
    pv = p[valid]
    # exon index = last boundary <= pv (boundaries are tx_starts, ascending).
    idx = np.searchsorted(tx_starts, pv, side="right") - 1
    off = pv - tx_starts[idx]
    gen = np.where(minus, e_arr[idx] - 1 - off, s_arr[idx] + off)
    out[valid] = gen
    return out


def _isub(a: Sequence[Tuple[int, int]], b: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Interval subtraction A \\ B (genomic, half-open)."""
    res: List[Tuple[int, int]] = []
    for s, e in a:
        cur = [(s, e)]
        for bs, be in b:
            nxt: List[Tuple[int, int]] = []
            for cs, ce in cur:
                if be <= cs or bs >= ce:
                    nxt.append((cs, ce))
                else:
                    if cs < bs:
                        nxt.append((cs, bs))
                    if be < ce:
                        nxt.append((be, ce))
            cur = nxt
        res.extend(cur)
    return [iv for iv in res if iv[1] > iv[0]]


def _merge(segs: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Merge overlapping/adjacent intervals."""
    out: List[Tuple[int, int]] = []
    for s, e in sorted(segs):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def _internal_bounds(c: TranscriptCandidate) -> List[int]:
    """Flat internal donor/acceptor coordinates (introns[1:-1] flattened)."""
    flat: List[int] = []
    for s, e in sorted(c.intron_chain.introns):
        flat += [s, e]
    return flat[1:-1]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def internal_diff_regions(
    class_cands: List[TranscriptCandidate],
    flank: int = 2,
) -> List[Tuple[int, int]]:
    """Internal-junction discrimination windows for one wobble class.

    For every pair of candidates in the class, take the union of their exon
    symmetric difference, restrict to the internal-intron-boundary span
    ``[min(internal_bounds), max(internal_bounds)]`` (so 5'/3' and UTR
    differences are excluded), pad each sliver by ``flank``, and merge.

    Args:
        class_cands: Candidates in a single wobble-tolerant class.
        flank: bp padding added on each side of every diff sliver.

    Returns:
        Sorted, merged list of (g_start, g_end) genomic windows. Empty when the
        class has < 2 candidates or no internal disagreement.
    """
    if len(class_cands) < 2:
        return []

    # Internal-bound span pooled across the class (candidates share intron count).
    all_bounds: List[int] = []
    for c in class_cands:
        all_bounds.extend(_internal_bounds(c))
    if not all_bounds:
        return []
    lo_b, hi_b = min(all_bounds), max(all_bounds)

    slivers: List[Tuple[int, int]] = []
    for i in range(len(class_cands)):
        ex_i = _build_exons(class_cands[i])
        for j in range(i + 1, len(class_cands)):
            ex_j = _build_exons(class_cands[j])
            sym = _isub(ex_i, ex_j) + _isub(ex_j, ex_i)
            for s, e in sym:
                if s >= lo_b and e <= hi_b:
                    slivers.append((s, e))
    if not slivers:
        return []

    padded = [(s - flank, e + flank) for s, e in slivers]
    return _merge(padded)


# ---------------------------------------------------------------------------
# Diff-region coverage gate (m2_em): genomic spans that bracket each wobbling
# junction, plus a read-straddle test. Unlike internal_diff_regions (the narrow
# sym-diff sliver) and class_junction_window_set (the per-candidate +-K scoring
# window), these spans run donor->acceptor across the whole wobbling intron so the
# coverage test asks "did the read's signal actually traverse this junction?".
# ---------------------------------------------------------------------------


def wobble_diff_spans(
    class_cands: List[TranscriptCandidate],
    flank: int = 2,
) -> List[Tuple[int, int]]:
    """Genomic spans bracketing each wobbling junction in a tie.

    Collect every intron across the tie, cluster by genomic overlap, and emit a
    flank-padded ``[min intron-start, max intron-end]`` span for each cluster that
    is NOT unanimous — i.e. the candidates disagree there, either because the
    intron boundaries wobble (donor/acceptor shift, intron retention) or because
    some candidate has no intron in the cluster at all (exon skip). Shared
    (unanimous) junctions are skipped.

    Each returned span runs from the donor side to the acceptor side of the
    wobbling intron (the intronic middle has no events by construction), so a read
    "covers" the span only when its events straddle it (see :func:`read_straddles`).

    Args:
        class_cands: Candidates in a single wobble-tolerant tie.
        flank: bp padding added on each side of every span.

    Returns:
        Sorted, merged list of (g_start, g_end) genomic spans. Empty when the tie
        has < 2 candidates or no internal disagreement.
    """
    if len(class_cands) < 2:
        return []

    per_cand = [set(c.intron_chain.introns) for c in class_cands]
    union: List[Tuple[int, int]] = sorted({iv for s in per_cand for iv in s})
    if not union:
        return []

    # Cluster union introns by genomic overlap (one cluster per logical junction).
    clusters: List[List[Tuple[int, int]]] = []
    cur: List[Tuple[int, int]] = [union[0]]
    cur_hi = union[0][1]
    for s, e in union[1:]:
        if s < cur_hi:  # overlaps the running cluster
            cur.append((s, e))
            cur_hi = max(cur_hi, e)
        else:
            clusters.append(cur)
            cur = [(s, e)]
            cur_hi = e
    clusters.append(cur)

    n = len(class_cands)
    spans: List[Tuple[int, int]] = []
    for cluster in clusters:
        cluster_set = set(cluster)
        distinct = len(cluster_set)
        contributing = sum(1 for cs in per_cand if cs & cluster_set)
        if distinct > 1 or contributing < n:  # not unanimous => wobbling
            lo = min(s for s, _ in cluster) - flank
            hi = max(e for _, e in cluster) + flank
            spans.append((lo, hi))
    return _merge(spans)


def read_straddles(event_gpos: Sequence[int], lo: int, hi: int) -> bool:
    """True iff the read's genomic event positions bracket the span ``[lo, hi]``.

    Straddling (>=1 event at <= lo AND >=1 event at >= hi) means the read's signal
    traverses the wobbling intron rather than clipping one side of it.
    """
    return any(g <= lo for g in event_gpos) and any(g >= hi for g in event_gpos)


def event_genomic_positions(res: dict, cand: TranscriptCandidate) -> List[int]:
    """Genomic positions of a read's eventalign events against ``cand``.

    Projects every event's transcript-frame position back to genomic via
    :func:`_tx2genome` (same projection :func:`_mean_nll_in_gset` uses). Used by the
    diff-region coverage gate to test whether the read straddles a wobbling junction.
    """
    pos = res.get("position", [])
    if len(pos) == 0:
        return []
    gen = tx2genome_array(cand, np.asarray(pos, dtype=np.int64))
    return gen[gen >= 0].tolist()


# ---------------------------------------------------------------------------
# NEW transcript-frame two-sided junction window (validated 99.2% pairwise;
# see experiments/_gt_zval_window_compare.py). Each candidate contributes K
# tx-bp on EACH side of its OWN wobbling junction; the genomic positions are
# UNIONed across the class so the disagreement sliver (one's exon = other's
# intron) lands inside the window -- that asymmetry IS the discriminator.
# ---------------------------------------------------------------------------


def _tx_junctions(c: TranscriptCandidate) -> List[Tuple[int, Tuple[int, int]]]:
    """Per transcript-consecutive exon pair: (J_tx, (intron_lo, intron_hi)).

    J_tx = transcript offset where the downstream exon (in tx 5'->3' order)
    begins; intron = genomic gap between the two exons (strand-agnostic).
    """
    em = _exon_tx_map(c)
    out: List[Tuple[int, Tuple[int, int]]] = []
    for k in range(len(em) - 1):
        s0, e0, _ = em[k]
        s1, e1, t1 = em[k + 1]
        lo = min(e0, e1)
        hi = max(s0, s1)
        out.append((t1, (min(lo, hi), max(lo, hi))))
    return out


def _nearest_junction(c: TranscriptCandidate, alo: int, ahi: int) -> Optional[int]:
    """Tx offset of the junction whose intron is closest to genomic [alo, ahi]."""
    js = _tx_junctions(c)
    if not js:
        return None
    ac = 0.5 * (alo + ahi)

    def dist(intron: Tuple[int, int]) -> int:
        ilo, ihi = intron
        if ihi < alo:
            return alo - ihi
        if ilo > ahi:
            return ilo - ahi
        return 0

    return min(js, key=lambda jt: (dist(jt[1]), abs(0.5 * (jt[1][0] + jt[1][1]) - ac)))[0]


def class_junction_window_set(
    class_cands: List[TranscriptCandidate],
    flank: int = 2,
    k: int = 10,
) -> set:
    """Genomic position set: K tx-bp each side of every wobbling junction.

    Uses :func:`internal_diff_regions` only to locate the genomic anchor span(s)
    of the disagreement; then for each candidate takes ``[J-K, J+K)`` in its own
    transcript frame around the junction nearest each anchor, projects to genomic
    via :func:`_tx2genome`, and unions across all candidates. Empty when the
    class has no internal discrimination region.
    """
    regions = internal_diff_regions(class_cands, flank=flank)
    if not regions:
        return set()
    gset: set = set()
    for (alo, ahi) in regions:
        for c in class_cands:
            j = _nearest_junction(c, alo, ahi)
            if j is None:
                continue
            n = len(c.sequence)
            for t in range(max(0, j - k), min(n, j + k)):
                g = _tx2genome(c, t)
                if g is not None:
                    gset.add(g)
    return gset


def _mean_nll_in_gset(
    res: dict,
    cand: TranscriptCandidate,
    gset: set,
) -> Tuple[float, int]:
    """Mean per-event NLL over events whose genomic position is in ``gset``."""
    pos = res["position"]
    z, sd = _zrecords(res)
    if len(pos) == 0:
        return float("nan"), 0
    gen = tx2genome_array(cand, np.asarray(pos, dtype=np.int64))
    total = 0.0
    n = 0
    for i in range(len(pos)):
        if not np.isfinite(z[i]):
            continue
        g = gen[i]
        if g < 0 or g not in gset:  # -1 == _tx2genome None (out of range)
            continue
        total += 0.5 * z[i] * z[i] + math.log(sd[i])
        n += 1
    if n == 0:
        return float("nan"), 0
    return total / n, n


def _zrecords(res: dict) -> Tuple[np.ndarray, np.ndarray]:
    """Per-event (z, model_stdv) from an align_read_variants record.

    ``z = (event_level_mean - model_mean) / model_stdv``. For non-HMM
    eventalign every event has ``model_kmer == reference_kmer`` and a positive
    stdv, so the per-read kmer->(mean,stdv) recovery map is a defensive no-op
    (kept for parity with the HMM-capable blueprint).
    """
    rk = res["reference_kmer"]
    mk = res["model_kmer"]
    lev = np.asarray(res["event_level_mean"], dtype=float)
    mm = np.asarray(res["model_mean"], dtype=float)
    msd = np.asarray(res["model_stdv"], dtype=float)
    n = len(rk)

    kmap: Dict[str, Tuple[float, float]] = {}
    for i in range(n):
        if mk[i] != "NNNNN" and msd[i] > 0:
            kmap.setdefault(rk[i], (mm[i], msd[i]))

    z = np.full(n, np.nan)
    sd = np.full(n, np.nan)
    for i in range(n):
        if mk[i] != "NNNNN" and msd[i] > 0:
            m, s = mm[i], msd[i]
        elif rk[i] in kmap:
            m, s = kmap[rk[i]]
        else:
            continue
        z[i] = (lev[i] - m) / s
        sd[i] = s
    return z, sd


def _mean_nll_in_window(
    res: dict,
    cand: TranscriptCandidate,
    windows: List[Tuple[int, int]],
) -> Tuple[float, int]:
    """Mean per-event NLL over events whose genomic position falls in any window."""
    pos = res["position"]
    z, sd = _zrecords(res)
    total = 0.0
    n = 0
    for i in range(len(pos)):
        if not np.isfinite(z[i]):
            continue
        gen = _tx2genome(cand, int(pos[i]))
        if gen is None:
            continue
        if not any(w0 <= gen < w1 for w0, w1 in windows):
            continue
        total += 0.5 * z[i] * z[i] + math.log(sd[i])
        n += 1
    if n == 0:
        return float("nan"), 0
    return total / n, n


def read_cand_mean_nll(
    read_id: str,
    read_seq: str,
    cand: TranscriptCandidate,
    windows: List[Tuple[int, int]],
    krill_aligner,
    mappy_aligner,
    sig_path: str,
    pore: str,
    gset: Optional[set] = None,
    use_gpu: bool = True,
    num_thread: int = 0,
) -> Tuple[float, int]:
    """Per-event mean NLL of one read against one candidate over ``windows``.

    Pipeline: mappy best hit on ``cand`` (rejecting single indels > cap via
    :func:`score_hit`) -> slice ``cand.sequence[r_st:r_en]`` -> krill
    ``align_read_variants(start=r_st)`` -> project each event back to genomic
    via :func:`_tx2genome`, keep those inside ``windows``, average
    ``0.5·z² + log(model_stdv)``.

    Returns ``(mean_nll, n_events)``; ``(nan, 0)`` when there is no usable hit,
    eventalign fails, or no event falls in the window.
    """
    import krill

    # mappy best hit by reconstructed AS (rejects > M1_MAX_INDEL_BP single indel)
    best_hit = None
    best_as = None
    for h in mappy_aligner.map(read_seq):
        v = score_hit(h)
        if v is None:
            continue
        if best_as is None or v > best_as:
            best_as, best_hit = v, h
    if best_hit is None:
        return float("nan"), 0

    try:
        recs = krill.align_read_variants(
            sig_path,
            read_id,
            {cand.candidate_id: cand.sequence[best_hit.r_st:best_hit.r_en]},
            pore=pore,
            use_gpu=use_gpu,
            num_thread=num_thread,
            aligner=krill_aligner,
            start=best_hit.r_st,
        )
    except Exception as e:  # noqa: BLE001 - krill raises broad errors
        logger.debug("align_read_variants failed for %s/%s: %s", read_id, cand.candidate_id, e)
        return float("nan"), 0

    res = {x.get("variant_label"): x for x in recs}.get(cand.candidate_id)
    if res is None or res.get("status", -1) != 0:
        return float("nan"), 0
    if gset is not None:
        return _mean_nll_in_gset(res, cand, gset)
    return _mean_nll_in_window(res, cand, windows)


# ---------------------------------------------------------------------------
# Production entry point (Stage M2-0): confine M2 to its ONLY validated niche —
# a post-hoc 2-way (small) tiebreak on an M1-already-consolidated tie set. This
# is the simultaneously-best-AS set that ``_quant_argmax_keep`` would otherwise
# blindly 1/K-split; M2 picks the single best wobble sibling instead. NOT a
# full read x candidate matrix (that scatters and loses in every prior ablation).
# ---------------------------------------------------------------------------


def m2_resolve_tie(
    read_id: str,
    read_seq: str,
    tied_cands: List[TranscriptCandidate],
    sig_path: str,
    pore: str = "rna002",
    junction_k: int = 10,
    flank: int = 2,
    krill_aligner=None,
    mappy_aligners: Optional[List] = None,
    return_scored: bool = False,
    use_gpu: bool = True,
    num_thread: int = 0,
) -> Tuple:
    """Resolve an M1 tie with the validated junction-window mean-NLL metric.

    Builds the transcript-frame two-sided junction discrimination window for the
    (small) tie set, scores each tied candidate's per-event mean NLL (non-HMM
    eventalign) over that window, and returns the lowest-NLL candidate plus the
    NLL margin to the runner-up. Lower NLL = the read's junction signal fits this
    candidate's claimed exon structure better.

    Args:
        read_id: Read identifier (for krill signal lookup).
        read_seq: Read query sequence.
        tied_cands: The simultaneously-best-AS candidates to disambiguate.
        sig_path: blow5/slow5 signal path for krill.
        pore: krill pore model.
        junction_k: tx-bp on each side of the wobbling junction (10 = SIRV sweet
            spot, validated 99.2% pairwise).
        flank: diff-window padding (bp).
        krill_aligner: optional pre-built non-HMM ``krill.Aligner`` (built lazily
            when None; pass one to reuse across reads in a loop).
        mappy_aligners: optional list aligned with ``tied_cands`` of pre-built
            ``mappy.Aligner`` (built lazily when None).
        return_scored: when True, append a third element to the returned tuple:
            the sorted (best-first) list of LOCAL indices into ``tied_cands``
            that eventalign could score in the window. Enables the caller's
            score-gated fallback split (assign the read only to scored siblings;
            empty list -> nothing scored -> caller keeps its full default split).

    Returns:
        ``(best_local_idx, margin)`` (or ``(best_local_idx, margin, scored_idxs)``
        when ``return_scored``):
          * ``best_local_idx`` indexes ``tied_cands``; ``None`` when M2 cannot
            decide (fewer than 2 candidates, no discrimination window, or no
            candidate has signal in the window) -> caller should fall back to its
            default (1/K split).
          * ``margin`` = runner-up NLL - best NLL (>= 0); ``inf`` when exactly one
            tied candidate has window signal; ``0.0`` when ``best_local_idx`` is
            None. Larger = more confident; gate with ``margin >= threshold``.
          * ``scored_idxs`` (only with ``return_scored``): local indices of the
            tied candidates with finite window NLL, sorted best (lowest NLL)
            first. Empty when none scored.
    """
    if len(tied_cands) < 2 or not read_seq:
        return (None, 0.0, []) if return_scored else (None, 0.0)

    gset = class_junction_window_set(tied_cands, flank=flank, k=junction_k)
    if not gset:
        return (None, 0.0, []) if return_scored else (None, 0.0)

    if krill_aligner is None:
        try:
            import krill

            from fin.scoring.krill_aligner import make_krill_aligner
        except Exception as exc:  # krill not importable -> abstain (no tiebreak)
            logger.warning("m2_resolve_tie: krill import failed (%s); skipping", exc)
            return (None, 0.0, []) if return_scored else (None, 0.0)

        krill_aligner, use_gpu = make_krill_aligner(
            krill, pore, use_gpu, hmm_confidence=False, num_thread=num_thread
        )
        if krill_aligner is None:
            return (None, 0.0, []) if return_scored else (None, 0.0)
    if mappy_aligners is None:
        import mappy

        preset = get_m1_preset()
        mappy_aligners = [
            mappy.Aligner(seq=c.sequence, preset=preset) if c.sequence else None
            for c in tied_cands
        ]

    scored: List[Tuple[int, float]] = []
    for idx, cand in enumerate(tied_cands):
        aln = mappy_aligners[idx] if idx < len(mappy_aligners) else None
        if aln is None:
            continue
        nll, n_ev = read_cand_mean_nll(
            read_id, read_seq, cand, [], krill_aligner, aln, sig_path, pore,
            gset=gset, use_gpu=use_gpu, num_thread=num_thread,
        )
        if n_ev > 0 and math.isfinite(nll):
            scored.append((idx, nll))

    if not scored:
        return (None, 0.0, []) if return_scored else (None, 0.0)
    scored.sort(key=lambda t: t[1])
    best_idx = scored[0][0]
    margin = float("inf") if len(scored) == 1 else scored[1][1] - scored[0][1]
    if return_scored:
        return best_idx, margin, [idx for idx, _ in scored]
    return best_idx, margin
