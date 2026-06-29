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


def support_gate_drops(cands, ties_by_read, nlls_by_read, tie_ok: bool = True):
    """Columns to drop by the M2/M1 read-support gate.

    A multi-exon candidate (any source except fusion; mono-exon exempt) is KEPT
    iff it earns >=1 read's support:
      (a) it is some read's M1 SOLE best-AS — the read's AS tie set is exactly
          this one candidate (unique mappy-AS winner), OR
      (b) it is some read's M2-best — lowest junction-NLL among that read's
          scored tie candidates (``tie_ok`` True counts ties for the min;
          False requires a strict unique minimum).
    Otherwise the candidate is dropped. A jitter-corrupted GTF junction wins
    neither (it loses the M2 signal contest to the true-junction candidate and
    never holds a read's sole AS), whereas a genuine isoform always earns one or
    the other — so this drops corrupted-annotation shadows without an annotation
    oracle.

    Args:
        cands: list of TranscriptCandidate (column order = list/index order).
        ties_by_read: read index -> list of candidate columns in its AS tie set.
        nlls_by_read: read index -> {candidate column -> junction NLL}.
        tie_ok: count candidates tied for the lowest NLL as M2-best (recall-safe).

    Returns:
        set of candidate columns to drop.
    """
    # No support evidence at all (e.g. no signal backend AND the diff-cover gate
    # off, so neither ties nor NLLs were populated): the gate has nothing to judge
    # on, so it must NOT drop anything (dropping every candidate would be wrong).
    if not ties_by_read and not nlls_by_read:
        return set()
    supported: set = set()
    for tie in ties_by_read.values():
        if tie is not None and len(tie) == 1:
            supported.add(tie[0])
    for nlls in nlls_by_read.values():
        if not nlls:
            continue
        best = min(nlls.values())
        winners = [j for j, v in nlls.items() if v <= best]
        if tie_ok:
            supported.update(winners)
        elif len(winners) == 1:
            supported.add(winners[0])
    drops: set = set()
    for j, c in enumerate(cands):
        if c.source == "fusion":
            continue
        if len(c.intron_chain.introns) < 1:
            continue  # mono-exon exempt (no junction to validate)
        if j not in supported:
            drops.add(j)
    return drops


def structural_wobble_clusters(cands, bp: int, cassette_max_exon_bp: int = 0):
    """Union-find clusters of candidates that differ only by junction wobble.

    Two multi-exon candidates join iff same chrom/strand and EITHER:
      (a) same intron count and every intron boundary within ``bp`` of the other's
          (classic wobble), OR
      (b) intron counts differ by 1 AND the chains align as a CASSETTE-skip pair:
          the longer chain has one extra exon (< ``cassette_max_exon_bp``) flanked
          by two introns whose outer donor/acceptor match the shorter chain's
          spanning intron within ``bp``, with every other intron matching within
          ``bp``. Handles minimap2's small-exon-skip failure mode where one read
          population aligned through the cassette exon and another collapsed it
          into one long intron. Disabled when ``cassette_max_exon_bp <= 0``.

    Single-exon candidates are never clustered. ALL sources (gtf/novel/fusion)
    participate, so a high-abundance true isoform (often a GTF passthrough) anchors
    its +-bp novel shadows. Pure structure — no GTF "rightness" judgement.

    Args:
        cands: list of TranscriptCandidate (column order = list order).
        bp: per-junction tolerance (each donor & acceptor within +-bp).
        cassette_max_exon_bp: max length (bp) of an extra cassette exon that
            triggers a K vs K-1 cluster join (0 disables the cassette extension).

    Returns:
        dict root_col -> [member cols] (only the connected-components map; callers
        filter size>=2).
    """
    n = len(cands)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    # bucket by (chrom, strand, intron_count) so we only compare plausible pairs.
    from collections import defaultdict
    buckets = defaultdict(list)
    introns = []
    for j, c in enumerate(cands):
        iv = tuple(sorted(c.intron_chain.introns))
        introns.append(iv)
        if len(iv) >= 1:  # multi-exon only
            buckets[(c.chrom, c.strand, len(iv))].append(j)

    def wobble(a, b):
        ia, ib = introns[a], introns[b]
        return all(abs(s1 - s2) <= bp and abs(e1 - e2) <= bp
                   for (s1, e1), (s2, e2) in zip(ia, ib))

    def cassette(a_idx, b_idx):
        """True iff ``a`` (K introns) is a cassette-skip sibling of ``b`` (K-1).

        Aligns the chains: there must be exactly one position ``i`` in ``a``
        where ``a[i]`` and ``a[i+1]`` together replace ``b[i]`` (the long span),
        with the in-between exon shorter than ``cassette_max_exon_bp``, and
        every other intron within ``bp``.
        """
        A, B = introns[a_idx], introns[b_idx]
        for i in range(len(A) - 1):
            # prefix [0..i): A[k] matches B[k] within wobble
            if any(abs(A[k][0] - B[k][0]) > bp or abs(A[k][1] - B[k][1]) > bp
                   for k in range(i)):
                continue
            # cassette pivot: outer ends of (A[i], A[i+1]) match B[i]
            if abs(A[i][0] - B[i][0]) > bp:
                continue
            if abs(A[i+1][1] - B[i][1]) > bp:
                continue
            # suffix: A[i+2+k] matches B[i+1+k]
            rest = len(B) - i - 1
            if any(abs(A[i+2+k][0] - B[i+1+k][0]) > bp or
                   abs(A[i+2+k][1] - B[i+1+k][1]) > bp
                   for k in range(rest)):
                continue
            cas_size = A[i+1][0] - A[i][1]
            if 0 < cas_size < cassette_max_exon_bp:
                return True
        return False

    # (a) same-intron-count wobble within each bucket
    for cols in buckets.values():
        for i in range(len(cols)):
            for k in range(i + 1, len(cols)):
                a, b = cols[i], cols[k]
                if find(a) != find(b) and wobble(a, b):
                    parent[find(a)] = find(b)

    # (b) cassette-skip across adjacent intron-count buckets at the same locus
    if cassette_max_exon_bp > 0:
        by_loc = defaultdict(dict)
        for (chrom, strand, k), cols in buckets.items():
            by_loc[(chrom, strand)][k] = cols
        for counts in by_loc.values():
            for k, A_cols in counts.items():
                B_cols = counts.get(k - 1)
                if not B_cols:
                    continue
                for a in A_cols:
                    for b in B_cols:
                        if find(a) != find(b) and cassette(a, b):
                            parent[find(a)] = find(b)

    out = defaultdict(list)
    for col in {c for cols in buckets.values() for c in cols}:
        out[find(col)].append(col)
    return out


def wobble_shadow_drops(cands, em_ab, clusters, frac, gtf_drop_enabled=True,
                        observed_jct=None, jct_tol=2, gtf_min_jct_reads=1):
    """Columns to drop as wobble shadows within structural clusters.

    Within each cluster (>=2 members) the highest-EM-abundance candidate is the
    anchor. A non-anchor sibling whose EM abundance is below ``frac`` * anchor is a
    candidate shadow, resolved by SOURCE:

      - NOVEL siblings are dropped on abundance alone (an unannotated low-support
        near-duplicate of the anchor is a read-mapping artefact — the original
        recheck behaviour).
      - GTF siblings are dropped only when ``gtf_drop_enabled`` AND a direct
        read-support guard fires: the sibling's DISTINGUISHING junction(s) — those
        not within ``jct_tol`` bp of any anchor junction — carry fewer than
        ``gtf_min_jct_reads`` directly-observed reads (``observed_jct``). Rationale:
        an annotated isoform is a documented hypothesis, so abundance alone must NOT
        delete it; only the data totally failing to traverse its specific junction
        (a jittered / phantom passthrough) justifies the drop. A genuine annotated
        isoform a few bp from a stronger neighbour (e.g. SIRV606: 419 reads on its
        own donor) keeps its reads and survives. This is symmetric in the anchor's
        source — a GTF is judged by ITS OWN read support, not by whether the anchor
        is novel.
      - FUSION siblings are never dropped.

    When ``observed_jct`` is None the read-support guard cannot be evaluated, so GTF
    siblings are conservatively KEPT (never dropped) — recall-safe fallback.

    Args:
        cands: list of TranscriptCandidate (column order = em_ab index order).
        em_ab: per-candidate EM abundance (column sum of R), index-aligned to cands.
        clusters: dict root_col -> [member cols] from structural_wobble_clusters.
        frac: shadow threshold (sibling dropped when em_ab/anchor_ab < frac).
        gtf_drop_enabled: allow GTF siblings to be dropped at all.
        observed_jct: optional STRAND-KEYED dict {strand: {(donor, acceptor): n_reads}}
            of intron junctions directly observed in the interval's read CIGARs.
            Strand-keyed so an antisense read at the same coordinates cannot falsely
            support a candidate on the other strand.
        jct_tol: bp tolerance matching a candidate junction to an observed one.
        gtf_min_jct_reads: a GTF sibling's distinguishing junction must have FEWER
            than this many observed reads to be droppable (1 == drop only zero-support).

    Returns:
        set of candidate columns to drop.
    """
    def jct_support(j, strand):
        counter = observed_jct.get(strand, {}) if observed_jct else {}
        d, a = j
        return sum(n for (dd, aa), n in counter.items()
                   if abs(dd - d) <= jct_tol and abs(aa - a) <= jct_tol)

    drops: set = set()
    for cols in clusters.values():
        if len(cols) < 2:
            continue
        anchor_j = max(cols, key=lambda j: em_ab[j])
        anchor_ab = em_ab[anchor_j]
        if anchor_ab <= 0.0:
            continue
        anchor_introns = cands[anchor_j].intron_chain.introns
        for j in cols:
            if j == anchor_j:
                continue
            src = cands[j].source
            if src == "fusion":
                continue
            if (em_ab[j] / anchor_ab) >= frac:
                continue
            if src == "novel":
                drops.add(j)
                continue
            # GTF sibling: needs the gate AND the direct read-support guard.
            if not gtf_drop_enabled:
                continue
            if observed_jct is None:
                continue  # cannot judge support -> keep (recall-safe)
            cj = cands[j].intron_chain.introns
            diff = [x for x in cj
                    if not any(abs(x[0] - y[0]) <= jct_tol and abs(x[1] - y[1]) <= jct_tol
                               for y in anchor_introns)]
            if not diff:
                # GTF shares every junction with the anchor (same intron chain,
                # differing only in 5'/3' ends — a real alt-TSS/TES isoform whose
                # junctions are all read-supported via the anchor). No distinguishing
                # junction => no failed-support evidence => keep (recall-safe).
                continue
            strand = cands[j].strand
            if min(jct_support(x, strand) for x in diff) < gtf_min_jct_reads:
                drops.add(j)
    return drops


def containment_shadow_drops(cands, em_ab, *, tol_bp, min_ratio, exclude=None):
    """Map NOVEL 5'-truncation containment shadows to the longer parent to fold into.

    A candidate ``a`` is a 5'-truncation shadow of a longer candidate ``b`` iff
    they share the same chrom/strand and ``a``'s intron chain is a STRICT, exact
    3' SUFFIX of ``b``'s chain (i.e. ``a`` has fewer introns, all equal to ``b``'s
    downstream-most introns), their 3' termini agree within ``tol_bp``, and ``a``'s
    5' end is interior to ``b`` (``a`` starts downstream). Strand sets which genomic
    end is "3'":

      - "+" strand: 3' = high coordinate. Shared suffix = ``b``'s LAST introns;
        3' terminus = ``end``; ``a.start > b.start``.
      - "-" strand: 3' = low coordinate. Shared suffix = ``b``'s FIRST introns;
        3' terminus = ``start``; ``a.end < b.end``.
      - ".": unstranded -> never folded (conservative).

    Only NOVEL candidates are foldable shadows; a parent may be gtf or novel (a gtf
    parent absorbs a novel truncation, but a gtf candidate is never itself folded
    away). A shadow is emitted only when ``em_ab[a] <= em_ab[b] * min_ratio`` (the
    shadow is the minor member). Among eligible parents the highest-EM-abundance one
    is chosen, ties broken by a RUN-STABLE structural key (start, end, intron chain,
    source: gtf<novel) -- NEVER by ``candidate_id`` (novel ids are random uuids).

    The match is exact on the shared suffix introns (no per-junction wobble), so a
    genuine exon-skipping or alternative-internal-splicing isoform -- whose chain is
    NOT an exact suffix -- is never folded. NOTE (recall caveat): an exact 3'-suffix
    match also fits a genuine low-abundance alternative-TSS isoform that begins at a
    downstream exon; this function cannot distinguish that from a 5'-truncation
    artifact. It is therefore a precision heuristic, not recall-safe by construction
    -- the caller keeps it default-off.

    Chained containment (a in b in c) is resolved to the TERMINAL longest parent, so
    every returned value is a candidate that is not itself a shadow. Candidates in
    ``exclude`` (already dropped by another lever) are eligible neither as shadow nor
    as parent, so a fold never targets a doomed candidate.

    Args:
        cands: list of TranscriptCandidate (column order = em_ab index order).
        em_ab: per-candidate EM abundance (column sum of R), index-aligned to cands.
        tol_bp: bp tolerance for matching the shadow's 3' terminus to the parent's.
        min_ratio: fold only when em_ab[shadow] <= em_ab[parent] * min_ratio.
        exclude: optional set of columns to skip entirely (e.g. already dropped).

    Returns:
        dict {shadow_col: parent_col}, parent_col resolved to the terminal longest
        parent. Empty when nothing folds.
    """
    exclude = exclude or set()
    from collections import defaultdict

    chains = [tuple(sorted(c.intron_chain.introns)) for c in cands]
    buckets = defaultdict(list)
    for j, c in enumerate(cands):
        if j in exclude:
            continue
        if len(chains[j]) >= 1:  # multi-exon only
            buckets[(c.chrom, c.strand)].append(j)

    def is_trunc(a, b):
        """True iff ``a`` is a pure 5'-truncation strictly contained in ``b``."""
        ca, cb = chains[a], chains[b]
        if len(ca) >= len(cb):
            return False
        strand = cands[a].strand
        if strand == "+":
            if ca != cb[-len(ca):]:
                return False
            if abs(cands[a].end - cands[b].end) > tol_bp:
                return False
            return cands[a].start > cands[b].start
        if strand == "-":
            if ca != cb[:len(ca)]:
                return False
            if abs(cands[a].start - cands[b].start) > tol_bp:
                return False
            return cands[a].end < cands[b].end
        return False  # unstranded: never fold

    def pkey(j):
        c = cands[j]
        return (c.start, c.end, chains[j], 0 if c.source == "gtf" else 1)

    raw: dict = {}  # shadow -> immediate parent
    for cols in buckets.values():
        for a in cols:
            if cands[a].source != "novel":
                continue  # only novels are foldable shadows
            parents = [
                b for b in cols
                if b != a and is_trunc(a, b) and em_ab[a] <= em_ab[b] * min_ratio
            ]
            if not parents:
                continue
            # highest EM abundance wins; ties -> smallest structural key.
            raw[a] = min(parents, key=lambda b: (-em_ab[b], pkey(b)))

    # Resolve chains (a->b->c) to the terminal longest parent.
    def terminal(x):
        seen = set()
        while x in raw and x not in seen:
            seen.add(x)
            x = raw[x]
        return x

    return {s: terminal(p) for s, p in raw.items() if terminal(p) != s}


def junction_support_drops(cands, observed_jct, *, min_reads, tol):
    """Columns to drop: NOVEL multi-exon candidates with an under-supported junction.

    Every junction (intron) of a NOVEL multi-exon candidate must be directly
    spliced by ``>= min_reads`` reads observed in the interval — i.e. a novel
    junction has to be carried by at least N independent reads, not just one. The
    support count for a junction (donor, acceptor) is the number of observed read
    introns within ``tol`` bp of BOTH boundaries (``observed_jct`` is the
    STRAND-KEYED ``{strand: {(donor, acceptor): n_reads}}`` map built from primary
    -read CIGARs). If ANY of a candidate's junctions falls below ``min_reads`` the
    whole candidate is dropped.

    GTF-passthrough (``source="gtf"``), fusion, and single-exon (mono, no
    junctions) candidates are EXEMPT. ``min_reads <= 1`` (a junction trivially has
    >= 1 read once assembled) or an empty/None ``observed_jct`` disables the gate
    (returns an empty set — byte-identical / recall-safe).

    RECALL NOTE: a genuine low-abundance novel isoform whose junction is seen by a
    single read is dropped — this is the intended precision/recall trade ("a novel
    junction needs >= 2 reads"). Keep ``min_reads`` small and tune on real data.
    """
    if min_reads <= 1 or not observed_jct:
        return set()

    def support(d, a, strand):
        counter = observed_jct.get(strand, {})
        return sum(n for (dd, aa), n in counter.items()
                   if abs(dd - d) <= tol and abs(aa - a) <= tol)

    drops: set = set()
    for j, c in enumerate(cands):
        if c.source != "novel":
            continue
        introns = c.intron_chain.introns
        if len(introns) < 1:
            continue  # mono: no junctions to support
        strand = c.strand
        if any(support(d, a, strand) < min_reads for (d, a) in introns):
            drops.add(j)
    return drops


def dominant_junction_drops(cands, observed_jct, *, min_reads, window, tol=2):
    """Columns to drop: NOVEL candidates with a weak OR non-dominant junction.

    The "junction-first" gate: a candidate junction (donor, acceptor) is VALID
    iff BOTH
      - it has ``>= min_reads`` directly-observed reads within ``tol`` bp
        (``observed_jct`` = strand-keyed ``{strand: {(d,a): n}}`` from CIGARs), AND
      - it is LOCALLY DOMINANT: no DIFFERENT observed junction within ``window`` bp
        (beyond the ``tol`` "same-junction" zone) carries STRICTLY more reads.
    A NOVEL multi-exon candidate is dropped if ANY of its junctions is invalid.

    The dominance clause is what a pure read-count gate lacks: a multi-read wobble
    shadow (its junction mis-placed a few bp off the true site) has read support
    but LOSES to the stronger true junction a few bp away -> dropped. A genuine
    alternative splice site that is the local winner (or ties) survives.

    GTF/fusion/mono exempt. ``min_reads <= 0`` or empty ``observed_jct`` disables
    (returns an empty set — byte-identical).

    RECALL NOTE: two genuine, comparably-supported alternative splice sites closer
    than ``window`` bp could demote the weaker real one. Keep ``window`` modest and
    tune on real data; this is the precision/recall knob.
    """
    if min_reads <= 0 or not observed_jct:
        return set()

    def support(d, a, strand):
        c = observed_jct.get(strand, {})
        return sum(n for (dd, aa), n in c.items()
                   if abs(dd - d) <= tol and abs(aa - a) <= tol)

    def competitor(d, a, strand):
        c = observed_jct.get(strand, {})
        best = 0
        for (dd, aa), n in c.items():
            near = abs(dd - d) <= window and abs(aa - a) <= window
            same = abs(dd - d) <= tol and abs(aa - a) <= tol
            if near and not same and n > best:
                best = n
        return best

    drops: set = set()
    for j, c in enumerate(cands):
        if c.source != "novel":
            continue
        introns = c.intron_chain.introns
        if len(introns) < 1:
            continue
        strand = c.strand
        for (d, a) in introns:
            s = support(d, a, strand)
            if s < min_reads or s < competitor(d, a, strand):
                drops.add(j)
                break
    return drops


def gtf_guard_needed(cands, em_ab, clusters, frac):
    """True iff any cluster has a non-anchor GTF sibling below ``frac`` * anchor.

    Cheap pre-check (no read I/O) so the caller only builds the observed-junction
    map — which reopens the BAM — when a clustered low-abundance GTF candidate
    could actually be dropped by the read-support guard. Mirrors the anchor/frac
    semantics of :func:`wobble_shadow_drops`.
    """
    for cols in clusters.values():
        if len(cols) < 2:
            continue
        anchor_j = max(cols, key=lambda j: em_ab[j])
        anchor_ab = em_ab[anchor_j]
        if anchor_ab <= 0.0:
            continue
        for j in cols:
            if j == anchor_j:
                continue
            src = cands[j].source
            if src in ("novel", "fusion"):
                continue
            if (em_ab[j] / anchor_ab) < frac:
                return True
    return False


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
