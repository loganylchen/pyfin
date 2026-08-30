"""Is a contained (nested) short transcript real, or a 5'-degradation artifact?

THE PROBLEM. In nanopore direct RNA (dRNA), sequencing proceeds 3'->5': the
poly(A) tail and 3' end are present in essentially every read, while the 5'
end is frequently missing because the RNA was degraded or the read terminated
early. So when a short transcript S is fully contained in a long transcript L,
a read that "supports S" can be indistinguishable from a read of L that simply
stopped early. If S shares L's splice chain and TES and differs ONLY by an
internal TSS, and if the degradation process is allowed an arbitrary profile,
the two explanations are **mathematically unidentifiable from single reads**.
This module therefore never forces a binary answer: every verdict is one of
``supported`` / ``unsupported`` / ``unidentifiable``.

THE STATISTIC. Raw 5'-end density is the wrong object: it decays with depth
simply because fewer reads reach further 5'. The right object is the
**conditional termination hazard** in transcript coordinates, walking in the
sequencing direction 3' -> 5':

    offset d      = distance from L's 5' end along the SPLICED transcript
    at_risk(d)    = reads that actually reached d  (5' end offset <= d)
    ends_at(d)    = reads whose 5' end lies in the bin at d
    hazard(d)     = ends_at(d) / at_risk(d)

Under a degradation-only null this hazard is a smooth, slowly varying function
of position (plus the terminal mass at d = 0, which is L's own TSS). A genuine
internal TSS at d0 injects an EXTRA spike into the hazard at d0, because reads
transcribed from S can only begin at d0 -- they are not the tail of a smooth
decay.

IDENTIFIABILITY LADDER (reported per candidate, hardest last):
  1. ``own_junction``  - S has a splice junction or exon L does not: direct
     junction-specific reads settle it; the hazard test is not needed.
  2. ``own_tes``       - S has its own 3' end: dRNA's RELIABLE end plus poly(A)
     anchors the read partition first, and the TSS test then runs within that
     3'-anchored group.
  3. ``tss_only``      - same chain, same TES, internal TSS only: the hardest
     case; only the degradation-null hazard test applies, and abstention is a
     legitimate and frequent outcome.

LEAKAGE RULES. Every input here is observable at inference time. Reference
annotation, gffcompare class codes, and expression oracles are used ONLY as
labels in the offline study, never as features. Cross-sample recurrence is a
cohort-only feature and is absent (not zero) in single-sample runs.
"""
from __future__ import annotations

import bisect
import logging
import math
import random
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

VERDICT_SUPPORTED = "supported"
VERDICT_UNSUPPORTED = "unsupported"
VERDICT_UNIDENTIFIABLE = "unidentifiable"


# --------------------------------------------------------------------------
# transcript-coordinate geometry
# --------------------------------------------------------------------------
def spliced_length(exons: Sequence[Tuple[int, int]]) -> int:
    return sum(e - s for s, e in exons)


def genomic_to_offset(
    pos: int, exons: Sequence[Tuple[int, int]], strand: str
) -> Optional[int]:
    """Spliced distance from the transcript's 5' end to genomic ``pos``.

    Returns None when ``pos`` falls outside the exons (intronic or beyond the
    model): an unmappable position must abstain, never silently clamp.
    """
    ex = sorted(exons)
    if not ex or pos < ex[0][0] or pos > ex[-1][1]:
        return None
    acc = 0
    for s, e in ex:
        if s <= pos <= e:
            fwd = acc + (pos - s)
            total = spliced_length(ex)
            return fwd if strand == "+" else total - fwd
        if pos < s:
            return None
        acc += e - s
    return None


def read_five_prime_offset(
    span: Tuple[int, int], exons: Sequence[Tuple[int, int]], strand: str
) -> Optional[int]:
    """Offset of a read's genomic 5' end in the model's transcript frame."""
    start, end = span
    five = start if strand == "+" else end
    off = genomic_to_offset(five, exons, strand)
    if off is not None:
        return off
    # 5' end outside the model (soft-clip / overhang): clamp to the nearest
    # in-model boundary only when the read clearly starts beyond the model's
    # 5' terminus, which is the "reached the very start" case.
    ex = sorted(exons)
    if strand == "+" and start < ex[0][0]:
        return 0
    if strand == "-" and end > ex[-1][1]:
        return 0
    return None


# --------------------------------------------------------------------------
# hazard profile
# --------------------------------------------------------------------------
@dataclass
class HazardProfile:
    """Binned conditional termination hazard along a transcript."""

    bin_bp: int
    ends: Dict[int, int] = field(default_factory=dict)      # bin -> terminations
    at_risk: Dict[int, int] = field(default_factory=dict)   # bin -> reached
    n_reads: int = 0

    def hazard(self, b: int) -> Optional[float]:
        risk = self.at_risk.get(b, 0)
        if risk <= 0:
            return None
        return self.ends.get(b, 0) / risk


def build_hazard_profile(
    offsets: Sequence[int], tx_len: int, *, bin_bp: int = 25
) -> HazardProfile:
    """Bin read 5' offsets into ends/at-risk counts.

    ``at_risk(b)`` counts reads that actually reached bin ``b`` -- i.e. whose
    5' offset is <= the bin's upper edge -- so the hazard is not confounded by
    the trivial fact that fewer reads reach further 5'.
    """
    prof = HazardProfile(bin_bp=bin_bp, n_reads=len(offsets))
    if not offsets:
        return prof
    srt = sorted(offsets)
    nbins = max(1, tx_len // bin_bp + 1)
    for b in range(nbins):
        lo, hi = b * bin_bp, (b + 1) * bin_bp - 1
        ends = bisect.bisect_right(srt, hi) - bisect.bisect_left(srt, lo)
        reached = bisect.bisect_right(srt, hi)
        if reached:
            prof.ends[b] = ends
            prof.at_risk[b] = reached
    return prof


def local_background_hazard(
    offsets: Sequence[int],
    tss_offset: int,
    *,
    bin_bp: int = 25,
    neighbourhood_bins: int = 12,
    exclude_bins: int = 1,
    floor: float = 1e-4,
) -> Optional[float]:
    """Robust degradation hazard from the candidate's OWN neighbourhood.

    Measured on SIRV, a single pooled hazard is not a safe null: real dRNA
    degradation is position-dependent (secondary structure, modifications,
    motor stalling), so degradation HOTSPOTS in single-isoform transcripts
    were called as TSS peaks at a 52% rate. The neighbourhood median asks the
    sharper question -- is this bin anomalous relative to the local
    degradation regime it sits in -- which a broad hotspot cannot satisfy but
    a genuine sharp TSS can.

    Returns None when the neighbourhood has too few usable bins to estimate.
    """
    if not offsets:
        return None
    srt = sorted(offsets)
    b0 = int(tss_offset) // bin_bp
    rates: List[float] = []
    for b in range(max(0, b0 - neighbourhood_bins), b0 + neighbourhood_bins + 1):
        if abs(b - b0) <= exclude_bins:
            continue          # never let the candidate peak set its own null
        if b < 1:
            continue          # bin 0 is the model's own TSS, not degradation
        lo, hi = b * bin_bp, (b + 1) * bin_bp - 1
        reached = bisect.bisect_right(srt, hi)
        if reached <= 0:
            continue
        ends = reached - bisect.bisect_left(srt, lo)
        rates.append(ends / reached)
    if len(rates) < 4:
        return None
    rates.sort()
    n = len(rates)
    median = (rates[n // 2] if n % 2 else 0.5 * (rates[n // 2 - 1] + rates[n // 2]))
    return max(median, floor)


def pooled_background_hazard(
    profiles: Iterable[HazardProfile], *, exclude_first_bins: int = 1
) -> float:
    """One pooled degradation hazard from loci with no internal TSS.

    Bins near offset 0 are excluded because that terminal mass is the model's
    OWN TSS, not degradation. A single pooled rate is deliberately simple: it
    is the conservative null (a real internal TSS must beat the average
    degradation rate), and it cannot overfit a locus it is being tested on.
    """
    ends = risk = 0
    for p in profiles:
        for b, r in p.at_risk.items():
            if b < exclude_first_bins:
                continue
            ends += p.ends.get(b, 0)
            risk += r
    if risk <= 0:
        return 0.0
    return ends / risk


# --------------------------------------------------------------------------
# the test
# --------------------------------------------------------------------------
@dataclass
class TssEvidence:
    """Per-candidate verdict with the numbers that produced it."""

    candidate_id: str
    parent_id: str
    tss_offset: int
    identifiability: str
    verdict: str
    reason: str
    n_reads_locus: int = 0
    n_at_risk: int = 0
    n_peak: int = 0
    n_upstream_of_tss: int = 0
    peak_fraction: float = 0.0
    background_hazard: float = 0.0
    expected_peak: float = 0.0
    effect_size: float = 0.0
    mixture_pi: float = 0.0
    llr: float = 0.0
    # None = no calibrated probability (structural rule or abstention). Never
    # NaN: these rows are serialised to JSON, and NaN is not valid JSON.
    p_value: Optional[float] = None
    q_value: Optional[float] = None
    # Within-parent Bonferroni-adjusted p for THIS candidate. The parent-level
    # q alone cannot certify an individual alternative: one strong alternative
    # would otherwise carry its weak siblings.
    p_within: Optional[float] = None
    selection_corrected: bool = False
    peak_mad: float = -1.0
    parent_unobserved: bool = False

    def as_row(self) -> dict:
        return {k: getattr(self, k) for k in (
            "candidate_id", "parent_id", "tss_offset", "identifiability",
            "verdict", "reason", "n_reads_locus", "n_at_risk", "n_peak",
            "n_upstream_of_tss", "peak_fraction", "background_hazard",
            "expected_peak", "effect_size", "mixture_pi", "llr", "p_value",
            "q_value", "p_within", "peak_mad", "parent_unobserved",
            "selection_corrected",
        )}


def apply_grouped_fdr(
    evidences: Sequence[TssEvidence], *, alpha: float = 0.05
) -> None:
    """Two-level FDR over whatever unit ``parent_id`` names.

    The grouping unit is deliberately NOT fixed by this function: production
    passes the parent *candidate model* id (de novo assembly has no gene
    annotation), while the SIRV study passes the GTF ``gene_id``. The caller
    owns that choice, so the name is neutral; describe the unit in the
    caller's own reporting rather than assuming a biological locus here.

    Three corrections are needed and they are different:

    1. Selection: each candidate's p-value is already Bonferroni-adjusted for
       the search space it was found in (see ``evaluate_internal_tss``).
    2. Within a parent, several alternative starts may be tested. Each
       candidate gets ``p_within = min(1, p * k)`` for the ``k`` alternatives
       tested on that parent, and the parent-level p-value is the minimum of
       those, so one parent contributes exactly one test upward.
    3. Across parents, Benjamini-Yekutieli rather than BH. BH assumes
       independence or positive dependence; alternative starts share reads,
       coverage and the degradation background, so the dependence structure is
       not guaranteed. BY is valid under ARBITRARY dependence at the cost of a
       log factor.

    A candidate stays ``supported`` only when BOTH its own ``p_within`` and its
    parent's ``q_value`` clear ``alpha``: the parent q alone would let one
    strong alternative certify its weak siblings. Structural ``own_tes``
    decisions carry ``p_value=None`` and never join this family.
    """
    tested = [e for e in evidences
              if e.verdict != VERDICT_UNIDENTIFIABLE
              and e.p_value is not None
              and math.isfinite(e.p_value)
              and e.selection_corrected]
    if not tested:
        return
    by_locus: Dict[str, List[TssEvidence]] = {}
    for e in tested:
        by_locus.setdefault(e.parent_id, []).append(e)

    loci = sorted(by_locus)
    # within-parent: Bonferroni over the alternatives actually tested, applied
    # to EVERY candidate (not just the winner)
    for lid in loci:
        k = len(by_locus[lid])
        for x in by_locus[lid]:
            x.p_within = min(1.0, x.p_value * k)
    locus_p = {lid: min(x.p_within for x in by_locus[lid]) for lid in loci}
    m = len(loci)
    harmonic = sum(1.0 / k for k in range(1, m + 1))   # BY penalty C(m)
    order = sorted(loci, key=lambda lid: locus_p[lid])
    running = 1.0
    locus_q: Dict[str, float] = {}
    for rank in range(m - 1, -1, -1):
        lid = order[rank]
        q = locus_p[lid] * m * harmonic / (rank + 1)
        running = min(running, q)
        locus_q[lid] = min(1.0, running)

    for lid in loci:
        q = locus_q[lid]
        for e in by_locus[lid]:
            e.q_value = q
            if e.verdict != VERDICT_SUPPORTED:
                continue
            if q > alpha:
                e.verdict = VERDICT_UNSUPPORTED
                e.reason = (f"parent_q={q:.4f}>fdr={alpha} "
                            f"(BY over {m} parent models)")
            elif e.p_within is not None and e.p_within > alpha:
                # the parent has a real signal, but not from THIS alternative
                e.verdict = VERDICT_UNSUPPORTED
                e.reason = (f"p_within={e.p_within:.4f}>alpha={alpha} "
                            f"(parent passed on another alternative)")


def _binom_logpmf(k: int, n: int, p: float) -> float:
    if n <= 0:
        return 0.0
    p = min(max(p, 1e-12), 1 - 1e-12)
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
            + k * math.log(p) + (n - k) * math.log1p(-p))


def evaluate_internal_tss(
    *,
    candidate_id: str,
    parent_id: str,
    offsets: Sequence[int],
    tss_offset: int,
    background_hazard: float,
    identifiability: str = "tss_only",
    bin_bp: int = 25,
    min_peak_reads: int = 3,
    min_peak_fraction: float = 0.10,
    min_at_risk: int = 10,
    n_bootstrap: int = 2000,
    alpha: float = 0.05,
    seed: int = 17,
    use_local_background: bool = True,
    neighbourhood_bins: int = 12,
    max_peak_mad: float = 12.0,
    n_eligible_bins: Optional[int] = None,
) -> TssEvidence:
    # NOTE: `n_eligible_bins` defaults to 1 (offset given a priori). Callers
    # that SEARCHED for the peak must pass the number of bins searched.
    """Degradation-null vs TSS-mixture test at one candidate TSS offset.

    H0: every read is a (possibly degraded) read of the long parent, so a
        termination in the candidate bin happens at the pooled background
        hazard.
    H1: an extra population begins exactly at ``tss_offset``, adding mass to
        that bin on top of the background.

    The p-value is an EMPIRICAL bootstrap under H0 (binomial resampling at the
    background hazard), not an asymptotic chi-square: the statistic is a
    boundary case (pi >= 0) where the asymptotic null is wrong.
    """
    ev = TssEvidence(
        candidate_id=candidate_id, parent_id=parent_id,
        tss_offset=int(tss_offset), identifiability=identifiability,
        verdict=VERDICT_UNIDENTIFIABLE, reason="not_evaluated",
        n_reads_locus=len(offsets), background_hazard=float(background_hazard),
    )
    if not offsets:
        ev.reason = "no_reads"
        return ev

    b0 = int(tss_offset) // bin_bp
    lo, hi = b0 * bin_bp, (b0 + 1) * bin_bp - 1
    srt = sorted(offsets)
    # How many distinct bins could this peak have been selected from? Used to
    # build the max-statistic null below. Defaults to the observed span of the
    # locus, which is the search space a peak-finder actually scans.
    # Search-space size for the selection correction. The DEFAULT is 1,
    # i.e. "this offset was specified in advance, no search was performed".
    # A caller that DISCOVERED the offset by scanning must declare how many
    # bins it scanned; deriving it from the observed read span instead would
    # scale the penalty with sequencing depth, which is not the search space.
    n_eligible_bins = 1 if n_eligible_bins is None else max(
        int(n_eligible_bins), 1)
    n_peak = bisect.bisect_right(srt, hi) - bisect.bisect_left(srt, lo)
    n_at_risk = bisect.bisect_right(srt, hi)
    # Reads strictly 5' of the candidate TSS can only come from the long
    # parent: they are positive proof the parent exists, and they are what
    # makes the contrast possible at all.
    n_up = bisect.bisect_left(srt, lo)

    # Prefer the candidate's own neighbourhood as the null: a global pooled
    # hazard cannot represent position-dependent degradation hotspots.
    if use_local_background:
        local = local_background_hazard(
            offsets, tss_offset, bin_bp=bin_bp,
            neighbourhood_bins=neighbourhood_bins,
        )
        if local is not None:
            background_hazard = max(background_hazard, local)
            ev.background_hazard = background_hazard

    ev.n_peak, ev.n_at_risk, ev.n_upstream_of_tss = n_peak, n_at_risk, n_up
    ev.peak_fraction = n_peak / n_at_risk if n_at_risk else 0.0
    ev.expected_peak = background_hazard * n_at_risk

    # Peak sharpness: a genuine TSS is a tight mode, a degradation hotspot is
    # broad. MAD of the in-bin 5' ends about the candidate offset.
    in_bin = [o for o in srt if lo <= o <= hi]
    if in_bin:
        devs = sorted(abs(o - int(tss_offset)) for o in in_bin)
        k = len(devs)
        ev.peak_mad = float(
            devs[k // 2] if k % 2 else 0.5 * (devs[k // 2 - 1] + devs[k // 2])
        )

    # --- abstention gates: too little evidence to distinguish anything -----
    if n_at_risk < min_at_risk:
        ev.reason = f"insufficient_depth(at_risk={n_at_risk}<{min_at_risk})"
        return ev
    if background_hazard <= 0.0:
        ev.reason = "no_background_model"
        return ev
    if n_up == 0:
        # No read extends 5' of the candidate TSS. The degradation explanation
        # REQUIRES the long parent to be present and then truncated, but here
        # the parent's 5' region is unobserved entirely, and degradation is a
        # smooth process that would not deposit every read in one bin. A sharp
        # dominant start with no upstream evidence is therefore support FOR a
        # transcript beginning at d0 (and silence about the parent), not an
        # unresolvable tie. Weak/diffuse starts still abstain.
        ev.parent_unobserved = True
        # Sharpness is checked BEFORE the parent-unobserved shortcut, else a
        # broad hump with no upstream reads bypasses the guard entirely.
        if ev.peak_mad > max_peak_mad:
            ev.verdict = VERDICT_UNSUPPORTED
            ev.reason = f"peak_too_broad(mad={ev.peak_mad:.1f}>{max_peak_mad})"
            return ev
        if n_peak >= min_peak_reads and (n_peak / n_at_risk) >= 0.5:
            ev.verdict = VERDICT_SUPPORTED
            # NOTE: this supports the SHORT model's existence only. It says
            # nothing about whether the long parent exists -- the parent's 5'
            # region is simply unobserved here.
            ev.reason = "dominant_start_no_upstream_parent_unobserved"
            ev.peak_fraction = n_peak / n_at_risk
            ev.effect_size = float("inf") if background_hazard <= 0 else (
                (n_peak - background_hazard * n_at_risk)
                / math.sqrt(max(background_hazard * n_at_risk, 1e-9))
            )
            ev.mixture_pi = 1.0
            # STRUCTURAL call, not a calibrated probability: None keeps it out
            # of the BH family (see apply_bh_fdr) and out of JSON as null.
            ev.p_value = None
        else:
            ev.reason = "no_upstream_reads_and_diffuse_start"
        return ev

    # --- mixture fit + LLR ------------------------------------------------
    ll0 = _binom_logpmf(n_peak, n_at_risk, background_hazard)
    best_pi, best_ll = 0.0, ll0
    for i in range(1, 101):
        pi = i / 100.0
        p1 = background_hazard + pi * (1.0 - background_hazard)
        ll = _binom_logpmf(n_peak, n_at_risk, p1)
        if ll > best_ll:
            best_ll, best_pi = ll, pi
    ev.mixture_pi = best_pi
    ev.llr = 2.0 * (best_ll - ll0)
    ev.effect_size = (
        (n_peak - ev.expected_peak) / math.sqrt(max(ev.expected_peak, 1e-9))
    )

    # --- empirical bootstrap p-value under H0 -----------------------------
    rng = random.Random(seed)
    ge = 0
    for _ in range(n_bootstrap):
        sim = sum(1 for _ in range(n_at_risk)
                  if rng.random() < background_hazard)
        if sim >= n_peak:
            ge += 1
    p_pointwise = (ge + 1) / (n_bootstrap + 1)

    # SELECTION CORRECTION. The candidate bin was not chosen at random: it was
    # proposed because it already looked like a peak, then tested on the very
    # reads that made it look that way, so `p_pointwise` is anti-conservative.
    #
    # A Bonferroni adjustment over the search space is used rather than a
    # simulated max-statistic. Simulating the maximum would require the risk
    # set at EVERY candidate bin; reusing this bin's `n_at_risk` for all of
    # them understates the null's spread and is itself anti-conservative.
    # Bonferroni needs only the search-space SIZE and is conservative for any
    # dependence structure between bins.
    n_bins = max(int(n_eligible_bins), 1)
    ev.p_value = min(1.0, p_pointwise * n_bins)
    ev.selection_corrected = True

    # --- verdict ----------------------------------------------------------
    if n_peak < min_peak_reads:
        ev.verdict = VERDICT_UNSUPPORTED
        ev.reason = f"peak_reads({n_peak})<{min_peak_reads}"
    elif ev.peak_mad > max_peak_mad:
        ev.verdict = VERDICT_UNSUPPORTED
        ev.reason = f"peak_too_broad(mad={ev.peak_mad:.1f}>{max_peak_mad})"
    elif ev.peak_fraction < min_peak_fraction:
        ev.verdict = VERDICT_UNSUPPORTED
        ev.reason = f"peak_fraction({ev.peak_fraction:.3f})<{min_peak_fraction}"
    elif ev.p_value > alpha:
        ev.verdict = VERDICT_UNSUPPORTED
        ev.reason = f"p={ev.p_value:.4f}>alpha={alpha}"
    else:
        ev.verdict = VERDICT_SUPPORTED
        ev.reason = "peak_exceeds_degradation_background"
    return ev


def evaluate_tes_support(
    *,
    candidate_id: str,
    parent_id: str,
    read_three_prime_offsets: Sequence[int],
    candidate_tes_offset: int,
    parent_tes_offset: int,
    window_bp: int = 25,
    min_reads: int = 3,
    min_fraction: float = 0.10,
    polya_supported: Optional[int] = None,
) -> TssEvidence:
    """Rung-2 decision: a contained model that differs by its 3' END.

    dRNA sequences 3'->5', so the 3' end is the RELIABLE end: essentially
    every read carries it and poly(A) can confirm it. When the contained
    model's start coincides with the parent's (an `own_tes` candidate whose
    TSS offset is inside the first bin), the 5' hazard test is asking the
    wrong question -- there is no internal TSS to detect. The evidence that
    settles such a model is how many reads TERMINATE at its own 3' end rather
    than the parent's.

    Unlike 5' truncation, reads do not lose their 3' end to degradation, so a
    distinct 3' cluster is positive evidence rather than an artifact.
    """
    ev = TssEvidence(
        candidate_id=candidate_id, parent_id=parent_id,
        tss_offset=int(candidate_tes_offset), identifiability="own_tes",
        verdict=VERDICT_UNIDENTIFIABLE, reason="not_evaluated",
        n_reads_locus=len(read_three_prime_offsets),
    )
    if not read_three_prime_offsets:
        ev.reason = "no_reads"
        return ev
    at_cand = sum(1 for o in read_three_prime_offsets
                  if abs(o - candidate_tes_offset) <= window_bp)
    at_parent = sum(1 for o in read_three_prime_offsets
                    if abs(o - parent_tes_offset) <= window_bp)
    total = len(read_three_prime_offsets)
    ev.n_peak = at_cand
    ev.n_at_risk = total
    ev.n_upstream_of_tss = at_parent
    ev.peak_fraction = at_cand / total if total else 0.0
    if abs(candidate_tes_offset - parent_tes_offset) <= window_bp:
        ev.reason = "candidate_and_parent_share_a_3prime_end"
        return ev
    if total < min_reads * 2:
        ev.reason = f"insufficient_depth(reads={total})"
        return ev
    if polya_supported is not None:
        ev.mixture_pi = polya_supported / at_cand if at_cand else 0.0
    if at_cand < min_reads:
        ev.verdict = VERDICT_UNSUPPORTED
        ev.reason = f"tes_reads({at_cand})<{min_reads}"
    elif ev.peak_fraction < min_fraction:
        ev.verdict = VERDICT_UNSUPPORTED
        ev.reason = f"tes_fraction({ev.peak_fraction:.3f})<{min_fraction}"
    else:
        ev.verdict = VERDICT_SUPPORTED
        ev.reason = "distinct_3prime_cluster_on_the_reliable_end"
    # STRUCTURAL threshold rule on the reliable 3' end, not a calibrated
    # probability: None excludes it from the BH family and serialises as null.
    ev.p_value = None
    # Ranking score WITHIN this rung only -- not comparable with the hazard
    # effect size of the `tss_only` rung.
    expected = total * (2 * window_bp) / max(
        abs(candidate_tes_offset - parent_tes_offset) + 2 * window_bp, 1)
    ev.expected_peak = expected
    ev.effect_size = (at_cand - expected) / math.sqrt(max(expected, 1e-9))
    return ev


def classify_identifiability(
    short_exons: Sequence[Tuple[int, int]],
    long_exons: Sequence[Tuple[int, int]],
    strand: str,
) -> str:
    """Which rung of the identifiability ladder a contained candidate sits on."""
    def chain(ex):
        ex = sorted(ex)
        return tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))

    sc, lc = chain(short_exons), chain(long_exons)
    if not set(sc) <= set(lc):
        return "own_junction"
    s, l = sorted(short_exons), sorted(long_exons)
    s_tes = s[-1][1] if strand == "+" else s[0][0]
    l_tes = l[-1][1] if strand == "+" else l[0][0]
    if s_tes != l_tes:
        return "own_tes"
    return "tss_only"


__all__ = [
    "VERDICT_SUPPORTED", "VERDICT_UNSUPPORTED", "VERDICT_UNIDENTIFIABLE",
    "HazardProfile", "TssEvidence", "build_hazard_profile",
    "pooled_background_hazard", "local_background_hazard",
    "genomic_to_offset", "read_five_prime_offset", 
    "spliced_length", "evaluate_internal_tss", "evaluate_tes_support",
    "classify_identifiability", "apply_grouped_fdr",
]
