"""5'-TSS short-isoform recovery (post-EM).

Generation folds every exact contiguous sub-chain into its container, so a genuine
SHORT isoform that is an exact 5'-truncation of a longer one disappears as a candidate
and its reads are pooled into the long member. This module brings it back AFTER the
per-cluster EM has confirmed the long member: the folded sub-chains are dormant SHADOW
proposals (``GenCandidate.folded``); each proposes a 5' transcription start (its
truncation boundary). We test each proposal against the long member's read 5'-end
distribution and activate only those that show a SHARP peak -- a real TSS -- rather than
the smooth ramp of dRNA 5'-degradation.

The test is deliberately conservative and signal-free (structure + read-end pileup only,
per the validated 5'-TSS-pileup finding). Direct peak-excess support only; a full
degradation/abundance model is deferred. Only 5' TSS peaks for now (dRNA's dominant
gap); 3' polyA recovery is a later addition.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple


@dataclass
class RecoveredPeak:
    """An accepted 5'-TSS proposal: the peak position and its read support."""
    pos: int            # genomic 5' TSS position (the proposal that passed)
    support: float      # EM-weighted reads whose 5' end falls in +-window of pos
    excess: float       # support minus the local degradation background (peak excess)
    oe: float           # observed / expected (enrichment over flank background)
    concentration: float  # fraction of +-(flank_mult*w) reads that lie in +-w


def _weighted_in(ends: Sequence[Tuple[int, float]], center: int, half: int) -> float:
    """Total weight of 5'-end positions within ``+-half`` of ``center``."""
    tot = 0.0
    for pos, w in ends:
        if abs(pos - center) <= half:
            tot += w
    return tot


def evaluate_tss_peak(
    ends: Sequence[Tuple[int, float]],
    pos: int,
    *,
    window: int = 20,
    flank_mult: int = 3,
    min_oe: float = 4.0,
    min_concentration: float = 0.6,
    min_reads: float = 8.0,
    min_frac: float = 0.03,
    total_mass: float = 0.0,
) -> RecoveredPeak | None:
    """Peak-vs-ramp test for one proposed 5'-TSS position.

    ``ends`` are (genomic 5'-end position, EM weight) for the containment family's
    reads. A proposal at ``pos`` is accepted only if the read-end pileup in ``+-window``
    is a sharp peak over the local degradation background (estimated from the flanks
    ``window < |d| <= flank_mult*window``), meeting ALL of:
      * support (weight in +-window) >= max(``min_reads``, ``min_frac`` * family mass),
      * enrichment O/E >= ``min_oe`` (observed vs flank-density extrapolation),
      * concentration >= ``min_concentration`` (share of +-(flank_mult*w) inside +-w).
    Returns a :class:`RecoveredPeak` (with peak-excess support) or None if rejected.

    Flank-based background is direction-agnostic (a real TSS stands above the local
    ramp on either side), so no strand handling is needed here.
    """
    w_in = _weighted_in(ends, pos, window)
    w_wide = _weighted_in(ends, pos, flank_mult * window)
    w_flank = w_wide - w_in
    flank_bp = 2 * (flank_mult - 1) * window            # total flank width (bp)
    bg_density = (w_flank / flank_bp) if flank_bp > 0 else 0.0
    expected = bg_density * (2 * window)                # background reads in +-window
    oe = (w_in / expected) if expected > 0 else float("inf")
    concentration = (w_in / w_wide) if w_wide > 0 else 0.0
    support_floor = max(min_reads, min_frac * total_mass)
    if w_in >= support_floor and oe >= min_oe and concentration >= min_concentration:
        return RecoveredPeak(
            pos=pos, support=w_in, excess=max(0.0, w_in - expected),
            oe=oe, concentration=concentration)
    return None


def recover_5p_peaks(
    ends: Sequence[Tuple[int, float]],
    proposals: Sequence[int],
    *,
    window: int = 20,
    flank_mult: int = 3,
    min_oe: float = 4.0,
    min_concentration: float = 0.6,
    min_reads: float = 8.0,
    min_frac: float = 0.03,
    merge_bp: int = 20,
) -> List[RecoveredPeak]:
    """Test every candidate 5'-TSS ``proposals`` against the read-end pileup ``ends``.

    ``proposals`` are candidate genomic TSS positions (from folded shadow boundaries).
    Returns the accepted peaks, de-duplicated: proposals within ``merge_bp`` collapse
    to the single strongest (highest support) peak, so a nested short isoform is not
    counted under several near-identical proposals. Deterministic (proposals are sorted;
    ties break on position).
    """
    total_mass = sum(w for _, w in ends)
    hits: List[RecoveredPeak] = []
    for p in sorted(set(proposals)):
        r = evaluate_tss_peak(
            ends, p, window=window, flank_mult=flank_mult, min_oe=min_oe,
            min_concentration=min_concentration, min_reads=min_reads,
            min_frac=min_frac, total_mass=total_mass)
        if r is not None:
            hits.append(r)
    # merge near-duplicate accepted peaks (keep the strongest, then lowest pos).
    hits.sort(key=lambda r: (-r.support, r.pos))
    kept: List[RecoveredPeak] = []
    for r in hits:
        if all(abs(r.pos - k.pos) > merge_bp for k in kept):
            kept.append(r)
    kept.sort(key=lambda r: r.pos)
    return kept
