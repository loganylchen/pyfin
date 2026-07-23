"""Per-cluster read assignment (the "compute inside a cluster" stage of chain-cluster
generation).

Given ONE generation cluster -- its member candidate chains and the reads pooled into
it -- decide, WITHOUT ever comparing a read to a candidate in another cluster:

  1. Which member each read supports (read assignment), and
  2. Which members earn enough independent support to survive as real isoforms.

The policy here is pure/deterministic and signal-agnostic; the two mechanisms it needs
-- mappy alignment scores (M1) and the summed-LLR junction-signal tiebreak (M2) -- are
injected by the caller (the runner), so this module has no krill/mappy import and is
unit-testable on synthetic inputs.

Assignment builds an M3 read x member compatibility matrix, then EM quantifies it:
  * M1 within-cluster: each read is scored ONLY against this cluster's members; its
    best-AS tie set is computed from those scores.
  * tie == 1  -> the read is COMPATIBLE with that one member (a certain anchor).
  * tie >= 2:
      - the read STRADDLES a junction that differs among the tied members AND M2
        (summed-LLR) returns a STRICT best (positive margin -- not an NLL tie)
        -> compatible with only that signal winner (a certain anchor). No margin
        threshold: any strict winner counts.
      - otherwise (no straddle, or M2 is an NLL tie / abstains) -> the read stays
        AMBIGUOUS: compatible with ALL its tied members (initial 1/k each).
  * EM: certain reads (compat set size 1) anchor member abundance every iteration;
    ambiguous reads redistribute among their tied members in proportion to the
    learned abundances (RSEM/kallisto-style). No hard "main peak" routing.
  * A member SURVIVES iff its EM abundance >= ``min_support``.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Set, Tuple

Chain = Tuple[Tuple[int, int], ...]

# M2 resolver contract: given (read_id, [tied member indices]) return
# (winner_member_idx, margin) where winner_member_idx is one of the passed indices and
# margin is the summed-LLR gap to the runner-up (> 0 == strict winner), or None to
# abstain. A non-positive margin is treated as "no strict winner" (read stays ambiguous).
M2Resolver = Callable[[str, List[int]], Optional[Tuple[int, float]]]


@dataclass
class ClusterAssignment:
    """Result of assigning one cluster's reads to its members."""
    # read_id -> {member_idx: weight}; weights for a read sum to 1.0 (empty = dropped).
    weights: Dict[str, Dict[int, float]] = field(default_factory=dict)
    survivors: Set[int] = field(default_factory=set)
    # diagnostics
    n_unique: int = 0
    n_m2_called: int = 0
    n_m2_confident: int = 0
    n_containment: int = 0
    n_deferred: int = 0


def _differing_coords(members: Sequence[Chain], tied: Sequence[int]) -> Set[int]:
    """Donor/acceptor coords that are NOT shared by every tied member."""
    per = [set(members[i]) for i in tied]
    common: Set[Tuple[int, int]] = set(per[0])
    for s in per[1:]:
        common &= s
    coords: Set[int] = set()
    for s in per:
        for (a, b) in (s - common):
            coords.add(a)
            coords.add(b)
    return coords


def _straddles_difference(
    span: Tuple[int, int], members: Sequence[Chain], tied: Sequence[int]
) -> bool:
    """True iff the read's genomic span covers a whole intron that differs among the
    tied members -- i.e. the read actually crosses a discriminating junction."""
    rst, ren = span
    per = [set(members[i]) for i in tied]
    common: Set[Tuple[int, int]] = set(per[0])
    for s in per[1:]:
        common &= s
    for s in per:
        for (a, b) in (s - common):
            if rst <= a and b <= ren:
                return True
    return False


def _em_quantify(
    compat: Dict[str, List[int]],
    n_members: int,
    max_iter: int,
    tol: float,
) -> Tuple[List[float], Dict[str, Dict[int, float]]]:
    """EM over read compatibility sets. Deterministic (fixed 1/k init, fixed order).

    Each read's mass is confined to its compatibility set and split by the current
    member abundances; a size-1 (certain) set always contributes 1.0 to its member.
    Returns (per-member abundance, per-read responsibilities).
    """
    abundance = [0.0] * n_members
    for cs in compat.values():
        w = 1.0 / len(cs)
        for m in cs:
            abundance[m] += w
    for _ in range(max_iter):
        nxt = [0.0] * n_members
        for cs in compat.values():
            tot = 0.0
            for m in cs:
                tot += abundance[m]
            if tot <= 0.0:                       # all-zero abundance -> uniform
                w = 1.0 / len(cs)
                for m in cs:
                    nxt[m] += w
            else:
                for m in cs:
                    nxt[m] += abundance[m] / tot
        delta = 0.0
        for m in range(n_members):
            delta += abs(nxt[m] - abundance[m])
        abundance = nxt
        if delta < tol:
            break
    weights: Dict[str, Dict[int, float]] = {}
    for rid, cs in compat.items():
        tot = 0.0
        for m in cs:
            tot += abundance[m]
        if tot <= 0.0:
            w = 1.0 / len(cs)
            weights[rid] = {m: w for m in cs}
        else:
            weights[rid] = {m: abundance[m] / tot for m in cs}
    return abundance, weights


def assign_cluster_reads(
    members: Sequence[Chain],
    as_rows: Dict[str, Dict[int, float]],
    spans: Dict[str, Tuple[int, int]],
    *,
    m2_resolve: Optional[M2Resolver] = None,
    min_support: float = 1.0,
    m1_tie_margin: float = 1e-9,
    em_max_iter: int = 200,
    em_tol: float = 1e-9,
) -> ClusterAssignment:
    """Assign one cluster's reads to its members (M3 compat matrix -> EM) + survivors.

    Args:
        members: member intron chains (index == member id used everywhere else).
        as_rows: read_id -> {member_idx: mappy AS} for members THIS read aligned to
            (already restricted to this cluster's members by the caller).
        spans: read_id -> (genomic_start, genomic_end) of the read's alignment.
        m2_resolve: callback for straddling ties; returns (winner_member_idx, margin)
            or None to abstain. A STRICT winner (margin > 0) makes the read certain;
            an NLL tie / abstain leaves it ambiguous. When None, straddling ties stay
            ambiguous.
        min_support: a member survives iff its EM abundance >= this.
        m1_tie_margin: AS margin for the best-AS tie set -- a read is unique-best only
            when its top member beats the runner-up by MORE than this many AS points;
            within the margin the members are TIED (go to M2 / ambiguous 1/k). A 2-3 bp
            wobble shifts AS by only ~<=20, so a margin here stops a near-tie read from
            spuriously anchoring a wobble shadow (measured: wobble AS-gaps concentrate
            <=20, genuinely-different candidates differ by hundreds).
        em_max_iter / em_tol: EM iteration cap / L1-convergence tolerance.

    Returns:
        ClusterAssignment with per-read EM responsibilities and the surviving members.
    """
    res = ClusterAssignment()
    n = len(members)

    # --- Build the M3 compatibility sets (which members each read can come from). ---
    compat: Dict[str, List[int]] = {}
    for rid, row in as_rows.items():
        if not row:
            continue
        best = max(row.values())
        if best <= 0:
            continue
        tie = sorted(j for j, v in row.items() if v >= best - m1_tie_margin)
        if len(tie) == 1:                        # M1 certain
            compat[rid] = [tie[0]]
            res.n_unique += 1
            continue
        span = spans.get(rid)
        straddles = span is not None and _straddles_difference(span, members, tie)
        decided = None
        if straddles and m2_resolve is not None:
            res.n_m2_called += 1
            out = m2_resolve(rid, tie)
            # strict best: a finite winner with a POSITIVE margin (not an NLL tie).
            if out is not None and out[0] is not None and out[1] > 0:
                decided = out[0]
        if decided is not None:                  # M2 certain
            compat[rid] = [decided]
            res.n_m2_confident += 1
        else:                                    # ambiguous -> 1/k over the tie
            compat[rid] = list(tie)
            if straddles:
                res.n_deferred += 1
            else:
                res.n_containment += 1

    # --- EM: certain reads anchor abundance; ambiguous reads redistribute by it. ----
    abundance, weights = _em_quantify(compat, n, em_max_iter, em_tol)
    res.weights = weights
    res.survivors = {m for m in range(n) if abundance[m] >= min_support}
    return res
