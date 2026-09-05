"""Calibrated null for the existence statistic: centred molecule-level bootstrap,
plus the sample-level hierarchical error control (spec B4 / B6).

Why not chi-square. ``T_C`` is a cross-fitted, boundary-constrained (pi >= 0)
mixture contrast over reads that share candidates, on a candidate that was
itself selected -- none of the asymptotic conditions hold. The only null we can
defend is one that re-runs the whole procedure on resampled MOLECULES.

Resampling unit. A molecule is resampled as a WHOLE row (every channel, every
candidate's weight): resampling channels, junctions, or already-computed fold
deltas separately would destroy the within-read correlation the statistic
depends on. Each replicate re-assigns folds, re-fits H1 and H0 per fold, and
re-scores held-out reads -- the same code path as the observed statistic.

The null. H0 says C generates none of the reads: within a read, C's weight is
EXCHANGEABLE with an alternative's -- the read has no reason to prefer one
label over the other. Each replicate therefore resamples whole molecules
(pairs bootstrap) AND, per read with probability 1/2, swaps C's weight with
a uniformly chosen finite alternative's weight (a within-read label
permutation; the read's BEST alternative is not an exchangeable partner), then
re-runs the full fold / fit / score path to get ``T^(b)_null``. A real
generator's systematic per-read advantage is destroyed by the swap; a
candidate whose column duplicates an alternative's is unchanged by it, scores
0 with zero spread, and is reported degenerate (unresolved), which is the
correct verdict for "adds nothing". Then

    p+ = (1 + #{T^(b)_null >= T^obs}) / (B + 1)      (adding C helps)
    p- = (1 + #{T^(b)_null <= T^obs}) / (B + 1)      (adding C hurts)

Why not centre on T_obs: a redundant candidate has NEGATIVE expected
cross-fitted T (the held-out cost of a superfluous parameter). A centred
null puts its mass at T_obs and is therefore anti-conservative in the "harm"
direction (measured: 0.15 at nominal 0.05). A duplicate-column null has zero
spread and is worse. The permutation null carries the bias and the spread by
construction. It is written into the artifact as
``label_permutation_molecule_bootstrap``; a parametric null would need a
read/CIGAR/signal generator we do not have.

Resolution. With ``(b+1)/(B+1)`` p-values the smallest attainable p is
``1/(B+1)``; the correction must be reachable, so B grows with the family:
``B = 2^ceil(log2 max(4096, B_min + 1)) - 1`` with ``B_min + 1 >= max(G*H_G,
max_g m_g) / alpha_d``. A p that cannot cross its threshold is reported as
``untestable_resolution``, never rounded.

Error control (B6). Per direction: a component omnibus max-statistic (same
replicates), Benjamini-Yekutieli across components (arbitrary dependence),
Holm within a passing component. Both layers are reported; the artifact does
not claim a flat candidate-level FDR.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

NULL_METHOD = "label_permutation_molecule_bootstrap"
DIRECTIONAL_ALPHA = 0.025          # existence alpha 0.05 split over the two directions


# --------------------------------------------------------------------------
# vectorised mixture machinery over a component weight matrix
# --------------------------------------------------------------------------
def fit_mixture(W: np.ndarray, *, max_iter: int = 200, tol: float = 1e-9) -> Tuple[np.ndarray, bool]:
    """EM for mixture weights over columns of ``W`` (log weights, -inf = absent).

    Exact for the fixed-kernel finite mixture of Part A8: the only parameters
    are the weights. Rows with no finite entry are ignored.
    """
    n, k = W.shape
    theta = np.full(k, 1.0 / k)
    if n == 0 or k == 0:
        return theta, True
    finite = np.isfinite(W)
    use = finite.any(axis=1)
    if not use.any():
        return theta, True
    Wu = W[use]
    fu = finite[use]
    rowmax = np.where(fu, Wu, -np.inf).max(axis=1, keepdims=True)
    E = np.where(fu, np.exp(Wu - rowmax), 0.0)
    for _ in range(max_iter):
        R = E * theta[None, :]
        z = R.sum(axis=1, keepdims=True)
        ok = z[:, 0] > 0
        if not ok.any():
            return theta, True
        R = R[ok] / z[ok]
        new = R.mean(axis=0)
        delta = float(np.abs(new - theta).max())
        theta = new
        if delta < tol:
            return theta, True
    return theta, False


def predictive_ll(W: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Per-row log sum_j theta_j exp(W_ij); -inf where nothing is finite."""
    finite = np.isfinite(W) & (theta[None, :] > 0)
    if W.shape[0] == 0:
        return np.zeros(0)
    rowmax = np.where(finite, W, -np.inf).max(axis=1, keepdims=True)
    safe = np.where(np.isfinite(rowmax), rowmax, 0.0)
    E = np.where(finite, np.exp(W - safe) * theta[None, :], 0.0)
    s = E.sum(axis=1)
    out = np.full(W.shape[0], -np.inf)
    ok = s > 0
    out[ok] = np.log(s[ok]) + safe[ok, 0]
    return out


def balanced_folds(n: int, max_folds: int, rng: np.random.Generator) -> Tuple[np.ndarray, int]:
    """Round-robin folds over a random order: every fold non-empty, K=min(max,n)."""
    if n < 2:
        return np.zeros(n, dtype=int), 0
    k = min(max_folds, n)
    order = rng.permutation(n)
    folds = np.empty(n, dtype=int)
    folds[order] = np.arange(n) % k
    return folds, k


def crossfit_statistics(
    W: np.ndarray, elig: np.ndarray, folds: np.ndarray, k_folds: int,
    test_cols: Sequence[int],
) -> Dict[int, float]:
    """``T_C`` for every column in ``test_cols`` under one fold assignment.

    Per fold: H1 weights fitted on the training rows over all columns; for each
    C, H0 weights fitted over the columns without C; the held-out contribution
    of a read is ``ll1 - ll0`` when the read is eligible, has C in its support,
    and has at least one alternative in its support.
    """
    n, k = W.shape
    T = {j: 0.0 for j in test_cols}
    if k_folds < 2:
        return T
    for f in range(k_folds):
        test = folds == f
        train = ~test
        if not test.any() or not train.any():
            continue
        theta1, _ = fit_mixture(W[train])
        ll1_all = predictive_ll(W[test], theta1)
        Wt = W[test]
        el = elig[test]
        for j in test_cols:
            others = [c for c in range(k) if c != j]
            theta0, _ = fit_mixture(W[train][:, others])
            ll0 = predictive_ll(Wt[:, others], theta0)
            hasC = np.isfinite(Wt[:, j])
            n_sup = np.isfinite(Wt).sum(axis=1)
            use = el & hasC & (n_sup >= 2) & np.isfinite(ll1_all) & np.isfinite(ll0)
            if use.any():
                T[j] += float((ll1_all[use] - ll0[use]).sum())
    return {j: 2.0 * v for j, v in T.items()}


# --------------------------------------------------------------------------
# bootstrap
# --------------------------------------------------------------------------
@dataclass
class BootstrapResult:
    t_obs: Dict[int, float]
    p_plus: Dict[int, float]
    p_minus: Dict[int, float]
    exceed_plus: Dict[int, int]
    exceed_minus: Dict[int, int]
    component_p_plus: float
    component_p_minus: float
    replicates: int
    degenerate: bool


def required_replicates(n_components: int, max_component_size: int,
                        alpha_d: float = DIRECTIONAL_ALPHA) -> int:
    """Smallest B making the strictest threshold attainable (spec B4.3)."""
    g = max(n_components, 1)
    h_g = sum(1.0 / i for i in range(1, g + 1))
    b_min_plus_1 = max(g * h_g / alpha_d, max(max_component_size, 1) / alpha_d)
    target = max(4096.0, b_min_plus_1)
    return int(2 ** math.ceil(math.log2(target))) - 1


def _fit_mixture_batched(E: np.ndarray, train: np.ndarray, *, max_iter: int = 60,
                         tol: float = 1e-6, init: Optional[np.ndarray] = None) -> np.ndarray:
    """EM over a batch: ``E`` (B, n, k) = exp-normalised weights (0 = absent),
    ``train`` (B, n) row mask. Returns theta (B, k). Columns absent from every
    training row get weight 0, which is how H0 (candidate removed) is expressed.
    ``init`` warm-starts (H0 from the H1 fit with the removed column zeroed)."""
    Bn, n, k = E.shape
    theta = np.full((Bn, k), 1.0 / k) if init is None else init.copy()
    tm = train[:, :, None].astype(float)
    n_train = np.maximum((train & (E.sum(axis=2) > 0)).sum(axis=1), 1)[:, None]
    for _ in range(max_iter):
        R = E * theta[:, None, :]
        z = R.sum(axis=2, keepdims=True)
        R = np.where(z > 0, R / np.where(z > 0, z, 1.0), 0.0)
        new = (R * tm).sum(axis=1) / n_train
        delta = np.abs(new - theta).max()
        theta = new
        if delta < tol:
            break
    return theta


def _predictive_batched(E: np.ndarray, rowmax: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """(B, n) log-likelihood per row: log sum_j theta_j E_ij + rowmax."""
    s = (E * theta[:, None, :]).sum(axis=2)
    out = np.full(s.shape, -np.inf)
    ok = s > 0
    out[ok] = np.log(s[ok]) + rowmax[ok]
    return out


def crossfit_batched(W: np.ndarray, elig: np.ndarray, folds: np.ndarray, k_folds: int,
                     test_cols: Sequence[int]) -> Dict[int, np.ndarray]:
    """Batched ``T_C``: ``W`` (B, n, k) log weights, ``folds`` (B, n). Returns
    {col: T (B,)}. The observed data is the batch of size 1."""
    Bn, n, k = W.shape
    finite = np.isfinite(W)
    rowmax = np.where(finite, W, -np.inf).max(axis=2)
    safe = np.where(np.isfinite(rowmax), rowmax, 0.0)
    E = np.where(finite, np.exp(W - safe[:, :, None]), 0.0)
    n_sup = finite.sum(axis=2)
    T = {j: np.zeros(Bn) for j in test_cols}
    if k_folds < 2:
        return T
    for f in range(k_folds):
        test = folds == f
        train = ~test
        theta1 = _fit_mixture_batched(E, train)
        ll1 = _predictive_batched(E, safe, theta1)
        for j in test_cols:
            E0 = E.copy(); E0[:, :, j] = 0.0
            init = theta1.copy(); init[:, j] = 0.0
            init = init / np.maximum(init.sum(axis=1, keepdims=True), 1e-12)
            theta0 = _fit_mixture_batched(E0, train, init=init)
            ll0 = _predictive_batched(E0, safe, theta0)
            use = test & elig & finite[:, :, j] & (n_sup >= 2) & np.isfinite(ll1) & np.isfinite(ll0)
            contrib = np.where(use, ll1 - ll0, 0.0)
            T[j] += contrib.sum(axis=1)
    return {j: 2.0 * v for j, v in T.items()}


def bootstrap_component(
    W: np.ndarray, elig: np.ndarray, test_cols: Sequence[int], *,
    max_folds: int, replicates: int, seed: int, component_tag: str,
    batch: int = 256,
) -> BootstrapResult:
    """Within-read label-permutation bootstrap of ``T_C`` for every testable column.

    The observed statistic and every replicate go through the SAME fold /
    fit / score path. Replicate ``b`` resamples whole rows with replacement,
    draws a fresh balanced fold assignment, and -- for each tested column C --
    swaps, per read and with probability 1/2, C's weight with a random finite
    alternative's, so ``T^(b)`` is the score of a candidate no read prefers
    systematically.
    """
    n, k = W.shape
    rng = np.random.default_rng([seed, hash(component_tag) & 0xFFFFFFFF])
    test_cols = list(test_cols)
    folds, kf = balanced_folds(n, max_folds, rng)
    t_obs_arr = crossfit_batched(W[None], elig[None], folds[None], kf, test_cols)
    t_obs = {j: float(v[0]) for j, v in t_obs_arr.items()}
    empty = BootstrapResult(t_obs=t_obs, p_plus={}, p_minus={}, exceed_plus={j: 0 for j in test_cols},
                            exceed_minus={j: 0 for j in test_cols}, component_p_plus=1.0,
                            component_p_minus=1.0, replicates=0, degenerate=True)
    if n < 2 or kf < 2 or k < 2 or replicates <= 0:
        return empty
    ex_p = {j: 0 for j in test_cols}
    ex_m = {j: 0 for j in test_cols}
    comp_obs_plus = max((max(v, 0.0) for v in t_obs.values()), default=0.0)
    comp_obs_minus = max((max(-v, 0.0) for v in t_obs.values()), default=0.0)
    comp_ex_plus = comp_ex_minus = 0
    spread = 0.0
    done = 0
    while done < replicates:
        m = min(batch, replicates - done)
        idx = rng.integers(0, n, size=(m, n))
        Wb = W[idx]                                  # (m, n, k) whole rows
        eb = elig[idx]
        fb = np.empty((m, n), dtype=int)
        for b in range(m):
            fb[b], _ = balanced_folds(n, max_folds, rng)
        zp = np.zeros(m); zm = np.zeros(m)
        for j in test_cols:
            Wn = Wb.copy()
            # a RANDOM finite alternative per (replicate, read): the exchangeable
            # partner. (The read's BEST alternative is not exchangeable with C --
            # swapping with it biases the null upward; measured p- 0.26 at 0.05.)
            others = np.asarray([c for c in range(k) if c != j])
            fin_o = np.isfinite(Wn[:, :, others])                      # (m, n, k-1)
            has_alt = fin_o.any(axis=2) & np.isfinite(Wn[:, :, j])
            u = rng.random((m, n, 1))
            cum = np.cumsum(fin_o, axis=2)
            tot = cum[:, :, -1:]
            pick = np.argmax(cum > u * tot, axis=2)                     # uniform among finite
            alt_col = others[pick]
            swap = has_alt & (rng.random((m, n)) < 0.5)
            bi, ri = np.nonzero(swap)
            ai = alt_col[bi, ri]
            cvals = Wn[bi, ri, j].copy()
            Wn[bi, ri, j] = Wn[bi, ri, ai]
            Wn[bi, ri, ai] = cvals
            t_b = crossfit_batched(Wn, eb, fb, kf, [j])[j]
            spread += float(np.abs(t_b - t_b.mean()).sum())
            ex_p[j] += int((t_b >= t_obs[j]).sum())
            ex_m[j] += int((t_b <= t_obs[j]).sum())
            zp = np.maximum(zp, np.maximum(t_b, 0.0))
            zm = np.maximum(zm, np.maximum(-t_b, 0.0))
        comp_ex_plus += int((zp >= comp_obs_plus).sum())
        comp_ex_minus += int((zm >= comp_obs_minus).sum())
        done += m
    if spread == 0.0:
        return empty                              # the null has no spread at all
    B1 = replicates + 1
    return BootstrapResult(
        t_obs=t_obs,
        p_plus={j: (1 + ex_p[j]) / B1 for j in test_cols},
        p_minus={j: (1 + ex_m[j]) / B1 for j in test_cols},
        exceed_plus=ex_p, exceed_minus=ex_m,
        component_p_plus=(1 + comp_ex_plus) / B1,
        component_p_minus=(1 + comp_ex_minus) / B1,
        replicates=replicates, degenerate=False,
    )


# --------------------------------------------------------------------------
# hierarchical error control
# --------------------------------------------------------------------------
def by_adjust(pvals: Sequence[float]) -> List[float]:
    """Benjamini-Yekutieli q-values (valid under arbitrary dependence)."""
    m = len(pvals)
    if m == 0:
        return []
    c_m = sum(1.0 / i for i in range(1, m + 1))
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [1.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        val = min(prev, pvals[i] * m * c_m / rank)
        q[i] = val
        prev = val
    return q


def holm_adjust(pvals: Sequence[float]) -> List[float]:
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    adj = [1.0] * m
    running = 0.0
    for rank, i in enumerate(order):
        val = min(1.0, (m - rank) * pvals[i])
        running = max(running, val)
        adj[i] = running
    return adj
