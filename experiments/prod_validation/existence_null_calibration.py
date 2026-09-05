#!/usr/bin/env python3
"""Monte-Carlo calibration of the existence bootstrap (spec D3 audit gate).

Under H0 the candidate adds nothing: its column is exchangeable with an
alternative's. The one-sided p+ must then be (super-)uniform: the fraction of
datasets with p <= a must not exceed a beyond Monte-Carlo error. Under H1 the
candidate fits every read better by `effect` nats; power must rise with n and
effect. n=1 and n=2 components must be degenerate / low-power, never
"significant by construction".
"""
import argparse, time
import numpy as np
from fin.analysis.existence_bootstrap import bootstrap_component


def h0_dataset(rng, n, k, sparsity, noise):
    """Exchangeable H0: every column, INCLUDING the tested one, is an independent
    draw of the same law (same noise, same sparsity). A candidate that is
    always present while the others are sparse is not H0 -- it explains reads
    nobody else can, which is H1."""
    W = rng.normal(0, noise, (n, k))
    W[rng.random((n, k)) < sparsity] = -np.inf
    # guarantee at least one finite column per row so rows are not empty
    empty = ~np.isfinite(W).any(axis=1)
    W[empty, rng.integers(0, k, size=int(empty.sum()))] = 0.0
    return W


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", type=int, default=400)
    ap.add_argument("--B", type=int, default=511)
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--k", type=int, default=3)
    a = ap.parse_args()
    rng = np.random.default_rng(20260905)
    t0 = time.time()
    for label, sparsity, noise in (("dense_lownoise", 0.0, 0.3), ("sparse_noisy", 0.5, 1.0)):
        ps, pm = [], []
        for s in range(a.sets):
            W = h0_dataset(rng, a.n, a.k, sparsity, noise)
            r = bootstrap_component(W, np.ones(a.n, bool), [0], max_folds=5,
                                    replicates=a.B, seed=s, component_tag=label)
            if r.degenerate:
                continue
            ps.append(r.p_plus[0]); pm.append(r.p_minus[0])
        ps, pm = np.array(ps), np.array(pm)
        print(f"H0 {label}: n_sets={len(ps)} B={a.B}")
        for lvl in (0.01, 0.025, 0.05, 0.10, 0.25, 0.5):
            se = (lvl * (1 - lvl) / len(ps)) ** 0.5
            print(f"   P(p+<={lvl:5.3f}) = {(ps<=lvl).mean():.3f}   P(p-<={lvl:5.3f}) = {(pm<=lvl).mean():.3f}   (nominal {lvl}, MC se {se:.3f})")
    print(f"[{time.time()-t0:.0f}s]")
    print("H1 power (p+ <= 0.025), 100 sets each:")
    for n in (6, 12, 30, 80):
        for eff in (0.5, 1.0, 2.0):
            hits = 0; deg = 0
            for s in range(100):
                W = np.zeros((n, 2)); W[:, 1] = -eff + rng.normal(0, 0.3, n)
                r = bootstrap_component(W, np.ones(n, bool), [0], max_folds=5,
                                        replicates=a.B, seed=s, component_tag="h1")
                if r.degenerate: deg += 1; continue
                hits += r.p_plus[0] <= 0.025
            print(f"   n={n:3d} effect={eff:.1f}: power {hits/100:.2f}  degenerate {deg}")
    for n in (1, 2, 3):
        W = np.zeros((n, 2)); W[:, 1] = -5.0
        r = bootstrap_component(W, np.ones(n, bool), [0], max_folds=5, replicates=a.B, seed=0, component_tag="tiny")
        print(f"   n={n}: degenerate={r.degenerate} p+={r.p_plus.get(0)}")


if __name__ == "__main__":
    main()
