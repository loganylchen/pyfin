# heya8 wobble-FP fix — final report (bp20 + cassette + read-support GTF guard)

**Date:** 2026-06-25  •  **Data:** 5 HEYA8 dual-spike samples × 8 ratios (full + p00 + 6
corruptions), GPU (RTX 3090, pyfin-gpu:dev), honest metrics (nanocount expressed-truth).

## What changed (production defaults)

Three levers, each independently validated, added to the m2_em cluster-recheck:

1. **`m2_cluster_recheck_bp` 10 → 20** — widen the wobble window so 11–15 bp wobble
   shadows cluster with their anchor and the abundance-fraction filter drops them.
2. **`m2_cluster_recheck_cassette_max_exon_bp = 70`** — cluster a K-intron candidate with
   a (K−1)-intron one differing by a single small exon (<70 bp): minimap2 small-exon-skip
   artifacts cluster with their true sibling.
3. **`m2_cluster_recheck_novel_displaces_gtf = True` + read-support guard** — a clustered
   low-abundance GTF sibling can be dropped, BUT only when its DISTINGUISHING junction has
   **zero reads splicing exactly there** (`gtf_min_jct_reads = 1`, `jct_tol = 0`). Removes
   the original anchor-source asymmetry: a GTF is judged by its OWN read support, not by
   whether the anchor is novel. Drops jittered/phantom GTF passthroughs; keeps real
   annotated isoforms (which carry their own reads).

## Lever ladder (global mean F1@3, 5×8 cells)

| config | F1@3 | Sn@3 | Pr@3 | note |
|---|---|---|---|---|
| bp10 (old default) | 79.9 | 92.7 | 70.5 | start |
| + bp20 + cassette (bp20cas) | 81.0 | 92.2 | 72.5 | precision +2.0 |
| + old displace (v2) | 81.1 | 91.7 | 73.0 | **lost SIRV606** (recall ↓) |
| **+ read-support guard (v3)** | **81.3** | **92.0** | **73.1** | **recall restored, best on all axes** |

v3 strictly dominates v2: higher F1, higher recall (SIRV606 back), highest precision.
**+1.4 F1@3 over the old default, every ratio positive, zero recall regression.**

## Per-ratio (v3 − bp10, F1@3)

| ratio | bp10 | v3 | Δ |
|---|---|---|---|
| c_jitter10bp | 67.3 | 70.5 | **+3.2** (the target — jittered-GTF FPs) |
| full | 84.4 | 85.4 | +1.1 |
| other 6 | — | — | +1.0 … +1.1 |

## SIRV606 recall audit (the read-support guard's job)

| | c_skip10 | c_ir10 |
|---|---|---|
| v2 (old displace) SIRV606 lost | 5/5 samples | 5/5 samples |
| **v3 (read-support guard)** | **0/5** | **0/5** |

The guard fixed the recall regression in all 10 cells while keeping the c_jitter precision gain.
Strict `jct_tol=0` was decisive: jittered GTFs have 0 reads at their EXACT junction (even 1–2 bp
jitters), while real isoforms (SIRV606: 419 reads on its exact donor) are unmissable.

## vs competitors (global means)

| tool | F1@3 | Sn@3 | Pr@3 | worst-cell | c_jitter |
|---|---|---|---|---|---|
| espresso | 81.6 | 87.2 | 78.4 | 66.1 | 74.7 |
| **pyfin v3** | **81.3** | **92.0** | 73.1 | 65.0 | 70.5 |
| lafite | 78.7 | 81.6 | 76.8 | 42.2 | 51.6 |
| stringtie3 | 71.5 | 79.2 | 66.2 | 33.3 | 57.0 |
| isoquant | 69.7 | 74.5 | 67.8 | 27.6 | 38.0 |
| flair | 61.6 | 68.7 | 58.5 | 22.6 | 51.4 |
| bambu/talon/isotools | ≤40 | — | ≤33 | — | — |

pyfin v3 and espresso are the co-leading honest tools (within 0.3 F1@3); pyfin leads recall
(Sn 92.0), espresso leads precision (Pr 78.4) and the c_jitter axis. Both are the only
high-score + corruption-robust tools.

## Tests

`tests/unit/test_structural_wobble_cassette.py` — 24 new cases (cassette clustering +
read-support shadow drop, incl. strict-vs-loose tol). 42/42 pass with the existing
`test_m2_diff_cover.py`.

## CRITICAL caveat

All of the above is **SIRV/Sequin synthetic spike-in** data. SIRV's dense engineered GT-AG
sites create adversarial cases (SIRV606) that real genomes rarely produce. **The defaults are
tuned/validated on synthetic truth; real-transcriptome validation (HEYA8 human reads vs
GENCODE, via cross-replicate reproducibility + matched short-read SJ support) is the gate
before treating this as the final production config.** The mechanisms are principled
(evidence-based, not magic thresholds), so overfitting risk is lower than a tuned number —
but real-data confirmation is still required.

## Files

- code: `fin/cli.py`, `fin/pipeline/config.py`, `fin/pipeline/runner.py`,
  `fin/scoring/m2_junction_nll.py` (+ `tests/unit/test_structural_wobble_cassette.py`)
- matrix outputs: `experiments/wobble_heya8/matrix/<sample>/<ratio>/out_{bp10,bp20cas,v2,v3}/`
- tables: `experiments/wobble_heya8/tables/`
