#!/usr/bin/env python3
"""Reproducible T3/T1 gap decomposition backing PERFORMANCE_STRATEGY.md §3.

Compares a PyFIN run against a same-contract competitor rescore at the
gffcompare-tmap level. Counting semantics, stated once and used everywhere:

* "query row"  = one line of `gc.assembly.gtf.tmap` (one emitted transcript).
  Several rows may match the same reference transcript.
* "truth ID"   = version-stripped reference transcript ID
  (`profile_sweep.normalized_ref_id`), the unit of recall.
* Expressed truth sets come from the sample's NanoCount table at est_count
  thresholds 1.0 (T1) and 3.0 (T3) — identical to the production scorer.
* Official honest F1 always comes from `gc.stats` + `profile_sweep`; nothing
  here re-derives it. This script only decomposes *who matches what* and what
  the non-exact budget is made of.

Usage:
    python3 t3_gap_decomposition.py \
        --pyfin-dir <run dir with gc.assembly.gtf.tmap, scores.tsv, assembly.gtf> \
        --competitor-tmap <competitor gc.assembly.gtf.tmap> \
        --sample gencode_p00
"""
from __future__ import annotations

import argparse
import csv
import re
import statistics as st  # noqa: F401  (used for per-stratum medians)
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import profile_sweep as ps  # noqa: E402


def load_tmap(path: Path):
    """Yield (ref_id_normalized_or_None, class_code, num_exons, qry_id)."""
    rows = []
    with open(path) as handle:
        handle.readline()
        for line in handle:
            x = line.rstrip("\n").split("\t")
            ref = ps.normalized_ref_id(x[1]) if x[1] and x[1] != "-" else None
            rows.append((ref, x[2], int(x[5]) if x[5].isdigit() else 0, x[4]))
    return rows


def intron_chains(gtf: Path):
    exons = defaultdict(list)
    meta = {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[2] != "exon":
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        exons[tid].append((int(f[3]) - 1, int(f[4])))
        meta[tid] = (f[0], f[6])
    chains = {}
    for tid, ex in exons.items():
        ex = sorted(ex)
        chains[tid] = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
    return chains, meta


def is_strict_subchain(a, b):
    if len(a) == 0 or len(a) >= len(b):
        return False
    return any(b[i:i + len(a)] == a for i in range(len(b) - len(a) + 1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pyfin-dir", required=True, type=Path)
    ap.add_argument("--competitor-tmap", required=True, type=Path)
    ap.add_argument("--sample", default="gencode_p00")
    args = ap.parse_args()

    ds = ps.real_dataset(args.sample)
    t1 = ps.expressed_truth(ds.nanocount, threshold=1.0)
    t3 = ps.expressed_truth(ds.nanocount, threshold=3.0)
    print(f"truth IDs: T1={len(t1)} T3={len(t3)}")

    pf = load_tmap(args.pyfin_dir / "gc.assembly.gtf.tmap")
    cp = load_tmap(args.competitor_tmap)

    def matched_ids(rows):
        return {r for r, c, _, _ in rows if c == "=" and r}

    pm, cm = matched_ids(pf), matched_ids(cp)
    for label, rows, m in (("pyfin", pf, pm), ("competitor", cp, cm)):
        eq = sum(1 for _, c, _, _ in rows if c == "=")
        print(f"{label}: rows={len(rows)} exact_rows={eq} "
              f"uniq_matched_ids={len(m)} T3_matched_ids={len(m & t3)} "
              f"T1_matched_ids={len(m & t1)}")

    p3, c3 = pm & t3, cm & t3
    print(f"T3 truth-ID overlap: both={len(p3 & c3)} "
          f"pyfin_only={len(p3 - c3)} competitor_only={len(c3 - p3)}")

    near = Counter()
    by_ref = defaultdict(list)
    for r, c, _, _ in pf:
        if r:
            by_ref[r].append(c)
    order = "=ckmnjeosxiypru"
    for ref in c3 - p3:
        codes = by_ref.get(ref)
        near["absent" if not codes else "pyfin_has_" + min(
            codes, key=lambda c: order.index(c) if c in order else 99)] += 1
    print(f"competitor-only truths, pyfin context: {dict(near.most_common())}")

    # Query-row strata for the PyFIN run (row-level, NOT truth-ID-level),
    # with per-stratum abundance/read medians from scores.tsv.
    scores_by_id = {r["candidate_id"]: r for r in csv.DictReader(
        open(args.pyfin_dir / "scores.tsv"), delimiter="\t")}
    strata = Counter()
    uniq = defaultdict(set)
    stratum_rows = defaultdict(list)
    for r, c, nex, q in pf:
        if c == "=":
            key = ("exact_T3" if r in t3 else
                   "exact_T1_only" if r in t1 else "exact_below_T1")
            strata[key] += 1
            uniq[key].add(r)
        else:
            key = f"nonexact_{c}_{'mono' if nex <= 1 else 'multi'}"
            strata[key] += 1
        if q in scores_by_id:
            stratum_rows[key].append(scores_by_id[q])
    print("\npyfin query-row strata (rows / unique IDs / med_abundance / med_reads):")
    for k in sorted(strata):
        u = len(uniq[k]) if k in uniq else ""
        g = stratum_rows.get(k, [])
        med_ab = st.median(float(x["abundance"]) for x in g) if g else float("nan")
        med_rd = st.median(int(x["num_reads"]) for x in g) if g else float("nan")
        print(f"  {k:26} {strata[k]:>6} {str(u):>5} {med_ab:>8.1f} {med_rd:>6}")

    # Observable sibling containment among pyfin outputs.
    chains, meta = intron_chains(args.pyfin_dir / "assembly.gtf")
    by_key = defaultdict(list)
    for tid in chains:
        by_key[meta[tid]].append(tid)
    contained = set()
    for _, tids in by_key.items():
        multi = [(t, chains[t]) for t in tids if chains[t]]
        for t, ch in multi:
            if any(t != u and is_strict_subchain(ch, ch2) for u, ch2 in multi):
                contained.add(t)
    qcls = {q: (c, r) for r, c, _, q in pf}
    cross = Counter()
    for tid in chains:
        c, r = qcls.get(tid, ("?", None))
        cat = ("exact_T3" if c == "=" and r in t3 else
               "exact_subT3" if c == "=" else f"nonexact_{c}")
        obs = ("mono" if not chains[tid] else
               "sub_of_sibling" if tid in contained else "standalone")
        cross[(cat, obs)] += 1
    print("\nsibling-containment cross-tab (category, observability): ")
    for (cat, obs), n in sorted(cross.items()):
        print(f"  {cat:20} {obs:14} {n}")

    # Abundance-threshold separability (the negative result).
    ranked = sorted(
        (float(scores_by_id[q]["abundance"]), qcls.get(q, ("?", None))[0] == "=")
        for q in scores_by_id if q in qcls)
    for cut in (500, 1000, 1500, 2000):
        lost = Counter(t for _, t in ranked[:cut])
        print(f"cut lowest {cut:>4} by abundance: "
              f"removes nonexact={lost[False]:>4} exact_rows={lost[True]:>4}")


if __name__ == "__main__":
    main()
