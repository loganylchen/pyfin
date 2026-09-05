#!/usr/bin/env python3
"""Score assembled transcripts against the Sequin CONCENTRATION truth.

SIRV4 gives us only a structural axis (gffcompare Sn/Pr). Sequin ships an
independent molar ladder (``RNAsequins_MixA_v1.tsv``, attomol/ul), so an
assembly can be judged on a second, orthogonal axis: given that you emitted a
model for a truth transcript, did you get its ABUNDANCE right?

Two numbers are always reported together, because either alone is gameable:

* ``matched_*``  — correlation over the truth transcripts the assembly actually
  recovered (class code ``=``). This is quantification quality CONDITIONAL on
  recovery. A tool that emits only the ten most abundant transcripts scores
  well here, which is why it is never reported alone.
* ``complete_*`` — the same statistics over EVERY truth transcript, with an
  unrecovered transcript entered at abundance 0. This is the joint
  discovery+quantification metric: missing a transcript is a quantification
  error of its full concentration.

Statistics mirror the benchmark's own ``quant_metrics.tsv`` so PyFIN is
comparable to the NanoCount rows already in the tree: Spearman rho on raw
values, Pearson r on log10(x+1), and NRMSE on log10(x+1) after scaling the
query to the truth's total (assemblers report TPM, the truth is attomol/ul, so
only the SHAPE of the ladder is comparable -- never the units).

Usage:
    python3 sequin_quant_score.py --runs _sequingrid/gv_off --tag pyfin_gv_off \
        --out _sequingrid/quant_score.tsv
"""
from __future__ import annotations

import argparse
import collections
import csv
import math
import re
from pathlib import Path

BENCH_REFS = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/refs")
CONCENTRATION = BENCH_REFS / "concentration" / "RNAsequins_MixA_v1.tsv"

COLS = ["sample", "cond", "tool", "status", "n_truth", "n_query", "n_matched",
        "matched_spearman", "matched_pearson_log10p1", "matched_nrmse_log10p1",
        "complete_spearman", "complete_pearson_log10p1", "complete_nrmse_log10p1"]


def load_concentration(path: Path) -> dict:
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["transcript_id"]] = float(row["concentration_attomol_per_ul"])
    return out


def _ranks(xs):
    """Average ranks (ties shared), so Spearman is correct on tied zeros."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def _pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sxx = sum((x - mx) ** 2 for x in xs)
    syy = sum((y - my) ** 2 for y in ys)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def _spearman(xs, ys):
    return _pearson(_ranks(xs), _ranks(ys))


def _nrmse_log10p1(truth, query):
    """RMSE on log10(x+1) after scaling query to the truth's total, normalised
    by the truth's log-range. Scaling first is what makes TPM comparable to a
    molar ladder at all: only the shape carries information."""
    tq = sum(query)
    tt = sum(truth)
    if tq <= 0 or tt <= 0:
        return float("nan")
    scale = tt / tq
    lt = [math.log10(v + 1.0) for v in truth]
    lq = [math.log10(v * scale + 1.0) for v in query]
    rng = max(lt) - min(lt)
    if rng <= 0:
        return float("nan")
    mse = sum((a - b) ** 2 for a, b in zip(lt, lq)) / len(lt)
    return math.sqrt(mse) / rng


def query_abundance(cell: Path) -> dict:
    """query transcript_id -> abundance (TPM preferred, else EM abundance)."""
    scores = cell / "out" / "scores.tsv"
    out = {}
    if not scores.exists():
        return out
    with open(scores) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            v = row.get("tpm") or row.get("abundance") or "0"
            try:
                out[row["candidate_id"]] = float(v)
            except ValueError:
                continue
    return out


def matched_pairs(cell: Path):
    """truth transcript_id -> query transcript_id, from gffcompare's tmap.

    Only class code ``=`` (exact intron chain) counts as recovery: this is the
    same criterion the structural score uses, so the two axes agree on what
    "found it" means. When several query models match one truth transcript,
    their abundances are SUMMED -- splitting one transcript's reads across two
    emitted models is a quantification error, not a reason to pick the best.
    """
    tmap = None
    for cand in (cell / "gc" / "cmp.assembly.gtf.tmap",
                 cell / "out" / "cmp.assembly.gtf.tmap"):
        if cand.exists():
            tmap = cand
            break
    if tmap is None:
        return None
    pairs = collections.defaultdict(list)
    with open(tmap) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("class_code") == "=":
                pairs[row["ref_id"]].append(row["qry_id"])
    return pairs


def score_cell(cell: Path, sample: str, cond: str, tool: str, conc: dict) -> dict:
    row = {"sample": sample, "cond": cond, "tool": tool, "status": "?",
           "n_truth": len(conc)}
    qa = query_abundance(cell)
    pairs = matched_pairs(cell)
    if pairs is None:
        row["status"] = "no_tmap"
        return row
    row["n_query"] = len(qa)
    m_t, m_q = [], []
    c_t, c_q = [], []
    n_matched = 0
    for tid, amount in sorted(conc.items()):
        qids = pairs.get(tid, [])
        val = sum(qa.get(q, 0.0) for q in qids)
        if qids:
            n_matched += 1
            m_t.append(amount)
            m_q.append(val)
        c_t.append(amount)
        c_q.append(val)
    row["n_matched"] = n_matched
    if n_matched >= 2:
        row["matched_spearman"] = _spearman(m_t, m_q)
        row["matched_pearson_log10p1"] = _pearson(
            [math.log10(v + 1) for v in m_t], [math.log10(v + 1) for v in m_q])
        row["matched_nrmse_log10p1"] = _nrmse_log10p1(m_t, m_q)
    row["complete_spearman"] = _spearman(c_t, c_q)
    row["complete_pearson_log10p1"] = _pearson(
        [math.log10(v + 1) for v in c_t], [math.log10(v + 1) for v in c_q])
    row["complete_nrmse_log10p1"] = _nrmse_log10p1(c_t, c_q)
    row["status"] = "ok"
    return row


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True,
                    help="run root holding <sample>/<cond>/ cells")
    ap.add_argument("--tag", required=True, help="tool label for the rows")
    ap.add_argument("--out", required=True)
    ap.add_argument("--concentration", default=str(CONCENTRATION))
    a = ap.parse_args()

    conc = load_concentration(Path(a.concentration))
    root = Path(a.runs)
    rows = []
    for cell in sorted(root.glob("*/*")):
        if not (cell / "out").is_dir():
            continue
        rows.append(score_cell(cell, cell.parent.name, cell.name, a.tag, conc))
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["sample"], r["cond"])):
            w.writerow(r)
    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"{ok}/{len(rows)} cells scored -> {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
