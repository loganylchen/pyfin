#!/usr/bin/env python3
"""Quantification accuracy against the Sequin spike-in's KNOWN concentrations.

Every benchmark this project runs so far scores STRUCTURE (gffcompare `=`).
The Sequin mix has real quantitative ground truth -- attomol/uL per transcript
-- so abundance can be scored directly instead of against a pseudo-truth.

For each (sample, condition, tool) it reports, over the transcripts the tool
emitted and gffcompare matched `=` to a truth transcript:

  Spearman  rank correlation with the known concentration (scale-free)
  Pearson   on log10, after normalizing both sides to parts-per-million
  MARD      median absolute relative difference, |est - true| / ((est+true)/2)
            on the ppm scale (the standard quantification metric)

and the same three over every EXPRESSED truth transcript (the project's
expression contract: NanoCount TPM > 0 in that library), with a missed
transcript scored as zero -- that version mixes detection and quantification,
which is what a user actually experiences, so both are given. The denominator
is the expressed set, not the whole mix: most of Mix A sits below this
sequencing depth, and scoring undetectable transcripts as misses saturates
MARD at its ceiling of 2.0 for every tool including the oracle.

Abundance sources: pyfin `scores.tsv` (tpm), stringtie3's GTF `TPM` attribute,
NanoCount's `tpm` (a REFERENCE-based quantifier given the true annotation --
the oracle baseline, not a competitor).
"""
import argparse, collections, csv, math, re, statistics
from pathlib import Path

BENCH = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results"
             "/sequin_v1_full_sweep_ratio_grid_v2")
CONC = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/refs"
            "/concentration/RNAsequins_MixA_v1.tsv")
WORK = Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_sequingrid")


def spearman(xs, ys):
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(ranks(xs), ranks(ys))


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0 or sy == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy)


def ppm(d):
    tot = sum(d.values())
    return {k: 1e6 * v / tot for k, v in d.items()} if tot > 0 else dict(d)


def mard(est, true):
    """Median |e-t| / ((e+t)/2) over pairs where at least one side is > 0."""
    vals = [abs(e - t) / ((e + t) / 2) for e, t in zip(est, true) if e + t > 0]
    return statistics.median(vals) if vals else float("nan")


def load_conc():
    out = {}
    with open(CONC) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["transcript_id"]] = float(r["concentration_attomol_per_ul"])
    return out


def tmap_matches(path):
    """ref_id -> every query id matched `=` to it.

    More than one emitted model can carry a reference transcript's intron chain
    (different ends), and those models SPLIT that transcript's reads, so their
    abundances are summed rather than one being picked -- picking one would
    under-report exactly where a tool fragments a transcript.
    """
    out = collections.defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["class_code"] == "=":
                out[r["ref_id"]].append(r["qry_id"])
    return dict(out)


def pyfin_abundance(cell_out, column="tpm"):
    """``tpm`` is our reported value: abundance / length_kb, renormalised.

    ``abundance`` is the raw probability-weighted READ COUNT. On dRNA each
    molecule yields at most one read, so the count is already proportional to
    MOLAR abundance -- which is what the Sequin truth is (attomol/uL) -- and
    dividing by length turns a molar quantity into a mass-like one. Both are
    scored so the difference is measured rather than assumed.
    """
    d = {}
    with open(Path(cell_out) / "scores.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[r["candidate_id"]] = float(r[column])
    return d


def gtf_tpm(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            t = re.search(r'transcript_id "([^"]+)"', f[8])
            v = re.search(r'TPM "([^"]+)"', f[8])
            if t and v:
                d[t.group(1)] = float(v.group(1))
    return d


def nanocount_tpm(path):
    d = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d[r["transcript_name"]] = float(r["tpm"])
    return d


def evaluate(name, est_by_qry, ref2qry, conc, expressed, common=None):
    """est_by_qry: query id -> estimate. Returns the metric row.

    ``common``: when given, the conditional columns are computed over THIS set
    instead of the tool's own detections. Comparing per-tool sets is comparing
    different problems -- an easy set and a hard set renormalised separately --
    so a cross-tool claim needs the common set; the own-set columns describe
    each tool internally and are kept for that.
    """
    est_ref = {}
    for ref, qrys in ref2qry.items():
        if ref not in conc:
            continue
        vals = [est_by_qry[q] for q in (qrys if isinstance(qrys, list) else [qrys])
                if q in est_by_qry]
        if vals:
            est_ref[ref] = sum(vals)
    if len(est_ref) < 3:
        return None
    # unmatched mass: abundance this tool put on models that matched no truth
    # transcript. Dropping it silently would let a tool assign mass to wrong
    # structures for free, so it is reported (as a fraction of its total).
    matched_q = {q for qs in ref2qry.values() for q in (qs if isinstance(qs, list) else [qs])}
    tot = sum(est_by_qry.values())
    unmatched = (sum(v for q, v in est_by_qry.items() if q not in matched_q) / tot
                 if tot > 0 else float("nan"))
    # conditional: transcripts we detected (or the common set, when given)
    keys = sorted(k for k in (common if common is not None else est_ref)
                  if k in conc and (common is None or k in est_ref or True))
    if len(keys) < 3:
        return None
    e = ppm({k: est_ref.get(k, 0.0) for k in keys})
    t = ppm({k: conc[k] for k in keys})
    ex, tx = [e[k] for k in keys], [t[k] for k in keys]
    sp = spearman(ex, tx)
    pe = pearson([math.log10(v + 1) for v in ex], [math.log10(v + 1) for v in tx])
    md = mard(ex, tx)
    # unconditional: every EXPRESSED transcript, missed ones = 0
    all_keys = sorted(k for k in conc if k in expressed)
    if len(all_keys) < 3:
        all_keys = keys
    e2 = ppm({k: est_ref.get(k, 0.0) for k in all_keys})
    t2 = ppm({k: conc[k] for k in all_keys})
    ex2, tx2 = [e2[k] for k in all_keys], [t2[k] for k in all_keys]
    return dict(tool=name, n_detected=len(est_ref), n_scored=len(keys),
                n_truth=len(all_keys), unmatched_mass_frac=round(unmatched, 4),
                spearman=sp, pearson_log=pe, mard=md,
                spearman_all=spearman(ex2, tx2),
                pearson_log_all=pearson([math.log10(v + 1) for v in ex2],
                                        [math.log10(v + 1) for v in tx2]),
                mard_all=mard(ex2, tx2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="+", required=True,
                    help="LABEL=<dir> pyfin run roots (dir/<sample>/<cond>/out)")
    ap.add_argument("--conds", default="full,p50,p00,c_jitter50_10bp")
    ap.add_argument("--samples", default="")
    ap.add_argument("--contract-cells", default=str(WORK / "expression_contract" / "cells"),
                    help="where the expression-contract gffcompare cells live")
    ap.add_argument("--out", default="")
    ap.add_argument("--column", default="tpm", choices=["tpm", "abundance", "num_reads"],
                    help="which pyfin value to score: tpm (count/length_kb), "
                         "abundance (soft read count) or num_reads (hard count)")
    ap.add_argument("--common", action="store_true",
                    help="print the COMMON-set conditional columns (both are written to --out)")
    a = ap.parse_args()
    conc = load_conc()
    conds = [c for c in a.conds.split(",") if c]
    samples = ([s for s in a.samples.split(",") if s]
               or sorted(p.name for p in BENCH.glob("SGNex_*") if p.is_dir()))
    rows = []
    pending = collections.defaultdict(list)   # (sample, cond) -> [(name, est_by_qry, ref2qry)]
    for sample in samples:
        # the expression contract's expressed set for this library
        ncp = BENCH / "_ref" / sample / "expression" / "nanocount.tsv"
        expressed = set()
        if ncp.exists():
            with open(ncp) as fh:
                expressed = {r["transcript_name"] for r in csv.DictReader(fh, delimiter="\t")
                             if float(r["tpm"]) > 0.0}
        for cond in conds:
            # pyfin arms: use the expression-contract tmap for that label
            for spec in a.arms:
                label, _, root = spec.partition("=")
                cell = Path(root) / sample / cond / "out"
                tmap = Path(a.contract_cells) / sample / cond / label / "cmp.query.gtf.tmap"
                if not (cell / "scores.tsv").exists() or not tmap.exists():
                    continue
                pending[(sample, cond, frozenset(expressed))].append(
                    (label, pyfin_abundance(cell, a.column), tmap_matches(tmap)))
            asm = BENCH / sample / cond / "assembly"
            st3 = asm / "stringtie3.gtf"
            if st3.exists():
                tm = asm / "cmp.stringtie3.gtf.tmap"
                if tm.exists():
                    pending[(sample, cond, frozenset(expressed))].append(
                        ("stringtie3", gtf_tpm(st3), tmap_matches(tm)))
            nc = asm / "nanocount.tsv"
            if nc.exists():
                d = nanocount_tpm(nc)
                pending[(sample, cond, frozenset(expressed))].append(
                    ("nanocount(oracle)", d, {k: [k] for k in d}))
    # the common set for a cell = truth transcripts EVERY non-oracle tool
    # matched, so the conditional columns compare the same problem
    for (sample, cond, expressed), items in sorted(pending.items(), key=lambda x: x[0][:2]):
        sets = []
        for name, est, ref2q in items:
            if name == "nanocount(oracle)":
                continue
            sets.append({ref for ref, qs in ref2q.items()
                         if ref in conc and any(q in est for q in
                                                (qs if isinstance(qs, list) else [qs]))})
        common = set.intersection(*sets) if sets else set()
        for name, est, ref2q in items:
            r = evaluate(name, est, ref2q, conc, set(expressed))
            if r:
                rows.append(dict(sample=sample, cond=cond, common=0, **r))
            rc = evaluate(name, est, ref2q, conc, set(expressed), common=common)
            if rc:
                rows.append(dict(sample=sample, cond=cond, common=1, **rc))
    if not rows:
        raise SystemExit("no cells scored")
    if a.out:
        with open(a.out, "w") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)
    print(f"{'cond':18s} {'tool':20s} {'n_det':>6s} {'rho':>6s} {'r_log':>6s} {'MARD':>6s}   "
          f"{'rho_all':>7s} {'r_log_all':>9s} {'MARD_all':>8s}")
    for cond in conds:
        for tool in sorted({r["tool"] for r in rows}):
            g = [r for r in rows if r["cond"] == cond and r["tool"] == tool
                 and r["common"] == (1 if a.common else 0)]
            if not g:
                continue
            m = lambda k: statistics.mean(r[k] for r in g if r[k] == r[k])
            print(f"{cond:18s} {tool:20s} {statistics.mean(r['n_detected'] for r in g):6.0f} "
                  f"{m('spearman'):6.3f} {m('pearson_log'):6.3f} {m('mard'):6.3f}   "
                  f"{m('spearman_all'):7.3f} {m('pearson_log_all'):9.3f} {m('mard_all'):8.3f}")
    print(f"\n({len(samples)} libraries; rho/r/MARD over the transcripts each tool DETECTED, "
          f"*_all over every EXPRESSED mix transcript with misses scored 0)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
