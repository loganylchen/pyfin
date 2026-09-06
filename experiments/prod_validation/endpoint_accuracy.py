#!/usr/bin/env python3
"""How close are the EMITTED transcript ends to the truth's ends?

gffcompare's transcript-level match is a statement about intron chains: two
models with the same chain and 500 bp different 3' ends are both `=`. So the
structural F1 used everywhere else in this benchmark is BLIND to the extent,
which is exactly what the endpoint channels re-estimate. This scores that axis
directly: for every query that gffcompare matched to a truth transcript (`=`),
the strand-aware distance from the query's TSS/TES to that transcript's.
"""
import argparse, collections, csv, re, statistics as st
from pathlib import Path


def ends_by_tx(gtf):
    """transcript_id -> (chrom, strand, tss, tes) from transcript lines."""
    out = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m:
                continue
            s, e, strand = int(f[3]) - 1, int(f[4]), f[6]
            out[m.group(1)] = (f[0], strand, s if strand == "+" else e,
                               e if strand == "+" else s)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True,
                    help="dir of expression_contract/cells/<sample>/<cond>/<tool>")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--arms", required=True, help="comma list of tool labels")
    ap.add_argument("--runs", required=True,
                    help="comma list LABEL=run_root, aligned with --arms")
    a = ap.parse_args()
    truth = ends_by_tx(a.truth)
    roots = dict(s.split("=", 1) for s in a.runs.split(","))
    rows = collections.defaultdict(lambda: collections.defaultdict(list))
    for arm in a.arms.split(","):
        root = Path(roots[arm])
        for tmap in sorted(Path(a.cells).glob(f"*/*/{arm}/cmp.query.gtf.tmap")):
            cond = tmap.parents[1].name
            sample = tmap.parents[2].name
            q = ends_by_tx(root / sample / cond / "out" / "assembly.gtf")
            for r in csv.DictReader(open(tmap), delimiter="\t"):
                if r["class_code"] != "=" or r["ref_id"] not in truth:
                    continue
                qe, te = q.get(r["qry_id"]), truth[r["ref_id"]]
                if qe is None or qe[1] != te[1]:
                    continue
                sign = 1 if te[1] == "+" else -1
                rows[cond][arm].append((abs(qe[2] - te[2]), abs(qe[3] - te[3]),
                                        sign * (qe[3] - te[3])))
    arms = a.arms.split(",")
    print(f"{'cond':18s} {'arm':8s} {'n':>6s} {'medTSS':>7s} {'medTES':>7s} "
          f"{'TES<=20':>8s} {'TSS<=20':>8s}")
    for cond in sorted(rows):
        for arm in arms:
            v = rows[cond].get(arm)
            if not v:
                continue
            tss = [x[0] for x in v]; tes = [x[1] for x in v]
            print(f"{cond:18s} {arm:8s} {len(v):6d} {st.median(tss):7.0f} "
                  f"{st.median(tes):7.0f} "
                  f"{100*sum(1 for x in tes if x <= 20)/len(tes):7.1f}% "
                  f"{100*sum(1 for x in tss if x <= 20)/len(tss):7.1f}%")
    print()
    for arm in arms:
        v = [x for c in rows for x in rows[c].get(arm, [])]
        if not v:
            continue
        tss = [x[0] for x in v]; tes = [x[1] for x in v]
        print(f"ALL {arm:8s} n={len(v):6d} medTSS={st.median(tss):5.0f} "
              f"medTES={st.median(tes):5.0f} "
              f"TES<=20 {100*sum(1 for x in tes if x <= 20)/len(tes):.1f}% "
              f"TSS<=20 {100*sum(1 for x in tss if x <= 20)/len(tss):.1f}%")


if __name__ == "__main__":
    main()
