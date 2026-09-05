#!/usr/bin/env python3
"""Per-channel discriminative power of the existence audit, by corruption type.

For every annotated candidate in the audit: is the SUPPLIED chain the reference
chain (true) or a corrupted one? Reports, per condition, the median structural
and 3'/polyA channel contributions and how often each is positive, plus the
rank-AUC of true-vs-corrupted for the full statistic and for each channel.
Evaluation-only."""
import argparse, collections, csv, glob, re, statistics as st


def chains(gtf):
    ex = collections.defaultdict(list)
    for l in open(gtf):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if m: ex[m.group(1)].append((f[0], int(f[3]), int(f[4])))
    out = {}
    for t, e in ex.items():
        e.sort(key=lambda x: x[1])
        out[t] = (e[0][0], e[0][1], e[-1][2], tuple((a[2], b[1] - 1) for a, b in zip(e, e[1:])))
    return out


def auc(pos, neg):
    if not pos or not neg: return float("nan")
    w = sum((1.0 if x > y else 0.5 if x == y else 0.0) for x in pos for y in neg)
    return w / (len(pos) * len(neg))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True); ap.add_argument("--ref-root", required=True)
    ap.add_argument("--truth", required=True)
    a = ap.parse_args()
    truth = chains(a.truth)
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for f in sorted(glob.glob(f"{a.runs}/*/*/out/candidate_existence.tsv")):
        parts = f.split("/"); sample, cond = parts[-4], parts[-3]
        try: anno = chains(f"{a.ref_root}/{sample}/{cond}/annotation.gtf")
        except OSError: continue
        for r in csv.DictReader(open(f), delimiter="\t"):
            if r["identifiable"] != "true": continue
            cid = r["candidate_id"]; sup, tru = anno.get(cid), truth.get(cid)
            if sup is None or tru is None: continue
            # corruption type: chain differs, or same chain but endpoints differ (trunc)
            if sup[3] != tru[3]: kind = "corrupted_chain"
            elif (sup[1], sup[2]) != (tru[1], tru[2]): kind = "corrupted_ends"
            else: kind = "true"
            full = float(r["delta_predictive_ll"])
            struct = float(r["structural_llr"]) if r["structural_llr"] else 0.0
            tes = float(r["tes_polya_llr"]) if r["tes_polya_llr"] else 0.0
            rt = int(r.get("n_tes_readthrough") or 0)
            by[cond][kind].append((full, struct, tes, rt))
    print(f"{'cond':18s} {'kind':16s} {'n':>4s} {'full med':>9s} {'struct':>8s} {'tes':>8s} {'tes>0':>6s} {'readthru':>9s}")
    for cond in sorted(by):
        for kind in ("true", "corrupted_chain", "corrupted_ends"):
            v = by[cond].get(kind)
            if not v: continue
            print(f"{cond:18s} {kind:16s} {len(v):4d} {st.median(x[0] for x in v):9.2f} "
                  f"{st.median(x[1] for x in v):8.2f} {st.median(x[2] for x in v):8.2f} "
                  f"{sum(1 for x in v if x[2] > 1e-9):6d} {sum(1 for x in v if x[3] > 0):9d}")
        t = by[cond].get("true", []); c = by[cond].get("corrupted_chain", []) + by[cond].get("corrupted_ends", [])
        if t and c:
            print(f"{'':18s} AUC true>corrupted  full {auc([x[0] for x in t],[x[0] for x in c]):.3f}"
                  f"  structural {auc([x[1] for x in t],[x[1] for x in c]):.3f}"
                  f"  tes {auc([x[2] for x in t],[x[2] for x in c]):.3f}")


if __name__ == "__main__":
    main()
