#!/usr/bin/env python3
"""Dissect the pyfin-vs-stringtie3 `=` gap on one human cell.

For the three sets of expressed truth transcripts -- matched by stringtie3
only, by pyfin only, by both -- report what read evidence exists in the BAM
(exact full-chain reads, 3'-suffix reads), what NanoCount assigned, what
stringtie3 reported as coverage, and which class code the OTHER tool gave at
that reference. The point is to say where the 554 stringtie3-only matches
would have to come from in our pipeline, before proposing any change.
"""
import argparse, collections, csv, re, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from human_recall_decomposition import chains_from_gtf, read_chains


def eq_set(tracking):
    out = set()
    with open(tracking) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[3] == "=":
                out.add(f[2].split("|")[1])
    return out


def ref_codes(tracking):
    d = collections.defaultdict(set)
    with open(tracking) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[2] != "-":
                d[f[2].split("|")[1]].add(f[3])
    return d


def st3_cov(tracking):
    """ref tx -> max stringtie3 cov among its '=' query models."""
    out = {}
    with open(tracking) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[3] != "=":
                continue
            t = f[2].split("|")[1]
            for q in f[4:]:
                parts = q.split("|")
                if len(parts) >= 5:
                    try:
                        out[t] = max(out.get(t, 0.0), float(parts[3]))
                    except ValueError:
                        pass
    return out


def q(v, p):
    v = sorted(v)
    return v[int(p * (len(v) - 1))] if v else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--pyfin-tracking", required=True)
    ap.add_argument("--st3-tracking", required=True)
    ap.add_argument("--nanocount", required=True)
    ap.add_argument("--out", required=True, help="per-transcript TSV")
    ap.add_argument("--chroms", default="chr1,chr19,chr21")
    a = ap.parse_args()
    chroms = set(a.chroms.split(","))
    truth = chains_from_gtf(a.truth, chroms=chroms)
    anno = chains_from_gtf(a.annotation, chroms=chroms)
    anno_keys = {(c, s, ch) for c, s, ch in anno.values()}
    reads = read_chains(a.bam, chroms)
    by_locus = collections.defaultdict(list)
    for (c, s, ch), k in reads.items():
        if ch:
            by_locus[(c, s)].append((ch, k))
    nc = {r["transcript_name"]: float(r["est_count"])
          for r in csv.DictReader(open(a.nanocount), delimiter="\t")}
    pe, se = eq_set(a.pyfin_tracking), eq_set(a.st3_tracking)
    pc, sc = ref_codes(a.pyfin_tracking), ref_codes(a.st3_tracking)
    cov = st3_cov(a.st3_tracking)

    def suffix_support(c, s, ch):
        tot = 0
        for rch, k in by_locus.get((c, s), ()):
            if not rch or len(rch) > len(ch):
                continue
            seg = ch[len(ch) - len(rch):] if s == "+" else ch[:len(rch)]
            if seg == rch:
                tot += k
        return tot

    rows = []
    for t, (c, s, ch) in truth.items():
        grp = ("both" if t in pe and t in se else "st3_only" if t in se
               else "pyfin_only" if t in pe else "neither")
        exact = reads.get((c, s, ch), 0)
        suf = suffix_support(c, s, ch) if ch else 0
        rows.append(dict(tx=t, group=grp, chrom=c, strand=s, n_introns=len(ch),
                         in_annotation=int((c, s, ch) in anno_keys),
                         exact_reads=exact, suffix_reads=suf,
                         nanocount=nc.get(t, float("nan")),
                         st3_cov=cov.get(t, float("nan")),
                         pyfin_codes=",".join(sorted(pc.get(t, ()))) or "-",
                         st3_codes=",".join(sorted(sc.get(t, ()))) or "-"))
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    for grp in ("st3_only", "pyfin_only", "both", "neither"):
        g = [r for r in rows if r["group"] == grp]
        if not g:
            continue
        n = len(g)
        mono = sum(1 for r in g if r["n_introns"] == 0)
        ex0 = sum(1 for r in g if r["exact_reads"] == 0)
        ex1 = sum(1 for r in g if r["exact_reads"] == 1)
        ex2 = sum(1 for r in g if r["exact_reads"] >= 2)
        suf_only = sum(1 for r in g if r["exact_reads"] == 0 and r["suffix_reads"] > 0)
        none = sum(1 for r in g if r["exact_reads"] == 0 and r["suffix_reads"] == 0
                   and r["n_introns"] > 0)
        print(f"\n== {grp}: n={n}  mono={mono}  in_annotation={sum(r['in_annotation'] for r in g)}")
        print(f"   exact full-chain reads: 0:{ex0} ({100*ex0/n:.0f}%)  1:{ex1}  >=2:{ex2}")
        print(f"   of exact=0 multi-exon: suffix-only {suf_only}, no read at all {none}")
        print(f"   nanocount est_count median {q([r['nanocount'] for r in g], .5):.2f}, "
              f"<1: {sum(1 for r in g if r['nanocount'] < 1)}")
        cv = [r["st3_cov"] for r in g if r["st3_cov"] == r["st3_cov"]]
        if cv:
            print(f"   stringtie3 cov: median {q(cv, .5):.2f}  q25 {q(cv, .25):.2f}  "
                  f"q75 {q(cv, .75):.2f}  <1: {sum(1 for x in cv if x < 1)}  <2: {sum(1 for x in cv if x < 2)}")
        for lab, key in (("pyfin codes", "pyfin_codes"), ("st3 codes", "st3_codes")):
            cc = collections.Counter(r[key] for r in g)
            print(f"   {lab}: {cc.most_common(8)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
