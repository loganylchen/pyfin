#!/usr/bin/env python3
"""What would admitting gate-dropped annotated candidates buy, by evidence tier?

Joins the selection audit of one cell with the BAM's read chains and the
expression truth: for every annotated candidate the m2 support gate dropped,
count exact full-chain reads and 3'-suffix reads, and whether the transcript
is expressed (truth). Reports, per evidence tier, how many candidates a rule
admitting that tier would add and what fraction of them are truth -- an
offline precision estimate for complete-chain rescue at K=1,2,3 and for the
suffix-only tail. EM redistribution and downstream floors are NOT replayed,
so this bounds, not measures, an end-to-end arm.
"""
import argparse, collections, csv, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from human_recall_decomposition import chains_from_gtf, read_chains


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit-merged", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--chroms", default="chr1,chr19,chr21")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    chroms = set(a.chroms.split(","))
    anno = chains_from_gtf(a.annotation, chroms=chroms)
    truth = set(chains_from_gtf(a.truth, chroms=chroms))
    reads = read_chains(a.bam, chroms)
    by_locus = collections.defaultdict(list)
    for (c, s, ch), k in reads.items():
        if ch:
            by_locus[(c, s)].append((ch, k))

    def suffix_support(c, s, ch):
        tot = 0
        for rch, k in by_locus.get((c, s), ()):
            if not rch or len(rch) > len(ch):
                continue
            seg = ch[len(ch) - len(rch):] if s == "+" else ch[:len(rch)]
            if seg == rch:
                tot += k
        return tot

    rows = list(csv.DictReader(open(a.audit_merged), delimiter="\t"))
    discovered = {r["candidate_id"]: r for r in rows if r["stage"] == "discovered"}
    dropped = {r["candidate_id"]: r for r in rows
               if r["stage"] == "interval_cascade" and r["action"] == "drop"
               and r["reason"] == "m2_support_gate"}
    final = {r["candidate_id"] for r in rows if r["stage"] == "final"}
    out = []
    for cid, r in dropped.items():
        if cid not in anno:
            continue  # novel
        c, s, ch = anno[cid]
        exact = reads.get((c, s, ch), 0)
        suf = suffix_support(c, s, ch) if ch else 0
        out.append(dict(candidate_id=cid, n_introns=len(ch), exact_reads=exact,
                        suffix_reads=suf, disc_reads=discovered.get(cid, {}).get("n_reads", ""),
                        truth=int(cid in truth)))
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader(); w.writerows(out)

    n_final_truth = sum(1 for cid in final if cid in truth)
    print(f"annotated candidates dropped by m2_support_gate: {len(out)}  "
          f"(truth {sum(o['truth'] for o in out)}); final emitted {len(final)} (annotated truth in final {n_final_truth})")
    tiers = [("exact>=3", lambda o: o["exact_reads"] >= 3),
             ("exact==2", lambda o: o["exact_reads"] == 2),
             ("exact==1", lambda o: o["exact_reads"] == 1),
             ("exact==0 & suffix>=3", lambda o: o["exact_reads"] == 0 and o["suffix_reads"] >= 3),
             ("exact==0 & suffix 1-2", lambda o: o["exact_reads"] == 0 and 0 < o["suffix_reads"] < 3),
             ("no chain read (multi-exon)", lambda o: o["n_introns"] > 0 and o["exact_reads"] == 0 and o["suffix_reads"] == 0),
             ("mono-exon", lambda o: o["n_introns"] == 0)]
    print(f"{'tier':28s} {'n':>6s} {'truth':>6s} {'prec':>6s}")
    for name, f in tiers:
        g = [o for o in out if f(o)]
        t = sum(o["truth"] for o in g)
        print(f"{name:28s} {len(g):6d} {t:6d} {100*t/max(1,len(g)):5.1f}%")
    print("\ncumulative admission (multi-exon only), F1 break-even needs marginal precision > F1/2 ~ 22%:")
    cum_n = cum_t = 0
    for name, f in tiers[:5]:
        g = [o for o in out if f(o) and o["n_introns"] > 0]
        cum_n += len(g); cum_t += sum(o["truth"] for o in g)
        print(f"  + {name:26s} cum n={cum_n:5d} cum truth={cum_t:5d} cum prec={100*cum_t/max(1,cum_n):5.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
