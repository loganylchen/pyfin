#!/usr/bin/env python3
"""How precise is the annotation's own TSS/TES, measured from the reads?

The extent-contradiction test asks whether a read's end is consistent with an
annotated end under a kernel fitted to READ-END measurement noise (b_narrow
~1.9 bp on human). That is the wrong yardstick for the question "is the
ANNOTATION's end right": annotation ends carry their own error, and if that
error is tens of bases then a read cluster 15 bp away is not evidence against
the annotation, it is inside the annotation's own uncertainty.

This estimates the annotation's end precision WITHOUT using the contradiction
test, so it is independent of the mechanism under review. For every annotated
multi-exon transcript whose intron chain is carried exactly by >= K primary
reads (so the reads are unambiguously that transcript's, not a nested
sibling's), take the MODE of those reads' 5' ends (and 3' ends) -- the modal
position of a cluster, not the mean, so degradation tails do not move it -- and
report the signed offset of that mode from the annotated end.

Two selection controls are applied, because the naive version measures dRNA
degradation instead of annotation error:
  * the 5' offset is reported both over all such transcripts and over the
    subset whose read 5' ends are TIGHT (>= `frac` of them within `tol`), where
    a shared start is more likely to be a real TSS than a degradation stop;
  * transcripts whose gene has another annotated TSS within `tol` of the mode
    are reported separately, since there the read may belong to the sibling.
"""
import argparse, collections, csv, re, statistics, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from human_recall_decomposition import chains_from_gtf


def modal_position(vals, tol):
    """Position of the densest cluster: the value with the most neighbours
    within +-tol, ties to the smallest. Returns (mode, n_in_cluster)."""
    best, best_n = None, -1
    for v in sorted(vals):
        n = sum(1 for x in vals if abs(x - v) <= tol)
        if n > best_n:
            best, best_n = v, n
    return best, best_n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--chroms", default="chr1,chr19,chr21")
    ap.add_argument("--min-reads", type=int, default=3)
    ap.add_argument("--tol", type=int, default=10)
    ap.add_argument("--frac", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sys.path.insert(0, "/SSD/logan/dev/pyfin")
    import pysam

    chroms = set(a.chroms.split(","))
    anno = chains_from_gtf(a.annotation, chroms=chroms)
    # annotated ends and the gene's other ends
    ex = collections.defaultdict(list)
    meta = {}
    with open(a.annotation) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon" or f[0] not in chroms:
                continue
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            g = re.search(r'gene_id "([^"]+)"', f[8])
            ex[t].append((int(f[3]) - 1, int(f[4])))
            meta[t] = (f[0], f[6], g.group(1) if g else "")
    ends = {}
    gene_tss, gene_tes = collections.defaultdict(set), collections.defaultdict(set)
    for t, e in ex.items():
        e.sort()
        c, s, g = meta[t]
        tss = e[0][0] if s == "+" else e[-1][1]
        tes = e[-1][1] if s == "+" else e[0][0]
        ends[t] = (c, s, g, tss, tes)
        gene_tss[(c, s, g)].add(tss)
        gene_tes[(c, s, g)].add(tes)

    # reads: chain -> [(start, end)] for MAPQ>=1 primaries
    want = {}
    for t, (c, s, ch) in anno.items():
        if ch:
            want.setdefault((c, s, ch), []).append(t)
    by_chain = collections.defaultdict(list)
    bam = pysam.AlignmentFile(a.bam, "rb")
    for aln in bam.fetch():
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if aln.mapping_quality < 1 or aln.reference_name not in chroms:
            continue
        pos, chain = aln.reference_start, []
        for op, ln in aln.cigartuples or ():
            if op in (0, 7, 8, 2):
                pos += ln
            elif op == 3:
                chain.append((pos, pos + ln))
                pos += ln
        if not chain:
            continue
        key = (aln.reference_name, "-" if aln.is_reverse else "+", tuple(chain))
        if key in want:
            by_chain[key].append((aln.reference_start, aln.reference_end))
    bam.close()

    rows = []
    for key, spans in by_chain.items():
        if len(spans) < a.min_reads:
            continue
        c, s, ch = key
        for t in want[key]:
            _, _, g, tss, tes = ends[t]
            starts = [x[0] if s == "+" else x[1] for x in spans]
            stops = [x[1] if s == "+" else x[0] for x in spans]
            m5, n5 = modal_position(starts, a.tol)
            m3, n3 = modal_position(stops, a.tol)
            d5 = (m5 - tss) if s == "+" else (tss - m5)
            d3 = (m3 - tes) if s == "+" else (tes - m3)
            sib5 = min((abs(m5 - x) for x in gene_tss[(c, s, g)] if x != tss), default=-1)
            sib3 = min((abs(m3 - x) for x in gene_tes[(c, s, g)] if x != tes), default=-1)
            rows.append(dict(tx=t, n_reads=len(spans), tight5=int(n5 >= a.frac * len(spans)),
                             tight3=int(n3 >= a.frac * len(spans)),
                             d5=d5, d3=d3, sib5=sib5, sib3=sib3))
    if not rows:
        raise SystemExit("no annotated chain had enough exact reads")
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)

    def summ(name, vals):
        if not vals:
            print(f"  {name:44s} n=0")
            return
        v = sorted(abs(x) for x in vals)
        sgn = sorted(vals)
        q = lambda p: v[int(p * (len(v) - 1))]
        print(f"  {name:44s} n={len(v):5d}  |d| median {q(.5):4.0f}  p75 {q(.75):5.0f}  "
              f"p90 {q(.90):6.0f}  <=2bp {100*sum(1 for x in v if x<=2)/len(v):4.0f}%  "
              f"<=10 {100*sum(1 for x in v if x<=10)/len(v):4.0f}%  <=25 {100*sum(1 for x in v if x<=25)/len(v):4.0f}%  "
              f"signed median {statistics.median(sgn):+.0f}")

    print(f"annotated chains with >= {a.min_reads} exact primary reads: {len(rows)}")
    print("\noffset of the READ-CLUSTER MODE from the annotated end "
          "(negative 5' = cluster upstream of the annotation):")
    summ("5' all", [r["d5"] for r in rows])
    summ("5' tight clusters only", [r["d5"] for r in rows if r["tight5"]])
    summ("5' tight AND no sibling TSS within 25bp",
         [r["d5"] for r in rows if r["tight5"] and not (0 <= r["sib5"] <= 25)])
    summ("3' all", [r["d3"] for r in rows])
    summ("3' tight clusters only", [r["d3"] for r in rows if r["tight3"]])
    summ("3' tight AND no sibling TES within 25bp",
         [r["d3"] for r in rows if r["tight3"] and not (0 <= r["sib3"] <= 25)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
