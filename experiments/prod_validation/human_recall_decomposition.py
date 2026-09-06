#!/usr/bin/env python3
"""Where does human recall go?

On the human slice PyFIN recalls 25% of the expressed transcripts while the
annotation itself contains 61% of them. This decomposes the gap into things
that are fixable and things that are not, by asking of every EXPRESSED truth
transcript, in order:

  in_annotation   - is its structure in the input annotation at all?
  has_reads       - does any primary read's CIGAR chain match its intron chain
                    exactly (a read that could have generated it)?
  full_chain>=K   - do >= K reads attest the COMPLETE chain?
  emitted         - did we output a model gffcompare matched to it?

A transcript with no read that carries its chain is a coverage limit, not a
pipeline defect; one with K attesting reads that we did not emit is ours to fix.
"""
import argparse, collections, csv, re
from pathlib import Path

import pysam


def chains_from_gtf(path, keep=None, chroms=None):
    """transcript_id -> (chrom, strand, intron chain tuple)."""
    exons = collections.defaultdict(list)
    meta = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            if chroms and f[0] not in chroms:
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m or (keep is not None and m.group(1) not in keep):
                continue
            t = m.group(1)
            exons[t].append((int(f[3]) - 1, int(f[4])))
            meta[t] = (f[0], f[6])
    out = {}
    for t, ex in exons.items():
        ex.sort()
        out[t] = (meta[t][0], meta[t][1],
                  tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)))
    return out


def read_chains(bam_path, chroms):
    """(chrom, strand, chain) -> number of primary reads carrying it exactly."""
    counts = collections.Counter()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam.fetch():
            if (r.is_unmapped or r.is_secondary or r.is_supplementary
                    or r.reference_name not in chroms):
                continue
            pos, introns = r.reference_start, []
            for op, ln in (r.cigartuples or []):
                if op in (0, 2, 7, 8):
                    pos += ln
                elif op == 3:
                    introns.append((pos, pos + ln)); pos += ln
            counts[(r.reference_name, "-" if r.is_reverse else "+",
                    tuple(introns))] += 1
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--tmap", required=True, help="gffcompare tmap of the run")
    ap.add_argument("--chroms", default="chr1,chr19,chr21")
    ap.add_argument("--k", type=int, default=2)
    a = ap.parse_args()
    chroms = set(a.chroms.split(","))
    truth = chains_from_gtf(a.truth, chroms=chroms)
    anno = chains_from_gtf(a.annotation, chroms=chroms)
    anno_keys = {(c, s, ch) for c, s, ch in anno.values()}
    reads = read_chains(a.bam, chroms)
    matched = {r["ref_id"] for r in csv.DictReader(open(a.tmap), delimiter="\t")
               if r["class_code"] == "="}

    # 5'-degraded dRNA reads carry a 3' SUFFIX of the transcript's chain: index
    # every observed read chain by its (chrom, strand) so a truth chain can ask
    # how many reads carry a contiguous suffix of it.
    by_locus = collections.defaultdict(list)
    for (c, s, ch), k in reads.items():
        if ch:
            by_locus[(c, s)].append((ch, k))

    def suffix_support(c, s, ch):
        """reads whose chain is a contiguous 3' suffix of ``ch`` (strand-aware)."""
        if not ch:
            return 0, 0
        tot = best = 0
        for rch, k in by_locus.get((c, s), ()):
            if not rch or len(rch) > len(ch):
                continue
            # on '+' the 3' end is the chain's tail; on '-' it is its head
            seg = ch[len(ch) - len(rch):] if s == "+" else ch[:len(rch)]
            if seg == rch:
                tot += k
                best = max(best, len(rch))
        return tot, best

    n = collections.Counter()
    for t, (c, s, ch) in truth.items():
        n["expressed"] += 1
        key = (c, s, ch)
        in_anno = key in anno_keys
        nr = reads.get(key, 0)
        emitted = t in matched
        n["in_annotation"] += in_anno
        n["has_read_chain"] += nr > 0
        n[f"reads>={a.k}"] += nr >= a.k
        n["emitted"] += emitted
        if nr >= a.k and not emitted:
            n["LOST: attested but not emitted"] += 1
        if nr == 0 and emitted:
            n["emitted without an exact read chain"] += 1
        if nr == 0 and not emitted:
            n["no read carries the chain"] += 1
        if 0 < nr < a.k and not emitted:
            n[f"only 1..{a.k-1} reads"] += 1
        if nr == 0 and ch:
            sup, longest = suffix_support(c, s, ch)
            if sup >= a.k:
                n[f"  of those: >={a.k} reads carry a 3' SUFFIX"] += 1
                if longest >= max(1, len(ch) - 1):
                    n["    suffix missing only the 5'-most junction"] += 1
            elif sup > 0:
                n["  of those: 1 suffix read"] += 1
            else:
                n["  of those: no suffix read either"] += 1
        elif nr == 0 and not ch:
            n["  of those: mono-exon truth"] += 1
    tot = n["expressed"]
    for k in ("expressed", "in_annotation", "has_read_chain", f"reads>={a.k}",
              "emitted", "LOST: attested but not emitted",
              f"only 1..{a.k-1} reads", "no read carries the chain",
              f"  of those: >={a.k} reads carry a 3' SUFFIX",
              "    suffix missing only the 5'-most junction",
              "  of those: 1 suffix read", "  of those: no suffix read either",
              "  of those: mono-exon truth",
              "emitted without an exact read chain"):
        print(f"{k:38s} {n[k]:7d}  {100*n[k]/max(tot,1):5.1f}%")


if __name__ == "__main__":
    main()
