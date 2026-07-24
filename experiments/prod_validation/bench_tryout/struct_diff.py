#!/usr/bin/env python3
"""Structural byte-identity gate for the refactor: compare two pyfin GTFs by the
MULTISET of canonical transcript records, ignoring the random novel UUIDs
(gene_id / transcript_id). Two GTFs are "structurally identical" iff every
transcript maps 1:1 on (chrom, strand, source, num_reads, abundance, exon coords).

Usage: struct_diff.py A.gtf B.gtf
Exit 0 = identical, 1 = differences (with a sample printed).
"""
import re
import sys
from collections import Counter

ATTR = re.compile(r'(\w+) "([^"]*)"')


def parse(path):
    """Return Counter{canonical_transcript_tuple: n}."""
    tx = {}          # tid -> [chrom, strand, source, num_reads, abundance]
    exons = {}       # tid -> list[(start,end)]
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or "\t" not in line:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _src, feat, start, end, _sc, strand, _fr, attrs = f[:9]
            a = dict(ATTR.findall(attrs))
            tid = a.get("transcript_id", "")
            if feat == "transcript":
                tx[tid] = (
                    chrom, strand,
                    a.get("transcript_source", ""),
                    a.get("num_reads", ""),
                    a.get("abundance", ""),
                )
            elif feat == "exon":
                exons.setdefault(tid, []).append((int(start), int(end)))
    recs = Counter()
    for tid, base in tx.items():
        ex = tuple(sorted(exons.get(tid, [])))
        recs[base + (ex,)] += 1
    return recs, len(tx)


def main():
    a, b = sys.argv[1], sys.argv[2]
    ca, na = parse(a)
    cb, nb = parse(b)
    print(f"A={a}  transcripts={na}  distinct_records={len(ca)}")
    print(f"B={b}  transcripts={nb}  distinct_records={len(cb)}")
    only_a = ca - cb
    only_b = cb - ca
    da = sum(only_a.values())
    db = sum(only_b.values())
    if da == 0 and db == 0:
        print("STRUCTURALLY IDENTICAL (multiset match, ignoring novel UUIDs)")
        return 0
    print(f"DIFFERENT: only-in-A={da}  only-in-B={db}")
    for tag, c in (("A", only_a), ("B", only_b)):
        for rec, n in list(c.items())[:8]:
            print(f"  only-{tag} x{n}: chrom={rec[0]} strand={rec[1]} src={rec[2]} "
                  f"nreads={rec[3]} ab={rec[4]} exons={rec[5][:4]}{'...' if len(rec[5])>4 else ''}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
