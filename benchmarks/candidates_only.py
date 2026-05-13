"""Run only interval discovery + candidate discovery; print summary.

Usage:
    python benchmarks/candidates_only.py --bam X --genome Y [--gtf Z] --out out.tsv
"""
from __future__ import annotations

import argparse
import logging
import sys
from collections import Counter

from fin.candidates.discovery import discover_candidates
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_fasta import FASTAReader
from fin.io.io_gtf import GTFReader

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
logger = logging.getLogger(__name__)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", default=None)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threshold", type=int, default=24)
    ap.add_argument("--canonical-search-bp", type=int, default=0,
                    help="Scan +/- N bp around each junction for canonical "
                         "GT-AG and emit chain alternatives. 0 disables.")
    ap.add_argument("--max-chains-per-read", type=int, default=16,
                    help="Cap on per-read alternative chains during canonical "
                         "search.")
    args = ap.parse_args()

    gtf_reader = None
    if args.gtf:
        gtf_reader = GTFReader(args.gtf)
        gtf_reader.open()
        gtf_reader.parse()

    fasta = FASTAReader(args.genome)
    genome = {rec.id: rec.sequence for rec in fasta.get_records()}

    result = generate_isolated_intervals(
        args.bam, gtf_path=args.gtf, max_gap=200, max_reads=None
    )
    intervals = result["intervals"]
    logger.info("Generated %d intervals", len(intervals))

    rows = []
    source_counts = Counter()
    total = 0
    for i, iv in enumerate(intervals):
        chrom_seq = genome.get(iv.chrom, "")
        cs = discover_candidates(
            interval=iv,
            bam_path=args.bam,
            gtf_reader=gtf_reader,
            genome_fasta=chrom_seq,
            threshold=args.threshold,
            min_novel_reads=1,
            canonical_search_bp=args.canonical_search_bp,
            max_chains_per_read=args.max_chains_per_read,
        )
        for c in cs.candidates:
            source_counts[c.source] += 1
            total += 1
            rows.append((
                c.candidate_id,
                c.source,
                c.chrom,
                c.strand,
                c.start,
                c.end,
                c.three_prime_pos,
                len(c.intron_chain.introns),
                ";".join(f"{a}-{b}" for a, b in c.intron_chain.introns),
                len(c.supporting_read_ids),
            ))

    with open(args.out, "w") as fh:
        fh.write("candidate_id\tsource\tchrom\tstrand\tstart\tend\tthree_prime\tn_introns\tintron_chain\tn_reads\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    logger.info("Wrote %d candidates to %s", total, args.out)
    logger.info("By source: %s", dict(source_counts))
    logger.info("Intervals: %d", len(intervals))


if __name__ == "__main__":
    main()
