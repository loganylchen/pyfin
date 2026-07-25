#!/usr/bin/env python3
"""Compare the GENERATION candidate set of the two clustering primitives -- "up to this step",
before EM/selection -- so the comparison isolates clustering+fold from the (read_chains-tuned)
downstream.

Runs the REAL production generation (`discover_candidates`, chain_cluster=True) twice per interval:
  read_chains: legacy cluster_read_chains (fold_span_guard on, fold_monoexon off, matching the
               --mono-resolve-post-em end-to-end config)
  families   : cluster_families + collapse (the redesign)
and writes each mode's post-fold candidate set as a GTF. gffcompare each vs truth to get the
generation-level Sn (recall ceiling) + Pr (precision floor) of each clustering primitive.

Usage: measure_gen_compare.py <bam> <genome.fa> <out_read_chains.gtf> <out_families.gtf>
"""
import sys

import pysam

from fin.candidates.discovery import _exons_from_chain, discover_candidates
from fin.io.interval_manager import generate_isolated_intervals

BAM, GENOME, OUT_RC, OUT_FAM = sys.argv[1:5]


def write_tx(fh, chrom, strand, exons, tid):
    exons = sorted(exons)
    st, en = exons[0][0], exons[-1][1]
    attr = f'gene_id "{tid}"; transcript_id "{tid}";'
    fh.write(f"{chrom}\tgen\ttranscript\t{st + 1}\t{en}\t.\t{strand}\t.\t{attr}\n")
    for i, (s, e) in enumerate(exons, 1):
        fh.write(f'{chrom}\tgen\texon\t{s + 1}\t{e}\t.\t{strand}\t.\t{attr} exon_number "{i}";\n')


def emit(fh, cs, tag, counter):
    for c in cs.candidates:
        if c.source != "novel":            # p00: all novel; skip any GTF passthrough
            continue
        exons = _exons_from_chain(c.start, c.end, c.intron_chain)
        counter[0] += 1
        write_tx(fh, c.chrom, c.strand, exons, f"{tag}_{counter[0]}")


def run_mode(interval, chrom_seq, clustering):
    return discover_candidates(
        interval=interval, bam_path=BAM, gtf_reader=None, genome_fasta=chrom_seq,
        threshold=24, min_novel_reads=1, chain_cluster=True, clustering=clustering,
        chain_cluster_wobble_bp=6, chain_cluster_cassette_max_exon_bp=70,
        chain_cluster_fold_monoexon=False,        # match --mono-resolve-post-em (mono deferred)
        chain_cluster_fold_span_guard=True)       # read_chains production default


def main():
    intervals = generate_isolated_intervals(BAM)["intervals"]
    fa = pysam.FastaFile(GENOME)
    n_rc = [0]
    n_fam = [0]
    seq_cache = {}
    with open(OUT_RC, "w") as f_rc, open(OUT_FAM, "w") as f_fam:
        for iv in intervals:
            if iv.chrom not in seq_cache:
                seq_cache.clear()
                try:
                    seq_cache[iv.chrom] = fa.fetch(iv.chrom)
                except (KeyError, ValueError):
                    seq_cache[iv.chrom] = ""
            chrom_seq = seq_cache[iv.chrom]
            emit(f_rc, run_mode(iv, chrom_seq, "read_chains"), "rc", n_rc)
            emit(f_fam, run_mode(iv, chrom_seq, "families"), "fam", n_fam)
    sys.stderr.write(f"[gen_compare] read_chains novel={n_rc[0]}  families novel={n_fam[0]}\n")


if __name__ == "__main__":
    main()
