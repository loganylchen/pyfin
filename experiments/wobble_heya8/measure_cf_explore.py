#!/usr/bin/env python3
"""Structural recall/precision ceiling of the clustering-only step (cluster_families).

Emit EVERY family variant (+ optional per-locus mono) as a transcript -- NO fold, NO
selection, NO EM, NO quant. This is "stop at cluster_families and keep everything": the
recall ceiling (nothing folded/dropped) and the precision floor (every wobble variant kept)
of the new clustering. gffcompare the resulting GTF vs truth to score it.

Reads are handled like production generation: primary only, strand-filtered into the
interval, fusion-like (>=250bp terminal soft-clip) reads excluded.

Usage: measure_cf.py <bam> <out.gtf> [--mono] [--wobble N]
"""
import sys

import pysam

from fin.candidates.chain_cluster import cluster_families
from fin.candidates.explore import explore_family_paths
from fin.candidates.dataclasses import IntronChain
from fin.candidates.discovery import _exons_from_chain, _spatial_read_clusters
from fin.candidates.intron_chains import extract_intron_chain
from fin.io.interval_manager import generate_isolated_intervals

BAM = sys.argv[1]
OUT = sys.argv[2]
EMIT_MONO = "--mono" in sys.argv[3:]
WOBBLE = 6
if "--wobble" in sys.argv:
    WOBBLE = int(sys.argv[sys.argv.index("--wobble") + 1])
SOFT_CLIP = 250


def write_tx(fh, chrom, strand, exons, tid):
    exons = sorted(exons)
    st, en = exons[0][0], exons[-1][1]
    attr = f'gene_id "{tid}"; transcript_id "{tid}";'
    fh.write(f"{chrom}\tcf\ttranscript\t{st + 1}\t{en}\t.\t{strand}\t.\t{attr}\n")
    for i, (s, e) in enumerate(exons, 1):
        fh.write(f'{chrom}\tcf\texon\t{s + 1}\t{e}\t.\t{strand}\t.\t{attr} exon_number "{i}";\n')


def main():
    intervals = generate_isolated_intervals(BAM)["intervals"]
    bam = pysam.AlignmentFile(BAM, "rb")
    n_var = n_mono = tidx = 0
    with open(OUT, "w") as fh:
        for iv in intervals:
            strand = iv.strand or "+"
            read_chains = []
            spans = {}
            region = f"{iv.chrom}:{max(iv.start + 1, 1)}-{iv.end}"
            try:
                it = bam.fetch(region=region)
            except (ValueError, KeyError):
                continue
            for aln in it:
                if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                    continue
                rstrand = "-" if aln.is_reverse else "+"
                if iv.strand is not None and rstrand != iv.strand:
                    continue
                ct = aln.cigartuples
                if not ct:
                    continue
                # exclude fusion-like reads (>=250bp terminal soft-clip), like generation
                if (ct[0][0] == 4 and ct[0][1] >= SOFT_CLIP) or \
                   (ct[-1][0] == 4 and ct[-1][1] >= SOFT_CLIP):
                    continue
                q = aln.query_name
                chain = extract_intron_chain(ct, aln.reference_start)
                read_chains.append(({"query_name": q}, chain))
                spans[q] = (aln.reference_start, aln.reference_end)
            if not read_chains:
                continue
            cr = cluster_families(read_chains, wobble_bp=WOBBLE, cassette_max_exon_bp=70)
            for fam in cr.families:
                for v in fam.variants:
                    reads = fam.variant_reads.get(v, set())
                    st = [spans[r][0] for r in reads if r in spans]
                    en = [spans[r][1] for r in reads if r in spans]
                    if not st:
                        continue
                    exons = _exons_from_chain(min(st), max(en), IntronChain(introns=v))
                    tidx += 1
                    n_var += 1
                    write_tx(fh, iv.chrom, strand, exons, f"cf_{tidx}")
            # explore: emit assembled cross-read full-length candidates (additive)
            for ep in explore_family_paths(fam):
                st = [spans[r][0] for r in ep.reads if r in spans]
                en = [spans[r][1] for r in ep.reads if r in spans]
                if not st:
                    continue
                exons = _exons_from_chain(min(st), max(en), IntronChain(introns=ep.chain))
                tidx += 1
                n_var += 1
                write_tx(fh, iv.chrom, strand, exons, f"cfx_{tidx}")
            if EMIT_MONO:
                for grp in _spatial_read_clusters(cr.mono_reads, spans):
                    st = [spans[r][0] for r in grp]
                    en = [spans[r][1] for r in grp]
                    if not st:
                        continue
                    tidx += 1
                    n_mono += 1
                    write_tx(fh, iv.chrom, strand, [(min(st), max(en))], f"cfmono_{tidx}")
    sys.stderr.write(
        f"[measure_cf_explore] wobble={WOBBLE} variants={n_var} mono={n_mono} "
        f"total={n_var + n_mono}\n")


if __name__ == "__main__":
    main()
