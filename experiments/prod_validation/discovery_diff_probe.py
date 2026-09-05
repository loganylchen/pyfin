#!/usr/bin/env python3
"""Directly diff discover_candidates() emissions with and without a corrupt GTF.

For each target chromosome/strand, call the production discovery (families
clustering, production knobs) once with gtf_reader=None and once with the
condition annotation, then report whether each target TRUE chain appears among
emitted candidate chains. Pinpoints the exact stage where truncation-corrupt
annotation suppresses read-derived true chains. Read-only probe.
"""
from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, "/SSD/logan/dev/pyfin")
from fin.candidates.discovery import discover_candidates  # noqa: E402
from fin.io.interval_manager import GenomicInterval  # noqa: E402
from fin.io.io_gtf import GTFReader  # noqa: E402

BENCH = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/"
             "results/sirv4_full_sweep_ratio_grid_v2")
W = Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid")
S = "SGNex_H9_directRNA_replicate3_run1"
BAM = BENCH / S / "align" / f"{S}.sorted.bam"
FA = W / "_ref" / "SIRV_Set4.fa"

TARGETS = {
    "c_trunc3p50": ["SIRV204", "SIRV307", "SIRV409", "SIRV410", "SIRV501",
                    "SIRV502", "SIRV506", "SIRV507"],
    "c_trunc5p50": ["SIRV307", "SIRV409", "SIRV511", "SIRV611", "SIRV616",
                    "SIRV702", "SIRV704"],
    "c_merge10": ["SIRV204", "SIRV405", "SIRV507"],
}


def parse_truth(path):
    ex = defaultdict(list)
    meta = {}
    for line in Path(path).open():
        f = line.rstrip().split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        ex[m.group(1)].append((int(f[3]), int(f[4])))
        meta[m.group(1)] = (f[0], f[6])
    out = {}
    for tid, e in ex.items():
        e = sorted(e)
        chain = tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1))
        out[tid] = (*meta[tid], chain)
    return out


def load_fasta(path):
    seqs = {}
    name = None
    buf = []
    for line in Path(path).open():
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].split()[0]
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def emitted_chains(gtf_reader, chrom, strand, seq):
    iv = GenomicInterval(chrom=chrom, start=0, end=len(seq), strand=strand)
    cs = discover_candidates(
        iv, str(BAM), gtf_reader, seq,
        chain_cluster=True, chain_cluster_wobble_bp=6,
        chain_cluster_cassette_max_exon_bp=70,
        chain_cluster_fold_monoexon=False,
        chain_cluster_fold_span_guard=True, clustering="families",
    )
    chains = defaultdict(list)
    for c in cs.candidates:
        chains[tuple(c.intron_chain.introns)].append(
            (c.candidate_id, c.source, len(c.supporting_read_ids)))
    return chains


def main() -> int:
    truth = parse_truth(W / "_ref" / "SIRV_Set4.gtf")
    seqs = load_fasta(FA)
    for cond, tids in TARGETS.items():
        anno = BENCH / "_ref" / S / cond / "annotation.gtf"
        reader = GTFReader(str(anno))
        combos = sorted({truth[t][:2] for t in tids})
        for chrom, strand in combos:
            wanted = [t for t in tids
                      if truth[t][0] == chrom and truth[t][1] == strand]
            base = emitted_chains(None, chrom, strand, seqs[chrom])
            wgtf = emitted_chains(reader, chrom, strand, seqs[chrom])
            for t in wanted:
                ch = truth[t][2]
                b = base.get(ch)
                g = wgtf.get(ch)
                print("%-11s %-8s %s%-2s p00:%-28s gtf:%s" % (
                    cond, t, chrom, strand,
                    b[:1] if b else "ABSENT", g[:1] if g else "ABSENT"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
