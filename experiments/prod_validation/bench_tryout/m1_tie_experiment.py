#!/usr/bin/env python3
"""On NEW chain-cluster candidates, apply ONLY the M1-tie support rule (drop any
candidate that is in NO read's best-AS tie set) and compare to no gate at all.

M1 = production mappy AS (get_m1_preset + score_hit). Per interval: build a mappy
aligner per candidate, align every read, each read's best-AS tie = candidates within
1e-9 of its max AS. A candidate kept iff it appears in >=1 read's tie set. NO EM, NO
other gates. Writes NEW_nogate.gtf (all novel) and NEW_m1tie.gtf (M1-tie survivors).
"""
import os, sys
import numpy as np
import mappy
import pysam

REPO = "/SSD/logan/dev/pyfin"
sys.path.insert(0, REPO)
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_gtf import GTFReader
from fin.candidates.discovery import discover_candidates, _exons_from_chain
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S = "SGNex_H9_directRNA_replicate2_run2"
BAM = f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
GENOME = f"{B}/refs/host/GRCh38.primary_assembly.genome.fa"
GTF = f"{B}/results/gencode_full_sweep/_ref/p00/annotation.gtf"
OUT = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/cluster_cmp"
LIMIT = int(os.environ.get("LIMIT", "0"))
os.makedirs(OUT, exist_ok=True)

fa = pysam.FastaFile(GENOME)
_cc = {}
def chrom_seq(c):
    if c not in _cc:
        try: _cc[c] = fa.fetch(c)
        except Exception: _cc[c] = ""
    return _cc[c]

gtf_reader = GTFReader(GTF); gtf_reader.open(); gtf_reader.parse()
intervals = generate_isolated_intervals(BAM, gtf_path=GTF, max_gap=0)["intervals"]
if LIMIT: intervals = intervals[:LIMIT]
sys.stderr.write(f"intervals: {len(intervals)}\n")
preset = get_m1_preset()

nogate, m1tie = [], []
for n, iv in enumerate(intervals):
    cs = chrom_seq(iv.chrom)
    cset = discover_candidates(
        iv, BAM, gtf_reader, cs, threshold=24, min_novel_reads=1,
        chain_cluster=True, canonical_search_bp=0)
    cand_list = list(cset.candidates)
    if not cand_list:
        continue
    read_seqs = getattr(cset, "read_sequences", {}) or {}
    reads = [(r, s) for r, s in read_seqs.items() if s]
    # M1: mappy AS matrix (reads x candidates)
    aligners = [mappy.Aligner(seq=c.sequence, preset=preset) if c.sequence else None
                for c in cand_list]
    in_tie = [False] * len(cand_list)
    for rid, seq in reads:
        best = -1e18
        row = [-1e18] * len(cand_list)
        for j, aln in enumerate(aligners):
            if aln is None:
                continue
            b = None
            for h in aln.map(seq):
                v = score_hit(h)
                if v is not None and (b is None or v > b):
                    b = v
            if b is not None:
                row[j] = b
                if b > best:
                    best = b
        if best <= 0:
            continue
        for j in range(len(cand_list)):
            if row[j] > 0 and row[j] >= best - 1e-9:
                in_tie[j] = True
    for j, c in enumerate(cand_list):
        if c.source != "novel":
            continue
        rec = (c.chrom, c.strand, c.start, c.end, c.intron_chain)
        nogate.append(rec)
        if in_tie[j]:
            m1tie.append(rec)
    if (n + 1) % 2000 == 0:
        sys.stderr.write(f"  {n+1}/{len(intervals)} nogate={len(nogate)} m1tie={len(m1tie)}\n")


def gtf_lines(cands, method):
    out = []
    for i, (chrom, strand, start, end, chain) in enumerate(cands):
        exons = _exons_from_chain(start, end, chain)
        if not exons:
            continue
        tid = f"{method}_{i}"
        attr = f'gene_id "{tid}"; transcript_id "{tid}";'
        out.append(f"{chrom}\tpyfin\ttranscript\t{start+1}\t{end}\t.\t{strand}\t.\t{attr}\n")
        for es, ee in exons:
            out.append(f"{chrom}\tpyfin\texon\t{es+1}\t{ee}\t.\t{strand}\t.\t{attr}\n")
    return out


SUF = os.environ.get("OUT_SUFFIX", "")
for method, cands in [("NEW_nogate" + SUF, nogate), ("NEW_m1tie" + SUF, m1tie)]:
    with open(f"{OUT}/{method}.gtf", "w") as fh:
        fh.writelines(gtf_lines(cands, method))
    sys.stderr.write(f"{method}: {len(cands)} candidates\n")
sys.stderr.write(f"GEN DONE (M1_MAX_INDEL_BP={os.environ.get('M1_MAX_INDEL_BP','50')})\n")
