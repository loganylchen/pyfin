#!/usr/bin/env python3
"""On NEW chain-cluster candidates + M1 (mappy AS), measure the per-read best-AS tie
set size: size 1 = unique best (assigned 100%, no M2 needed); size >=2 = ambiguous
(needs M2 to discriminate wobble siblings). Reports the distribution and, for the tied
reads, how many are within a signal-resolvable wobble cluster.
"""
import os, sys
from collections import Counter
import mappy, pysam

REPO = "/SSD/logan/dev/pyfin"
sys.path.insert(0, REPO)
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_gtf import GTFReader
from fin.candidates.discovery import discover_candidates
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S = "SGNex_H9_directRNA_replicate2_run2"
BAM = f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
GENOME = f"{B}/refs/host/GRCh38.primary_assembly.genome.fa"
GTF = f"{B}/results/gencode_full_sweep/_ref/p00/annotation.gtf"
LIMIT = int(os.environ.get("LIMIT", "0"))

fa = pysam.FastaFile(GENOME); _cc = {}
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

n_reads = 0        # reads with >=1 candidate at AS>0 (kept by M1)
n_noalign = 0      # reads with no candidate at AS>0 (dropped by M1 keep)
tie_hist = Counter()
for n, iv in enumerate(intervals):
    cs = chrom_seq(iv.chrom)
    cset = discover_candidates(iv, BAM, gtf_reader, cs, threshold=24,
                               min_novel_reads=1, chain_cluster=True, canonical_search_bp=0)
    cand_list = list(cset.candidates)
    if not cand_list:
        continue
    read_seqs = getattr(cset, "read_sequences", {}) or {}
    aligners = [mappy.Aligner(seq=c.sequence, preset=preset) if c.sequence else None
                for c in cand_list]
    for rid, seq in read_seqs.items():
        if not seq:
            continue
        best = -1e18; row = [-1e18] * len(cand_list)
        for j, aln in enumerate(aligners):
            if aln is None: continue
            b = None
            for h in aln.map(seq):
                v = score_hit(h)
                if v is not None and (b is None or v > b): b = v
            if b is not None:
                row[j] = b
                if b > best: best = b
        if best <= 0:
            n_noalign += 1
            continue
        tie = sum(1 for x in row if x > 0 and x >= best - 1e-9)
        n_reads += 1
        tie_hist[tie] += 1
    if (n + 1) % 2000 == 0:
        sys.stderr.write(f"  {n+1}/{len(intervals)} reads={n_reads}\n")

uniq = tie_hist.get(1, 0)
tied = n_reads - uniq
print(f"\n=== per-read best-AS tie size (NEW candidates) ===")
print(f"reads kept by M1 (AS>0): {n_reads}   dropped (no align): {n_noalign}")
print(f"  UNIQUE best (tie=1, 100% assigned, NO M2): {uniq}  ({100*uniq/n_reads if n_reads else 0:.1f}%)")
print(f"  TIED (>=2, needs M2):                      {tied}  ({100*tied/n_reads if n_reads else 0:.1f}%)")
print(f"\ntie-size histogram:")
for k in sorted(tie_hist):
    print(f"  tie={k:2}: {tie_hist[k]:8} ({100*tie_hist[k]/n_reads:.1f}%)")
