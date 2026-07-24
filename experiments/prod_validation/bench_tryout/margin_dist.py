#!/usr/bin/env python3
"""For NEW chain-cluster candidates' TIED reads, compute the M2 margin
(2nd-best junction-NLL - best junction-NLL, over the distinguishing window) via
krill eventalign, and report its distribution. Answers: how separated are the
best vs 2nd-best candidate for a typical tied read -> what margin threshold means
"distinguishable". Runs on LIMIT intervals for a representative sample.
"""
import os, sys
from collections import Counter
import numpy as np
import mappy, pysam

REPO = "/SSD/logan/dev/pyfin"
sys.path.insert(0, REPO)
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_gtf import GTFReader
from fin.candidates.discovery import discover_candidates
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit
from fin.scoring.m2_junction_nll import m2_resolve_tie
from fin.scoring.krill_aligner import make_krill_aligner, krill_thread_count
import krill

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S = "SGNex_H9_directRNA_replicate2_run2"
BAM = f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
GENOME = f"{B}/refs/host/GRCh38.primary_assembly.genome.fa"
GTF = f"{B}/results/gencode_full_sweep/_ref/p00/annotation.gtf"
SIGNAL = "/autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/blow5/nanopore.blow5"
LIMIT = int(os.environ.get("LIMIT", "3000"))
USE_GPU = os.environ.get("USE_GPU", "1") == "1"

fa = pysam.FastaFile(GENOME); _cc = {}
def chrom_seq(c):
    if c not in _cc:
        try: _cc[c] = fa.fetch(c)
        except Exception: _cc[c] = ""
    return _cc[c]

gtf_reader = GTFReader(GTF); gtf_reader.open(); gtf_reader.parse()
intervals = generate_isolated_intervals(BAM, gtf_path=GTF, max_gap=0)["intervals"]
if LIMIT: intervals = intervals[:LIMIT]
sys.stderr.write(f"intervals: {len(intervals)}  gpu={USE_GPU}\n")
preset = get_m1_preset()
nthr = krill_thread_count()
krill_aligner, eff_gpu = make_krill_aligner(krill, "rna002", USE_GPU,
                                            hmm_confidence=False, num_thread=nthr)

margins = []          # finite margins (best decidable)
n_tied = 0            # tied reads seen
n_scored = 0          # tied reads where >=2 candidates got a window NLL
n_nowindow = 0        # tied reads with no discrimination window / <2 scored
for n, iv in enumerate(intervals):
    cs = chrom_seq(iv.chrom)
    cset = discover_candidates(iv, BAM, gtf_reader, cs, threshold=24,
                               min_novel_reads=1, chain_cluster=True, canonical_search_bp=0)
    cand_list = list(cset.candidates)
    if len(cand_list) < 2:
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
        if best <= 0: continue
        tie = [j for j in range(len(cand_list)) if row[j] > 0 and row[j] >= best - 1e-9]
        if len(tie) < 2: continue
        n_tied += 1
        best_idx, margin, scored = m2_resolve_tie(
            rid, seq, [cand_list[j] for j in tie], SIGNAL, pore="rna002",
            krill_aligner=krill_aligner, mappy_aligners=[aligners[j] for j in tie],
            return_scored=True, use_gpu=eff_gpu, num_thread=nthr)
        if best_idx is None or len(scored) < 2:
            n_nowindow += 1
            continue
        n_scored += 1
        if margin != float("inf") and np.isfinite(margin):
            margins.append(float(margin))
    if (n + 1) % 200 == 0 and margins:
        mm = np.array(margins)
        sys.stderr.write(
            f"  [{n+1}/{len(intervals)}] tied={n_tied} scored={n_scored} "
            f"nowin={n_nowindow}({100*n_nowindow/n_tied:.0f}%) | "
            f"margin median={np.median(mm):.3f} "
            f">=0.1:{100*np.mean(mm>=0.1):.0f}% >=0.5:{100*np.mean(mm>=0.5):.0f}% "
            f">=1:{100*np.mean(mm>=1):.0f}%\n")
        sys.stderr.flush()

m = np.array(margins) if margins else np.array([0.0])
print(f"\n=== M2 margin (2nd-best NLL - best NLL) on TIED reads ===")
print(f"tied reads: {n_tied}   scored (>=2 window NLL): {n_scored}   "
      f"no-window/<2-scored: {n_nowindow} ({100*n_nowindow/n_tied if n_tied else 0:.0f}%)")
print(f"finite margins collected: {len(margins)}")
print(f"  median={np.median(m):.3f}  mean={np.mean(m):.3f}  "
      f"p25={np.percentile(m,25):.3f}  p75={np.percentile(m,75):.3f}  max={m.max():.2f}")
for thr in (0.05, 0.1, 0.2, 0.5, 1.0, 2.0):
    frac = 100 * np.mean(m >= thr)
    print(f"  fraction of scored ties with margin >= {thr:<4}: {frac:.0f}%")
