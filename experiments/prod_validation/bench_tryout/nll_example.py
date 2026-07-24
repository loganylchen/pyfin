#!/usr/bin/env python3
"""Dump CONCRETE examples of a tied read whose per-candidate junction-window NLL is
nearly identical (margin ~ 0) vs clearly separated (margin >= 0.5), so we can judge
whether the signal genuinely can't tell wobble siblings apart. For each tie read that
eventalign can score on >=2 candidates, print: read id, the discrimination window
(genomic span), and for every tied candidate its intron chain, mean NLL, and #events.
Early-exits once enough examples of each kind are collected.
"""
import os, sys
import numpy as np
import mappy, pysam

REPO = "/SSD/logan/dev/pyfin"
sys.path.insert(0, REPO)
from fin.io.interval_manager import generate_isolated_intervals
from fin.io.io_gtf import GTFReader
from fin.candidates.discovery import discover_candidates
from fin.scoring.mappy_preset import get_m1_preset
from fin.scoring.mappy_score import score_hit
from fin.scoring.m2_junction_nll import class_junction_window_set, read_cand_mean_nll
from fin.scoring.krill_aligner import make_krill_aligner, krill_thread_count
import krill

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S = "SGNex_H9_directRNA_replicate2_run2"
BAM = f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
GENOME = f"{B}/refs/host/GRCh38.primary_assembly.genome.fa"
GTF = f"{B}/results/gencode_full_sweep/_ref/p00/annotation.gtf"
SIGNAL = "/autofs/NAS25_Shared/public_data/NanoporeDRS/human/sgnexdata/data/H9directRNAreplicate2run2/blow5/nanopore.blow5"
LIMIT = int(os.environ.get("LIMIT", "800"))
USE_GPU = os.environ.get("USE_GPU", "1") == "1"
N_NEAR = int(os.environ.get("N_NEAR", "12"))     # margin < 0.05 examples wanted
N_CLEAR = int(os.environ.get("N_CLEAR", "6"))    # margin >= 0.5 examples wanted

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


def fmt_chain(chain):
    return " ".join(f"({s},{e})" for s, e in chain)


def diff_junctions(chains):
    """junction coords that are NOT shared by all candidates (the discriminators)."""
    common = set(chains[0])
    for c in chains[1:]:
        common &= set(c)
    return [sorted(set(c) - common) for c in chains]


def emit(iv, rid, gset, tie_cands, tie_as, scored):
    sj = sorted(gset)
    span = f"{iv.chrom}:{sj[0]}-{sj[-1]}" if sj else "(empty)"
    chains = [tuple(c.intron_chain.introns) for c in tie_cands]
    diffs = diff_junctions(chains)
    nll_by_idx = dict(scored)
    print(f"\n--- read {rid}  window {span} ({len(gset)} gpos)  tie={len(tie_cands)} ---")
    order = sorted(range(len(tie_cands)), key=lambda i: nll_by_idx.get(i, 9e9))
    for rank, i in enumerate(order):
        c = tie_cands[i]
        nll = nll_by_idx.get(i)
        tag = "BEST" if rank == 0 else f"#{rank+1}"
        nlls = f"{nll:.4f}" if nll is not None else "  n/a "
        print(f"  {tag:5} NLL={nlls}  AS={tie_as[i]:.0f}  "
              f"diffJ={fmt_chain(diffs[i]) or '(none)'}")
        print(f"        chain: {fmt_chain(chains[i])}")


near, clear = [], []
for n, iv in enumerate(intervals):
    if len(near) >= N_NEAR and len(clear) >= N_CLEAR:
        break
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
        tie_cands = [cand_list[j] for j in tie]
        tie_as = [row[j] for j in tie]
        tie_aln = [aligners[j] for j in tie]
        gset = class_junction_window_set(tie_cands, flank=2, k=10)
        if not gset:
            continue
        scored = []
        for idx, cand in enumerate(tie_cands):
            nll, nev = read_cand_mean_nll(rid, seq, cand, [], krill_aligner,
                                          tie_aln[idx], SIGNAL, "rna002", gset=gset,
                                          use_gpu=eff_gpu, num_thread=nthr)
            if nev > 0 and np.isfinite(nll):
                scored.append((idx, float(nll)))
        if len(scored) < 2:
            continue
        scored.sort(key=lambda t: t[1])
        margin = scored[1][1] - scored[0][1]
        if margin < 0.05 and len(near) < N_NEAR:
            near.append(1); print(f"\n[NEAR-TIE margin={margin:.4f}]", end="")
            emit(iv, rid, gset, tie_cands, tie_as, scored)
        elif margin >= 0.5 and len(clear) < N_CLEAR:
            clear.append(1); print(f"\n[CLEAR   margin={margin:.4f}]", end="")
            emit(iv, rid, gset, tie_cands, tie_as, scored)
        if len(near) >= N_NEAR and len(clear) >= N_CLEAR:
            break
    if (n + 1) % 50 == 0:
        sys.stderr.write(f"  {n+1}/{len(intervals)} near={len(near)} clear={len(clear)}\n")

sys.stderr.write(f"DONE near={len(near)} clear={len(clear)}\n")
