#!/usr/bin/env python3
"""Compare candidate GENERATION only (no gates, no EM): OLD production clustering
(canonical-search + group-by-3'+exact-chain + collapse) vs NEW chain-cluster
(group-by-chain, collapse EXACT sub-chains, cluster wobble/cassette/containment keep
members, NO canonical search). Candidates ARE the final result -> gffcompare =/j/c.

Same intervals + reads for both (generate_isolated_intervals); NOVEL candidates only.
Writes two GTFs and runs gffcompare vs gencode truth.
"""
import os, sys, subprocess, re
import pysam

REPO = "/SSD/logan/dev/pyfin"
sys.path.insert(0, REPO)
from fin.io.interval_manager import generate_isolated_intervals, is_fusion_read
from fin.io.io_gtf import GTFReader
from fin.io.io_bam import BamReader
from fin.candidates.discovery import discover_candidates, _exons_from_chain
from fin.candidates.intron_chains import extract_intron_chain
from fin.candidates.dataclasses import IntronChain
from fin.candidates.chain_cluster import cluster_read_chains

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S = "SGNex_H9_directRNA_replicate2_run2"
BAM = f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
GENOME = f"{B}/refs/host/GRCh38.primary_assembly.genome.fa"
GTF = f"{B}/results/gencode_full_sweep/_ref/p00/annotation.gtf"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
OUT = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/cluster_cmp"
IMG = "quay.io/biocontainers/gffcompare:0.12.10--h9948957_0"
os.makedirs(OUT, exist_ok=True)

fa = pysam.FastaFile(GENOME)
_chrom_cache = {}
def chrom_seq(c):
    if c not in _chrom_cache:
        try: _chrom_cache[c] = fa.fetch(c)
        except Exception: _chrom_cache[c] = ""
    return _chrom_cache[c]

LIMIT = int(os.environ.get("LIMIT", "0"))
gtf_reader = GTFReader(GTF); gtf_reader.open(); gtf_reader.parse()
intervals = generate_isolated_intervals(BAM, gtf_path=GTF, max_gap=0)["intervals"]
if LIMIT:
    intervals = intervals[:LIMIT]
sys.stderr.write(f"intervals: {len(intervals)}\n")


def fetch_read_chains(interval):
    """Replicate discovery's read fetch: strand-filter, skip fusion; return
    [(read_dict, IntronChain)] and {query_name: (ref_start, ref_end)}."""
    rc, spans = [], {}
    with BamReader(BAM) as bam:
        for aln in bam.fetch(reference=interval.chrom, start=interval.start, end=interval.end):
            rd = bam.alignment_to_dict(aln)
            if not rd or not rd.get("is_mapped"):
                continue
            if interval.strand is not None:
                if ("+" if rd.get("is_forward") else "-") != interval.strand:
                    continue
            if is_fusion_read(rd):
                continue
            cig = rd.get("cigartuples"); rs = rd.get("reference_start")
            if cig is None or rs is None:
                continue
            rc.append((rd, extract_intron_chain(cig, rs)))
            q = rd.get("query_name")
            if q:
                spans[q] = (rs, rd.get("reference_end", rs))
    return rc, spans


def gtf_lines(cands, method):
    """cands: list of (chrom, strand, start, end, IntronChain). Yield GTF lines."""
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


old_cands, new_cands = [], []
for n, iv in enumerate(intervals):
    cs = chrom_seq(iv.chrom)
    # OLD: production discovery (canonical_search=4), novel only.
    try:
        cset = discover_candidates(
            iv, BAM, gtf_reader, cs, threshold=24, min_novel_reads=1,
            canonical_search_bp=4, max_chains_per_read=16,
            canonical_motifs=("GT-AG", "GC-AG", "AT-AC"))
        for c in cset.candidates:
            if c.source == "novel":
                old_cands.append((c.chrom, c.strand, c.start, c.end, c.intron_chain))
    except Exception as e:
        sys.stderr.write(f"OLD fail {iv.region_string}: {e}\n")
    # NEW: chain-cluster (no canonical search), all members; span = union of reads.
    try:
        rc, spans = fetch_read_chains(iv)
        strand = iv.strand or "+"
        for cl in cluster_read_chains(rc):
            for m in cl.members:
                rs = [spans[r][0] for r in m.read_ids if r in spans]
                re_ = [spans[r][1] for r in m.read_ids if r in spans]
                if not rs:
                    continue
                new_cands.append((iv.chrom, strand, min(rs), max(re_), m.chain))
    except Exception as e:
        sys.stderr.write(f"NEW fail {iv.region_string}: {e}\n")
    if (n + 1) % 2000 == 0:
        sys.stderr.write(f"  {n+1}/{len(intervals)} old={len(old_cands)} new={len(new_cands)}\n")

for method, cands in [("OLD", old_cands), ("NEW", new_cands)]:
    p = f"{OUT}/{method}.gtf"
    with open(p, "w") as fh:
        fh.writelines(gtf_lines(cands, method))
    sys.stderr.write(f"{method}: {len(cands)} candidates -> {p}\n")
sys.stderr.write("GEN DONE — now score on host: gffcompare OLD.gtf/NEW.gtf vs truth\n")
