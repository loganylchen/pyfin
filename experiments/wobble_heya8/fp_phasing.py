#!/usr/bin/env python3
"""For each FP candidate chain, measure read-level phasing support from the BAM.

full_support      = reads whose intron chain *contains* the FP chain as a contiguous
                    run (all junctions within tol) AND span 5'..3' of it
junction_support  = per-junction #reads carrying that exact junction
pair_support      = per-adjacent-junction-PAIR #reads carrying both consecutive junctions
                    -> a pair with 0 reads means the chain is an unobserved recombination

Also lists overlapping truth transcripts and their full-support, for contrast.
"""
import sys, re
import pysam

BAM = sys.argv[1]; FPGTF = sys.argv[2]; TRUTH = sys.argv[3]
FP_IDS = sys.argv[4].split(",")
TOL = 8

def gtf_chains(path):
    tx = {}
    for line in open(path):
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "exon": continue
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if not m: continue
        d = tx.setdefault(m.group(1), {"chrom": p[0], "strand": p[6], "exons": []})
        d["exons"].append((int(p[3])-1, int(p[4])))
    out = {}
    for t, d in tx.items():
        ex = sorted(d["exons"])
        intr = tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))
        out[t] = {"chrom": d["chrom"], "strand": d["strand"], "introns": intr,
                  "start": ex[0][0], "end": ex[-1][1]}
    return out

fp = gtf_chains(FPGTF)
truth = gtf_chains(TRUTH)
bam = pysam.AlignmentFile(BAM, "rb")

def read_introns(aln):
    intr = []; pos = aln.reference_start
    for op, ln in aln.cigartuples or []:
        if op in (0,7,8,2): pos += ln
        elif op == 3: intr.append((pos, pos+ln)); pos += ln
    return intr, aln.reference_start, pos

def jeq(a, b): return abs(a[0]-b[0]) <= TOL and abs(a[1]-b[1]) <= TOL

def analyze(name, g, chrom):
    introns = g["introns"]
    njun = len(introns)
    jcount = [0]*njun
    paircount = [0]*max(njun-1,1)
    full = 0
    total_overlap = 0
    for aln in bam.fetch(chrom, g["start"], g["end"]):
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary: continue
        rintr, rs, re_ = read_introns(aln)
        if not rintr: continue
        total_overlap += 1
        # which of the candidate junctions does this read carry?
        carried = [any(jeq(rj, cj) for rj in rintr) for cj in introns]
        for i, c in enumerate(carried):
            if c: jcount[i] += 1
        for i in range(njun-1):
            if carried[i] and carried[i+1]: paircount[i] += 1
        # full support: read spans whole candidate & carries every junction
        if all(carried) and rs <= g["start"]+TOL and re_ >= g["end"]-TOL:
            full += 1
    return njun, full, jcount, paircount, total_overlap

print(f"{'cand':<16}{'njun':>5}{'fullchain':>10}{'minJun':>8}{'minPair':>8}{'overlapRds':>11}")
for fid in FP_IDS:
    if fid not in fp:
        print(f"{fid}: not found"); continue
    g = fp[fid]; chrom = g["chrom"]
    njun, full, jc, pc, tot = analyze(fid, g, chrom)
    print(f"{fid:<16}{njun:>5}{full:>10}{min(jc) if jc else 0:>8}{min(pc) if pc else 0:>8}{tot:>11}")
    print(f"    per-junction support: {jc}")
    print(f"    per-pair    support: {pc}")
    # nearest truth siblings (same chrom, overlapping)
    sibs = [(t, td) for t, td in truth.items() if td["chrom"]==chrom
            and not (td["end"] < g["start"] or td["start"] > g["end"])]
    for t, td in sorted(sibs, key=lambda x:-len(x[1]["introns"]))[:3]:
        nj, fu, jcs, pcs, to = analyze(t, td, chrom)
        print(f"    TRUTH {t:<12} njun={nj} fullchain={fu} minPair={min(pcs) if pcs else 0}")
