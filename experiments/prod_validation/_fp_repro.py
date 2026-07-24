#!/usr/bin/env python3
"""Are FPs systematic (tool-intrinsic) or stochastic (per-run noise)?

For each dataset/tool at the FULL clean annotation:
  - FP = a predicted transcript whose gffcompare class_code != '='
         (i.e. it does NOT exactly reproduce a truth transcript).
  - Signature = structural fingerprint independent of transcript_id:
        multi-exon: (chrom, strand, tuple(intron coords))
        mono-exon : (chrom, strand, round(start,-2), round(end,-2))
  - Across the N replicates of that dataset, we ask how often the SAME FP
    signature recurs.

Metrics per tool:
  mean_fp        = mean # FP signatures per replicate
  recur_ge2_frac = fraction of distinct FP signatures seen in >=2 replicates
  shared_all     = # FP signatures present in ALL replicates
  mean_jaccard   = mean pairwise Jaccard of FP signature sets across replicates

High recurrence => FP is a deterministic property of the tool+data (systematic).
Low recurrence  => FP is run-to-run stochastic noise.
"""
import os, glob, re, itertools
from collections import defaultdict, Counter

DS_ROOT = "/autofs/mnemosyne4_SSD/logan/dev/pyfin/experiments/prod_validation"
DATASETS = ["sirv4", "heya8", "sequin"]
TOOLS = ["pyfin_prod", "bambu", "espresso", "flair", "isoquant",
         "isotools", "lafite", "stringtie3", "talon"]
RATIO = "full"


def gtf_path(dsdir, s, tool):
    if tool == "pyfin_prod":
        p = os.path.join(dsdir, s, RATIO, "pyfin.gtf")
        if os.path.exists(p):
            return p
    return os.path.join(dsdir, s, RATIO, "scoring", f"gc_{tool}.input.gtf")


def class_codes(dsdir, s, tool):
    """qry_id -> class_code from the tmap."""
    g = glob.glob(os.path.join(dsdir, s, RATIO, "scoring", f"gc_{tool}.*.tmap"))
    cc = {}
    if not g:
        return cc
    with open(g[0]) as fh:
        h = next(fh).rstrip("\n").split("\t")
        ci = {x: i for i, x in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            cc[f[ci["qry_id"]]] = f[ci["class_code"]]
    return cc


def signatures(gtf, cc):
    """return set of structural signatures for FP transcripts (class != '=')."""
    exons = defaultdict(list)   # tid -> [(start,end)]
    meta = {}                   # tid -> (chrom, strand)
    if not os.path.exists(gtf):
        return set(), 0
    for ln in open(gtf):
        if ln.startswith("#"):
            continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if not m:
            continue
        tid = m.group(1)
        exons[tid].append((int(p[3]), int(p[4])))
        meta[tid] = (p[0], p[6])
    sigs = set()
    ntx = 0
    for tid, ex in exons.items():
        ntx += 1
        code = cc.get(tid, "u")   # unmatched if absent from tmap
        if code == "=":
            continue              # true positive, skip
        ex.sort()
        chrom, strand = meta[tid]
        if len(ex) == 1:
            sig = (chrom, strand, round(ex[0][0], -2), round(ex[0][1], -2))
        else:
            introns = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
            sig = (chrom, strand, introns)
        sigs.add(sig)
    return sigs, ntx


print(f"{'dataset':8s} {'tool':11s} | {'mean_fp':>7s} {'recur>=2':>8s} {'shared_all':>10s} {'jaccard':>7s}")
print("-" * 64)
rows = []
for ds in DATASETS:
    dsdir = os.path.join(DS_ROOT, ds)
    samples = sorted(os.path.basename(os.path.dirname(p))
                     for p in glob.glob(os.path.join(dsdir, "*", "stage")))
    for tool in TOOLS:
        persample = []
        for s in samples:
            cc = class_codes(dsdir, s, tool)
            sigs, _ = signatures(gtf_path(dsdir, s, tool), cc)
            persample.append(sigs)
        persample = [x for x in persample if x is not None]
        nrep = len(persample)
        if nrep < 2:
            continue
        counter = Counter()
        for st in persample:
            for sig in st:
                counter[sig] += 1
        distinct = len(counter)
        mean_fp = sum(len(st) for st in persample) / nrep
        recur2 = sum(1 for c in counter.values() if c >= 2)
        recur2_frac = 100 * recur2 / distinct if distinct else 0.0
        shared_all = sum(1 for c in counter.values() if c == nrep)
        # mean pairwise jaccard
        js = []
        for a, b in itertools.combinations(persample, 2):
            u = len(a | b)
            js.append(len(a & b) / u if u else 1.0)
        mj = 100 * sum(js) / len(js) if js else 0.0
        rows.append((ds, tool, mean_fp, recur2_frac, shared_all, mj))
        print(f"{ds:8s} {tool:11s} | {mean_fp:7.1f} {recur2_frac:7.1f}% {shared_all:10d} {mj:6.1f}%")
    print()

out = os.path.join(DS_ROOT, "_fp_repro.tsv")
with open(out, "w") as fh:
    fh.write("dataset\ttool\tmean_fp\trecur_ge2_pct\tshared_all\tmean_jaccard_pct\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in
                 (r[0], r[1], f"{r[2]:.2f}", f"{r[3]:.2f}", r[4], f"{r[5]:.2f}")) + "\n")
print("wrote", out)
