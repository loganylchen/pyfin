"""How informative is krill's polyA call on these datasets? Evaluation-only.
For each sample: fraction of reads with a PASS call, and for reads assigned to
truth transcripts, how the read's genomic 3' end sits relative to the annotated
TES (is the polyA read a reliable TES witness?)."""
import sys, collections, gzip, re, statistics as st
sys.path.insert(0, "/SSD/logan/dev/pyfin")
from fin.scoring.polya import compute_polya
import pysam

ds, sample, bam, fq, blow5, truth = sys.argv[1:7]
seqs = {}
with gzip.open(fq, "rt") as fh:
    while True:
        h = fh.readline()
        if not h: break
        s = fh.readline().rstrip(); fh.readline(); fh.readline()
        seqs[h[1:].split()[0]] = s
pm = compute_polya(seqs, blow5, pore="rna002", use_gpu=False)
n = len(seqs); scored = len(pm)
passed = {r for r, (L, qc) in pm.items() if qc == "PASS" and L is not None and L > 0}
print(f"{ds} {sample}: reads={n} scored={scored} PASS={len(passed)} ({100*len(passed)/max(n,1):.1f}%)")
lens = sorted(L for r, (L, qc) in pm.items() if r in passed)
if lens: print(f"  polyA length PASS: q10 {lens[len(lens)//10]:.0f} median {lens[len(lens)//2]:.0f} q90 {lens[9*len(lens)//10]:.0f}")
# TES witness quality: read 3' end vs annotated TES of the transcript it overlaps best
tes = {}
for l in open(truth):
    if l.startswith("#"): continue
    f = l.split("\t")
    if len(f) < 9 or f[2] != "transcript": continue
    t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    tes[t] = (f[0], int(f[3]), int(f[4]), f[6])
by_chrom = collections.defaultdict(list)
for t, (c, s, e, strand) in tes.items(): by_chrom[c].append((s, e, strand, t))
d_pass, d_fail = [], []
with pysam.AlignmentFile(bam) as b:
    for r in b.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        cands = [(min(r.reference_end, e) - max(r.reference_start, s), s, e, strand)
                 for (s, e, strand, t) in by_chrom.get(r.reference_name, [])]
        if not cands: continue
        ov, s, e, strand = max(cands)
        if ov <= 0: continue
        read3 = r.reference_end if strand == "+" else r.reference_start
        anno3 = e if strand == "+" else s
        d = (read3 - anno3) if strand == "+" else (anno3 - read3)   # >0 = read runs PAST the TES
        (d_pass if r.query_name in passed else d_fail).append(d)
def summ(v):
    if not v: return "n=0"
    v = sorted(v); n = len(v)
    within = sum(1 for x in v if abs(x) <= 25)
    return f"n={n} median {v[n//2]:+.0f} |d|<=25bp {100*within/n:.0f}%  readthrough>50bp {100*sum(1 for x in v if x>50)/n:.0f}%  short>50bp {100*sum(1 for x in v if x<-50)/n:.0f}%"
print("  read 3' end vs annotated TES, polyA PASS reads :", summ(d_pass))
print("  read 3' end vs annotated TES, other reads      :", summ(d_fail))
