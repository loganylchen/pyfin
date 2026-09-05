"""What does the 3'-endpoint distribution of polyA PASS reads look like around
the NEAREST annotated TES on the same strand (not the best-overlap transcript)?
This is the empirical g_s / b_s the endpoint model must fit within-sample."""
import sys, collections, gzip, re
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
passed = {r for r, (L, qc) in pm.items() if qc == "PASS" and L and L > 0}
tes = collections.defaultdict(list)   # (chrom,strand) -> [tes positions]
for l in open(truth):
    if l.startswith("#"): continue
    f = l.split("\t")
    if len(f) < 9 or f[2] != "transcript": continue
    tes[(f[0], f[6])].append(int(f[4]) if f[6] == "+" else int(f[3]))
d_pass, d_other = [], []
with pysam.AlignmentFile(bam) as b:
    for r in b.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        strand = "-" if r.is_reverse else "+"
        ends = tes.get((r.reference_name, strand))
        if not ends: continue
        read3 = r.reference_end if strand == "+" else r.reference_start
        # signed distance to the nearest annotated TES: >0 = past it (readthrough)
        d = min(((read3 - t) if strand == "+" else (t - read3) for t in ends), key=abs)
        (d_pass if r.query_name in passed else d_other).append(d)
def hist(v, name):
    v = sorted(v); n = len(v)
    bins = [(-10**9,-100),(-100,-50),(-50,-25),(-25,-10),(-10,-3),(-3,3),(3,10),(10,25),(25,50),(50,100),(100,10**9)]
    print(f"  {name} n={n}:", " ".join(f"[{lo},{hi}):{100*sum(1 for x in v if lo<=x<hi)/n:.1f}%" for lo,hi in bins))
    a = [abs(x) for x in v]; a.sort()
    print(f"     |d| median {a[n//2]}  q75 {a[3*n//4]}  q90 {a[9*n//10]}  MAD-of-d {sorted(abs(x - v[n//2]) for x in v)[n//2]}")
print(f"{ds} {sample}  PASS={len(passed)}/{len(seqs)}")
hist(d_pass, "polyA PASS: read3' - nearest same-strand TES")
hist(d_other, "other      : read3' - nearest same-strand TES")
