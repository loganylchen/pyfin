"""Honest Sn / Pr / F1 that the annotation-passthrough trick cannot fool.

Expressed-truth (the recall denominator) = TRUTH transcripts that are actually
expressed in this sample. Preferred measure: **nanocount est_count** (reads
aligned to the truth transcriptome, annotation-independent; the sweep ships a
per-condition nanocount.tsv). Fallback: chain-matching read count from the BAM.

At expression threshold k:
  expr_k  = TRUTH transcripts with est_count (or chain-reads) >= k
  TP_k    = expr_k transcripts the tool recovered (gffcompare class '=' in refmap)
  Sn_k    = TP_k / |expr_k|                 (recall over EXPRESSED truth only)
  Pr_k    = TP_k / (tool's TOTAL reported transcripts)
            -> every phantom (0-read annotation copy) and shadow counts against Pr.
  F1_k    = 2*Sn*Pr/(Sn+Pr)

Args: <workdir> <truth.gtf> <bam> [nanocount.tsv]
  4th arg present -> expressed-truth from nanocount est_count (recommended).
  4th arg absent  -> legacy chain-match read-count proxy from the BAM.
workdir holds <tool>.gtf and gc_<tool>.<tool>.gtf.refmap per tool.
"""
import sys, re, glob, os
import pysam
from collections import defaultdict

W, TRUTH, BAM = sys.argv[1], sys.argv[2], sys.argv[3]
NANOCOUNT = sys.argv[4] if len(sys.argv) > 4 else None

def chains(path):
    tx = {}
    for line in open(path):
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "exon": continue
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if not m: continue
        t = tx.setdefault(m.group(1), {"chrom": p[0], "strand": p[6], "exons": []})
        t["exons"].append((int(p[3]) - 1, int(p[4])))
    out = []
    for tid, t in tx.items():
        ex = sorted(t["exons"])
        intr = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        out.append({"tid": tid, "chrom": t["chrom"], "strand": t["strand"], "chain": intr,
                    "3p": ex[-1][1] if t["strand"] == "+" else ex[0][0],
                    "start": ex[0][0], "end": ex[-1][1]})
    return out

def rchain(aln):
    intr = []; pos = aln.reference_start
    for op, ln in aln.cigartuples:
        if op in (0, 7, 8, 2): pos += ln
        elif op == 3: intr.append((pos, pos + ln)); pos += ln
    return tuple(intr), aln.reference_start, pos

def close(c1, c2, tol=6):
    if len(c1) != len(c2): return False
    return all(abs(a - c) <= tol and abs(b - d) <= tol for (a, b), (c, d) in zip(c1, c2))

bam = pysam.AlignmentFile(BAM, "rb")
def support(g):
    seen = set()
    try: it = bam.fetch(g["chrom"], g["start"], g["end"])
    except ValueError: return 0
    for aln in it:
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary: continue
        ch, rs, re_ = rchain(aln)
        r3 = re_ if g["strand"] == "+" else rs
        if g["chain"]:
            if close(ch, g["chain"]) and abs(r3 - g["3p"]) <= 50: seen.add(aln.query_name)
        elif len(ch) == 0 and abs(r3 - g["3p"]) <= 50 and rs >= g["start"] - 100 and re_ <= g["end"] + 100:
            seen.add(aln.query_name)
    return len(seen)

gt = chains(TRUTH)
gt_tids = {g["tid"] for g in gt}

if NANOCOUNT:
    # Expressed-truth from nanocount est_count (annotation-independent). Truth
    # transcripts absent from nanocount.tsv have zero assigned reads -> est_count 0.
    est = {t: 0.0 for t in gt_tids}
    fh = open(NANOCOUNT); next(fh, None)
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if len(p) >= 3 and p[0] in est:
            try: est[p[0]] = float(p[2])  # est_count column
            except ValueError: pass
    expr1 = {t for t, v in est.items() if v >= 1}
    expr3 = {t for t, v in est.items() if v >= 3}
    EXPR_SRC = "nanocount est_count"
else:
    sup = {g["tid"]: support(g) for g in gt}
    expr1 = {t for t, n in sup.items() if n >= 1}
    expr3 = {t for t, n in sup.items() if n >= 3}
    EXPR_SRC = "chain-match reads"

def matched(refmap):
    out = set()
    with open(refmap) as fh:
        next(fh, None)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2] == "=": out.add(f[1])
    return out

tools = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(W, "*.gtf"))
               if os.path.basename(p) != "truth.gtf" and not os.path.basename(p).startswith("gc_"))
print(f"truth={len(gt)}  expressed>=1={len(expr1)}  expressed>=3={len(expr3)}  "
      f"[{EXPR_SRC}]\n")
hdr = f"{'tool':<12}|{'Sn@1':>6}{'Pr@1':>6}{'F1@1':>6} |{'Sn@3':>6}{'Pr@3':>6}{'F1@3':>6} |{'output':>7}"
print(hdr); print("-" * len(hdr))
def f1(s, p): return 2 * s * p / (s + p) if s + p else 0.0
res = []
for t in tools:
    rm = glob.glob(os.path.join(W, f"gc_{t}.*.refmap"))
    if not rm: continue
    mt = matched(rm[0])
    nout = len(chains(os.path.join(W, f"{t}.gtf")))
    m1, m3 = len(mt & expr1), len(mt & expr3)
    sn1, pr1 = 100 * m1 / max(len(expr1), 1), 100 * m1 / max(nout, 1)
    sn3, pr3 = 100 * m3 / max(len(expr3), 1), 100 * m3 / max(nout, 1)
    res.append((t, sn1, pr1, f1(sn1, pr1), sn3, pr3, f1(sn3, pr3), nout))
for r in sorted(res, key=lambda x: -x[6]):
    print(f"{r[0]:<12}|{r[1]:>6.1f}{r[2]:>6.1f}{r[3]:>6.1f} |{r[4]:>6.1f}{r[5]:>6.1f}{r[6]:>6.1f} |{r[7]:>7}")
