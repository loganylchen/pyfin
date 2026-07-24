#!/usr/bin/env python3
"""Test the 5'-end-distribution hypothesis:

A short isoform that GRAPH assembly destroyed (matched '=' in OFF p00 but NOT in
GRAPH p00) which is 5'-internal to a surviving LONGER transcript is the exact
degradation-lookalike case. If the short isoform is REAL (its own TSS), reads
should show a PEAK of 5'-ends at its annotated TSS, rising above the smooth
degradation background. If it were pure degradation of the long transcript, the
5'-ends would be a smooth ramp with no peak there.

We measure, per such locus, the enrichment of read-5'-ends in a +/-TOL window at
the lost isoform's TSS vs the smooth local background (median per-window count
across the transcript body). Enrichment >> 1 across loci confirms the signal.
"""
import os, re, sys, statistics
import pysam

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
GTF = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
BAM = (f"{B}/results/gencode_full_sweep/SGNex_H9_directRNA_replicate2_run2/"
       "align/SGNex_H9_directRNA_replicate2_run2.sorted.bam")
S = "SGNex_H9_directRNA_replicate2_run2"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
ENST = re.compile(r"(ENST\d+)")
TOL = 25           # +/- bp window counted as "at the TSS"
BIN = 50           # bp bin for background histogram
MIN_READS = 20     # skip low-depth loci (peak not detectable)
MAX_LOCI = 400     # cap for speed


def matched_eq(track):
    s = set()
    for ln in open(track):
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 4 and c[3] == "=":
            m = ENST.search(c[2])
            if m:
                s.add(m.group(1))
    return s


# --- truth transcripts: exons -> (chr, strand, TSS, TES, intron chain) ---
tx = {}
for ln in open(GTF):
    if ln.startswith("#"):
        continue
    f = ln.split("\t")
    if len(f) < 9 or f[2] != "exon":
        continue
    m = ENST.search(f[8])
    if not m:
        continue
    tx.setdefault(m.group(1), [f[0], f[6], []])[2].append((int(f[3]) - 1, int(f[4])))

TX = {}
for e, (chrom, strand, exons) in tx.items():
    exons.sort()
    introns = tuple((exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1))
    tss = exons[0][0] if strand == "+" else exons[-1][1]
    tes = exons[-1][1] if strand == "+" else exons[0][0]
    span = (exons[0][0], exons[-1][1])
    TX[e] = dict(chrom=chrom, strand=strand, tss=tss, tes=tes,
                 introns=introns, span=span)

off = matched_eq(f"{HERE}/gc_full/{S}__p00__pyfin.tracking")
grf = matched_eq(f"{HERE}/gc_graph/{S}__p00__pyfin.tracking")
lost = off - grf
present = grf
print(f"OFF matched= {len(off)}  GRAPH matched= {len(grf)}  LOST {len(lost)}")

# index ALL truth transcripts by 3'-terminal intron for fast long-partner lookup
# (a real contained short isoform is defined by truth structure, independent of
#  whether the long partner happened to survive in this run)
by_term = {}
for e in TX:
    t = TX.get(e)
    if not t or not t["introns"]:
        continue
    term = t["introns"][-1] if t["strand"] == "+" else t["introns"][0]
    by_term.setdefault((t["chrom"], t["strand"], term), []).append(e)


def is_5p_internal(shortE):
    """Return a surviving longer transcript that `shortE` is a 5'-truncated
    contiguous 3'-suffix of (its TSS internal to the long one), else None."""
    s = TX[shortE]
    if not s["introns"]:
        return None
    term = s["introns"][-1] if s["strand"] == "+" else s["introns"][0]
    for e in by_term.get((s["chrom"], s["strand"], term), []):
        L = TX[e]
        if len(L["introns"]) <= len(s["introns"]):
            continue
        if s["strand"] == "+":
            if L["introns"][-len(s["introns"]):] == s["introns"] and L["tss"] < s["tss"]:
                return e
        else:
            if L["introns"][:len(s["introns"])] == s["introns"] and L["tss"] > s["tss"]:
                return e
    return None


internal = [(e, is_5p_internal(e)) for e in lost]
internal = [(e, L) for e, L in internal if L]
print(f"LOST that are 5'-internal to a surviving longer transcript: {len(internal)}")

bam = pysam.AlignmentFile(BAM, "rb")
enrich = []
fracs = []
examples = []
for shortE, longE in internal[:MAX_LOCI]:
    s = TX[shortE]
    lo, hi = TX[longE]["span"]
    if hi - lo > 200000:
        continue
    fivep = []
    for r in bam.fetch(s["chrom"], max(0, lo), hi):
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue
        p = r.reference_start if s["strand"] == "+" else r.reference_end
        if p is None:
            continue
        if lo <= p <= hi:
            fivep.append(p)
    if len(fivep) < MIN_READS:
        continue
    # count at the short-isoform TSS window
    at = sum(1 for p in fivep if abs(p - s["tss"]) <= TOL)
    # smooth background: histogram in BIN bins across the body, median density
    nb = max(1, (hi - lo) // BIN)
    hist = [0] * nb
    for p in fivep:
        b = min(nb - 1, (p - lo) // BIN)
        hist[b] += 1
    # scale median per-bin count to the TOL window width
    med_bin = statistics.median(hist) if hist else 0
    bg = med_bin * (2 * TOL + 1) / BIN
    e = at / bg if bg > 0 else (float("inf") if at > 0 else 0)
    enrich.append(e)
    fracs.append(at / len(fivep))
    examples.append((shortE, longE, len(fivep), at, round(bg, 1), round(e, 1), s, lo, hi, fivep))

bam.close()
enrich_f = [e for e in enrich if e != float("inf")]
print(f"\nLoci with enough depth tested: {len(enrich)}")
if enrich:
    fin = enrich_f or [0]
    print(f"TSS 5'-end enrichment vs smooth background (>1 = real TSS peak):")
    print(f"  median={statistics.median(fin):.2f}  mean={statistics.mean(fin):.2f}  "
          f"max={max(fin):.1f}")
    for thr in (1.5, 2, 3, 5):
        frac = 100 * sum(1 for e in enrich if e >= thr) / len(enrich)
        print(f"  fraction of loci with enrichment >= {thr}x : {frac:.0f}%")
    print(f"\nRaw pile-up: fraction of a locus's reads whose 5'-end is within "
          f"+/-{TOL}bp of the lost-isoform TSS:")
    fracs.sort(reverse=True)
    print("  " + "  ".join(f"{f*100:.0f}%" for f in fracs))
    print(f"  median={statistics.median(fracs)*100:.0f}%  "
          f"loci with >=40% of reads piled at TSS: "
          f"{100*sum(1 for f in fracs if f>=0.4)/len(fracs):.0f}%")

# a few concrete ASCII histograms (highest enrichment examples)
examples.sort(key=lambda x: -(x[5] if x[5] != float("inf") else 1e9))
print("\n=== example loci (read 5'-end histogram; ^ = lost short-isoform TSS) ===")
for shortE, longE, n, at, bg, e, s, lo, hi, fivep in examples[:4]:
    nb = min(60, max(10, (hi - lo) // BIN))
    step = (hi - lo) / nb
    hist = [0] * nb
    for p in fivep:
        hist[min(nb - 1, int((p - lo) / step))] += 1
    tss_bin = min(nb - 1, int((s["tss"] - lo) / step))
    mx = max(hist) or 1
    print(f"\n{shortE} (short, strand {s['strand']}) inside {longE} | "
          f"reads={n} at_TSS={at} bg~{bg} enrich={e}x")
    for i, h in enumerate(hist):
        bar = "#" * int(30 * h / mx)
        mark = " <== short TSS" if i == tss_bin else ""
        print(f"  {'^' if i==tss_bin else ' '}{bar:<30}{h:>4}{mark}")
