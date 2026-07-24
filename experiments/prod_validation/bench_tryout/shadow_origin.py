#!/usr/bin/env python3
"""Are the wobble shadows from ALIGNMENT ambiguity (fixable by read-consensus snap)
or from genuinely different reads? For expressed truth transcripts that pyfin matched
'=', take each junction and look at how reads (minimap CIGAR introns) place it. If a
strong majority sit at ONE coordinate (the mode) and a minority scatter +/- a few bp,
the shadows are alignment wobble and snapping every read to the per-junction mode
recovers the truth. Report the read-coordinate concentration at the mode."""
import re, statistics
import pysam
from collections import defaultdict, Counter

B = "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
TRUTH = f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
S = "SGNex_H9_directRNA_replicate2_run2"
HERE = "/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
BAM = f"{B}/results/gencode_full_sweep/{S}/align/{S}.sorted.bam"
NC = f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
PYF_TRK = f"{HERE}/gc_frontier/{S}__p00__pyfin.tracking"
ENST = re.compile(r"(ENST\d+)")
WIN = 15   # bp window to gather reads' version of a junction

est = {}
for ln in open(NC):
    p = ln.rstrip("\n").split("\t")
    if len(p) >= 3:
        m = ENST.search(p[0])
        if m:
            try: est[m.group(1)] = float(p[2])
            except ValueError: pass
expr3 = {k for k, v in est.items() if v >= 3}
# transcripts pyfin matched '=' (has candidates) & expressed -> real loci with shadows
matched = set()
for ln in open(PYF_TRK):
    c = ln.rstrip("\n").split("\t")
    if len(c) >= 4 and c[3] == "=":
        m = ENST.search(c[2])
        if m: matched.add(m.group(1))
target = matched & expr3

# truth junctions for target transcripts
tx = {}
for ln in open(TRUTH):
    if ln.startswith("#"): continue
    f = ln.split("\t")
    if len(f) < 9 or f[2] != "exon": continue
    m = ENST.search(f[8])
    if not m or m.group(1) not in target: continue
    tx.setdefault(m.group(1), [f[0], []])[1].append((int(f[3]) - 1, int(f[4])))
junctions = []   # (chrom, donor, acceptor)
for e, (chrom, ex) in tx.items():
    ex.sort()
    for i in range(len(ex) - 1):
        junctions.append((chrom, ex[i][1], ex[i + 1][0]))
# de-dup, cap for speed
junctions = list({j for j in junctions})[:4000]

def read_introns(r):
    out = []; pos = r.reference_start
    for op, ln in r.cigartuples or []:
        if op in (0, 2, 7, 8): pos += ln
        elif op == 3: out.append((pos, pos + ln)); pos += ln
    return out

bam = pysam.AlignmentFile(BAM, "rb")
mode_fracs = []; spreads = []; offsets = Counter()
analyzed = 0
for chrom, td, ta in junctions:
    # gather reads whose CIGAR has an intron with donor within WIN of td (acceptor free)
    coords = []
    for r in bam.fetch(chrom, max(0, td - 200), ta + 200):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        for (d, a) in read_introns(r):
            if abs(d - td) <= WIN and abs(a - ta) <= WIN:
                coords.append((d, a))
    if len(coords) < 5:
        continue
    analyzed += 1
    cc = Counter(coords)
    mode, mode_n = cc.most_common(1)[0]
    mode_fracs.append(mode_n / len(coords))
    # spread = stdev of donor offset from true
    doff = [d - td for d, a in coords]
    spreads.append(statistics.pstdev(doff) if len(doff) > 1 else 0)
    offsets[mode[0] - td] += 1
bam.close()

mode_fracs.sort()
print(f"analyzed {analyzed} expressed-TP junctions (>=5 spanning reads)")
print(f"\nread-coordinate concentration at the dominant mode:")
print(f"  median fraction of reads AT the mode coordinate: {statistics.median(mode_fracs)*100:.0f}%")
print(f"  fraction of junctions where mode carries >=80% of reads: "
      f"{100*sum(1 for f in mode_fracs if f>=0.8)/len(mode_fracs):.0f}%")
print(f"  fraction where mode carries >=90%: "
      f"{100*sum(1 for f in mode_fracs if f>=0.9)/len(mode_fracs):.0f}%")
print(f"  median donor-coordinate spread (stdev, bp): {statistics.median(spreads):.1f}")
print(f"\ndoes the mode == the TRUE (annotated) coordinate? offset of mode from truth:")
tot = sum(offsets.values())
for off in sorted(offsets, key=lambda x: -offsets[x])[:7]:
    print(f"  mode offset {off:+d}bp from truth: {offsets[off]} ({100*offsets[off]/tot:.0f}%)")
print("\nINTERPRETATION: high mode-fraction + mode-at-truth => shadows are ALIGNMENT wobble")
print("(minority of reads mis-placed); snapping every read to the per-junction mode recovers truth.")
