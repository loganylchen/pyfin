#!/usr/bin/env python3
"""Is part of our spike-in junction wobble an ALIGNER artifact?

minimap2 scores a GT donor higher when the base immediately after it is A or G
(true for ~91% of human introns, which is where the prior comes from). Both our
spike-in and human grids were aligned with that default (`-ax splice -uf -k14`,
no `--splice-flank=no`). Measured on our own references:

    human   93.9% of GT donors are followed by A/G   (prior HOLDS)
    SIRV4   62.3%                                    (prior VIOLATED)
    HEYA8   50.9%                                    (prior is noise)

So the prior is right for the real target and wrong for the two grids most of
this project's junction machinery was tuned on. This measures the CONSEQUENCE
without realigning: for each annotated junction, what fraction of the reads
that cross it place it EXACTLY, split by whether its post-GT base satisfies the
prior. If the prior is manufacturing wobble, the non-A/G junctions should place
worse -- and the gap should be large on the spike-ins and absent on human.

A within-alignment comparison, so aligner version, preset and read chemistry
are held fixed by construction; only the sequence context differs.
"""
import argparse, collections, re, sys


def truth_junctions(gtf, fasta, chroms=None):
    """{(chrom, strand, donor, acceptor): flank_is_AG}"""
    ex, meta = collections.defaultdict(list), {}
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        if chroms and f[0] not in chroms:
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        ex[m.group(1)].append((int(f[3]) - 1, int(f[4])))
        meta[m.group(1)] = (f[0], f[6])
    rc = lambda x: x[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))
    out = {}
    for t, e in ex.items():
        e.sort()
        c, s = meta[t]
        for i in range(len(e) - 1):
            d, acc = e[i][1], e[i + 1][0]
            try:
                if s == "+":
                    donor = fasta.fetch(c, d, d + 2).upper()
                    nxt = fasta.fetch(c, d + 2, d + 3).upper()
                else:
                    donor = rc(fasta.fetch(c, acc - 2, acc).upper())
                    nxt = rc(fasta.fetch(c, acc - 3, acc - 2).upper())
            except Exception:
                continue
            if donor != "GT":
                continue
            out[(c, s, d, acc)] = nxt in ("A", "G")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--chroms", default="")
    ap.add_argument("--tol", type=int, default=10,
                    help="a read junction within this many bp of a truth "
                         "junction on BOTH boundaries is counted as an attempt "
                         "at it (wobble), beyond this it is a different event")
    ap.add_argument("--label", default="grid")
    a = ap.parse_args()
    import pysam

    fasta = pysam.FastaFile(a.genome)
    chroms = set(a.chroms.split(",")) if a.chroms else None
    truth = truth_junctions(a.gtf, fasta, chroms)
    fasta.close()
    # bucket truth by (chrom, donor // BIN) so the per-read-junction lookup is
    # O(1) instead of a scan over every annotated junction -- on human that scan
    # is 25k junctions per read junction and does not finish.
    BIN = 256
    index = collections.defaultdict(list)
    for (c, s, d, acc), ag in truth.items():
        for b in {(d - a.tol) // BIN, (d + a.tol) // BIN}:
            index[(c, b)].append((d, acc, ag))

    exact = collections.Counter()
    near = collections.Counter()
    bam = pysam.AlignmentFile(a.bam, "rb")
    for aln in bam.fetch(until_eof=False) if bam.has_index() else bam:
        if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
            continue
        if aln.mapping_quality < 1:
            continue
        if chroms and aln.reference_name not in chroms:
            continue
        pos = aln.reference_start
        for op, ln in aln.cigartuples or ():
            if op in (0, 7, 8, 2):
                pos += ln
            elif op == 3:
                d, acc = pos, pos + ln
                pos += ln
                # nearest truth junction within tol on BOTH boundaries
                best = None
                for td, tacc, ag in index.get((aln.reference_name, d // BIN), ()):
                    if abs(td - d) <= a.tol and abs(tacc - acc) <= a.tol:
                        off = abs(td - d) + abs(tacc - acc)
                        if best is None or off < best[0]:
                            best = (off, ag)
                if best is None:
                    continue
                off, ag = best
                near[ag] += 1
                if off == 0:
                    exact[ag] += 1
    bam.close()

    print(f"{a.label}: read junctions landing within {a.tol}bp of an annotated "
          f"GT junction, by whether that junction satisfies minimap2's "
          f"post-GT A/G prior\n")
    print(f"{'post-GT base':<16}{'read junctions':>16}{'placed EXACTLY':>16}{'exact rate':>12}")
    for ag, name in ((True, "A/G (prior ok)"), (False, "C/T (violates)")):
        n, e = near[ag], exact[ag]
        if n:
            print(f"{name:<16}{n:>16d}{e:>16d}{100*e/n:>11.1f}%")
    if near[True] and near[False]:
        gap = 100 * exact[True] / near[True] - 100 * exact[False] / near[False]
        print(f"\ngap = {gap:+.1f} pp. A large positive gap means the aligner "
              f"places junctions worse exactly where its own prior is wrong, "
              f"i.e. some of the wobble is manufactured.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
