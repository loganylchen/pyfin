#!/usr/bin/env python3
"""Freeze the K562 BCR-ABL1 positive control as a compact test fixture.

Takes the 98 real records overlapping the BCR and ABL1 loci (SG-NEx K562
direct RNA rep4 run1, GRCh38, minimap2 2.30 `-ax splice -uf -k14`), their
unmodified full read sequences, every BCR and ABL1 transcript of GENCODE v44,
and two 300-kb reference slices. Coordinates are shifted onto the slices;
nothing else is changed -- SA tags (including the summarised `175M723D`
CIGAR that motivated locator reconciliation), hard clips, indels, tails and
MAPQ are the originals. Whole transcripts are kept so exon ranks (e14a2) are
computed on the real exon order.
"""
import re, subprocess, sys
from pathlib import Path

import pysam

SRC = Path("/SSD/logan/dev/pyfin/experiments/fusion/_k562_bcrabl/fixture")
DST = Path("/SSD/logan/dev/pyfin/tests/data/fusion_k562_bcr_abl1")
SLICES = {"chr22": (23100000, 23400000), "chr9": (130600000, 130900000)}


def shift_cigar_position(chrom, pos0):
    lo, hi = SLICES[chrom]
    assert lo <= pos0 < hi, (chrom, pos0)
    return pos0 - lo


def main():
    DST.mkdir(parents=True, exist_ok=True)
    fasta = pysam.FastaFile(str(SRC / "chr9_chr22.fa"))
    with (DST / "reference.fa").open("w") as out:
        for chrom, (lo, hi) in SLICES.items():
            seq = fasta.fetch(chrom, lo, hi)
            out.write(f">{chrom}\n")
            for i in range(0, len(seq), 80):
                out.write(seq[i:i + 80] + "\n")
    subprocess.run(["samtools", "faidx", str(DST / "reference.fa")], check=True)
    with (DST / "annotation.gtf").open("w") as out:
        for line in (SRC / "chr9_chr22.gtf").open():
            if 'gene_name "BCR"' not in line and 'gene_name "ABL1"' not in line:
                continue
            f = line.rstrip("\n").split("\t")
            lo, hi = SLICES[f[0]]
            start, end = int(f[3]), int(f[4])
            assert lo < start and end <= hi, line[:80]
            f[3], f[4] = str(start - lo), str(end - lo)
            out.write("\t".join(f) + "\n")
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": c, "LN": hi - lo} for c, (lo, hi) in SLICES.items()]}
    names = list(SLICES)
    dropped = []
    with pysam.AlignmentFile(str(SRC / "bcr_abl1.bam"), "rb") as reader, \
            pysam.AlignmentFile(str(DST / "reads.unsorted.bam"), "wb", header=header) as writer:
        for record in reader:
            # a read whose SA partner lies outside the two slices is a chimera
            # with somewhere else, not part of this positive control; it is
            # left out entirely rather than stripped of its partner (which
            # would silently turn a chimera into a clipped normal read)
            sa = record.get_tag("SA") if record.has_tag("SA") else ""
            if any(entry.split(",")[0] not in SLICES for entry in sa.rstrip(";").split(";") if entry):
                dropped.append(record.query_name)
                continue
            new = pysam.AlignedSegment(writer.header)
            new.query_name = record.query_name
            new.query_sequence = record.query_sequence
            new.query_qualities = record.query_qualities
            new.flag = record.flag
            new.reference_id = names.index(record.reference_name)
            new.reference_start = shift_cigar_position(record.reference_name, record.reference_start)
            new.mapping_quality = record.mapping_quality
            new.cigarstring = record.cigarstring
            tags = []
            for tag, value in record.get_tags():
                if tag == "SA":
                    entries = []
                    for entry in value.rstrip(";").split(";"):
                        chrom, pos, strand, cigar, mapq, nm = entry.split(",")
                        entries.append(",".join((chrom, str(shift_cigar_position(chrom, int(pos) - 1) + 1),
                                                 strand, cigar, mapq, nm)))
                    value = ";".join(entries) + ";"
                tags.append((tag, value))
            new.set_tags(tags)
            writer.write(new)
    subprocess.run(["samtools", "sort", "-o", str(DST / "reads.bam"), str(DST / "reads.unsorted.bam")], check=True)
    (DST / "reads.unsorted.bam").unlink()
    subprocess.run(["samtools", "index", str(DST / "reads.bam")], check=True)
    kept = set()
    with pysam.AlignmentFile(str(DST / "reads.bam"), "rb") as reader:
        kept = {record.query_name for record in reader}
    with (SRC / "reads.fastq").open() as fin, (DST / "reads.fastq").open("w") as fout:
        while True:
            block = [fin.readline() for _ in range(4)]
            if not block[0]:
                break
            if block[0][1:].split()[0] in kept:
                fout.writelines(block)
    (DST / "README.md").write_text(
        "K562 (SG-NEx directRNA replicate4 run1) records overlapping BCR (chr22) and ABL1 (chr9),\n"
        "GRCh38 coordinates shifted onto two 300-kb slices (chr22:23,100,000-23,400,000 and\n"
        "chr9:130,600,000-130,900,000). Every BCR/ABL1 GENCODE v44 transcript is kept whole so exon\n"
        "ranks are real. Reads whose SA partner lies outside the slices were left out entirely:\n"
        + "".join(f"  {name}\n" for name in sorted(set(dropped)))
        + "Built by experiments/fusion/build_bcr_abl1_fixture.py.\n")
    print("done", sorted(p.name for p in DST.iterdir()), "dropped", sorted(set(dropped)))


if __name__ == "__main__":
    sys.exit(main())
