"""Read-only final-assembly junction nominations, never fusion read votes."""

from __future__ import annotations

import csv
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from functools import cached_property

from fin.fusion.correction import junction_locus_key, path_junctions


@dataclass(frozen=True)
class AssemblyJunctions:
    junctions: tuple = ()
    min_clean_reads: int = 3

    @cached_property
    def _indexes(self):
        introns, donors, acceptors = {}, {}, {}
        for chrom, strand, donor, acceptor, clean in self.junctions:
            locus = chrom, strand
            introns.setdefault(locus, {})[donor, acceptor] = clean
            donor_counts = donors.setdefault(locus, {})
            acceptor_counts = acceptors.setdefault(locus, {})
            donor_counts[donor] = max(donor_counts.get(donor, 0), clean)
            acceptor_counts[acceptor] = max(acceptor_counts.get(acceptor, 0), clean)
        return introns, donors, acceptors

    @cached_property
    def _sites(self):
        return tuple({locus: tuple(sorted(counts)) for locus, counts in index.items()}
                     for index in self._indexes[1:])

    def _nearby(self, arm, position, side, tolerance):
        sites = self._sites[side].get((arm.chrom, arm.strand), ())
        return sites[bisect_left(sites, position - tolerance):bisect_right(sites, position + tolerance)]

    def registry(self, paths, *, tolerance=3):
        registry = {}
        for path in paths:
            for kind, _, _, left, right in path_junctions(path):
                observed = left.exit_position, right.entry_position
                positions = set()
                donors = self._nearby(left, observed[0], 0, tolerance)
                acceptors = self._nearby(right, observed[1], 1, tolerance)
                if kind == "intron":
                    introns = self._indexes[0].get((left.chrom, left.strand), {})
                    positions.update((donor, acceptor) for donor in donors for acceptor in acceptors
                                     if (donor, acceptor) in introns)
                else:
                    positions.update((donor, acceptor) for donor in set(donors) | {observed[0]}
                                     for acceptor in set(acceptors) | {observed[1]})
                nearby = {position for position in positions
                          if max(abs(value - baseline) for value, baseline in zip(position, observed)) <= tolerance}
                registry.setdefault(junction_locus_key(kind, left, right), set()).update(nearby)
        return {key: tuple(sorted(values)) for key, values in registry.items()}

    def annotate(self, path):
        result = []
        for kind, arm_index, block_index, left, right in path_junctions(path):
            introns, donors, acceptors = self._indexes
            left_locus, right_locus = (left.chrom, left.strand), (right.chrom, right.strand)
            donor_support = donors.get(left_locus, {}).get(left.exit_position, 0)
            acceptor_support = acceptors.get(right_locus, {}).get(right.entry_position, 0)
            exact = (introns.get(left_locus, {}).get((left.exit_position, right.entry_position), 0)
                     if kind == "intron" else 0)
            result.append(dict(kind=kind, arm_index=arm_index, block_index=block_index,
                               exact_intron_clean_reads=exact, donor_clean_reads=donor_support,
                               acceptor_clean_reads=acceptor_support,
                               interpretation="auxiliary_structure_not_signal_or_fusion_support"))
        return tuple(result)


def load_assembly_junctions(gtf_path, support_path, *, min_clean_reads=3):
    """Require final-model membership plus clean reads; overlapping rows use max.

    Aggregated counts cannot exclude the current molecule, so they nominate and
    annotate only. They never independently validate a fusion or a correction.
    """
    from fin.io.io_gtf import GTFReader

    if min_clean_reads < 2:
        raise ValueError("assembly junction evidence needs multiple clean reads")
    final = set()
    with GTFReader(str(gtf_path)) as reader:
        reader.parse()
        for transcript in reader.transcripts.values():
            if transcript.strand not in ("+", "-"):
                continue
            exons = sorted(set(transcript.exons))
            final.update((transcript.chrom, transcript.strand, previous[1], following[0])
                         for previous, following in zip(exons, exons[1:]) if previous[1] < following[0])
    counts = {}
    with open(support_path) as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        if not {"chrom", "strand", "donor", "acceptor", "raw_reads", "indel_clean_reads"} <= set(rows.fieldnames or ()):
            raise ValueError("invalid assembly junction support header")
        for row in rows:
            key = row["chrom"], row["strand"], int(row["donor"]), int(row["acceptor"])
            clean, raw = int(row["indel_clean_reads"]), int(row["raw_reads"])
            if not 0 <= clean <= raw or not 0 <= key[2] < key[3]:
                raise ValueError("invalid assembly junction support counts or coordinates")
            if key in final and clean >= min_clean_reads:
                counts[key] = max(counts.get(key, 0), clean)
    return AssemblyJunctions(tuple((chrom, strand, start if strand == "+" else end,
                                    end if strand == "+" else start, clean)
                                   for (chrom, strand, start, end), clean in sorted(counts.items())), min_clean_reads)
