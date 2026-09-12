"""Strand-aware gene labels and ordered fusion arms without assembly mutation."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, replace
from typing import Iterable, Tuple

from fin.fusion.paths import ReadPath
from fin.fusion.records import AlignmentBlock


def _merge_intervals(intervals):
    merged = []
    for start, end in sorted(set(intervals)):
        if start < 0 or end <= start:
            raise ValueError("invalid annotation exon")
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return tuple(merged)


@dataclass(frozen=True)
class GeneModel:
    gene_id: str
    chrom: str
    strand: str
    exons: Tuple[Tuple[int, int], ...]

    def __post_init__(self):
        if not self.gene_id or not self.chrom or self.strand not in ("+", "-"):
            raise ValueError("gene model needs identity, chromosome and strand")
        if not self.exons or self.exons != _merge_intervals(self.exons):
            raise ValueError("gene model exons must be merged and sorted")


class GeneIndex:
    """Exonic overlap index; annotation nominates identity, never read support."""

    def __init__(self, genes: Iterable[GeneModel] = (), *, bin_size: int = 16384):
        if bin_size <= 0:
            raise ValueError("gene-index bin size must be positive")
        self.bin_size = bin_size
        self.genes = tuple(sorted(set(genes), key=lambda gene: (
            gene.chrom, gene.strand, gene.gene_id, gene.exons,
        )))
        self._bins = defaultdict(set)
        for index, gene in enumerate(self.genes):
            for start, end in gene.exons:
                for bucket in range(start // bin_size, (end - 1) // bin_size + 1):
                    self._bins[(gene.chrom, gene.strand, bucket)].add(index)

    @classmethod
    def from_transcripts(cls, transcripts):
        grouped = defaultdict(list)
        for transcript in transcripts:
            if transcript.gene_id and transcript.strand in ("+", "-"):
                grouped[(transcript.gene_id, transcript.chrom, transcript.strand)].extend(transcript.exons)
        return cls(GeneModel(gene_id, chrom, strand, _merge_intervals(exons))
                   for (gene_id, chrom, strand), exons in grouped.items() if exons)

    def labels(self, chrom, strand, blocks, *, min_overlap_bases=25, min_fraction=0.5):
        if min_overlap_bases <= 0 or not 0 < min_fraction <= 1:
            raise ValueError("invalid gene-overlap thresholds")
        intervals = _merge_intervals((block.reference_start, block.reference_end) for block in blocks)
        if not intervals:
            return ()
        candidates = set()
        for start, end in intervals:
            for bucket in range(start // self.bin_size, (end - 1) // self.bin_size + 1):
                candidates.update(self._bins.get((chrom, strand, bucket), ()))
        span = sum(end - start for start, end in intervals)
        supported = set()
        for index in candidates:
            gene = self.genes[index]
            overlap = sum(max(0, min(end, gene_end) - max(start, gene_start))
                          for start, end in intervals for gene_start, gene_end in gene.exons)
            if overlap >= min_overlap_bases and overlap / span >= min_fraction:
                supported.add(gene.gene_id)
        return tuple(sorted(supported))


@dataclass(frozen=True)
class FusionArm:
    chrom: str
    strand: str
    blocks: Tuple[AlignmentBlock, ...]
    gene_ids: Tuple[str, ...]

    @property
    def entry_position(self):
        return self.blocks[0].reference_start if self.strand == "+" else self.blocks[0].reference_end

    @property
    def exit_position(self):
        return self.blocks[-1].reference_end if self.strand == "+" else self.blocks[-1].reference_start

    @property
    def query_start(self):
        return self.blocks[0].query_start

    @property
    def query_end(self):
        return self.blocks[-1].query_end

    @property
    def exons(self):
        return tuple((block.reference_start, block.reference_end) for block in self.blocks)


@dataclass(frozen=True)
class AnnotatedReadPath:
    read_id: str
    arms: Tuple[FusionArm, ...]
    event_class: str
    identity_status: str

    @property
    def junctions(self):
        return tuple(((previous.chrom, previous.exit_position, previous.strand),
                      (following.chrom, following.entry_position, following.strand))
                     for previous, following in zip(self.arms, self.arms[1:]))

    @property
    def gene_path(self):
        return tuple(arm.gene_ids for arm in self.arms)


def _collinear(previous, following):
    if previous.chrom != following.chrom or previous.strand != following.strand:
        return False
    if previous.query_end > following.query_start:
        return False
    if previous.strand == "+":
        return previous.exit_position <= following.entry_position
    return previous.exit_position >= following.entry_position


def annotate_read_path(path: ReadPath, genes: GeneIndex, *, readthrough_distance=100000):
    """Split a spliced alignment at annotated gene transitions; retain ambiguity."""
    if readthrough_distance < 0:
        raise ValueError("readthrough distance must be nonnegative")
    arms = []
    for segment in path.segments:
        block_labels = [genes.labels(segment.chrom, segment.strand, (block,), min_overlap_bases=1)
                        for block in segment.blocks]
        segment_groups = []
        for block, labels in zip(segment.blocks, block_labels):
            if segment_groups and labels == segment_groups[-1][1]:
                segment_groups[-1][0].append(block)
            else:
                segment_groups.append(([block], labels))
        for blocks, labels in segment_groups:
            arm = FusionArm(segment.chrom, segment.strand, tuple(blocks), labels)
            if arms and len(labels) == 1 and arms[-1].gene_ids == labels and _collinear(arms[-1], arm):
                arms[-1] = replace(arms[-1], blocks=arms[-1].blocks + arm.blocks)
            else:
                arms.append(arm)
    identities = [arm.gene_ids for arm in arms]
    identity_status = ("ambiguous" if any(len(labels) > 1 for labels in identities)
                       else "unannotated" if any(not labels for labels in identities)
                       else "resolved")
    if len(arms) == 1:
        event_class = "ordinary_alignment"
    elif identity_status == "ambiguous":
        event_class = "ambiguous_gene_path"
    elif any(previous.chrom != following.chrom for previous, following in zip(arms, arms[1:])):
        event_class = "interchromosomal_candidate"
    elif all(_collinear(previous, following)
             and abs(previous.exit_position - following.entry_position) <= readthrough_distance
             for previous, following in zip(arms, arms[1:])):
        event_class = "readthrough_or_splicing_candidate"
    elif identity_status == "resolved" and len({labels[0] for labels in identities}) == 1:
        event_class = "intragenic_rearrangement_candidate"
    else:
        event_class = "intrachromosomal_candidate"
    return AnnotatedReadPath(path.read_id, tuple(arms), event_class, identity_status)
