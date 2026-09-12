"""Alignment evidence in full, forward-query coordinates for fusion analysis."""

from __future__ import annotations

import re
from dataclasses import dataclass, replace
from typing import Iterable, Optional, Tuple


_CIGAR_CODES = {symbol: index for index, symbol in enumerate("MIDNSHP=X")}
_QUERY_CODES = frozenset((0, 1, 4, 5, 7, 8))
_ALIGNED_CODES = frozenset((0, 7, 8))


@dataclass(frozen=True)
class AlignmentBlock:
    reference_start: int
    reference_end: int
    query_start: int
    query_end: int
    aligned_bases: int

    def __post_init__(self):
        if not 0 <= self.reference_start < self.reference_end:
            raise ValueError("alignment block has invalid reference coordinates")
        if not 0 <= self.query_start < self.query_end:
            raise ValueError("alignment block has invalid query coordinates")
        if self.aligned_bases <= 0:
            raise ValueError("alignment block must contain aligned query bases")


@dataclass(frozen=True)
class AlignmentSegment:
    read_id: str
    chrom: str
    strand: str
    read_length: int
    blocks: Tuple[AlignmentBlock, ...]
    mapping_quality: int
    edit_distance: Optional[int] = None
    sources: Tuple[str, ...] = ()
    secondary: bool = False

    def __post_init__(self):
        if not self.read_id or not self.chrom or self.strand not in ("+", "-"):
            raise ValueError("alignment segment requires read, chromosome and strand")
        if not self.blocks or self.read_length <= 0:
            raise ValueError("alignment segment requires blocks and full read length")
        if not 0 <= self.mapping_quality <= 255:
            raise ValueError("mapping quality must be in [0, 255]")
        if self.edit_distance is not None and self.edit_distance < 0:
            raise ValueError("edit distance must be nonnegative")
        if any(block.query_end > self.read_length for block in self.blocks):
            raise ValueError("alignment block extends past full read length")
        for previous, following in zip(self.blocks, self.blocks[1:]):
            if following.query_start < previous.query_end:
                raise ValueError("blocks must follow forward-query order")
            reference_ordered = (
                previous.reference_end <= following.reference_start
                if self.strand == "+"
                else following.reference_end <= previous.reference_start
            )
            if not reference_ordered:
                raise ValueError("blocks must follow the segment's genomic orientation")

    @property
    def query_start(self) -> int:
        return self.blocks[0].query_start

    @property
    def query_end(self) -> int:
        return self.blocks[-1].query_end

    @property
    def aligned_bases(self) -> int:
        return sum(block.aligned_bases for block in self.blocks)

    @property
    def effective_mapping_quality(self) -> int:
        return 0 if self.mapping_quality == 255 else self.mapping_quality

    @property
    def exons(self) -> Tuple[Tuple[int, int], ...]:
        return tuple((block.reference_start, block.reference_end) for block in self.blocks)

    @property
    def entry_position(self) -> int:
        block = self.blocks[0]
        return block.reference_start if self.strand == "+" else block.reference_end

    @property
    def exit_position(self) -> int:
        block = self.blocks[-1]
        return block.reference_end if self.strand == "+" else block.reference_start

    @property
    def geometry_key(self) -> tuple:
        return (
            self.read_id, self.chrom, self.strand, self.read_length,
            tuple((block.reference_start, block.reference_end,
                   block.query_start, block.query_end) for block in self.blocks),
        )


def parse_cigar(cigar: str) -> Tuple[Tuple[int, int], ...]:
    """Parse SAM CIGAR without accepting partial strings or internal clips."""
    tokens = re.findall(r"([1-9][0-9]*)([MIDNSHP=X])", cigar)
    if not tokens or "".join(length + symbol for length, symbol in tokens) != cigar:
        raise ValueError(f"invalid CIGAR: {cigar!r}")
    operations = "".join(symbol for _, symbol in tokens)
    if not re.fullmatch(r"H?S?[MIDNP=X]+S?H?", operations):
        raise ValueError(f"invalid CIGAR clipping: {cigar!r}")
    if not any(symbol in "M=X" for _, symbol in tokens):
        raise ValueError("CIGAR has no aligned query bases")
    return tuple((_CIGAR_CODES[symbol], int(length)) for length, symbol in tokens)


def segment_from_cigar(
    read_id: str,
    chrom: str,
    reference_start: int,
    strand: str,
    cigar: str,
    mapping_quality: int,
    *,
    edit_distance: Optional[int] = None,
    source: str = "bam",
    secondary: bool = False,
) -> AlignmentSegment:
    """Include hard clips in the molecule axis, then reverse that axis if needed."""
    if reference_start < 0 or strand not in ("+", "-"):
        raise ValueError("invalid reference start or strand")
    operations = parse_cigar(cigar)
    read_length = sum(length for code, length in operations if code in _QUERY_CODES)
    core = [index for index, (code, _) in enumerate(operations) if code not in (4, 5)]
    query_cursor = sum(length for _, length in operations[:core[0]])
    reference_cursor = reference_start
    block_query_start = query_cursor
    block_reference_start = reference_cursor
    aligned_bases = 0
    blocks = []
    for code, length in operations[core[0]:core[-1] + 1]:
        if code in _ALIGNED_CODES:
            query_cursor += length
            reference_cursor += length
            aligned_bases += length
        elif code == 1:
            query_cursor += length
        elif code == 2:
            reference_cursor += length
        elif code == 3:
            blocks.append(AlignmentBlock(
                block_reference_start, reference_cursor,
                block_query_start, query_cursor, aligned_bases,
            ))
            reference_cursor += length
            block_reference_start = reference_cursor
            block_query_start = query_cursor
            aligned_bases = 0
    blocks.append(AlignmentBlock(
        block_reference_start, reference_cursor,
        block_query_start, query_cursor, aligned_bases,
    ))
    if strand == "-":
        blocks = [AlignmentBlock(
            block.reference_start, block.reference_end,
            read_length - block.query_end, read_length - block.query_start,
            block.aligned_bases,
        ) for block in reversed(blocks)]
    return AlignmentSegment(
        read_id, chrom, strand, read_length, tuple(blocks), mapping_quality,
        edit_distance, (source,), secondary,
    )


def parse_sa_tag(read_id: str, value: str) -> Tuple[AlignmentSegment, ...]:
    """Parse every SA entry atomically; a malformed alternative is not ignored."""
    if not value:
        return ()
    entries = value[:-1].split(";") if value.endswith(";") else value.split(";")
    segments = []
    for entry in entries:
        fields = entry.split(",")
        if len(fields) != 6:
            raise ValueError(f"invalid SA entry for {read_id!r}: expected six fields")
        chrom, position, strand, cigar, quality, distance = fields
        if not chrom or chrom == "*" or any(character.isspace() for character in chrom):
            raise ValueError("invalid SA reference name")
        if not all(re.fullmatch(r"[0-9]+", number) for number in (position, quality, distance)):
            raise ValueError("invalid SA numeric field")
        segments.append(segment_from_cigar(
            read_id, chrom, int(position) - 1, strand, cigar, int(quality),
            edit_distance=int(distance), source="sa",
        ))
    return tuple(segments)


def segments_from_record(record, *, include_sa: bool = True) -> Tuple[AlignmentSegment, ...]:
    """Read one pysam alignment without changing its sequence or flags."""
    if record.is_unmapped:
        return ()
    source = ("bam_secondary" if record.is_secondary else
              "bam_supplementary" if record.is_supplementary else "bam_primary")
    segment = segment_from_cigar(
        record.query_name, record.reference_name, record.reference_start,
        "-" if record.is_reverse else "+", record.cigarstring, record.mapping_quality,
        edit_distance=int(record.get_tag("NM")) if record.has_tag("NM") else None,
        source=source, secondary=bool(record.is_secondary),
    )
    segments = [segment]
    if include_sa and record.has_tag("SA"):
        segments.extend(parse_sa_tag(record.query_name, record.get_tag("SA")))
    if any(other.read_length != segment.read_length for other in segments):
        raise ValueError(f"inconsistent full query lengths for {record.query_name!r}")
    return tuple(segments)


def deduplicate_segments(segments: Iterable[AlignmentSegment]) -> Tuple[AlignmentSegment, ...]:
    """Collapse SA/BAM copies without inflating confidence or molecule support."""
    groups = {}
    for segment in segments:
        groups.setdefault(segment.geometry_key, []).append(segment)
    output = []
    for key in sorted(groups):
        copies = groups[key]
        quality = min(copies, key=lambda item: item.effective_mapping_quality).mapping_quality
        distances = [item.edit_distance for item in copies if item.edit_distance is not None]
        sources = tuple(sorted({source for item in copies for source in item.sources}))
        representative = min(copies, key=lambda item: (item.aligned_bases, item.sources))
        output.append(replace(
            representative, mapping_quality=quality,
            edit_distance=max(distances) if distances else None,
            sources=sources, secondary=all(item.secondary for item in copies),
        ))
    return tuple(output)
