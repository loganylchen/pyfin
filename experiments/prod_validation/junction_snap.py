#!/usr/bin/env python3
"""Offline evaluation of annotation-free, read-supported junction snapping."""
from __future__ import annotations

import argparse
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path

import pysam

from fin.candidates.intron_chains import extract_intron_chain


TX_ID = re.compile(r'transcript_id "([^"]+)')
NUM_READS = re.compile(r'num_reads "([0-9.]+)')
ABUNDANCE = re.compile(r'abundance "([0-9.]+)')


@dataclass
class Transcript:
    transcript_line: list[str]
    exons: list[list[str]]

    @property
    def transcript_id(self) -> str:
        match = TX_ID.search(self.transcript_line[8])
        return match.group(1) if match else ""

    @property
    def score(self) -> tuple[float, float, str]:
        reads = NUM_READS.search(self.transcript_line[8])
        abundance = ABUNDANCE.search(self.transcript_line[8])
        return (
            float(reads.group(1)) if reads else 0.0,
            float(abundance.group(1)) if abundance else 0.0,
            self.transcript_id,
        )


def observed_junctions(bam_path: Path) -> dict[tuple[str, str], Counter]:
    observed: dict[tuple[str, str], Counter] = defaultdict(Counter)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or not read.cigartuples
            ):
                continue
            strand = "-" if read.is_reverse else "+"
            chain = extract_intron_chain(read.cigartuples, read.reference_start)
            for intron in chain.introns:
                observed[(read.reference_name, strand)][intron] += 1
    return observed


def read_transcripts(path: Path) -> list[Transcript]:
    transcripts: dict[str, Transcript] = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in {"transcript", "exon"}:
                continue
            match = TX_ID.search(fields[8])
            if not match:
                continue
            transcript_id = match.group(1)
            if fields[2] == "transcript":
                transcripts[transcript_id] = Transcript(fields, [])
            elif transcript_id in transcripts:
                transcripts[transcript_id].exons.append(fields)
    return list(transcripts.values())


def snap_transcript(
    transcript: Transcript,
    support: Counter,
    tolerance: int,
    min_support: int,
    min_ratio: float,
) -> tuple[Transcript, int]:
    exons = sorted((fields.copy() for fields in transcript.exons), key=lambda x: int(x[3]))
    snapped = 0
    for index in range(len(exons) - 1):
        current = (int(exons[index][4]), int(exons[index + 1][3]) - 1)
        nearby = [
            (count, intron)
            for intron, count in support.items()
            if abs(intron[0] - current[0]) <= tolerance
            and abs(intron[1] - current[1]) <= tolerance
        ]
        if not nearby:
            continue
        nearby.sort(
            key=lambda item: (
                -item[0],
                abs(item[1][0] - current[0]) + abs(item[1][1] - current[1]),
                item[1],
            )
        )
        best_count, best = nearby[0]
        current_count = support[current]
        required = max(min_support, math.floor(current_count * min_ratio) + 1)
        if best == current or best_count < required:
            continue
        if best[0] >= best[1]:
            continue
        exons[index][4] = str(best[0])
        exons[index + 1][3] = str(best[1] + 1)
        snapped += 1
    return Transcript(transcript.transcript_line.copy(), exons), snapped


def structural_key(transcript: Transcript) -> tuple:
    exons = tuple(
        (int(fields[3]), int(fields[4]))
        for fields in sorted(transcript.exons, key=lambda x: int(x[3]))
    )
    return (
        transcript.transcript_line[0],
        transcript.transcript_line[6],
        exons,
    )


def write_transcripts(path: Path, transcripts: list[Transcript]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as handle:
        for transcript in sorted(
            transcripts,
            key=lambda tx: (
                tx.transcript_line[0],
                int(tx.transcript_line[3]),
                int(tx.transcript_line[4]),
                tx.transcript_id,
            ),
        ):
            handle.write("\t".join(transcript.transcript_line) + "\n")
            for number, exon in enumerate(
                sorted(transcript.exons, key=lambda x: int(x[3])), start=1
            ):
                handle.write("\t".join(exon) + "\n")
    temporary.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--input-gtf", required=True, type=Path)
    parser.add_argument("--output-gtf", required=True, type=Path)
    parser.add_argument("--tolerance", type=int, required=True)
    parser.add_argument("--min-support", type=int, default=2)
    parser.add_argument("--min-ratio", type=float, default=1.0)
    args = parser.parse_args()

    observed = observed_junctions(args.bam)
    snapped_transcripts = 0
    snapped_junctions = 0
    representatives: dict[tuple, Transcript] = {}
    for transcript in read_transcripts(args.input_gtf):
        key = (transcript.transcript_line[0], transcript.transcript_line[6])
        snapped, count = snap_transcript(
            transcript,
            observed.get(key, Counter()),
            args.tolerance,
            args.min_support,
            args.min_ratio,
        )
        snapped_transcripts += count > 0
        snapped_junctions += count
        structure = structural_key(snapped)
        current = representatives.get(structure)
        if current is None or snapped.score > current.score:
            representatives[structure] = snapped

    write_transcripts(args.output_gtf, list(representatives.values()))
    print(
        f"input={len(read_transcripts(args.input_gtf))} "
        f"output={len(representatives)} snapped_transcripts={snapped_transcripts} "
        f"snapped_junctions={snapped_junctions}"
    )


if __name__ == "__main__":
    main()
