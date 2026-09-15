"""Independent, complete BAM scans for fusion evidence and molecule sequences."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional, Tuple

from fin.fusion.annotation import GeneIndex, annotate_read_path
from fin.fusion.paths import ReadPath
from fin.fusion.records import AlignmentSegment, deduplicate_segments, segments_from_record


class FusionInputError(RuntimeError):
    """No fusion calls may be published from an incomplete or changing input."""


@dataclass(frozen=True)
class MoleculeEvidence:
    read_id: str
    segments: Tuple[AlignmentSegment, ...]
    sequence: Optional[str]
    status: str
    reason: str


@dataclass(frozen=True)
class FusionInput:
    molecules: Tuple[MoleculeEvidence, ...]
    record_counts: dict
    complete: bool


def _snapshot(path):
    stat = Path(path).stat()
    return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns


def _is_nomination(record, genes, min_soft_clip):
    """Which route, if any, nominates this record. Returns "" or the route."""
    if record.is_unmapped:
        return ""
    if record.is_supplementary or record.has_tag("SA"):
        return "supplementary_or_sa"
    if any(code in (4, 5) and length >= min_soft_clip for code, length in record.cigartuples or ()):
        return "soft_clip"
    if genes.genes and any(code == 3 for code, _ in record.cigartuples or ()):
        segment, = segments_from_record(record, include_sa=False)
        path = annotate_read_path(ReadPath(record.query_name, (segment,)), genes)
        if len(path.arms) < 2:
            return ""
        # A GENE TRANSITION, not merely a change of label set. `labels()`
        # attaches a gene to a block on ONE base of exonic overlap, so on a
        # real annotation a read that simply enters an overlapping gene's
        # territory ({A} -> {A,B}) splits into arms and would be nominated:
        # measured at 6.1% of records on chr9+chr22 of a real library, which
        # alone exceeds the candidate cap. A transition requires each side to
        # carry a gene the other lacks, and both sides to be annotated at all
        # -- leaving annotated territory ({A} -> {}) is novel structure, not a
        # gene change. This route only ever sees ONE collinear segment, so it
        # can only find adjacent-gene readthrough; a genuine distant fusion
        # arrives with a supplementary alignment or a large clip. Tightening
        # it takes 6.09% -> 0.19% on that measurement.
        #
        # Over ALL pairs of arms, not adjacent ones. A read crossing an
        # intermediate arm -- {A} -> {} -> {B} through unannotated sequence, or
        # {A} -> {A,B} -> {B} through the overlap of two genes -- connects two
        # genes without any single adjacent pair saying so ({A} -> {A,B} is a
        # subset step, {A,B} -> {B} a superset step). The arm count is small,
        # and a path that merely enters or leaves one gene's territory still
        # yields no qualifying pair.
        sets = [set(arm.gene_ids) for arm in path.arms]
        if any(left and right and (left - right) and (right - left)
               for index, left in enumerate(sets) for right in sets[index + 1:]):
            return "gene_transition"
    return ""


def scan_fusion_bam(
    bam_path,
    genes: Optional[GeneIndex] = None,
    *,
    min_soft_clip=250,
    max_candidate_reads=100000,
    max_segments_per_read=64,
) -> FusionInput:
    """Nominate molecules, then retrieve all records including normal competitors.

    Both scans are sequential and require no BAM index. The read cap fails the
    analysis rather than selecting a silent prefix. Bad per-read alignment data
    is unresolved; a file-level failure invalidates the entire scan.
    """
    if min_soft_clip <= 0 or min(max_candidate_reads, max_segments_per_read) <= 0:
        raise ValueError("fusion input limits must be positive")
    import logging

    import pysam

    log = logging.getLogger(__name__)
    genes = genes if genes is not None else GeneIndex()
    counts = Counter()
    nominated = set()
    invalid = {}
    try:
        identity = _snapshot(bam_path)
        with pysam.AlignmentFile(str(bam_path), "rb") as reader:
            for record in reader.fetch(until_eof=True):
                counts["alignment_records"] += 1
                if counts["alignment_records"] % 200000 == 0:
                    log.info("fusion scan pass 1: %d records, %d nominated",
                             counts["alignment_records"], len(nominated))
                if not record.is_secondary and not record.is_supplementary:
                    counts["primary_records"] += 1
                if record.is_unmapped:
                    counts["unmapped_records"] += 1
                    continue
                read_id = record.query_name
                if not read_id:
                    raise FusionInputError("mapped BAM record has no read identity")
                try:
                    route = _is_nomination(record, genes, min_soft_clip)
                except (ValueError, TypeError) as exc:
                    route = "invalid_alignment"
                    invalid[read_id] = f"invalid_alignment:{type(exc).__name__}"
                if route:
                    counts["nominated_by_" + route] += 1
                    if read_id not in nominated:
                        counts["nominated_first_by_" + route] += 1
                    nominated.add(read_id)
                    if len(nominated) > max_candidate_reads:
                        raise FusionInputError("candidate read limit exceeded; no complete fusion analysis")
        if _snapshot(bam_path) != identity:
            raise FusionInputError("BAM changed during fusion nomination")
        counts["nominated_molecules"] = len(nominated)
        segments = defaultdict(dict)
        sequences = {}
        if nominated:
            with pysam.AlignmentFile(str(bam_path), "rb") as reader:
                for record in reader.fetch(until_eof=True):
                    read_id = record.query_name
                    if read_id not in nominated or record.is_unmapped:
                        continue
                    counts["candidate_alignment_records"] += 1
                    if read_id in invalid:
                        continue
                    try:
                        parsed = segments_from_record(record)
                    except (ValueError, TypeError) as exc:
                        invalid[read_id] = f"invalid_alignment:{type(exc).__name__}"
                        segments.pop(read_id, None)
                        sequences.pop(read_id, None)
                        continue
                    stored = segments[read_id]
                    for segment in parsed:
                        if segment.geometry_key in stored:
                            stored[segment.geometry_key] = deduplicate_segments((
                                stored[segment.geometry_key], segment,
                            ))[0]
                        else:
                            stored[segment.geometry_key] = segment
                    if len(stored) > max_segments_per_read:
                        invalid[read_id] = "segment_limit"
                        segments.pop(read_id, None)
                        sequences.pop(read_id, None)
                        continue
                    lengths = {segment.read_length for segment in stored.values()}
                    if len(lengths) > 1:
                        invalid[read_id] = "inconsistent_read_length"
                        segments.pop(read_id, None)
                        sequences.pop(read_id, None)
                        continue
                    full_length = next(iter(lengths))
                    if not record.is_secondary and not any(code == 5 for code, _ in record.cigartuples or ()):
                        sequence = record.get_forward_sequence()
                        if sequence is not None:
                            sequence = sequence.upper()
                            if len(sequence) != full_length:
                                invalid[read_id] = "sequence_length_mismatch"
                            elif read_id in sequences and sequences[read_id] != sequence:
                                invalid[read_id] = "conflicting_forward_sequences"
                            else:
                                sequences[read_id] = sequence
            if _snapshot(bam_path) != identity:
                raise FusionInputError("BAM changed during fusion evidence collection")
    except FusionInputError:
        raise
    except (OSError, ValueError) as exc:
        raise FusionInputError(f"fusion BAM scan incomplete: {exc}") from exc
    molecules = []
    for read_id in sorted(nominated):
        if read_id in invalid:
            molecules.append(MoleculeEvidence(read_id, (), None, "unresolved", invalid[read_id]))
            counts["invalid_molecules"] += 1
        else:
            evidence = tuple(sorted(segments[read_id].values(), key=lambda segment: segment.geometry_key))
            sequence = sequences.get(read_id)
            molecules.append(MoleculeEvidence(
                read_id, evidence, sequence, "ready" if sequence else "needs_sequence",
                "complete_alignment_scan" if sequence else "full_sequence_unavailable",
            ))
            counts["ready_molecules" if sequence else "molecules_missing_sequence"] += 1
    return FusionInput(tuple(molecules), dict(sorted(counts.items())), True)


def supply_fusion_sequences(evidence: FusionInput, sequence_path) -> FusionInput:
    """Recover hard-clipped molecules from FASTQ/FASTA without hiding conflicts."""
    import pysam

    if not evidence.complete:
        raise FusionInputError("cannot supplement an incomplete BAM scan")
    wanted = {molecule.read_id: molecule for molecule in evidence.molecules}
    supplied = {}
    invalid = {}
    counts = Counter(evidence.record_counts)
    try:
        identity = _snapshot(sequence_path)
        with pysam.FastxFile(str(sequence_path)) as reader:
            for record in reader:
                counts["sequence_records_scanned"] += 1
                if record.name not in wanted:
                    continue
                sequence = record.sequence.upper().replace("U", "T")
                if record.name in supplied and supplied[record.name] != sequence:
                    invalid[record.name] = "conflicting_sequence_records"
                supplied[record.name] = sequence
        if _snapshot(sequence_path) != identity:
            raise FusionInputError("sequence input changed during fusion scan")
    except FusionInputError:
        raise
    except (OSError, ValueError) as exc:
        raise FusionInputError(f"fusion sequence scan incomplete: {exc}") from exc
    molecules = []
    for molecule in evidence.molecules:
        sequence = supplied.get(molecule.read_id)
        if molecule.status == "unresolved" or sequence is None:
            molecules.append(molecule)
            continue
        reason = invalid.get(molecule.read_id)
        if not reason and (not sequence or not molecule.segments
                           or any(segment.read_length != len(sequence) for segment in molecule.segments)):
            reason = "sequence_length_mismatch"
        if (not reason and molecule.sequence is not None
                and molecule.sequence.upper().replace("U", "T") != sequence):
            reason = "bam_sequence_conflict"
        if reason:
            molecules.append(replace(molecule, status="unresolved", reason=reason, sequence=None))
            counts["sequence_conflict_molecules"] += 1
        else:
            molecules.append(replace(molecule, sequence=sequence, status="ready",
                                     reason="sequence_input_validated"))
            counts["sequence_validated_molecules"] += 1
            if molecule.sequence is None:
                counts["sequence_recovered_molecules"] += 1
    counts["ready_molecules"] = sum(molecule.status == "ready" for molecule in molecules)
    counts["invalid_molecules"] = sum(molecule.status == "unresolved" for molecule in molecules)
    counts["molecules_missing_sequence"] = sum(molecule.status == "needs_sequence" for molecule in molecules)
    return FusionInput(tuple(molecules), dict(sorted(counts.items())), True)
