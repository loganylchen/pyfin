"""Fused-transcript models on top of molecule-level fusion chains.

The chain layer (`chains.py`) answers "which ordered gene path and breakpoint
tuple do complete molecules support?". It does not say WHICH TRANSCRIPT of
each gene took part, does not extend a 5'-degraded direct-RNA read to a
transcript start, and keys its identity on observed coordinate modes. This
module adds the transcript level without touching any of that:

* a transcript-preserving annotation index (``TranscriptIndex``) -- the gene
  index deliberately merges isoforms, so it cannot resolve them;
* per-arm resolution against annotated transcripts on the arm's OBSERVED
  introns, with arm-relative classes (``arm_fsm``, ``arm_ism_5p`` ...) that
  keep splice identity separate from endpoint completeness;
* one ``FusedTranscriptModel`` per distinct observed structure, identified by
  its ordered exons -- never by support counts, transcript ids or chain ids;
* truncation-aware assignment: a molecule that observed a 5' suffix of a model
  is COMPATIBLE with it, and gets a unique vote only if it is compatible with
  exactly one model. Compatible annotation never becomes observation;
* labelled 5'/3' completion hypotheses (``annotated_extension``) as separate
  products from the evidence-constrained model.

Everything here reads the finished ``FusionAnalysis``; it never feeds back
into chain support, correction or signal adjudication.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from dataclasses import dataclass, field, replace
from typing import Dict, Iterable, List, Mapping, Optional, Tuple

from fin.fusion.annotation import AnnotatedReadPath, FusionArm
from fin.utils.sequences import reverse_complement

Exon = Tuple[int, int]


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


# --------------------------------------------------------------------------
# annotation, isoform-preserving
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class TranscriptModel:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    exons: Tuple[Exon, ...]  # genomic order, 0-based half-open, sorted

    def __post_init__(self):
        if not self.exons or self.strand not in ("+", "-"):
            raise ValueError("transcript model needs exons and a strand")
        if any(end <= start for start, end in self.exons):
            raise ValueError("transcript exon must be non-empty")
        if list(self.exons) != sorted(self.exons):
            raise ValueError("transcript exons must be sorted")
        if any(following[0] < previous[1] for previous, following in zip(self.exons, self.exons[1:])):
            raise ValueError("transcript exons must not overlap")

    @property
    def introns(self):
        """Genomic introns as (end of previous exon, start of next)."""
        return tuple((previous[1], following[0]) for previous, following in zip(self.exons, self.exons[1:]))

    @property
    def tss(self):
        return self.exons[0][0] if self.strand == "+" else self.exons[-1][1]

    @property
    def tes(self):
        return self.exons[-1][1] if self.strand == "+" else self.exons[0][0]

    def exon_rank(self, exon):
        """1-based rank in transcript orientation of a genomic exon."""
        index = self.exons.index(exon)
        return index + 1 if self.strand == "+" else len(self.exons) - index


class TranscriptIndex:
    """Bin index over annotated transcripts that keeps every isoform intact."""

    def __init__(self, transcripts: Iterable[TranscriptModel] = (), *, bin_size: int = 16384):
        if bin_size <= 0:
            raise ValueError("transcript-index bin size must be positive")
        self.bin_size = bin_size
        self.transcripts = tuple(sorted(set(transcripts), key=lambda t: (t.chrom, t.strand, t.gene_id, t.transcript_id)))
        self._bins: Dict[Tuple[str, str, int], set] = defaultdict(set)
        for index, transcript in enumerate(self.transcripts):
            for start, end in transcript.exons:
                for bucket in range(start // bin_size, (end - 1) // bin_size + 1):
                    self._bins[(transcript.chrom, transcript.strand, bucket)].add(index)

    @classmethod
    def from_transcripts(cls, transcripts):
        models = []
        for transcript in transcripts:
            if not transcript.gene_id or transcript.strand not in ("+", "-") or not transcript.exons:
                continue
            attributes = getattr(transcript, "attributes", None) or {}
            models.append(TranscriptModel(
                transcript.transcript_id, transcript.gene_id,
                str(attributes.get("gene_name") or ""), transcript.chrom, transcript.strand,
                tuple(sorted((int(start), int(end)) for start, end in transcript.exons)),
            ))
        return cls(models)

    def overlapping(self, chrom, strand, exons):
        candidates = set()
        for start, end in exons:
            for bucket in range(start // self.bin_size, (end - 1) // self.bin_size + 1):
                candidates.update(self._bins.get((chrom, strand, bucket), ()))
        return tuple(self.transcripts[index] for index in sorted(candidates))


# --------------------------------------------------------------------------
# per-arm resolution
# --------------------------------------------------------------------------

ARM_FIRST, ARM_MIDDLE, ARM_LAST = "first", "middle", "last"


def _oriented(exons, strand):
    """Exons in transcript orientation (5' -> 3')."""
    return tuple(exons) if strand == "+" else tuple(reversed(tuple(exons)))


def _retained_slice(transcript: TranscriptModel, role, entry, exit_):
    """The part of an annotated transcript a fusion arm could retain.

    ``entry``/``exit_`` are the arm's incoming/outgoing junction-facing
    positions (genomic coordinates, half-open convention of the arm blocks).
    Returns exons in GENOMIC order, or None when the breakpoint does not fall
    inside/at the boundary of this transcript's exonic sequence.
    """
    plus = transcript.strand == "+"
    if role == ARM_FIRST:
        lo, hi = (transcript.tss, exit_) if plus else (exit_, transcript.tss)
    elif role == ARM_LAST:
        lo, hi = (entry, transcript.tes) if plus else (transcript.tes, entry)
    else:
        lo, hi = (entry, exit_) if plus else (exit_, entry)
    if hi <= lo:
        return None
    kept = []
    for start, end in transcript.exons:
        s, e = max(start, lo), min(end, hi)
        if e > s:
            kept.append((s, e))
    if not kept:
        return None
    # the breakpoint-facing boundary must be exonic in this transcript, i.e.
    # the slice must actually reach it (a breakpoint inside an intron of T
    # means T cannot be the retained isoform)
    if role in (ARM_FIRST, ARM_MIDDLE):
        reach = kept[-1][1] if plus else kept[0][0]
        if reach != exit_:
            return None
    if role in (ARM_MIDDLE, ARM_LAST):
        reach = kept[0][0] if plus else kept[-1][1]
        if reach != entry:
            return None
    return tuple(kept)


@dataclass(frozen=True)
class ArmResolution:
    role: str
    match_class: str
    compatible_transcript_ids: Tuple[str, ...]
    gene_ids: Tuple[str, ...]
    gene_names: Tuple[str, ...]
    # annotated completion in GENOMIC order for the best-supported class, or ()
    completion_exons: Tuple[Exon, ...]
    completion_transcript_ids: Tuple[str, ...]
    # exon ranks (transcript orientation, 1-based) of the exon that FACES the
    # outgoing breakpoint (exit) and the one that RECEIVES the incoming one
    # (entry), per compatible transcript -- what "e13a2" is made of
    breakpoint_exon_ranks: Tuple[Tuple[str, int], ...]
    entry_exon_ranks: Tuple[Tuple[str, int], ...] = ()


def genomic_exons(arm: FusionArm):
    """A FusionArm's blocks follow TRANSCRIPT order (descending on minus);
    every comparison against annotation below is done in ascending genomic
    order, so this is the one place the conversion happens."""
    return tuple(sorted((block.reference_start, block.reference_end) for block in arm.blocks))


def introns_of(exons_ascending):
    return tuple((p[1], f[0]) for p, f in zip(exons_ascending, exons_ascending[1:]))


def _annotation_intron_overlaps_exon(transcript, exons_ascending):
    """An annotated intron overlapping an observed exon means the molecule
    aligned continuously across a splice site of that transcript: incompatible.
    Abutting at an exon boundary (half-open) is not an overlap."""
    for start, end in exons_ascending:
        for i_start, i_end in transcript.introns:
            if i_start < end and start < i_end:
                return True
    return False


def resolve_arm(chrom, strand, exons_ascending, role, index: TranscriptIndex, gene_ids=()) -> ArmResolution:
    """Match one arm (ascending genomic exons) against annotated isoforms.

    Compatibility (all required): every observed intron is an intron of the
    transcript; no annotated intron overlaps an observed exon; the
    breakpoint-facing boundary is exonic in the transcript; every observed
    exon lies inside the retained annotated slice. Missing annotated introns
    are tolerated ONLY outside the arm's observed coverage.

    ALL compatible transcripts are kept and reported; ``match_class`` is the
    tightest class among them (display), and a completion is offered only
    when every compatible transcript retains the same exons -- a ranking is
    never used as evidence against a longer compatible isoform.
    """
    exons = tuple(exons_ascending)
    observed = introns_of(exons)
    plus = strand == "+"
    entry = exons[0][0] if plus else exons[-1][1]
    exit_ = exons[-1][1] if plus else exons[0][0]
    cover_lo, cover_hi = exons[0][0], exons[-1][1]
    candidates = index.overlapping(chrom, strand, exons) if index.transcripts else ()
    compatible: List[Tuple[TranscriptModel, Tuple[Exon, ...], str]] = []
    for transcript in candidates:
        if gene_ids and transcript.gene_id not in gene_ids:
            continue
        introns = set(transcript.introns)
        if any(intron not in introns for intron in observed):
            continue
        if _annotation_intron_overlaps_exon(transcript, exons):
            continue
        retained = _retained_slice(transcript, role, entry, exit_)
        if retained is None:
            continue
        if not all(any(s <= start and end <= e for s, e in retained) for start, end in exons):
            continue
        missing = [i for i in introns_of(retained) if i not in set(observed)]
        if not missing:
            klass = "arm_fsm"
        else:
            before = any(i[1] <= cover_lo for i in missing)
            after = any(i[0] >= cover_hi for i in missing)
            five_missing, three_missing = (before, after) if plus else (after, before)
            klass = ("arm_ism_both" if five_missing and three_missing else
                     "arm_ism_5p" if five_missing else "arm_ism_3p")
        compatible.append((transcript, retained, klass))
    if not compatible:
        gene_transcripts = [t for t in candidates if not gene_ids or t.gene_id in gene_ids]
        if not candidates and not gene_ids:
            klass = "unannotated_locus"
        elif observed and gene_transcripts and any(
                intron not in {i for t in gene_transcripts for i in t.introns} for intron in observed):
            klass = "novel_junction"
        elif observed:
            klass = "novel_combination_of_known_junctions"
        else:
            klass = "single_exon_unresolved"
        return ArmResolution(role, klass, (), tuple(gene_ids), (), (), (), ())
    compatible.sort(key=lambda item: item[0].transcript_id)
    order = {"arm_fsm": 0, "arm_ism_3p": 1, "arm_ism_5p": 1, "arm_ism_both": 2}
    tightest = min(order[k] for _, _, k in compatible)
    if not observed and len(exons) == 1:
        match_class = "single_exon_compatible"
    else:
        classes = {k for _, _, k in compatible if order[k] == tightest}
        match_class = next(iter(classes)) if len(classes) == 1 else "arm_ism_mixed"
    structures = {r for _, r, _ in compatible}
    completion = next(iter(structures)) if len(structures) == 1 else ()
    exit_ranks, entry_ranks = [], []
    for transcript, retained, _ in compatible:
        if role in (ARM_FIRST, ARM_MIDDLE):
            facing = retained[-1] if plus else retained[0]
            full = next((e for e in transcript.exons if e[0] <= facing[0] and facing[1] <= e[1]), None)
            if full is not None:
                exit_ranks.append((transcript.transcript_id, transcript.exon_rank(full)))
        if role in (ARM_MIDDLE, ARM_LAST):
            receiving = retained[0] if plus else retained[-1]
            full = next((e for e in transcript.exons if e[0] <= receiving[0] and receiving[1] <= e[1]), None)
            if full is not None:
                entry_ranks.append((transcript.transcript_id, transcript.exon_rank(full)))
    return ArmResolution(
        role, match_class, tuple(t.transcript_id for t, _, _ in compatible),
        tuple(sorted({t.gene_id for t, _, _ in compatible})),
        tuple(sorted({t.gene_name for t, _, _ in compatible if t.gene_name})),
        completion, tuple(t.transcript_id for t, _, _ in compatible) if completion else (),
        tuple(exit_ranks), tuple(entry_ranks),
    )


# --------------------------------------------------------------------------
# models
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class ModelArm:
    chrom: str
    strand: str
    exons: Tuple[Exon, ...]            # ASCENDING genomic order, observed (junction-corrected) coordinates
    resolution: ArmResolution


@dataclass(frozen=True)
class FusedTranscriptModel:
    model_id: str
    arms: Tuple[ModelArm, ...]
    breakpoints: Tuple[Tuple[str, int, str, str, int, str], ...]  # (chromL, exitL, strandL, chromR, entryR, strandR)
    primary_chain_id: str
    source_chain_ids: Tuple[str, ...]
    unique_read_ids: Tuple[str, ...]          # compatible with THIS reported model only
    compatible_read_ids: Tuple[str, ...]      # compatible with this and other reported models
    partial_read_ids: Tuple[str, ...]         # molecules covering only some arms (chain partial support)
    chain_edge_support: Tuple[int, ...]       # per edge, from the parent chain (includes full-support molecules)
    coordinate_status: str
    event_class: str
    alternative_group_id: str
    five_prime_status: str
    three_prime_status: str
    five_prime_observed_read_ids: Tuple[str, ...] = ()   # molecules whose start reaches the model's 5' extent
    three_prime_observed_read_ids: Tuple[str, ...] = ()
    # other placements this model's OWN supporting molecules hold as
    # alternatives: complete exon structures (JSON strings of arms), not just
    # breakpoints -- an internal acceptor 3 bp away leaves the fusion edge
    # unchanged and would be invisible in a breakpoint list
    coordinate_alternatives: Tuple[str, ...] = ()
    # molecules whose own decisive placement differs from the envelope by no
    # more than the coordinate tolerance: NOT members, NOT support -- reported
    # so a reader can see them, with their placements (JSON strings of arms)
    nearby_read_ids: Tuple[str, ...] = ()
    nearby_observed_placements: Tuple[str, ...] = ()
    # "single_molecule_spanned" when at least one member reaches both ends of
    # the envelope, else "aggregate": the envelope is a compatibility union
    # of pieces observed on different molecules, not one observed transcript
    envelope_status: str = "aggregate"

    @property
    def gene_path(self):
        return tuple(arm.resolution.gene_ids or () for arm in self.arms)

    @property
    def structure_key(self):
        return structure_key_of(self.arms)


@dataclass(frozen=True)
class CollapsedStructure:
    """An observed structure that is a strict 5' truncation of a reported
    model: not silently gone, recorded with the models it could be."""
    splice_summary: str
    read_ids: Tuple[str, ...]
    compatible_model_ids: Tuple[str, ...]
    interpretation: str = ("5'-degraded observation of a reported model OR a genuine shorter isoform "
                           "with a downstream start; direct RNA cannot distinguish these")


def structure_key_of(arms):
    """Identity of a fused transcript: ordered arms with explicit exons."""
    return tuple((arm.chrom, arm.strand, tuple(arm.exons)) for arm in arms)


@dataclass(frozen=True)
class _Observation:
    """One selected path of one molecule, in the ascending-exon convention."""
    read_id: str
    arms: Tuple[Tuple[str, str, Tuple[Exon, ...], Tuple[str, ...]], ...]  # chrom, strand, exons asc, gene_ids
    event_class: str

    @property
    def splice_key(self):
        """Everything observed EXCEPT the outer 5' and 3' ends."""
        parts = []
        n = len(self.arms)
        for index, (chrom, strand, exons, _) in enumerate(self.arms):
            plus = strand == "+"
            entry = (exons[0][0] if plus else exons[-1][1]) if index > 0 else None
            exit_ = (exons[-1][1] if plus else exons[0][0]) if index < n - 1 else None
            parts.append((chrom, strand, entry, introns_of(exons), exit_))
        return tuple(parts)


def _observation(read_id, path: AnnotatedReadPath):
    return _Observation(read_id, tuple((arm.chrom, arm.strand, genomic_exons(arm), tuple(arm.gene_ids))
                                       for arm in path.arms), path.event_class)


def _close(x, y, tolerance):
    return abs(x - y) <= tolerance


def _same_exons(a, b, tolerance):
    return len(a) == len(b) and all(_close(s1, s2, tolerance) and _close(e1, e2, tolerance)
                                    for (s1, e1), (s2, e2) in zip(a, b))


def _suffix_fit(short_exons, long_exons, strand, tolerance=0):
    """`short` is a 5' truncation of `long` (transcript orientation): its
    oriented exons are a suffix of the longer's, exact (within `tolerance`)
    for junction-bounded blocks, and its 5'-most exon lies inside the
    corresponding longer exon while ending at the same 3'-facing boundary."""
    sa, sb = _oriented(short_exons, strand), _oriented(long_exons, strand)
    if len(sa) > len(sb):
        return False
    tail = sb[len(sb) - len(sa):]
    if not _same_exons(sa[1:], tail[1:], tolerance):
        return False
    head, ref = sa[0], tail[0]
    if strand == "+":
        return _close(head[1], ref[1], tolerance) and ref[0] <= head[0] + tolerance
    return _close(head[0], ref[0], tolerance) and head[1] <= ref[1] + tolerance


def _prefix_fit(short_exons, long_exons, strand, tolerance=0):
    """`short` is a 3' truncation of `long`: its oriented exons are a prefix
    of the longer's, exact (within `tolerance`) except that its 3'-most exon
    may end early inside the corresponding longer exon."""
    sa, sb = _oriented(short_exons, strand), _oriented(long_exons, strand)
    if len(sa) > len(sb):
        return False
    head = sb[:len(sa)]
    if not _same_exons(sa[:-1], head[:-1], tolerance):
        return False
    last, ref = sa[-1], head[-1]
    if strand == "+":
        return _close(last[0], ref[0], tolerance) and last[1] <= ref[1] + tolerance
    return _close(last[1], ref[1], tolerance) and ref[0] <= last[0] + tolerance


def _fits_nearby(inner_arms, outer_arms, tolerance):
    """`inner` fits inside `outer` under the truncation relation with every
    coordinate compared within `tolerance` (used only to fold placements the
    evidence could not separate; see _merge_unresolved_alternatives)."""
    if len(inner_arms) != len(outer_arms):
        return False
    n = len(outer_arms)
    for k, (x, y) in enumerate(zip(inner_arms, outer_arms)):
        if (x.chrom, x.strand) != (y.chrom, y.strand):
            return False
        if 0 < k < n - 1:
            if not _same_exons(x.exons, y.exons, tolerance):
                return False
        elif k == 0:
            exit_x = x.exons[-1][1] if x.strand == "+" else x.exons[0][0]
            exit_y = y.exons[-1][1] if y.strand == "+" else y.exons[0][0]
            if not _close(exit_x, exit_y, tolerance) or not _suffix_fit(x.exons, y.exons, x.strand, tolerance):
                return False
        else:
            entry_x = x.exons[0][0] if x.strand == "+" else x.exons[-1][1]
            entry_y = y.exons[0][0] if y.strand == "+" else y.exons[-1][1]
            if not _close(entry_x, entry_y, tolerance) or not _prefix_fit(x.exons, y.exons, x.strand, tolerance):
                return False
    return True


def _fits_within(shorter: _Observation, envelope_arms):
    """True when `shorter` observed the same fusion and internal splicing as
    the envelope but less of the outer arms: a 5' truncation of the FIRST
    arm (direct-RNA degradation) and/or a 3' truncation of the LAST arm (a
    molecule that did not reach the end). Every middle arm and every
    breakpoint must be identical. Judged on THIS molecule's own coordinates."""
    if len(shorter.arms) != len(envelope_arms):
        return False
    n = len(envelope_arms)
    for k, ((chrom, strand, exons, _), env) in enumerate(zip(shorter.arms, envelope_arms)):
        if (chrom, strand) != (env.chrom, env.strand):
            return False
        if 0 < k < n - 1:
            if exons != env.exons:
                return False
            continue
        first, last = k == 0, k == n - 1
        if first and last:
            # a single-arm path is never a fusion; unreachable, but keep the
            # rule total: both ends may be truncated
            ok = any(_suffix_fit(exons, env.exons, strand) and _prefix_fit(exons, env.exons, strand)
                     for _ in (0,))
        elif first:
            # exit boundary must agree (a junction-facing coordinate)
            exit_s = exons[-1][1] if strand == "+" else exons[0][0]
            exit_e = env.exons[-1][1] if strand == "+" else env.exons[0][0]
            ok = exit_s == exit_e and _suffix_fit(exons, env.exons, strand)
        else:
            entry_s = exons[0][0] if strand == "+" else exons[-1][1]
            entry_e = env.exons[0][0] if strand == "+" else env.exons[-1][1]
            ok = entry_s == entry_e and _prefix_fit(exons, env.exons, strand)
        if not ok:
            return False
    return True


def _union_envelope(a_arms, b_arms):
    """The aggregate envelope of two compatible structures, or None.

    Two structures with identical middle arms and breakpoints, whose first
    arms stand in a 5'-truncation relation and whose last arms stand in a
    3'-truncation relation, describe one fused transcript observed to
    different extents; the envelope takes the longer end on each side. No
    single molecule need have observed the whole envelope -- that is why the
    model records which molecules reach each end.
    """
    if len(a_arms) != len(b_arms):
        return None
    n = len(a_arms)
    out = []
    for k, (x, y) in enumerate(zip(a_arms, b_arms)):
        if (x.chrom, x.strand) != (y.chrom, y.strand):
            return None
        if 0 < k < n - 1:
            if x.exons != y.exons:
                return None
            out.append(x)
            continue
        exons = None
        if k == 0:
            exit_x = x.exons[-1][1] if x.strand == "+" else x.exons[0][0]
            exit_y = y.exons[-1][1] if y.strand == "+" else y.exons[0][0]
            if exit_x != exit_y:
                return None
            if _suffix_fit(x.exons, y.exons, x.strand):
                exons = y.exons
            elif _suffix_fit(y.exons, x.exons, x.strand):
                exons = x.exons
        if k == n - 1 and exons is None:
            entry_x = x.exons[0][0] if x.strand == "+" else x.exons[-1][1]
            entry_y = y.exons[0][0] if y.strand == "+" else y.exons[-1][1]
            if entry_x != entry_y:
                return None
            if _prefix_fit(x.exons, y.exons, x.strand):
                exons = y.exons
            elif _prefix_fit(y.exons, x.exons, x.strand):
                exons = x.exons
        if exons is None:
            return None
        out.append(ModelArm(x.chrom, x.strand, tuple(exons), x.resolution))
    return tuple(out)


def _breakpoints_of(arms):
    out = []
    for l, r in zip(arms, arms[1:]):
        exit_ = l.exons[-1][1] if l.strand == "+" else l.exons[0][0]
        entry = r.exons[0][0] if r.strand == "+" else r.exons[-1][1]
        out.append((l.chrom, exit_, l.strand, r.chrom, entry, r.strand))
    return tuple(out)


def _end_statuses(arms):
    """(five, three) annotation-completeness, judged on THESE arms.

    A regenerated envelope reaches further than any member did, so the members'
    statuses cannot be carried over: the extended end may no longer coincide
    with an annotated boundary, and claiming it does would be a false claim of
    observation.
    """
    five, three = arms[0].resolution, arms[-1].resolution
    return (("observed_boundary_5p_completeness_unknown" if not five.completion_exons
             else "annotated_boundary_observed" if five.completion_exons == arms[0].exons
             else "annotated_extension_available"),
            ("observed_end_completeness_unknown" if not three.completion_exons
             else "annotated_boundary_observed" if three.completion_exons == arms[-1].exons
             else "annotated_extension_available"))


def _model_id(arms):
    return "FTM_" + hashlib.sha256(_json(structure_key_of(arms)).encode()).hexdigest()[:16]


def placement_json(arms):
    """A complete placement: every arm's chrom, strand and exons."""
    return _json([[arm.chrom, arm.strand, [list(e) for e in arm.exons]] for arm in arms])


def _envelope(observations, index: TranscriptIndex):
    """Model arms for one splice structure: exons shared by every observation
    of it, with the outer exons widened to the maximal observed extent."""
    first = observations[0]
    n = len(first.arms)
    arms = []
    for k, (chrom, strand, exons, gene_ids) in enumerate(first.arms):
        exons = list(exons)
        if k == 0:
            if strand == "+":
                exons[0] = (min(o.arms[0][2][0][0] for o in observations), exons[0][1])
            else:
                exons[-1] = (exons[-1][0], max(o.arms[0][2][-1][1] for o in observations))
        if k == n - 1:
            if strand == "+":
                exons[-1] = (exons[-1][0], max(o.arms[-1][2][-1][1] for o in observations))
            else:
                exons[0] = (min(o.arms[-1][2][0][0] for o in observations), exons[0][1])
        role = ARM_FIRST if k == 0 else ARM_LAST if k == n - 1 else ARM_MIDDLE
        exons = tuple(exons)
        arms.append(ModelArm(chrom, strand, exons, resolve_arm(chrom, strand, exons, role, index, gene_ids)))
    return tuple(arms)


def _resolved(arms, index, gene_ids_per_arm):
    n = len(arms)
    out = []
    for k, arm in enumerate(arms):
        role = ARM_FIRST if k == 0 else ARM_LAST if k == n - 1 else ARM_MIDDLE
        out.append(ModelArm(arm.chrom, arm.strand, arm.exons,
                            resolve_arm(arm.chrom, arm.strand, arm.exons, role, index, gene_ids_per_arm[k])))
    return tuple(out)


def build_fused_transcript_models(analysis, transcript_index: TranscriptIndex,
                                  genome: Optional[Mapping[str, str]] = None):
    """Turn molecule-level fusion evidence into explicit fused-transcript models.

    Returns ``(models, collapsed)``.

    Observed structures (splice keys: everything except the outer 5'/3'
    ends) are collected from every chain's full-support molecules -- the
    chain layer keys on exact topology, so five molecules of one fusion that
    lost different numbers of 5' exons to degradation are five chains of
    support one there. Here they are one transcript: structures whose
    middle arms and breakpoints agree and whose outer arms stand in a
    5'-truncation / 3'-truncation relation are merged into an aggregate
    ENVELOPE that takes the longer extent on each side. No single molecule
    need have observed the whole envelope; the model records which molecules
    reach each end.

    A model is reported when at least ``settings.min_support`` distinct
    molecules fit inside its envelope (their complete observed path, judged
    on their own coordinates -- a molecule may hold several selected paths
    and each is judged separately; votes are per molecule). A molecule is
    UNIQUE to a model when it fits that model alone among reported models,
    otherwise COMPATIBLE on each it fits. One molecule, at most one unique
    vote. Structures that fit only inside other envelopes are returned as
    collapsed, with the models they could be, never silently dropped.

    Same-sequence placements of one junction (microhomology) merge into one
    record that lists the other placements; unrelated events never merge.
    """
    reads = {read.read_id: read for read in analysis.reads}
    chain_by_id = {chain.chain_id: chain for chain in analysis.chains}
    min_support = max(1, int(getattr(analysis.settings, "min_support", 2)))
    structures: Dict[tuple, dict] = {}
    for chain in analysis.chains:
        for read_id in chain.full_support_read_ids:
            audit = reads.get(read_id)
            if audit is None:
                continue
            for path in audit.selected:
                if len(path.arms) < 2:
                    continue
                observation = _observation(read_id, path)
                entry = structures.setdefault(observation.splice_key, {"obs": [], "chains": set()})
                entry["obs"].append(observation)
                entry["chains"].add(chain.chain_id)
    if not structures:
        return (), ()
    keys = sorted(structures, key=_json)
    for key in keys:
        structures[key]["obs"].sort(key=lambda o: (o.read_id, _json(o.arms)))
        structures[key]["arms"] = _envelope(structures[key]["obs"], transcript_index)
    # aggregate envelopes: merge structures that describe one transcript
    # observed to different extents; deterministic (sorted keys, first fit)
    envelopes: List[dict] = [{"arms": structures[k]["arms"], "keys": [k]} for k in keys]
    merged = True
    while merged:
        merged = False
        for a in range(len(envelopes)):
            for b in range(a + 1, len(envelopes)):
                union = _union_envelope(envelopes[a]["arms"], envelopes[b]["arms"])
                if union is not None:
                    envelopes[a] = {"arms": union, "keys": envelopes[a]["keys"] + envelopes[b]["keys"]}
                    del envelopes[b]
                    merged = True
                    break
            if merged:
                break
    for env in envelopes:
        gene_ids = tuple(structures[env["keys"][0]]["obs"][0].arms[k][3] for k in range(len(env["arms"])))
        env["arms"] = _resolved(env["arms"], transcript_index, gene_ids)
        env["id"] = _model_id(env["arms"])
    # membership per OBSERVATION, judged on its own coordinates, EXACT only:
    # coordinate proximity is not evidence, and a molecule whose corrector
    # settled on one placement must not acquire another because it is close.
    # Molecules whose decisive placement lies within `breakpoint_tolerance`
    # of an envelope are recorded on it as `nearby` -- visible, not counted.
    tolerance = int(getattr(analysis.settings, "breakpoint_tolerance", 6))
    fits: Dict[Tuple[tuple, int], set] = {}
    nearby: Dict[str, dict] = defaultdict(dict)   # env id -> read id -> placement json
    for key in keys:
        for position, obs in enumerate(structures[key]["obs"]):
            exact = {e["id"] for e in envelopes if _fits_within(obs, e["arms"])}
            fits[(key, position)] = exact
            as_arms = tuple(ModelArm(c, st, ex, None) for c, st, ex, _ in obs.arms)
            for e in envelopes:
                if e["id"] not in exact and _fits_nearby(as_arms, e["arms"], tolerance):
                    nearby[e["id"]][obs.read_id] = _json([[c, st, list(map(list, ex))] for c, st, ex, _ in obs.arms])
    molecules_of: Dict[str, set] = defaultdict(set)
    for key in keys:
        for position, obs in enumerate(structures[key]["obs"]):
            for env_id in fits[(key, position)]:
                molecules_of[env_id].add(obs.read_id)
    eligible = {e["id"] for e in envelopes if len(molecules_of[e["id"]]) >= min_support}
    read_models: Dict[str, set] = defaultdict(set)
    for key in keys:
        for position, obs in enumerate(structures[key]["obs"]):
            read_models[obs.read_id].update(m for m in fits[(key, position)] if m in eligible)
    models: List[FusedTranscriptModel] = []
    env_by_id = {e["id"]: e for e in envelopes}
    observations_of: Dict[str, list] = {}
    for env in envelopes:
        if env["id"] not in eligible:
            continue
        arms = env["arms"]
        member_keys = env["keys"]
        obs = [o for k in member_keys for o in structures[k]["obs"]
               if env["id"] in fits[(k, structures[k]["obs"].index(o))]]
        unique = sorted(r for r, ms in read_models.items() if ms == {env["id"]})
        compatible = sorted(r for r, ms in read_models.items() if env["id"] in ms and len(ms) > 1)
        chains = sorted({c for k in member_keys for c in structures[k]["chains"]})
        # display parent: the chain with the most full support, then id
        primary = max(chains, key=lambda c: (len(chain_by_id[c].full_support_read_ids), c))
        chain = chain_by_id[primary]
        first, last = arms[0], arms[-1]
        five_reach = first.exons[0][0] if first.strand == "+" else first.exons[-1][1]
        three_reach = last.exons[-1][1] if last.strand == "+" else last.exons[0][0]
        five_obs = sorted({o.read_id for o in obs if (o.arms[0][2][0][0] if first.strand == "+" else o.arms[0][2][-1][1]) == five_reach})
        three_obs = sorted({o.read_id for o in obs if (o.arms[-1][2][-1][1] if last.strand == "+" else o.arms[-1][2][0][0]) == three_reach})
        five_status, three_status = _end_statuses(arms)
        siblings = {m for r in compatible for m in read_models[r]} | {env["id"]}
        group = "CAG_" + hashlib.sha256(_json(sorted(siblings)).encode()).hexdigest()[:12]
        partial = sorted({r for c in chains for r in chain_by_id[c].partial_support_read_ids})
        near_ids = {r for r in nearby.get(env["id"], {}) if r not in unique and r not in compatible}
        models.append(FusedTranscriptModel(
            env["id"], arms, _breakpoints_of(arms), primary, tuple(chains),
            tuple(unique), tuple(compatible), tuple(partial),
            tuple(len(ids) for ids in chain.edge_support_read_ids),
            chain.coordinate_status, structures[member_keys[0]]["obs"][0].event_class, group,
            five_status, three_status, tuple(five_obs), tuple(three_obs),
            nearby_read_ids=tuple(sorted(near_ids)),
            nearby_observed_placements=tuple(sorted(nearby[env["id"]][r] for r in near_ids)),
            envelope_status=("single_molecule_spanned" if set(five_obs) & set(three_obs) else "aggregate"),
        ))
        observations_of[env["id"]] = obs
    final_id: Dict[str, str] = {m.model_id: m.model_id for m in models}
    # the specific claim first (identical sequence: microhomology), then the
    # weaker one (placements the SAME molecules hold as alternatives). After
    # each pass aliases are resolved transitively and memberships remapped,
    # so a later pass can retire an earlier survivor without dangling ids.
    if genome is not None:
        models, final_id = _merge_same_sequence_placements(models, read_models, genome, final_id)
        _canonicalise(final_id, read_models)
    models, final_id = _merge_unresolved_alternatives(models, read_models, final_id, tolerance,
                                                     transcript_index, genome)
    _canonicalise(final_id, read_models)
    models = [_recompute_endpoints(m, observations_of, final_id) for m in models]
    models.sort(key=lambda m: (-len(m.unique_read_ids), -len(m.compatible_read_ids), m.model_id))
    emitted = {m.model_id for m in models}
    assert set(final_id.values()) <= emitted, "collapsed audit would reference a removed model"
    collapsed = []
    reported_keys = {k for e in envelopes if e["id"] in eligible for k in e["keys"]
                     if structures[k]["arms"] == e["arms"] or _model_id(structures[k]["arms"]) == e["id"]}
    for key in keys:
        if key in reported_keys:
            continue
        targets = sorted({final_id[m] for position in range(len(structures[key]["obs"]))
                          for m in fits[(key, position)] if m in final_id})
        if not targets:
            # fits no reported model: an unsupported structure, already
            # visible at chain level as insufficient support; nothing was
            # folded, so there is nothing to audit here
            continue
        summary = ";".join(f"{c}:{s}:" + ",".join(f"{a}-{b}" for a, b in e) for c, s, e, _ in structures[key]["obs"][0].arms)
        collapsed.append(CollapsedStructure(summary, tuple(sorted({o.read_id for o in structures[key]["obs"]})), tuple(targets)))
    return tuple(models), tuple(collapsed)


def _nearby(a: FusedTranscriptModel, b: FusedTranscriptModel, tolerance):
    """One structure whose exact placements the sequence corrector could not
    settle: the two could form one envelope (middle arms equal, first arms in
    a 5'-truncation relation either way, last arms in a 3'-truncation
    relation either way) with every coordinate within `tolerance`. Judged per
    end independently, so a member that is 5'-longer but 3'-shorter than
    another still counts as the same population."""
    x_arms, y_arms = a.arms, b.arms
    if len(x_arms) != len(y_arms):
        return False
    n = len(x_arms)
    for k, (x, y) in enumerate(zip(x_arms, y_arms)):
        if (x.chrom, x.strand) != (y.chrom, y.strand):
            return False
        if 0 < k < n - 1:
            if not _same_exons(x.exons, y.exons, tolerance):
                return False
            continue
        if k == 0:
            exit_x = x.exons[-1][1] if x.strand == "+" else x.exons[0][0]
            exit_y = y.exons[-1][1] if y.strand == "+" else y.exons[0][0]
            if not _close(exit_x, exit_y, tolerance):
                return False
            if not (_suffix_fit(x.exons, y.exons, x.strand, tolerance) or _suffix_fit(y.exons, x.exons, x.strand, tolerance)):
                return False
        if k == n - 1:
            entry_x = x.exons[0][0] if x.strand == "+" else x.exons[-1][1]
            entry_y = y.exons[0][0] if y.strand == "+" else y.exons[-1][1]
            if not _close(entry_x, entry_y, tolerance):
                return False
            if not (_prefix_fit(x.exons, y.exons, x.strand, tolerance) or _prefix_fit(y.exons, x.exons, x.strand, tolerance)):
                return False
    return True


def _annotation_agreement(model: FusedTranscriptModel):
    return sum(1 for arm in model.arms if arm.resolution.compatible_transcript_ids)


def _canonicalise(final_id, read_models):
    """Resolve alias chains (old -> intermediate -> survivor) to their root
    and rewrite memberships in surviving ids."""
    for key in list(final_id):
        target = final_id[key]
        seen = {key}
        while final_id.get(target, target) != target and target not in seen:
            seen.add(target)
            target = final_id[target]
        final_id[key] = target
    for read_id, ms in list(read_models.items()):
        read_models[read_id] = {final_id.get(m, m) for m in ms}


def _recompute_endpoints(model, observations_of, final_id):
    """Endpoint-reach fields must be judged against the FINAL envelope: after
    a merge, a member of a retired model does not reach the survivor's ends
    just because it reached its own."""
    obs = [o for mid, items in observations_of.items() if final_id.get(mid, mid) == model.model_id for o in items]
    first, last = model.arms[0], model.arms[-1]
    five_reach = first.exons[0][0] if first.strand == "+" else first.exons[-1][1]
    three_reach = last.exons[-1][1] if last.strand == "+" else last.exons[0][0]
    reaches_five = lambda o: (o.arms[0][2][0][0] if first.strand == "+" else o.arms[0][2][-1][1]) == five_reach
    reaches_three = lambda o: (o.arms[-1][2][-1][1] if last.strand == "+" else o.arms[-1][2][0][0]) == three_reach
    five = sorted({o.read_id for o in obs if reaches_five(o)})
    three = sorted({o.read_id for o in obs if reaches_three(o)})
    # spanning is a property of ONE observation: a molecule holding a
    # 5'-reaching path and a separate 3'-reaching path has not observed the
    # envelope end to end, and the intersection of the two per-end lists
    # would claim it had
    spanned = any(reaches_five(o) and reaches_three(o) for o in obs)
    return replace(model, five_prime_observed_read_ids=tuple(five), three_prime_observed_read_ids=tuple(three),
                   envelope_status="single_molecule_spanned" if spanned else "aggregate")


def _free_extent(arm, first):
    """How far the arm reaches at its FREE (non-junction-facing) end, as a
    sort key: more exons first, then further out in transcript orientation."""
    plus = arm.strand == "+"
    if first:
        reach = arm.exons[0][0] if plus else arm.exons[-1][1]
        return (len(arm.exons), -reach if plus else reach)
    reach = arm.exons[-1][1] if plus else arm.exons[0][0]
    return (len(arm.exons), reach if plus else -reach)


def _extend_outer(lead_arm, candidates, first, tolerance):
    """The furthest-reaching compatible extent at the arm's FREE end, carrying
    the record's junction-facing boundary.

    Taking one member's arm wholesale would silently discard an extension
    another member observed (a member can be 5'-longer and 3'-shorter than
    the record at once). Only the free end moves; the junction-facing
    coordinate stays the record's, and the other members' complete placements
    remain in ``coordinate_alternatives``.
    """
    best = lead_arm
    for candidate in candidates:
        if _free_extent(candidate, first) <= _free_extent(best, first):
            continue
        compatible = (_suffix_fit(best.exons, candidate.exons, candidate.strand, tolerance) if first
                      else _prefix_fit(best.exons, candidate.exons, candidate.strand, tolerance))
        if compatible:
            best = candidate
    if best.exons == lead_arm.exons:
        return lead_arm
    exons = list(best.exons)
    plus = best.strand == "+"
    # the junction-facing boundary is the RECORD's: exit for a first arm,
    # entry for a last arm; on minus these sit at the opposite genomic end
    if first:
        if plus:
            exons[-1] = (exons[-1][0], lead_arm.exons[-1][1])
        else:
            exons[0] = (lead_arm.exons[0][0], exons[0][1])
    else:
        if plus:
            exons[0] = (lead_arm.exons[0][0], exons[0][1])
        else:
            exons[-1] = (exons[-1][0], lead_arm.exons[-1][1])
    return ModelArm(best.chrom, best.strand, tuple(exons), best.resolution)


def _aggregate_arms(members, lead, tolerance, index):
    """The record's internal placement with every member's observed extent."""
    n = len(lead.arms)
    arms = list(lead.arms)
    arms[0] = _extend_outer(lead.arms[0], [m.arms[0] for m in members], True, tolerance)
    arms[-1] = _extend_outer(lead.arms[-1], [m.arms[-1] for m in members], False, tolerance)
    if arms[0].exons == lead.arms[0].exons and arms[-1].exons == lead.arms[-1].exons:
        return tuple(lead.arms)
    if index is None:
        return tuple(arms)
    gene_ids = tuple(arm.resolution.gene_ids for arm in lead.arms)
    return _resolved(tuple(arms), index, gene_ids)


def _merge_unresolved_alternatives(models, read_models, final_id, tolerance, index=None, genome=None):
    """Fold placements the evidence itself could not separate.

    A molecule whose corrector left two placements of an internal intron a
    few bases apart carries BOTH as selected paths; each lands in its own
    envelope. Two eligible models are one population's unresolved placement
    when (a) at least one molecule is an EXACT member of both -- it holds
    both placements itself; proximity alone is never evidence -- and (b) the
    whole group stays within `tolerance` at every corresponding coordinate,
    pairwise across ALL members (a bounded diameter, not single linkage, so
    a bridge molecule cannot chain two distinguishable populations). The
    record is the placement annotation agrees with on the most arms, then
    the longest; the others are listed as ``coordinate_alternatives``. The
    same-sequence status survives only if every member spells the same
    sequence; otherwise the group is ``unresolved_nearby_placements``.
    """
    remaining = list(models)
    merged = []
    while remaining:
        head = remaining.pop(0)
        members = [head]
        changed = True
        while changed:
            changed = False
            rest = []
            for other in remaining:
                shared = any(head_m.model_id in ms and other.model_id in ms
                             for ms in read_models.values() for head_m in members)
                within = all(_nearby(m, other, tolerance) for m in members)
                if shared and within:
                    members.append(other)
                    changed = True
                else:
                    rest.append(other)
            remaining = rest
        if len(members) == 1:
            merged.append(head)
            continue
        members.sort(key=lambda m: (-_annotation_agreement(m),
                                    -sum(e - s for arm in m.arms for s, e in arm.exons),
                                    m.breakpoints, m.model_id))
        lead = members[0]
        member_ids = {m.model_id for m in members}
        arms = _aggregate_arms(members, lead, tolerance, index)
        model_id = _model_id(arms)
        for m in members:
            final_id[m.model_id] = model_id
        unique = sorted(r for r, ms in read_models.items() if ms and ms <= member_ids)
        compatible = sorted(r for r, ms in read_models.items() if ms & member_ids and not ms <= member_ids)
        record_sequence = model_sequence(arms, genome)[0] if genome is not None else ""
        same_sequence = (genome is not None and record_sequence != ""
                         and {model_sequence(m.arms, genome)[0] for m in members} == {record_sequence})
        five_status, three_status = _end_statuses(arms)
        # every member's placement, and every alternative any member already
        # carried, minus the one that became the record
        alternatives = ({placement_json(m.arms) for m in members}
                        | {q for m in members for q in m.coordinate_alternatives}) - {placement_json(arms)}
        merged.append(FusedTranscriptModel(
            model_id, arms, _breakpoints_of(arms), lead.primary_chain_id,
            tuple(sorted({c for m in members for c in m.source_chain_ids})),
            tuple(unique), tuple(compatible),
            tuple(sorted({r for m in members for r in m.partial_read_ids})), lead.chain_edge_support,
            "unresolved_microhomology_same_sequence" if same_sequence else "unresolved_nearby_placements",
            lead.event_class, lead.alternative_group_id, five_status, three_status,
            lead.five_prime_observed_read_ids, lead.three_prime_observed_read_ids,   # recomputed after the pass
            tuple(sorted(alternatives)),
            tuple(sorted({r for m in members for r in m.nearby_read_ids} - set(unique) - set(compatible))),
            tuple(sorted({q for m in members for q in m.nearby_observed_placements})),
            lead.envelope_status,
        ))
    return merged, final_id


def _junction_free_structure(arms):
    """Every exon coordinate EXCEPT the fusion-junction-facing boundaries.

    Two models are alternative placements of one fusion junction only if they
    agree on everything else: internal introns, outer ends, and the far side
    of every junction-facing exon. Anything else that happens to spell the
    same sequence is a different structure and must not be merged."""
    out = []
    n = len(arms)
    for k, arm in enumerate(arms):
        exons = list(arm.exons)
        plus = arm.strand == "+"
        if k < n - 1:   # this arm's EXIT boundary faces a junction
            if plus:
                exons[-1] = (exons[-1][0], None)
            else:
                exons[0] = (None, exons[0][1])
        if k > 0:       # this arm's ENTRY boundary faces a junction
            if plus:
                exons[0] = (None, exons[0][1])
            else:
                exons[-1] = (exons[-1][0], None)
        out.append((arm.chrom, arm.strand, tuple(exons)))
    return tuple(out)


def _same_family(a: FusedTranscriptModel, b: FusedTranscriptModel, tolerance=20):
    """Alternative placements of ONE fusion junction: identical structure away
    from the junction-facing boundaries, and every breakpoint coordinate within
    `tolerance`. Unrelated events, and same-gene structures with different
    internal exons, never merge even when they spell the same sequence."""
    if len(a.arms) != len(b.arms):
        return False
    for x, y in zip(a.arms, b.arms):
        if x.resolution.gene_ids != y.resolution.gene_ids:
            return False
    if _junction_free_structure(a.arms) != _junction_free_structure(b.arms):
        return False
    for p, q in zip(a.breakpoints, b.breakpoints):
        if p[0] != q[0] or p[3] != q[3] or abs(p[1] - q[1]) > tolerance or abs(p[4] - q[4]) > tolerance:
            return False
    return True


def _merge_same_sequence_placements(models, read_models, genome, final_id):
    """Collapse same-family models whose reference-projected sequence is
    identical (microhomology placements). The record keeps the placement
    annotation agrees with on the most arms; the others are listed. Support is
    re-derived: unique when every reported model the molecule fits is in the
    merged group, compatible when it fits the group and something outside."""
    sequences = {m.model_id: model_sequence(m.arms, genome)[0] for m in models}
    remaining = list(models)
    merged = []
    while remaining:
        head = remaining.pop(0)
        members = [head]
        rest = []
        for other in remaining:
            if sequences[other.model_id] and sequences[other.model_id] == sequences[head.model_id] and _same_family(head, other):
                members.append(other)
            else:
                rest.append(other)
        remaining = rest
        if len(members) == 1:
            merged.append(head)
            continue
        members.sort(key=lambda m: (-sum(1 for arm in m.arms if arm.resolution.compatible_transcript_ids), m.breakpoints))
        member_keys = {m.model_id for m in members}
        for m in members[1:]:
            final_id[m.model_id] = members[0].model_id
        unique = sorted(r for r, ms in read_models.items() if ms and ms <= member_keys)
        compatible = sorted(r for r, ms in read_models.items() if ms & member_keys and not ms <= member_keys)
        lead = members[0]
        merged.append(FusedTranscriptModel(
            lead.model_id, lead.arms, lead.breakpoints, lead.primary_chain_id,
            tuple(sorted({c for m in members for c in m.source_chain_ids})),
            tuple(unique), tuple(compatible),
            tuple(sorted({r for m in members for r in m.partial_read_ids})), lead.chain_edge_support,
            "unresolved_microhomology_same_sequence", lead.event_class, lead.alternative_group_id,
            lead.five_prime_status, lead.three_prime_status,
            tuple(sorted({r for m in members for r in m.five_prime_observed_read_ids})),
            tuple(sorted({r for m in members for r in m.three_prime_observed_read_ids})),
            tuple(sorted(({placement_json(m.arms) for m in members}
                          | {q for m in members for q in m.coordinate_alternatives})
                         - {placement_json(lead.arms)})),
            tuple(sorted({r for m in members for r in m.nearby_read_ids} - set(unique) - set(compatible))),
            tuple(sorted({q for m in members for q in m.nearby_observed_placements})),
            lead.envelope_status,
        ))
    return merged, final_id


# --------------------------------------------------------------------------
# sequences
# --------------------------------------------------------------------------


def _slice(genome, chrom, start, end):
    if chrom not in genome or start < 0 or end <= start:
        return ""
    sequence = genome[chrom]
    if end > len(sequence):
        return ""
    return sequence[start:end].upper()


def model_sequence(arms: Iterable[ModelArm], genome: Mapping[str, str], *, exons_override=None):
    """Reference-projected sequence of ordered arms; returns (sequence, exon rows).

    Each row: (arm_index, exon_index_in_arm, chrom, strand, start0, end0, tx_start0, tx_end0).
    ``arm.exons`` are ASCENDING genomic; a minus-strand arm is emitted in
    transcript order (descending) with each exon reverse-complemented. Read
    substitutions/insertions are NOT reconstructed.
    """
    pieces = []
    rows = []
    offset = 0
    for arm_index, arm in enumerate(arms, start=1):
        exons = exons_override[arm_index - 1] if exons_override else arm.exons
        for exon_index, (start, end) in enumerate(_oriented(exons, arm.strand), start=1):
            base = _slice(genome, arm.chrom, start, end)
            if not base:
                return "", ()
            seq = base if arm.strand == "+" else reverse_complement(base)
            rows.append((arm_index, exon_index, arm.chrom, arm.strand, start, end, offset, offset + len(seq)))
            offset += len(seq)
            pieces.append(seq)
    return "".join(pieces), tuple(rows)


def completion_exons(model: FusedTranscriptModel):
    """Annotated completion of the outer arms, or None when unavailable.

    Returns per-arm exon tuples in genomic order: the first arm extended to
    the compatible transcript's TSS, the last to its TES, middle arms as
    observed. Only defined when each outer arm's compatible transcripts agree
    on the retained structure (``completion_exons`` non-empty).
    """
    first, last = model.arms[0].resolution, model.arms[-1].resolution
    if not first.completion_exons or not last.completion_exons:
        return None
    exons = [arm.exons for arm in model.arms]
    exons[0] = first.completion_exons
    exons[-1] = last.completion_exons
    return tuple(exons)
