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
from dataclasses import dataclass, field
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
    # other placements of the SAME sequence within the same fusion family
    # (microhomology): the sequence record is one, the coordinates are several
    coordinate_alternatives: Tuple[Tuple[Tuple[str, int, str, str, int, str], ...], ...] = ()

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


def _fits_within(shorter: _Observation, envelope_arms):
    """True when `shorter` observed the same fusion and splicing as the
    envelope but fewer 5'-most exons of the FIRST arm (direct-RNA 5'
    degradation), with its 5'-most exon lying INSIDE the envelope's
    corresponding exon and ending at the same 3'-facing boundary. Outer 3'
    ends are free; every other arm must carry identical splicing/breakpoints.
    Judged on THIS molecule's own coordinates, never on a representative."""
    if len(shorter.arms) != len(envelope_arms):
        return False
    env = _Observation("", tuple((a.chrom, a.strand, a.exons, ()) for a in envelope_arms), "")
    if shorter.splice_key[1:] != env.splice_key[1:]:
        return False
    (chrom, strand, exons, _), first = shorter.arms[0], envelope_arms[0]
    if (chrom, strand) != (first.chrom, first.strand):
        return False
    sa, sb = _oriented(exons, strand), _oriented(first.exons, strand)
    if len(sa) > len(sb):
        return False
    tail = sb[len(sb) - len(sa):]
    if sa[1:] != tail[1:]:
        return False
    head, ref = sa[0], tail[0]
    if strand == "+":
        return head[1] == ref[1] and ref[0] <= head[0]
    return head[0] == ref[0] and head[1] <= ref[1]


def _breakpoints_of(arms):
    out = []
    for l, r in zip(arms, arms[1:]):
        exit_ = l.exons[-1][1] if l.strand == "+" else l.exons[0][0]
        entry = r.exons[0][0] if r.strand == "+" else r.exons[-1][1]
        out.append((l.chrom, exit_, l.strand, r.chrom, entry, r.strand))
    return tuple(out)


def _model_id(arms):
    return "FTM_" + hashlib.sha256(_json(structure_key_of(arms)).encode()).hexdigest()[:16]


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


def build_fused_transcript_models(analysis, transcript_index: TranscriptIndex,
                                  genome: Optional[Mapping[str, str]] = None):
    """Turn supported chains into explicit fused-transcript models.

    Returns ``(models, collapsed)``.

    Observed structures (splice keys) are collected from EVERY chain's full-
    support molecules: a molecule in an insufficient-support chain may be a
    5'-degraded view of a supported model and then belongs to it as compatible
    evidence. Only structures seen in a SUPPORTED chain may become models.

    Maximality is decided among eligible structures and PER MOLECULE: a
    structure is collapsed (not reported) only when every one of its
    observations fits inside some eligible longer structure's envelope. A
    single observation that cannot fit -- it saw sequence the longer structure
    splices out -- keeps its structure as a model. Collapsed structures are
    returned, never silently dropped.

    Votes: a molecule is UNIQUE to a reported model when it is compatible with
    that model alone among reported models ("unique among reported models",
    not "observed the whole model"); otherwise it is COMPATIBLE on each model
    it fits. One molecule, at most one unique vote.

    Same-sequence alternative placements (microhomology) within one fusion
    family are merged into one record; unrelated structures that happen to
    spell the same sequence are not.
    """
    reads = {read.read_id: read for read in analysis.reads}
    chain_by_id = {chain.chain_id: chain for chain in analysis.chains}
    structures: Dict[tuple, dict] = {}
    for chain in analysis.chains:
        is_supported = chain.status.startswith("supported")
        for read_id in chain.full_support_read_ids:
            audit = reads.get(read_id)
            if audit is None:
                continue
            for path in audit.selected:
                if len(path.arms) < 2:
                    continue
                observation = _observation(read_id, path)
                entry = structures.setdefault(observation.splice_key, {"obs": [], "chains": set(), "supported": False})
                entry["obs"].append(observation)
                entry["chains"].add(chain.chain_id)
                entry["supported"] = entry["supported"] or is_supported
    if not structures:
        return (), ()
    keys = sorted(structures, key=_json)
    for key in keys:
        structures[key]["obs"].sort(key=lambda o: (o.read_id, _json(o.arms)))
        structures[key]["arms"] = _envelope(structures[key]["obs"], transcript_index)
    eligible = [k for k in keys if structures[k]["supported"]]
    # per OBSERVATION (a molecule may hold several selected paths under one
    # splice key -- e.g. two outer-start alternatives -- and each must be
    # judged on its own blocks): eligible structures other than its own that
    # it fits inside. Keyed by (key, position), never by read id alone.
    fits_other: Dict[Tuple[tuple, int], set] = {}
    for key in keys:
        for position, obs in enumerate(structures[key]["obs"]):
            fits_other[(key, position)] = {
                other for other in eligible
                if other != key and _fits_within(obs, structures[other]["arms"])}
    collapsed_keys = {
        key for key in keys
        if structures[key]["obs"] and all(fits_other[(key, position)]
                                          for position in range(len(structures[key]["obs"])))}
    model_keys = [k for k in eligible if k not in collapsed_keys]
    # per molecule: the reported models it is compatible with (union over its
    # observations -- votes are per molecule, compatibility per observation)
    read_models: Dict[str, set] = defaultdict(set)
    for key in keys:
        for position, obs in enumerate(structures[key]["obs"]):
            if key in model_keys:
                read_models[obs.read_id].add(key)
            read_models[obs.read_id].update(m for m in fits_other[(key, position)] if m in model_keys)
    models: List[FusedTranscriptModel] = []
    for key in model_keys:
        arms = structures[key]["arms"]
        obs = structures[key]["obs"]
        unique = sorted(r for r, ms in read_models.items() if ms == {key})
        compatible = sorted(r for r, ms in read_models.items() if key in ms and len(ms) > 1)
        chains = set(structures[key]["chains"])
        for other in collapsed_keys:
            if any(key in fits_other[(other, position)] for position in range(len(structures[other]["obs"]))):
                chains.update(structures[other]["chains"])
        primary = sorted(c for c in structures[key]["chains"] if chain_by_id[c].status.startswith("supported"))[0]
        chain = chain_by_id[primary]
        first, last = arms[0], arms[-1]
        five_reach = first.exons[0][0] if first.strand == "+" else first.exons[-1][1]
        three_reach = last.exons[-1][1] if last.strand == "+" else last.exons[0][0]
        five_obs = sorted({o.read_id for o in obs if (o.arms[0][2][0][0] if first.strand == "+" else o.arms[0][2][-1][1]) == five_reach})
        three_obs = sorted({o.read_id for o in obs if (o.arms[-1][2][-1][1] if last.strand == "+" else o.arms[-1][2][0][0]) == three_reach})
        five, three = first.resolution, last.resolution
        five_status = ("observed_boundary_5p_completeness_unknown" if not five.completion_exons
                       else "annotated_boundary_observed" if five.completion_exons == first.exons
                       else "annotated_extension_available")
        three_status = ("observed_end_completeness_unknown" if not three.completion_exons
                        else "annotated_boundary_observed" if three.completion_exons == last.exons
                        else "annotated_extension_available")
        siblings = {m for r in compatible for m in read_models[r]} | {key}
        group = "CAG_" + hashlib.sha256(_json(sorted(_json(m) for m in siblings)).encode()).hexdigest()[:12]
        models.append(FusedTranscriptModel(
            _model_id(arms), arms, _breakpoints_of(arms), primary, tuple(sorted(chains)),
            tuple(unique), tuple(compatible), tuple(chain.partial_support_read_ids),
            tuple(len(ids) for ids in chain.edge_support_read_ids),
            chain.coordinate_status, obs[0].event_class, group, five_status, three_status,
            tuple(five_obs), tuple(three_obs),
        ))
    # pre-merge structure key -> FINAL emitted model id, so the collapsed
    # audit can never point at a record the merge step removed
    final_id: Dict[tuple, str] = {key: _model_id(structures[key]["arms"]) for key in model_keys}
    if genome is not None:
        models, final_id = _merge_same_sequence_placements(models, read_models, model_keys, genome, final_id)
    models.sort(key=lambda m: (-len(m.unique_read_ids), -len(m.compatible_read_ids), m.model_id))
    emitted = {m.model_id for m in models}
    assert set(final_id.values()) <= emitted, "collapsed audit would reference a removed model"
    collapsed = []
    for key in sorted(collapsed_keys, key=_json):
        targets = sorted({final_id[m] for position in range(len(structures[key]["obs"]))
                          for m in fits_other[(key, position)] if m in final_id})
        summary = ";".join(f"{c}:{s}:" + ",".join(f"{a}-{b}" for a, b in e) for c, s, e, _ in structures[key]["obs"][0].arms)
        collapsed.append(CollapsedStructure(summary, tuple(sorted({o.read_id for o in structures[key]["obs"]})), tuple(targets)))
    return tuple(models), tuple(collapsed)


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


def _merge_same_sequence_placements(models, read_models, model_keys, genome, final_id):
    """Collapse same-family models whose reference-projected sequence is
    identical (microhomology placements). The record keeps the placement
    annotation agrees with on the most arms; the others are listed. Support is
    re-derived: unique when every reported model the molecule fits is in the
    merged group, compatible when it fits the group and something outside."""
    key_of = {m.structure_key: key for m, key in zip(models, model_keys)}
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
        member_keys = {key_of[m.structure_key] for m in members}
        for m in members[1:]:
            final_id[key_of[m.structure_key]] = members[0].model_id
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
            tuple(m.breakpoints for m in members[1:]),
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
