"""Sample-level fusion adjudication, without ordinary assembly candidates or EM."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional, Tuple

from fin.fusion.annotation import AnnotatedReadPath, GeneIndex, annotate_read_path
from fin.fusion.chains import FusionChain, build_fusion_chains, coordinate_key, topology_key
from fin.fusion.correction import PathCorrection, annotation_junction_registry, correct_read_path, correction_alternatives, observed_junction_registry
from fin.fusion.paths import find_read_paths
from fin.fusion.validation import TranscriptValidation, validate_stitched_path


@dataclass(frozen=True)
class FusionSettings:
    min_support: int = 2
    breakpoint_tolerance: int = 6
    correction_tolerance: int = 3
    min_mapping_quality: int = 20
    min_arm_bases: int = 25
    junction_anchor: int = 12
    junction_flank: int = 35
    min_identity: float = 0.7
    min_read_coverage: float = 0.9
    min_score_margin: float = 6.0
    max_query_gap: int = 29
    max_query_overlap: int = 10
    max_candidate_reads: int = 100000
    max_segments: int = 128
    max_paths: int = 64
    max_hits: int = 32
    max_windows: int = 64
    max_junction_candidates: int = 128
    min_clip_bases: int = 100

    def validate(self):
        import math

        positive = (self.min_support, self.min_arm_bases, self.junction_anchor,
                    self.max_candidate_reads, self.max_segments, self.max_paths,
                    self.max_hits, self.max_windows, self.max_junction_candidates, self.min_clip_bases)
        nonnegative = (self.breakpoint_tolerance, self.correction_tolerance,
                       self.max_query_gap, self.max_query_overlap, self.min_score_margin)
        if (min(positive) <= 0 or min(nonnegative) < 0
                or not all(math.isfinite(value) for value in (*positive, *nonnegative,
                                                             self.min_identity, self.min_read_coverage))
                or not 0 < self.min_identity <= 1 or not 0 < self.min_read_coverage <= 1
                or not 0 <= self.min_mapping_quality <= 254
                or self.junction_flank < self.junction_anchor):
            raise ValueError("invalid independent fusion settings")


@dataclass(frozen=True)
class EvaluatedPath:
    observed: AnnotatedReadPath
    proposed: AnnotatedReadPath
    correction: Optional[PathCorrection]
    validation: TranscriptValidation
    min_mapping_quality: int


@dataclass(frozen=True)
class FusionReadAudit:
    read_id: str
    status: str
    reason: str
    path_count: int = 0
    evaluated: Tuple[EvaluatedPath, ...] = ()
    selected: Tuple[AnnotatedReadPath, ...] = ()
    search_windows: Tuple[Tuple[int, int], ...] = ()
    corrections: Tuple[PathCorrection, ...] = ()
    signal: Tuple[dict, ...] = ()
    assembly_junctions: Tuple[dict, ...] = ()


@dataclass(frozen=True)
class FusionAnalysis:
    chains: Tuple[FusionChain, ...]
    reads: Tuple[FusionReadAudit, ...]
    record_counts: dict
    settings: FusionSettings
    annotation_status: str
    signal_status: str = "not_scored"


def _geometry(path):
    return (topology_key(path), tuple(tuple((block.reference_start, block.reference_end,
                                           block.query_start, block.query_end) for block in arm.blocks)
                                     for arm in path.arms))


def analyze_fusion_input(evidence, genome, aligner, genes=None, settings=None, *, annotation_transcripts=(),
                         assembly_junctions=None):
    """Adjudicate complete molecules before counting any sample-level support."""
    if not evidence.complete:
        raise ValueError("cannot analyze incomplete fusion evidence")
    settings = settings or FusionSettings()
    settings.validate()
    genes = genes if genes is not None else GeneIndex()
    audits = {}
    searches = {}
    molecules = {molecule.read_id: molecule for molecule in evidence.molecules}
    if len(molecules) != len(evidence.molecules):
        raise ValueError("duplicate molecule IDs in fusion input")
    for molecule in evidence.molecules:
        if molecule.status != "ready":
            audits[molecule.read_id] = FusionReadAudit(molecule.read_id, "unresolved", molecule.reason)
            continue
        mapped = aligner.map_molecule(molecule.read_id, molecule.sequence, molecule.segments)
        if mapped.status != "complete":
            audits[molecule.read_id] = FusionReadAudit(molecule.read_id, "unresolved", mapped.reason,
                                                      search_windows=mapped.searched_windows)
            continue
        paths = find_read_paths(
            molecule.segments + mapped.segments, min_mapping_quality=0,
            min_anchor_bases=settings.min_arm_bases, max_query_gap=settings.max_query_gap,
            max_query_overlap=settings.max_query_overlap, min_read_coverage=settings.min_read_coverage,
            max_segments=settings.max_segments, max_paths=settings.max_paths,
        )
        if not paths.paths:
            audits[molecule.read_id] = FusionReadAudit(molecule.read_id, "unresolved", paths.reason,
                                                      search_windows=mapped.searched_windows)
            continue
        annotated = {}
        for path in paths.paths:
            candidate = annotate_read_path(path, genes)
            quality = min(segment.effective_mapping_quality for segment in path.segments)
            key = _geometry(candidate)
            if key not in annotated or quality < annotated[key][1]:
                annotated[key] = candidate, quality
        searches[molecule.read_id] = tuple(annotated[key] for key in sorted(annotated)), mapped.searched_windows
    registry = observed_junction_registry(candidate for candidates, _ in searches.values()
                                           for candidate, _ in candidates)
    auxiliary = {}
    if assembly_junctions is not None:
        auxiliary = assembly_junctions.registry(
            (candidate for candidates, _ in searches.values() for candidate, _ in candidates),
            tolerance=settings.correction_tolerance,
        )
    annotation_registry = annotation_junction_registry(
        (candidate for candidates, _ in searches.values() for candidate, _ in candidates),
        annotation_transcripts, tolerance=settings.correction_tolerance,
    )
    representatives = []
    for read_id in sorted(searches):
        candidates, windows = searches[read_id]
        molecule = molecules[read_id]
        evaluated = []
        failures = []
        corrections = []
        for candidate, quality in candidates:
            correction = None
            if len(candidate.arms) < 2:
                proposals = (candidate,)
            else:
                correction = correct_read_path(
                    candidate, molecule.sequence, genome, supported_registry=registry, assembly_registry=auxiliary,
                    annotation_registry=annotation_registry,
                    tolerance=settings.correction_tolerance, flank=settings.junction_flank,
                    min_anchor=settings.junction_anchor, min_identity=settings.min_identity,
                    min_margin=settings.min_score_margin, max_candidates=settings.max_junction_candidates,
                )
                corrections.append(correction)
                proposals, reason = correction_alternatives(
                    correction, min_margin=settings.min_score_margin, max_paths=settings.max_paths,
                )
                if not proposals:
                    failures.append(reason)
            for proposal in proposals:
                identity_changed = genes.genes and any(
                    genes.labels(arm.chrom, arm.strand, (block,), min_overlap_bases=1) != arm.gene_ids
                    for arm in proposal.arms for block in arm.blocks
                )
                if identity_changed:
                    validation = TranscriptValidation("unresolved", "corrected_gene_identity_changed")
                else:
                    validation = validate_stitched_path(
                        proposal, molecule.sequence, genome, min_anchor=settings.junction_anchor,
                        min_identity=settings.min_identity, min_coverage=settings.min_read_coverage,
                        max_hits=settings.max_hits,
                    )
                evaluated.append(EvaluatedPath(candidate, proposal, correction, validation, quality))
                if validation.status == "unresolved":
                    failures.append(validation.reason)
        evaluated = tuple(evaluated)
        technical = sorted(set(failures) - {"no_anchored_hypothesis", "baseline_not_anchored",
                                           "joint_template_not_anchored", "incompatible_joint_edits"})
        valid = [item for item in evaluated if item.validation.status == "validated"]
        if technical or not valid:
            audits[read_id] = FusionReadAudit(read_id, "unresolved", technical[0] if technical else
                                             "no_validated_complete_transcript", len(candidates), evaluated,
                                             search_windows=windows, corrections=tuple(corrections))
            continue
        best = max(item.validation.score for item in valid)
        competitive = [item for item in valid if item.validation.score == best
                       or best - item.validation.score < settings.min_score_margin]
        selected = tuple({_geometry(item.proposed): item.proposed for item in competitive}.values())
        selected = tuple(sorted(selected, key=_geometry))
        topologies = {topology_key(path) for path in selected}
        coordinates = [coordinate_key(path) for path in selected]
        if any(len(path.arms) < 2 for path in selected):
            status = "unresolved" if any(len(path.arms) >= 2 for path in selected) else "not_fusion"
            reason = "normal_alignment_competitor"
        elif len(topologies) > 1:
            status, reason = "unresolved", "competing_gene_or_splice_paths"
        elif any(path.identity_status == "ambiguous" for path in selected):
            status, reason = "unresolved", "ambiguous_gene_identity"
        elif any(max(values) - min(values) > settings.breakpoint_tolerance for values in zip(*coordinates)):
            status, reason = "unresolved", "competing_distant_coordinates"
        elif any(item.min_mapping_quality < settings.min_mapping_quality for item in competitive):
            status, reason = "unresolved", "low_mapping_quality"
        else:
            status = "supported_path"
            reason = "coordinate_ambiguity" if len(set(coordinates)) > 1 else "whole_chain_validated"
            representatives.append(selected[0])
        audits[read_id] = FusionReadAudit(read_id, status, reason, len(candidates), evaluated, selected, windows,
                                         tuple(corrections))
    if assembly_junctions is not None:
        audits = {read_id: replace(read, assembly_junctions=tuple(
            {"alternative_index": index, "path_coordinates": coordinate_key(path),
             "junctions": assembly_junctions.annotate(path)}
            for index, path in enumerate(read.selected))) for read_id, read in audits.items()}
    return FusionAnalysis(finalize_chains(representatives, audits, settings),
                          tuple(audits[read_id] for read_id in sorted(audits)),
                          evidence.record_counts, settings, "provided" if genes.genes else "unavailable_loci_only")


def finalize_chains(representatives, audits, settings):
    chains = build_fusion_chains(representatives, min_full_support=settings.min_support,
                                 breakpoint_tolerance=settings.breakpoint_tolerance)
    output = []
    for chain in chains:
        ambiguous = any(audits[read_id].reason == "coordinate_ambiguity" for read_id in chain.full_support_read_ids)
        if ambiguous:
            chain = replace(chain, coordinate_status="unresolved_sequence_alternatives")
        if chain.status == "supported_candidate":
            chain = replace(chain, status="supported_rna_fusion" if chain.coordinate_status == "exact_observed"
                            else "supported_chain_unresolved_coordinates")
        if chain.coordinate_status == "exact_observed":
            signal_resolved = any(audits[read_id].reason == "signal_resolved_coordinates"
                                  for read_id in chain.full_support_read_ids)
            chain = replace(chain, coordinate_status="signal_resolved" if signal_resolved else "sequence_resolved")
        output.append(chain)
    return tuple(output)
