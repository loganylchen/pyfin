"""File-level independent fusion analysis with its own output lifecycle."""

from __future__ import annotations

from dataclasses import asdict
from contextlib import nullcontext
import hashlib
import logging
import time
from importlib.metadata import version
from pathlib import Path

from fin.fusion.analysis import FusionSettings, analyze_fusion_input
from fin.fusion.annotation import GeneIndex
from fin.fusion.input import (FusionInputError, SNAPSHOT_FIELDS, _snapshot, scan_fusion_bam,
                              supply_fusion_sequences)
from fin.fusion.output import FusionPublication, write_fusion_outputs
from fin.fusion.realignment import FusionGenomeAligner

logger = logging.getLogger(__name__)


def run_fusion_analysis(bam_path, reference_path, output_dir, *, gtf_path=None,
                        sequence_path=None, settings=None, threads=1, publication=None,
                        signal_path=None, signal_settings=None, assembly_gtf_path=None,
                        assembly_support_path=None, assembly_min_clean_reads=3,
                        realign_workers=1, genome_cache_bases=4_000_000_000):
    """Never instantiate an assembly runner or touch its candidates or caches."""
    from fin.io.io_gtf import GTFReader
    from fin.io.lazy_genome import open_genome

    settings = settings or FusionSettings()
    context = FusionPublication(output_dir) if publication is None else nullcontext(publication)
    with context as publication:
        settings.validate()
        if (assembly_gtf_path is None) != (assembly_support_path is None):
            raise ValueError("assembly junction evidence requires both final GTF and support TSV")
        inputs = {name: str(Path(path).resolve()) for name, path in
                  (("bam", bam_path), ("reference", reference_path), ("gtf", gtf_path), ("sequence", sequence_path),
                   ("signal", signal_path), ("assembly_gtf", assembly_gtf_path),
                   ("assembly_junction_support", assembly_support_path))
                  if path is not None}
        snapshots = {name: _snapshot(path) for name, path in inputs.items()}
        assembly_junctions = None
        if assembly_gtf_path is not None:
            from fin.fusion.assembly_evidence import load_assembly_junctions

            assembly_junctions = load_assembly_junctions(
                assembly_gtf_path, assembly_support_path, min_clean_reads=assembly_min_clean_reads,
            )
        # a real library is ~1e6 records and a genome index is minutes: a run
        # with no progress is indistinguishable from a hung one. The clock
        # starts HERE, before the annotation is read, because it used to start
        # after -- which reported "annotation indexed [0s]" for a GENCODE parse
        # that is tens of seconds of real work, and left every later phase
        # quoting an elapsed time the wall clock disagreed with.
        started = time.monotonic()
        elapsed = lambda: f"{time.monotonic() - started:.0f}s"
        genes = GeneIndex()
        transcripts = ()
        if gtf_path is not None:
            with GTFReader(str(gtf_path)) as reader:
                reader.parse()
                transcripts = tuple(reader.transcripts.values())
                genes = GeneIndex.from_transcripts(transcripts)
            logger.info("fusion: annotation indexed, %d genes [%s]", len(genes.genes), elapsed())
        logger.info("fusion: scanning %s (min clip %d bp) [%s]", bam_path, settings.min_clip_bases, elapsed())
        evidence = scan_fusion_bam(bam_path, genes, min_soft_clip=settings.min_clip_bases,
                                   max_candidate_reads=settings.max_candidate_reads,
                                   max_segments_per_read=settings.max_segments)
        logger.info("fusion: %s [%s]", " ".join(f"{k}={v}" for k, v in sorted(evidence.record_counts.items())),
                    elapsed())
        if sequence_path is not None:
            logger.info("fusion: recovering read sequences from %s [%s]", sequence_path, elapsed())
            evidence = supply_fusion_sequences(evidence, sequence_path)
        logger.info("fusion: building the genome index [%s]", elapsed())
        aligner = FusionGenomeAligner(
            reference_path, max_hits=settings.max_hits, max_windows=settings.max_windows,
            max_segments=settings.max_segments, min_anchor=settings.min_arm_bases,
            min_identity=settings.min_identity, threads=threads,
        )
        # The fusion adjudication walks candidate molecules in read-id order,
        # so the chromosomes it needs arrive in effectively random order --
        # nothing like the assembly pipeline's one-locus-at-a-time sweep that
        # the 2-chromosome default was chosen for. Measured on a whole-genome
        # run under that default: 499 misses, 497 of them evicting the entry
        # the next miss wanted, 66.6 Gb fetched from a 3.1 Gb genome, 246 s.
        # Bound by SIZE instead, so the same number is right for a spike-in
        # reference and for a genome larger than the budget.
        genome = open_genome(str(reference_path), cache_chroms=None,
                             cache_bases=genome_cache_bases)
        logger.info("fusion: index ready; adjudicating %d molecules on %d worker(s) [%s]",
                    sum(1 for m in evidence.molecules if m.status == "ready"), realign_workers, elapsed())
        try:
            analysis = analyze_fusion_input(evidence, genome, aligner, genes, settings,
                                            annotation_transcripts=transcripts, assembly_junctions=assembly_junctions,
                                            realign_workers=realign_workers)
            if signal_path is not None:
                from fin.fusion.signal import FusionSignalSettings, validate_fusion_signal

                signal_settings = signal_settings or FusionSignalSettings()
                analysis = validate_fusion_signal(analysis, evidence, genome, signal_path, signal_settings)
            logger.info("fusion: %d chains from %d molecules [%s]",
                        len(analysis.chains), len(analysis.reads), elapsed())
            write_fusion_outputs(analysis, genome, publication.work_dir)
            # Fused-transcript level, computed from the finished molecule-level
            # analysis and written BESIDE the chain artifacts: which isoform of
            # each gene, the conventional exon-rank junction name, labelled
            # 5'/3' completion, and molecule assignment with unique vs
            # compatible votes. Nothing here feeds back into chains/signal.
            from fin.fusion.transcript_models import TranscriptIndex, build_fused_transcript_models
            from fin.fusion.transcript_output import write_transcript_outputs

            transcript_index = TranscriptIndex.from_transcripts(transcripts)
            models, collapsed = build_fused_transcript_models(analysis, transcript_index, genome)
            write_transcript_outputs(models, collapsed, genome, publication.work_dir)
            logger.info("fusion: %d fused-transcript models (%d multi-gene), %d collapsed [%s]",
                        len(models), sum(1 for m in models if len(m.arms) > 2), len(collapsed), elapsed())
            drifted = {name: (snapshot, _snapshot(inputs[name])) for name, snapshot in snapshots.items()
                       if _snapshot(inputs[name]) != snapshot}
            if drifted:
                fields = SNAPSHOT_FIELDS
                detail = "; ".join(
                    f"{name} [{inputs[name]}] " + ", ".join(
                        f"{field} {was}->{now}" for field, was, now in zip(fields, before, after) if was != now)
                    for name, (before, after) in sorted(drifted.items()))
                raise FusionInputError("fusion input changed before publication: " + detail)
            source_root = Path(__file__).resolve().parents[2]
            source_digest = hashlib.sha256()
            for source in sorted((source_root / "fin").rglob("*.py")):
                source_digest.update(str(source.relative_to(source_root)).encode())
                source_digest.update(b"\0" + source.read_bytes() + b"\0")
            publication.complete({
                "settings": asdict(settings), "inputs": inputs, "input_snapshots": snapshots,
                "realign_workers": realign_workers,
                "record_counts": analysis.record_counts, "annotation_status": analysis.annotation_status,
                "fused_transcript_models": {
                    "count": len(models),
                    "unique_support_total": sum(len(m.unique_read_ids) for m in models),
                    "multi_gene_models": sum(1 for m in models if len(m.arms) > 2),
                    "collapsed_structures": len(collapsed),
                    "interpretation": ("evidence-constrained models in fusion_models.fa; annotated completion "
                                       "hypotheses in fusion_models.annotated.fa are NOT observations"),
                },
                "signal_status": analysis.signal_status, "mappy_version": version("mappy"),
                "signal_settings": asdict(signal_settings) if signal_path is not None else None,
                "signal_method": "krill_local_summed_llr" if signal_path is not None else None,
                "assembly_junction_evidence": {
                    "status": "provided" if assembly_junctions is not None else "not_provided",
                    "junction_count": len(assembly_junctions.junctions) if assembly_junctions is not None else 0,
                    "min_clean_reads": assembly_min_clean_reads,
                    "interpretation": "final_model_and_clean_reads; nominations_only; not_independent_signal_support",
                },
                "source_root": str(source_root), "source_sha256": source_digest.hexdigest(),
                "assembly_isolation": "independent_analysis_no_assembly_candidates",
                "interpretation": "RNA fusion evidence, not proof of a genomic rearrangement",
                "bedpe_coordinates": "0-based half-open junction-facing exonic bases; representative coordinates",
            })
        finally:
            close = getattr(genome, "close", None)
            if close is not None:
                close()
    return analysis
