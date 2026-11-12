"""Pipeline orchestrator: interval -> candidates -> scoring -> EM -> quantification."""

from __future__ import annotations

import logging
import math
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pysam

from fin.analysis.abundance_refit import (
    ResponsibilityLedger,
    annotate_selection_metadata,
    build_responsibility_ledger,
    refit_survivor_abundance,
    write_abundance_refit_diagnostics,
)
from fin.analysis.assignments import em_with_coherence
from fin.analysis.quantification import (
    QuantResult,
    _exons_from_candidate,
    aggregate_across_intervals,
    compute_fulllen_frac,
    polya5p_drops,
    quantify_transcripts,
)
from fin.candidates.dataclasses import CandidateSet
from fin.candidates.discovery import discover_candidates, merge_fusion_candidates
from fin.io.interval_manager import GenomicInterval, generate_isolated_intervals
from fin.pipeline.assignment import Assigner
from fin.pipeline.config import PipelineConfig
from fin.pipeline.finalize import finalize_outputs, write_unfiltered_diagnostic
from fin.pipeline.junction_snap import apply_junction_snap
from fin.pipeline.selection import (
    append_selection_audit_global,
    audit_rows_for_candidates,
    audit_rows_for_outcomes,
    canonical_gate_select,
    enforce_gtf_positive_after_refit,
    junction_dominance_select,
    select_global,
    select_m2_interval,
    unify_extent_pairs,
    write_selection_audit_shard,
)
from fin.scoring.em_inputs import build_em_matrices
from fin.scoring.krill_tiebreak import krill_tiebreak
from fin.scoring.mappy_distance import compute_mappy_distance

logger = logging.getLogger(__name__)

# Pure tie-break junction-NLL M2-EM constants (SIRV-validated experiment values).
# K (transcript-bp each side of the wobbling junction) comes from the configurable
# m2_tiebreak_junction_k; FLANK is the diff-window padding and PAD the penalty for
# an unscorable tie cell (a tied candidate the eventalign window could not score).
_M2_EM_FLANK = 2
_M2_EM_PAD = 1.0


IntervalOutput = Tuple[List[QuantResult], Optional[ResponsibilityLedger]]


def _usable_polya_call(value):
    """``[length, qc]`` when ``value`` is a usable polyA call, else None.

    A call is a pair whose length is absent or a finite number and whose qc is
    absent or a string. Anything else (a NaN length, a wrong shape) is not a
    failed call but a broken one: it is dropped at the source so it can never be
    written into the sample evidence, where it would invalidate the whole
    artifact -- read ends included -- for every later consumer."""
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        return None
    length, qc = value
    if not (qc is None or isinstance(qc, str)):
        return None
    if length is None:
        return [None, qc]
    if isinstance(length, bool) or not isinstance(length, (int, float)):
        return None
    if isinstance(length, float) and not math.isfinite(length):
        return None
    if not -2.0 ** 63 <= length <= 2.0 ** 63:
        return None
    return [length, qc]


def input_fingerprint(path) -> str:
    """Fingerprint of an INPUT file, used to decide whether derived evidence
    (read endpoints, polyA calls) may be reused for it.

    Default (fast): ``path:size:mtime:ctime:digest(head 64KiB + tail 64KiB)``.
    This DETECTS a file regenerated or rewritten at the same path -- the normal
    way an input changes -- because such a write moves the timestamps and almost
    always the size or the header/EOF bytes. It does NOT PROVE content identity:
    a same-length in-place edit of the interior that also restores the
    timestamps is indistinguishable, and on filesystems where ctime tracks mtime
    the restored-timestamp case is invisible too. Only a full hash proves
    identity, and that costs the scan this check exists to skip.

    Set ``PYFIN_EVIDENCE_STRICT_HASH=1`` to digest the whole file instead, for
    setups where inputs may be edited in place behind the pipeline. The cost is
    one sequential read of the input per run (still well under the per-record
    decode it guards)."""
    import hashlib
    import os

    if not path:
        return ""
    p = Path(path)
    strict = os.environ.get("PYFIN_EVIDENCE_STRICT_HASH", "") not in ("", "0")
    try:
        st = p.stat()
        h = hashlib.blake2b(digest_size=8)
        with p.open("rb") as fh:
            if strict:
                for chunk in iter(lambda: fh.read(8 << 20), b""):
                    h.update(chunk)
            else:
                edge = 64 * 1024
                h.update(fh.read(edge))
                if st.st_size > edge:
                    fh.seek(max(st.st_size - edge, edge))
                    h.update(fh.read(edge))
        return (f"{p.resolve()}:{st.st_size}:{int(st.st_mtime_ns)}"
                f":{int(st.st_ctime_ns)}:{'full' if strict else 'edge'}:{h.hexdigest()}")
    except OSError:
        return str(p)


def _order_interval_outputs(
    outputs: List[IntervalOutput], refit_effective: bool
) -> List[IntervalOutput]:
    """Canonicalize aggregation order only for the refit-enabled path."""
    if not refit_effective:
        return outputs
    if any(ledger is None for _, ledger in outputs):
        raise RuntimeError(
            "post-selection refit requires a responsibility ledger for every "
            "quantified interval"
        )
    # Genomic lexical order is sufficient; interval_key includes strand and is
    # unique for the isolated interval construction.
    return sorted(outputs, key=lambda pair: pair[1].interval_key)


class PipelineRunner:
    """Orchestrates the full pyfin pipeline."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self._gtf_reader = None
        self._genome_fasta = None
        self._signal_reader = None
        # Genome-wide mappy aligner for fusion soft-clip arm re-alignment; built
        # lazily on first fusion interval and cached (indexing is expensive).
        self._fusion_genome_aligner = None
        self._fusion_aligner_built = False
        # m2_em assignment layer (mappy AS -> tie_nll -> d_tx -> EM -> quantify).
        self._assigner = Assigner(config)

    def setup(self):
        """Open file handles and load references.

        All signal scoring is in-memory krill (no f5c CLI / external-tool
        validation / f5c index build required).
        """
        # Load GTF
        if self.config.gtf_path:
            self._gtf_reader = self._open_gtf_reader()

        # Load genome FASTA (full sequence per chrom)
        if self.config.genome_fasta_path:
            self._genome_fasta = self._load_genome_fasta(self.config.genome_fasta_path)

        # Open signal reader
        if self.config.signal_path:
            self._signal_reader = self._open_signal_reader()

        # Create work directory
        Path(self.config.work_dir).mkdir(parents=True, exist_ok=True)

        logger.info("Pipeline setup complete")

    def run(self) -> Dict[str, QuantResult]:
        """Run the full pipeline across all intervals.

        Returns:
            Aggregated quantification results.
        """
        # A runner may be reused after its inputs changed: re-fingerprint them
        # before anything derived from the old ones is consulted.
        self._reset_input_snapshot()
        # Generate intervals
        if getattr(self.config, "selection_audit", False):
            # a stale shard or global file left by an earlier run in the same
            # work dir would be read as this run's record, so a failure to
            # remove it is fatal, not ignorable
            import shutil as _shutil
            _shard_dir = Path(self.config.work_dir) / "selection_audit"
            _shutil.rmtree(_shard_dir, ignore_errors=True)
            (Path(self.config.work_dir) / "selection_audit_global.tsv").unlink(
                missing_ok=True)
            if _shard_dir.exists():
                raise RuntimeError(
                    f"could not clear the stale selection audit at {_shard_dir}; "
                    "refusing to append this run's records to another run's"
                )
        # Lend the annotation parsed in setup() instead of parsing it again.
        # This is AFTER _reset_input_snapshot(), so the reader lent here is the
        # one whose fingerprint describes the bytes this run is using.
        result = generate_isolated_intervals(
            self.config.bam_path,
            gtf_path=self.config.gtf_path,
            max_gap=self.config.max_gap,
            max_reads=self.config.max_reads,
            gtf_reader=self._gtf_reader,
        )
        intervals = result["intervals"]
        logger.info("Generated %d intervals", len(intervals))

        # Sample-level evidence prepass for the existence audit: one krill polyA
        # pass over the sample and the strand-aware 3' ends, cached to the work
        # dir so every worker judges intervals against the SAME versioned model.
        if (getattr(self.config, "candidate_existence_mode", "off") != "off"
                or getattr(self.config, "assignment_channels", "none") != "none"):
            self._prepare_sample_evidence()

        # Process each interval (serial by default; opt-in process parallelism).
        interval_outputs: List[IntervalOutput] = []
        if self.config.threads <= 1:
            for i, interval in enumerate(intervals):
                logger.info("Processing interval %d/%d: %s", i + 1, len(intervals), interval.region_string)
                quant = self.process_interval(interval)
                if quant is not None:
                    interval_outputs.append(quant)
        else:
            from fin.pipeline.parallel import run_parallel

            log_level = logging.getLevelName(logging.getLogger("fin").getEffectiveLevel())
            interval_outputs = run_parallel(
                self.config,
                intervals,
                self.config.threads,
                self.config.gpu_workers,
                log_level,
            )

        refit_effective = bool(getattr(
            self.config, "post_selection_refit_effective", False
        ))
        interval_outputs = _order_interval_outputs(
            interval_outputs, refit_effective
        )
        all_quant_results: List[List[QuantResult]] = []
        responsibility_ledgers: List[ResponsibilityLedger] = []
        for results, ledger in interval_outputs:
            if results:
                all_quant_results.append(results)
            if ledger is not None:
                responsibility_ledgers.append(ledger)

        # An annotated model's extent pair (annotation vs read-extent sibling)
        # that split across intervals is unified BEFORE aggregation, so its
        # reads and soft mass are deduplicated under one identity.
        extent_redirects, extent_outcomes = unify_extent_pairs(all_quant_results)
        if getattr(self.config, "selection_audit", False):
            append_selection_audit_global(
                self.config.work_dir,
                audit_rows_for_outcomes("*", "extent_unification", extent_outcomes))
        if extent_redirects:
            logger.info("Unified %d annotated extent pair(s) across intervals",
                        len(extent_redirects))
        # Aggregate across intervals. Responsibility ledgers remain separate so
        # the historical aggregation/selection evidence is not changed by refit.
        aggregated = aggregate_across_intervals(all_quant_results)
        return self._finalize_and_write(
            aggregated,
            output_gtf=self.config.output_gtf,
            output_tsv=self.config.output_tsv,
            responsibility_ledgers=responsibility_ledgers,
            extent_redirects=extent_redirects,
        )

    def _finalize_and_write(
        self,
        aggregated: Dict[str, QuantResult],
        output_gtf: Optional[str],
        output_tsv: Optional[str],
        responsibility_ledgers: Optional[List[ResponsibilityLedger]] = None,
        extent_redirects: Optional[Dict[str, str]] = None,
    ) -> Dict[str, QuantResult]:
        """Finalize: diagnostic dump -> GLOBAL selection -> gene-id resolution + writers.

        Thin orchestrator wiring three independent layers: the pre-filter diagnostic and the
        output writers live in fin.pipeline.finalize; the GLOBAL selection cascade lives in
        fin.pipeline.selection. Kept as a method so tests can call this seam directly."""
        # P2 containment finalization: merge worker shards, apply the
        # grouped-FDR pass to hazard rows, execute the CORRECTED require-mode
        # hazard actions (drop demoted/contradicted shorts, mark supported
        # ones for the soft-mass exemption), and write the deterministic
        # containment_evidence.tsv. Runs BEFORE the global cascade so the
        # cascade sees final verdicts. audit mode only writes the artifact.
        if getattr(self.config, "containment_evidence_mode", "off") != "off":
            from fin.analysis.containment_evidence import (
                finalize_containment,
            )

            aggregated, c_dropped, c_path = finalize_containment(
                self.config.work_dir, aggregated,
                getattr(self.config, "containment_evidence_mode", "off"),
            )
            if c_dropped:
                logger.info(
                    "Containment evidence (FDR-corrected) dropped %d "
                    "contained model(s): %s",
                    len(c_dropped), ",".join(c_dropped[:8]),
                )
            logger.info("Containment evidence written: %s", c_path)
        if getattr(self.config, "junction_election_mode", "off") != "off":
            from fin.analysis.junction_election import merge_election_shards

            je_path = merge_election_shards(self.config.work_dir)
            logger.info("Junction election written: %s", je_path)
        from fin.analysis.junction_support import merge_support_shards

        js_path = merge_support_shards(self.config.work_dir)
        if js_path:
            logger.info("Junction support written: %s", js_path)
        if getattr(self.config, "candidate_existence_mode", "off") != "off":
            from fin.analysis.candidate_existence import merge_existence_shards

            ce_path = merge_existence_shards(self.config.work_dir)
            if ce_path:
                logger.info("Candidate existence audit written: %s", ce_path)
        if getattr(self.config, "folded_tss_recovery_mode", "off") != "off":
            from fin.analysis.tss_shadow_recovery import merge_recovery_shards

            tr_path = merge_recovery_shards(self.config.work_dir)
            if tr_path:
                logger.info("Folded-TSS shadow recovery written: %s", tr_path)
        if getattr(self.config, "junction_correction_mode", "off") != "off":
            from fin.analysis.junction_correction import merge_correction_shards

            jc_path = merge_correction_shards(self.config.work_dir)
            logger.info("Junction correction written: %s", jc_path)
        write_unfiltered_diagnostic(
            self.config, aggregated, output_tsv, self._gtf_reader
        )
        aggregated, _global_outcomes = select_global(
            self.config, aggregated, self._apply_polya5p_filter,
        )
        if getattr(self.config, "selection_audit", False):
            append_selection_audit_global(
                self.config.work_dir,
                audit_rows_for_outcomes("*", "global_cascade", _global_outcomes))
        want_evidence = bool(getattr(self.config, "candidate_evidence", False))
        want_ranking = (
            getattr(self.config, "ranking_mode", "off") == "filter"
        )
        if want_evidence or want_ranking:
            evidence_rows = self._compute_candidate_evidence_rows(
                aggregated, strict=want_ranking,
            )
            if want_evidence and evidence_rows is not None:
                from fin.analysis.candidate_evidence import (
                    write_candidate_evidence,
                )
                write_candidate_evidence(
                    Path(self.config.work_dir) / "candidate_evidence.tsv",
                    evidence_rows,
                )
            if want_ranking:
                if evidence_rows is None:
                    # An explicitly requested filter must never silently
                    # degrade into a no-op baseline run.
                    raise RuntimeError(
                        "ranking_mode=filter requires complete candidate "
                        "evidence, but the evidence computation failed"
                    )
                from fin.analysis.candidate_ranking import ranking_filter

                aggregated, _scores, dropped = ranking_filter(
                    aggregated,
                    evidence_rows,
                    threshold=getattr(
                        self.config, "ranking_threshold", None
                    ),
                )
                if dropped:
                    logger.info(
                        "Ranking filter removed %d candidates before "
                        "finalization", len(dropped),
                    )
        refit_effective = bool(getattr(
            self.config, "post_selection_refit_effective", False
        ))
        endpoint_routes: dict = {}
        endpoint_primary: dict = {}
        if refit_effective and getattr(self.config, "endpoint_refine", False):
            from fin.analysis.endpoint_refine import (
                apply_endpoint_splits,
                plan_endpoint_splits,
            )

            bam_ev = getattr(self, "_ranking_bam_evidence", None)
            if bam_ev is None:
                from fin.analysis.candidate_evidence import (
                    collect_ranking_bam_evidence,
                )
                bam_ev = collect_ranking_bam_evidence(self.config.bam_path)
                self._ranking_bam_evidence = bam_ev
            if not bam_ev.complete:
                raise RuntimeError(
                    "endpoint_refine requires a complete BAM end-evidence "
                    f"scan ({bam_ev.error})"
                )
            polya_read_ids = self._endpoint_polya_read_ids()
            tss_mode = getattr(self.config, "tss_evidence_mode", "off")
            background = 0.0
            if tss_mode != "off":
                background = self._estimate_degradation_background(
                    aggregated, bam_ev.read_ends,
                    bin_bp=int(getattr(self.config, "endpoint_window_bp", 25)),
                )
            plan = plan_endpoint_splits(
                aggregated,
                bam_ev.read_ends,
                tss_evidence_mode=tss_mode,
                background_hazard=background,
                window_bp=int(getattr(self.config, "endpoint_window_bp", 25)),
                min_end_reads=int(getattr(self.config, "endpoint_min_reads", 3)),
                min_pair_frac=float(
                    getattr(self.config, "endpoint_min_pair_frac", 0.15)
                ),
                max_splits=int(getattr(self.config, "endpoint_max_splits", 2)),
                polya_read_ids=polya_read_ids,
            )
            import json as _json
            diag_path = Path(self.config.work_dir) / "endpoint_refine.json"
            tmp = diag_path.with_suffix(".json.tmp")
            tmp.write_text(_json.dumps({
                "schema_version": 1,
                "polya_available": polya_read_ids is not None,
                "n_polya_reads": (
                    len(polya_read_ids) if polya_read_ids is not None else None
                ),
                "tss_evidence_mode": tss_mode,
                "degradation_background_hazard": background,
                "models_split": len(plan.replacements),
                "endpoint_states": plan.n_splits,
                "details": plan.details,
                "tss_verdicts": plan.tss_verdicts,
            }, indent=2, sort_keys=True) + "\n")
            tmp.replace(diag_path)
            if plan.replacements:
                aggregated = apply_endpoint_splits(aggregated, plan)
                endpoint_routes = plan.read_routes
                endpoint_primary = plan.primary
                logger.info(
                    "EndpointRefine split %d models into %d endpoint states "
                    "(polyA %s)",
                    len(plan.replacements), plan.n_splits,
                    "available" if polya_read_ids is not None else "unavailable",
                )
        if refit_effective:
            aggregated, snap_redirects = apply_junction_snap(
                self.config,
                aggregated,
                getattr(self, "_genome_fasta", None),
                return_redirects=True,
            )
            ledgers = responsibility_ledgers or []
            if aggregated and not ledgers:
                raise RuntimeError(
                    "post-selection refit is enabled but no responsibility "
                    "ledgers reached global finalization"
                )
            # Cross-interval extent-pair folds (unify_extent_pairs) are id -> id
            # redirects exactly like junction-snap merges: a read whose interval
            # ledger still names the folded loser is re-quantified on the winner.
            if extent_redirects:
                snap_redirects = {**extent_redirects, **dict(snap_redirects or {})}
            aggregated, diagnostics = refit_survivor_abundance(
                aggregated, ledgers, snap_redirects=snap_redirects,
                split_routes=endpoint_routes, split_primary=endpoint_primary,
                length_normalize=bool(getattr(
                    self.config, "tpm_length_normalize", True)),
            )
            # P0 strict-positive GTF gate, post-refit enforcement: the refit
            # can zero a GTF survivor that passed the pre-refit A4 gate; the
            # row carries no mass or reads, so removal conserves every
            # accounting invariant. Runs BEFORE the diagnostics write so the
            # JSON always describes the emitted survivor set.
            aggregated, zero_gtf_ids = enforce_gtf_positive_after_refit(
                aggregated, diagnostics,
                enabled=bool(getattr(
                    self.config, "require_gtf_positive_mass", False
                )),
            )
            if zero_gtf_ids:
                logger.info(
                    "Strict-positive GTF gate removed %d zero-mass GTF "
                    "survivor(s) after refit (gtf_zero_mass): %s",
                    len(zero_gtf_ids), ",".join(zero_gtf_ids[:8]),
                )
            if getattr(self.config, "selection_audit", False):
                append_selection_audit_global(
                    self.config.work_dir,
                    [("*", "junction_snap", loser, "", "", "fold", "snap_merge", winner)
                     for loser, winner in sorted(dict(snap_redirects or {}).items())
                     if loser not in (extent_redirects or {})]
                    + [("*", "post_refit", cid, "", "", "drop", "gtf_zero_mass", "")
                       for cid in zero_gtf_ids])
            write_abundance_refit_diagnostics(
                Path(self.config.work_dir) / "abundance_refit.json",
                diagnostics,
            )
        else:
            aggregated = apply_junction_snap(
                self.config, aggregated, getattr(self, "_genome_fasta", None)
            )
            # Defensive sweep on the non-refit path: A4 already applied the
            # strict boundary pre-snap and snapping only merges mass, so this
            # is expected to be a no-op; it guarantees the flag's output
            # contract (no zero-mass GTF survivor) regardless of path.
            aggregated, zero_gtf_ids = enforce_gtf_positive_after_refit(
                aggregated, None,
                enabled=bool(getattr(
                    self.config, "require_gtf_positive_mass", False
                )),
            )
            if zero_gtf_ids:
                logger.info(
                    "Strict-positive GTF gate removed %d zero-mass GTF "
                    "survivor(s) after snap (gtf_zero_mass): %s",
                    len(zero_gtf_ids), ",".join(zero_gtf_ids[:8]),
                )
            if getattr(self.config, "selection_audit", False):
                # snap merges are not returned on this path; only the removals
                # are recorded here
                append_selection_audit_global(
                    self.config.work_dir,
                    [("*", "post_snap", cid, "", "", "drop", "gtf_zero_mass", "")
                     for cid in zero_gtf_ids])
            if (
                getattr(self.config, "post_selection_refit", False)
                and getattr(
                    self.config, "post_selection_refit_disable_reason", None
                )
            ):
                write_abundance_refit_diagnostics(
                    Path(self.config.work_dir) / "abundance_refit.json",
                    {
                        "schema_version": 1,
                        "mode": "post_selection_survivor_renormalization",
                        "effective": False,
                        "disable_reason": self.config.post_selection_refit_disable_reason,
                        "structural_identity": "not_applicable",
                    },
                )
        if getattr(self.config, "selection_audit", False):
            append_selection_audit_global(
                self.config.work_dir,
                [("*", "final", cid, getattr(qr, "source", ""),
                  getattr(qr, "num_assigned_reads", ""), "keep", "", "")
                 for cid, qr in sorted(aggregated.items())])
        return finalize_outputs(
            self.config, aggregated, output_gtf, output_tsv, self._gtf_reader,
        )

    def _estimate_degradation_background(
        self, aggregated: Dict[str, QuantResult],
        read_ends: Mapping[str, tuple], *, bin_bp: int = 25,
    ) -> float:
        """Pooled 5'-termination hazard from models with no contained sibling.

        The null must describe degradation, so it is estimated only from
        models that no other surviving model is nested inside: a locus that
        already holds a shorter sibling would contaminate the background with
        that sibling's genuine TSS.
        """
        from fin.analysis.tss_evidence import (
            build_hazard_profile,
            genomic_to_offset,
            pooled_background_hazard,
            spliced_length,
        )

        def chain(ex):
            ex = sorted(ex)
            return tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))

        by_key: Dict[tuple, list] = {}
        for cid, qr in aggregated.items():
            by_key.setdefault((qr.chrom, qr.strand), []).append(cid)
        has_nested = set()
        for _key, ids in by_key.items():
            for a in ids:
                ca = chain(aggregated[a].exons)
                for b in ids:
                    if a == b:
                        continue
                    qb = aggregated[b]
                    qa = aggregated[a]
                    if qb.start <= qa.start and qa.end <= qb.end and \
                            (qa.end - qa.start) < (qb.end - qb.start):
                        has_nested.add(b)
        profiles = []
        for cid, qr in aggregated.items():
            if cid in has_nested or len(qr.exons) < 2:
                continue
            offs = []
            plus = qr.strand == "+"
            for rid in getattr(qr, "assigned_read_ids", ()) or ():
                sp = read_ends.get(rid)
                if sp is None:
                    continue
                r5 = sp[0] if plus else sp[1]
                o = genomic_to_offset(r5, qr.exons, qr.strand)
                if o is not None:
                    offs.append(o)
            if len(offs) >= 20:
                profiles.append(build_hazard_profile(
                    offs, spliced_length(qr.exons), bin_bp=bin_bp))
        bg = pooled_background_hazard(profiles)
        logger.info(
            "TSS evidence: degradation background %.6f per %dbp bin from %d "
            "nesting-free models", bg, bin_bp, len(profiles),
        )
        return bg

    def _endpoint_polya_read_ids(self) -> Optional[set]:
        """Poly(A)-confident read IDs for EndpointRefine TES support.

        Reuses the same krill whole-read polyA pass and confidence predicate
        as the polyA+5' filter. Returns None (documented degraded mode: TES
        modes rely on end sharpness alone) when no signal file is configured
        or krill yields nothing; EndpointRefine logs which mode ran.
        """
        if not getattr(self.config, "signal_path", None):
            return None
        try:
            from fin.analysis.quantification import polya_read_passes
            from fin.scoring.polya import compute_polya

            read_seqs, _read_ends = self._fetch_read_seqs_and_ends()
            polya_map = compute_polya(
                read_seqs,
                self.config.signal_path,
                pore=self.config.krill_pore,
                use_gpu=self.config.use_gpu,
            )
            if not polya_map:
                return None
            min_len = float(getattr(self.config, "min_polya_length", 15.0))
            return {
                rid for rid in polya_map
                if polya_read_passes(rid, polya_map, min_len)
            }
        except Exception:
            logger.exception(
                "endpoint polyA pass failed; TES support degrades to "
                "end-sharpness only"
            )
            return None

    def _compute_candidate_evidence_rows(
        self, aggregated: Dict[str, QuantResult], *, strict: bool = False
    ) -> Optional[list]:
        """Observable-feature rows for the post-selection survivor set.

        Computed BEFORE junction snapping and the abundance refit, so rows
        show exactly the feature view the selection-stage scorer sees. One
        whole-BAM pass supplies junction support and read ends together.

        ``strict=False`` (pure audit): failures and incomplete BAM scans are
        logged and tolerated (partial evidence may still be written for
        inspection). ``strict=True`` (the ranking filter): an incomplete BAM
        scan or any failure returns None so the caller refuses to filter -
        undercounted junction support must never depress ranking scores.
        """
        try:
            from fin.analysis.candidate_evidence import (
                collect_ranking_bam_evidence,
                compute_candidate_evidence,
            )
            from fin.candidates.canonical import parse_motifs

            bam_ev = collect_ranking_bam_evidence(self.config.bam_path)
            self._ranking_bam_evidence = bam_ev
            if strict and not bam_ev.complete:
                logger.error(
                    "ranking evidence BAM scan incomplete (%s); refusing to "
                    "rank on partial junction support", bam_ev.error,
                )
                return None
            return compute_candidate_evidence(
                aggregated,
                junction_support=bam_ev.junction_support,
                read_ends=bam_ev.read_ends,
                genome=self._genome_fasta,
                canonical_motifs=parse_motifs(
                    getattr(self.config, "canonical_motifs", ()) or ()
                ),
                end_window_bp=int(getattr(self.config, "fulllen_window_bp", 25)),
            )
        except Exception:
            logger.exception("candidate-evidence computation failed")
            return None

    def process_interval(
        self, interval: GenomicInterval
    ) -> Optional[Tuple[List[QuantResult], Optional[ResponsibilityLedger]]]:
        """Discover candidates and quantify a single interval via ``quant_mode``.

        Returns:
            ``(results, ledger)`` for a quantified interval. The ledger is
            populated only when post-selection refit is effectively enabled.
            Returns ``None`` when the interval has no candidates or no reads.
        """
        work_dir = Path(self.config.work_dir) / interval.region_string.replace(":", "_").replace("-", "_")
        work_dir.mkdir(parents=True, exist_ok=True)

        # Get genome sequence for this chromosome
        chrom_seq = ""
        if self._genome_fasta and interval.chrom in self._genome_fasta:
            chrom_seq = self._genome_fasta[interval.chrom]

        # --- Phase 1: Candidate discovery ---
        candidate_set = discover_candidates(
            interval=interval,
            bam_path=self.config.bam_path,
            gtf_reader=self._gtf_reader,
            genome_fasta=chrom_seq,
            threshold=self.config.three_prime_threshold,
            min_novel_reads=self.config.min_novel_reads,
            chain_cluster=getattr(self.config, "chain_cluster_discovery", False),
            chain_cluster_wobble_bp=getattr(self.config, "chain_cluster_wobble_bp", 6),
            chain_cluster_cassette_max_exon_bp=getattr(
                self.config, "chain_cluster_cassette_max_exon_bp", 70),
            chain_cluster_fold_monoexon=(
                getattr(self.config, "chain_cluster_fold_monoexon", False)
                # Defer the generation-time mono fold ONLY when post-EM resolution will
                # actually run (mono_resolve_post_em is m2_em-only). For argmax/m1_em/cluster
                # there is no post-EM mono step, so keep the generation fold to suppress FPs.
                and not (getattr(self.config, "mono_resolve_post_em", False)
                         and getattr(self.config, "quant_mode", "m2_em") == "m2_em")),
            chain_cluster_fold_span_guard=getattr(
                self.config, "chain_cluster_fold_span_guard", False),
            clustering=getattr(self.config, "clustering", "read_chains"),
            # a per-read x per-model scan: computed only when the support-gate
            # attestation will read it, so the default path pays nothing
            collect_compatible_depths=bool(
                getattr(self.config, "m2_support_gate_attest_reads", 0) > 0),
            read_extent_mode_tol_bp=int(getattr(
                self.config, "read_extent_mode_tol_bp", 0)),
            junction_correction=self._make_junction_correction(interval, chrom_seq),
            junction_correction_apply=self._junction_correction_apply(),
            read_extent_contradicted=self._read_extent_contradicted(),
            canonical_search_bp=self.config.canonical_search_bp,
            max_chains_per_read=self.config.max_chains_per_read,
            canonical_motifs=self.config.canonical_motifs,
            denovo_wobble_tol=getattr(self.config, "denovo_wobble_tol", 0),
            denovo_wobble_shadow_ratio=getattr(
                self.config, "denovo_wobble_shadow_ratio", 1.0),
            denovo_graph=getattr(self.config, "denovo_graph", False),
            denovo_graph_tol=getattr(self.config, "denovo_graph_tol", 6),
            denovo_graph_min_edge_reads=getattr(
                self.config, "denovo_graph_min_edge_reads", 2),
            denovo_graph_tss_brake=getattr(
                self.config, "denovo_graph_tss_brake", True),
            denovo_graph_tss_tol=getattr(self.config, "denovo_graph_tss_tol", 20),
            denovo_graph_tss_min_reads=getattr(
                self.config, "denovo_graph_tss_min_reads", 3),
            denovo_graph_tss_frac=getattr(
                self.config, "denovo_graph_tss_frac", 0.4),
        )

        # Selection audit (opt-in): every candidate's fate through the gates,
        # written as one shard per interval at the end of this method. Rows are
        # collected here; nothing below reads them.
        audit = bool(getattr(self.config, "selection_audit", False))
        audit_key = interval.region_string
        audit_rows: List[tuple] = []
        if audit:
            audit_rows += audit_rows_for_candidates(
                audit_key, "discovered", candidate_set.candidates)

        # Always-on record of the interval junction registry (raw / indel-clean
        # support per junction). Evidence only; no gate reads this file.
        if getattr(candidate_set, "junction_support", None):
            from fin.analysis.junction_support import write_support_shard

            write_support_shard(self.config.work_dir, interval, candidate_set.junction_support)

        # --- Phase 1.5: Fusion candidate augmentation (optional) ---
        if self.config.fusion_enabled:
            candidate_set = self._augment_with_fusion_candidates(
                candidate_set, interval
            )

        # --- Phase 1.6: Canonical-motif gate (Stage B) ---
        # Drop NOVEL multi-exon candidates whose intron chain isn't all-canonical.
        # GTF-passthrough and fusion candidates are exempt (annotated junctions
        # are trusted); mono candidates trivially pass.
        if getattr(self.config, "canonical_gate", False):
            gate_outcomes = self._apply_canonical_gate(candidate_set, chrom_seq)
            if audit:
                audit_rows += audit_rows_for_outcomes(
                    audit_key, "canonical_gate", gate_outcomes or [])

        # --- Phase 1.7: Junction-dominance gate (PRE-EM, junction-first) ---
        # Drop NOVEL multi-exon candidates with a junction that is weak
        # (< min_reads observed) or non-dominant (a different observed junction
        # within window bp carries strictly more reads). Removes multi-read wobble
        # shadows before EM so they never compete for reads. gtf/fusion/mono
        # exempt. OFF (default) -> no drops (byte-identical).
        if getattr(self.config, "junction_dominance_filter", False):
            gate_outcomes = self._apply_junction_dominance_gate(candidate_set, interval)
            if audit:
                audit_rows += audit_rows_for_outcomes(
                    audit_key, "junction_dominance_gate", gate_outcomes or [])

        if candidate_set.num_candidates == 0:
            logger.info("No candidates for interval %s", interval.region_string)
            if audit:
                write_selection_audit_shard(self.config.work_dir, interval, audit_rows)
            return None

        if not candidate_set.read_ids:
            logger.info("No reads for interval %s", interval.region_string)
            if audit:
                audit_rows.append((audit_key, "interval", "*", "", 0,
                                   "drop", "no_reads_in_interval", ""))
                write_selection_audit_shard(self.config.work_dir, interval, audit_rows)
            return None

        # Sort once; canonical read axis for every distance matrix.
        read_ids = sorted(candidate_set.read_ids)
        candidate_ids = candidate_set.candidate_ids()

        # Quantification engine dispatch (quant_mode). All three modes are
        # krill-only (no f5c CLI).
        quant_mode = getattr(self.config, "quant_mode", "m2_em")
        ledger: Optional[ResponsibilityLedger] = None
        self._audit_interval_rows = [] if audit else None
        if quant_mode == "argmax":
            results = self._quant_argmax_keep(candidate_set, read_ids, interval)
        elif quant_mode == "m1_em":
            results = self._quant_m1_em(candidate_set, read_ids, interval)
        elif quant_mode == "m2_em":
            results, ledger = self._quant_m2_em(
                candidate_set, read_ids, interval
            )
        elif quant_mode == "cluster":
            results = self._quant_cluster(candidate_set, read_ids, interval)
        else:
            raise ValueError(f"unknown quant_mode: {quant_mode!r}")
        if audit:
            audit_rows += self._audit_interval_rows or []
            self._audit_interval_rows = None
            audit_rows += [(audit_key, "interval_survivor", qr.candidate_id,
                            getattr(qr, "source", ""), getattr(qr, "num_assigned_reads", ""),
                            "keep", "", "") for qr in (results or [])]
            write_selection_audit_shard(self.config.work_dir, interval, audit_rows)
        if not results and ledger is None:
            return None

        # Full-length end-coherence: compute fulllen_frac per candidate over its
        # assigned reads (the non-circular population; the fulllen METRIC itself
        # uses BAM spans only — no signal). Gated by the same switches as the
        # drop in _finalize_and_write so a filter-disabled run pays no BAM-fetch
        # cost and leaves the -1.0 sentinel untouched.
        _filter_on = getattr(self.config, "enable_score_filter", True)
        if (
            results
            and _filter_on
            and getattr(self.config, "min_fulllen_fraction", 0.0) > 0.0
        ):
            self._annotate_fulllen_frac(results, interval)
        return results, ledger

    def _quant_argmax_first(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> List[QuantResult]:
        """M1-first hard argmin assignment (production default; no signal, no EM).

        Each read is aligned against every candidate with the configured mappy
        preset (``MAPPY_PRESET``, default map-ont); the reconstructed map-ont AS
        (``score_hit``, single-indel > cap → rejected) is taken as the cell
        score. The read is assigned to its single best-AS candidate; AS TIES are
        broken by lowest candidate index — which is the implicit GTF prior, since
        ``discover_candidates`` returns ``gtf_passthrough + collapsed_novel +
        fusion`` (GTF first). Abundance = integer read count.

        This mirrors the ablation harness ``_build_m1`` + ``_quant_tie(mode=
        "first")`` so production and the harness produce the same assignment for
        the same candidate set.
        """
        import mappy

        from fin.scoring.mappy_preset import get_m1_preset
        from fin.scoring.mappy_score import score_hit

        cand_list = list(candidate_set.candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        n_c = len(cand_list)

        preset = get_m1_preset()
        aligners = [
            mappy.Aligner(seq=c.sequence, preset=preset) if c.sequence else None
            for c in cand_list
        ]

        counts = [0] * n_c
        assigned: List[List[str]] = [[] for _ in range(n_c)]
        for rid in read_ids:
            seq = read_sequences.get(rid)
            if not seq:
                continue
            best_as = None      # highest reconstructed AS over all candidates
            best_j = -1
            for j, aln in enumerate(aligners):
                if aln is None:
                    continue
                cell = None
                for h in aln.map(seq):
                    v = score_hit(h)
                    if v is None:
                        continue
                    if cell is None or v > cell:
                        cell = v
                if cell is None:
                    continue
                # strict '>' keeps the LOWEST-index candidate on an AS tie.
                if best_as is None or cell > best_as:
                    best_as = cell
                    best_j = j
            if best_j < 0:
                continue
            counts[best_j] += 1
            assigned[best_j].append(rid)

        quant_results: List[QuantResult] = []
        for j, cand in enumerate(cand_list):
            has = counts[j] > 0
            qr = QuantResult(
                candidate_id=cand.candidate_id,
                abundance=float(counts[j]),
                confidence=1.0 if has else 0.0,
                num_assigned_reads=len(assigned[j]),
                source=cand.source,
                chrom=cand.chrom,
                strand=cand.strand,
                start=cand.start,
                end=cand.end,
                exons=_exons_from_candidate(cand),
                family_id=cand.family_id,
                complete_chain_weakest_indel_clean_depth=getattr(
                    cand, "complete_chain_weakest_indel_clean_depth", 0),
                assigned_read_ids=tuple(assigned[j]),
                breakpoint_left=cand.breakpoint_left,
                breakpoint_right=cand.breakpoint_right,
                fusion_junction=cand.fusion_junction,
            )
            qr.max_R = 1.0 if has else 0.0
            quant_results.append(qr)

        logger.info(
            "R1 M1-first interval %s: %d reads -> %d candidates",
            interval.region_string, len(read_ids), len(cand_list),
        )
        return quant_results

    def _quant_argmax_keep(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> List[QuantResult]:
        """M1-keep SPLIT assignment (production default; no signal, no EM).

        Each read is aligned against every candidate with the configured mappy
        preset; the reconstructed AS (``score_hit``) is the cell score. The read
        is assigned to ALL of its simultaneously-best-AS candidates (within a
        1e-9 tolerance). Abundance is split 1/K across the K tied winners so the
        per-read mass is conserved (== 1.0 per read total).

        Mirrors the ablation harness ``_build_m1`` + ``_quant_tie(mode="split")``
        (config "M1-split"): production and harness produce the same split
        assignment for the same candidate set.
        """
        import mappy

        from fin.scoring.mappy_preset import get_m1_preset
        from fin.scoring.mappy_score import score_hit

        cand_list = list(candidate_set.candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        n_c = len(cand_list)

        preset = get_m1_preset()
        aligners = [
            mappy.Aligner(seq=c.sequence, preset=preset) if c.sequence else None
            for c in cand_list
        ]

        # Optional M2 tie resolution (Stage M2-1): when a read is simultaneously
        # best-AS across >=2 candidates, the validated junction-window mean-NLL
        # metric picks the single true wobble sibling and (if confident) takes the
        # read's FULL mass instead of the 1/K split. Lazily build a shared non-HMM
        # krill aligner; auto-skip the whole leg if signal is absent.
        m2_krill = None
        m2_gpu = False
        m2_on = bool(getattr(self.config, "m2_tiebreak", False)) and bool(
            self.config.signal_path
        )
        m2_threads = 0
        if m2_on:
            try:
                import krill

                from fin.scoring.krill_aligner import (
                    krill_thread_count,
                    make_krill_aligner,
                )

                m2_threads = krill_thread_count()
                m2_krill, m2_gpu = make_krill_aligner(
                    krill, self.config.krill_pore, self.config.use_gpu,
                    hmm_confidence=False, num_thread=m2_threads,
                )
            except Exception as exc:  # krill not importable -> keep 1/K split
                logger.warning("M2 tiebreak disabled (krill import failed): %s", exc)
                m2_krill = None
            if m2_krill is None:  # signal stack unavailable -> keep 1/K split
                logger.warning("M2 tiebreak disabled (krill init failed)")
                m2_on = False

        counts = [0.0] * n_c
        assigned: List[List[str]] = [[] for _ in range(n_c)]
        # Track the max single-read weight any candidate receives (== max_R).
        max_weight = [0.0] * n_c
        n_m2_override = 0
        for rid in read_ids:
            seq = read_sequences.get(rid)
            if not seq:
                continue
            cells: List[Optional[float]] = [None] * n_c
            best_as = None
            for j, aln in enumerate(aligners):
                if aln is None:
                    continue
                cell = None
                for h in aln.map(seq):
                    v = score_hit(h)
                    if v is None:
                        continue
                    if cell is None or v > cell:
                        cell = v
                if cell is None:
                    continue
                cells[j] = cell
                if best_as is None or cell > best_as:
                    best_as = cell
            if best_as is None:
                continue
            tied = [
                j for j in range(n_c)
                if cells[j] is not None and cells[j] >= best_as - 1e-9
            ]
            if not tied:
                continue
            # M2 single-winner override on a genuine (>=2) tie set.
            if m2_on and len(tied) >= 2:
                from fin.scoring.m2_junction_nll import m2_resolve_tie

                tied_cands = [cand_list[j] for j in tied]
                tied_aligners = [aligners[j] for j in tied]
                best_local, margin, scored_local = m2_resolve_tie(
                    rid, seq, tied_cands, self.config.signal_path,
                    pore=self.config.krill_pore,
                    junction_k=self.config.m2_tiebreak_junction_k,
                    krill_aligner=m2_krill, mappy_aligners=tied_aligners,
                    return_scored=True, use_gpu=m2_gpu, num_thread=m2_threads,
                )
                if best_local is not None and margin >= self.config.m2_tiebreak_margin:
                    j = tied[best_local]
                    counts[j] += 1.0
                    assigned[j].append(rid)
                    if 1.0 > max_weight[j]:
                        max_weight[j] = 1.0
                    n_m2_override += 1
                    continue
                # Score-gated fallback split (mk1): restrict the 1/K split to the
                # candidates eventalign could score. If none scored, keep the full
                # tied set (== M1). Recall-safe: the read is never dropped.
                if (
                    getattr(self.config, "m2_tie_scoregate_split", False)
                    and scored_local
                ):
                    tied = [tied[k] for k in scored_local]
            w = 1.0 / len(tied)
            for j in tied:
                counts[j] += w
                assigned[j].append(rid)
                if w > max_weight[j]:
                    max_weight[j] = w

        quant_results: List[QuantResult] = []
        for j, cand in enumerate(cand_list):
            has = counts[j] > 0.0
            qr = QuantResult(
                candidate_id=cand.candidate_id,
                abundance=float(counts[j]),
                confidence=1.0 if has else 0.0,
                num_assigned_reads=len(assigned[j]),
                source=cand.source,
                chrom=cand.chrom,
                strand=cand.strand,
                start=cand.start,
                end=cand.end,
                exons=_exons_from_candidate(cand),
                family_id=cand.family_id,
                complete_chain_weakest_indel_clean_depth=getattr(
                    cand, "complete_chain_weakest_indel_clean_depth", 0),
                assigned_read_ids=tuple(assigned[j]),
                breakpoint_left=cand.breakpoint_left,
                breakpoint_right=cand.breakpoint_right,
                fusion_junction=cand.fusion_junction,
            )
            qr.max_R = max_weight[j]
            quant_results.append(qr)

        logger.info(
            "R1 M1-keep(split) interval %s: %d reads -> %d candidates "
            "(M2 overrides=%d)",
            interval.region_string, len(read_ids), len(cand_list), n_m2_override,
        )
        return quant_results

    def _quant_cluster(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> List[QuantResult]:
        """Within-cluster read assignment (quant_mode="cluster").

        Each generation cluster (``CandidateSet.clusters``) is scored in isolation:
        a read is aligned ONLY against that cluster's members (M1), straddling ties
        are resolved by the summed-LLR junction signal (M2), and containment/
        truncation ties fall to the cluster's main peak. Per-cluster assignment is
        delegated to :func:`fin.pipeline.cluster_quant.assign_cluster_reads` (pure/
        signal-agnostic); the M2 tiebreak is injected via ``m2_resolve_tie``. Cluster
        results are folded back into per-candidate abundance/survival and emitted as
        QuantResults (mirrors ``_quant_argmax_keep``).
        """
        import mappy

        from fin.pipeline.cluster_quant import assign_cluster_reads
        from fin.scoring.m2_junction_nll import m2_resolve_tie
        from fin.scoring.mappy_preset import get_m1_preset
        from fin.scoring.mappy_score import score_hit

        cand_list = list(candidate_set.candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        n_c = len(cand_list)

        # candidate_set.clusters holds candidate_ids per generation cluster. Translate
        # to CURRENT list indices (ids dropped by post-discovery gates are skipped).
        idx_of = {c.candidate_id: j for j, c in enumerate(cand_list)}
        if candidate_set.clusters is None:
            # Not chain-cluster mode: run unscoped over the whole candidate list.
            logger.warning(
                "quant_mode=cluster: interval %s has no cluster grouping "
                "(clusters=None); treating all %d candidates as one cluster",
                interval.region_string, n_c,
            )
            clusters = [list(range(n_c))]
        else:
            clusters = [
                [idx_of[cid] for cid in ids if cid in idx_of]
                for ids in candidate_set.clusters
            ]
            # Candidates added AFTER discovery (e.g. fusion augmentation) are in
            # cand_list but not in any generation cluster; group them into ONE extra
            # cluster so they are scored (never silently zero) AND a read shared by
            # several of them competes once (one singleton each would let the same
            # read take weight 1.0 in every singleton -> inflated abundance).
            referenced = {j for ml in clusters for j in ml}
            extra = [j for j in range(n_c) if j not in referenced]
            if extra:
                clusters.append(extra)

        # Read genomic spans over the interval (primary mapped alignments only;
        # first occurrence per read id wins). Mirrors _annotate_fulllen_frac.
        spans: Dict[str, Tuple[int, int]] = {}
        with pysam.AlignmentFile(self.config.bam_path, "rb") as bam:
            for r in bam.fetch(interval.chrom, interval.start, interval.end):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                rid = r.query_name
                if rid is None or rid in spans:
                    continue
                rs = r.reference_start
                re = r.reference_end
                if rs is None or re is None:
                    continue
                spans[rid] = (int(rs), int(re))

        # One shared non-HMM krill aligner for the whole interval; when signal is
        # absent m2_resolve stays None (straddling ties defer, which is fine).
        krill_aligner = None
        eff_gpu = False
        krill_threads = 0
        if self.config.signal_path and getattr(self.config, "cluster_use_m2", True):
            try:
                import krill

                from fin.scoring.krill_aligner import (
                    krill_thread_count,
                    make_krill_aligner,
                )

                krill_threads = krill_thread_count()
                krill_aligner, eff_gpu = make_krill_aligner(
                    krill, self.config.krill_pore, self.config.use_gpu,
                    hmm_confidence=False, num_thread=krill_threads,
                )
            except Exception as exc:  # krill unavailable -> defer straddling ties
                logger.warning(
                    "quant_mode=cluster: krill unavailable, M2 tiebreak off: %s", exc)
                krill_aligner = None

        preset = get_m1_preset()
        read_id_set = set(read_ids)

        # Accumulators over ALL candidates (global index j).
        counts = [0.0] * n_c
        assigned: List[List[str]] = [[] for _ in range(n_c)]
        max_weight = [0.0] * n_c
        survivor = [False] * n_c

        # Genome chromosome sequence for 5'-TSS short-isoform recovery (empty ->
        # recovery skipped for this interval).
        genome_seq = (
            self._genome_fasta[interval.chrom]
            if self._genome_fasta and interval.chrom in self._genome_fasta
            else ""
        )
        shadows_map = candidate_set.shadows or {}

        tot_unique = tot_containment = tot_m2_confident = tot_deferred = 0
        recovered_total = 0
        for members_idx in clusters:
            if not members_idx:
                continue
            member_candidates = [cand_list[i] for i in members_idx]
            local_aligners = [
                mappy.Aligner(seq=c.sequence, preset=preset) if c.sequence else None
                for c in member_candidates
            ]
            member_chains = [
                tuple(cand_list[i].intron_chain.introns) for i in members_idx
            ]

            # Cluster reads = union of members' supporting reads, restricted to this
            # interval's read set and reads that carry a sequence.
            cluster_read_ids: set = set()
            for i in members_idx:
                cluster_read_ids.update(cand_list[i].supporting_read_ids)
            cluster_read_ids &= read_id_set
            cluster_read_ids = {
                rid for rid in cluster_read_ids if read_sequences.get(rid)
            }

            # M1 within-cluster: best AS per (read, member), local member index.
            as_rows: Dict[str, Dict[int, float]] = {}
            for rid in cluster_read_ids:
                seq = read_sequences.get(rid)
                if not seq:
                    continue
                row: Dict[int, float] = {}
                for local_j, aln in enumerate(local_aligners):
                    if aln is None:
                        continue
                    cell = None
                    for h in aln.map(seq):
                        v = score_hit(h)
                        if v is None:
                            continue
                        if cell is None or v > cell:
                            cell = v
                    if cell is not None:
                        row[local_j] = cell
                if row:
                    as_rows[rid] = row

            # M2 resolver for straddling ties (summed-LLR); None when no krill.
            m2_resolve = None
            if krill_aligner is not None:
                def m2_resolve(rid, tie_local, _mc=member_candidates,
                               _la=local_aligners):
                    tied_cands = [_mc[j] for j in tie_local]
                    tied_aligners = [_la[j] for j in tie_local]
                    # Honor the configured tight-window metric so the
                    # experimental sqrt_count variant behaves consistently in
                    # cluster mode; anything else keeps the validated summed
                    # contrast this resolver was built around.
                    _metric = (
                        self.config.m2_metric
                        if getattr(self.config, "m2_metric", "summed_llr")
                        in ("summed_llr", "sqrt_count_mean_llr")
                        else "summed_llr"
                    )
                    out = m2_resolve_tie(
                        rid, read_sequences[rid], tied_cands,
                        self.config.signal_path, pore=self.config.krill_pore,
                        krill_aligner=krill_aligner, mappy_aligners=tied_aligners,
                        use_gpu=eff_gpu, num_thread=krill_threads,
                        metric=_metric,
                    )
                    best, margin = out
                    if best is None:
                        return None
                    return (tie_local[best], margin)

            a = assign_cluster_reads(
                member_chains, as_rows,
                {rid: spans[rid] for rid in as_rows if rid in spans},
                m2_resolve=m2_resolve,
                min_support=self.config.cluster_min_support,
                m1_tie_margin=self.config.cluster_m1_tie_margin,
            )

            # Fold local results back to GLOBAL candidate indices.
            for rid, row in a.weights.items():
                for local_j, w in row.items():
                    gi = members_idx[local_j]
                    counts[gi] += w
                    assigned[gi].append(rid)
                    if w > max_weight[gi]:
                        max_weight[gi] = w
            for local_j in a.survivors:
                survivor[members_idx[local_j]] = True

            tot_unique += a.n_unique
            tot_containment += a.n_containment
            tot_m2_confident += a.n_m2_confident
            tot_deferred += a.n_deferred

            # 5'-TSS short-isoform recovery: each EM-confirmed survivor member may
            # have folded SHADOW sub-chains (exact 5'-truncations pooled into it at
            # generation). Re-test each shadow's truncation boundary against the
            # member's read-5'-end pileup; a SHARP peak (real TSS, not a degradation
            # ramp) re-activates the short isoform as a novel candidate. Mass is
            # conserved: the recovered peak-excess is subtracted from the parent.
            if genome_seq:
                for lj in sorted(a.survivors):
                    gj = members_idx[lj]
                    cand = cand_list[gj]
                    cid = cand.candidate_id
                    strand = cand.strand
                    shadows = shadows_map.get(cid)
                    if not shadows:
                        continue

                    # 5' end of a read (strand-aware): ref_start for +, ref_end for -.
                    def _p5(rid, _s=strand):
                        rs, re = spans[rid]
                        return rs if _s == "+" else re

                    # Member's EM-weighted assigned reads (this member's column).
                    member_reads = {
                        rid: row[lj]
                        for rid, row in a.weights.items()
                        if lj in row
                    }
                    ends = [
                        (_p5(rid), w)
                        for rid, w in member_reads.items()
                        if rid in spans
                    ]
                    if len(ends) < 8:
                        continue

                    # Build TSS proposals (median 5' end of each shadow's reads).
                    proposals: List[int] = []
                    prop_shadow: Dict[int, Tuple[tuple, tuple]] = {}
                    for schain, sreads in shadows:
                        s_ends = sorted(_p5(rid) for rid in sreads if rid in spans)
                        if not s_ends:
                            continue
                        tss = s_ends[len(s_ends) // 2]
                        prev = prop_shadow.get(tss)
                        if prev is None or len(sreads) > len(prev[1]):
                            prop_shadow[tss] = (schain, sreads)
                            if prev is None:
                                proposals.append(tss)

                    if not proposals:
                        continue

                    from fin.candidates.isoform_recovery import recover_5p_peaks

                    peaks = recover_5p_peaks(ends, list(proposals))
                    for peak in peaks:
                        # Nearest proposal within 20 bp.
                        near = min(
                            prop_shadow.keys(),
                            key=lambda p: abs(p - peak.pos),
                        )
                        if abs(near - peak.pos) > 20:
                            continue
                        schain, sreads = prop_shadow[near]

                        # Short-isoform genomic span: 5' end -> TSS, 3' end = parent's.
                        if strand == "+":
                            start, end = peak.pos, cand.end
                        else:
                            start, end = cand.start, peak.pos
                        if end <= start:
                            continue

                        from fin.candidates.discovery import (
                            _build_spliced_sequence,
                            _generate_novel_id,
                        )
                        from fin.candidates.dataclasses import (
                            IntronChain,
                            TranscriptCandidate,
                        )

                        short_chain = IntronChain(introns=schain)
                        seq = _build_spliced_sequence(
                            genome_seq, start, end, short_chain, strand,
                        )
                        if not seq:
                            continue

                        new_id = _generate_novel_id(
                            cand.chrom, strand, start, end, short_chain
                        )
                        new_cand = TranscriptCandidate(
                            candidate_id=new_id,
                            intron_chain=short_chain,
                            three_prime_pos=(end if strand == "+" else start),
                            sequence=seq,
                            source="novel",
                            supporting_read_ids=set(sreads),
                            chrom=cand.chrom,
                            strand=strand,
                            start=start,
                            end=end,
                            family_id=cand.family_id,
                        )
                        # Append in lockstep with the parallel accumulators so the
                        # emission loop (which enumerates cand_list) picks it up.
                        cand_list.append(new_cand)
                        counts.append(peak.excess)
                        assigned.append(list(sreads))
                        max_weight.append(1.0)
                        survivor.append(True)
                        # Conserve mass: remove the recovered excess from the parent.
                        counts[gj] = max(0.0, counts[gj] - peak.excess)
                        recovered_total += 1

        if recovered_total:
            logger.info(
                "R1 cluster interval %s: recovered %d 5'-TSS short isoforms",
                interval.region_string, recovered_total,
            )

        quant_results: List[QuantResult] = []
        for j, cand in enumerate(cand_list):
            # Honor the per-cluster survival gate: a member whose total assigned
            # weight is below cluster_min_support is NOT a survivor -> zero it so the
            # threshold is a real gate (not overridden here to confidence 1.0).
            kept = survivor[j]
            qr = QuantResult(
                candidate_id=cand.candidate_id,
                abundance=float(counts[j]) if kept else 0.0,
                confidence=1.0 if kept else 0.0,
                num_assigned_reads=len(assigned[j]) if kept else 0,
                source=cand.source,
                chrom=cand.chrom,
                strand=cand.strand,
                start=cand.start,
                end=cand.end,
                exons=_exons_from_candidate(cand),
                family_id=cand.family_id,
                complete_chain_weakest_indel_clean_depth=getattr(
                    cand, "complete_chain_weakest_indel_clean_depth", 0),
                assigned_read_ids=tuple(assigned[j]) if kept else (),
                breakpoint_left=cand.breakpoint_left,
                breakpoint_right=cand.breakpoint_right,
                fusion_junction=cand.fusion_junction,
            )
            qr.max_R = max_weight[j]
            quant_results.append(qr)

        logger.info(
            "R1 cluster interval %s: %d clusters, %d reads -> %d candidates "
            "(unique=%d containment=%d m2_confident=%d deferred=%d)",
            interval.region_string, len(clusters), len(read_ids), len(cand_list),
            tot_unique, tot_containment, tot_m2_confident, tot_deferred,
        )
        return quant_results

    def _annotate_fulllen_frac(
        self,
        results: List[QuantResult],
        interval: GenomicInterval,
    ) -> None:
        """Store fulllen_frac on NOVEL multi-exon QuantResults (in place).

        full-length read support: the fraction of a candidate's argmax-assigned
        reads whose primary genomic 5' AND 3' alignment ends both fall within
        ``fulllen_window_bp`` of the candidate's genomic 5'/3' ends. Signal-free
        (BAM primary spans only). A single fetch over the interval builds the
        read-end map shared by every candidate. Candidates with fewer than
        ``fulllen_min_reads`` reads carrying a span keep the -1.0 sentinel
        (unreachable -> never dropped). gtf/fusion/mono candidates are skipped.
        """
        novel_multi = [
            qr for qr in results
            if qr.source == "novel" and len(qr.exons) >= 3
        ]
        if not novel_multi:
            return
        window = getattr(self.config, "fulllen_window_bp", 25)
        min_reads = getattr(self.config, "fulllen_min_reads", 4)
        # Single BAM fetch over the interval -> genomic read-end map (primary
        # mapped alignments only; first occurrence per read id wins).
        read_ends: Dict[str, Tuple[int, int]] = {}
        with pysam.AlignmentFile(self.config.bam_path, "rb") as bam:
            for r in bam.fetch(interval.chrom, interval.start, interval.end):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                rid = r.query_name
                if rid is None or rid in read_ends:
                    continue
                rs = r.reference_start
                re = r.reference_end
                if rs is None or re is None:
                    continue
                read_ends[rid] = (int(rs), int(re))
        for qr in novel_multi:
            qr.fulllen_frac = compute_fulllen_frac(
                qr, read_ends, window, min_reads
            )

    # ---- sample-level evidence for the existence test -----------------------
    _SAMPLE_EVIDENCE_NAME = "_sample_evidence.json"      # legacy basename only

    def _invalidate_sample_caches(self) -> None:
        """Drop every in-memory view of the sample evidence. Called whenever the
        file on disk is (re)written, so a retry after a failed polyA pass does
        not keep serving the abstaining maps to the intervals."""
        for attr in ("_sample_evidence_cache", "_nuisance_evidence_cache",
                     "_read_extent_hook_cache"):
            if hasattr(self, attr):
                delattr(self, attr)

    def _open_gtf_reader(self):
        """Parse the configured annotation and remember its fingerprint, so a
        later run of this runner can tell that the annotation changed."""
        from fin.io.io_gtf import GTFReader

        before = input_fingerprint(self.config.gtf_path)
        reader = GTFReader(self.config.gtf_path)
        reader.open()
        reader.parse()
        if input_fingerprint(self.config.gtf_path) != before:
            reader.close()
            raise RuntimeError(
                "the annotation changed while it was being parsed "
                f"({self.config.gtf_path}); rerun on a stable input")
        # the fingerprint describes exactly the bytes this reader holds
        self._gtf_fingerprint = before
        return reader

    def _close_jc_context(self) -> None:
        """Drop the junction-correction context and close the slow5 handle it
        opened (its ``reader`` is a closure over that handle, not a file object)."""
        s5 = getattr(self, "_jc_slow5", None)
        if s5 is not None and hasattr(s5, "close"):
            try:
                s5.close()
            except Exception:  # noqa: BLE001 - closing a stale handle is best effort
                logger.debug("junction-correction slow5 close failed", exc_info=True)
        self._jc_slow5 = None
        if hasattr(self, "_jc_ctx"):
            del self._jc_ctx

    def _reset_input_snapshot(self) -> None:
        """Re-fingerprint the inputs at a run boundary. A runner may be reused
        after its BAM / signal / annotation changed; the identity must then be
        recomputed and EVERY view derived from the old inputs dropped -- the
        endpoint evidence, the annotation TSS/TES lists, and the junction
        correction context (which holds an aligner and an open signal reader)."""
        if hasattr(self, "_evidence_identity_cache"):
            del self._evidence_identity_cache
        self.config.sample_evidence_path = ""
        self._invalidate_sample_caches()
        for attr in ("_annotated_tss_cache", "_annotated_tes_cache"):
            if hasattr(self, attr):
                delattr(self, attr)
        self._close_jc_context()
        # An annotation can change between runs of one runner (a new --gtf, or
        # the same path rewritten). The endpoint channels score reads against
        # its TSS/TES, so a stale reader would silently apply the old
        # coordinates: reopen when the file's fingerprint moved.
        gtf_fp = input_fingerprint(getattr(self.config, "gtf_path", ""))
        if gtf_fp != getattr(self, "_gtf_fingerprint", None):
            if self._gtf_reader is not None:
                try:
                    self._gtf_reader.close()
                except Exception:  # noqa: BLE001 - closing a stale reader is best effort
                    logger.debug("stale GTF reader close failed", exc_info=True)
                self._gtf_reader = None
            self._gtf_fingerprint = None
            if getattr(self.config, "gtf_path", ""):
                # _open_gtf_reader records the fingerprint only after a
                # successful parse, so a failed reopen is retried next time
                # instead of being remembered as current
                self._gtf_reader = self._open_gtf_reader()
            else:
                self._gtf_fingerprint = gtf_fp
        elif getattr(self.config, "gtf_path", "") and self._gtf_reader is None:
            self._gtf_reader = self._open_gtf_reader()      # a previous open failed

    def _input_fingerprints(self) -> dict:
        """Fresh fingerprints of the files the evidence is derived from."""
        return {
            "bam": input_fingerprint(self.config.bam_path),
            "signal": input_fingerprint(self.config.signal_path),
        }

    def _evidence_identity(self):
        """The input fingerprint the evidence for THIS run must carry."""
        cached = getattr(self, "_evidence_identity_cache", None)
        if cached is None:
            cached = self._evidence_identity_cache = {
                "bam": input_fingerprint(self.config.bam_path),
                "signal": input_fingerprint(self.config.signal_path),
                "pore": str(self.config.krill_pore),
                # the polyA caller runs on GPU or CPU; the calls are the
                # evidence, so the backend belongs to their identity. This is
                # the REQUESTED backend: krill may fall back to CPU when a GPU
                # is unavailable, which this cannot see.
                "gpu": str(bool(self.config.use_gpu)),
            }
        return cached

    def _run_token(self, new_attempt: bool = False) -> str:
        """Token naming this runner's disposable files. A pid alone is not
        exclusive: one process can hold two runners, and a RETRY within one
        runner must not adopt the artifact of its own earlier failed attempt --
        so each publication attempt takes a fresh token."""
        token = getattr(self, "_run_token_value", None)
        if token is None or new_attempt:
            import uuid

            token = self._run_token_value = f"{os.getpid()}_{uuid.uuid4().hex[:8]}"
        return token

    def _evidence_path(self, identity=None, run_local=False, new_attempt=False):
        """Work-dir evidence file, ADDRESSED BY ITS INPUT IDENTITY.

        Two runs may share a work directory (the benchmark harness does). With
        one fixed filename they overwrite each other, and a consumer either
        reads endpoints derived from another sample or -- once that is detected
        -- silently loses its own. Deriving the name from the identity makes the
        file immutable for a given input set: a different input is a different
        file, so there is nothing to race over.

        ``run_local=True`` adds this process's id: evidence whose polyA pass
        FAILED is valid only for the run that produced it, and must never
        replace another run's successful evidence for the same inputs.

        Once the parent prepass has chosen a file it is recorded in
        ``config.sample_evidence_path``; workers (which receive the config) then
        read that exact file instead of fingerprinting the inputs again."""
        import hashlib
        import os

        recorded = getattr(self.config, "sample_evidence_path", "")
        if recorded and identity is None and not run_local:
            return Path(recorded)
        ident = identity if identity is not None else self._evidence_identity()
        key = hashlib.blake2b(
            "|".join(f"{k}={ident[k]}" for k in sorted(ident)).encode("utf-8"),
            digest_size=8).hexdigest()
        suffix = f"_run{self._run_token(new_attempt)}" if run_local else ""
        return Path(self.config.work_dir) / f"_sample_evidence_{key}{suffix}.json"

    def _prepare_sample_evidence(self) -> None:
        """Strand-aware 3'/5' ends for every primary read (always) plus krill
        polyA calls (only with ``--signal``), written once to
        ``<work_dir>/_sample_evidence.json``. Idempotent."""
        import hashlib
        import json
        import os
        import shutil
        import tempfile

        identity = self._evidence_identity()
        out = self._evidence_path(identity)
        self.config.sample_evidence_path = str(out)
        # A run WITH --signal whose krill pass failed must not leave a file that
        # looks like a valid "no polyA in this sample" answer: it is recorded as
        # unavailable, which is a cache MISS for a later run that can retry.
        _wanted_polya_status = "ok" if self.config.signal_path else "no_signal"
        # "partial"/"unavailable" artifacts are run-local by construction, so a
        # prepass only ever accepts a complete one.
        if out.exists():
            if self._read_sample_evidence(out, identity, _wanted_polya_status) is not None:
                return
            logger.info("Sample evidence in work dir is unusable (version/inputs/shape); rebuilding")
        Path(self.config.work_dir).mkdir(parents=True, exist_ok=True)
        # The evidence depends only on the sample (BAM + signal + pore), not on
        # the annotation or any flag, so a shared cache (PYFIN_SAMPLE_EVIDENCE_CACHE)
        # lets every run of the same sample reuse one krill pass.
        shared_dir = os.environ.get("PYFIN_SAMPLE_EVIDENCE_CACHE")
        shared = None
        if shared_dir:
            key = hashlib.blake2b("|".join(
                [f"{k}={identity[k]}" for k in sorted(identity)]
                + ["sample-evidence-v3"],
            ).encode("utf-8"), digest_size=8).hexdigest()
            shared = Path(shared_dir) / f"sample_evidence_{key}.json"
            if shared.exists():
                # Copy FIRST, then validate the copy: validating the shared path
                # and re-reading it afterwards would accept one payload and
                # install another if a concurrent writer replaced it between the
                # two opens.
                fd, ctmp = tempfile.mkstemp(prefix=out.name + ".", suffix=".tmp",
                                            dir=str(out.parent))
                os.close(fd)
                try:
                    # The entry can be evicted or become unreadable between the
                    # existence check and the copy: a cache miss, not a failure.
                    shutil.copyfile(shared, ctmp)
                    ok = self._read_sample_evidence(
                        ctmp, identity, _wanted_polya_status) is not None
                except OSError as exc:  # noqa: BLE001 - the cache is an optimisation
                    logger.info("Sample evidence cache unreadable (%s); rebuilding", exc)
                    ok = False
                try:
                    if ok:
                        # create-only, exactly like publication: importing must
                        # not replace an artifact another run may be reading
                        try:
                            os.link(ctmp, out)
                            imported = True
                        except FileExistsError:
                            imported = self._read_sample_evidence(
                                out, identity, _wanted_polya_status) is not None
                            if not imported:
                                logger.warning(
                                    "Sample evidence at %s is unusable and is not "
                                    "replaced; using a run-local copy", out)
                                out = self._evidence_path(
                                    identity, run_local=True, new_attempt=True)
                                self.config.sample_evidence_path = str(out)
                                os.link(ctmp, out)      # attempt-unique: cannot exist
                                self._run_local_evidence = getattr(
                                    self, "_run_local_evidence", set())
                                self._run_local_evidence.add(str(out))
                                imported = True
                        if imported:
                            logger.info("Sample evidence reused from cache: %s", shared)
                            self._invalidate_sample_caches()
                            return
                    logger.info("Sample evidence cache entry unusable; rebuilding: %s", shared)
                finally:
                    try:
                        if os.path.exists(ctmp):
                            os.unlink(ctmp)
                    except OSError:
                        logger.debug("could not remove evidence temp %s", ctmp, exc_info=True)
        polya: Dict[str, list] = {}
        ends: Dict[str, list] = {}
        # The strand-aware ends are BAM evidence and are ALWAYS collected: the
        # 5' start channel needs nothing else, so it must keep working without
        # a signal file. Read sequences are gathered only when krill can use
        # them; without --signal every polyA call is MISSING (abstention).
        want_polya = bool(self.config.signal_path)
        read_seqs: Dict[str, str] = {}
        with pysam.AlignmentFile(self.config.bam_path, "rb") as bam:
            for r in bam.fetch():
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                rid = r.query_name
                if rid is None or rid in ends:
                    continue
                if r.reference_start is None or r.reference_end is None:
                    continue
                strand = "-" if r.is_reverse else "+"
                three = int(r.reference_end) if strand == "+" else int(r.reference_start)
                five = int(r.reference_start) if strand == "+" else int(r.reference_end)
                ends[rid] = [r.reference_name, strand, three, five]
                if want_polya and r.query_sequence:
                    read_seqs[rid] = r.query_sequence
        polya_status = _wanted_polya_status
        if want_polya:  # noqa: SIM102 - status decided inside
            from fin.scoring.polya import compute_polya

            pm = compute_polya(read_seqs, self.config.signal_path,
                               pore=self.config.krill_pore, use_gpu=self.config.use_gpu)
            # Validate at the SOURCE: a call is (length, qc) with a finite
            # length. A record that is neither (NaN, a wrong shape) is not a
            # failed call, it is not a call -- it is dropped here so it can
            # never be published and later invalidate the whole artifact,
            # read ends included. Dropping a call is abstention, the same as
            # a read the caller never scored.
            polya, n_bad = {}, 0
            for rid, v in pm.items():
                call = _usable_polya_call(v)
                if call is None:
                    n_bad += 1
                    continue
                polya[rid] = call
            if n_bad:
                # a partial result is this run's evidence, never the answer
                # another run will inherit for these inputs
                polya_status = "partial"
                logger.warning(
                    "polyA caller returned %d unusable call(s); those reads "
                    "abstain on the 3' channel and this evidence stays local", n_bad)
            if read_seqs and not polya:
                # compute_polya returns {} when the krill backend or the batch
                # failed. That is a broken run, not "this sample has no polyA":
                # record it so the next run retries instead of inheriting a
                # permanently abstaining 3' channel.
                polya_status = "unavailable"
                logger.warning(
                    "polyA calling produced no calls for %d reads; the 3' channel "
                    "will abstain and this evidence is not cached as final",
                    len(read_seqs))
        # The identity was taken before collection. If an input was replaced
        # while this run was reading it, the evidence would describe a mixture
        # of two files under one name: refuse to publish it. Checked AFTER the
        # polyA pass so the signal file is bracketed too.
        after = self._input_fingerprints()
        if after != {k: identity[k] for k in after}:
            raise RuntimeError(
                "input files changed while the sample evidence was being "
                "collected (BAM/signal replaced mid-run); rerun on stable inputs")
        # Serialise ONCE and publish the same bytes everywhere: re-reading the
        # work file to fill the shared cache would publish whatever another run
        # sharing this work dir had replaced it with.
        payload = json.dumps({
            "polya": polya, "ends": ends, "version": "sample-evidence-v3",
            "identity": identity, "polya_status": polya_status,
            # completeness manifest: a truncated-but-parseable artifact (records
            # dropped, or both maps emptied) is otherwise indistinguishable from
            # a sample that genuinely has no ends and no calls
            "n_ends": len(ends), "n_polya": len(polya)})
        if polya_status in ("unavailable", "partial"):
            # This run could not produce a complete polyA result. Its evidence
            # is still usable (the read ends are there, and the calls it did
            # get) but it is NOT the answer for these inputs, so it lives under
            # a name of its own and is never shared.
            out = self._evidence_path(identity, run_local=True, new_attempt=True)
            self.config.sample_evidence_path = str(out)
            self._publish_evidence(out, payload, run_local=True)
        else:
            # The identity-addressed artifact is IMMUTABLE: it is created once
            # and never replaced. Whoever writes it first wins; a later run that
            # finds it there uses it (its inputs are the same by construction)
            # instead of overwriting a file other processes may be reading.
            existing = self._publish_evidence(out, payload)
            if existing is not None:
                # Another run published first. Its payload describes the same
                # inputs, so it is the artifact -- unless it does not parse, in
                # which case this run uses its own copy rather than replacing a
                # file other processes may be reading right now.
                if self._read_sample_evidence(existing, identity, polya_status) is not None:
                    logger.info("Sample evidence already published by another run: %s", out)
                else:
                    logger.warning(
                        "Sample evidence at %s is unusable and is not replaced "
                        "(another run may be reading it); using a run-local copy. "
                        "Delete it to let a later run republish.", existing)
                    out = self._evidence_path(identity, run_local=True, new_attempt=True)
                    self.config.sample_evidence_path = str(out)
                    self._publish_evidence(out, payload, run_local=True)
            if shared is not None:
                try:
                    shared.parent.mkdir(parents=True, exist_ok=True)
                    self._publish_evidence(shared, payload)
                except OSError as exc:  # noqa: BLE001 - the cache is an optimisation
                    logger.warning("sample evidence cache write failed: %s", exc)
        self._invalidate_sample_caches()
        logger.info("Sample evidence written: %s (%d polyA calls, %d ends)",
                    out, len(polya), len(ends))

    def _publish_evidence(self, dest, payload: str, run_local: bool = False):
        """Create ``dest`` with ``payload``, WITHOUT replacing an existing file.

        The evidence of one input set is immutable: several runs may derive it
        concurrently and their consumers may be reading it, so publication is a
        create (link from a writer-unique temp), never a replace. Losing that
        race is the normal case, not an error -- the winner's payload describes
        the same inputs. Returns the existing file's path when the race was lost,
        else None. ``run_local`` files belong to this run alone and are
        registered for removal at cleanup."""
        import tempfile

        dest = Path(dest)
        fd, tmp = tempfile.mkstemp(prefix=dest.name + ".", suffix=".tmp",
                                   dir=str(dest.parent))
        try:
            with os.fdopen(fd, "w") as fh:
                fh.write(payload)
            try:
                os.link(tmp, dest)
            except FileExistsError:
                if run_local:
                    # the name carries this attempt's token, so it cannot
                    # legitimately exist: refuse rather than adopt it
                    raise RuntimeError(
                        f"run-local sample evidence {dest} already exists") from None
                return dest
            if run_local:
                # register BEFORE the fallible temp cleanup: the file exists
                # now, and cleanup must know about it even if unlinking the
                # temp raises
                self._run_local_evidence = getattr(self, "_run_local_evidence", set())
                self._run_local_evidence.add(str(dest))
        finally:
            # the temp is ours either way: a failed write must not leave it
            try:
                if os.path.exists(tmp):
                    os.unlink(tmp)
            except OSError:
                logger.debug("could not remove evidence temp %s", tmp, exc_info=True)
        return None

    def _nuisance_evidence(self):
        """The sample's per-read endpoint distances, derived once per process.
        Both fit sites (the discovery contradiction hook and the per-interval
        assignment context) consume this; only the EM fit is per interval."""
        cached = getattr(self, "_nuisance_evidence_cache", "unset")
        if cached != "unset":
            return cached
        self._nuisance_evidence_cache = None
        ev = self._sample_evidence()
        if ev is not None:
            from fin.analysis.sample_nuisance import sample_nuisance_evidence

            polya_map, read_ends, read_starts = ev
            want_tss = getattr(self.config, "assignment_channels", "none") == "ends"
            self._nuisance_evidence_cache = sample_nuisance_evidence(
                polya_map, read_ends, claimed_tes=self._annotated_tes(),
                read_starts=read_starts if want_tss else None,
                claimed_tss=self._annotated_tss() if want_tss else None,
                accept_adapter=bool(getattr(
                    self.config, "polya_end_accept_adapter", False)),
                qc_conditional=bool(getattr(
                    self.config, "endpoint_qc_conditional", False)))
        return self._nuisance_evidence_cache

    def _read_extent_contradicted(self):
        """Callable(gtf_candidate, read_ids, alternatives) -> {"tes"/"tss":
        [contradicting read ids]} under the sample's endpoint kernels (empty =
        consistent); only those reads define the re-estimated end. Applied to single-exon annotation
        matched by overlap and to multi-exon annotated keys (their wobble-
        equivalent reads). None when the endpoint channels are off."""
        cached = getattr(self, "_read_extent_hook_cache", "unset")
        if cached != "unset":
            return cached                       # sample-level: identical per interval
        self._read_extent_hook_cache = None
        if getattr(self.config, "assignment_channels", "none") == "none":
            return None
        ev = self._sample_evidence()
        evidence = self._nuisance_evidence()
        if ev is None or evidence is None:
            return None
        from fin.analysis.molecule_evidence import polya_end_evidential
        from fin.analysis.sample_nuisance import fit_from_evidence

        polya_map, read_ends, read_starts = ev
        want_tss = getattr(self.config, "assignment_channels", "none") == "ends"
        # This hook stays on the POOLED kernel and the name-based admission
        # predicate even under --endpoint-qc-conditional. Its test is an
        # ABSOLUTE one -- narrow_posterior(d) < 0.5 -- and that threshold is not
        # invariant to the fitted eta: a class fitted as LESS reliable has a
        # lower narrow posterior at every distance, so it would produce MORE
        # contradictions and be protected by FEWER annotated alternatives, which
        # inverts the intended meaning. The class-conditional widening therefore
        # applies to assignment, whose question is a within-read contrast where a
        # low-confidence class correctly moves the answer less.
        # NOTE: no interval exclusion here, unlike `_endpoint_context`. This is
        # a sample-level hook consulted during discovery, before intervals are
        # judged, so every eligible read contributes to the kernel it is then
        # tested against (a read with no end coordinate never reaches the fit; a
        # read on a locus with no annotation is skipped by the annotated-distance
        # branch but still participates in the split-half fallback).
        # Pre-existing; the leave-one-interval-out contract covers assignment
        # only.
        nuisance = fit_from_evidence(evidence)
        endpoint = nuisance.endpoint if nuisance.polya_available else None
        start = nuisance.start if want_tss else None
        min_tss_share = float(getattr(self.config, "tss_contradiction_min_share", 0.0))
        if start is not None and getattr(self.config, "tss_contradiction_auto_gate", False):
            # Does THIS run's 5' end resolve annotated starts at all? Where it
            # does, an outlier start cluster is evidence; where it does not, it
            # is dispersion. See config for the measured bimodality.
            from fin.analysis.sample_nuisance import annotated_tss_tracking

            anno_tss = self._annotated_tss()
            frac, n_loci = annotated_tss_tracking(
                getattr(evidence, "starts_by_locus", None) or {},
                anno_tss or {},
                track_bp=int(getattr(self.config, "tss_contradiction_track_bp", 25)),
                min_mode_share=float(getattr(
                    self.config, "tss_contradiction_mode_share", 0.5)),
            )
            need = float(getattr(self.config, "tss_contradiction_track_frac", 0.89))
            if n_loci == 0:
                logger.info(
                    "5' extent contradiction: no annotated TSS could be scored; "
                    "leaving it ENABLED (unknown is not evidence)")
            elif frac < need:
                # NOT switched off: switched to the per-LOCUS form of the same
                # question. A sample whose 5' ends are imprecise still has loci
                # where the reads cluster on one start, and those are exactly
                # the loci a truncated annotation needs. Turning the mechanism
                # off wholesale costs 3.60 F1 on a 5'-truncated human
                # annotation; demanding a dominant local start instead costs
                # 0.28 and keeps most of the gain.
                gate_share = float(getattr(
                    self.config, "tss_contradiction_gate_share", 0.5))
                logger.info(
                    "5' extent contradiction RESTRICTED to dominant local "
                    "starts (share >= %.2f): this run resolves a 5' end at only "
                    "%.3f of %d annotated loci (need %.2f)",
                    gate_share, frac, n_loci, need)
                min_tss_share = max(min_tss_share, gate_share)
            else:
                logger.info(
                    "5' extent contradiction enabled: this run resolves a 5' end "
                    "at %.3f of %d annotated loci (need %.2f)", frac, n_loci, need)
        min_start_reads = int(getattr(self.config, "tss_contradiction_min_reads", 3))
        # Magnitude floors: a read end this close to the annotated end never
        # contradicts it, whatever the kernel says. See config for the
        # measurement that motivates an asymmetric 5'/3' treatment.
        min_tss_bp = int(getattr(self.config, "tss_contradiction_min_bp", 0))
        min_tes_bp = int(getattr(self.config, "tes_contradiction_min_bp", 0))
        accept_adapter = bool(getattr(self.config, "polya_end_accept_adapter", False))

        def _tes(c):
            return c.end if getattr(c, "strand", "+") == "+" else c.start

        def _tss(c):
            return c.start if getattr(c, "strand", "+") == "+" else c.end

        def contradicted(gc, read_ids, alternatives=()) -> Dict[str, List[str]]:
            """{"tes": [read ids], "tss": [read ids]} -- per annotated end of
            ``gc``, the reads that contradict it and that no alternative
            annotated model explains (a read landing within the narrow kernel
            of another annotated TES/TSS came from that model). Only those
            reads may define the re-estimated end. Empty (falsy) = consistent."""
            strand = getattr(gc, "strand", "+")
            t_c, s_c = _tes(gc), _tss(gc)
            alt_tes = [_tes(a) for a in alternatives]
            alt_tss = [_tss(a) for a in alternatives]
            out: Dict[str, List[str]] = {}
            for rid in sorted(read_ids):
                if endpoint is not None:
                    loc = read_ends.get(rid)
                    if (polya_end_evidential(polya_map.get(rid), accept_adapter)
                            and loc is not None):
                        e_r = loc[2]
                        d = (e_r - t_c) if strand == "+" else (t_c - e_r)
                        if (d > min_tes_bp and endpoint.narrow_posterior(d) < 0.5
                                and not any(endpoint.narrow_posterior(abs(e_r - t)) >= 0.5
                                            for t in alt_tes)):
                            out.setdefault("tes", []).append(rid)   # polyA read past the claimed end
                if start is not None and getattr(start, "available", False):
                    loc5 = read_starts.get(rid)
                    if loc5 is not None:
                        s_r = loc5[2]
                        u = (s_c - s_r) if strand == "+" else (s_r - s_c)
                        if (u > min_tss_bp and start.narrow_posterior(u) < 0.5
                                and not any(start.narrow_posterior(abs(s_r - t)) >= 0.5
                                            for t in alt_tss)
                                # a start is evidence about a TRANSCRIPT's start
                                # only when other molecules end there too
                                and evidence.start_support(loc5[0], loc5[1], s_r)
                                >= min_start_reads
                                # ... and be the DOMINANT start of its locus, not
                                # one shoulder of a 5'-degradation smear
                                and (min_tss_share <= 0.0
                                     or evidence.start_share(loc5[0], loc5[1], s_r)
                                     >= min_tss_share)):
                            out.setdefault("tss", []).append(rid)   # read begins upstream
            return out

        self._read_extent_hook_cache = contradicted
        return contradicted

    def _endpoint_context(self, candidate_set):
        """(polya states by read, endpoint model, read 3' ends) for the assignment's
        3'/polyA channel, or None when the channel is off / unavailable. The
        endpoint model is fitted leave-one-interval-out against the input
        annotation's claimed TES coordinates (see sample_nuisance)."""
        if getattr(self.config, "assignment_channels", "none") == "none":
            return None
        ev = self._sample_evidence()
        evidence = self._nuisance_evidence()
        if ev is None or evidence is None:
            return None
        from fin.analysis.sample_nuisance import fit_from_evidence

        polya_map, read_ends, read_starts = ev
        interval_reads = set(getattr(candidate_set, "read_ids", set()) or ())
        channels = getattr(self.config, "assignment_channels", "none")
        nuisance = fit_from_evidence(evidence, exclude_reads=interval_reads)
        qc_conditional = bool(getattr(self.config, "endpoint_qc_conditional", False))
        ctx = {"polya": polya_map, "ends": read_ends,
               "accept_adapter": bool(getattr(
                   self.config, "polya_end_accept_adapter", False)),
               "qc_conditional": qc_conditional,
               "endpoint_by_qc": nuisance.endpoint_by_qc if qc_conditional else None,
               "endpoint": nuisance.endpoint if nuisance.polya_available else None,
               "starts": read_starts if channels == "ends" else None,
               "start": nuisance.start if channels == "ends" else None}
        if (ctx["endpoint"] is None and not ctx["endpoint_by_qc"]
                and ctx["start"] is None):
            return None
        return ctx

    def _annotated_tss(self):
        """(chrom, strand) -> TSS positions of the INPUT annotation (None without one)."""
        cached = getattr(self, "_annotated_tss_cache", "unset")
        if cached != "unset":
            return cached
        self._annotated_tss_cache = None
        reader = getattr(self, "_gtf_reader", None)
        if reader is None:
            return None
        out: Dict[Tuple[str, str], list] = {}
        try:
            for tx in (reader.transcripts.values() if hasattr(reader, "transcripts") else []):
                strand = getattr(tx, "strand", "+")
                tss = tx.start if strand == "+" else tx.end
                out.setdefault((tx.chrom, strand), []).append(int(tss))
        except Exception:  # noqa: BLE001
            return None
        self._annotated_tss_cache = out or None
        return self._annotated_tss_cache

    def _annotated_tes(self):
        """(chrom, strand) -> TES positions of the INPUT annotation (None without one)."""
        cached = getattr(self, "_annotated_tes_cache", "unset")
        if cached != "unset":
            return cached
        self._annotated_tes_cache = None
        reader = getattr(self, "_gtf_reader", None)
        if reader is None:
            return None
        out: Dict[Tuple[str, str], list] = {}
        try:
            txs = reader.transcripts.values() if hasattr(reader, "transcripts") else []
            for tx in txs:
                strand = getattr(tx, "strand", "+")
                tes = tx.end if strand == "+" else tx.start
                out.setdefault((tx.chrom, strand), []).append(int(tes))
        except Exception:  # noqa: BLE001 - nuisance must degrade, never abort
            return None
        self._annotated_tes_cache = out or None
        return self._annotated_tes_cache

    @staticmethod
    def _read_sample_evidence(path, identity=None, polya_status=None):
        """Parse an evidence file, or None when it is unusable.

        Derived evidence read back from disk is INPUT, not a promise: a
        truncated write, a payload from another sample, or a record of the
        wrong shape must rebuild rather than crash or silently disable a
        channel. ``identity`` / ``polya_status``, when given, must match."""
        import json

        try:
            with Path(path).open() as fh:
                d = json.load(fh)
        except (OSError, ValueError):
            return None
        if not isinstance(d, dict):
            return None
        if d.get("version") != "sample-evidence-v3":
            return None
        if identity is not None and d.get("identity") != identity:
            return None
        if polya_status is not None and d.get("polya_status") != polya_status:
            return None
        ends, polya = d.get("ends"), d.get("polya")
        if not isinstance(ends, dict) or not isinstance(polya, dict):
            return None
        # the writer recorded how many records it wrote; anything else is a
        # damaged file, not a sample with fewer reads
        if (d.get("n_ends") != len(ends)) or (d.get("n_polya") != len(polya)):
            return None

        # Every consumer eventually converts these to float (polya_state, the
        # endpoint kernels); a value outside the double range would raise there,
        # far from the file that caused it. Reject it here as an unusable
        # payload -- genomic coordinates and polyA lengths are far below this.
        _NUMBER_MAX = 2.0 ** 63

        def _finite_number(x) -> bool:
            """A JSON number usable as a coordinate/length. ``math.isnan`` on a
            huge int raises OverflowError, so ints are checked as ints."""
            if isinstance(x, bool):
                return False                     # bool is an int subclass
            if isinstance(x, int):
                return -_NUMBER_MAX <= x <= _NUMBER_MAX
            return (isinstance(x, float) and not (math.isnan(x) or math.isinf(x))
                    and -_NUMBER_MAX <= x <= _NUMBER_MAX)

        def _coord(x):
            if not _finite_number(x):
                return None
            if isinstance(x, float) and x != int(x):
                return None                      # a coordinate is a whole number
            return int(x)

        read_ends, read_starts, polya_map = {}, {}, {}
        for rid, v in ends.items():
            # [chrom, strand, three_prime, five_prime]: a record of any other
            # shape is a corrupt file, not a read without a 5' end -- accepting
            # it would silently disable a channel or fabricate a coordinate.
            if (not isinstance(rid, str) or not isinstance(v, list) or len(v) != 4
                    or not isinstance(v[0], str) or v[1] not in ("+", "-")):
                return None
            three, five = _coord(v[2]), _coord(v[3])
            if three is None or five is None:
                return None
            read_ends[rid] = (v[0], v[1], three)
            read_starts[rid] = (v[0], v[1], five)
        for rid, v in polya.items():
            # [length, qc]: length None or a finite number, qc None or a string
            if (not isinstance(rid, str) or not isinstance(v, list) or len(v) != 2
                    or not (v[1] is None or isinstance(v[1], str))):
                return None
            if v[0] is not None and not _finite_number(v[0]):
                return None
            polya_map[rid] = (v[0], v[1])
        return polya_map, read_ends, read_starts

    def _sample_evidence(self):
        """Load this run's sample evidence (None when the prepass did not run or
        the file on disk does not belong to this run's inputs).

        The identity is re-checked HERE, not only at write time: a work
        directory can be shared, and a consumer must never score reads against
        endpoints derived from another sample."""
        cached = getattr(self, "_sample_evidence_cache", "unset")
        if cached != "unset":
            return cached
        recorded = getattr(self.config, "sample_evidence_path", "")
        path = Path(recorded) if recorded else self._evidence_path()
        identity = None if recorded else self._evidence_identity()
        parsed = self._read_sample_evidence(path, identity) if path.exists() else None
        if parsed is None:
            if path.exists():
                # Addressed by identity, so this is corruption, not a foreign
                # payload. Say so: abstaining in silence would look like a
                # sample with no polyA calls and no read ends.
                logger.warning(
                    "Sample evidence at %s is unusable; the endpoint channels "
                    "will abstain for this run (delete it and rerun)", path)
            elif recorded:
                # The prepass prepared this file and it is gone: a worker must
                # not silently score its intervals without endpoint evidence.
                logger.warning(
                    "Prepared sample evidence %s is missing; the endpoint "
                    "channels will abstain for this run", path)
        self._sample_evidence_cache = parsed
        return self._sample_evidence_cache

    def _fetch_read_seqs_and_ends(
        self,
    ) -> Tuple[Dict[str, str], Dict[str, Tuple[int, int]]]:
        """Single BAM pass -> ({rid: query_sequence}, {rid: (ref_start, ref_end)}).

        Primary mapped alignments only; first occurrence per read id wins. The
        sequences feed krill whole-read polyA; the genomic spans feed 5'
        proximity (mirrors ``_annotate_fulllen_frac``'s read-end map).
        """
        read_seqs: Dict[str, str] = {}
        read_ends: Dict[str, Tuple[int, int]] = {}
        with pysam.AlignmentFile(self.config.bam_path, "rb") as bam:
            for r in bam.fetch():
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                rid = r.query_name
                if rid is None or rid in read_ends:
                    continue
                rs, re = r.reference_start, r.reference_end
                if rs is None or re is None:
                    continue
                read_ends[rid] = (int(rs), int(re))
                if r.query_sequence:
                    read_seqs[rid] = r.query_sequence
        return read_seqs, read_ends

    def _apply_polya5p_filter(
        self, aggregated: Dict[str, QuantResult]
    ) -> Dict[str, QuantResult]:
        """Drop candidates failing the polyA + 5'-proximity gate.

        Runs a krill whole-read polyA pass over all reads, then removes any
        gated candidate with fewer than ``min_polya5p_reads`` reads that both
        have a confident polyA tail and map 5'-flush to the candidate. Fusion
        candidates are always exempt; GTF candidates are exempt by default
        (``polya5p_exempt_gtf``), so only novel candidates are gated. No-ops if
        krill returns nothing.
        """
        from fin.scoring.polya import compute_polya

        read_seqs, read_ends = self._fetch_read_seqs_and_ends()
        polya_map = compute_polya(
            read_seqs,
            self.config.signal_path,
            pore=self.config.krill_pore,
            use_gpu=self.config.use_gpu,
        )
        if not polya_map:
            logger.warning(
                "polyA+5' filter: krill returned no polyA estimates; skipping"
            )
            return aggregated

        drop_ids = polya5p_drops(
            aggregated,
            polya_map,
            read_ends,
            window=self.config.polya5p_window_bp,
            min_polya_len=self.config.min_polya_length,
            min_reads=self.config.min_polya5p_reads,
            exempt_gtf=self.config.polya5p_exempt_gtf,
        )
        if not drop_ids:
            return aggregated

        n_gtf = sum(
            1 for cid in drop_ids if aggregated[cid].source == "gtf"
        )
        n_novel = len(drop_ids) - n_gtf
        aggregated = {
            cid: qr for cid, qr in aggregated.items() if cid not in drop_ids
        }
        logger.info(
            "Dropped %d gtf + %d novel transcripts failing polyA+5' (>= %d reads)",
            n_gtf,
            n_novel,
            self.config.min_polya5p_reads,
        )
        return aggregated

    def _adaptive_sigma(self, dist_read_to_tx: np.ndarray) -> float:
        """Data-adaptive EM sigma: median per-read range of the read×tx distance,
        clipped to [em_sigma_min, em_sigma_max]. Falls back to em_sigma.

        With per-event-normalized krill distances the absolute scale is dataset-
        dependent (~0.5-5 nats); a fixed sigma=1.0 often collapses R to one-hot.
        """
        n_reads, n_tx = dist_read_to_tx.shape
        if n_tx >= 2 and n_reads > 0:
            d_max = dist_read_to_tx.max(axis=1)
            d_min = dist_read_to_tx.min(axis=1)
            adaptive = float(np.median(d_max - d_min))
            return float(np.clip(
                adaptive if adaptive > 0 else self.config.em_sigma,
                getattr(self.config, "em_sigma_min", 0.05),
                getattr(self.config, "em_sigma_max", 50.0),
            ))
        return self.config.em_sigma

    def _eff_lengths(self, cand_list) -> Optional[np.ndarray]:
        """Per-candidate spliced effective length in cand_list column order, for
        Salmon-style abundance-feedback length normalization. Returns None unless
        abundance_feedback AND abundance_length_norm are both on. Lengths are
        floored at 1.0 so a degenerate (zero-length) candidate cannot divide-by-0.
        """
        if not (self.config.abundance_feedback and self.config.abundance_length_norm):
            return None
        lengths = np.array(
            [
                max(1.0, float(sum(e - s for s, e in _exons_from_candidate(c))))
                for c in cand_list
            ],
            dtype=float,
        )
        return lengths

    def _quant_m1_em(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> Optional[List[QuantResult]]:
        """quant_mode='m1_em': EM seeded by the M1 mappy AS-gap distance (β=0,
        no signal coherence). Pure-alignment soft assignment over krill-free
        mappy distances."""
        from fin.ablation.mappy_argmax import mappy_multimap_responsibilities

        cand_list = list(candidate_set.candidates)
        read_sequences = getattr(candidate_set, "read_sequences", {}) or {}
        reads_iter = [(rid, read_sequences.get(rid, "")) for rid in read_ids]
        reads_iter = [(rid, seq) for rid, seq in reads_iter if seq]
        R_mm, kept_read_ids = mappy_multimap_responsibilities(reads_iter, cand_list)
        if R_mm.size == 0:
            return []

        read_seqs = {rid: seq for rid, seq in reads_iter}
        n_reads_em = len(kept_read_ids)
        n_cands_em = len(cand_list)
        max_iter_em = (
            self.config.em_max_iter_override
            if self.config.em_max_iter_override is not None
            else self.config.em_max_iter
        )

        m1 = compute_mappy_distance(read_seqs, cand_list, kept_read_ids)
        m2_dummy = np.zeros((n_reads_em, n_cands_em), dtype=np.float32)
        dist_read_to_tx, _, _ = build_em_matrices("m1", m1, m2_dummy)
        dist_read_to_read = np.zeros((n_reads_em, n_reads_em), dtype=np.float32)
        R, hard_assignments, _ = em_with_coherence(
            dist_read_to_tx=dist_read_to_tx,
            dist_read_to_read=dist_read_to_read,
            sigma=self.config.em_sigma,
            beta=0.0,
            max_iter=max_iter_em,
            tol=self.config.em_tol,
            verbose=False,
            use_gpu=self.config.use_gpu,
            abundance_feedback=self.config.abundance_feedback,
            abundance_length_norm=self.config.abundance_length_norm,
            eff_lengths=self._eff_lengths(cand_list),
        )

        if self.config.krill_tiebreak:
            R = krill_tiebreak(
                R=R, read_ids=kept_read_ids, read_seqs=read_seqs,
                candidates=cand_list, signal_path=self.config.signal_path,
                pore=self.config.krill_pore,
                ambig_threshold=self.config.tiebreak_ambig_threshold,
                use_gpu=self.config.use_gpu,
            )
        # Same explicit tie-break policy as the production m2_em path, applied
        # UNCONDITIONALLY: leaving it inside the krill branch would mean the
        # default (krill off) still resolved exact ties by column index, i.e. by
        # emission order, which is the behaviour this rule exists to remove.
        from fin.pipeline.assignment import hard_assign

        hard_assignments, _, _ = hard_assign(R, cand_list)

        quant_results = quantify_transcripts(
            R, hard_assignments, cand_list, kept_read_ids
        )
        for j, qr in enumerate(quant_results):
            qr.max_R = float(R[:, j].max()) if R.shape[0] > 0 else 0.0
        logger.info(
            "m1_em interval %s: %d reads -> %d candidates",
            interval.region_string, len(read_ids), len(cand_list),
        )
        return quant_results

    def _tie_nll(self, kept_read_ids, read_seqs, cand_list, aligners, raw):
        """Per-read best-AS tie junction-NLL. Delegates to the assignment layer
        (fin.pipeline.assignment.tie_nll); kept as a thin method so tests can patch
        this mock seam. Returns (nlls_by_read, ties_by_read, n_ties, n_refined,
        cover_by_read)."""
        from fin.pipeline.assignment import tie_nll
        return tie_nll(
            self.config, kept_read_ids, read_seqs, cand_list, aligners, raw)

    def _observed_junctions(self, interval: GenomicInterval):
        """Strand-keyed {strand: Counter{(donor, acceptor): n_reads}} of intron
        junctions directly observed in the interval's primary-read CIGARs.

        Delegates to ``evidence.compute_observed_junctions`` (the Evidence layer) with a
        single-entry per-interval memo: this map is consumed by up to three gates within one
        interval (cluster-recheck GTF guard, Lever-2 per-junction support gate,
        junction-dominance / guided support gate), which used to each re-open the BAM;
        computing it once here is byte-identical. Returns None on empty/unreadable BAM so every
        support gate self-disables (fail-open). See fin/pipeline/evidence.py for the caveat.
        """
        from fin.pipeline.evidence import compute_observed_junctions

        # Memo key = region only: the map fetches chrom:start-end and buckets by READ
        # strand, so it is independent of the interval's own strand (two same-region
        # intervals share it). region_string is present on every interval incl. test stubs.
        key = interval.region_string
        memo = getattr(self, "_obs_junc_memo", None)
        if memo is not None and memo[0] == key:
            return memo[1]
        value = compute_observed_junctions(self.config.bam_path, interval)
        self._obs_junc_memo = (key, value)
        return value

    def _quant_m2_em(
        self,
        candidate_set: CandidateSet,
        read_ids: List[str],
        interval: GenomicInterval,
    ) -> Tuple[List[QuantResult], Optional[ResponsibilityLedger]]:
        """quant_mode='m2_em' (production default): pure tie-break junction-NLL EM.

        M1/AS is used ONLY as a hard selector — it picks each read's best-AS tie
        set and masks mappability — and the per-event junction NLL is the SOLE
        graded distance over that tie set (m2_resolve_tie semantics, NOT the dense
        read×candidate matrix, which collapses). Production uses no read×read
        coherence term. All signal scoring is in-memory krill; there is no f5c
        CLI path.
        """
        endpoint_ctx = self._endpoint_context(candidate_set)
        qo = self._assigner.assign(
            candidate_set, read_ids,
            tie_nll_fn=self._tie_nll, eff_lengths_fn=self._eff_lengths,
            endpoint_context=endpoint_ctx,
        )
        if qo is None:
            return [], None
        # Post-EM folded-TSS recovery (opt-in). A `supported` shadow is added
        # as a candidate and the WHOLE assignment/EM is rerun once, so the
        # second pass alone decides read ownership and abundance; the ledger
        # below is built only from that final pass. Never recursive.
        recovery_mode = getattr(self.config, "folded_tss_recovery_mode", "off")
        if recovery_mode != "off" and getattr(candidate_set, "shadows", None):
            from fin.analysis.tss_shadow_recovery import (
                recover_folded_tss_candidates, write_recovery_shard,
            )

            genome_seq = (
                self._genome_fasta[interval.chrom]
                if self._genome_fasta and interval.chrom in self._genome_fasta
                else ""
            )
            rec = recover_folded_tss_candidates(
                candidate_set, qo, interval, genome_seq,
                indels_by_read=getattr(candidate_set, "indels_by_read", None),
            )
            write_recovery_shard(self.config.work_dir, interval, rec.rows)
            if recovery_mode == "require" and rec.support_merges:
                for cid, reads in rec.support_merges.items():
                    ex = candidate_set.get_candidate(cid)
                    if ex is not None:
                        ex.supporting_read_ids.update(reads)
                        ex.raw_supporting_read_ids.update(reads)
            if rec.candidates and recovery_mode == "require":
                candidate_set.candidates.extend(rec.candidates)
                if candidate_set.clusters is not None:
                    for cand in rec.candidates:
                        parent_id = rec.parent_of.get(cand.candidate_id)
                        for ids in candidate_set.clusters:
                            if parent_id in ids:
                                ids.append(cand.candidate_id)
                                break
                        else:
                            candidate_set.clusters.append([cand.candidate_id])
                logger.info(
                    "m2_em interval %s: recovered %d folded 5'-TSS shadow(s); "
                    "rerunning assignment",
                    interval.region_string, len(rec.candidates),
                )
                qo2 = self._assigner.assign(
                    candidate_set, read_ids,
                    tie_nll_fn=self._tie_nll, eff_lengths_fn=self._eff_lengths,
                    start_likelihood_context=rec.start_models,
                    endpoint_context=endpoint_ctx,
                )
                if qo2 is not None:
                    qo = qo2
        # Unified existence audit: interval rows only (the FDR family is the
        # whole SAMPLE, so it is applied at finalization, never here).
        if getattr(self.config, "candidate_existence_mode", "off") != "off":
            from fin.analysis.candidate_existence import (
                build_conflict_components, evaluate_candidate_existence,
                write_existence_shard,
            )

            comps = build_conflict_components(qo, candidate_set)
            if comps:
                molecules, nuisance = None, None
                ev = self._sample_evidence()
                if ev is not None:
                    from fin.analysis.molecule_evidence import build_molecule_evidence
                    from fin.analysis.sample_nuisance import fit_sample_nuisance

                    polya_map, read_ends, _read_starts = ev
                    molecules = build_molecule_evidence(qo, candidate_set, polya_map=polya_map)
                    # Leave-one-interval-out: the interval under test contributes
                    # nothing to the nuisance it is judged against.
                    nuisance = fit_sample_nuisance(
                        polya_map, read_ends, exclude_reads=set(qo.kept_read_ids),
                        claimed_tes=self._annotated_tes())
                write_existence_shard(
                    self.config.work_dir, interval,
                    evaluate_candidate_existence(
                        qo, candidate_set, comps, interval,
                        sample_id=str(getattr(self.config, "bam_path", "")),
                        molecules=molecules, nuisance=nuisance,
                        bootstrap_replicates=int(getattr(
                            self.config, "candidate_existence_bootstrap", 4095)),
                    ),
                )
        ledger = None
        if getattr(self.config, "post_selection_refit_effective", False):
            ledger = build_responsibility_ledger(
                qo.R,
                qo.kept_read_ids,
                [candidate.candidate_id for candidate in qo.cand_list],
                input_read_ids=read_ids,
                interval_key=(
                    f"{interval.chrom}:{interval.start}-{interval.end}:"
                    f"{interval.strand}"
                ),
            )
        quant_results, outcomes = select_m2_interval(
            self.config, qo, candidate_set, interval, self._observed_junctions,
        )
        if getattr(self, "_audit_interval_rows", None) is not None:
            key = interval.region_string
            self._audit_interval_rows += audit_rows_for_candidates(
                key, "entered_em", qo.cand_list)
            self._audit_interval_rows += audit_rows_for_outcomes(
                key, "interval_cascade", outcomes)
        if ledger is not None:
            annotate_selection_metadata(
                ledger,
                candidates=qo.cand_list,
                read_ids=qo.kept_read_ids,
                hard_assignments=qo.hard_assignments,
                surviving_results=quant_results,
                outcomes=outcomes,
                mono_resolve_applied=bool(getattr(
                    self.config, "mono_resolve_post_em", False
                )),
            )
        logger.info(
            "m2_em interval %s: %d reads -> %d candidates "
            "(ties=%d refined=%d)",
            interval.region_string, len(read_ids), len(qo.cand_list),
            qo.n_ties, qo.n_refined,
        )
        if logger.isEnabledFor(logging.DEBUG):
            import resource

            # ru_maxrss is this PROCESS's high-water mark (kB on Linux), so it
            # only ever rises: the interval that moved it is the one to look at.
            logger.debug(
                "m2_em interval %s: process peak RSS %.2f GB after this interval",
                interval.region_string,
                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0,
            )
        return quant_results, ledger

    def _junction_correction_context(self):
        """Lazily build the per-WORKER krill aligner + blow5 reader for step 5a.

        Each parallel worker constructs exactly one PipelineRunner, so this is
        created once per worker and reused across intervals -- never per contest.
        Returns (aligner, krill_module, signal_reader) or None when unavailable.
        """
        if getattr(self, "_jc_ctx", "unset") != "unset":
            return self._jc_ctx
        self._jc_ctx = None
        try:
            import krill
            import pyslow5

            from fin.scoring.krill_aligner import make_krill_aligner
        except Exception as exc:  # noqa: BLE001 - optional signal stack
            logger.warning("junction correction: signal stack unavailable: %s", exc)
            return None
        aligner, _eff_gpu = make_krill_aligner(
            krill, self.config.krill_pore, self.config.use_gpu,
            hmm_confidence=False, num_thread=0,
        )
        if aligner is None:
            logger.warning("junction correction: krill aligner init failed")
            return None
        s5 = pyslow5.Open(self.config.signal_path, "r")
        self._jc_slow5 = s5          # owned here; closed by reset/cleanup
        cache: dict = {}

        def reader(read_id):
            if read_id in cache:
                return cache[read_id]
            try:
                rec = s5.get_read(read_id, pA=False)
            except Exception:  # noqa: BLE001 - a missing read must abstain
                rec = None
            out = None if rec is None else dict(
                signal=rec["signal"].astype("float32"),
                digitisation=float(rec["digitisation"]), offset=float(rec["offset"]),
                range=float(rec["range"]), sample_rate=float(rec["sampling_rate"]))
            if len(cache) < 4096:
                cache[read_id] = out
            return out

        self._jc_ctx = (aligner, krill, reader)
        return self._jc_ctx

    def _make_junction_correction(self, interval, chrom_seq):
        """Callback handed to discover_candidates; None unless the mode is on."""
        if getattr(self.config, "junction_correction_mode", "off") == "off":
            return None
        ctx = self._junction_correction_context()
        if ctx is None:
            return None
        aligner, krill_mod, reader = ctx

        def _ballot(read_chains, indels_by_read):
            from fin.analysis.junction_correction import (
                correct_read_chains, write_correction_tsv,
            )

            res = correct_read_chains(
                read_chains, indels_by_read, config=self.config, interval=interval,
                genome_seq=chrom_seq, signal_reader=reader, aligner=aligner,
                krill_mod=krill_mod,
            )
            if res.rows:
                shard = (Path(self.config.work_dir) / "_jcorrection_shards"
                         / f"{interval.chrom}_{interval.start}_{interval.end}_"
                           f"{interval.strand or 'u'}.tsv")
                write_correction_tsv(shard, res.rows)
            return res

        return _ballot

    def _junction_correction_apply(self):
        """Instantiation policy object for discovery; None in audit mode (record only)."""
        if getattr(self.config, "junction_correction_mode", "off") != "require":
            return None

        class _Apply:
            min_reads = int(getattr(self.config, "junction_correction_new_chain_min_reads", 2))

        return _Apply()

    def _get_fusion_genome_aligner(self):
        """Lazily build and cache the genome-wide mappy aligner for fusion arms."""
        if not self._fusion_aligner_built:
            from fin.fusion import build_genome_aligner

            self._fusion_aligner_built = True
            if self.config.genome_fasta_path:
                self._fusion_genome_aligner = build_genome_aligner(
                    self.config.genome_fasta_path
                )
        return self._fusion_genome_aligner

    def _augment_with_fusion_candidates(
        self, candidate_set: CandidateSet, interval: GenomicInterval
    ) -> CandidateSet:
        """Run the read-driven fusion sub-pipeline and merge into candidate_set.

        Collects chimeric (soft-clip) reads in the interval, re-aligns their
        soft-clip arms to find partners (Stage F1), infers per-arm splice
        structure (F2), and stitches cross-breakpoint fusion candidates (F3).
        """
        from fin.candidates.canonical import parse_motifs
        from fin.fusion import detect_fusion_candidates
        from fin.io.io_bam import BamReader

        aligner = self._get_fusion_genome_aligner()
        if aligner is None:
            return candidate_set

        with BamReader(self.config.bam_path) as bam:
            read_dicts = bam.get_reads_in_region(interval.region_string)
        if not read_dicts:
            return candidate_set

        fusion_cands = detect_fusion_candidates(
            read_dicts,
            aligner,
            self._genome_fasta or {},
            gtf_reader=self._gtf_reader,
            motif_set=parse_motifs(self.config.canonical_motifs),
            max_internal_gap_bp=self.config.fusion_max_internal_gap_bp,
            max_dist=self.config.fusion_max_dist,
            search_bp=self.config.canonical_search_bp,
            max_chains_per_read=self.config.max_chains_per_read,
            min_support=self.config.fusion_min_support,
        )
        if not fusion_cands:
            return candidate_set

        return merge_fusion_candidates(candidate_set, fusion_cands)

    def _apply_canonical_gate(
        self, candidate_set: CandidateSet, chrom_seq: str
    ) -> None:
        """Stage B canonical-junction gate. Delegates to the selection layer
        (fin.pipeline.selection.canonical_gate_select); kept as a method so tests can
        call this seam directly. Mutates ``candidate_set.candidates`` in place."""
        return canonical_gate_select(self.config, candidate_set, chrom_seq)

    def _apply_junction_dominance_gate(
        self, candidate_set: CandidateSet, interval: GenomicInterval
    ) -> None:
        """Pre-EM junction-dominance gate. Delegates to the selection layer
        (fin.pipeline.selection.junction_dominance_select); kept as a method so tests can
        call this seam directly. Mutates ``candidate_set.candidates`` in place."""
        return junction_dominance_select(
            self.config, candidate_set, interval, self._observed_junctions
        )

    def cleanup(self):
        """Close file handles and drop this run's disposable artifacts."""
        # evidence from a failed polyA pass: valid only for this run (a retry
        # may have moved on to the shared artifact, so every file this run
        # created is tracked, not just the current one)
        for local in sorted(getattr(self, "_run_local_evidence", ())):  # this runner's own
            try:
                Path(local).unlink()
            except OSError:
                pass
        self._run_local_evidence = set()
        self._close_jc_context()
        if self._gtf_reader:
            self._gtf_reader.close()
        if self._signal_reader:
            self._signal_reader.close()
        genome = getattr(self, "_genome_fasta", None)
        if genome is not None and hasattr(genome, "close"):
            genome.close()

    def _load_genome_fasta(self, path: str):
        """Open the genome as a mapping (lazy indexed access by default).

        Memory attribution showed the eager whole-genome dict dominated
        worker RSS (~3.1 GB per spawn worker); the lazy mapping keeps the
        same ``chrom -> sequence`` contract while holding only a bounded
        chromosome cache per process. ``lazy_genome=False`` restores the
        historical eager dict.
        """
        from fin.io.lazy_genome import open_genome

        return open_genome(
            path,
            lazy=bool(getattr(self.config, "lazy_genome", True)),
            cache_chroms=int(getattr(self.config, "genome_cache_chroms", 2)),
        )

    def _open_signal_reader(self):
        """Open the appropriate signal reader."""
        if self.config.signal_format == "pod5":
            from fin.io.io_pod5 import Pod5Reader

            reader = Pod5Reader(self.config.signal_path)
            reader.open()
            return reader
        else:
            from fin.io.io_slow5 import Slow5Reader

            reader = Slow5Reader(self.config.signal_path)
            reader.open()
            return reader
