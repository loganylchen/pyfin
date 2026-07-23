"""Finalize layer: gene-id resolution and output writers.

The last stage of the pipeline, after GLOBAL selection (fin.pipeline.selection.select_global). This
module is intentionally PURE-ish: it resolves gene_ids from the annotation and writes the GTF / TSV /
BEDPE outputs, plus the optional pre-filter diagnostic TSV. No candidate is dropped here — selection
is done by the time these run. Writers are the existing standalone ``fin.io`` functions; this module
just wires them to the aggregated result and the config.

Kept separate from selection so "which candidates survive" (selection.py) and "how survivors are
written out" (finalize.py) are independent concerns — a change to a writer cannot change the result
set, and a change to a filter cannot change output formatting. ``EndpointRefine`` (5'/3'-endpoint
polishing) is the natural future hook here; it is a no-op stub for now.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, Optional

logger = logging.getLogger(__name__)


def write_unfiltered_diagnostic(
    config, aggregated: Dict[str, "object"], output_tsv: Optional[str], gtf_reader
) -> None:
    """Write the pre-filter scoring TSV (diagnostic) when ``write_unfiltered_scores`` is set.

    Runs BEFORE the GLOBAL selection cascade so downstream FN-root-cause analysis can distinguish
    "candidate never reached EM" from "candidate dropped by a filter". Resolves gene_ids on the
    unfiltered snapshot too (so the diagnostic TSV has no empty gene_id column). No-op unless both
    ``write_unfiltered_scores`` and ``output_tsv`` are set. Byte-identical move of the former
    ``_finalize_and_write`` diagnostic block.
    """
    if not (getattr(config, "write_unfiltered_scores", False) and output_tsv):
        return
    from fin.io.io_tsv import write_scoring_tsv

    unfiltered_path = str(Path(output_tsv).with_suffix(".unfiltered.tsv"))
    transcript_lengths_uf = {
        cid: sum(end - start for start, end in qr.exons) if qr.exons else 0
        for cid, qr in aggregated.items()
    }
    # Resolve gene_ids on the unfiltered snapshot too so downstream
    # consumers don't see empty gene_id columns.
    for cid, qr in aggregated.items():
        if qr.source == "gtf" and gtf_reader:
            tx = gtf_reader.get_transcript(cid)
            if tx and not qr.gene_id:
                qr.gene_id = tx.gene_id
        if not qr.gene_id:
            qr.gene_id = qr.candidate_id
    write_scoring_tsv(aggregated, transcript_lengths_uf, unfiltered_path)
    logger.info(
        "Wrote unfiltered scoring TSV (pre-filter): %s (n=%d)",
        unfiltered_path,
        len(aggregated),
    )


def finalize_outputs(
    config,
    aggregated: Dict[str, "object"],
    output_gtf: Optional[str],
    output_tsv: Optional[str],
    gtf_reader,
) -> Dict[str, "object"]:
    """Resolve gene_ids and write the GTF / TSV / BEDPE outputs; return ``aggregated`` unchanged.

    Byte-identical move of the runner's former ``_finalize_and_write`` tail (post-selection): the
    gene-id resolution loop, the three optional writers, and the final "Pipeline complete" log. No
    candidate is dropped here.
    """
    # Resolve gene_ids from GTF annotation
    for cid, qr in aggregated.items():
        if qr.source == "gtf" and gtf_reader:
            tx = gtf_reader.get_transcript(cid)
            if tx:
                qr.gene_id = tx.gene_id
        if not qr.gene_id:
            qr.gene_id = qr.candidate_id

    # Write GTF output
    if output_gtf:
        from fin.io.io_gtf import write_gtf

        write_gtf(aggregated, output_gtf)
        logger.info("Wrote GTF output: %s", output_gtf)

    # US-013: Additional output writers
    if output_tsv:
        from fin.io.io_tsv import write_scoring_tsv

        transcript_lengths = {
            cid: sum(end - start for start, end in qr.exons) if qr.exons else 0
            for cid, qr in aggregated.items()
        }
        write_scoring_tsv(aggregated, transcript_lengths, output_tsv)
        logger.info("Wrote TSV output: %s", output_tsv)

    if config.fusion_enabled and config.output_bedpe:
        from fin.io.io_bedpe import write_fusion_bedpe

        write_fusion_bedpe(aggregated, config.output_bedpe)
        logger.info("Wrote BEDPE output: %s", config.output_bedpe)

    logger.info(
        "Pipeline complete: %d transcripts quantified", len(aggregated)
    )
    return aggregated
