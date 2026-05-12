"""TSV writer for scoring output from transcript quantification."""

from __future__ import annotations

import logging
from typing import Dict

from fin.analysis.quantification import QuantResult, compute_tpm

logger = logging.getLogger(__name__)

_HEADER = (
    "candidate_id\tgene_id\tchrom\tstrand\tstart\tend\tsource\t"
    "abundance\tconfidence\tcoherence_score\tdiscrimination_score\t"
    "combined_score\tnum_reads\ttpm\tmax_R\n"
)


def write_scoring_tsv(
    results: Dict[str, QuantResult],
    transcript_lengths: Dict[str, int],
    path: str,
) -> None:
    """Write quantification scoring results to a TSV file.

    Args:
        results: Dict mapping candidate_id -> QuantResult.
        transcript_lengths: Dict mapping candidate_id -> spliced length in bp.
        path: Output file path.
    """
    tpm_values = compute_tpm(results, transcript_lengths)

    with open(path, "w") as f:
        f.write(_HEADER)
        for cid in sorted(results):
            qr = results[cid]
            tpm = tpm_values.get(cid, 0.0)
            max_r_val = getattr(qr, "max_R", 0.0)
            row = "\t".join([
                qr.candidate_id,
                qr.gene_id,
                qr.chrom,
                qr.strand,
                str(qr.start),
                str(qr.end),
                qr.source,
                f"{qr.abundance:.4f}",
                f"{qr.confidence:.4f}",
                f"{qr.coherence_score:.4f}",
                f"{qr.discrimination_score:.4f}",
                f"{qr.combined_score:.4f}",
                str(qr.num_assigned_reads),
                f"{tpm:.4f}",
                f"{max_r_val:.4f}",
            ]) + "\n"
            f.write(row)

    logger.info("Wrote %d scoring rows to %s", len(results), path)
