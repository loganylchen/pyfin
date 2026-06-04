"""I/O helpers for ablation results: write per-row TSV summary tables."""

from __future__ import annotations

import csv
import logging
from pathlib import Path
from typing import List

from fin.ablation.runner import AblationResult

logger = logging.getLogger(__name__)


def write_ablation_summary(
    results: List[AblationResult],
    output_dir: str,
    filename: str = "ablation_summary.tsv",
) -> str:
    """Write a single TSV with one row per (ablation_row, transcript).

    Columns:
        row_id, label, candidate_id, source, abundance, num_assigned_reads,
        confidence, max_R, coherence_score, discrimination_score, combined_score

    Args:
        results: List of AblationResult, one per ablation row.
        output_dir: Directory to write the TSV.
        filename: Output filename (default: ablation_summary.tsv).

    Returns:
        Absolute path of the written file.
    """
    out_path = Path(output_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "row_id",
        "label",
        "candidate_id",
        "source",
        "abundance",
        "num_assigned_reads",
        "confidence",
        "max_R",
        "coherence_score",
        "discrimination_score",
        "combined_score",
    ]

    total_rows = 0
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for ablation_result in results:
            for cid, qr in sorted(ablation_result.aggregated.items()):
                writer.writerow({
                    "row_id": ablation_result.row_id,
                    "label": ablation_result.label,
                    "candidate_id": cid,
                    "source": qr.source,
                    "abundance": f"{qr.abundance:.4f}",
                    "num_assigned_reads": qr.num_assigned_reads,
                    "confidence": f"{qr.confidence:.4f}",
                    "max_R": f"{qr.max_R:.4f}",
                    "coherence_score": f"{qr.coherence_score:.4f}",
                    "discrimination_score": f"{qr.discrimination_score:.4f}",
                    "combined_score": f"{qr.combined_score:.4f}",
                })
                total_rows += 1

    logger.info(
        "Wrote ablation summary: %s (%d transcript-rows across %d rows)",
        out_path,
        total_rows,
        len(results),
    )
    return str(out_path)


def write_per_row_tsv(
    ablation_result: AblationResult,
    output_dir: str,
) -> str:
    """Write a per-row TSV for a single AblationResult.

    Filename is ``<row_id>_<label>.tsv`` inside output_dir.

    Args:
        ablation_result: Single AblationResult to write.
        output_dir: Directory to write the TSV.

    Returns:
        Absolute path of the written file.
    """
    fname = f"{ablation_result.row_id}_{ablation_result.label}.tsv"
    out_path = Path(output_dir) / fname
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "candidate_id",
        "source",
        "abundance",
        "num_assigned_reads",
        "confidence",
        "max_R",
        "coherence_score",
        "discrimination_score",
        "combined_score",
    ]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for cid, qr in sorted(ablation_result.aggregated.items()):
            writer.writerow({
                "candidate_id": cid,
                "source": qr.source,
                "abundance": f"{qr.abundance:.4f}",
                "num_assigned_reads": qr.num_assigned_reads,
                "confidence": f"{qr.confidence:.4f}",
                "max_R": f"{qr.max_R:.4f}",
                "coherence_score": f"{qr.coherence_score:.4f}",
                "discrimination_score": f"{qr.discrimination_score:.4f}",
                "combined_score": f"{qr.combined_score:.4f}",
            })

    logger.info(
        "Wrote %s: %d transcripts", out_path, len(ablation_result.aggregated)
    )
    return str(out_path)
