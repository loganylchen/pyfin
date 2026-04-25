"""BEDPE writer for fusion transcript quantification results."""

from __future__ import annotations

import logging
from typing import Dict, Iterable, Union

from fin.analysis.quantification import QuantResult

logger = logging.getLogger(__name__)


def write_fusion_bedpe(
    results: Union[Dict[str, QuantResult], Iterable[QuantResult]],
    path: str,
) -> None:
    """Write fusion QuantResult entries to a BEDPE file.

    Only entries with source=='fusion' are written. Entries missing
    breakpoint_left or breakpoint_right are skipped with a warning.

    BEDPE columns (tab-separated, no header):
        chromA  startA  endA  chromB  startB  endB  name  score  strandA  strandB

    Score is int(round(combined_score * 1000)) clamped to [0, 1000].
    Rows are sorted by (chromA, startA) for determinism.

    Args:
        results: Dict mapping candidate_id -> QuantResult, or any iterable
                 of QuantResult objects.
        path: Output file path.
    """
    if isinstance(results, dict):
        items: Iterable[QuantResult] = results.values()
    else:
        items = results

    rows = []
    for qr in items:
        if qr.source != "fusion":
            continue

        if qr.breakpoint_left is None or qr.breakpoint_right is None:
            logger.warning(
                "Skipping fusion entry '%s': missing breakpoint(s) "
                "(breakpoint_left=%r, breakpoint_right=%r)",
                qr.candidate_id,
                qr.breakpoint_left,
                qr.breakpoint_right,
            )
            continue

        chrom_a, pos_a, strand_a = qr.breakpoint_left
        chrom_b, pos_b, strand_b = qr.breakpoint_right

        start_a = pos_a
        end_a = pos_a + 1
        start_b = pos_b
        end_b = pos_b + 1

        raw_score = int(round(qr.combined_score * 1000))
        score = max(0, min(1000, raw_score))

        rows.append((chrom_a, start_a, end_a, chrom_b, start_b, end_b,
                     qr.candidate_id, score, strand_a, strand_b))

    rows.sort(key=lambda r: (r[0], r[1]))

    with open(path, "w") as fh:
        for row in rows:
            fh.write(
                f"{row[0]}\t{row[1]}\t{row[2]}\t"
                f"{row[3]}\t{row[4]}\t{row[5]}\t"
                f"{row[6]}\t{row[7]}\t{row[8]}\t{row[9]}\n"
            )
