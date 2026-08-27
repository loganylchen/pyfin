#!/usr/bin/env python
"""Ablation benchmark: run the active quantification rows on an interval set.

Requires a real dataset (BAM + signal + optional GTF). Writes per-row TSVs
and a combined summary to --out-dir. Intended to be run on a 5-row SIRV or
SG-NEx subset; not a unit test.

CLI example:
    python benchmarks/ablation_5rows.py \\
        --bam       /data/sirv/reads.bam \\
        --fastq     /data/sirv/reads.fastq \\
        --signal    /data/sirv/reads.slow5 \\
        --gtf       /data/sirv/sirv.gtf \\
        --genome    /data/sirv/genome.fa \\
        --out-dir   /tmp/ablation_5rows \\
        --max-reads 500
"""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from fin.pipeline.config import PipelineConfig

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("ablation_5rows")


def _build_base_config(args) -> "PipelineConfig":
    from fin.pipeline.config import PipelineConfig

    return PipelineConfig(
        bam_path=args.bam,
        gtf_path=args.gtf or None,
        genome_fasta_path=args.genome,
        fastq_path=args.fastq,
        signal_path=args.signal,
        signal_format=args.signal_format,
        work_dir=str(Path(args.out_dir) / "work"),
        use_gpu=not args.no_gpu,
        max_reads=args.max_reads,
        em_sigma=1.0,
        em_max_iter=200,
        use_prior=True,
        score_alpha=0.5,
        signal_normalize=True,
        persist_R_matrix=False,
        # Score filter thresholds — only active when enable_score_filter=True (R5)
        min_abundance=0.5,
        # Ablation-specific overrides set per row below
        quant_mode="m2_em",
        enable_score_filter=True,
    )


def _row_config(base, row_cfg) -> "PipelineConfig":
    """Return a shallow copy of base with ablation row overrides applied."""
    import dataclasses

    overrides = {
        "quant_mode": row_cfg.quant_mode,
        "em_max_iter_override": row_cfg.em_max_iter_override,
        "enable_score_filter": row_cfg.enable_score_filter,
    }
    return dataclasses.replace(base, **overrides)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--signal", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", default=None)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--signal-format", default="slow5")
    ap.add_argument("--max-reads", type=int, default=None)
    ap.add_argument("--rows", nargs="+", default=None,
                    help="Subset of row IDs to run, e.g. --rows R1 R3")
    ap.add_argument("--no-gpu", action="store_true")
    args = ap.parse_args(argv)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    from fin.ablation.io import write_ablation_summary, write_per_row_tsv
    from fin.ablation.runner import ABLATION_ROWS, AblationResult
    from fin.pipeline.runner import PipelineRunner

    base_config = _build_base_config(args)

    rows_to_run = ABLATION_ROWS
    if args.rows:
        rows_to_run = [r for r in ABLATION_ROWS if r.row_id in args.rows]
        if not rows_to_run:
            logger.error("No matching rows for: %s", args.rows)
            sys.exit(1)

    all_results: list[AblationResult] = []

    for row_cfg in rows_to_run:
        logger.info("=== Running ablation row %s: %s ===", row_cfg.row_id, row_cfg.label)
        cfg = _row_config(base_config, row_cfg)

        runner = PipelineRunner(cfg)
        runner.setup()

        t0 = time.perf_counter()
        aggregated = runner.run()
        elapsed = time.perf_counter() - t0
        runner.cleanup()

        logger.info(
            "Row %s done in %.1fs: %d transcripts quantified",
            row_cfg.row_id,
            elapsed,
            len(aggregated),
        )

        ablation_result = AblationResult(
            row_id=row_cfg.row_id,
            label=row_cfg.label,
            aggregated=aggregated,
        )
        all_results.append(ablation_result)

        tsv_path = write_per_row_tsv(ablation_result, str(out_dir))
        logger.info("Wrote: %s", tsv_path)

        # Write per-row GTF for intron-chain evaluation against ground truth
        from fin.io.io_gtf import write_gtf
        gtf_path = str(out_dir / f"{row_cfg.row_id}_{row_cfg.label}.gtf")
        write_gtf(aggregated, gtf_path)
        logger.info("Wrote GTF: %s", gtf_path)

    summary_path = write_ablation_summary(all_results, str(out_dir))
    logger.info("Summary written: %s", summary_path)

    # Print quick stats table to stdout
    header = f"{'row_id':<6} {'label':<28} {'n_transcripts':>14} {'n_novel':>10}"
    print(header)
    print("-" * len(header))
    for ar in all_results:
        n_total = len(ar.aggregated)
        n_novel = sum(1 for qr in ar.aggregated.values() if qr.source not in ("gtf", "fusion"))
        print(f"{ar.row_id:<6} {ar.label:<28} {n_total:>14} {n_novel:>10}")


if __name__ == "__main__":
    main()
