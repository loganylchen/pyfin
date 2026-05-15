#!/usr/bin/env python
"""EM matrix-subset ablation: 7 subsets × {with-GTF, no-GTF} = 14 runs.

Drives `PipelineRunner` once per cell of the 7×2 grid (em_matrix_subset
∈ {m1, m2, m3, m1+m2, m1+m3, m2+m3, all} × {with-GTF, no-GTF}). Each run
writes its own GTF + scoring TSV under ``--out-dir`` and contributes a
row to ``--out-dir/summary.tsv``.

Output layout (under --out-dir):
    {subset_safe}_{mode}/work/                # per-run pyfin work dir
    {subset_safe}_{mode}.gtf                  # assembled GTF (intron-chain eval input)
    {subset_safe}_{mode}.scoring.tsv          # per-candidate scores
    summary.tsv                               # one row per (subset, mode) cell

Evaluate with ``benchmarks/em_subset_evaluate.py``.

CLI example:
    python benchmarks/em_subset_ablation.py \\
        --bam       /data/sirv/reads.bam \\
        --fastq     /data/sirv/reads.fastq \\
        --signal    /data/sirv/reads.slow5 \\
        --gtf       /data/sirv/sirv.gtf \\
        --genome    /data/sirv/genome.fa \\
        --out-dir   /tmp/em_subset \\
        --subsets   m2+m3 all
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import logging
import sys
import time
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("em_subset_ablation")


ALL_SUBSETS = ("m1", "m2", "m3", "m1+m2", "m1+m3", "m2+m3", "all")
ALL_MODES = ("with-gtf", "no-gtf")


def _safe(subset: str) -> str:
    """Filesystem-safe subset tag (m1+m3 -> m1_m3)."""
    return subset.replace("+", "_")


def _build_config(args, subset: str, mode: str, out_dir: Path):
    """Build a PipelineConfig for one cell of the 7x2 grid."""
    from fin.pipeline.config import PipelineConfig

    cell_tag = f"{_safe(subset)}_{mode}"
    cell_work = out_dir / cell_tag / "work"
    cell_work.mkdir(parents=True, exist_ok=True)

    gtf_in = args.gtf if mode == "with-gtf" else None

    return PipelineConfig(
        bam_path=args.bam,
        gtf_path=gtf_in,
        genome_fasta_path=args.genome,
        fastq_path=args.fastq,
        signal_path=args.signal,
        signal_format=args.signal_format,
        work_dir=str(cell_work),
        use_gpu=not args.no_gpu,
        max_reads=args.max_reads,
        em_sigma=args.em_sigma,
        em_beta=args.em_beta,
        em_max_iter=args.em_max_iter,
        em_tol=1e-4,
        use_prior=True,
        score_alpha=0.5,
        signal_normalize=True,
        persist_R_matrix=False,
        m4_source=args.m4_source,
        em_matrix_subset=subset,
        # Filters off for assembly-only evaluation (mirrors P/R-curve sweep).
        enable_score_filter=args.enable_score_filter,
        min_abundance=args.min_abundance,
        min_max_r=0.0,
        min_novel_combined_score=0.0,
        output_gtf=str(out_dir / f"{cell_tag}.gtf"),
        output_tsv=str(out_dir / f"{cell_tag}.scoring.tsv"),
    )


def _summarize(aggregated) -> dict:
    """Count totals/novel/gtf from a Dict[cid, QuantResult]."""
    n_total = len(aggregated)
    n_gtf = sum(1 for qr in aggregated.values() if qr.source == "gtf")
    n_novel = sum(1 for qr in aggregated.values() if qr.source == "novel")
    n_fusion = sum(1 for qr in aggregated.values() if qr.source == "fusion")
    return {
        "n_transcripts": n_total,
        "n_gtf": n_gtf,
        "n_novel": n_novel,
        "n_fusion": n_fusion,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--signal", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--gtf", default=None,
                    help="Reference GTF (used for with-gtf mode; required if --modes includes with-gtf)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--signal-format", default="slow5")
    ap.add_argument("--max-reads", type=int, default=None)
    ap.add_argument("--subsets", nargs="+", default=list(ALL_SUBSETS),
                    choices=list(ALL_SUBSETS),
                    help="Which em_matrix_subset values to sweep")
    ap.add_argument("--modes", nargs="+", default=list(ALL_MODES),
                    choices=list(ALL_MODES),
                    help="Which annotation modes to sweep")
    ap.add_argument("--m4-source", default="whole_read",
                    choices=("whole_read", "diff_region", "none"))
    ap.add_argument("--em-sigma", type=float, default=1.0)
    ap.add_argument("--em-beta", type=float, default=0.5)
    ap.add_argument("--em-max-iter", type=int, default=200)
    ap.add_argument("--enable-score-filter", action="store_true",
                    help="Apply min_abundance/min_max_r filters (default: off, raw assembly)")
    ap.add_argument("--min-abundance", type=float, default=0.0)
    ap.add_argument("--no-gpu", action="store_true")
    ap.add_argument("--skip-existing", action="store_true",
                    help="Skip cells whose summary row + GTF already exist")
    args = ap.parse_args(argv)

    if "with-gtf" in args.modes and not args.gtf:
        ap.error("--gtf is required when --modes includes with-gtf")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "summary.tsv"

    from fin.pipeline.runner import PipelineRunner

    # Fieldnames for summary.tsv
    fieldnames = [
        "subset", "mode", "elapsed_s",
        "n_transcripts", "n_gtf", "n_novel", "n_fusion",
        "gtf_out", "tsv_out",
    ]

    existing_rows: list[dict] = []
    completed_cells: set[tuple] = set()
    if summary_path.exists() and args.skip_existing:
        with open(summary_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                existing_rows.append(row)
                completed_cells.add((row["subset"], row["mode"]))

    new_rows: list[dict] = []
    grid = [(s, m) for s in args.subsets for m in args.modes]
    logger.info("Running %d cells: subsets=%s modes=%s",
                len(grid), args.subsets, args.modes)

    for subset, mode in grid:
        cell_tag = f"{_safe(subset)}_{mode}"
        if (subset, mode) in completed_cells:
            logger.info("[skip] %s — already present in summary.tsv", cell_tag)
            continue

        logger.info("=== Cell %s (subset=%s mode=%s) ===", cell_tag, subset, mode)
        cfg = _build_config(args, subset, mode, out_dir)

        runner = PipelineRunner(cfg)
        runner.setup()
        t0 = time.perf_counter()
        try:
            aggregated = runner.run()
        finally:
            elapsed = time.perf_counter() - t0
            runner.cleanup()

        counts = _summarize(aggregated)
        row = {
            "subset": subset,
            "mode": mode,
            "elapsed_s": f"{elapsed:.1f}",
            **{k: counts[k] for k in ("n_transcripts", "n_gtf", "n_novel", "n_fusion")},
            "gtf_out": cfg.output_gtf,
            "tsv_out": cfg.output_tsv,
        }
        new_rows.append(row)
        logger.info("Cell %s done in %.1fs: %s",
                    cell_tag, elapsed, {k: counts[k] for k in counts})

        # Persist summary incrementally so failures don't lose work.
        with open(summary_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            w.writeheader()
            for r in existing_rows + new_rows:
                w.writerow(r)

    logger.info("Summary written: %s (%d rows)", summary_path,
                len(existing_rows) + len(new_rows))

    # Print a compact table to stdout
    header = f"{'subset':<8} {'mode':<10} {'n_total':>8} {'n_gtf':>6} {'n_novel':>8} {'time(s)':>8}"
    print(header)
    print("-" * len(header))
    for r in existing_rows + new_rows:
        print(f"{r['subset']:<8} {r['mode']:<10} {r['n_transcripts']:>8} "
              f"{r['n_gtf']:>6} {r['n_novel']:>8} {r['elapsed_s']:>8}")


if __name__ == "__main__":
    main()
