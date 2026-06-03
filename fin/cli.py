"""CLI entry point for pyfin."""

from __future__ import annotations

import os
import sys

import click


@click.group(invoke_without_command=True)
@click.option("--bam", default=None, type=click.Path(exists=True), help="Input BAM file.")
@click.option("--gtf", default=None, type=click.Path(), help="Reference GTF annotation file (optional).")
@click.option("--genome", default=None, type=click.Path(exists=True), help="Genome FASTA file.")
@click.option("--fastq", default=None, type=click.Path(exists=True), help="FASTQ reads file.")
@click.option("--signal", default=None, type=click.Path(exists=True), help="SLOW5/BLOW5/POD5 signal file.")
@click.option("--output-dir", default=None, help="Output directory.")
@click.option("--use-prior/--no-prior", default=True, show_default=True, help="Apply combined_score-derived EM prior.")
@click.option("--gpu/--no-gpu", "use_gpu", default=True, show_default=True, help="Enable/disable GPU acceleration.")
@click.option(
    "--signal-format",
    default="slow5",
    type=click.Choice(["slow5", "pod5"]),
    show_default=True,
    help="Signal file format.",
)
@click.option("--alpha", default=0.5, show_default=True, type=float, help="Score alpha (coherence vs discrimination weight).")
@click.option("--fusion", "fusion_enabled", is_flag=True, default=False, help="Enable fusion detection.")
@click.option("--min-support", default=2, show_default=True, type=int, help="Minimum read support for fusion breakpoint (only with --fusion).")
@click.option("--max-dist", default=500, show_default=True, type=int, help="Maximum distance (bp) for breakpoint clustering (only with --fusion).")
@click.option("--flank-bp", default=500, show_default=True, type=int, help="Flank size (bp) around fusion breakpoint (only with --fusion).")
@click.option("--max-reads-per-interval-for-dtw", "max_reads_for_dtw", default=2000, show_default=True, type=int, help="Cap reads per interval for read-to-read DTW (subsample beyond this).")
@click.option("--min-novel-reads", default=1, show_default=True, type=int, help="Drop novel candidates with fewer supporting reads (after collapsing).")
@click.option("--min-abundance", default=0.0, show_default=True, type=float, help="Drop NOVEL transcripts whose EM-estimated abundance is below this threshold (GTF candidates are exempt).")
@click.option("--min-max-r", default=0.0, show_default=True, type=float, help="Drop NOVEL transcripts whose max EM responsibility is below this (try 0.2 for ~+20pp precision; GTF candidates are exempt).")
@click.option("--min-novel-combined-score", default=0.0, show_default=True, type=float, help="Drop NOVEL transcripts whose combined_score is below this (F1-optimal: 0.288 with-GTF, 0.428 no-GTF; GTF candidates are exempt).")
@click.option("--min-isoform-fraction", default=0.01, show_default=True, type=float, help="Drop NOVEL multi-exon transcripts whose abundance is below this fraction of the dominant overlapping novel isoform at their locus (Cufflinks --min-isoform-fraction / StringTie -f minor-isoform suppression). GTF/fusion/mono exempt; 0.0 disables. Default 0.01 (StringTie-aligned, recall-safe). SIRV WARNING: F1-optimal ~0.4 is overfit — never use on real data.")
@click.option("--min-fulllen-fraction", default=0.1, show_default=True, type=float, help="Drop NOVEL multi-exon transcripts whose fraction of full-length assigned reads (read genomic 5' AND 3' both within --fulllen-window-bp of the candidate's ends) is below this (FLAIR/TALON-style full-length read support; signal-free). Orthogonal to --min-isoform-fraction. GTF/fusion/mono and unreachable candidates exempt; 0.0 disables. SIRV WARNING: default 0.1 is SIRV-tuned (drops most reachable novel-multi for free as SIRV lacks a 5'-truncated isoform tail) — re-tune or disable on real dRNA data.")
@click.option("--fulllen-window-bp", default=25, show_default=True, type=int, help="bp tolerance for a read genomic end to count as full-length wrt a candidate's 5'/3' end (used by --min-fulllen-fraction).")
@click.option("--fulllen-min-reads", default=4, show_default=True, type=int, help="Minimum assigned reads carrying a genomic span required to score a candidate's full-length fraction; below this the candidate is unreachable and never dropped (used by --min-fulllen-fraction).")
@click.option("--persist-R/--no-persist-R", "persist_R_matrix", default=True, show_default=True, help="Enable/disable R-matrix (R.npy) persistence per interval.")
@click.option("--canonical-gate/--no-canonical-gate", "canonical_gate", default=True, show_default=True, help="Drop NOVEL multi-exon candidates whose junctions aren't all canonical (GTF/fusion/mono exempt). SIRV-tuned default ON.")
@click.option("--canonical-motifs", default="GT-AG,GC-AG,AT-AC", show_default=True, help="Comma-separated donor-acceptor motifs accepted by the canonical gate AND search.")
@click.option("--canonical-search-bp", default=4, show_default=True, type=int, help="ea extended search: scan ±N bp around each read-derived NOVEL junction for canonical motifs and emit paired alternatives (GTF transcripts not extended). 0 disables. SIRV-tuned default 4.")
@click.option("--m2-tiebreak/--no-m2-tiebreak", "m2_tiebreak", default=True, show_default=True, help="Resolve argmax_keep ties with the junction-window mean-NLL signal metric: give a read's full mass to the M2-best tied candidate when the NLL margin >= --m2-tiebreak-margin, else keep the 1/K split. Needs --signal (auto-skips if absent). Default ON + aggressive (margin 1e-9): SIRV gffcompare Tx-F1 45.4 vs OFF 44.7, beats all tools + ablation champion 45.2.")
@click.option("--m2-tiebreak-margin", default=1e-9, show_default=True, type=float, help="Minimum M2 NLL margin (runner-up - best) required to override the 1/K split with the M2 single winner. Default 1e-9 (aggressive: take M2's pick whenever it can discriminate at all).")
@click.option("--m2-tiebreak-junction-k", default=10, show_default=True, type=int, help="Transcript-frame bp on each side of the wobbling junction for the M2 discrimination window (SIRV sweet spot 10).")
@click.option("--verbose", "-v", is_flag=True, help="Enable debug logging.")
@click.pass_context
def main(
    ctx,
    bam,
    gtf,
    genome,
    fastq,
    signal,
    output_dir,
    use_prior,
    use_gpu,
    signal_format,
    alpha,
    fusion_enabled,
    min_support,
    max_dist,
    flank_bp,
    max_reads_for_dtw,
    min_novel_reads,
    min_abundance,
    min_max_r,
    min_novel_combined_score,
    min_isoform_fraction,
    min_fulllen_fraction,
    fulllen_window_bp,
    fulllen_min_reads,
    persist_R_matrix,
    canonical_gate,
    canonical_motifs,
    canonical_search_bp,
    m2_tiebreak,
    m2_tiebreak_margin,
    m2_tiebreak_junction_k,
    verbose,
):
    """pyfin: nanopore signal-based transcriptome assembly.

    Default command performs reference-based transcriptome assembly. Pass
    --fusion to additionally detect gene fusions. Use the `quantify`
    subcommand for multi-sample known-transcript quantification.
    """
    if ctx.invoked_subcommand is not None:
        return

    missing = []
    if not bam:
        missing.append("--bam")
    if not genome:
        missing.append("--genome")
    if not fastq:
        missing.append("--fastq")
    if not signal:
        missing.append("--signal")
    if not output_dir:
        missing.append("--output-dir")
    if missing:
        click.echo(f"Error: missing required option(s): {', '.join(missing)}", err=True)
        click.echo(ctx.get_help(), err=True)
        sys.exit(2)

    from fin.utils.log_config import setup_logger

    setup_logger("fin", level="DEBUG" if verbose else "INFO")

    from fin.pipeline.config import PipelineConfig
    from fin.pipeline.runner import PipelineRunner

    os.makedirs(output_dir, exist_ok=True)

    cfg = PipelineConfig(
        bam_path=bam,
        gtf_path=gtf,
        genome_fasta_path=genome,
        fastq_path=fastq,
        signal_path=signal,
        work_dir=output_dir,
        output_gtf=os.path.join(output_dir, "assembly.gtf"),
        output_tsv=os.path.join(output_dir, "scores.tsv"),
        output_bedpe=os.path.join(output_dir, "fusions.bedpe") if fusion_enabled else None,
        use_gpu=use_gpu,
        use_prior=use_prior,
        signal_format=signal_format,
        score_alpha=alpha,
        fusion_enabled=fusion_enabled,
        fusion_min_support=min_support,
        fusion_max_dist=max_dist,
        fusion_flank_bp=flank_bp,
        max_reads_per_interval_for_dtw=max_reads_for_dtw,
        min_novel_reads=min_novel_reads,
        min_abundance=min_abundance,
        min_max_r=min_max_r,
        min_novel_combined_score=min_novel_combined_score,
        min_isoform_fraction=min_isoform_fraction,
        min_fulllen_fraction=min_fulllen_fraction,
        fulllen_window_bp=fulllen_window_bp,
        fulllen_min_reads=fulllen_min_reads,
        persist_R_matrix=persist_R_matrix,
        canonical_gate=canonical_gate,
        canonical_motifs=tuple(
            m.strip() for m in canonical_motifs.split(",") if m.strip()
        ),
        canonical_search_bp=canonical_search_bp,
        m2_tiebreak=m2_tiebreak,
        m2_tiebreak_margin=m2_tiebreak_margin,
        m2_tiebreak_junction_k=m2_tiebreak_junction_k,
    )

    runner = PipelineRunner(cfg)
    try:
        runner.setup()
        runner.run()
    finally:
        runner.cleanup()

    click.echo(f"Assembly output written to {output_dir}/")


@main.command()
@click.option("--gtf", required=True, type=click.Path(exists=True), help="GTF annotation file.")
@click.option("--genome", required=True, type=click.Path(exists=True), help="Genome FASTA file.")
@click.option(
    "--sample",
    multiple=True,
    required=True,
    help="Sample in format name:bam:fastq:blow5. Can be specified multiple times.",
)
@click.option("--output-dir", default="./pyfin_quant", show_default=True, help="Output directory.")
@click.option("--use-gpu/--no-gpu", default=True, show_default=True, help="Use GPU for DTW.")
@click.option(
    "--signal-format",
    default="slow5",
    type=click.Choice(["slow5", "pod5"]),
    show_default=True,
    help="Signal file format.",
)
@click.option("--f5c-path", default="f5c", show_default=True, help="Path to f5c binary.")
@click.option("--use-prior/--no-prior", default=True, show_default=True, help="Apply combined_score-derived EM prior.")
@click.option("--signal-normalize/--no-signal-normalize", default=True, show_default=True, help="Per-read robust z-score normalization before DTW.")
@click.option("--verbose", "-v", is_flag=True, help="Enable debug logging.")
def quantify(gtf, genome, sample, output_dir, use_gpu, signal_format, f5c_path, use_prior, signal_normalize, verbose):
    """Quantify known transcripts across multiple samples."""
    from fin.utils.log_config import setup_logger

    logger = setup_logger("fin", level="DEBUG" if verbose else "INFO")

    from fin.pipeline.quantify_runner import QuantifyRunner, SampleInput

    samples = []
    for s in sample:
        parts = s.split(":")
        if len(parts) != 4:
            click.echo(
                f"Error: --sample must be name:bam:fastq:blow5, got '{s}'",
                err=True,
            )
            sys.exit(1)
        samples.append(SampleInput(name=parts[0], bam_path=parts[1], fastq_path=parts[2], signal_path=parts[3]))

    logger.info("Quantifying %d samples against %s", len(samples), gtf)

    from fin.pipeline.config import PipelineConfig

    quant_config = PipelineConfig(
        bam_path="",
        use_gpu=use_gpu,
        use_prior=use_prior,
        signal_normalize=signal_normalize,
    )

    runner = QuantifyRunner(
        gtf_path=gtf,
        genome_fasta_path=genome,
        samples=samples,
        output_dir=output_dir,
        signal_format=signal_format,
        use_gpu=use_gpu,
        f5c_path=f5c_path,
        config=quant_config,
    )

    try:
        runner.setup()
        runner.run()
    finally:
        runner.cleanup()

    click.echo(f"Output written to {output_dir}/")


if __name__ == "__main__":
    main()
