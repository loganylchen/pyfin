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
@click.option("--no-gpu", "use_gpu", is_flag=True, default=True, flag_value=False, help="Disable GPU acceleration.")
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
@click.option("--verbose", "-v", is_flag=True, help="Enable debug logging.")
def quantify(gtf, genome, sample, output_dir, use_gpu, signal_format, f5c_path, use_prior, verbose):
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
