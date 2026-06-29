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
@click.option("--gpu/--no-gpu", "use_gpu", default=True, show_default=True, help="Enable/disable GPU acceleration.")
@click.option(
    "--signal-format",
    default="slow5",
    type=click.Choice(["slow5", "pod5"]),
    show_default=True,
    help="Signal file format.",
)
@click.option("--fusion", "fusion_enabled", is_flag=True, default=False, help="Enable fusion detection.")
@click.option("--min-support", default=2, show_default=True, type=int, help="Minimum read support for fusion breakpoint (only with --fusion).")
@click.option("--max-dist", default=500, show_default=True, type=int, help="Maximum distance (bp) for breakpoint clustering (only with --fusion).")
@click.option("--flank-bp", default=500, show_default=True, type=int, help="Flank size (bp) around fusion breakpoint (only with --fusion).")
@click.option("--fusion-max-internal-gap-bp", default=30, show_default=True, type=int, help="Adapter-bridged chimera guard (only with --fusion): drop a chimeric read whose two arms are separated by an internal stretch of read sequence mapping to NEITHER arm of >= this many bp (typically a nanopore ONT internal adapter — a false fusion, cf. DeepChopper). 0 disables.")
@click.option("--max-reads-per-interval-for-dtw", "max_reads_for_dtw", default=2000, show_default=True, type=int, help="Cap reads per interval for read-to-read DTW (subsample beyond this).")
@click.option("--min-novel-reads", default=1, show_default=True, type=int, help="Drop novel candidates with fewer supporting reads (after collapsing).")
@click.option("--min-abundance", default=3.0, show_default=True, type=float, help="Drop NOVEL transcripts whose EM-estimated abundance is below this threshold (GTF uses its own lighter floor --min-gtf-abundance unless --floor-gtf-abundance; fusion always exempt). Default 3 reproduces the SIRV gffcompare T>=3 NOVEL operating point. SIRV-tuned WARNING: a nonzero floor guts recall of genuine low-abundance isoforms — on real data set 0 to disable.")
@click.option("--min-gtf-abundance", default=1.0, show_default=True, type=float, help="Abundance floor for GTF-annotated transcripts on soft EM abundance (default 1: a GTF candidate must accrue at least one read of EM soft-mass to be kept, so pyfin is not a copy-annotation tool). Independent of --min-abundance (which gates novel); fusion always exempt. With --floor-gtf-abundance the effective GTF floor becomes max(this, --min-abundance). Set 0 to disable and keep every GTF candidate in an expressed locus. SIRV-tuned: on real data a nonzero floor can drop genuine low-abundance annotated isoforms.")
@click.option("--floor-gtf-abundance/--no-floor-gtf-abundance", "floor_gtf_abundance", default=False, show_default=True, help="Raise the GTF abundance floor to the NOVEL floor: GTF must clear max(--min-gtf-abundance, --min-abundance) (never lowers the explicit GTF floor). Default OFF: GTF uses its own lighter --min-gtf-abundance. Turn ON to reproduce the SIRV gffcompare T>=3 WITH-GTF operating point (e.g. T=3 -> Sn 43.2 / Pr 91.6). SIRV-tuned: on real data this drops genuine low-abundance annotated isoforms — leave OFF unless benchmarking. Fusion stays exempt regardless.")
@click.option("--min-isoform-fraction", default=0.01, show_default=True, type=float, help="Drop NOVEL multi-exon transcripts whose abundance is below this fraction of the dominant overlapping novel isoform at their locus (Cufflinks --min-isoform-fraction / StringTie -f minor-isoform suppression). GTF/fusion/mono exempt; 0.0 disables. Default 0.01 (StringTie-aligned, recall-safe). SIRV WARNING: F1-optimal ~0.4 is overfit — never use on real data.")
@click.option("--min-fulllen-fraction", default=0.1, show_default=True, type=float, help="Drop NOVEL multi-exon transcripts whose fraction of full-length assigned reads (read genomic 5' AND 3' both within --fulllen-window-bp of the candidate's ends) is below this (FLAIR/TALON-style full-length read support; signal-free). Orthogonal to --min-isoform-fraction. GTF/fusion/mono and unreachable candidates exempt; 0.0 disables. SIRV WARNING: default 0.1 is SIRV-tuned (drops most reachable novel-multi for free as SIRV lacks a 5'-truncated isoform tail) — re-tune or disable on real dRNA data.")
@click.option("--max-soft-mass-ratio", default=2.0, show_default=True, type=float, help="Drop NOVEL multi-exon transcripts whose EM soft abundance / hard argmax read count is >= this (0 hard reads always dropped). Real isoforms deposit ~1 soft mass per hard read (ratio ~1); wobble shadows borrow fractional soft crumbs from a high-abundance structural near-copy and show an inflated ratio — catches the high-relative-abundance shadows the fraction/cluster-recheck levers miss. Pure EM evidence (no GTF). GTF/fusion/mono exempt; 0.0 disables. SIRV-tuned default 2.0 (24-cell F1@3 up-or-equal, recall held); re-tune on real dRNA.")
@click.option("--fulllen-window-bp", default=25, show_default=True, type=int, help="bp tolerance for a read genomic end to count as full-length wrt a candidate's 5'/3' end (used by --min-fulllen-fraction).")
@click.option("--fulllen-min-reads", default=4, show_default=True, type=int, help="Minimum assigned reads carrying a genomic span required to score a candidate's full-length fraction; below this the candidate is unreachable and never dropped (used by --min-fulllen-fraction).")
@click.option("--min-polya5p-reads", default=1, show_default=True, type=int, help="Drop a candidate unless >= N of its assigned reads BOTH have a krill whole-read polyA tail (qc PASS & length > --min-polya-length) AND map with their genomic 5' end within --polya5p-window-bp of the candidate's 5' end. By default only novel candidates are gated (GTF and fusion exempt; pass --no-polya5p-exempt-gtf to gate GTF too). Needs --signal; adds a krill polyA pass to every run; 0 disables. SIRV WARNING: ON by default (=1); lifts no-GTF Tx-F1 — re-tune or disable on real dRNA.")
@click.option("--polya5p-window-bp", default=25, show_default=True, type=int, help="bp tolerance for a read's genomic 5' end to count toward --min-polya5p-reads.")
@click.option("--min-polya-length", default=10.0, show_default=True, type=float, help="Minimum krill polya_length (with polya_qc PASS) for a read to support a candidate under --min-polya5p-reads.")
@click.option("--polya5p-exempt-gtf/--no-polya5p-exempt-gtf", "polya5p_exempt_gtf", default=True, show_default=True, help="Exempt GTF-sourced candidates from the --min-polya5p-reads filter (like fusion); only novel candidates are gated. Default ON: validated as the corrected-recall-optimal config on SIRV (dRNA reads are often 5'-truncated/lack polyA, so gating GTF drops genuine annotated transcripts). Pass --no-polya5p-exempt-gtf to gate GTF candidates too.")
@click.option("--persist-R/--no-persist-R", "persist_R_matrix", default=True, show_default=True, help="Enable/disable R-matrix (R.npy) persistence per interval.")
@click.option("--canonical-gate/--no-canonical-gate", "canonical_gate", default=True, show_default=True, help="Drop NOVEL multi-exon candidates whose junctions aren't all canonical (GTF/fusion/mono exempt). SIRV-tuned default ON.")
@click.option("--canonical-motifs", default="GT-AG,GC-AG,AT-AC", show_default=True, help="Comma-separated donor-acceptor motifs accepted by the canonical gate AND search.")
@click.option("--canonical-search-bp", default=4, show_default=True, type=int, help="ea extended search: scan ±N bp around each read-derived NOVEL junction for canonical motifs and emit paired alternatives (GTF transcripts not extended). 0 disables. SIRV-tuned default 4.")
@click.option("--m2-tiebreak/--no-m2-tiebreak", "m2_tiebreak", default=True, show_default=True, help="Resolve argmax_keep ties with the junction-window mean-NLL signal metric: give a read's full mass to the M2-best tied candidate when the NLL margin >= --m2-tiebreak-margin, else keep the 1/K split. Needs --signal (auto-skips if absent). Default ON + aggressive (margin 1e-9): SIRV gffcompare Tx-F1 45.4 vs OFF 44.7, beats all tools + ablation champion 45.2.")
@click.option("--m2-tiebreak-margin", default=1e-9, show_default=True, type=float, help="Minimum M2 NLL margin (runner-up - best) required to override the 1/K split with the M2 single winner. Default 1e-9 (aggressive: take M2's pick whenever it can discriminate at all).")
@click.option("--m2-tiebreak-junction-k", default=10, show_default=True, type=int, help="Transcript-frame bp on each side of the wobbling junction for the M2 discrimination window (SIRV sweet spot 10).")
@click.option("--m2-diff-cover-gate/--no-m2-diff-cover-gate", "m2_diff_cover_gate", default=True, show_default=True, help="(--quant-mode m2_em only) Diff-region coverage gate. For a read in a >=2 best-AS tie: if the M2 NLL margin >= --m2-diff-cover-margin hard-assign its full mass to the best candidate (and, if it straddles every wobbling junction donor->acceptor, contribute that vote to the locus isoform-ratio prior); otherwise (margin below threshold, whether or not it covers) the read is ambiguous and its tie mass is redistributed in proportion to the prior learned from covered+distinguishing reads (flat 1/K only when the tie has no covered prior). No read is dropped (recall-safe). Default ON; --no-m2-diff-cover-gate reverts to the soft NLL-graded d_tx. SIRV/dense-locus precision lever; re-tune the margin on real dRNA.")
@click.option("--m2-diff-cover-margin", default=0.5, show_default=True, type=float, help="(--m2-diff-cover-gate) Minimum M2 NLL margin (runner-up - best) to HARD-assign a tied read to its lowest-NLL candidate. Below this the read is redistributed by the covered-read prior (flat 1/K when no prior exists). SIRV-tuned; sweep on real data.")
@click.option("--m2-cluster-recheck/--no-m2-cluster-recheck", "m2_cluster_recheck", default=True, show_default=True, help="(--quant-mode m2_em only) After EM, cluster ALL multi-exon candidates by structure (same intron count + every junction within --m2-cluster-recheck-bp); within a cluster the highest-abundance candidate anchors (often the true isoform / a GTF passthrough) and a NOVEL sibling whose abundance < --m2-cluster-recheck-fraction of the anchor is dropped as a wobble shadow (GTF/fusion never dropped). Pure abundance evidence — GTF only joins the abundance race, never used as a correctness oracle, so it stays robust on corrupted/absent annotation. Default ON; --no-m2-cluster-recheck disables (no drops).")
@click.option("--m2-cluster-recheck-bp", default=20, show_default=True, type=int, help="(--m2-cluster-recheck) Per-junction wobble tolerance (bp): two candidates cluster iff same intron count and every donor/acceptor within this many bp. Default 20 (validated on the 5-sample heya8 + sirv4 full sweep: 40/40 cells precision up vs bp=10).")
@click.option("--m2-cluster-recheck-fraction", default=0.15, show_default=True, type=float, help="(--m2-cluster-recheck) Relative-abundance threshold; a novel cluster sibling is a shadow when its EM abundance is below this fraction of the cluster's highest-abundance anchor. 0 falls back to --min-isoform-fraction. SIRV-tuned; sweep on real data.")
@click.option("--m2-cluster-recheck-cassette-max-exon-bp", default=70, show_default=True, type=int, help="(--m2-cluster-recheck) Extra cassette-skip equivalence: cluster a K-intron candidate with a (K-1)-intron candidate if they align with one extra exon shorter than this many bp on the K side, all other junctions within --m2-cluster-recheck-bp. Targets minimap2's small-exon-skip artifacts (one read pop got the exon, another collapsed it into one long intron). 0 disables. Default 70 (fires on heya8 p00 and c_jitter10bp where bp alone leaves cassette FPs).")
@click.option("--m2-cluster-recheck-novel-displaces-gtf/--no-m2-cluster-recheck-novel-displaces-gtf", "m2_cluster_recheck_novel_displaces_gtf", default=True, show_default=True, help="(--m2-cluster-recheck) Allow a clustered low-support GTF sibling (below --m2-cluster-recheck-fraction of the anchor) to be dropped as a wobble shadow, gated by a direct read-support guard (see --m2-cluster-recheck-gtf-min-jct-reads). Targets jittered/mis-annotated GTF passthroughs whose junction no read traverses. Judged by the GTF's OWN read support, not the anchor source. Default ON; --no-... restores the legacy 'GTF never dropped' behaviour.")
@click.option("--m2-cluster-recheck-gtf-min-jct-reads", default=1, show_default=True, type=int, help="(--m2-cluster-recheck) A clustered low-abundance GTF sibling is dropped only if its distinguishing junction has FEWER than this many reads splicing exactly there. 1 = drop only zero-exact-support (jittered/phantom) GTFs; a real annotated isoform keeps its own reads and survives.")
@click.option("--m2-cluster-recheck-jct-tol", default=0, show_default=True, type=int, help="(--m2-cluster-recheck) bp tolerance matching a GTF candidate's junction to an observed read junction for the read-support guard. Default 0 (strict exact match): validated on heya8 c_jitter10bp — a loose tol lets 1-2bp jitters borrow read support from the true neighbouring junction.")
@click.option("--m2-support-gate/--no-m2-support-gate", "m2_support_gate", default=True, show_default=True, help="(--quant-mode m2_em only) Keep a multi-exon candidate (GTF or novel; fusion/mono exempt) only if it earns >=1 read's support: it is some read's M1 SOLE best-AS, OR some read's M2-best (lowest junction-NLL in its tie). Drops corrupted-annotation junctions that win neither. Default ON (tie-accept): lifts c_jitter F1@3 with zero recall loss; --no-m2-support-gate disables.")
@click.option("--m2-support-gate-tie/--no-m2-support-gate-tie", "m2_support_gate_tie", default=True, show_default=True, help="(--m2-support-gate) Count a candidate tied for the lowest M2 NLL as M2-best (recall-safer). --no-... requires a strict unique M2 win.")
@click.option("--containment-collapse/--no-containment-collapse", "containment_collapse", default=False, show_default=True, help="(--quant-mode m2_em only) Lever 1: after EM, fold a NOVEL candidate whose intron chain is a pure 3' SUFFIX of a longer candidate (a 5'-truncation shadow: same downstream junctions, 3' terminus within --containment-3p-tol-bp, 5' end interior, EM abundance <= parent * --containment-min-abundance-ratio) into that longer parent (reads + soft mass reassigned, shadow dropped). gtf candidates are never folded away. STRICT suffix match -> exon-skip/alt-3'-end isoforms never folded. WARNING: cannot distinguish 5'-truncation from a genuine low-abundance alt-TSS isoform -> NOT recall-safe; default OFF, enable only after real-data validation.")
@click.option("--containment-3p-tol-bp", default=20, show_default=True, type=int, help="(--containment-collapse) bp tolerance for the shadow's 3' terminus matching the parent's 3' terminus.")
@click.option("--containment-min-abundance-ratio", default=1.0, show_default=True, type=float, help="(--containment-collapse) Fold a shadow only when its EM abundance <= parent's * this ratio (shadow must be the minor member). 1.0 = fold any shadow at or below the parent.")
@click.option("--drop-mono-exon-novel/--no-drop-mono-exon-novel", "drop_mono_exon_novel", default=False, show_default=True, help="Lever 3: drop NOVEL single-exon (mono) candidates with weak support — hard read count < --min-mono-exon-reads OR genomic length < --min-mono-exon-length. Suppresses single-exon de novo noise (like IsoQuant's ONT default) WITHOUT a blanket drop: a high-support/long real intronless gene survives. gtf/fusion/multi-exon exempt. Needs at least one threshold > 0 to fire. Default OFF.")
@click.option("--min-mono-exon-reads", default=0, show_default=True, type=int, help="(--drop-mono-exon-novel) Min hard reads for a novel mono candidate to survive. 0 disables this threshold.")
@click.option("--min-mono-exon-length", default=0, show_default=True, type=int, help="(--drop-mono-exon-novel) Min genomic length (bp) for a novel mono candidate to survive. 0 disables this threshold.")
@click.option("--novel-junction-min-reads", default=0, show_default=True, type=int, help="(--quant-mode m2_em only) Lever 2: drop a NOVEL multi-exon candidate if ANY of its junctions is spliced by fewer than N directly-observed reads (primary-read CIGAR introns, strand-keyed, within --novel-junction-reads-tol bp). I.e. a novel junction must be carried by >= N reads, not just 1. gtf/fusion/mono exempt. <=1 disables (byte-identical). Default 0; set 2 to require >=2 reads/junction.")
@click.option("--novel-junction-reads-tol", default=2, show_default=True, type=int, help="(--novel-junction-min-reads) bp tolerance matching a candidate junction to an observed read junction.")
@click.option("--junction-dominance-filter/--no-junction-dominance-filter", "junction_dominance_filter", default=False, show_default=True, help="(--quant-mode m2_em only) Junction-first PRE-EM gate: drop a NOVEL multi-exon candidate if ANY junction has < --junction-dominance-min-reads observed reads OR is not locally dominant (a different observed junction within --junction-dominance-window-bp carries strictly more reads). Removes multi-read wobble shadows (they lose to the stronger true junction nearby) BEFORE EM, so shadows never compete for reads. Pure mapping evidence, no snap. gtf/fusion/mono exempt. Default OFF.")
@click.option("--junction-dominance-min-reads", default=2, show_default=True, type=int, help="(--junction-dominance-filter) Min observed reads for a novel junction to survive.")
@click.option("--junction-dominance-window-bp", default=20, show_default=True, type=int, help="(--junction-dominance-filter) Neighborhood (bp) within which a stronger DIFFERENT junction dominates/demotes this one.")
@click.option("--junction-dominance-tol-bp", default=2, show_default=True, type=int, help="(--junction-dominance-filter) bp tolerance treating an observed junction as the SAME as the candidate junction (vs a competing neighbor).")
@click.option(
    "--quant-mode",
    default="m2_em",
    type=click.Choice(["argmax", "m1_em", "m2_em"]),
    show_default=True,
    help="Quantification engine. 'm2_em' (production default): EM seeded by the PURE tie-break junction-NLL M2 distance (M1/AS selects each read's best-AS tie set + mappability mask; per-event junction NLL is the sole graded distance over that tie set), no M3 unless --m3-coherence. 'argmax': mappy AS argmax + M2 krill junction tiebreak, hard counts, no EM. 'm1_em': EM seeded by M1 mappy distance (beta=0). All signal scoring is in-memory krill.",
)
@click.option(
    "--m3-coherence/--no-m3-coherence",
    "m3_coherence",
    default=False,
    show_default=True,
    help="Add the M3 read×read junction-window krill DTW coherence to the EM quant modes (beta=em_beta). OFF by default because the pairwise DTW is expensive; the m2_em default runs the pure tie-break junction-NLL EM with no DTW. SIRV with-GTF Tx-F1 59.7 (M3 on, beta=1) vs 59.2 (off) — high-precision, recall-neutral; the payoff is the real-dRNA bet.",
)
@click.option("--abundance-feedback/--no-abundance-feedback", default=False, show_default=True, help="RSEM/Salmon-style iterative abundance feedback in the EM: each M-step re-estimates per-transcript abundance theta and biases a read shared between candidates toward the more abundant one (a 0.5/0.5 split migrating toward 0.99/0.01 once theta diverges, e.g. 100:1). Only affects the EM quant modes --quant-mode m1_em|m2_em (m2_em is the default; no effect under 'argmax'). Experimental; OFF by default.")
@click.option("--abundance-length-norm/--no-abundance-length-norm", default=False, show_default=True, help="With --abundance-feedback, divide abundance counts by per-transcript spliced effective length before forming theta (Salmon effective-length normalization). Experimental; OFF by default.")
@click.option("--threads", default=1, show_default=True, type=int, help="Interval-level worker processes. The pipeline is CPU-bound serial Python, so prefer many light workers; each worker is pinned to 1 BLAS/krill thread so total threads stay ~N. NOTE: the genome FASTA is loaded per worker (N× memory). 1 keeps the serial path.")
@click.option("--gpu-workers", default=0, show_default=True, type=int, help="Of --threads workers, how many hold a GPU context (live CUDA contexts <= G; VRAM bound = G× per-context). 0 = all workers CPU-only. Forced to 0 with --no-gpu. Over-provisioned GPU workers auto-fall back to CPU on OOM.")
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
    use_gpu,
    signal_format,
    fusion_enabled,
    min_support,
    max_dist,
    flank_bp,
    fusion_max_internal_gap_bp,
    max_reads_for_dtw,
    min_novel_reads,
    min_abundance,
    min_gtf_abundance,
    floor_gtf_abundance,
    min_isoform_fraction,
    max_soft_mass_ratio,
    min_fulllen_fraction,
    fulllen_window_bp,
    fulllen_min_reads,
    min_polya5p_reads,
    polya5p_window_bp,
    min_polya_length,
    polya5p_exempt_gtf,
    persist_R_matrix,
    canonical_gate,
    canonical_motifs,
    canonical_search_bp,
    m2_tiebreak,
    m2_tiebreak_margin,
    m2_tiebreak_junction_k,
    m2_diff_cover_gate,
    m2_diff_cover_margin,
    m2_cluster_recheck,
    m2_cluster_recheck_bp,
    m2_cluster_recheck_fraction,
    m2_cluster_recheck_cassette_max_exon_bp,
    m2_cluster_recheck_novel_displaces_gtf,
    m2_cluster_recheck_gtf_min_jct_reads,
    m2_cluster_recheck_jct_tol,
    m2_support_gate,
    m2_support_gate_tie,
    containment_collapse,
    containment_3p_tol_bp,
    containment_min_abundance_ratio,
    drop_mono_exon_novel,
    min_mono_exon_reads,
    min_mono_exon_length,
    novel_junction_min_reads,
    novel_junction_reads_tol,
    junction_dominance_filter,
    junction_dominance_min_reads,
    junction_dominance_window_bp,
    junction_dominance_tol_bp,
    quant_mode,
    m3_coherence,
    abundance_feedback,
    abundance_length_norm,
    threads,
    gpu_workers,
    verbose,
):
    """pyfin: nanopore signal-based transcriptome assembly.

    Performs reference-based transcriptome assembly with EM-based abundance/TPM
    quantification baked into the output. Pass --fusion to additionally detect
    gene fusions.
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

    if threads < 1:
        raise click.BadParameter("must be >= 1", param_hint="--threads")
    if gpu_workers < 0:
        raise click.BadParameter("must be >= 0", param_hint="--gpu-workers")
    if not use_gpu:
        gpu_workers = 0
    if gpu_workers > threads:
        raise click.BadParameter(
            f"must be <= --threads ({threads})", param_hint="--gpu-workers"
        )

    # Pin per-process thread pools to 1 BEFORE importing numpy (BLAS reads these at
    # import). Workers spawn-inherit os.environ, so total threads stay ~--threads.
    if threads > 1:
        for _var in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "KRILL_THREADS",
        ):
            os.environ.setdefault(_var, "1")

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
        signal_format=signal_format,
        fusion_enabled=fusion_enabled,
        fusion_min_support=min_support,
        fusion_max_dist=max_dist,
        fusion_flank_bp=flank_bp,
        fusion_max_internal_gap_bp=fusion_max_internal_gap_bp,
        max_reads_per_interval_for_dtw=max_reads_for_dtw,
        min_novel_reads=min_novel_reads,
        min_abundance=min_abundance,
        min_gtf_abundance=min_gtf_abundance,
        floor_gtf_abundance=floor_gtf_abundance,
        min_isoform_fraction=min_isoform_fraction,
        max_soft_mass_ratio=max_soft_mass_ratio,
        min_fulllen_fraction=min_fulllen_fraction,
        fulllen_window_bp=fulllen_window_bp,
        fulllen_min_reads=fulllen_min_reads,
        min_polya5p_reads=min_polya5p_reads,
        polya5p_window_bp=polya5p_window_bp,
        min_polya_length=min_polya_length,
        polya5p_exempt_gtf=polya5p_exempt_gtf,
        persist_R_matrix=persist_R_matrix,
        canonical_gate=canonical_gate,
        canonical_motifs=tuple(
            m.strip() for m in canonical_motifs.split(",") if m.strip()
        ),
        canonical_search_bp=canonical_search_bp,
        m2_tiebreak=m2_tiebreak,
        m2_tiebreak_margin=m2_tiebreak_margin,
        m2_tiebreak_junction_k=m2_tiebreak_junction_k,
        m2_diff_cover_gate=m2_diff_cover_gate,
        m2_diff_cover_margin=m2_diff_cover_margin,
        m2_cluster_recheck=m2_cluster_recheck,
        m2_cluster_recheck_bp=m2_cluster_recheck_bp,
        m2_cluster_recheck_fraction=m2_cluster_recheck_fraction,
        m2_cluster_recheck_cassette_max_exon_bp=m2_cluster_recheck_cassette_max_exon_bp,
        m2_cluster_recheck_novel_displaces_gtf=m2_cluster_recheck_novel_displaces_gtf,
        m2_cluster_recheck_gtf_min_jct_reads=m2_cluster_recheck_gtf_min_jct_reads,
        m2_cluster_recheck_jct_tol=m2_cluster_recheck_jct_tol,
        m2_support_gate=m2_support_gate,
        m2_support_gate_tie=m2_support_gate_tie,
        containment_collapse=containment_collapse,
        containment_3p_tol_bp=containment_3p_tol_bp,
        containment_min_abundance_ratio=containment_min_abundance_ratio,
        drop_mono_exon_novel=drop_mono_exon_novel,
        min_mono_exon_reads=min_mono_exon_reads,
        min_mono_exon_length=min_mono_exon_length,
        novel_junction_min_reads=novel_junction_min_reads,
        novel_junction_reads_tol=novel_junction_reads_tol,
        junction_dominance_filter=junction_dominance_filter,
        junction_dominance_min_reads=junction_dominance_min_reads,
        junction_dominance_window_bp=junction_dominance_window_bp,
        junction_dominance_tol_bp=junction_dominance_tol_bp,
        quant_mode=quant_mode,
        m3_coherence=m3_coherence,
        abundance_feedback=abundance_feedback,
        abundance_length_norm=abundance_length_norm,
        threads=threads,
        gpu_workers=gpu_workers,
    )

    runner = PipelineRunner(cfg)
    try:
        runner.setup()
        runner.run()
    finally:
        runner.cleanup()

    click.echo(f"Assembly output written to {output_dir}/")


if __name__ == "__main__":
    main()
