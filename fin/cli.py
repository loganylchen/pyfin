"""CLI entry point for pyfin."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import click


def _source_sha256(source_root: Path) -> str:
    """Hash the loaded Python source even when git is absent in a container."""
    digest = hashlib.sha256()
    for path in sorted((source_root / "fin").rglob("*.py")):
        digest.update(str(path.relative_to(source_root)).encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def _read_git_head(source_root: Path) -> str | None:
    """Resolve a loose Git HEAD without requiring the git executable."""
    git_dir = source_root / ".git"
    if git_dir.is_file():
        line = git_dir.read_text(errors="ignore").strip()
        if line.startswith("gitdir:"):
            git_dir = (source_root / line.split(":", 1)[1].strip()).resolve()
    head = git_dir / "HEAD"
    if not head.is_file():
        return None
    value = head.read_text(errors="ignore").strip()
    if not value.startswith("ref:"):
        return value or None
    ref = value.split(None, 1)[1]
    loose = git_dir / ref
    if loose.is_file():
        return loose.read_text(errors="ignore").strip() or None
    packed = git_dir / "packed-refs"
    if packed.is_file():
        for line in packed.read_text(errors="ignore").splitlines():
            if line and not line.startswith(("#", "^")):
                commit, name = line.split(" ", 1)
                if name == ref:
                    return commit
    return None


def _write_run_manifest(cfg, output_dir: str) -> None:
    """Persist the resolved scientific configuration before execution."""
    source_root = Path(__file__).resolve().parents[1]
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=source_root,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
        dirty = bool(subprocess.check_output(
            ["git", "status", "--porcelain", "--untracked-files=no"],
            cwd=source_root,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip())
    except (OSError, subprocess.CalledProcessError):
        commit = _read_git_head(source_root)
        dirty = None
    payload = {
        "schema_version": 1,
        "profile": cfg.profile,
        "profile_overrides": list(cfg.profile_overrides),
        "git_commit": commit,
        "git_dirty": dirty,
        "source_root": str(source_root),
        "source_sha256": _source_sha256(source_root),
        "result_environment": {
            key: os.environ.get(key)
            for key in (
                "KRILL_THREADS",
                "M1_MAX_INDEL_BP",
                "MAPPY_PRESET",
                "MAPPY_R1_MIN_AS",
                "PYFIN_CODE_ROOT",
            )
            if os.environ.get(key) is not None
        },
        "config": asdict(cfg),
    }
    path = Path(output_dir) / "run_manifest.json"
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


@click.group(invoke_without_command=True)
@click.option("--bam", default=None, type=click.Path(exists=True), help="Input BAM file.")
@click.option("--gtf", default=None, type=click.Path(), help="Reference GTF annotation file (optional).")
@click.option("--genome", default=None, type=click.Path(exists=True), help="Genome FASTA file.")
@click.option("--fastq", default=None, type=click.Path(exists=True), help="FASTQ reads file.")
@click.option("--signal", default=None, type=click.Path(exists=True), help="SLOW5/BLOW5/POD5 signal file.")
@click.option("--output-dir", default=None, help="Output directory.")
@click.option(
    "--profile",
    type=click.Choice(["real-drna", "real-drna-precision", "sirv", "custom"]),
    default="real-drna",
    show_default=True,
    help="Resolved scientific operating point. Explicit CLI options override profile values. Use 'sirv' for synthetic benchmarking, balanced 'real-drna' for biological samples, 'real-drna-precision' for the validated higher-F1/lower-recall support floor, or 'custom' for raw option defaults.",
)
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
@click.option("--min-abundance", default=3.0, show_default=True, type=float, help="NOVEL soft-abundance floor (GTF has --min-gtf-abundance; fusion exempt). Named profiles resolve SIRV=3 inclusive and real-dRNA=1 strict; explicit values override either profile.")
@click.option("--strict-novel-abundance-floor/--inclusive-novel-abundance-floor", "strict_novel_abundance_floor", default=False, show_default=True, help="Boundary semantics for --min-abundance. Strict requires abundance > floor; inclusive keeps abundance == floor. The real-dRNA profile uses strict >1 read-equivalent; SIRV uses inclusive >=3.")
@click.option("--min-gtf-abundance", default=1.0, show_default=True, type=float, help="Abundance floor for GTF-annotated transcripts on soft EM abundance (default 1: a GTF candidate must accrue at least one read of EM soft-mass to be kept, so pyfin is not a copy-annotation tool). Independent of --min-abundance (which gates novel); fusion always exempt. With --floor-gtf-abundance the effective GTF floor becomes max(this, --min-abundance). Set 0 to disable and keep every GTF candidate in an expressed locus. SIRV-tuned: on real data a nonzero floor can drop genuine low-abundance annotated isoforms.")
@click.option("--floor-gtf-abundance/--no-floor-gtf-abundance", "floor_gtf_abundance", default=False, show_default=True, help="Raise the GTF abundance floor to the NOVEL floor: GTF must clear max(--min-gtf-abundance, --min-abundance) (never lowers the explicit GTF floor). The real-dRNA/custom default is OFF; the optimized SIRV profile turns it ON for the expressed-truth operating point (e.g. T=3 -> Sn 43.2 / Pr 91.6). SIRV-tuned: on real data this drops genuine low-abundance annotated isoforms — leave OFF unless benchmarking. Fusion stays exempt regardless.")
@click.option("--min-isoform-fraction", default=0.01, show_default=True, type=float, help="Drop NOVEL multi-exon transcripts whose abundance is below this fraction of the dominant transcript in their selected locus (Cufflinks/StringTie minor-isoform suppression). GTF/fusion/mono exempt; 0.0 disables.")
@click.option("--isoform-fraction-locus", type=click.Choice(["family", "overlap"]), default="family", show_default=True, help="Denominator for --min-isoform-fraction. 'family' uses the persisted discovery splice family and falls back to overlap only for family-less candidates; 'overlap' restores the historical same-strand genomic-overlap behavior.")
@click.option("--post-selection-refit/--no-post-selection-refit", default=False, show_default=True, help="After the final survivor set and junction snapping are fixed, renormalize each assignable read over surviving candidates and report selection-orphaned mass. Named profiles enable this; custom/raw defaults leave historical model-filtered abundance unchanged.")
@click.option("--min-fulllen-fraction", default=0.1, show_default=True, type=float, help="Drop NOVEL multi-exon transcripts whose fraction of full-length assigned reads (read genomic 5' AND 3' both within --fulllen-window-bp of the candidate's ends) is below this (FLAIR/TALON-style full-length read support; signal-free). Orthogonal to --min-isoform-fraction. GTF/fusion/mono and unreachable candidates exempt; 0.0 disables. The SIRV profile resolves 0.1; real-dRNA resolves 0 because genuine 5'-truncated isoforms would otherwise lose recall.")
@click.option("--max-soft-mass-ratio", default=2.0, show_default=True, type=float, help="Drop NOVEL multi-exon transcripts whose EM soft abundance / hard argmax read count is >= this (0 hard reads always dropped). Real isoforms deposit ~1 soft mass per hard read (ratio ~1); wobble shadows borrow fractional soft crumbs from a high-abundance structural near-copy and show an inflated ratio — catches the high-relative-abundance shadows the fraction/cluster-recheck levers miss. Pure EM evidence (no GTF). GTF/fusion/mono exempt; 0.0 disables. SIRV-tuned default 2.0 (24-cell F1@3 up-or-equal, recall held); re-tune on real dRNA.")
@click.option("--fulllen-window-bp", default=25, show_default=True, type=int, help="bp tolerance for a read genomic end to count as full-length wrt a candidate's 5'/3' end (used by --min-fulllen-fraction).")
@click.option("--fulllen-min-reads", default=4, show_default=True, type=int, help="Minimum assigned reads carrying a genomic span required to score a candidate's full-length fraction; below this the candidate is unreachable and never dropped (used by --min-fulllen-fraction).")
@click.option("--min-polya5p-reads", default=1, show_default=True, type=int, help="Drop a candidate unless >= N of its assigned reads BOTH have a krill whole-read polyA tail (qc PASS & length > --min-polya-length) AND map with their genomic 5' end within --polya5p-window-bp of the candidate's 5' end. By default only novel candidates are gated (GTF and fusion exempt; pass --no-polya5p-exempt-gtf to gate GTF too). Needs --signal; adds a krill polyA pass to every run; 0 disables. Raw custom default is 1; both optimized named profiles resolve 0 after the current honest-F1 sweep. Explicitly set 1 for the historical gate.")
@click.option("--polya5p-window-bp", default=25, show_default=True, type=int, help="bp tolerance for a read's genomic 5' end to count toward --min-polya5p-reads.")
@click.option("--min-polya-length", default=10.0, show_default=True, type=float, help="Minimum krill polya_length (with polya_qc PASS) for a read to support a candidate under --min-polya5p-reads.")
@click.option("--polya5p-exempt-gtf/--no-polya5p-exempt-gtf", "polya5p_exempt_gtf", default=True, show_default=True, help="Exempt GTF-sourced candidates from the --min-polya5p-reads filter (like fusion); only novel candidates are gated. Default ON: validated as the corrected-recall-optimal config on SIRV (dRNA reads are often 5'-truncated/lack polyA, so gating GTF drops genuine annotated transcripts). Pass --no-polya5p-exempt-gtf to gate GTF candidates too.")
@click.option("--persist-R/--no-persist-R", "persist_R_matrix", default=True, show_default=True, help="Enable/disable R-matrix (R.npy) persistence per interval.")
@click.option("--canonical-gate/--no-canonical-gate", "canonical_gate", default=True, show_default=True, help="Drop NOVEL multi-exon candidates whose junctions aren't all canonical (GTF/fusion/mono exempt). SIRV-tuned default ON.")
@click.option("--canonical-motifs", default="GT-AG,GC-AG,AT-AC", show_default=True, help="Comma-separated donor-acceptor motifs accepted by the canonical gate AND search.")
@click.option("--chain-cluster-discovery/--no-chain-cluster-discovery", "chain_cluster_discovery", default=True, show_default=True, help="PRODUCTION DEFAULT generation. Group reads by intron chain (3' IGNORED), collapse EXACT sub-chains, cluster wobble/cassette/containment keeping members; NO canonical-search shadows. Replaces the legacy 3'+exact-chain grouping + canonical expansion. p00 no-gates: raw candidates 175k->37k, c -73%, j -88%, structural precision 10->26%, at ~20% distinct-truth recall cost (92% from the exact sub-chain collapse). --no- reverts to the legacy discovery.")
@click.option("--clustering", type=click.Choice(["read_chains", "families"]), default="families", show_default=True, help="(--chain-cluster-discovery) generation clustering primitive. 'families' = clustering redesign (PRODUCTION DEFAULT): cluster_families (grouping only) + collapse (span-guarded exact-subchain fold); mono deferred to a bucket resolved post-EM. Emits the same multi-exon candidate set as the legacy path (EM cluster scoping may differ in rare subchain-bridge cases; mono handling differs by design) and beats the old default end-to-end on real human p00 (recall +0.4, Pr +1.9). 'read_chains' = legacy cluster_read_chains (fold+cluster inline); add --no-mono-resolve-post-em for the exact old behavior.")
@click.option("--chain-cluster-wobble-bp", default=6, show_default=True, type=int, help="(--chain-cluster-discovery) per-junction bp tolerance for the wobble/cassette/containment cluster joins.")
@click.option("--chain-cluster-cassette-max-exon-bp", default=70, show_default=True, type=int, help="(--chain-cluster-discovery) max skipped-exon length (bp) for a cassette (K vs K-1) cluster join. 0 disables cassette joins.")
@click.option("--chain-cluster-fold-monoexon/--no-chain-cluster-fold-monoexon", "chain_cluster_fold_monoexon", default=True, show_default=True, help="(--chain-cluster-discovery) Generation-time mono fold. ONLY takes effect when post-EM mono resolution is NOT running (i.e. --no-mono-resolve-post-em, or a non-m2_em quant mode); under the production default (--mono-resolve-post-em) it is superseded by the post-EM resolution and this flag is inert. Fold a single-exon read whose aligned span lies wholly inside one exon of a multi-exon candidate INTO that candidate (a 5'/3' degradation fragment) instead of emitting a standalone single-exon candidate. Reads inside an intron (a different gene) or contained in no multi-exon candidate stay separate. Suppresses single-exon truncation FPs (p00 m2_em: out -400, `=` +91, Pr 66.0->69.9). --no- reverts to standalone mono candidates.")
@click.option("--chain-cluster-fold-span-guard/--no-chain-cluster-fold-span-guard", "chain_cluster_fold_span_guard", default=True, show_default=True, help="(--chain-cluster-discovery) PRODUCTION DEFAULT. Only fold a read into an exact-sub-chain container if its aligned span does NOT run exonically across one of the container's EXTRA introns. A read that spans such an intron (retained-intron / alternative isoform) contradicts it and is kept as its own candidate instead of being absorbed. --no- reverts to the unconditional exact-sub-chain fold.")
@click.option("--mono-resolve-post-em/--no-mono-resolve-post-em", "mono_resolve_post_em", default=True, show_default=True, help="(quant_mode=m2_em) PRODUCTION DEFAULT (pairs with --clustering families). Defer single-exon read resolution to AFTER the multi-exon EM: generation stops folding mono reads into multi (chain-cluster-fold-monoexon forced off); post-EM each surviving mono candidate's reads are re-resolved against SURVIVING multi candidates by strict strand-aware exonic containment (1 cover -> fold; several -> highest-EM-abundance; uncovered kept, mono dropped if < --mono-resolve-min-reads). On real human p00 this beats the old generation-time fold_monoexon (recall +0.4, Pr +1.9). --no- disables post-EM resolution, falling back to the generation-time fold (--chain-cluster-fold-monoexon, on by default) -- or standalone mono if that is also off.")
@click.option("--mono-resolve-min-reads", default=2, show_default=True, type=int, help="(--mono-resolve-post-em) uncovered reads a mono candidate must retain to survive as a single-exon transcript.")
@click.option("--mono-resolve-slop-bp", default=10, show_default=True, type=int, help="(--mono-resolve-post-em) terminal boundary slop (bp) for the mono read exonic-containment test.")
@click.option("--canonical-search-bp", default=4, show_default=True, type=int, help="ea extended search: scan ±N bp around each read-derived NOVEL junction for canonical motifs and emit paired alternatives (GTF transcripts not extended). 0 disables. SIRV-tuned default 4. NOTE: ignored when --chain-cluster-discovery is on.")
@click.option("--m2-tiebreak/--no-m2-tiebreak", "m2_tiebreak", default=True, show_default=True, help="Resolve argmax_keep ties with the junction-window mean-NLL signal metric: give a read's full mass to the M2-best tied candidate when the NLL margin >= --m2-tiebreak-margin, else keep the 1/K split. Needs --signal (auto-skips if absent). Default ON + aggressive (margin 1e-9): SIRV gffcompare Tx-F1 45.4 vs OFF 44.7, beats all tools + ablation champion 45.2.")
@click.option("--m2-tiebreak-margin", default=1e-9, show_default=True, type=float, help="Minimum M2 NLL margin (runner-up - best) required to override the 1/K split with the M2 single winner. Default 1e-9 (aggressive: take M2's pick whenever it can discriminate at all).")
@click.option("--m2-tiebreak-junction-k", default=10, show_default=True, type=int, help="Transcript-frame bp on each side of the wobbling junction for the M2 discrimination window (SIRV sweet spot 10).")
@click.option(
    "--m2-metric",
    type=click.Choice(["off", "mean", "summed_llr", "auto"]),
    default="mean",
    show_default=True,
    help="(--quant-mode m2_em) Signal refinement for exact M1 ties. 'mean' is the legacy wide-window mean NLL; 'summed_llr' uses tight differing-junction windows and an undivided NLL sum; 'off' keeps the M1 tie; 'auto' selects mean with a usable GTF and summed_llr when unguided. Profile-controlled unless explicitly set.",
)
@click.option("--m2-summed-llr-margin", default=2.0, show_default=True, type=float, help="(--m2-metric summed_llr) Minimum summed-NLL runner-up gap for a hard winner. Sum has a different scale from the legacy mean margin.")
@click.option("--m2-summed-llr-flank", default=6, show_default=True, type=int, help="(--m2-metric summed_llr) Genomic bp flank around each differing junction boundary in the tight scoring window.")
@click.option("--m2-diff-cover-gate/--no-m2-diff-cover-gate", "m2_diff_cover_gate", default=True, show_default=True, help="(--quant-mode m2_em only) Diff-region coverage gate. For a read in a >=2 best-AS tie: if the configured M2 margin reaches its metric-specific threshold, hard-assign its full mass to the best candidate (and, if it straddles every wobbling junction donor->acceptor, contribute that vote to the locus isoform-ratio prior); otherwise the read is ambiguous and follows the covered-read prior (flat 1/K only when no prior exists). A valid contrast requires two scored hypotheses. No read is dropped. Profile-controlled.")
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
@click.option("--containment-cluster/--no-containment-cluster", "containment_cluster", default=True, show_default=True, help="(--quant-mode m2_em only, DEFAULT ON) Recall-SAFER generalisation of --containment-collapse (runs AFTER all structural/support gates, so a shadow is never folded into an already-dropped parent). Drop a NOVEL candidate whose intron chain is a contiguous SUB-CHAIN (within --containment-cluster-wobble-bp per junction) of a longer candidate — a truncation/exon-skip shadow the same-intron-count wobble cluster never groups — ONLY when it is a low-support shadow by BOTH EM abundance (<= parent * --containment-cluster-min-ab-ratio) AND supporting-read count (<= parent * --containment-cluster-min-read-ratio). The read-support guard keeps MOST genuine short/alt-TSS isoforms (on p00 a truncation shadow carries ~1 read vs a real short isoform's ~13), but a genuine low-fraction minor isoform can still be dropped — hence recall-safER not recall-safe. gtf/fusion never dropped. DEFAULT ON; --no-containment-cluster disables.")
@click.option("--containment-cluster-wobble-bp", default=6, show_default=True, type=int, help="(--containment-cluster) per-junction bp tolerance matching the sub-chain to the parent's window.")
@click.option("--containment-cluster-min-ab-ratio", default=0.3, show_default=True, type=float, help="(--containment-cluster) drop only when shadow EM abundance <= parent's * this ratio.")
@click.option("--containment-cluster-min-read-ratio", default=0.3, show_default=True, type=float, help="(--containment-cluster) drop only when shadow supporting-read count <= parent's * this ratio (the recall-safER guard: keeps most genuine short isoforms, but a real low-fraction minor isoform can still fall below the ratio).")
@click.option("--containment-cluster-max-shadow-reads", default=10, show_default=True, type=int, help="(--containment-cluster) absolute cap: never drop a shadow carrying more than N supporting reads regardless of the ratio (protects a genuine low-fraction isoform of a very-high-support parent). 0 disables (ratio-only).")
@click.option("--drop-mono-exon-novel/--no-drop-mono-exon-novel", "drop_mono_exon_novel", default=False, show_default=True, help="Lever 3: drop NOVEL single-exon (mono) candidates with weak support — hard read count < --min-mono-exon-reads OR genomic length < --min-mono-exon-length. Suppresses single-exon de novo noise (like IsoQuant's ONT default) WITHOUT a blanket drop: a high-support/long real intronless gene survives. gtf/fusion/multi-exon exempt. Needs at least one threshold > 0 to fire. Default OFF.")
@click.option("--min-mono-exon-reads", default=0, show_default=True, type=int, help="(--drop-mono-exon-novel) Min hard reads for a novel mono candidate to survive. 0 disables this threshold.")
@click.option("--min-mono-exon-length", default=0, show_default=True, type=int, help="(--drop-mono-exon-novel) Min genomic length (bp) for a novel mono candidate to survive. 0 disables this threshold.")
@click.option("--junction-snap/--no-junction-snap", "junction_snap", default=False, show_default=True, help="After global selection, correct NOVEL multi-exon junctions to a more strongly supported nearby primary-read CIGAR mode and merge structurally identical corrected models while preserving read/abundance mass. GTF/fusion exempt. Default off pending live validation.")
@click.option("--junction-snap-tolerance", default=6, show_default=True, type=int, help="(--junction-snap) Maximum donor and acceptor shift in bp.")
@click.option("--junction-snap-min-support", default=2, show_default=True, type=int, help="(--junction-snap) Minimum direct read support for the target junction mode.")
@click.option("--junction-snap-min-ratio", default=2.0, show_default=True, type=float, help="(--junction-snap) Target junction support must exceed current support times this ratio.")
@click.option("--novel-junction-min-reads", default=2, show_default=True, type=int, help="(--quant-mode m2_em only) Lever 2: drop a NOVEL multi-exon candidate if ANY of its junctions is spliced by fewer than N directly-observed reads (primary-read CIGAR introns, strand-keyed, within --novel-junction-reads-tol bp). I.e. a novel junction must be carried by >= N reads, not just 1. gtf/fusion/mono exempt. <=1 disables. DEFAULT 2 (production): on real human GENCODE it lifts precision at zero recall cost across every scenario (de novo Pr +3.8, ratios +0.8..+3.2, all corrupted-GTF +0.6..+0.9) and keeps pyfin #1 on SIRV4. Set 0 to disable.")
@click.option("--novel-junction-reads-tol", default=2, show_default=True, type=int, help="(--novel-junction-min-reads) bp tolerance matching a candidate junction to an observed read junction.")
@click.option("--guided-junction-min-reads", default=0, show_default=True, type=int, help="(--quant-mode m2_em only, EXPERIMENT) GUIDED junction-support gate: mirror of Lever 2 for GTF-passthrough candidates — drop a guided multi-exon candidate if ANY of its junctions is spliced by fewer than N directly-observed reads (primary-read CIGAR introns, within --guided-junction-reads-tol bp). Extends the coordinate-EXACT read-support check (Lever 2 is novel-only) to GTF junctions, so a jitter-corrupted annotation junction (coords > tol from the true site) is dropped instead of surviving the coordinate-inexact M2/M1 gate. Targets the c_jitter precision collapse, and trims GTF echo in loci that DO show observed splicing (a FULLY read-sparse locus fails open — echo there is unaffected). NOT recall-safe (also drops low-coverage annotated junctions). Default 0 (OFF) pending real-data validation.")
@click.option("--guided-junction-reads-tol", default=2, show_default=True, type=int, help="(--guided-junction-min-reads) bp tolerance matching a guided candidate junction to an observed read junction.")
@click.option("--denovo-wobble-tol", default=0, show_default=True, type=int, help="(EXPERIMENT) De-novo wobble-tolerant collapse: merge NOVEL candidates whose intron chains match within N bp per junction (and 3' within --three-prime-threshold) into the highest-read-support consensus, BEFORE they become separate candidates. Attacks pyfin's low de-novo structural precision (wobble shadows: a junction and its ±few-bp minimap variant currently survive as separate novel transcripts). Consensus = the reads' own mode, NOT annotation (no snap). 0 disables (byte-identical). isoquant uses ~6 for ONT.")
@click.option("--denovo-wobble-shadow-ratio", default=0.5, show_default=True, type=float, help="(--denovo-wobble-tol) Absorb a wobble-matching candidate into a rep only if it is a true SHADOW: len(cand.reads) <= ratio * len(rep.reads). Protects genuine close isoforms (real alt donor/acceptor a few bp apart) from being merged. 1.0 = merge any wobble-match. isoquant's ½-support bulge criterion ≈ 0.5.")
@click.option("--denovo-graph/--no-denovo-graph", "denovo_graph", default=False, show_default=True, help="(EXPERIMENT) De-novo intron-graph assembly: pool all reads' junctions, cluster wobbles to a read-count consensus (--denovo-graph-tol), build a read-adjacency graph, and EXTEND each read's chain through UNAMBIGUOUSLY-supported edges (>= --denovo-graph-min-edge-reads) into a maximal full-length chain — assembling truncated 3'-biased dRNA reads into full transcripts (attacks the measured #1 de-novo error: 'c' contained/truncated). Stops at genuine branch points (no fabricated combos). Reads grouped by extended chain; candidate span = union of its reads' extents, sequence from genome. NOT annotation-snap. Default OFF (byte-identical).")
@click.option("--denovo-graph-tol", default=6, show_default=True, type=int, help="(--denovo-graph) bp tolerance clustering wobbled junctions to a read-count consensus before building the graph. isoquant uses ~6 for ONT.")
@click.option("--denovo-graph-min-edge-reads", default=2, show_default=True, type=int, help="(--denovo-graph) Minimum reads supporting a junction->junction adjacency edge for it to be used when extending a chain. Higher = more conservative assembly.")
@click.option("--denovo-graph-tss-brake/--no-denovo-graph-tss-brake", "denovo_graph_tss_brake", default=True, show_default=True, help="(--denovo-graph) 5'-TSS brake: don't extend a chain 5'-ward past a real transcription start. dRNA truncation is 5'-ward, so a truncated read of a LONG transcript stops at a random 5' position while a COMPLETE read of a genuine short isoform stops at its real TSS. If >= --denovo-graph-tss-frac of the reads starting at a chain's 5' junction pile their 5'-end within --denovo-graph-tss-tol bp (and >= --denovo-graph-tss-min-reads), that's a real TSS and the short isoform is kept instead of merged into the longer one. Validated on p00: real short-isoform TSS carry 40-90% of read-5'-ends; degradation scatters. --no- disables (pure assembly, over-merges).")
@click.option("--denovo-graph-tss-tol", default=20, show_default=True, type=int, help="(--denovo-graph-tss-brake) bp window for clustering read-5'-ends into a TSS peak.")
@click.option("--denovo-graph-tss-min-reads", default=3, show_default=True, type=int, help="(--denovo-graph-tss-brake) Minimum reads piled at a 5' position to call it a real TSS.")
@click.option("--denovo-graph-tss-frac", default=0.4, show_default=True, type=float, help="(--denovo-graph-tss-brake) Minimum fraction of a chain's 5'-starting reads that must pile at one 5' position to call it a real TSS (peak vs degradation background).")
@click.option("--junction-dominance-filter/--no-junction-dominance-filter", "junction_dominance_filter", default=False, show_default=True, help="(--quant-mode m2_em only) Junction-first PRE-EM gate: drop a NOVEL multi-exon candidate if ANY junction has < --junction-dominance-min-reads observed reads OR is not locally dominant (a different observed junction within --junction-dominance-window-bp carries strictly more reads). Removes multi-read wobble shadows (they lose to the stronger true junction nearby) BEFORE EM, so shadows never compete for reads. Pure mapping evidence, no snap. gtf/fusion/mono exempt. Default OFF.")
@click.option("--junction-dominance-min-reads", default=2, show_default=True, type=int, help="(--junction-dominance-filter) Min observed reads for a novel junction to survive.")
@click.option("--junction-dominance-window-bp", default=20, show_default=True, type=int, help="(--junction-dominance-filter) Neighborhood (bp) within which a stronger DIFFERENT junction dominates/demotes this one.")
@click.option("--junction-dominance-tol-bp", default=2, show_default=True, type=int, help="(--junction-dominance-filter) bp tolerance treating an observed junction as the SAME as the candidate junction (vs a competing neighbor).")
@click.option(
    "--quant-mode",
    default="m2_em",
    type=click.Choice(["argmax", "m1_em", "m2_em", "cluster"]),
    show_default=True,
    help="Quantification engine. 'm2_em' (production default): EM seeded by the PURE tie-break junction-NLL M2 distance (M1/AS selects each read's best-AS tie set + mappability mask; per-event junction NLL is the sole graded distance over that tie set). 'argmax': mappy AS argmax + M2 krill junction tiebreak, hard counts, no EM. 'm1_em': EM seeded by M1 mappy distance (beta=0). 'cluster': assign reads to candidates WITHIN each generation cluster only (CandidateSet.clusters), M1 within-cluster + summed-LLR M2 tiebreak + main-peak containment (fin.pipeline.cluster_quant). All signal scoring is in-memory krill.",
)
@click.option(
    "--cluster-use-m2/--no-cluster-use-m2", "cluster_use_m2", default=True,
    show_default=True,
    help="(--quant-mode cluster) Run the summed-LLR M2 signal tiebreak on straddling M1 ties. With a wide --cluster-m1-tie-margin this fires on many ties (expensive eventalign) for little end-to-end gain; --no-cluster-use-m2 sends near-ties straight to the ambiguous 1/K -> EM path (much faster). Signal is still used by the polyA finalize gate.")
@click.option(
    "--cluster-m1-tie-margin", "cluster_m1_tie_margin", default=20.0, show_default=True, type=float,
    help="(--quant-mode cluster) AS margin for the within-cluster M1 best tie: a read is unique-best only when its top member beats the runner-up by MORE than this many AS points; within it the members tie (-> ambiguous 1/K -> EM). A 2-3bp wobble shifts AS by <=~20. Larger => more wobble near-ties merged by EM (fewer wobble FPs but risks merging the exact-coordinate variant); smaller => keep wobble variants separate.")
@click.option("--abundance-feedback/--no-abundance-feedback", default=False, show_default=True, help="RSEM/Salmon-style iterative abundance feedback in the EM: each M-step re-estimates per-transcript abundance theta and biases a read shared between candidates toward the more abundant one (a 0.5/0.5 split migrating toward 0.99/0.01 once theta diverges, e.g. 100:1). Only affects the EM quant modes --quant-mode m1_em|m2_em (m2_em is the default; no effect under 'argmax'). Experimental; OFF by default.")
@click.option("--abundance-length-norm/--no-abundance-length-norm", default=False, show_default=True, help="With --abundance-feedback, divide abundance counts by per-transcript spliced effective length before forming theta (Salmon effective-length normalization). Experimental; OFF by default.")
@click.option("--threads", default=1, show_default=True, type=int, help="Interval-level worker processes. The pipeline is CPU-bound serial Python, so prefer many light workers; each worker is pinned to 1 BLAS/krill thread so total threads stay ~N. NOTE: the genome FASTA is loaded per worker (N× memory). 1 keeps the serial path.")
@click.option("--gpu-workers", default=0, show_default=True, type=int, help="Of --threads workers, how many hold a GPU context (live CUDA contexts <= G; VRAM bound = G× per-context). 0 = all workers CPU-only. Forced to 0 with --no-gpu. Over-provisioned GPU workers auto-fall back to CPU on OOM.")
@click.option("--write-unfiltered-scores/--no-write-unfiltered-scores", "write_unfiltered_scores", default=False, show_default=True, help="Diagnostic: also write scores.unfiltered.tsv (post-EM, PRE post-EM-filter snapshot of every candidate) alongside scores.tsv. Lets FN root-cause analysis split 'candidate dropped by a post-EM filter' (present in unfiltered, absent in final) from 'never reached EM / not generated' (absent in both). No effect on the emitted GTF. OFF by default.")
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
    profile,
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
    strict_novel_abundance_floor,
    min_gtf_abundance,
    floor_gtf_abundance,
    min_isoform_fraction,
    isoform_fraction_locus,
    post_selection_refit,
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
    chain_cluster_discovery,
    clustering,
    chain_cluster_wobble_bp,
    chain_cluster_cassette_max_exon_bp,
    chain_cluster_fold_monoexon,
    chain_cluster_fold_span_guard,
    mono_resolve_post_em,
    mono_resolve_min_reads,
    mono_resolve_slop_bp,
    canonical_search_bp,
    m2_tiebreak,
    m2_tiebreak_margin,
    m2_tiebreak_junction_k,
    m2_metric,
    m2_summed_llr_margin,
    m2_summed_llr_flank,
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
    containment_cluster,
    containment_cluster_wobble_bp,
    containment_cluster_min_ab_ratio,
    containment_cluster_min_read_ratio,
    containment_cluster_max_shadow_reads,
    drop_mono_exon_novel,
    min_mono_exon_reads,
    min_mono_exon_length,
    junction_snap,
    junction_snap_tolerance,
    junction_snap_min_support,
    junction_snap_min_ratio,
    novel_junction_min_reads,
    novel_junction_reads_tol,
    guided_junction_min_reads,
    guided_junction_reads_tol,
    denovo_wobble_tol,
    denovo_wobble_shadow_ratio,
    denovo_graph,
    denovo_graph_tol,
    denovo_graph_min_edge_reads,
    denovo_graph_tss_brake,
    denovo_graph_tss_tol,
    denovo_graph_tss_min_reads,
    denovo_graph_tss_frac,
    junction_dominance_filter,
    junction_dominance_min_reads,
    junction_dominance_window_bp,
    junction_dominance_tol_bp,
    quant_mode,
    cluster_use_m2,
    cluster_m1_tie_margin,
    abundance_feedback,
    abundance_length_norm,
    threads,
    gpu_workers,
    write_unfiltered_scores,
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

    from click.core import ParameterSource

    from fin.pipeline.config import (
        PROFILE_FIELDS,
        PipelineConfig,
        resolve_m2_metric,
        resolve_profile_values,
    )
    from fin.pipeline.runner import PipelineRunner

    explicit_profile_fields = {
        name
        for name in PROFILE_FIELDS
        if name in ctx.params
        and ctx.get_parameter_source(name) is not ParameterSource.DEFAULT
    }
    resolved = resolve_profile_values(
        profile, ctx.params, explicit_fields=explicit_profile_fields
    )
    min_novel_reads = resolved["min_novel_reads"]
    min_abundance = resolved["min_abundance"]
    strict_novel_abundance_floor = resolved["strict_novel_abundance_floor"]
    min_gtf_abundance = resolved["min_gtf_abundance"]
    floor_gtf_abundance = resolved["floor_gtf_abundance"]
    min_isoform_fraction = resolved["min_isoform_fraction"]
    post_selection_refit = resolved["post_selection_refit"]
    max_soft_mass_ratio = resolved["max_soft_mass_ratio"]
    min_fulllen_fraction = resolved["min_fulllen_fraction"]
    min_polya5p_reads = resolved["min_polya5p_reads"]
    canonical_gate = resolved["canonical_gate"]
    novel_junction_min_reads = resolved["novel_junction_min_reads"]
    containment_cluster = resolved["containment_cluster"]
    drop_mono_exon_novel = resolved["drop_mono_exon_novel"]
    min_mono_exon_reads = resolved["min_mono_exon_reads"]
    junction_snap = resolved["junction_snap"]
    junction_snap_tolerance = resolved["junction_snap_tolerance"]
    junction_snap_min_support = resolved["junction_snap_min_support"]
    junction_snap_min_ratio = resolved["junction_snap_min_ratio"]
    m2_cluster_recheck = resolved["m2_cluster_recheck"]
    m2_metric, m2_metric_route = resolve_m2_metric(
        resolved["m2_metric"], gtf
    )
    m2_diff_cover_margin = resolved["m2_diff_cover_margin"]
    m2_summed_llr_margin = resolved["m2_summed_llr_margin"]
    m2_summed_llr_flank = resolved["m2_summed_llr_flank"]

    os.makedirs(output_dir, exist_ok=True)

    cfg = PipelineConfig(
        bam_path=bam,
        profile=profile,
        profile_overrides=tuple(sorted(explicit_profile_fields)),
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
        strict_novel_abundance_floor=strict_novel_abundance_floor,
        min_gtf_abundance=min_gtf_abundance,
        floor_gtf_abundance=floor_gtf_abundance,
        min_isoform_fraction=min_isoform_fraction,
        isoform_fraction_locus=isoform_fraction_locus,
        post_selection_refit=post_selection_refit,
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
        chain_cluster_discovery=chain_cluster_discovery,
        clustering=clustering,
        chain_cluster_wobble_bp=chain_cluster_wobble_bp,
        chain_cluster_cassette_max_exon_bp=chain_cluster_cassette_max_exon_bp,
        chain_cluster_fold_monoexon=chain_cluster_fold_monoexon,
        chain_cluster_fold_span_guard=chain_cluster_fold_span_guard,
        mono_resolve_post_em=mono_resolve_post_em,
        mono_resolve_min_reads=mono_resolve_min_reads,
        mono_resolve_slop_bp=mono_resolve_slop_bp,
        canonical_search_bp=canonical_search_bp,
        m2_tiebreak=m2_tiebreak,
        m2_tiebreak_margin=m2_tiebreak_margin,
        m2_tiebreak_junction_k=m2_tiebreak_junction_k,
        m2_metric=m2_metric,
        m2_metric_route=m2_metric_route,
        m2_summed_llr_margin=m2_summed_llr_margin,
        m2_summed_llr_flank=m2_summed_llr_flank,
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
        containment_cluster=containment_cluster,
        containment_cluster_wobble_bp=containment_cluster_wobble_bp,
        containment_cluster_min_ab_ratio=containment_cluster_min_ab_ratio,
        containment_cluster_min_read_ratio=containment_cluster_min_read_ratio,
        containment_cluster_max_shadow_reads=containment_cluster_max_shadow_reads,
        drop_mono_exon_novel=drop_mono_exon_novel,
        min_mono_exon_reads=min_mono_exon_reads,
        min_mono_exon_length=min_mono_exon_length,
        junction_snap=junction_snap,
        junction_snap_tolerance=junction_snap_tolerance,
        junction_snap_min_support=junction_snap_min_support,
        junction_snap_min_ratio=junction_snap_min_ratio,
        novel_junction_min_reads=novel_junction_min_reads,
        novel_junction_reads_tol=novel_junction_reads_tol,
        guided_junction_min_reads=guided_junction_min_reads,
        guided_junction_reads_tol=guided_junction_reads_tol,
        denovo_wobble_tol=denovo_wobble_tol,
        denovo_wobble_shadow_ratio=denovo_wobble_shadow_ratio,
        denovo_graph=denovo_graph,
        denovo_graph_tol=denovo_graph_tol,
        denovo_graph_min_edge_reads=denovo_graph_min_edge_reads,
        denovo_graph_tss_brake=denovo_graph_tss_brake,
        denovo_graph_tss_tol=denovo_graph_tss_tol,
        denovo_graph_tss_min_reads=denovo_graph_tss_min_reads,
        denovo_graph_tss_frac=denovo_graph_tss_frac,
        junction_dominance_filter=junction_dominance_filter,
        junction_dominance_min_reads=junction_dominance_min_reads,
        junction_dominance_window_bp=junction_dominance_window_bp,
        junction_dominance_tol_bp=junction_dominance_tol_bp,
        quant_mode=quant_mode,
        cluster_use_m2=cluster_use_m2,
        cluster_m1_tie_margin=cluster_m1_tie_margin,
        abundance_feedback=abundance_feedback,
        abundance_length_norm=abundance_length_norm,
        threads=threads,
        gpu_workers=gpu_workers,
        write_unfiltered_scores=write_unfiltered_scores,
    )

    try:
        cfg.validate()
    except (FileNotFoundError, ValueError) as exc:
        raise click.UsageError(str(exc)) from exc

    import logging

    logging.getLogger(__name__).info(
        "Resolved profile=%s overrides=%s m2_metric=%s(%s) abundance=%s%s/%s refit=%s full=%s polyA=%s",
        cfg.profile,
        ",".join(cfg.profile_overrides) or "none",
        cfg.m2_metric,
        cfg.m2_metric_route,
        ">" if cfg.strict_novel_abundance_floor else ">=",
        cfg.min_abundance,
        cfg.min_gtf_abundance,
        cfg.post_selection_refit_effective,
        cfg.min_fulllen_fraction,
        cfg.min_polya5p_reads,
    )
    _write_run_manifest(cfg, output_dir)

    runner = PipelineRunner(cfg)
    try:
        runner.setup()
        runner.run()
    finally:
        runner.cleanup()

    click.echo(f"Assembly output written to {output_dir}/")


if __name__ == "__main__":
    main()
