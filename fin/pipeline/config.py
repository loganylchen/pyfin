"""Pipeline configuration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Optional


@dataclass
class PipelineConfig:
    """Configuration for the full pyfin pipeline."""

    # Input files
    bam_path: str
    gtf_path: Optional[str] = None
    genome_fasta_path: str = ""
    fastq_path: str = ""
    signal_path: str = ""  # SLOW5/BLOW5 file
    signal_format: str = "slow5"  # "pod5" or "slow5"

    # Working directory
    work_dir: str = "./pyfin_work"

    # Candidate discovery
    three_prime_threshold: int = 24
    max_gap: int = 0
    min_novel_reads: int = 1          # A2: min supporting reads for a novel candidate (after collapsing)
    min_abundance: float = 0.0        # A4: drop quantified transcripts with abundance < this threshold
    min_max_r: float = 0.0            # Phase A T2: drop NOVEL transcripts whose max EM responsibility < this (Cohen's d=1.34)
    min_novel_combined_score: float = 0.0  # Step 3: drop NOVEL transcripts whose combined_score < this (F1-optimal: 0.288 with-GTF, 0.428 no-GTF; Cohen's d=0.70)

    # Minimum isoform fraction (locus-relative abundance) FILTER. Drop a NOVEL
    # multi-exon transcript whose abundance is below this fraction of the
    # dominant OVERLAPPING novel isoform at its locus. This is the standard
    # Cufflinks --min-isoform-fraction (-F, default 0.10) / StringTie -f
    # (default 0.01) minor-isoform suppression heuristic: the low-fraction tail
    # of a locus is enriched for incompletely-spliced pre-mRNA, RT/template-
    # switching artifacts, and assembly noise rather than genuine isoforms.
    # Applied as a post-quant candidate FILTER (never inside EM/assignment);
    # GTF-passthrough, fusion, and single-exon candidates are EXEMPT. Default
    # 0.01 (StringTie-aligned, recall-safe for long reads). Set 0.0 to disable.
    # SIRV WARNING: synthetic SIRV's F1-optimal value (~0.4) is 4-40x more
    # aggressive than any real tool because SIRV lacks a real low-abundance
    # isoform tail — NEVER use it on real data. Configurable; re-tune on real
    # transcriptomes.
    min_isoform_fraction: float = 0.01

    # Full-length end-coherence FILTER (FLAIR/TALON-style full-length read
    # support). Drop a NOVEL multi-exon (>=2 intron) transcript whose fraction
    # of full-length assigned reads is below `min_fulllen_fraction`. A read is
    # full-length wrt candidate C when its genomic 5' AND 3' alignment ends are
    # BOTH within `fulllen_window_bp` of C's genomic 5'/3' ends (strand-aware).
    # Signal-free: uses BAM primary-alignment spans only (no f5c/f5c_rna).
    # Applied as a post-quant candidate FILTER (never inside EM/assignment),
    # over the argmax assignment population where it is non-circular (discovery
    # groups support reads by 3' end, so it must NOT be applied at candidate
    # generation). It is ORTHOGONAL to min_isoform_fraction and stacks with it.
    # GTF-passthrough, fusion, and single-exon candidates are EXEMPT. Candidates
    # with fewer than `fulllen_min_reads` assigned reads carrying a genomic span
    # (unreachable), and candidates on the legacy EM path (fulllen never
    # computed), keep the -1.0 sentinel and are NEVER dropped. Set 0.0 to
    # disable. SIRV WARNING: the default 0.1 is SIRV-tuned (it drops ~74% of
    # reachable novel-multi candidates "for free" because synthetic SIRV lacks a
    # real 5'-truncated minor-isoform tail). On a real dRNA transcriptome
    # (sequenced 3'->5', genuine 5'-truncated isoforms) this same value would
    # cost recall — re-tune or disable on real data. Configurable; never
    # hardcode a SIRV-tuned threshold for production real-data runs.
    min_fulllen_fraction: float = 0.1
    fulllen_window_bp: int = 25          # bp tolerance for a read end to count as full-length
    fulllen_min_reads: int = 4           # min assigned reads-with-span to score (else unreachable)

    # Canonical-motif alternative expansion: scan ±N bp around each read's
    # CIGAR-derived junction for GT-AG and emit additional chains.
    # Stage C: ea extended canonical SEARCH. For each read-derived novel
    # junction, scan ±N bp for canonical (donor,acceptor) motifs (from
    # canonical_motifs — SAME set as the gate, "search what you filter") and emit
    # the PAIRED alternatives alongside the original chain. GTF-passthrough
    # transcripts are NOT extended (only read-derived novels go through search).
    # SIRV-tuned default ON (search_bp=4); net-neutral on SIRV (recall↑ cancels
    # precision↓). Set 0 to disable. Revisit for real transcriptomes.
    canonical_search_bp: int = 4           # 0 disables; 4 is the SIRV sweet spot
    max_chains_per_read: int = 16          # cap per-read alternative count (16 is SIRV sweet spot)

    # Canonical-motif GATE (Stage B): post-discovery FILTER that drops a NOVEL
    # multi-exon candidate if any internal junction's (donor,acceptor) motif is
    # not in canonical_motifs. This is the decisive precision lever on SIRV.
    # GTF-passthrough (source="gtf") and fusion (source="fusion") candidates are
    # EXEMPT — annotated junctions are trusted as-is. Single-exon (mono)
    # candidates have no internal junction and trivially pass.
    # SIRV WARNING: default ON + EXTENDED motif set tuned on synthetic SIRV;
    # revisit for real transcriptomes (toggle stays configurable).
    canonical_gate: bool = True
    canonical_motifs: tuple[str, ...] = ("GT-AG", "GC-AG", "AT-AC")

    # EM parameters
    em_sigma: float = 1.0
    em_beta: float = 0.5
    em_max_iter: int = 1000
    em_tol: float = 1e-4

    # R1 default path. enable_signal=False skips signal/EM and assigns reads by
    # mappy-AS argmax. This is now the PRODUCTION DEFAULT: the ablation champion
    # (M1-first) strictly dominates every M2/M3 signal variant on SIRV. SIRV
    # WARNING: tuned on synthetic SIRV; revisit for real transcriptomes (the
    # toggle stays configurable — set enable_signal=True for the legacy EM path).
    enable_signal: bool = False
    # R1 sub-variant. Only honored when enable_signal=False.
    #   "argmax_keep" : M1-keep SPLIT — each read assigned to ALL of its
    #                   simultaneously-best-AS candidates; abundance contribution
    #                   is split 1/K across the K tied winners (read mass
    #                   conserved). Mirrors ablation harness _quant_tie(mode=
    #                   "split") == config "M1-split". PRODUCTION DEFAULT (user
    #                   directive: "如果是同时最优，那么就都保留" with +1/K abundance).
    #   "argmax_first": M1-first HARD argmin — each read assigned to its single
    #                   best-AS candidate, ties broken by lowest candidate index
    #                   (== GTF prior: candidates are gtf_passthrough +
    #                   collapsed_novel + fusion). Integer read-count abundance.
    #                   Mirrors ablation harness _quant_tie(mode="first").
    #   "argmax_only" : mappy AS-weighted multimap → abundance = column sum (legacy R1)
    #   "argmax_em"   : mappy R_mm used as soft-init for EM (β=0, no signal),
    #                   then quantify_transcripts; aggregated with optional filter.
    r1_variant: Literal[
        "argmax_keep", "argmax_first", "argmax_only", "argmax_em"
    ] = "argmax_keep"
    # R2 uses em_max_iter_override=1 for single-step EM; None = use em_max_iter.
    em_max_iter_override: Optional[int] = None
    # R5 only: when False, min_abundance / min_max_r / min_novel_combined_score
    # are treated as 0 (single-switch filter, AC7).
    enable_score_filter: bool = True
    # m4 (read-to-read distance) source. "whole_read" preserves legacy
    # production behavior. "diff_region" uses the new intron-chain-derived
    # diff-region DTW (R4/R5 of the ablation). "none" forces a zero matrix
    # (no coherence contribution) regardless of em_beta.
    m4_source: Literal["whole_read", "diff_region", "none"] = "whole_read"

    # EM matrix subset selector. Controls which distance matrices feed EM.
    #   "m1"     : mappy distance only (read×tx); zeros for read×read; β=0
    #   "m2"     : eventalign distance only; zeros for read×read; β=0
    #   "m3"     : zeros for read×tx (uniform); m3 for coherence; β=em_beta
    #   "m1+m2"  : per-row z-score mean of M1/M2 → read×tx; zeros r×r; β=0
    #   "m1+m3"  : M1 → read×tx; M3 → read×read; β=em_beta
    #   "m2+m3"  : M2 → read×tx; M3 → read×read; β=em_beta  (DEFAULT, legacy)
    #   "all"    : zscore_mean(M1,M2) → read×tx; M3 → read×read; β=em_beta
    em_matrix_subset: Literal[
        "m1", "m2", "m3", "m1+m2", "m1+m3", "m2+m3", "all"
    ] = "m2+m3"
    # M2 normalization before EM. "none" = raw absolute distances (default,
    # legacy). "center" = per-row min subtraction (best→0, like M1).
    # "zscore" = per-row z-score shifted non-negative.
    m2_norm: Literal["none", "center", "zscore"] = "none"
    # M2 distance metric. "nll" = mean negative log-likelihood per event
    # (default, legacy). "gap" = weighted gap-run distance with flanking
    # context confidence — better discrimination for isoforms sharing exons.
    m2_metric: Literal["nll", "gap", "outlier"] = "nll"
    # M1+M2 fusion method for em_matrix_subset="m1+m2" or "all".
    m1m2_fusion: Literal["zscore_mean", "rank"] = "zscore_mean"
    # Post-hoc M2 validation: after M1-based EM, check if M2 agrees with
    # each read assignment. Discount abundance for candidates where M2
    # doesn't support the M1 assignment. Only effective when enable_signal
    # is True and em_matrix_subset uses M1 (not M2) for EM.
    m2_posthoc: bool = False
    m2_posthoc_top_k: int = 3  # M2 must rank the candidate in top-K
    # Signal tiebreaker: after M1 EM, resolve d=0 ties using diff-region
    # signal quality. Only effective when enable_signal=True.
    signal_tiebreak: bool = False
    tiebreak_ambig_threshold: float = 0.90
    # Fake-FASTQ tiebreak: binary fail/pass f5c test for d=0 ties.
    fake_fastq_tiebreak: bool = False
    # f5c_rna in-memory tiebreak (preferred over fake_fastq_tiebreak)
    f5c_rna_tiebreak: bool = False
    f5c_rna_pore: str = "rna002"
    # M2 tie resolution on the argmax_keep tie set (Stage M2-1). When a read is
    # simultaneously-best-AS across >=2 candidates, score the (small) tie set with
    # the validated junction-window mean-NLL metric (m2_resolve_tie); if the NLL
    # margin >= m2_tiebreak_margin, give the read's FULL mass to the M2-best
    # candidate instead of the 1/K split. Below the margin (or signal absent /
    # no discrimination window) -> keep the 1/K split. DEFAULT ON + AGGRESSIVE
    # (margin=1e-9: take M2's pick whenever it can discriminate at all). On the
    # competitor metric (gffcompare Tx-F1, NO-GTF) this peaks 45.4 vs M2-OFF 44.7,
    # beating all external tools AND the ablation M1-first champion 45.2. SIRV-tuned;
    # signal-absent auto-skips. All thresholds configurable (--no-m2-tiebreak reverts).
    m2_tiebreak: bool = True
    m2_tiebreak_junction_k: int = 10
    m2_tiebreak_margin: float = 1e-9
    # Scoring method for the R1 (enable_signal=False) path.
    # "mappy" = M1 mappy AS-gap distance (default).
    # "f5c_rna" = M2 f5c_rna match score distance (pure signal scoring).
    r1_scoring: Literal["mappy", "f5c_rna"] = "mappy"

    # EM prior / scoring
    score_alpha: float = 0.5          # weight for coherence vs discrimination in combined_score
    prior_weight_cap: float = 10.0    # max multiplicative boost from prior
    use_prior: bool = True            # apply combined_score-derived EM prior (set False for backward compat)

    # DTW
    use_gpu: bool = True
    max_reads_per_interval_for_dtw: int = 2000
    signal_normalize: bool = True   # per-read robust z-score before DTW

    # External tools
    f5c_path: str = "f5c"

    # Parallelism (future)
    num_workers: int = 1

    # Output
    output_gtf: Optional[str] = None
    output_tsv: Optional[str] = None        # scoring TSV output path
    output_bedpe: Optional[str] = None      # fusion BEDPE output path

    # Limits
    max_reads: Optional[int] = None

    # R-matrix persistence (T8: FP-by-EM data infrastructure)
    persist_R_matrix: bool = True  # write R.npy + R_meta.json per interval after EM

    # Diagnostic: write the per-candidate scoring TSV BEFORE the post-EM filters
    # (min_abundance / min_max_r / min_novel_combined_score) so downstream
    # FN-root-cause analysis can attribute drops to the correct filter rather
    # than mis-classifying filter drops as "missing_candidate". Path is
    # derived from output_tsv with the ".unfiltered.tsv" suffix.
    write_unfiltered_scores: bool = False

    # Fusion detection
    fusion_enabled: bool = False
    fusion_min_support: int = 2
    fusion_max_dist: int = 500     # bp window for breakpoint clustering
    fusion_flank_bp: int = 500     # bp to extract on each side of breakpoint for fusion sequence

    def validate(self):
        """Validate that required paths exist."""
        for attr in ("bam_path", "genome_fasta_path", "fastq_path", "signal_path"):
            val = getattr(self, attr)
            if val and not Path(val).exists():
                raise FileNotFoundError(f"{attr}: {val}")
