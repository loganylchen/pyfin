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
    min_abundance: float = 0.0        # A4: drop quantified NOVEL transcripts with abundance < this. Field default stays 0 (programmatic/quantify/ablation callers unchanged); the `fin` CLI defaults it to 3 (NOVEL only).
    floor_gtf_abundance: bool = False  # A4: when True the GTF floor is raised to max(min_gtf_abundance, min_abundance) (i.e. up to the NOVEL floor, never below the explicit GTF floor); default False = GTF uses its own lighter min_gtf_abundance. The `fin` CLI exposes this as --floor-gtf-abundance (OFF). Fusion is always exempt. Ablation path is separate and keeps GTF exempt regardless.
    min_gtf_abundance: float = 0.0    # A4: own (lighter) abundance floor for GTF transcripts, on soft EM abundance. Field default stays 0 (programmatic/quantify/ablation callers unchanged); the `fin` CLI defaults it to 1, so a GTF candidate whose EM soft-mass is below 1 read is dropped (avoids echoing annotation as a copy-tool). When floor_gtf_abundance is set the effective GTF floor becomes max(min_gtf_abundance, min_abundance). Fusion always exempt.
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
    # Signal-free: uses BAM primary-alignment spans only (no f5c CLI / krill).
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

    # polyA + 5'-proximity candidate-retention FILTER. Drop a candidate unless at
    # least `min_polya5p_reads` of its assigned reads BOTH (a) have a krill
    # whole-read polyA tail with polya_qc == "PASS" and polya_length >
    # `min_polya_length`, AND (b) map with their genomic 5' end within
    # `polya5p_window_bp` of the candidate's genomic 5' end (strand-aware).
    # Requires `signal_path`; it adds a krill whole-read eventalign (polya=True)
    # pass over all reads to every run, and no-ops gracefully when signal is
    # absent or krill produces nothing. Set `min_polya5p_reads = 0` to disable.
    # SIRV WARNING: SIRV-tuned and ON by default (=1). On SIRV no-GTF this lifts
    # gffcompare Tx-F1 (48.8 -> 49.5) by removing FP novels. Re-tune or disable
    # on real dRNA.
    #
    # `polya5p_exempt_gtf` (default True): GTF-sourced candidates are EXEMPT from
    # this filter, like fusion candidates — only `novel` candidates face the
    # polyA evidence bar. This is the validated optimal configuration: a 6-sample
    # x 6-ratio SIRV benchmark showed exempting GTF makes pyfin the corrected-recall
    # leader at every annotation-completeness ratio with precision only ever rising
    # vs gating, because dRNA reads are frequently 5'-truncated or lack a detectable
    # polyA tail and would otherwise drop genuine annotated transcripts. Set
    # `polya5p_exempt_gtf = False` to restore the old behavior of gating GTF too.
    min_polya5p_reads: int = 1
    polya5p_window_bp: int = 25          # bp 5' tolerance (mirrors fulllen_window_bp)
    min_polya_length: float = 10.0       # min krill polya_length for a read to support (qc PASS required)
    polya5p_exempt_gtf: bool = True      # exempt GTF candidates from polyA+5' (like fusion); False gates GTF too

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
    # M3 read×read coherence weight in em_with_coherence (only applied when
    # m3_coherence=True). Default 1.0 = the SIRV beta-sweep sweet spot for the
    # pure tie-break junction-NLL EM (with-GTF Tx-F1 59.7 at β=1 vs 59.2 M3-off);
    # β≥2 over-smooths. SIRV-tuned; revisit on real dRNA.
    em_beta: float = 1.0
    em_max_iter: int = 1000
    em_tol: float = 1e-4

    # RSEM/Salmon-style iterative abundance feedback inside the EM (assignments.
    # em_with_coherence). When ON, each M-step re-estimates a transcript-abundance
    # vector theta from the current responsibilities and biases the energy by
    # -sigma*log(theta_j), so a read shared between candidates is pulled toward the
    # more abundant one (a 0.5/0.5 split migrating toward 0.99/0.01 once theta
    # diverges, e.g. 100:1). Only affects the EM quant modes ("m1_em"/"m2_em");
    # the production "argmax" mode does no EM and is unaffected. OFF by default
    # (production output is byte-identical). abundance_length_norm additionally
    # divides theta counts by per-transcript effective length (Salmon-style);
    # requires the runner to supply eff_lengths. SIRV-experimental.
    abundance_feedback: bool = False
    abundance_length_norm: bool = False

    # Quantification engine selector (replaces the enable_signal/r1_variant/
    # r1_scoring triplet). The three modes are mutually exclusive:
    #   "argmax" : mappy AS argmax + M2 krill junction tiebreak, hard counts,
    #              NO EM. (M1-first; on SIRV peaks slightly above m2_em.)
    #   "m1_em"  : EM seeded by the M1 mappy AS-gap distance (beta=0, no signal
    #              coherence). Pure-alignment soft assignment.
    #   "m2_em"  : PRODUCTION DEFAULT. EM seeded by the PURE tie-break junction-NLL
    #              M2 distance: M1/AS picks each read's best-AS tie set + mappability
    #              mask, and the per-event junction-NLL is the sole graded distance
    #              over that tie set (m2_resolve_tie semantics, NOT the dense
    #              read×candidate matrix). Optionally adds the M3 read×read DTW
    #              coherence when m3_coherence=True (beta=em_beta). Chosen as the
    #              default for real-dRNA robustness (principled soft assignment);
    #              on SIRV it clears the competitor floor (Tx-F1 ~48.9 no-GTF /
    #              59.2 with-GTF, M3 off) just under argmax.
    # All signal scoring is krill (in-memory eventalign); the legacy f5c CLI
    # path is gone. SIRV WARNING: tuned on synthetic SIRV; revisit for real data.
    quant_mode: Literal["argmax", "m1_em", "m2_em"] = "m2_em"

    # R2 uses em_max_iter_override=1 for single-step EM; None = use em_max_iter.
    em_max_iter_override: Optional[int] = None
    # R5 only: when False, the post-EM abundance/fraction filters are treated as 0
    # (single-switch filter, AC7).
    enable_score_filter: bool = True
    # m4 (read-to-read distance) source. "whole_read" preserves legacy
    # production behavior. "diff_region" uses the new intron-chain-derived
    # diff-region DTW (R4/R5 of the ablation). "none" forces a zero matrix
    # (no coherence contribution) regardless of em_beta.
    m4_source: Literal["whole_read", "diff_region", "none"] = "whole_read"

    # Read×read junction-window DTW coherence (M3) in the EM quant modes. OFF by
    # default because the pairwise DTW is expensive; the production m2_em default
    # runs the pure tie-break junction-NLL EM with NO M3 (no DTW). When True the
    # runner builds the M3 read×read matrix and mixes it into em_with_coherence at
    # beta=em_beta. SIRV: with-GTF Tx-F1 59.7 (M3 on, β=1) vs 59.2 (off) — small,
    # high-precision, recall-neutral; the payoff is expected on real dRNA. Gates
    # M3 in the assembly runner (m2_em).
    m3_coherence: bool = False

    tiebreak_ambig_threshold: float = 0.90
    # krill in-memory tiebreak
    krill_tiebreak: bool = False
    krill_pore: str = "rna002"
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
    # M2-EM diff-region coverage gate (quant_mode="m2_em" only). Default ON. For a
    # read in a >=2 best-AS tie, only let it discriminate wobble siblings when its
    # eventalign signal STRADDLES the wobbling junction(s) (donor->acceptor span;
    # see scoring.m2_junction_nll.wobble_diff_spans). Per read x tie:
    #   * M2 margin (2nd-best NLL - best NLL) >= m2_diff_cover_margin -> HARD assign
    #     the read's full mass to the lowest-NLL candidate; if it also covered every
    #     diff region it contributes that vote to the locus isoform-ratio prior;
    #   * margin < threshold (whether or not the read covers) -> AMBIGUOUS: its tie
    #     mass is redistributed in proportion to the prior (the ratio learned from
    #     covered+distinguishing reads). If no tie candidate earned prior votes ->
    #     1/K flat (signal-dead loci are never starved). The margin is thus the sole
    #     decider of hard-assign vs ratio-follow.
    # No read is ever dropped (recall-safe). OFF reverts to the soft NLL-graded d_tx
    # skeleton (byte-identical to the prior behavior). SIRV/dense-locus precision
    # lever; on real dRNA re-tune the margin.
    m2_diff_cover_gate: bool = True
    m2_diff_cover_margin: float = 0.5
    # M2-EM cluster-internal wobble recheck (quant_mode="m2_em" only). Default ON.
    # AFTER the full EM, cluster ALL multi-exon candidates by STRUCTURE (same
    # chrom/strand, same intron count, every junction within m2_cluster_recheck_bp);
    # within a cluster the highest-EM-abundance candidate anchors (usually the true
    # isoform, often a GTF passthrough) and a NOVEL sibling whose abundance is below
    # m2_cluster_recheck_fraction of the anchor is a wobble shadow and is dropped.
    # GTF/fusion are never dropped (they only anchor). The judgement is pure
    # abundance: GTF participates in the abundance race but is NOT a correctness
    # oracle, so precision stays robust when the annotation is wrong; with no GTF a
    # real novel anchors its own cluster. OFF == no drops (byte-identical).
    # m2_cluster_recheck_fraction <= 0 falls back to min_isoform_fraction.
    # bp=10 / fraction=0.15 tuned on the SGNex heya8 + sirv4 12-condition matrix
    # (nanocount expressed-truth): lifts dense-locus F1@3 by ~6 pts toward the head
    # tools while holding full Sn@3=100 (no real isoform dropped); re-tune on real
    # dRNA if genuine sub-15%-abundance wobble isoforms are expected.
    m2_cluster_recheck: bool = True
    m2_cluster_recheck_bp: int = 10
    m2_cluster_recheck_fraction: float = 0.15
    # Score-gated fallback split (mk1). On a tie that M2 does NOT win
    # outright (margin < m2_tiebreak_margin), restrict the 1/K split to the
    # candidates eventalign could actually score in the junction window; the
    # candidates with NO window signal are excluded for THIS read. If NONE
    # scored, keep the full 1/K split over all tied candidates (== M1). The read
    # is never dropped, so this is recall-safe. On SIRV it is a no-op (only ~2
    # "mixed" reads; verified gffcompare == baseline); kept ON as the principled
    # per-read rule that may matter on real dRNA. Configurable.
    m2_tie_scoregate_split: bool = True

    # EM prior / scoring. VESTIGIAL: the combined_score-derived EM prior was driven
    # by the composite scorer, which has been removed; the assembly runner never
    # reads these. Retained as inert config defaults for backward compatibility.
    score_alpha: float = 0.5          # (inert) weight for coherence vs discrimination in combined_score
    prior_weight_cap: float = 10.0    # (inert) max multiplicative boost from prior
    use_prior: bool = True            # (inert) no code path applies a prior in assembly

    # DTW
    use_gpu: bool = True
    max_reads_per_interval_for_dtw: int = 2000
    signal_normalize: bool = True   # per-read robust z-score before DTW

    # Parallelism: interval-level process workers. threads=1 keeps the serial
    # path verbatim. gpu_workers of the threads workers hold a GPU context (VRAM
    # bound = G x per-context footprint); the rest are CPU-only. gpu_workers must
    # be 0 when use_gpu is False.
    threads: int = 1
    gpu_workers: int = 0

    # Output
    output_gtf: Optional[str] = None
    output_tsv: Optional[str] = None        # scoring TSV output path
    output_bedpe: Optional[str] = None      # fusion BEDPE output path

    # Limits
    max_reads: Optional[int] = None

    # R-matrix persistence (T8: FP-by-EM data infrastructure)
    persist_R_matrix: bool = True  # write R.npy + R_meta.json per interval after EM

    # Diagnostic: write the per-candidate scoring TSV BEFORE the post-EM filters
    # (min_abundance / isoform-fraction / full-length / polyA) so downstream
    # FN-root-cause analysis can attribute drops to the correct filter rather
    # than mis-classifying filter drops as "missing_candidate". Path is
    # derived from output_tsv with the ".unfiltered.tsv" suffix.
    write_unfiltered_scores: bool = False

    # Fusion detection
    fusion_enabled: bool = False
    fusion_min_support: int = 2
    fusion_max_dist: int = 500     # bp window for breakpoint clustering
    fusion_flank_bp: int = 500     # bp to extract on each side of breakpoint for fusion sequence
    # Adapter-bridged chimera guard (nanopore dRNA false fusions): a soft-clip
    # read whose two arms are separated by an internal stretch of read sequence
    # that maps to NEITHER arm (typically an ONT internal adapter, ~adapter
    # length) is a sequencing artifact, not a real fusion. Drop the read when the
    # internal unmapped gap (read coords, between the two aligned arms) is >= this
    # many bp. See DeepChopper (Nat. Commun. 2026). 0 disables the guard.
    fusion_max_internal_gap_bp: int = 30

    def validate(self):
        """Validate that required paths exist and parallelism knobs are coherent."""
        for attr in ("bam_path", "genome_fasta_path", "fastq_path", "signal_path"):
            val = getattr(self, attr)
            if val and not Path(val).exists():
                raise FileNotFoundError(f"{attr}: {val}")
        if self.threads < 1:
            raise ValueError(f"threads must be >= 1, got {self.threads}")
        if not 0 <= self.gpu_workers <= self.threads:
            raise ValueError(
                f"gpu_workers must be in [0, threads={self.threads}], got {self.gpu_workers}"
            )
        if self.gpu_workers and not self.use_gpu:
            raise ValueError("gpu_workers must be 0 when use_gpu is False")
