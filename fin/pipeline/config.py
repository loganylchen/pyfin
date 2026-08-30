"""Pipeline configuration."""

from __future__ import annotations

from dataclasses import dataclass
import logging
from pathlib import Path
from typing import Iterable, Literal, Mapping, Optional


logger = logging.getLogger(__name__)


# CLI operating points. Programmatic PipelineConfig construction keeps the
# dataclass defaults below; profiles are resolved by fin.cli before construction.
PIPELINE_PROFILES: dict[str, dict[str, object]] = {
    "sirv": {
        "min_abundance": 3.0,
        "strict_novel_abundance_floor": False,
        "min_gtf_abundance": 1.0,
        "floor_gtf_abundance": True,
        "min_isoform_fraction": 0.01,
        "post_selection_refit": True,
        "max_soft_mass_ratio": 2.0,
        "min_fulllen_fraction": 0.1,
        "min_polya5p_reads": 0,
        "canonical_gate": True,
        "novel_junction_min_reads": 2,
        "containment_cluster": True,
        "m2_cluster_recheck": True,
        "m2_metric": "auto",
        "m2_diff_cover_margin": 0.5,
        # Inactive under the mean default; these are the best SIRV summed
        # fallback values when a user explicitly switches only m2_metric.
        "m2_summed_llr_margin": 1.0,
        "m2_summed_llr_flank": 4,
    },
    "real-drna": {
        "min_abundance": 1.0,
        "strict_novel_abundance_floor": True,
        "min_gtf_abundance": 1.0,
        "floor_gtf_abundance": False,
        "min_isoform_fraction": 0.01,
        "post_selection_refit": True,
        "max_soft_mass_ratio": 0.0,
        "min_fulllen_fraction": 0.0,
        "min_polya5p_reads": 0,
        "drop_mono_exon_novel": True,
        "min_mono_exon_reads": 5,
        "junction_snap": True,
        "junction_snap_tolerance": 6,
        "junction_snap_min_support": 2,
        "junction_snap_min_ratio": 2.0,
        "canonical_gate": True,
        "novel_junction_min_reads": 2,
        "containment_cluster": True,
        "m2_cluster_recheck": True,
        "m2_metric": "summed_llr",
        "m2_diff_cover_margin": 0.5,
        "m2_summed_llr_margin": 1.0,
        "m2_summed_llr_flank": 8,
    },
    "custom": {},
}
PIPELINE_PROFILES["real-drna-precision"] = {
    **PIPELINE_PROFILES["real-drna"],
    "min_novel_reads": 2,
}
PROFILE_FIELDS = frozenset(
    key for values in PIPELINE_PROFILES.values() for key in values
)


def gtf_has_usable_guide(gtf_path: Optional[str], min_transcripts: int = 2) -> bool:
    """Use the benchmark convention: fewer than two GTF transcripts is unguided."""
    if not gtf_path:
        return False
    path = Path(gtf_path)
    if not path.exists():
        return False
    found = 0
    with path.open(errors="ignore") as handle:
        for line in handle:
            if "\ttranscript\t" in line:
                found += 1
                if found >= min_transcripts:
                    return True
    return False


def resolve_m2_metric(metric: str, gtf_path: Optional[str]) -> tuple[str, str]:
    """Resolve the profile-only auto sentinel to a concrete assignment metric."""
    if metric != "auto":
        return metric, "fixed"
    if gtf_has_usable_guide(gtf_path):
        return "mean", "auto-guided"
    return "summed_llr", "auto-unguided"


def resolve_profile_values(
    profile: str,
    values: Mapping[str, object],
    explicit_fields: Iterable[str] = (),
) -> dict[str, object]:
    """Overlay a named profile while preserving explicit caller values."""
    if profile not in PIPELINE_PROFILES:
        choices = ", ".join(sorted(PIPELINE_PROFILES))
        raise ValueError(f"unknown profile {profile!r}; expected one of: {choices}")
    explicit = set(explicit_fields)
    resolved = dict(values)
    for key, value in PIPELINE_PROFILES[profile].items():
        if key in resolved and key not in explicit:
            resolved[key] = value
    return resolved


@dataclass
class PipelineConfig:
    """Configuration for the full pyfin pipeline."""

    # Input files
    bam_path: str

    # Resolved CLI operating point. "custom" is the programmatic/back-compat
    # default; fin.cli records command-line fields that overrode a named profile.
    profile: str = "custom"
    profile_overrides: tuple[str, ...] = ()

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
    # NEW generation-side clustering (PRODUCTION DEFAULT). Group reads by intron chain
    # (3' ignored), collapse EXACT sub-chains, cluster wobble/cassette/containment
    # keeping members; NO canonical-search shadow generation. Replaces the legacy
    # 3'+exact-chain grouping + canonical expansion. On p00 (no gates, candidates = the
    # result) this cut raw candidates 175k->37k, c -73%, j -88%, structural precision
    # 10.1%->26.4% (2.6x), at a ~20% distinct-truth recall cost that is 92% attributable
    # to the exact sub-chain collapse (recoverable later by guarding the collapse).
    chain_cluster_discovery: bool = True
    clustering: str = "families"  # generation clustering primitive (--chain-cluster-discovery only). "families" = clustering redesign (PRODUCTION DEFAULT): cluster_families (grouping only, single-linkage) + collapse (explicit exact-subchain fold, span-guarded); mono reads deferred to a holding bucket resolved post-EM (mono_resolve_post_em). Emits the SAME multi-exon candidate SET as the legacy path (SIRV 250=250, gencode p00 51069=51069/82.4% FL; the EM cluster scoping may differ in rare exact-subchain-bridge cases, and mono handling differs by design -- deferred vs generation-folded) and the full pipeline BEATS the old default on real human p00 (recall +0.4, Pr +1.9). "read_chains" = legacy cluster_read_chains (fold+cluster inline); pair with --no-mono-resolve-post-em for the exact old behavior.
    chain_cluster_wobble_bp: int = 6
    chain_cluster_cassette_max_exon_bp: int = 70
    chain_cluster_fold_monoexon: bool = True  # generation-time mono fold: a single-exon read wholly inside a multi-exon candidate's exon folds into it (5'/3' degradation fragment) instead of a standalone mono candidate; intronic/uncontained mono stay separate. ONLY effective when post-EM mono resolution is NOT running (mono_resolve_post_em False, or non-m2_em quant mode); under the production default (mono_resolve_post_em=True) the post-EM resolution supersedes it and this flag is inert. --no- reverts to standalone mono candidates.
    chain_cluster_fold_span_guard: bool = True  # only fold a read into an exact-sub-chain container if its aligned span does NOT run exonically across one of the container's EXTRA introns; a read spanning such an intron (retained-intron / alternative isoform) is kept as its own candidate instead of being absorbed. PRODUCTION DEFAULT. --no- reverts to unconditional exact-sub-chain fold.
    # Post-EM mono-exon resolution (quant_mode="m2_em" only). When True, single-exon
    # reads are NOT folded into multi-exon candidates at generation (chain_cluster_fold_
    # monoexon is forced off); instead, AFTER the multi-exon EM + structural gates, each
    # surviving mono candidate's reads are re-resolved against the SURVIVING multi-exon
    # candidates by strict strand-aware exonic containment: a read wholly inside one
    # exon of exactly one surviving multi folds into it; inside several -> folds into the
    # highest-EM-abundance one; uncovered reads stay on the mono candidate, which is kept
    # only if it retains >= mono_resolve_min_reads (else dropped). Fragment reads follow
    # inferred multi-isoform abundance instead of a static generation-time container.
    # PRODUCTION DEFAULT (ON): pairs with clustering="families" (which defers mono to a bucket).
    # On real human gencode p00 this survivor-conditioned, cross-cluster resolution BEATS the old
    # generation-time fold_monoexon (recall 33.0->33.4, Pr 70.9->72.8) -- a mono read can be shared
    # across clusters, so deferring the host decision until after multi-exon EM/selection (against
    # SURVIVING candidates) is both more principled and better. When False, mono handling falls back
    # to the generation-time fold (chain_cluster_fold_monoexon, on by default) -- or standalone mono
    # candidates if that is also off. Only active for quant_mode="m2_em".
    mono_resolve_post_em: bool = True
    mono_resolve_min_reads: int = 2   # uncovered reads a mono candidate must retain to survive
    mono_resolve_slop_bp: int = 10    # terminal boundary slop for the exonic-containment test
    min_abundance: float = 0.0        # A4: NOVEL soft-abundance floor. Programmatic default stays 0; CLI profiles resolve SIRV=3 and real-dRNA=1.
    strict_novel_abundance_floor: bool = False  # False keeps abundance == floor (SIRV >=3); True requires abundance > floor (real-dRNA >1 read-equivalent). This encodes boundary semantics explicitly instead of an epsilon/magic threshold.
    floor_gtf_abundance: bool = False  # A4: when True the GTF floor is raised to max(min_gtf_abundance, min_abundance). Programmatic/real-dRNA default False; the optimized SIRV profile enables it. Fusion is always exempt.
    min_gtf_abundance: float = 0.0    # A4: own lighter GTF soft-abundance floor. Programmatic default stays 0; both named CLI profiles resolve it to 1. With floor_gtf_abundance, the effective floor is max(min_gtf_abundance, min_abundance). Fusion always exempt.
    # Minimum isoform fraction (locus-relative abundance) FILTER. Drop a NOVEL
    # multi-exon transcript whose abundance is below this fraction of the
    # dominant transcript in its persisted splice family. This is the standard
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
    # "family" uses the discovery splice family; candidates without a family
    # ID fall back to historical same-strand genomic overlap. "overlap" forces
    # the historical denominator for every candidate.
    isoform_fraction_locus: Literal["family", "overlap"] = "family"

    # Refit read responsibilities after GLOBAL selection and finalized junction
    # snapping. Named production profiles enable this; the raw/programmatic
    # default stays off for historical reproduction. The effective/disable
    # fields are resolved by validate() and serialized in the run manifest.
    post_selection_refit: bool = False
    post_selection_refit_effective: bool = False
    post_selection_refit_disable_reason: Optional[str] = None

    # Emit candidate_evidence.tsv: one row per post-selection survivor with
    # inference-time observable features (junction support, family share,
    # end support, canonicality, containment geometry). Read-only diagnostic:
    # never changes candidates, abundances, or ordinary outputs. Off by
    # default; the ranking work consumes this table offline.
    candidate_evidence: bool = False

    # Calibrated candidate ranking (fin/analysis/candidate_ranking.py).
    # "filter" scores every post-selection survivor with the frozen v1
    # logistic model over inference-time observable evidence and removes
    # NOVEL candidates below the frozen operating point BEFORE junction
    # snapping and the final abundance refit (GTF/fusion always exempt).
    # Trained on H9 r2r2 only; threshold frozen on the tuning frontier under
    # the hard T1-not-lower constraint and validated untouched on r3r1.
    # "off" preserves the historical selection exactly.
    ranking_mode: Literal["off", "filter"] = "off"
    # None -> the frozen model threshold; explicit values override (expert).
    ranking_threshold: Optional[float] = None

    # EXPERIMENTAL observability: also write per-interval machine-readable M2
    # contrast records (per-comparison event counts/values/margins plus one
    # aggregate line) as JSONL under <work_dir>/m2_contrasts/. Off by
    # default; log lines always exist.
    m2_contrast_stats_jsonl: bool = False

    # EndpointRefine (EXPERIMENTAL, off by default; requires the final
    # abundance refit because post-split requantification is mandatory).
    # Splits a novel multi-exon survivor into at most endpoint_max_splits
    # endpoint states when strand-aware read-end modes support distinct
    # (TSS, TES) pairs; interior-TSS (degradation-direction) modes need
    # 2x support, and every read is re-routed through the refit so mass
    # conservation is enforced by the existing invariants.
    endpoint_refine: bool = False
    endpoint_window_bp: int = 25
    endpoint_min_reads: int = 3
    endpoint_min_pair_frac: float = 0.15
    endpoint_max_splits: int = 2

    # TSS evidence (EXPERIMENTAL). Decides whether a CONTAINED shorter model
    # is a real transcript or a 5'-degradation artifact of the longer one, by
    # testing its start against the local conditional-termination-hazard
    # background (fin/analysis/tss_evidence.py). "audit" records verdicts
    # without changing output; "require" keeps only endpoint states whose
    # alternative TSS is `supported`. An `unidentifiable` verdict never drops
    # a model - insufficient evidence is not evidence of absence.
    tss_evidence_mode: Literal["off", "audit", "require"] = "off"

    # Genome access. lazy_genome=True (default) opens the FASTA through an
    # indexed lazy mapping holding at most genome_cache_chroms chromosomes
    # per process - measured to be the dominant worker-memory hotspot when
    # loaded eagerly (~3.1 GB x N workers). False restores the historical
    # eager whole-genome dict.
    lazy_genome: bool = True
    genome_cache_chroms: int = 2

    # Soft-mass / hard-read ratio ceiling. Drop a NOVEL multi-exon candidate
    # whose EM soft abundance (R.sum) divided by its hard argmax read count
    # (num_assigned_reads) is >= this value (a candidate with 0 hard reads is
    # always dropped). A real isoform deposits ~1 soft mass per hard read
    # (ratio ~1); a wobble shadow borrows fractional soft crumbs from a
    # high-abundance structural near-copy and shows an inflated ratio, so this
    # catches the HIGH-relative-abundance shadows the locus-fraction /
    # cluster-recheck levers miss. Pure EM evidence (no GTF). gtf/fusion/mono
    # exempt. 0 disables. Default 2.0 tuned on the SGNex heya8+sirv4 12-condition
    # nanocount matrix: ALL 24 cells F1@3 up-or-equal (mean +1.04), recall held
    # on every cell (true isoforms cluster at ratio ~1, a few reach ~1.4-1.6, so
    # 2.0 clears them). The SIRV profile keeps 2.0; the real-dRNA profile
    # disables this gate because its recall gain outweighed the tiny F1 change.
    max_soft_mass_ratio: float = 2.0

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
    # disable. SIRV WARNING: the SIRV-profile 0.1 is tuned to that synthetic
    # domain; the real-dRNA profile resolves this field to 0. It drops ~74% of
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
    # Historical SIRV runs enabled this at 1, but the current two-replicate
    # p00/full honest-F1 sweep selected 0; both named profiles now disable it.
    # It remains available as an explicit A/B lever.
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
    #              read×candidate matrix). Production uses no read×read coherence.
    #              Chosen for real-dRNA robustness (principled soft assignment).
    # All signal scoring is krill (in-memory eventalign); the legacy f5c CLI
    # path is gone. SIRV WARNING: tuned on synthetic SIRV; revisit for real data.
    quant_mode: Literal["argmax", "m1_em", "m2_em", "cluster"] = "m2_em"

    # quant_mode="cluster" only: assign reads to candidates WITHIN each generation
    # cluster (CandidateSet.clusters) using fin.pipeline.cluster_quant. A straddling
    # M1 tie is resolved by the summed-LLR junction signal only when |LLR| >=
    # cluster_llr_threshold; a member survives iff its total assigned weight >=
    # cluster_min_support. Defaults leave every other mode unchanged.
    cluster_llr_threshold: float = 2.0
    cluster_min_support: float = 1.0
    # M1 best-AS tie margin (AS points): a read is unique-best only when its top
    # member beats the runner-up by more than this; within it the members tie (-> M2 /
    # ambiguous 1/k -> EM). Stops a 2-3 bp wobble near-tie from anchoring a shadow
    # (measured: wobble AS-gaps concentrate <=20, real differences differ by hundreds).
    cluster_m1_tie_margin: float = 20.0
    # cluster mode: run the summed-LLR M2 signal tiebreak on straddling ties. With a
    # wide m1_tie_margin this fires on many ties (expensive eventalign) for little
    # end-to-end gain, so it can be turned off (near-ties then fall straight to the
    # ambiguous 1/k -> EM path). Signal is still used by the polyA finalize gate.
    cluster_use_m2: bool = True

    # R2 uses em_max_iter_override=1 for single-step EM; None = use em_max_iter.
    em_max_iter_override: Optional[int] = None
    # R5 only: when False, the post-EM abundance/fraction filters are treated as 0
    # (single-switch filter, AC7).
    enable_score_filter: bool = True
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
    # Metric used by the default m2_em tie scorer. "mean" is the legacy wide
    # class-window mean NLL. "summed_llr" uses tight differing-junction windows
    # and undivided per-event NLL, restricted to same-intron-count contrasts.
    # "off" retains the M1 exact tie and assigns it without signal refinement.
    m2_metric: Literal["off", "mean", "summed_llr", "sqrt_count_mean_llr", "auto"] = "mean"
    # "sqrt_count_mean_llr" (EXPERIMENTAL, never auto-routed): same tight
    # differing-junction footprint and same-intron-count guard as summed_llr,
    # but each candidate's MEAN window NLL is rescaled by sqrt(min effective
    # event count), so the decision margin is a z-like difference of means
    # instead of an undivided sum whose magnitude tracks event count. Gated by
    # the same m2_summed_llr_margin threshold. Promote only after paired live
    # validation; observation counters are logged for both tight metrics.
    m2_metric_route: str = "fixed"
    m2_summed_llr_margin: float = 2.0
    m2_summed_llr_flank: int = 6
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
    # bp=20 / fraction=0.15 from the 5-sample heya8 + sirv4 full sweep:
    # 40/40 cells precision up, 36/40 F1@3 up vs bp=10 (the prior SIRV-tuned default);
    # full Sn@3=100 preserved everywhere it was already 100. cassette_max_exon_bp=70
    # extends the cluster equivalence to K vs K-1 intron pairs differing by one small
    # exon (<70bp), so minimap2's small-exon-skip artifacts cluster with their true
    # sibling and the abundance-fraction filter drops the shadow; this lever fires
    # exactly on p00 / c_jitter10bp where bp=20 alone leaves cassette FPs.
    m2_cluster_recheck: bool = True
    m2_cluster_recheck_bp: int = 20
    m2_cluster_recheck_fraction: float = 0.15
    m2_cluster_recheck_cassette_max_exon_bp: int = 70
    # Allow a clustered low-support GTF sibling (below fraction*anchor) to be dropped
    # as a wobble shadow — BUT only when a direct read-support guard fires: the GTF's
    # distinguishing junction (vs the anchor) carries fewer than
    # m2_cluster_recheck_gtf_min_jct_reads reads splicing EXACTLY there
    # (m2_cluster_recheck_jct_tol bp tolerance). Rationale: an annotated isoform is a
    # documented hypothesis, so abundance alone must not delete it; only the data
    # totally failing to traverse its specific junction (a jittered/phantom passthrough,
    # 0 exact reads) justifies the drop. A genuine annotated isoform keeps its own reads
    # (e.g. SIRV606: 419 reads on its exact donor) and survives. Symmetric in the
    # anchor's source — a GTF is judged by its OWN read support. Default ON; set False
    # for the legacy "GTF never dropped" behaviour. jct_tol=0 (strict exact match):
    # validated on heya8 c_jitter10bp — all jittered-GTF FPs have 0 exact reads at their
    # junction while real isoforms have many; a loose tol lets small (1-2bp) jitters
    # borrow support from the true neighbouring junction.
    m2_cluster_recheck_novel_displaces_gtf: bool = True
    m2_cluster_recheck_gtf_min_jct_reads: int = 1
    m2_cluster_recheck_jct_tol: int = 0
    # Lever 1 — containment / 5'-truncation collapse (quant_mode="m2_em" only).
    # Post-EM, fold a NOVEL candidate whose intron chain is a pure 3' SUFFIX of a
    # longer candidate (a 5'-truncation shadow: same downstream junctions, 3'
    # terminus within containment_3p_tol_bp, 5' end interior) into that longer
    # parent: the shadow's reads + soft mass are reassigned to the parent and the
    # shadow is dropped. Parent may be gtf or novel; a gtf candidate is NEVER
    # folded away (only novels are droppable). The suffix match is STRICT (exact
    # intron equality on the shared suffix) so genuine exon-skipping / alt-3'-end
    # isoforms (different splicing) are never folded. WARNING: this rule cannot
    # distinguish a dRNA 5'-truncation artifact from a genuine low-abundance
    # alt-TSS isoform that starts at a downstream exon and shares the 3' end, so it
    # is NOT recall-safe by construction; default OFF, enable only after real-truth
    # honest-F1 validation. OFF (default) -> no folds (byte-identical).
    containment_collapse: bool = False
    containment_3p_tol_bp: int = 20
    # Fold the shadow only when its EM abundance <= parent's * this ratio (shadow
    # must be the minor member). 1.0 = fold any shadow at or below the parent.
    containment_min_abundance_ratio: float = 1.0
    # Containment-CLUSTER drop (recall-SAFER generalisation of containment_collapse):
    # drop a NOVEL candidate whose intron chain is a contiguous SUB-CHAIN (within
    # containment_cluster_wobble_bp per junction) of a longer candidate — the
    # truncation / exon-skip shadow the same-intron-count wobble cluster never groups
    # — ONLY when it is a low-support shadow by BOTH EM abundance (<= parent *
    # min_ab_ratio) AND supporting-read count (<= parent * min_read_ratio). The
    # read-support guard is what containment_collapse lacked: measured on p00 a
    # truncation shadow carries ~1 read while a genuine short/alt-TSS isoform carries
    # reads comparable to the parent (median 13), so requiring BOTH keeps MOST real short
    # isoforms — but a genuine low-fraction minor isoform (e.g. parent 100 reads, real
    # 20-read isoform <30% abundance) can still be dropped, so this is recall-safER, not
    # recall-safe. DEFAULT-ON (production): a targeted, read-guarded precision lever that
    # drops truncation/exon-skip mapping shadows the wobble cluster misses. Validated on
    # the gencode sweep (p00 honestF1 +1.0 at <=0.5 corrRec cost; c_jitter/full
    # non-regressing in the containment-only ablation). The absolute cap
    # (containment_cluster_max_shadow_reads) bounds the residual recall risk.
    containment_cluster: bool = True
    containment_cluster_wobble_bp: int = 6
    containment_cluster_min_ab_ratio: float = 0.3
    containment_cluster_min_read_ratio: float = 0.3
    # absolute cap: never drop a shadow carrying more than this many supporting reads,
    # regardless of the ratio (a real low-fraction isoform of a very-high-support parent
    # can exceed min_read_ratio*parent). 0 disables (ratio-only). 10 keeps the validated
    # benefit (shadows are ~1-3 reads) while never folding a >10-read candidate.
    containment_cluster_max_shadow_reads: int = 10
    # Lever 3 — mono-exon (single-exon) read-support gate (post-aggregate, in
    # _finalize_and_write). Drop a NOVEL single-exon candidate whose hard read
    # count < min_mono_exon_reads OR genomic length < min_mono_exon_length.
    # Suppresses single-exon de novo noise (IsoQuant drops novel unspliced by
    # default for ONT) WITHOUT a blanket hard drop: a high-support / long real
    # intronless gene (histone, many ncRNAs) survives. gtf/fusion/multi-exon are
    # EXEMPT (only novel mono are gated). Master switch drop_mono_exon_novel must
    # be True AND at least one threshold > 0 to fire; otherwise no drops
    # (byte-identical). Defaults OFF.
    drop_mono_exon_novel: bool = False
    min_mono_exon_reads: int = 0
    min_mono_exon_length: int = 0
    # Finalized-model junction consensus correction. Novel multi-exon models
    # may snap each intron to a more strongly supported exact CIGAR junction
    # within tolerance, then structurally identical corrected models merge with
    # their abundance/read mass preserved. GTF/fusion are exempt. Default off
    # until the live multi-sample validation promotes a profile setting.
    junction_snap: bool = False
    junction_snap_tolerance: int = 6
    junction_snap_min_support: int = 2
    junction_snap_min_ratio: float = 2.0
    # Lever 2 — per-junction read-support gate (quant_mode="m2_em" only). Drop a
    # NOVEL multi-exon candidate if ANY of its junctions is spliced by fewer than
    # novel_junction_min_reads directly-observed reads (intron junctions extracted
    # from primary-read CIGARs, strand-keyed, matched within novel_junction_reads
    # _tol bp). I.e. a novel junction must be carried by >= N independent reads,
    # not just 1. gtf/fusion/mono are EXEMPT (only novel multi-exon are gated).
    # <= 1 disables (every assembled junction trivially has >=1 read).
    # DEFAULT 2 (production, landed after the full off-vs-on sweep): a novel
    # junction must be carried by >= 2 reads. On real human GENCODE it lifts
    # precision at ZERO recall cost across EVERY scenario — de novo Pr +3.8,
    # p10 +3.2, p50 +1.5, p90/full +0.8, and all six corrupted-GTF ratios
    # (c_jitter/skip/flip/merge/spurious/ir) +0.6..+0.9 — and keeps pyfin #1 on
    # SIRV4 (F1@3 85.0, wash). NOTE: on SYNTHETIC de novo / aggressive-corruption
    # this trades ~0.4-0.8 F1@3 (recall) for precision, but the real-data arbiter
    # and #1 standing justify it. Set 0 to disable.
    novel_junction_min_reads: int = 2
    novel_junction_reads_tol: int = 2
    # DE-NOVO wobble-tolerant collapse (EXPERIMENT default-OFF). At discovery,
    # novel candidates are bucketed by EXACT intron chain, so a junction at
    # (1000,2000) and a minimap-wobbled (1002,1998) become SEPARATE candidates
    # (wobble shadows) — the measured cause of pyfin's low de-novo structural
    # precision (p00 65.6% vs isoquant 91.3%). When > 0, _collapse_candidates
    # merges novel candidates whose intron chains match within this many bp per
    # junction (AND 3' within three_prime_threshold) into the highest-read-support
    # representative (consensus of the reads' OWN mode — NOT annotation-snap). 0
    # disables (byte-identical). isoquant uses Δ=6 for ONT.
    denovo_wobble_tol: int = 0
    # Shadow-ratio guard for the wobble merge: absorb a wobble-matching novel
    # candidate into a rep ONLY if it is a true SHADOW — len(cand.reads) <=
    # ratio * len(rep.reads). Protects genuine CLOSE isoforms (real alt
    # donor/acceptor a few bp apart, e.g. NAGNAG) from being merged: a minimap
    # shadow has few reads, a real close isoform has comparable support. 1.0
    # merges any wobble-match (naive; measured to eat correct transcripts at
    # tol=6). isoquant's ½-support bulge criterion ≈ 0.5. Only used when
    # denovo_wobble_tol > 0.
    denovo_wobble_shadow_ratio: float = 0.5
    # DE-NOVO intron-graph assembly (EXPERIMENT default-OFF). Attacks the measured
    # #1 de-novo error: truncation (dRNA 3'-bias -> reads cover only 3' junctions
    # -> pyfin emits truncated candidates, gffcompare class 'c'). When ON, pool all
    # reads' junctions, cluster wobbles to a read-count consensus (denovo_graph_tol
    # bp), build a read-adjacency graph, and EXTEND each read's chain through
    # UNAMBIGUOUSLY-supported edges (>= denovo_graph_min_edge_reads) to a maximal
    # chain — assembling truncated partials into full-length transcripts, stopping
    # at genuine branch points (no fabricated wrong combos). Reads are then grouped
    # by the extended chain; the candidate span is the union of its reads' extents
    # (5'/3' ends) and its sequence is built from the genome. NOT annotation-snap
    # (consensus + edges from reads only). OFF => byte-identical.
    denovo_graph: bool = False
    denovo_graph_tol: int = 6
    denovo_graph_min_edge_reads: int = 2
    # 5'-TSS brake for --denovo-graph. dRNA truncation is 5'-ward, so a truncated
    # read of a LONG transcript stops at a random 5' position while a COMPLETE read
    # of a genuine short isoform stops at that isoform's real TSS. A chain whose
    # 5'-terminal junction sits at a TSS peak (>= tss_frac of its reads pile their
    # 5'-end within tss_tol bp, and >= tss_min_reads reads) is NOT extended 5'-ward,
    # so short isoforms contained inside a longer transcript survive instead of
    # being merged away. Validated on p00: real short-isoform TSS carry 40-90% of a
    # locus's read-5'-ends; degradation-only starts scatter and never cross frac.
    denovo_graph_tss_brake: bool = True
    denovo_graph_tss_tol: int = 20
    denovo_graph_tss_min_reads: int = 3
    denovo_graph_tss_frac: float = 0.4
    # The mirror of Lever 2 for GTF-passthrough (source="gtf") candidates: drop a
    # guided multi-exon candidate if ANY of its junctions is spliced by fewer than
    # guided_junction_min_reads directly-observed reads (primary-read CIGAR introns,
    # strand-keyed, within guided_junction_reads_tol bp of BOTH boundaries). Novel
    # candidates keep going through Lever 2; fusion/mono are exempt.
    # WHY: Lever 2 + the coordinate-EXACT gates are novel-only, so a jitter-corrupted
    # GTF junction (coords shifted >tol from the true site) is never checked against
    # read CIGARs — it only faces the coordinate-INEXACT M2/M1 support gate, which a
    # ±10bp shift survives (it still holds reads' mappy sole-AS / ties the NLL). That
    # is the measured c_jitter precision collapse (94.8→41.3). Requiring exact read
    # support for GTF junctions attacks that directly; it also trims GTF echo
    # (orphan) in loci with observed splicing — but a FULLY read-sparse locus fails
    # open (observed map is None -> no drop), so echo there is unaffected.
    # RECALL RISK: this also drops genuine but low-coverage annotated junctions —
    # the classic real-vs-synthetic sign-flip knob. Keep min_reads small, tol tight,
    # and validate on real data across samples before any default flip. 0 disables
    # (byte-identical / recall-safe).
    guided_junction_min_reads: int = 0
    guided_junction_reads_tol: int = 2
    # Junction-dominance gate (PRE-EM, in process_interval after the canonical
    # gate). The "junction-first" idea: before quantification, drop a NOVEL
    # multi-exon candidate if ANY of its junctions is either (a) supported by
    # fewer than junction_dominance_min_reads directly-observed reads, OR (b) NOT
    # locally dominant — i.e. a DIFFERENT observed junction within
    # junction_dominance_window_bp (beyond junction_dominance_tol_bp) carries
    # strictly more reads. Unlike a pure read-count gate this also removes
    # multi-read wobble shadows (they lose to the stronger true junction a few bp
    # away). Pure mapping evidence (CIGAR introns, strand-keyed), no snap; runs
    # before EM so shadows never compete for reads. gtf/fusion/mono exempt.
    # OFF (default) -> no drops (byte-identical).
    junction_dominance_filter: bool = False
    junction_dominance_min_reads: int = 2
    junction_dominance_window_bp: int = 20
    junction_dominance_tol_bp: int = 2
    # M2/M1 read-support gate (quant_mode="m2_em" only). A multi-exon candidate
    # (GTF or novel; fusion/mono exempt) is KEPT iff it earns >=1 read's support:
    #   (a) it is some read's M1 SOLE best-AS (the read's tie set is exactly this
    #       candidate -> unique mappy-AS winner), OR
    #   (b) it is some read's M2 best (lowest junction-NLL among that read's tie).
    # Otherwise it is dropped. Rationale: a jitter-corrupted GTF junction loses the
    # M2 signal contest to the true-junction candidate AND never wins a read's sole
    # M1-AS, so it earns no support; a genuine isoform always wins either its own
    # sole-AS reads or the M2 contest. Catches corrupted-annotation shadows the
    # abundance levers miss. OFF == no drops. m2_support_gate_tie: when True (default)
    # a candidate tied for the lowest M2 NLL also counts as M2-best (recall-safer);
    # False requires a strict unique M2 win.
    # Default ON (tie-accept) tuned on the SGNex heya8+sirv4 12-condition matrix:
    # lifts the corrupted-annotation (c_jitter) F1@3 with ZERO recall loss anywhere
    # (full Sn@3 stays 100). The strict (tie=False) variant drops more shadows but
    # costs full-condition recall (Sn@3 -> 98.8), so tie-accept is the default.
    m2_support_gate: bool = True
    m2_support_gate_tie: bool = True
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

    # R-matrix persistence (T8: FP-by-EM data infrastructure). INERT: no writer path
    # exists under fin/ (the R.npy/R_meta.json dump was never wired into the assembly
    # runner); retained for CLI/constructor back-compat. Do not delete.
    persist_R_matrix: bool = True  # (inert) no per-interval R dump is written

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

    @classmethod
    def from_profile(cls, profile: str, **kwargs) -> "PipelineConfig":
        """Construct a programmatic config from a named operating point.

        Explicit keyword arguments override profile values, matching the CLI.
        """
        if "profile" in kwargs:
            raise TypeError("profile must be passed as the first argument")
        config = cls(profile=profile, **kwargs)
        values = {key: getattr(config, key) for key in PROFILE_FIELDS}
        explicit = PROFILE_FIELDS.intersection(kwargs)
        resolved = resolve_profile_values(profile, values, explicit)
        for key, value in resolved.items():
            setattr(config, key, value)
        config.m2_metric, config.m2_metric_route = resolve_m2_metric(
            config.m2_metric, config.gtf_path
        )
        config.profile_overrides = tuple(sorted(explicit))
        return config

    def validate(self):
        """Validate required paths, profiles, metrics, and parallelism knobs."""
        for attr in (
            "bam_path",
            "gtf_path",
            "genome_fasta_path",
            "fastq_path",
            "signal_path",
        ):
            val = getattr(self, attr)
            if val and not Path(val).exists():
                raise FileNotFoundError(f"{attr}: {val}")
        if self.profile not in PIPELINE_PROFILES:
            raise ValueError(f"unknown profile: {self.profile!r}")
        if self.m2_metric == "auto":
            self.m2_metric, self.m2_metric_route = resolve_m2_metric(
                self.m2_metric, self.gtf_path
            )
        if self.isoform_fraction_locus not in ("family", "overlap"):
            raise ValueError(
                "isoform_fraction_locus must be 'family' or 'overlap'"
            )
        if getattr(self, "ranking_mode", "off") not in ("off", "filter"):
            raise ValueError(f"unknown ranking_mode: {self.ranking_mode!r}")
        if getattr(self, "tss_evidence_mode", "off") not in (
            "off", "audit", "require"
        ):
            raise ValueError(
                f"unknown tss_evidence_mode: {self.tss_evidence_mode!r}"
            )
        if getattr(self, "tss_evidence_mode", "off") != "off" and not getattr(
            self, "endpoint_refine", False
        ):
            raise ValueError(
                "tss_evidence_mode requires --endpoint-refine: the TSS test "
                "scores the alternative starts that EndpointRefine proposes"
            )
        if getattr(self, "endpoint_refine", False) and not (
            self.post_selection_refit and self.quant_mode == "m2_em"
            and not self.abundance_feedback
        ):
            raise ValueError(
                "endpoint_refine requires the post-selection refit "
                "(quant_mode=m2_em, post_selection_refit on, no "
                "abundance_feedback): post-split requantification is "
                "mandatory, splitting without it is not supported"
            )
        if self.m2_metric not in ("off", "mean", "summed_llr", "sqrt_count_mean_llr"):
            raise ValueError(f"unknown m2_metric: {self.m2_metric!r}")

        self.post_selection_refit_effective = False
        self.post_selection_refit_disable_reason = None
        if self.post_selection_refit:
            if self.quant_mode != "m2_em":
                self.post_selection_refit_disable_reason = (
                    f"unsupported quant_mode={self.quant_mode}; requires m2_em"
                )
                logger.warning(
                    "Post-selection abundance refit disabled: %s",
                    self.post_selection_refit_disable_reason,
                )
            elif self.abundance_feedback:
                profile_implied = (
                    self.profile != "custom"
                    and "post_selection_refit" not in self.profile_overrides
                )
                if not profile_implied:
                    raise ValueError(
                        "post_selection_refit is not equivalent to rerunning "
                        "abundance-feedback EM; disable one of them"
                    )
                self.post_selection_refit_disable_reason = (
                    "abundance_feedback requires a full EM rerun"
                )
                logger.warning(
                    "Post-selection abundance refit disabled: %s",
                    self.post_selection_refit_disable_reason,
                )
            else:
                self.post_selection_refit_effective = True

        # Inert-parameter transparency (review Medium 7): the production
        # m2_em path is a single beta=0 softmax with sigma hardcoded to 1.0,
        # so EM iteration/sigma knobs cannot change its output. Non-default
        # values are accepted for the other quant modes but must not be
        # silently ignored here.
        if self.quant_mode == "m2_em" and not self.abundance_feedback:
            inert = []
            if self.em_sigma != 1.0:
                inert.append(f"em_sigma={self.em_sigma} (m2_em fixes sigma=1.0)")
            if self.em_max_iter != 1000:
                inert.append(f"em_max_iter={self.em_max_iter}")
            if self.em_tol != 1e-4:
                inert.append(f"em_tol={self.em_tol}")
            if self.em_max_iter_override is not None:
                inert.append(
                    f"em_max_iter_override={self.em_max_iter_override}"
                )
            if inert:
                logger.warning(
                    "Inert under quant_mode=m2_em without abundance_feedback "
                    "(single beta=0 softmax; values have no effect): %s",
                    "; ".join(inert),
                )

        if self.m2_diff_cover_margin < 0 or self.m2_summed_llr_margin < 0:
            raise ValueError("M2 margins must be >= 0")
        if self.m2_summed_llr_flank < 1:
            raise ValueError("m2_summed_llr_flank must be >= 1")
        if self.junction_snap_tolerance < 1:
            raise ValueError("junction_snap_tolerance must be >= 1")
        if self.junction_snap_min_support < 1:
            raise ValueError("junction_snap_min_support must be >= 1")
        if self.junction_snap_min_ratio < 1.0:
            raise ValueError("junction_snap_min_ratio must be >= 1")
        if self.threads < 1:
            raise ValueError(f"threads must be >= 1, got {self.threads}")
        if not 0 <= self.gpu_workers <= self.threads:
            raise ValueError(
                f"gpu_workers must be in [0, threads={self.threads}], got {self.gpu_workers}"
            )
        if self.gpu_workers and not self.use_gpu:
            raise ValueError("gpu_workers must be 0 when use_gpu is False")
