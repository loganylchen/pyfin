"""Post-EM recovery of folded 5'-truncation shadows on the production m2_em path.

Discovery folds every EXACT sub-chain into its longest container
(``collapse``), keeping the folded chain and its reads as a dormant *shadow*
(``CandidateSet.shadows``). Whether such a shadow is a genuine short isoform
(its own TSS) or 5' degradation of the container is NOT decided at generation:
it is decided here, after a first assignment/EM pass, by testing the shadow's
proposed TSS against the parent's read-5'-end hazard profile
(:func:`fin.analysis.tss_evidence.evaluate_internal_tss`).

Contract (Codex gpt-5.6-sol spec, 2026-09-03):

- eligibility: parent present in the first-pass candidate axis with soft
  abundance > 0; shadow is a pure 5' truncation (exact suffix of the parent
  chain on ``+``, exact prefix on ``-``); >= 80 % of the shadow reads end
  within +-25 bp of the parent 3' end (else it is a potential own-TES model,
  not this consumer's business); >= 3 distinct shadow reads.
- TSS proposal = median strand-aware genomic 5' end of the shadow reads,
  mapped to the parent transcript frame.
- offsets under test = 5' offsets of the reads with first-pass responsibility
  ``R[read, parent] >= 0.05`` (read-level, unweighted, one row per read).
- background = pooled hazard of the parent's own profile (first bin
  excluded); the test itself also applies the local-neighbourhood floor.
- ``n_eligible_bins = parent_spliced_length // 25``: the shadow TSS was
  found by the same reads, so the selection correction must cover the whole
  searchable parent, not a single a-priori bin.
- only ``verdict == "supported"`` re-emits; the recovered candidate carries
  the shadow's own chain, a 5' endpoint at the median of the IN-BIN shadow
  reads, the parent's 3' endpoint, ``source="novel"``, the parent's
  ``family_id``, and ONLY the in-bin shadow reads as discovery support.
  Complete-chain depths are computed from the raw CIGAR indel records of
  those reads, never from EM responsibilities.
- the caller then appends the candidate(s) to the candidate set and reruns
  assignment/EM ONCE; abundance and read ownership are decided by that second
  pass (no manual mass transfer). No recursion.

Everything is recorded to ``tss_shadow_recovery.tsv`` (one row per shadow
examined) regardless of mode, so ``audit`` and ``require`` differ only in
whether the candidate is added.
"""
from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from fin.analysis.tss_evidence import (
    VERDICT_SUPPORTED,
    build_hazard_profile,
    evaluate_internal_tss,
    genomic_to_offset,
    pooled_background_hazard,
    read_five_prime_offset,
    spliced_length,
)
from fin.candidates.dataclasses import IntronChain, TranscriptCandidate

TSV_FIELDS = (
    "interval", "parent_id", "shadow_chain", "tss_genomic", "tss_offset",
    "n_shadow_reads", "n_at_risk", "n_peak", "n_upstream", "peak_fraction",
    "background_hazard", "peak_mad", "p_value", "verdict", "reason",
    "reemitted_candidate_id",
)
SHARD_DIR = "_tss_recovery_shards"
OUT_NAME = "tss_shadow_recovery.tsv"

BIN_BP = 25
MIN_SHADOW_READS = 3
TES_TOL_BP = 25
TES_AGREE_FRACTION = 0.80
MIN_RESPONSIBILITY = 0.05


@dataclass
class ShadowRecoveryRow:
    interval: str
    parent_id: str
    shadow_chain: str
    tss_genomic: int = -1
    tss_offset: int = -1
    n_shadow_reads: int = 0
    n_at_risk: int = 0
    n_peak: int = 0
    n_upstream: int = 0
    peak_fraction: float = 0.0
    background_hazard: float = 0.0
    peak_mad: float = -1.0
    p_value: Optional[float] = None
    verdict: str = "unidentifiable"
    reason: str = ""
    reemitted_candidate_id: str = ""

    def as_tsv(self) -> str:
        vals = []
        for k in TSV_FIELDS:
            v = getattr(self, k)
            if v is None:
                vals.append("")
            elif isinstance(v, float):
                vals.append(f"{v:.6g}")
            else:
                vals.append(str(v))
        return "\t".join(vals)


@dataclass
class ShadowRecoveryResult:
    candidates: List[TranscriptCandidate] = field(default_factory=list)
    # recovered candidate_id -> parent candidate_id
    parent_of: Dict[str, str] = field(default_factory=dict)
    rows: List[ShadowRecoveryRow] = field(default_factory=list)
    # candidate_id -> read ids to add to an ALREADY EXISTING candidate's discovery
    # support. Returned rather than applied: the evaluator never mutates the
    # candidate set, so `audit` stays byte-inert.
    support_merges: Dict[str, set] = field(default_factory=dict)
    # frozen read-5'-end start models for the second assignment pass, one per
    # recovered candidate: parent_id, recovered_id, tss_bin, bin_bp, hazard
    # (the parent's first-pass hazard of terminating in that bin, i.e. the
    # background the test beat), parent exons/strand and the read spans needed
    # to place a read's 5' end in the parent frame. Consumed by
    # fin.pipeline.assignment.apply_start_likelihood.
    start_models: List[dict] = field(default_factory=list)


def _chain_str(chain: Sequence[Tuple[int, int]]) -> str:
    return ";".join(f"{d}-{a}" for d, a in chain) or "."


def _is_pure_5p_truncation(shadow: Tuple, parent: Tuple, strand: str) -> bool:
    n = len(shadow)
    if n == 0 or n >= len(parent):
        return False
    if strand == "-":
        return tuple(parent[:n]) == tuple(shadow)
    return tuple(parent[len(parent) - n:]) == tuple(shadow)


def _median(vals: Sequence[int]) -> int:
    s = sorted(vals)
    return int(s[len(s) // 2])


def _exons(cand: TranscriptCandidate) -> List[Tuple[int, int]]:
    from fin.candidates.discovery import _exons_from_chain
    return _exons_from_chain(cand.start, cand.end, cand.intron_chain)


def recover_folded_tss_candidates(
    candidate_set,
    qo,
    interval,
    genome_seq: str,
    *,
    indels_by_read: Optional[Dict[str, tuple]] = None,
) -> ShadowRecoveryResult:
    """Test every folded shadow of every first-pass parent; build the
    recovered candidates for the ``supported`` ones. Pure: mutates neither
    ``candidate_set`` nor ``qo``."""
    out = ShadowRecoveryResult()
    shadows_map = getattr(candidate_set, "shadows", None) or {}
    spans = getattr(candidate_set, "read_spans", None) or {}
    if not shadows_map or qo is None:
        return out
    indels_by_read = indels_by_read or {}
    ivl_key = f"{interval.chrom}:{interval.start}-{interval.end}:{interval.strand}"

    col_of = {c.candidate_id: j for j, c in enumerate(qo.cand_list)}
    row_of = {rid: i for i, rid in enumerate(qo.kept_read_ids)}
    existing_ids = {c.candidate_id for c in candidate_set.candidates}
    # (chain, tss_bin, parent) -> recovered candidate, so several shadows that
    # propose the same thing yield ONE candidate with merged support.
    built: Dict[Tuple[tuple, int, str], TranscriptCandidate] = {}

    for parent_id in sorted(shadows_map):
        shadows = shadows_map[parent_id]
        parent = candidate_set.get_candidate(parent_id)
        j = col_of.get(parent_id)
        for schain, sreads in shadows:
            row = ShadowRecoveryRow(
                interval=ivl_key, parent_id=parent_id,
                shadow_chain=_chain_str(schain), n_shadow_reads=len(set(sreads)),
            )
            out.rows.append(row)
            if parent is None or j is None:
                row.reason = "parent_not_in_first_pass"
                continue
            parent_mass = float(qo.R[:, j].sum()) if qo.R.size else 0.0
            if parent_mass <= 0.0:
                row.reason = "parent_zero_abundance"
                continue
            strand = parent.strand
            pchain = tuple(parent.intron_chain.introns)
            if not _is_pure_5p_truncation(tuple(schain), pchain, strand):
                row.reason = "not_pure_5p_truncation"
                continue
            reads = sorted({r for r in set(sreads) if r in spans})
            if len(reads) < MIN_SHADOW_READS:
                row.reason = f"too_few_shadow_reads({len(reads)}<{MIN_SHADOW_READS})"
                continue
            parent_3p = parent.end if strand == "+" else parent.start
            n_tes = sum(
                1 for r in reads
                if abs((spans[r][1] if strand == "+" else spans[r][0]) - parent_3p)
                <= TES_TOL_BP
            )
            if n_tes < TES_AGREE_FRACTION * len(reads):
                row.reason = f"tes_disagrees_with_parent({n_tes}/{len(reads)})"
                continue
            p5 = {r: (spans[r][0] if strand == "+" else spans[r][1]) for r in reads}
            tss_genomic = _median(list(p5.values()))
            row.tss_genomic = tss_genomic
            exons = _exons(parent)
            tx_len = spliced_length(exons)
            tss_offset = genomic_to_offset(tss_genomic, exons, strand)
            if tss_offset is None:
                row.reason = "tss_outside_parent_exons"
                continue
            row.tss_offset = int(tss_offset)

            # first-pass parent reads (responsibility >= 0.05), one offset per read
            offsets: List[int] = []
            for rid, i in row_of.items():
                if qo.R[i, j] < MIN_RESPONSIBILITY or rid not in spans:
                    continue
                off = read_five_prime_offset(spans[rid], exons, strand)
                if off is not None:
                    offsets.append(int(off))
            if not offsets:
                row.reason = "no_parent_offsets"
                continue
            profile = build_hazard_profile(offsets, tx_len, bin_bp=BIN_BP)
            background = pooled_background_hazard([profile], exclude_first_bins=1)
            ev = evaluate_internal_tss(
                candidate_id=f"{parent_id}|{_chain_str(schain)}",
                parent_id=parent_id, offsets=offsets, tss_offset=int(tss_offset),
                background_hazard=background, bin_bp=BIN_BP, min_peak_reads=3,
                min_peak_fraction=0.10, min_at_risk=10, alpha=0.05,
                max_peak_mad=12.0, n_bootstrap=2000, use_local_background=True,
                neighbourhood_bins=12, n_eligible_bins=max(1, tx_len // BIN_BP),
            )
            row.n_at_risk, row.n_peak, row.n_upstream = ev.n_at_risk, ev.n_peak, ev.n_upstream_of_tss
            row.peak_fraction, row.background_hazard = ev.peak_fraction, ev.background_hazard
            row.peak_mad, row.p_value = ev.peak_mad, ev.p_value
            row.verdict, row.reason = ev.verdict, ev.reason
            if ev.verdict != VERDICT_SUPPORTED:
                continue

            # --- build the recovered candidate from the IN-BIN shadow reads ----
            b0 = int(tss_offset) // BIN_BP
            lo, hi = b0 * BIN_BP, (b0 + 1) * BIN_BP - 1
            in_bin = []
            for r in reads:
                off = read_five_prime_offset(spans[r], exons, strand)
                if off is not None and lo <= off <= hi:
                    in_bin.append(r)
            if not in_bin:
                row.verdict, row.reason = "unidentifiable", "no_shadow_reads_in_supported_bin"
                continue
            five = _median([p5[r] for r in in_bin])
            if strand == "+":
                start, end = five, parent.end
            else:
                start, end = parent.start, five
            key = (tuple(schain), b0, parent_id)
            prev = built.get(key)
            if prev is not None:
                prev.supporting_read_ids.update(in_bin)
                prev.raw_supporting_read_ids.update(in_bin)
                row.reemitted_candidate_id = prev.candidate_id
                continue
            from fin.candidates.discovery import (
                _build_spliced_sequence, _generate_novel_id, junction_clean_depths,
            )
            chain = IntronChain(introns=tuple(schain))
            if end <= start:
                row.verdict, row.reason = "unidentifiable", "empty_span"
                continue
            seq = _build_spliced_sequence(genome_seq, start, end, chain, strand)
            if not seq:
                row.verdict, row.reason = "unidentifiable", "no_sequence"
                continue
            cid = _generate_novel_id(parent.chrom, strand, start, end, chain)
            row.reemitted_candidate_id = cid
            if cid in existing_ids:
                # Already a candidate. The evaluator is PURE -- it must not touch
                # the candidate set, in audit mode or any other, or the audit
                # stops being output-inert (support counts feed the
                # containment-cluster and support gates). The merge is returned
                # instead, for the caller to apply only in `require` mode.
                out.support_merges.setdefault(cid, set()).update(in_bin)
                row.reason += ";existing_candidate"
                continue
            depths = junction_clean_depths(chain, in_bin, (), indels_by_read)
            cand = TranscriptCandidate(
                candidate_id=cid, intron_chain=chain,
                three_prime_pos=(end if strand == "+" else start), sequence=seq,
                source="novel", supporting_read_ids=set(in_bin),
                raw_supporting_read_ids=set(in_bin),
                chrom=parent.chrom, strand=strand, start=start, end=end,
                family_id=parent.family_id,
                complete_chain_junction_depths=tuple(depths),
                complete_chain_weakest_indel_clean_depth=(min(depths) if depths else 0),
            )
            built[key] = cand
            existing_ids.add(cid)
            out.candidates.append(cand)
            out.parent_of[cid] = parent_id
            out.start_models.append(dict(
                parent_id=parent_id, recovered_id=cid, tss_bin=b0, bin_bp=BIN_BP,
                hazard=float(ev.background_hazard), exons=list(exons), strand=strand,
                read_spans=spans,
            ))
    return out


def write_recovery_shard(work_dir, interval, rows: List[ShadowRecoveryRow]) -> None:
    if not rows:
        return
    d = Path(work_dir) / SHARD_DIR
    d.mkdir(parents=True, exist_ok=True)
    name = f"{interval.chrom}_{interval.start}_{interval.end}_{interval.strand or 'u'}.tsv"
    tmp = d / (name + ".tmp")
    with tmp.open("w") as fh:
        for r in rows:
            fh.write(r.as_tsv() + "\n")
    tmp.replace(d / name)


def merge_recovery_shards(work_dir) -> Optional[str]:
    wd = Path(work_dir)
    shard_dir = wd / SHARD_DIR
    rows = []
    if shard_dir.is_dir():
        for shard in sorted(shard_dir.glob("*.tsv")):
            with shard.open() as fh:
                rows.extend(line.rstrip("\n") for line in fh if line.strip())
    if not rows:
        if shard_dir.is_dir():
            shutil.rmtree(shard_dir)
        return None
    rows.sort(key=lambda l: (l.split("\t")[0], l.split("\t")[1], l.split("\t")[2]))
    out = wd / OUT_NAME
    tmp = out.with_suffix(".tsv.tmp")
    with tmp.open("w") as fh:
        fh.write("\t".join(TSV_FIELDS) + "\n")
        for line in rows:
            fh.write(line + "\n")
    tmp.replace(out)
    shutil.rmtree(shard_dir)
    return str(out)
