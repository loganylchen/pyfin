"""Containment decision layer: is a contained (nested) candidate real?

For every SHORT candidate whose structure nests inside a LONGER candidate of
the same interval (pairs discovered from inference-time candidates only,
never from truth), route the question through the degradation-aware evidence
ladder of :mod:`fin.analysis.tss_evidence`:

- distinct 3' end            -> :func:`evaluate_tes_support` (dRNA's reliable
                                end: a distinct TES cluster is positive
                                evidence, its absence negative);
- same TES, internal TSS     -> :func:`evaluate_internal_tss` (conditional
                                termination hazard against the candidate's
                                local degradation null);
- mono-exon short model      -> the same two routes over the parent's
                                transcript frame using read spans (chain
                                evidence is undefined for mono models).

Verdicts are tri-state (``supported`` / ``unsupported`` -> reported here as
``contradicted`` / ``unidentifiable``) and abstention always retains. A
short model may be CONTRADICTED only when the long parent is independently
supported: every parent junction absent from the short model must carry at
least ``PARENT_MIN_UNIQUE_CLEAN_READS`` indel-clean exact complete-chain
reads (see discovery's ``complete_chain_junction_depths``); without that,
the degradation explanation has no observed parent and the pair abstains.

Leakage rules: every input is observable at inference time (candidate
structures, read spans, EM assignments). NanoCount, gffcompare classes, and
truth GTFs never enter. ``p_value`` is populated only by the hazard rung;
TES/structural rules stay outside the probability family, mirroring the TSS
module contract. q-values are attached by the caller across the whole run
(audit metadata; behavioral verdicts use the per-test abstention-aware
machinery, not the FDR pass).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from fin.analysis.tss_evidence import (
    TssEvidence,
    evaluate_internal_tss,
    evaluate_tes_support,
    genomic_to_offset,
    spliced_length,
)

logger = logging.getLogger(__name__)

PARENT_MIN_UNIQUE_CLEAN_READS = 2
TES_WINDOW_BP = 25

VERDICT_SUPPORTED = "supported"
VERDICT_CONTRADICTED = "contradicted"
VERDICT_UNIDENTIFIABLE = "unidentifiable"


@dataclass
class ContainmentRow:
    """One short-in-long containment decision (audit artifact row)."""

    interval: str
    short_id: str
    parent_id: str
    short_source: str
    parent_source: str
    rung: str                      # own_tes | tss_only | equal_chain | mono_span
    verdict: str
    reason: str
    # structural_tes rows act at interval level (reliable 3' end, outside
    # the probability family); hazard rows act only after the global
    # grouped-FDR pass (p_within + BY q before any drop/exempt).
    test_family: str = ""          # structural_tes | hazard | ""
    action: str = "none"           # none | drop | exempt
    parent_unique_min_clean_depth: int = -1
    short_weakest_clean_depth: int = 0
    n_pool_reads: int = 0
    n_peak: int = 0
    n_at_risk: int = 0
    peak_fraction: float = 0.0
    background_hazard: float = 0.0
    p_value: Optional[float] = None
    q_value: Optional[float] = None

    def as_tsv_row(self) -> List[str]:
        def fmt(v):
            if v is None:
                return ""
            if isinstance(v, float):
                return f"{v:.6g}"
            return str(v)

        return [fmt(getattr(self, f)) for f in TSV_FIELDS]


TSV_FIELDS = [
    "interval", "short_id", "parent_id", "short_source", "parent_source",
    "rung", "verdict", "reason", "test_family", "action",
    "parent_unique_min_clean_depth", "short_weakest_clean_depth",
    "n_pool_reads", "n_peak", "n_at_risk", "peak_fraction",
    "background_hazard", "p_value", "q_value",
]


def _chain(cand) -> Tuple[Tuple[int, int], ...]:
    ic = getattr(cand, "intron_chain", None)
    return tuple(getattr(ic, "introns", ()) or ())


def _exons_from(cand) -> Tuple[Tuple[int, int], ...]:
    """0-based half-open exon blocks reconstructed from span + chain."""
    chain = _chain(cand)
    start, end = int(cand.start), int(cand.end)
    if not chain:
        return ((start, end),)
    exons = [(start, chain[0][0])]
    for i in range(1, len(chain)):
        exons.append((chain[i - 1][1], chain[i][0]))
    exons.append((chain[-1][1], end))
    return tuple(exons)


def _is_contig_subchain(short, long) -> Optional[int]:
    """Offset of ``short`` inside ``long`` as a contiguous block, else None."""
    n, m = len(short), len(long)
    if n == 0 or n >= m:
        return None
    for off in range(m - n + 1):
        if long[off:off + n] == short:
            return off
    return None


def _classify_pair(short, long) -> Optional[str]:
    """Containment rung for a candidate pair (None = not contained).

    Rungs: ``own_tes`` (distinct 3' end), ``tss_only`` (sub-chain, same
    TES, internal TSS), ``equal_chain`` (IDENTICAL chain, same TES,
    strictly interior 5' start - the central internal-TSS case), and
    ``mono_span`` (chainless short nested in the parent span).
    """
    cs, cl = _chain(short), _chain(long)
    strand = long.strand
    if not cl:
        return None                     # mono parent never a containment frame
    if not cs:
        # mono short model: span must nest inside the parent's genomic span
        if long.start <= short.start and short.end <= long.end:
            return "mono_span"
        return None
    equal = cs == cl
    off = None if equal else _is_contig_subchain(cs, cl)
    if not equal and off is None:
        return None
    tes_short = short.end if strand == "+" else short.start
    tes_long = long.end if strand == "+" else long.start
    exl = _exons_from(long)
    o_s = genomic_to_offset(
        tes_short - 1 if strand == "+" else tes_short, exl, strand)
    o_l = genomic_to_offset(
        tes_long - 1 if strand == "+" else tes_long, exl, strand)
    if o_s is None or o_l is None:
        return None
    if abs(o_s - o_l) > TES_WINDOW_BP:
        return "own_tes"
    if not equal:
        return "tss_only"
    # Equal chain + shared TES: containment is by the 5' terminus alone.
    g5s = short.start if strand == "+" else short.end - 1
    g5l = long.start if strand == "+" else long.end - 1
    d_s = genomic_to_offset(g5s, exl, strand)
    d_l = genomic_to_offset(g5l, exl, strand)
    if d_s is None or d_l is None:
        return None
    if d_s - d_l > TES_WINDOW_BP:       # strictly interior 5' start
        return "equal_chain"
    return None


def _parent_unique_min_depth(short, long) -> int:
    """Minimum indel-clean depth over the parent's junctions NOT shared with
    the short model (its unique structure). -1 = no evidence collected /
    no unique junction (equal-chain pairs use parent SPAN evidence instead:
    reads extending 5' beyond the short start - see evaluate_containment)."""
    depths = tuple(getattr(long, "complete_chain_junction_depths", ()) or ())
    cl = _chain(long)
    if not depths or len(depths) != len(cl):
        return -1
    cs = set(_chain(short))
    unique = [d for j, d in zip(cl, depths) if j not in cs]
    if not unique:
        return -1                        # no unique junction (same chain)
    return min(unique)


def _pool_offsets(short, long, quant_by_id, read_spans, strand):
    """5' and 3' read-end offsets (long-parent frame) of the pair's reads."""
    exl = _exons_from(long)
    tx_len = spliced_length(exl)
    ids: List[str] = []
    for cand in (short, long):
        qr = quant_by_id.get(cand.candidate_id)
        if qr is not None:
            ids.extend(qr.assigned_read_ids or ())
    off5: List[int] = []
    off3: List[int] = []
    seen = set()
    for rid in ids:
        if rid in seen or rid not in read_spans:
            continue
        seen.add(rid)
        s, e = read_spans[rid]
        five = s if strand == "+" else e - 1
        three = e - 1 if strand == "+" else s
        o5 = genomic_to_offset(five, exl, strand)
        if o5 is None:
            # read 5' beyond the parent's 5' terminus -> offset 0
            if (strand == "+" and s < exl[0][0]) or (
                strand == "-" and e > exl[-1][1]
            ):
                o5 = 0
        o3 = genomic_to_offset(three, exl, strand)
        if o5 is not None:
            off5.append(int(o5))
        if o3 is not None:
            off3.append(int(o3))
    return off5, off3, tx_len, len(seen)


def evaluate_containment(
    cand_list: Sequence,
    quant_results: Sequence,
    read_spans: Dict[str, Tuple[int, int]],
    *,
    interval_name: str = "",
) -> List[ContainmentRow]:
    """Evaluate every contained candidate pair of one interval.

    Deterministic: pairs are scanned in candidate-list order and each short
    model is judged against its LONGEST containing parent (most junctions,
    then longest spliced span, then lexicographic id).
    """
    quant_by_id = {q.candidate_id: q for q in quant_results}
    rows: List[ContainmentRow] = []
    multi = [c for c in cand_list if c.source != "fusion"]
    for short in multi:
        if short.candidate_id not in quant_by_id:
            continue
        best = None
        for long in multi:
            if long is short or long.strand != short.strand:
                continue
            if long.candidate_id not in quant_by_id:
                continue
            rung = _classify_pair(short, long)
            if rung is None:
                continue
            key = (len(_chain(long)), spliced_length(_exons_from(long)),
                   long.candidate_id)
            if best is None or key > best[0]:
                best = (key, long, rung)
        if best is None:
            continue
        _, long, rung = best
        strand = long.strand
        row = ContainmentRow(
            interval=interval_name, short_id=short.candidate_id,
            parent_id=long.candidate_id, short_source=short.source,
            parent_source=long.source, rung=rung,
            verdict=VERDICT_UNIDENTIFIABLE, reason="not_evaluated",
            short_weakest_clean_depth=int(getattr(
                short, "complete_chain_weakest_indel_clean_depth", 0)),
        )
        rows.append(row)
        pmin = _parent_unique_min_depth(short, long)
        row.parent_unique_min_clean_depth = pmin
        span_evidence_rungs = ("equal_chain", "mono_span")
        if rung not in span_evidence_rungs and (
            pmin < PARENT_MIN_UNIQUE_CLEAN_READS
        ):
            row.reason = (
                f"parent_unsupported(unique_clean_depth={pmin}<"
                f"{PARENT_MIN_UNIQUE_CLEAN_READS})"
            )
            continue
        off5, off3, tx_len, n_pool = _pool_offsets(
            short, long, quant_by_id, read_spans, strand)
        row.n_pool_reads = n_pool
        if not off5:
            row.reason = "no_mappable_reads"
            continue
        exl = _exons_from(long)
        if rung == "own_tes" or (rung == "mono_span" and _mono_tes_distinct(
                short, long, exl, strand)):
            g3 = short.end - 1 if strand == "+" else short.start
            dc = genomic_to_offset(g3, exl, strand)
            if dc is None:
                row.reason = "short_tes_outside_parent"
                continue
            row.test_family = "structural_tes"
            ev = evaluate_tes_support(
                candidate_id=short.candidate_id,
                parent_id=long.candidate_id,
                read_three_prime_offsets=off3,
                candidate_tes_offset=int(dc),
                parent_tes_offset=tx_len - 1,
                window_bp=TES_WINDOW_BP,
            )
        else:
            row.test_family = "hazard"
            g5 = short.start if strand == "+" else short.end - 1
            d0 = genomic_to_offset(g5, exl, strand)
            if d0 is None:
                row.reason = "short_tss_outside_parent"
                continue
            if d0 // 25 == 0:
                row.reason = "short_start_equals_parent_start"
                continue
            # Selection correction: a GTF-provided start is an a-priori
            # offset (one eligible bin). A discovered (novel) start was
            # chosen from the data, so the honest search space is every bin
            # of the parent frame.
            n_bins = (1 if short.source == "gtf" and not getattr(short, "extent_of", None)
                      else max(1, tx_len // 25))
            ev = evaluate_internal_tss(
                candidate_id=short.candidate_id,
                parent_id=long.candidate_id,
                offsets=off5,
                tss_offset=int(d0),
                background_hazard=0.0,
                identifiability=rung,
                n_eligible_bins=n_bins,
            )
            # Independent parent evidence for span-based rungs: the
            # degradation (parent) explanation requires reads extending 5'
            # beyond the short start. Without >= PARENT_MIN qualified
            # upstream reads a CONTRADICTION is unfounded -> abstain.
            # (A 'supported' spike stands on its own evidence.)
            if (
                rung in span_evidence_rungs
                and ev.verdict == "unsupported"
                and ev.n_upstream_of_tss < PARENT_MIN_UNIQUE_CLEAN_READS
            ):
                row.reason = (
                    f"parent_span_unsupported(upstream={ev.n_upstream_of_tss}"
                    f"<{PARENT_MIN_UNIQUE_CLEAN_READS})"
                )
                row.n_peak = ev.n_peak
                row.n_at_risk = ev.n_at_risk
                continue
        _absorb(row, ev)
    return rows


def _mono_tes_distinct(short, long, exl, strand) -> bool:
    g3 = short.end - 1 if strand == "+" else short.start
    dc = genomic_to_offset(g3, exl, strand)
    if dc is None:
        return False
    return abs(int(dc) - (spliced_length(exl) - 1)) > TES_WINDOW_BP


def _absorb(row: ContainmentRow, ev: TssEvidence) -> None:
    row.reason = ev.reason
    row.n_peak = ev.n_peak
    row.n_at_risk = ev.n_at_risk
    row.peak_fraction = round(ev.peak_fraction, 6)
    row.background_hazard = round(ev.background_hazard, 6)
    row.p_value = ev.p_value
    if ev.verdict == "supported":
        row.verdict = VERDICT_SUPPORTED
    elif ev.verdict == "unsupported":
        row.verdict = VERDICT_CONTRADICTED
    else:
        row.verdict = VERDICT_UNIDENTIFIABLE


def containment_interval_actions(
    rows: Sequence[ContainmentRow], mode: str
) -> Tuple[set, set]:
    """Interval-time (drop_ids, exempt_ids) under ``mode``.

    audit -> nothing. require -> STRUCTURAL TES rows act immediately
    (contradicted drops, supported exempts: the reliable 3' end sits outside
    the probability family). HAZARD rows never drop here - their final
    verdict needs the global grouped-FDR pass - but a provisionally
    supported hazard short is exempted from the competition gates so the
    corrected global decision still has a live candidate to rule on.
    Abstention never acts.
    """
    if mode != "require":
        return set(), set()
    drops = set()
    exempts = set()
    for r in rows:
        if r.test_family == "structural_tes":
            if r.verdict == VERDICT_CONTRADICTED:
                r.action = "drop"
                drops.add(r.short_id)
            elif r.verdict == VERDICT_SUPPORTED:
                r.action = "exempt"
                exempts.add(r.short_id)
        elif r.test_family == "hazard" and r.verdict == VERDICT_SUPPORTED:
            exempts.add(r.short_id)      # provisional; finalized globally
    return drops, exempts


def write_containment_tsv(path, rows: Sequence[ContainmentRow]) -> None:
    """Deterministic audit artifact (sorted by interval, short, parent)."""
    ordered = sorted(
        rows, key=lambda r: (r.interval, r.short_id, r.parent_id))
    with open(path, "w") as fh:
        fh.write("\t".join(TSV_FIELDS) + "\n")
        for r in ordered:
            fh.write("\t".join(r.as_tsv_row()) + "\n")


def finalize_containment(
    work_dir, aggregated: Dict[str, "object"], mode: str
) -> Tuple[Dict[str, "object"], List[str], Optional[str]]:
    """Global containment finalization: FDR-correct hazard rows, apply the
    corrected require-mode hazard actions, and write the merged artifact.

    Runs once in the parent process AFTER cross-interval aggregation and
    BEFORE the global selection cascade. Worker shards are re-sorted
    globally (completion order cannot perturb the artifact). Hazard rows
    get ``p_within`` + grouped Bonferroni-within-parent / BY-across-parents
    ``q_value``; a provisionally supported hazard short whose corrected
    verdict is demoted loses its exemption and (in require mode with parent
    evidence already established) is dropped; a still-supported one keeps
    the soft-mass exemption via ``QuantResult.containment_verdict``.
    Structural TES rows are untouched (outside the probability family;
    their actions already ran at interval level). ``audit`` never mutates
    ``aggregated``. Writes a header-only artifact when no pair exists.
    Returns (aggregated, dropped_ids, artifact_path).
    """
    import csv
    import shutil
    from pathlib import Path as _Path

    from fin.analysis.tss_evidence import apply_grouped_fdr

    shard_dir = _Path(work_dir) / "_containment_shards"
    rows: List[dict] = []
    if shard_dir.is_dir():
        for shard in sorted(shard_dir.glob("*.tsv")):
            with shard.open() as fh:
                rows.extend(csv.DictReader(fh, delimiter="\t"))
    rows.sort(key=lambda r: (r["interval"], r["short_id"], r["parent_id"]))
    stubs = []
    for r in rows:
        if r.get("test_family") == "hazard" and r.get("p_value"):
            stub = TssEvidence(
                candidate_id=r["short_id"], parent_id=r["parent_id"],
                tss_offset=0, identifiability=r["rung"],
                verdict=("unsupported" if r["verdict"] == "contradicted"
                         else r["verdict"]),
                reason=r["reason"],
                p_value=float(r["p_value"]),
                selection_corrected=True,
            )
            stubs.append((stub, r))
    if stubs:
        apply_grouped_fdr([s for s, _ in stubs])
        for stub, r in stubs:
            if stub.q_value is not None:
                r["q_value"] = f"{stub.q_value:.6g}"
            # FDR can only demote supported -> unsupported; propagate.
            if r["verdict"] == VERDICT_SUPPORTED and (
                stub.verdict == "unsupported"
            ):
                r["verdict"] = VERDICT_CONTRADICTED
                r["reason"] = stub.reason
    # A short model is tested once per PARENT, so it can own several rows with
    # different verdicts. Resolve them by explicit PRECEDENCE, never by row
    # order (a dict comprehension over the rows keeps whichever row happened to
    # sort last, which silently discarded `supported` exemptions).
    #
    #   supported > contradicted > unsupported > unidentifiable
    #
    # `supported` outranks `contradicted` deliberately: the two verdicts come
    # from comparisons against DIFFERENT parents, and one comparison producing
    # positive evidence for the short model's own existence is not erased by
    # another comparison failing to reject it against a different parent. This
    # is the recall-safe reading of the project's red line that insufficient
    # evidence is not evidence of absence.
    _PRECEDENCE = {VERDICT_SUPPORTED: 3, VERDICT_CONTRADICTED: 2,
                   "unsupported": 1, VERDICT_UNIDENTIFIABLE: 0}
    final_by_id: Dict[str, str] = {}
    for r in rows:
        if r.get("test_family") != "hazard":
            continue
        cid, v = r["short_id"], r["verdict"]
        cur = final_by_id.get(cid)
        if cur is None or _PRECEDENCE.get(v, -1) > _PRECEDENCE.get(cur, -1):
            final_by_id[cid] = v

    dropped: List[str] = []
    if mode == "require":
        for r in rows:
            if r.get("test_family") != "hazard":
                continue
            cid = r["short_id"]
            resolved = final_by_id.get(cid, r["verdict"])
            if resolved == VERDICT_CONTRADICTED:
                r["action"] = "drop"
                if cid in aggregated:
                    del aggregated[cid]
                    dropped.append(cid)
            elif resolved == VERDICT_SUPPORTED:
                r["action"] = "exempt"
        # Final verdicts onto surviving QuantResults (soft-mass exemption
        # consumes 'supported'); deterministic single pass.
        for cid, qr in aggregated.items():
            v = final_by_id.get(cid)
            if v is not None:
                qr.containment_verdict = v
    out = _Path(work_dir) / "containment_evidence.tsv"
    tmp = out.with_suffix(".tsv.tmp")
    with tmp.open("w") as fh:
        fh.write("\t".join(TSV_FIELDS) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(f, "") for f in TSV_FIELDS) + "\n")
    tmp.replace(out)
    if shard_dir.is_dir():
        shutil.rmtree(shard_dir)
    return aggregated, sorted(dropped), str(out)
