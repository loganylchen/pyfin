"""Inference-time candidate evidence table (ranking substrate).

Emits one row per post-selection survivor with ONLY features that are
observable at inference time - no truth annotation, no expression oracle.
This is the shared feature layer for candidate-ranking work: the offline
model fit and any future in-pipeline scorer consume the same functions, so
train-time and inference-time features cannot drift.

The computation is strictly read-only over the selected results: enabling it
changes no candidate, no abundance, and no output row; it only writes an
additional ``candidate_evidence.tsv`` next to the ordinary outputs.
"""
from __future__ import annotations

import logging
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple

from fin.candidates.canonical import _revcomp2

logger = logging.getLogger(__name__)

EVIDENCE_COLUMNS = (
    "candidate_id", "chrom", "strand", "start", "end", "source",
    "n_exons", "tx_length", "is_mono",
    "abundance", "num_reads", "soft_hard_ratio", "confidence", "max_R",
    "family_id", "family_size", "family_share", "family_rank",
    "family_dominant_share",
    "is_subchain_of_sibling", "is_superchain_of_sibling",
    "weakest_junction_support", "median_junction_support",
    "n_junctions_below3",
    "canonical_fraction",
    "end5_support_frac", "end3_support_frac", "fulllen_frac",
    "n_end_reads",
)


@dataclass
class RankingBamEvidence:
    """Single-pass whole-BAM evidence for the ranking feature layer.

    ``complete`` is True only when the BAM iterator was exhausted without an
    error. The ranking filter must refuse to run on an incomplete scan -
    undercounted junction support would silently depress scores; the pure
    audit path may still use the partial maps for inspection.
    """

    junction_support: Dict[Tuple[str, str], Counter] = None
    read_ends: Dict[str, Tuple[int, int]] = None
    n_primary: int = 0
    n_supplementary: int = 0
    n_secondary: int = 0
    complete: bool = False
    error: Optional[str] = None


def collect_ranking_bam_evidence(bam_path: str) -> RankingBamEvidence:
    """ONE whole-BAM pass -> junction support + primary read ends + counts.

    Replaces the former separate ``collect_junction_support`` +
    ``collect_primary_read_ends`` calls for ranking/evidence, halving BAM
    scan cost. Junction counting matches ``junction_snap.collect_junction_support``
    (primary records only, exact intron coordinates); read ends match the
    runner's full-length map (first primary occurrence per read id wins).
    """
    import pysam

    from fin.candidates.intron_chains import extract_intron_chain

    out = RankingBamEvidence(
        junction_support=defaultdict(Counter), read_ends={},
    )
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for record in bam.fetch(until_eof=True):
                if record.is_unmapped:
                    continue
                if record.is_secondary:
                    out.n_secondary += 1
                    continue
                if record.is_supplementary:
                    out.n_supplementary += 1
                    continue
                out.n_primary += 1
                rid = record.query_name
                start, end = record.reference_start, record.reference_end
                if rid is not None and rid not in out.read_ends \
                        and start is not None and end is not None:
                    out.read_ends[rid] = (int(start), int(end))
                if record.cigartuples:
                    strand = "-" if record.is_reverse else "+"
                    chain = extract_intron_chain(
                        record.cigartuples, record.reference_start
                    )
                    for intron in chain.introns:
                        out.junction_support[
                            (record.reference_name, strand)
                        ][intron] += 1
        out.complete = True
    except Exception as exc:
        out.error = f"{type(exc).__name__}: {exc}"
        logger.warning(
            "ranking BAM evidence scan incomplete (%s); the ranking filter "
            "must not run on this partial evidence", out.error,
        )
    return out


def collect_primary_read_ends(bam_path: str) -> Dict[str, Tuple[int, int]]:
    """One full-BAM pass -> {read_id: (reference_start, reference_end)}.

    Primary mapped alignments only; first occurrence per read id wins
    (mirrors the runner's full-length end map). Retained for callers that
    only need spans; the ranking path uses the combined single-pass
    :func:`collect_ranking_bam_evidence`.
    """
    import pysam

    ends: Dict[str, Tuple[int, int]] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for record in bam.fetch(until_eof=True):
            if record.is_unmapped or record.is_secondary or record.is_supplementary:
                continue
            rid = record.query_name
            if rid is None or rid in ends:
                continue
            start, end = record.reference_start, record.reference_end
            if start is None or end is None:
                continue
            ends[rid] = (int(start), int(end))
    return ends


def _introns(qr) -> Tuple[Tuple[int, int], ...]:
    exons = sorted(qr.exons)
    return tuple(
        (exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)
    )


def _junction_canonical(
    intron: Tuple[int, int], strand: str, genome: str,
    motif_set: frozenset,
) -> Optional[bool]:
    s, e = intron
    n = len(genome)
    if s < 0 or e > n or e - 2 < 0 or s + 2 > n:
        return None
    up = genome[s:s + 2].upper()
    dn = genome[e - 2:e].upper()
    donor, acceptor = (up, dn) if strand == "+" else (_revcomp2(dn), _revcomp2(up))
    return (donor, acceptor) in motif_set


def _is_strict_subchain(a: Sequence, b: Sequence) -> bool:
    if len(a) == 0 or len(a) >= len(b):
        return False
    return any(tuple(b[i:i + len(a)]) == tuple(a) for i in range(len(b) - len(a) + 1))


def compute_candidate_evidence(
    results: Mapping[str, object],
    *,
    junction_support: Optional[Mapping[Tuple[str, str], Mapping[Tuple[int, int], int]]] = None,
    read_ends: Optional[Mapping[str, Tuple[int, int]]] = None,
    genome: Optional[Mapping[str, str]] = None,
    canonical_motifs: frozenset = frozenset(),
    end_window_bp: int = 25,
) -> list[dict]:
    """Compute the observable feature row for every candidate in ``results``.

    Sentinel -1.0 marks "evidence source unavailable" (no BAM junction map,
    no genome, no read spans); a consumer must treat -1 as missing, never as
    a low value.
    """
    chains = {cid: _introns(qr) for cid, qr in results.items()}

    # Family grouping: persisted discovery family when present, else the
    # candidate stands alone (share 1.0 within itself).
    families: Dict[str, list] = defaultdict(list)
    for cid, qr in results.items():
        fam = getattr(qr, "family_id", None) or f"__solo__{cid}"
        families[fam].append(cid)

    # Locus geometry: sibling sub/superchain among survivors on same key.
    by_key: Dict[Tuple[str, str], list] = defaultdict(list)
    for cid, qr in results.items():
        by_key[(qr.chrom, qr.strand)].append(cid)

    rows = []
    for cid, qr in results.items():
        introns = chains[cid]
        n_exons = len(qr.exons)
        length = sum(e - s for s, e in qr.exons)
        abundance = float(qr.abundance)
        n_reads = int(getattr(qr, "num_assigned_reads", 0))
        ratio = abundance / n_reads if n_reads > 0 else -1.0

        fam = getattr(qr, "family_id", None) or f"__solo__{cid}"
        members = families[fam]
        fam_ab = [(float(results[m].abundance), m) for m in members]
        fam_total = sum(a for a, _ in fam_ab)
        share = abundance / fam_total if fam_total > 0 else -1.0
        # rank: 1 = dominant. Deterministic tie-break by candidate id.
        rank = 1 + sum(
            1 for a, m in fam_ab
            if m != cid and (a > abundance or (a == abundance and m < cid))
        )
        dominant = max(a for a, _ in fam_ab) / fam_total if fam_total > 0 else -1.0

        sub = sup = False
        if introns:
            for other in by_key[(qr.chrom, qr.strand)]:
                if other == cid:
                    continue
                oc = chains[other]
                if not oc:
                    continue
                if _is_strict_subchain(introns, oc):
                    sub = True
                if _is_strict_subchain(oc, introns):
                    sup = True
                if sub and sup:
                    break

        if introns and junction_support is not None:
            counter = junction_support.get((qr.chrom, qr.strand), {})
            support = sorted(int(counter.get(i, 0)) for i in introns)
            weakest = float(support[0])
            median = float(support[len(support) // 2])
            below3 = sum(1 for s in support if s < 3)
        else:
            weakest = median = -1.0
            below3 = -1

        canonical_fraction = -1.0
        if introns and genome and qr.chrom in genome and canonical_motifs:
            verdicts = [
                _junction_canonical(i, qr.strand, genome[qr.chrom], canonical_motifs)
                for i in introns
            ]
            known = [v for v in verdicts if v is not None]
            if known:
                canonical_fraction = sum(known) / len(known)

        end5 = end3 = full = -1.0
        n_span = 0
        assigned = tuple(getattr(qr, "assigned_read_ids", ()) or ())
        if read_ends is not None and assigned:
            five = qr.start if qr.strand == "+" else qr.end
            three = qr.end if qr.strand == "+" else qr.start
            hits5 = hits3 = hitsf = 0
            for rid in assigned:
                span = read_ends.get(rid)
                if span is None:
                    continue
                n_span += 1
                r5 = span[0] if qr.strand == "+" else span[1]
                r3 = span[1] if qr.strand == "+" else span[0]
                ok5 = abs(r5 - five) <= end_window_bp
                ok3 = abs(r3 - three) <= end_window_bp
                hits5 += ok5
                hits3 += ok3
                hitsf += ok5 and ok3
            if n_span > 0:
                end5 = hits5 / n_span
                end3 = hits3 / n_span
                full = hitsf / n_span

        rows.append({
            "candidate_id": cid,
            "chrom": qr.chrom,
            "strand": qr.strand,
            "start": qr.start,
            "end": qr.end,
            "source": qr.source,
            "n_exons": n_exons,
            "tx_length": length,
            "is_mono": int(not introns),
            "abundance": round(abundance, 4),
            "num_reads": n_reads,
            "soft_hard_ratio": round(ratio, 4),
            "confidence": round(float(getattr(qr, "confidence", 0.0)), 4),
            "max_R": round(float(getattr(qr, "max_R", 0.0)), 4),
            "family_id": fam,
            "family_size": len(members),
            "family_share": round(share, 6),
            "family_rank": rank,
            "family_dominant_share": round(dominant, 6),
            "is_subchain_of_sibling": int(sub),
            "is_superchain_of_sibling": int(sup),
            "weakest_junction_support": weakest,
            "median_junction_support": median,
            "n_junctions_below3": below3,
            "canonical_fraction": round(canonical_fraction, 4),
            "end5_support_frac": round(end5, 4),
            "end3_support_frac": round(end3, 4),
            "fulllen_frac": round(full, 4),
            "n_end_reads": n_span,
        })
    rows.sort(key=lambda r: (r["chrom"], r["start"], r["candidate_id"]))
    return rows


def write_candidate_evidence(path: str | Path, rows: Iterable[dict]) -> None:
    """Atomically write the evidence table."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.with_suffix(output.suffix + ".tmp")
    with tmp.open("w") as handle:
        handle.write("\t".join(EVIDENCE_COLUMNS) + "\n")
        for row in rows:
            handle.write(
                "\t".join(str(row[c]) for c in EVIDENCE_COLUMNS) + "\n"
            )
    tmp.replace(output)
    logger.info("Wrote candidate evidence: %s", output)
