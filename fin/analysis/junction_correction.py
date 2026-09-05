"""Per-read junction correction (step 5a), applied BEFORE clustering.

The problem this solves: within one cluster every read may carry the same
mis-assembled junction, so no candidate in that cluster has the correct
coordinates -- but another cluster in the SAME interval carries the correct
junction. Comparing only inside the cluster can never recover it.

The mechanism, per interval:

1. REGISTRY -- every distinct junction observed in the interval's primary
   MAPQ>=1 reads, with its raw read count, its indel-clean count (no indel
   within 10 bp of either boundary) and the ids of the supporting reads.
2. NOMINATION -- for a junction J of a read, EVERY other registry junction
   on the same strand with both coordinates within ``tol_bp`` is an
   alternative. There is no support threshold: the registry is a plain
   record and the signal decides.
3. BALLOT -- the read's own chain and, for each alternative J', the chain
   with J replaced by J' are stitched into full transcript sequences from
   the genome. All K hypotheses are eventaligned against the read's raw
   signal in ONE batched krill call with the SAME fixed initial scaling
   (krill >= 0.1.6 ``init_scale``/``init_shift``, taken from a first
   full-read run), so they are normalised identically. The scored window is
   the differing region plus at least ``flank_bp`` (default 4) of SHARED
   sequence on each side; a hypothesis that has no event in both shared
   flanks is unscorable and the read abstains. Over the union raw-sample
   window we sum the per-event Gaussian NLL under each hypothesis; the
   lowest wins, and the read is re-chained only when it beats the read's
   own junction by more than ``margin``.
4. PROPOSE -- correction is ADDITIVE. A read's raw chain and its discovery
   membership are never changed. A read yields at most ONE proposal: the
   single-junction substitution with the largest confidence, accepted only
   when the alternative beats the read's own junction by >= ``margin`` AND
   beats the runner-up alternative by >= ``runner_up_margin`` (so a K-way
   contest where two alternatives are indistinguishable abstains).
5. INSTANTIATE -- after the ordinary discovery (grouping, union-find,
   collapse, emission) has produced the baseline candidates from the RAW
   chains, proposals are aggregated by corrected full chain. A chain that
   already exists as a candidate only gains correction support; a chain that
   does not exist becomes a NEW ``source="novel"`` candidate when at least
   ``new_chain_min_reads`` distinct reads proposed exactly it. Baseline
   candidates are never removed or re-collapsed, so
   ``candidates_after ⊇ candidates_before``.

A corrected candidate then competes in assignment/EM and every selection gate
like any other; correction grants no exemption anywhere, and correction votes
never count as raw complete-chain attestation.

Opt-in (``junction_correction_mode``); ``off`` does nothing, ``audit``
computes and records everything but appends no candidate.
"""
from __future__ import annotations

import logging
import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)

Junction = Tuple[int, int]
Chain = Tuple[Junction, ...]

JC_FLANK_BP = 10        # indel-clean flank, matches discovery's junction_clean_depths
JC_SHARED_FLANK_BP = 4  # shared sequence required on EACH side of the differing region
JC_MIN_EVENTS = 3       # per hypothesis, inside the window
NNNNN_NLL = 6.0         # charge for an unmodelable event (matches m2_junction_nll)

TSV_FIELDS = (
    "interval", "read_id", "strand", "own_junction", "alt_junction", "shift_donor",
    "shift_acceptor", "own_raw", "own_clean", "alt_raw", "alt_clean", "own_nll",
    "alt_nll", "margin", "runner_up_margin", "n_events_own", "n_events_alt",
    "n_samples", "decision", "reason",
)


@dataclass
class JunctionRecord:
    raw: int = 0
    clean: int = 0
    read_ids: set = field(default_factory=set)


@dataclass
class CorrectionRow:
    interval: str
    read_id: str
    strand: str
    own_junction: str
    alt_junction: str
    shift_donor: int
    shift_acceptor: int
    own_raw: int
    own_clean: int
    alt_raw: int
    alt_clean: int
    own_nll: float = float("nan")
    alt_nll: float = float("nan")
    margin: float = float("nan")
    runner_up_margin: float = float("nan")
    n_events_own: int = 0
    n_events_alt: int = 0
    n_samples: int = 0
    decision: str = "abstain"
    reason: str = ""

    def as_row(self) -> List[str]:
        def f(x):
            return "" if x is None or (isinstance(x, float) and not np.isfinite(x)) else (
                f"{x:.4f}" if isinstance(x, float) else str(x))
        return [f(getattr(self, k)) for k in TSV_FIELDS]


@dataclass
class CorrectionProposal:
    """One read's single best single-junction substitution."""
    read_id: str
    original_chain: Chain
    corrected_chain: Chain
    changed_junction: Junction
    replacement: Junction
    own_margin: float
    runner_up_margin: float
    confidence: float
    reference_start: int
    reference_end: int


@dataclass
class JunctionCorrectionResult:
    proposals: List[CorrectionProposal] = field(default_factory=list)
    rows: List[CorrectionRow] = field(default_factory=list)


# --------------------------------------------------------------------------- registry
def build_registry(read_chains, indels_by_read) -> Dict[Junction, JunctionRecord]:
    """Interval-wide junction registry from the RAW read chains.

    ``read_chains`` is discovery's ``[(read_dict, IntronChain)]``; a read counts
    toward a junction's clean support when it has no indel within
    ``JC_FLANK_BP`` of either boundary. Reads absent from ``indels_by_read``
    (MAPQ-filtered at collection) contribute raw support only.
    """
    reg: Dict[Junction, JunctionRecord] = {}
    for rd, chain in read_chains:
        q = rd.get("query_name")
        if q is None:
            continue
        ind = indels_by_read.get(q)
        for j in tuple(chain.introns):
            rec = reg.setdefault(j, JunctionRecord())
            rec.raw += 1
            rec.read_ids.add(q)
            if ind is not None and all(
                not ((lo - JC_FLANK_BP <= j[0] <= hi + JC_FLANK_BP)
                     or (lo - JC_FLANK_BP <= j[1] <= hi + JC_FLANK_BP))
                for lo, hi in ind
            ):
                rec.clean += 1
    return reg


def nominate(junction: Junction, registry: Dict[Junction, JunctionRecord], *, tol_bp: int) -> List[Junction]:
    """Every OTHER registry junction with both coordinates within ``tol_bp``.

    No support threshold: the registry only records; the signal ballot decides.
    Strand is implicit (the registry is built per strand-separated interval).
    Deterministically sorted.
    """
    return sorted(
        k for k in registry
        if k != junction
        and abs(k[0] - junction[0]) <= tol_bp and abs(k[1] - junction[1]) <= tol_bp
    )


# --------------------------------------------------------------------------- hypotheses
def substitute(chain: Chain, old: Junction, new: Junction) -> Chain:
    return tuple(new if j == old else j for j in chain)


def chain_is_valid(start: int, end: int, chain: Chain) -> bool:
    """Ordered, non-overlapping introns strictly inside the read's span."""
    prev = start
    for d, a in chain:
        if not (prev < d < a < end):
            return False
        prev = a
    return True


def _exon_blocks(start: int, end: int, chain: Chain):
    blocks, cur = [], start
    for d, a in chain:
        blocks.append((cur, d))
        cur = a
    blocks.append((cur, end))
    return blocks


_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def build_transcript(genome_seq: str, start: int, end: int, chain: Chain, strand: str):
    """(spliced sequence, transcript index -> genomic position), transcript 5'->3'."""
    blocks = _exon_blocks(start, end, chain)
    seq = "".join(genome_seq[s:e] for s, e in blocks)
    gmap = np.concatenate([np.arange(s, e) for s, e in blocks]).astype(np.int64)
    if strand == "-":
        seq = seq.translate(_COMP)[::-1]
        gmap = gmap[::-1].copy()
    return seq, gmap


# --------------------------------------------------------------------------- ballot
def _events_in_genomic(res, gmap, g_lo: int, g_hi: int):
    """Indices of events whose k-mer centre maps into genomic [g_lo, g_hi]."""
    if res is None or res.get("status", -1) != 0 or len(res["position"]) == 0:
        return None
    pos = np.asarray(res["position"], np.int64)
    ok = (pos >= 0) & (pos < len(gmap))
    gp = np.where(ok, gmap[np.clip(pos, 0, len(gmap) - 1)], -1)
    return np.nonzero(ok & (gp >= g_lo) & (gp <= g_hi))[0]


def _sample_range(res, idx) -> Optional[Tuple[int, int]]:
    if idx is None or len(idx) == 0:
        return None
    s0 = np.asarray(res["start_idx"], np.int64)[idx]
    s1 = np.asarray(res["end_idx"], np.int64)[idx]
    return int(np.minimum(s0, s1).min()), int(np.maximum(s0, s1).max())


def _sum_nll_samples(res, lo_s: int, hi_s: int) -> Optional[Tuple[float, int]]:
    """Sum per-event Gaussian NLL over every event overlapping raw samples [lo_s, hi_s).

    sigma = model_stdv * sqrt(var) reproduces krill's own ``standardized_level``.
    An unmodelable event is charged NNNNN_NLL so a hypothesis cannot win by
    declaring the contested bases unexplainable.
    """
    if res is None or res.get("status", -1) != 0 or len(res["position"]) == 0:
        return None
    s0 = np.asarray(res["start_idx"], np.int64)
    s1 = np.asarray(res["end_idx"], np.int64)
    lo = np.minimum(s0, s1); hi = np.maximum(s0, s1)
    sel = (hi > lo_s) & (lo < hi_s)
    if not sel.any():
        return None
    mm = np.asarray(res["model_mean"], float)[sel]
    ms = np.asarray(res["model_stdv"], float)[sel]
    lev = np.asarray(res["event_level_mean"], float)[sel]
    var = float(res.get("var", 1.0)) or 1.0
    sig = ms * math.sqrt(max(var, 1e-12))
    good = np.isfinite(mm) & np.isfinite(lev) & np.isfinite(sig) & (sig > 0)
    nll = np.full(int(sel.sum()), NNNNN_NLL, float)
    if good.any():
        z = (lev[good] - mm[good]) / sig[good]
        v = 0.5 * z * z + np.log(sig[good])
        v[~np.isfinite(v)] = NNNNN_NLL
        nll[good] = v
    return float(nll.sum()), int(sel.sum())


def _diff_span(own_j: Junction, alts: Sequence[Junction], strand: str, flank: int):
    """Genomic scored span for one contest and the two SHARED exonic flanks.

    All hypotheses keep the upstream exon up to at least ``min(donors)`` and the
    downstream exon from ``max(acceptors)`` on, so those are the shared contexts;
    everything between ``min(donors)`` and ``max(acceptors)`` is where they differ
    (exonic in some, intronic in others -- intronic positions carry no events).
    Returns ``(g_lo, g_hi, up_flank, dn_flank)`` or None when nothing differs.
    """
    donors = {own_j[0]} | {a[0] for a in alts}
    accs = {own_j[1]} | {a[1] for a in alts}
    if len(donors) == 1 and len(accs) == 1:
        return None
    up_flank = (min(donors) - flank, min(donors) - 1)       # last shared exonic bases
    dn_flank = (max(accs), max(accs) + flank - 1)           # first shared exonic bases
    return up_flank[0], dn_flank[1], up_flank, dn_flank


def score_read_kway(
    aligner,
    read_id: str,
    signal_rec: dict,
    genome_seq: str,
    start: int,
    end: int,
    strand: str,
    own_chain: Chain,
    own_j: Junction,
    alts: Sequence[Junction],
    *,
    flank_bp: int = JC_SHARED_FLANK_BP,
):
    """Score one read against its own chain and K alternative chains (one junction
    swapped) in ONE batched krill call with fixed initial scaling.

    Returns ``(nlls, n_events, n_samples, reason)`` where ``nlls[0]`` is the read's
    own hypothesis and ``nlls[1:]`` follow ``alts`` (None where unscorable).
    """
    hyps = [own_chain]
    kept_alts = []
    for a in alts:
        ch = substitute(own_chain, own_j, a)
        if ch != own_chain and chain_is_valid(start, end, ch) and len(set(ch)) == len(ch):
            hyps.append(ch); kept_alts.append(a)
    if not kept_alts:
        return None, [], 0, "no_valid_alternative", []
    seqs, gmaps = [], []
    for ch in hyps:
        tx, gm = build_transcript(genome_seq, start, end, ch, strand)
        seqs.append(tx); gmaps.append(gm)
    if any(not t for t in seqs):
        return None, [], 0, "empty_sequence", kept_alts
    common = dict(signal=signal_rec["signal"], digitisation=signal_rec["digitisation"],
                  offset=signal_rec["offset"], range=signal_rec["range"],
                  sample_rate=signal_rec["sample_rate"])
    try:
        full = aligner.align([dict(read_id=read_id, sequence=seqs[0], start=0, **common)])[0]
        if full.get("status", -1) != 0:
            return None, [], 0, "full_align_failed", kept_alts
        isc, ish = float(full["init_scale"]), float(full["init_shift"])
        results = aligner.align([
            dict(read_id=read_id, sequence=t, start=0, init_scale=isc, init_shift=ish, **common)
            for t in seqs])
    except Exception as exc:  # noqa: BLE001 - a signal failure must abstain, never crash
        logger.debug("junction correction: krill failed for %s: %s", read_id, exc)
        return None, [], 0, "krill_error", kept_alts
    span = _diff_span(own_j, kept_alts, strand, flank_bp)
    if span is None:
        return None, [], 0, "no_difference", kept_alts
    g_lo, g_hi, lo_flank, hi_flank = span
    # union raw-sample window over every hypothesis' in-span events; every hypothesis
    # must have an event in BOTH shared flanks (otherwise the contrast is one-sided)
    lo_s = hi_s = None
    for res, gm in zip(results, gmaps):
        idx = _events_in_genomic(res, gm, g_lo, g_hi)
        if idx is None or len(idx) == 0:
            return None, [], 0, "no_events_in_window", kept_alts
        if (len(_events_in_genomic(res, gm, *lo_flank)) == 0
                or len(_events_in_genomic(res, gm, *hi_flank)) == 0):
            return None, [], 0, "flank_not_covered", kept_alts
        r = _sample_range(res, idx)
        lo_s = r[0] if lo_s is None else min(lo_s, r[0])
        hi_s = r[1] if hi_s is None else max(hi_s, r[1])
    if hi_s <= lo_s:
        return None, [], 0, "no_events_in_window", kept_alts
    nlls, nevs = [], []
    for res in results:
        sc = _sum_nll_samples(res, lo_s, hi_s)
        if sc is None or sc[1] < JC_MIN_EVENTS:
            return None, [], hi_s - lo_s, "too_few_events", kept_alts
        nlls.append(sc[0]); nevs.append(sc[1])
    return nlls, nevs, hi_s - lo_s, "", kept_alts


# --------------------------------------------------------------------------- driver
def require_krill_init_scaling(aligner) -> None:
    """Fail loudly when the installed krill cannot echo the initial scaling.

    Both hypotheses MUST be eventaligned with the same fixed init_scale/init_shift
    (krill >= 0.1.6). Without it each hypothesis re-estimates its own scaling from
    its own events and the two NLLs are not comparable -- a silent degradation the
    frozen validation protocol forbids.
    """
    import numpy as _np

    probe = aligner.align([dict(
        read_id="__krill_capability_probe__", sequence="ACGTTGCA" * 80, start=0,
        signal=_np.zeros(4000, dtype="float32"), digitisation=8192.0, offset=0.0,
        range=1500.0, sample_rate=4000.0)])
    if not probe or "init_scale" not in probe[0]:
        raise RuntimeError(
            "junction_correction_mode requires krill >= 0.1.6: the installed build "
            "does not report init_scale/init_shift, so the two junction hypotheses "
            "cannot be scored under one fixed initial scaling. Install the 0.1.6 "
            "wheel (see docs/install) or run with --junction-correction-mode off."
        )


def correct_read_chains(
    read_chains,
    indels_by_read,
    *,
    config,
    interval,
    genome_seq: str,
    signal_reader,
    aligner,
    krill_mod,
) -> JunctionCorrectionResult:
    """Step 5a ballots. Returns proposals + audit rows; never mutates ``read_chains``.

    Every read carrying a junction that has a registry neighbour within tolerance
    is scored (no cap, no support threshold); reads are processed in sorted id
    order over a frozen registry, so no decision can depend on another.
    """
    mode = getattr(config, "junction_correction_mode", "off")
    tol = int(getattr(config, "junction_correction_tol_bp", 15))
    margin = float(getattr(config, "junction_correction_margin", 10.0))
    ru_margin = float(getattr(config, "junction_correction_runner_up_margin", 10.0))
    flank = int(getattr(config, "junction_correction_flank_bp", JC_SHARED_FLANK_BP))
    result = JunctionCorrectionResult()
    rows = result.rows
    if mode == "off" or not read_chains:
        return result

    registry = build_registry(read_chains, indels_by_read)
    if not registry:
        return result
    noms = {j: nominate(j, registry, tol_bp=tol) for j in sorted(registry)}
    noms = {j: v for j, v in noms.items() if v}
    if not noms:
        return result
    require_krill_init_scaling(aligner)   # only when there is real work to do

    by_id = {}
    for rd, chain in read_chains:
        q = rd.get("query_name")
        if q is not None:
            by_id.setdefault(q, (rd, chain))
    todo = sorted(q for q, (rd, ch) in by_id.items()
                  if any(j in noms for j in tuple(ch.introns)))

    iname = interval.region_string
    for q in todo:
        rd, chain = by_id[q]
        rec = signal_reader(q)
        own_chain = tuple(chain.introns)
        strand = interval.strand or ("-" if rd.get("is_reverse") else "+")
        if rec is None:
            for j in own_chain:
                for alt in noms.get(j, ()):
                    rows.append(_row(iname, q, rd, registry, j, alt, reason="no_signal"))
            continue
        best_for_read = None          # (confidence, j, alt, own_m, ru_m, row)
        for j in own_chain:
            alts = noms.get(j)
            if not alts:
                continue
            nlls, nevs, nsamp, reason, kept = score_read_kway(
                aligner, q, rec, genome_seq, rd.get("reference_start"),
                rd.get("reference_end"), strand, own_chain, j, alts, flank_bp=flank)
            if nlls is None:
                for alt in alts:
                    rows.append(_row(iname, q, rd, registry, j, alt, reason=reason))
                continue
            own_nll = nlls[0]
            alt_nlls = sorted(((nlls[k + 1], alt) for k, alt in enumerate(kept)))
            best_nll, best_alt = alt_nlls[0]
            second = alt_nlls[1][0] if len(alt_nlls) > 1 else float("inf")
            own_m = own_nll - best_nll
            ru_m = second - best_nll
            for k, alt in enumerate(kept):
                r = _row(iname, q, rd, registry, j, alt)
                r.own_nll, r.alt_nll = own_nll, nlls[k + 1]
                r.margin = own_nll - nlls[k + 1]
                r.runner_up_margin = ru_m if alt == best_alt else float("nan")
                r.n_events_own, r.n_events_alt, r.n_samples = nevs[0], nevs[k + 1], nsamp
                if alt == best_alt:
                    if own_m < margin:
                        r.reason = "own_margin_below_threshold"
                    elif ru_m < ru_margin:
                        r.reason = "runner_up_margin_below_threshold"
                    else:
                        conf = min(own_m, ru_m)
                        if best_for_read is None or conf > best_for_read[0]:
                            if best_for_read is not None:
                                best_for_read[5].decision = "abstain"
                                best_for_read[5].reason = "not_best_edit_for_read"
                            best_for_read = (conf, j, alt, own_m, ru_m, r)
                        else:
                            r.reason = "not_best_edit_for_read"
                rows.append(r)
        if best_for_read is None:
            continue
        conf, j, alt, own_m, ru_m, r = best_for_read
        new_chain = substitute(own_chain, j, alt)
        if not chain_is_valid(rd.get("reference_start"), rd.get("reference_end"), new_chain):
            r.reason = "invalid_chain"
            continue
        r.decision = "propose"
        result.proposals.append(CorrectionProposal(
            read_id=q, original_chain=own_chain, corrected_chain=new_chain,
            changed_junction=j, replacement=alt, own_margin=own_m, runner_up_margin=ru_m,
            confidence=conf, reference_start=int(rd.get("reference_start")),
            reference_end=int(rd.get("reference_end"))))
    if result.proposals:
        logger.info("junction correction %s: %d proposal(s) from %d ballots",
                    iname, len(result.proposals), len(rows))
    return result


def aggregate_proposals(proposals: Sequence[CorrectionProposal], *, min_reads: int):
    """Group proposals by exact corrected chain -> {chain: sorted voter ids};
    only chains with >= ``min_reads`` DISTINCT voters are returned."""
    by_chain: Dict[Chain, set] = {}
    for p in proposals:
        by_chain.setdefault(p.corrected_chain, set()).add(p.read_id)
    return {ch: sorted(v) for ch, v in sorted(by_chain.items()) if len(v) >= min_reads}


def _row(iname, q, rd, registry, j, alt, reason="") -> CorrectionRow:
    ro = registry.get(j, JunctionRecord()); ra = registry.get(alt, JunctionRecord())
    return CorrectionRow(
        interval=iname, read_id=q, strand=("-" if rd.get("is_reverse") else "+"),
        own_junction=f"{j[0]}-{j[1]}", alt_junction=f"{alt[0]}-{alt[1]}",
        shift_donor=alt[0] - j[0], shift_acceptor=alt[1] - j[1],
        own_raw=ro.raw, own_clean=ro.clean, alt_raw=ra.raw, alt_clean=ra.clean,
        reason=reason)


# --------------------------------------------------------------------------- artifact
def write_correction_tsv(path, rows: Sequence[CorrectionRow]) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    body = sorted(r.as_row() for r in rows)
    tmp = p.with_suffix(p.suffix + ".tmp")
    with tmp.open("w") as fh:
        fh.write("\t".join(TSV_FIELDS) + "\n")
        for r in body:
            fh.write("\t".join(r) + "\n")
    tmp.replace(p)


def merge_correction_shards(work_dir) -> Optional[str]:
    """Merge per-interval worker shards into ``junction_correction.tsv``."""
    wd = Path(work_dir)
    shard_dir = wd / "_jcorrection_shards"
    out = wd / "junction_correction.tsv"
    body: List[str] = []
    if shard_dir.is_dir():
        for shard in sorted(shard_dir.glob("*.tsv")):
            with shard.open() as fh:
                next(fh, None)
                body.extend(line.rstrip("\n") for line in fh if line.strip())
    body.sort()
    tmp = out.with_suffix(".tsv.tmp")
    with tmp.open("w") as fh:
        fh.write("\t".join(TSV_FIELDS) + "\n")
        for line in body:
            fh.write(line + "\n")
    tmp.replace(out)
    if shard_dir.is_dir():
        shutil.rmtree(shard_dir)
    return str(out)
