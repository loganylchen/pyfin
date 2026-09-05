"""Signal-ballot junction election: which of several nearby junctions is real?

When alignment noise (or corrupt annotation) proposes several junction
variants within a few bp of each other, base-space evidence is inconclusive
by construction: the misaligned reads themselves "support" the spurious
variant. The dRNA raw signal is the orthogonal referee. This module holds a
K-WAY election per contested neighbourhood:

1. NOMINATION (cheap, cigar-level screen): a junction whose indel-clean
   attestation is <= JE_SCREEN_CLEAN_MAX while at least one competitor
   within +-JE_WINDOW_BP holds >= max(JE_SCREEN_COMPETITOR_MIN,
   JE_SCREEN_RATIO x) clean attestations is ACCUSED. Every qualifying
   competitor within the window joins the same election; accused variants
   sharing a competitor merge into one contest, so a five-variant
   neighbourhood is decided by one five-way vote, not five pairwise duels.
   Both sides must be carried by emitted candidates; GTF junctions
   participate exactly like read-derived ones - the signal decides, the
   annotation is only a nominee.
2. WINDOW PURITY: ballots are scored on the merged +-JE_FLANK_BP genomic
   windows around every nominated coordinate. Each nominee's representative
   (its highest-abundance carrier) must present the SAME exonic structure
   inside those windows apart from its own variant, otherwise the contest
   abstains (``window_impure``) - unrelated structural differences must not
   leak into the vote.
3. BALLOTS: every assigned read spanning the contested region votes with
   its raw signal - the krill eventalign per-event mean NLL
   (``0.5 z^2 + log sigma``) against each nominee's representative over the
   shared windows (the validated M2 wobble metric; SIRV wobble pairs,
   pooled per-read accuracy 86.2%). A ballot is valid only when at least
   two hypotheses are scored with >= JE_BALLOT_MIN_EVENTS events each and
   the best-vs-runner-up margin reaches ``m2_tiebreak_margin`` (the same
   validated production threshold M2 uses); anything weaker abstains.
4. VERDICT (tri-state, per accused variant): with >= JE_QUORUM valid
   ballots, a winner holding a >= JE_MAJORITY supermajority DEFEATS every
   accused variant whose own share is <= 1 - JE_MAJORITY; ``require`` mode
   then drops the carriers of defeated variants. ``upheld`` (the accused is
   itself the supermajority winner) and every insufficient outcome
   (``no_quorum`` / ``split`` / ``no_ballots`` / ``window_impure``) retain -
   abstention never removes a model. ``audit`` only writes the artifact and
   is byte-inert.

Ballot ablation modes (``junction_election_ballot``, development-only, for
proving that the raw signal - not the nomination screen - carries the
decision): ``signal`` (production), ``m1`` (mappy alignment-score vote),
``screen`` (defeat on nomination alone, no ballots), ``shuffle``
(signal scores computed, then deterministically permuted per read -
negative control).

Burden of proof (``junction_election_burden``): ``prosecution`` (default) -
the accused is retained unless the vote convicts it (abstention retains);
``defense`` - a NOMINATED junction already carries the misalignment
signature (indel-unclean, heavily outvoted in clean support), so it is
dropped unless the signal vote EXONERATES it (``upheld`` with quorum and
supermajority); every other outcome, including abstention, follows the
screen and drops. Measured on the grid: ballots are 97.6% accurate when
they decide but abstain 45% of the time (contested junk is read-poor), so
prosecution is the conservative mode and defense the aggressive one.

Leakage rules: every input is observable at inference time (read chains,
cigar indels, EM abundances, raw signal). Truth GTFs, gffcompare classes,
and NanoCount never enter. All thresholds are development-fit on the six
SIRV H9 samples and require independent validation before promotion.
"""

from __future__ import annotations

import collections
import hashlib
import logging
from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

# Development-fit election constants (six H9 SIRV samples; module docstring).
JE_WINDOW_BP = 15            # competitor search window per junction coordinate
JE_FLANK_BP = 8              # signal scoring flank around nominated coordinates
JE_SCREEN_CLEAN_MAX = 1      # accused junction: indel-clean support <= this
JE_SCREEN_COMPETITOR_MIN = 3  # competitor junction: clean support >= this
JE_SCREEN_RATIO = 3          # ... and >= ratio x accused clean support
JE_QUORUM = 5                # minimum valid ballots to call an election
JE_MAJORITY = 0.75           # supermajority fraction required to elect
JE_BALLOT_MIN_EVENTS = 3     # min eventalign events per scored hypothesis
JE_MAX_ELECTORATE = 30       # deterministic cap on ballots per contest
JE_MAX_CONTESTS = 50         # deterministic cap on contests per interval

VERDICT_DEFEATED = "defeated"
VERDICT_UPHELD = "upheld"


@dataclass
class ElectionRow:
    """One accused variant's outcome within a K-way contest (artifact row)."""

    interval: str
    contest_id: str              # deterministic id of the K-way contest
    accused_junction: str        # "donor-acceptor" of this accused variant
    nominees: str                # all K junctions, comma-joined
    accused_clean: int
    best_competitor_clean: int
    accused_ids: str             # candidate ids carrying the accused variant
    winner_junction: str = ""
    n_electorate: int = 0
    n_ballots: int = 0
    votes_accused: int = 0
    votes_winner: int = 0
    verdict: str = "no_ballots"  # defeated|upheld|split|no_quorum|no_ballots|window_impure
    action: str = "none"         # none | drop

    def as_tsv_row(self) -> List[str]:
        return [str(getattr(self, f)) for f in ELECTION_TSV_FIELDS]


ELECTION_TSV_FIELDS = [
    "interval", "contest_id", "accused_junction", "nominees",
    "accused_clean", "best_competitor_clean", "accused_ids",
    "winner_junction", "n_electorate", "n_ballots", "votes_accused",
    "votes_winner", "verdict", "action",
]


def _chain(cand) -> Tuple[Tuple[int, int], ...]:
    ic = getattr(cand, "intron_chain", None)
    return tuple(getattr(ic, "introns", ()) or ())


def _jstr(j) -> str:
    return f"{j[0]}-{j[1]}"


def _near(a, b) -> bool:
    return abs(a[0] - b[0]) <= JE_WINDOW_BP and abs(a[1] - b[1]) <= JE_WINDOW_BP


def nominate_contests(
    cand_list: Sequence,
    quant_by_id: Dict[str, "object"],
    junction_support: Dict[Tuple[int, int], Tuple[int, int]],
) -> Tuple[List[dict], int]:
    """Cigar-level screen producing K-way contests.

    Returns ``(contests, n_truncated)``. Each contest dict:
    ``{nominees: [junction,...], reps: {junction: candidate},
    carriers: {junction: [candidates]}, accused: {junction: [candidates]},
    clean: {junction: int}}``.
    Deterministic; capped at JE_MAX_CONTESTS with loud truncation reporting.
    """
    carriers: Dict[Tuple[int, int], List] = {}
    for c in cand_list:
        if getattr(c, "source", "") == "fusion":
            continue
        for j in _chain(c):
            carriers.setdefault(j, []).append(c)
    if not carriers or not junction_support:
        return [], 0

    def clean(j) -> int:
        return int(junction_support.get(j, (0, 0))[1])

    def abundance(cid) -> float:
        qr = quant_by_id.get(cid)
        return float(getattr(qr, "abundance", 0.0) or 0.0) if qr is not None else 0.0

    def rep(j):
        return max(
            carriers[j],
            key=lambda c: (abundance(c.candidate_id), c.candidate_id),
        )

    accused_js = []
    champ_by_accused: Dict[Tuple[int, int], List[Tuple[int, int]]] = {}
    for j in sorted(carriers):
        cj = clean(j)
        if cj > JE_SCREEN_CLEAN_MAX:
            continue
        bar = max(JE_SCREEN_COMPETITOR_MIN, JE_SCREEN_RATIO * max(cj, 1))
        champs = [
            k for k in sorted(carriers)
            if k != j and _near(k, j) and clean(k) >= bar
        ]
        if champs:
            accused_js.append(j)
            champ_by_accused[j] = champs

    if not accused_js:
        return [], 0

    # Merge accused variants that share any champion into one K-way contest.
    parent: Dict[Tuple[int, int], Tuple[int, int]] = {}

    def find(x):
        while parent.get(x, x) != x:
            parent[x] = parent.get(parent[x], parent[x])
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    by_champ: Dict[Tuple[int, int], List] = {}
    for j in accused_js:
        parent.setdefault(j, j)
        for k in champ_by_accused[j]:
            by_champ.setdefault(k, []).append(j)
    for k, js in by_champ.items():
        for other in js[1:]:
            union(js[0], other)

    groups: Dict[Tuple[int, int], List] = {}
    for j in accused_js:
        groups.setdefault(find(j), []).append(j)

    contests = []
    truncated = 0
    for root in sorted(groups):
        accused_group = sorted(groups[root])
        champs = sorted({
            k for j in accused_group for k in champ_by_accused[j]
        })
        nominees = sorted(set(accused_group) | set(champs))
        accused_map = {}
        for j in accused_group:
            champ_set = set(champs)
            acc = [
                c for c in sorted(
                    carriers[j], key=lambda c: c.candidate_id)
                if not (set(_chain(c)) & champ_set)
            ]
            if acc:
                accused_map[j] = acc
        if not accused_map:
            continue
        if len(contests) >= JE_MAX_CONTESTS:
            truncated += 1
            continue
        contests.append(dict(
            nominees=nominees,
            reps={j: rep(j) for j in nominees},
            # EVERY candidate carrying each nominated junction, not just the
            # representative: the electorate is drawn from all of them (a
            # junction's support is spread over its carriers, and taking only
            # the representative's reads under-counts a well-supported junction
            # into `no_quorum` -- which the defense burden then treats as a
            # conviction).
            carriers={j: list(carriers[j]) for j in nominees},
            accused=accused_map,
            clean={j: clean(j) for j in nominees},
        ))
    return contests, truncated


def _contest_windows(nominees) -> List[Tuple[int, int]]:
    """Merged +-flank genomic windows around every nominated coordinate."""
    spans = []
    donors = [j[0] for j in nominees]
    accs = [j[1] for j in nominees]
    spans.append((min(donors) - JE_FLANK_BP, max(donors) + JE_FLANK_BP))
    spans.append((min(accs) - JE_FLANK_BP, max(accs) + JE_FLANK_BP))
    spans.sort()
    merged = [spans[0]]
    for lo, hi in spans[1:]:
        if lo <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], hi))
        else:
            merged.append((lo, hi))
    return merged


def _exon_edges_in_windows(cand, windows, own_junction) -> Tuple[int, ...]:
    """Exon boundary coordinates inside ``windows`` excluding the
    candidate's own nominated junction pair (window-purity fingerprint)."""
    edges = set()
    for (d, a) in _chain(cand):
        for x in (d, a):
            if any(lo <= x <= hi for lo, hi in windows):
                edges.add(x)
    edges.discard(own_junction[0])
    edges.discard(own_junction[1])
    return tuple(sorted(edges))


# One krill aligner per WORKER PROCESS, not per contested interval. Building it
# is expensive (model load + GPU/CPU context) and `run_election` is called once
# per interval, so an uncached scorer rebuilt the whole signal stack for every
# contested locus. Keyed on the settings that actually change the aligner.
_SCORER_CACHE: Dict[tuple, object] = {}


def default_signal_scorer(config, num_thread: int = 0):
    """Production ballot scorer: krill eventalign mean NLL.

    Returns ``scorer(read_id, read_seq, cand, windows) -> (nll, n_events)``
    or None when the signal stack is unavailable (the election then abstains
    everywhere rather than inventing a verdict from weaker evidence).

    Process-cached: the krill aligner and the per-candidate mappy cache are
    reused across intervals within a worker.
    """
    key = (getattr(config, "krill_pore", None), bool(getattr(config, "use_gpu", False)),
           getattr(config, "signal_path", None), num_thread)
    if key in _SCORER_CACHE:
        return _SCORER_CACHE[key]
    scorer = _build_signal_scorer(config, num_thread)
    _SCORER_CACHE[key] = scorer
    return scorer


def _build_signal_scorer(config, num_thread: int = 0):
    try:
        import krill
        import mappy

        from fin.scoring.krill_aligner import make_krill_aligner
        from fin.scoring.m2_junction_nll import read_cand_mean_nll
        from fin.scoring.mappy_preset import get_m1_preset
    except Exception as exc:  # noqa: BLE001 - optional signal stack
        logger.warning("junction election: signal stack unavailable: %s", exc)
        return None
    krill_aligner, eff_gpu = make_krill_aligner(
        krill, config.krill_pore, config.use_gpu,
        hmm_confidence=False, num_thread=num_thread,
    )
    if krill_aligner is None:
        logger.warning("junction election: krill aligner init failed")
        return None
    preset = get_m1_preset()
    # Bounded: a worker sees many intervals, and candidate ids are unique per
    # interval, so an unbounded cache would grow with the whole run.
    mappy_cache: "collections.OrderedDict[tuple, object]" = collections.OrderedDict()
    MAX_ALIGNERS = 256

    def scorer(read_id, read_seq, cand, windows):
        # Key on the SEQUENCE, not the candidate id: ids are unique per interval
        # but recur across intervals/runs, so an id-only key would hand a second
        # input's candidate the aligner built from the first one's sequence.
        key = (cand.candidate_id, hashlib.blake2b(
            (cand.sequence or "").encode("ascii", "ignore"),
            digest_size=8).hexdigest())
        aln = mappy_cache.get(key)
        if aln is None:
            aln = mappy.Aligner(seq=cand.sequence, preset=preset)
            mappy_cache[key] = aln
            if len(mappy_cache) > MAX_ALIGNERS:
                mappy_cache.popitem(last=False)
        else:
            mappy_cache.move_to_end(key)
        return read_cand_mean_nll(
            read_id, read_seq, cand, list(windows), krill_aligner, aln,
            config.signal_path, config.krill_pore,
            use_gpu=eff_gpu, num_thread=num_thread,
        )

    return scorer


def m1_scorer():
    """Ablation ballot: mappy alignment score (base-space), NLL-oriented
    (lower is better), margin in AS units."""
    import mappy

    from fin.scoring.mappy_preset import get_m1_preset
    from fin.scoring.mappy_score import score_hit

    preset = get_m1_preset()
    cache: Dict[str, object] = {}

    def scorer(read_id, read_seq, cand, windows):
        aln = cache.get(cand.candidate_id)
        if aln is None:
            aln = mappy.Aligner(seq=cand.sequence, preset=preset)
            cache[cand.candidate_id] = aln
        best = None
        for h in aln.map(read_seq):
            v = score_hit(h)
            if v is not None and (best is None or v > best):
                best = v
        if best is None:
            return float("nan"), 0
        return -float(best), JE_BALLOT_MIN_EVENTS  # AS as pseudo-NLL

    return scorer


def _shuffle_permutation(read_id: str, k: int) -> List[int]:
    """Deterministic per-read permutation of k hypothesis scores (negative
    control: destroys the read<->hypothesis signal mapping)."""
    digest = hashlib.sha1(read_id.encode("utf-8")).digest()
    order = list(range(k))
    # Fisher-Yates driven by the digest bytes (k is tiny).
    for i in range(k - 1, 0, -1):
        j = digest[i % len(digest)] % (i + 1)
        order[i], order[j] = order[j], order[i]
    return order


def run_junction_election(
    config,
    cand_list: Sequence,
    quant_results: Sequence,
    candidate_set,
    interval_name: str = "",
    scorer: Optional[Callable] = None,
) -> Tuple[set, List[ElectionRow]]:
    """Hold every nominated K-way election for one interval.

    Returns ``(drop_ids, rows)``. ``drop_ids`` is empty unless
    ``junction_election_mode == "require"`` AND an accused variant is
    ``defeated``.
    """
    mode = getattr(config, "junction_election_mode", "off")
    ballot = getattr(config, "junction_election_ballot", "signal")
    burden = getattr(config, "junction_election_burden", "prosecution")
    defense = burden == "defense"
    junction_support = getattr(candidate_set, "junction_support", None) or {}
    quant_by_id = {q.candidate_id: q for q in quant_results}
    contests, truncated = nominate_contests(
        cand_list, quant_by_id, junction_support)
    if truncated:
        logger.warning(
            "junction election interval %s: %d contest(s) beyond the "
            "JE_MAX_CONTESTS=%d cap were NOT held (loud truncation; raise "
            "the cap if this occurs on real data)",
            interval_name, truncated, JE_MAX_CONTESTS,
        )
    if not contests:
        return set(), []
    shuffle = ballot == "shuffle"
    if scorer is None and ballot in ("signal", "shuffle"):
        scorer = default_signal_scorer(config)
    elif scorer is None and ballot == "m1":
        scorer = m1_scorer()
    margin = float(getattr(config, "m2_tiebreak_margin", 1e-9))
    read_spans = getattr(candidate_set, "read_spans", None) or {}
    read_seqs = getattr(candidate_set, "read_sequences", None) or {}
    rows: List[ElectionRow] = []
    drop_ids: set = set()
    for ci, con in enumerate(contests):
        nominees = con["nominees"]
        reps = con["reps"]
        accused_map = con["accused"]
        windows = _contest_windows(nominees)
        contest_id = f"{interval_name}#{ci}"
        best_comp_clean = max(
            (con["clean"][j] for j in nominees if j not in accused_map),
            default=0,
        )

        def make_row(j):
            return ElectionRow(
                interval=interval_name, contest_id=contest_id,
                accused_junction=_jstr(j),
                nominees=",".join(_jstr(x) for x in nominees),
                accused_clean=con["clean"][j],
                best_competitor_clean=best_comp_clean,
                accused_ids=",".join(
                    c.candidate_id for c in accused_map[j]),
            )

        contest_rows = {j: make_row(j) for j in sorted(accused_map)}
        rows.extend(contest_rows.values())

        def _screen_drop(reason: str) -> None:
            """defense burden: the screen stands unless ballots exonerate."""
            for j, r in contest_rows.items():
                r.verdict = reason
                if defense and mode == "require":
                    r.action = "drop"
                    drop_ids.update(c.candidate_id for c in accused_map[j])

        # Window purity: inside the ballot windows every representative must
        # agree on all exon boundaries other than its own nominated variant.
        fingerprints = {
            _exon_edges_in_windows(reps[j], windows, j) for j in nominees
        }
        if len(fingerprints) > 1:
            _screen_drop("window_impure")
            continue

        if ballot == "screen":
            # Ablation: the nomination screen acts alone (no ballots).
            for j, r in contest_rows.items():
                r.verdict = VERDICT_DEFEATED
                if mode == "require":
                    r.action = "drop"
                    drop_ids.update(c.candidate_id for c in accused_map[j])
            continue

        lo = min(w[0] for w in windows)
        hi = max(w[1] for w in windows)
        electorate = set()
        for j in nominees:
            for cand in con["carriers"].get(j, (reps[j],)):
                qr = quant_by_id.get(cand.candidate_id)
                for rid in (getattr(qr, "assigned_read_ids", ()) or ()):
                    sp = read_spans.get(rid)
                    if sp is not None and sp[0] <= lo and sp[1] >= hi:
                        electorate.add(rid)
        # Deterministic unbiased sampling: order by sha1(read_id).
        electorate = sorted(
            electorate,
            key=lambda r: hashlib.sha1(r.encode("utf-8")).hexdigest(),
        )[:JE_MAX_ELECTORATE]
        for r in contest_rows.values():
            r.n_electorate = len(electorate)
        if scorer is None or not electorate:
            _screen_drop("no_ballots")
            continue

        votes = {j: 0 for j in nominees}
        for rid in electorate:
            seq = read_seqs.get(rid)
            if not seq:
                continue
            scored = []
            for j in nominees:
                nll, n_ev = scorer(rid, seq, reps[j], windows)
                if n_ev >= JE_BALLOT_MIN_EVENTS and nll == nll:
                    scored.append((nll, j))
            if len(scored) < 2:
                continue
            if shuffle:
                vals = [s[0] for s in scored]
                order = _shuffle_permutation(rid, len(scored))
                scored = [(vals[order[i]], scored[i][1])
                          for i in range(len(scored))]
            scored.sort()
            if scored[1][0] - scored[0][0] < margin:
                continue                     # too close to call: abstain
            votes[scored[0][1]] += 1
        valid = sum(votes.values())
        winner = max(nominees, key=lambda j: (votes[j], j))
        for j, r in contest_rows.items():
            r.n_ballots = valid
            r.votes_accused = votes[j]
            r.votes_winner = votes[winner]
            r.winner_junction = _jstr(winner) if valid else ""
            upheld = (
                valid >= JE_QUORUM
                and winner == j
                and votes[winner] / valid >= JE_MAJORITY
            )
            convicted = (
                valid >= JE_QUORUM
                and winner != j
                and votes[winner] / valid >= JE_MAJORITY
                and votes[j] / valid <= 1 - JE_MAJORITY
            )
            if upheld:
                r.verdict = VERDICT_UPHELD
            elif convicted:
                r.verdict = VERDICT_DEFEATED
            elif valid < JE_QUORUM:
                r.verdict = "no_quorum"
            else:
                r.verdict = "split"
            should_drop = (
                r.verdict == VERDICT_DEFEATED
                if not defense
                else r.verdict != VERDICT_UPHELD
            )
            if mode == "require" and should_drop:
                r.action = "drop"
                drop_ids.update(c.candidate_id for c in accused_map[j])
    return drop_ids, rows


def write_election_tsv(path, rows: Sequence[ElectionRow]) -> None:
    """Deterministic shard/artifact writer."""
    ordered = sorted(
        rows, key=lambda r: (r.interval, r.contest_id, r.accused_junction))
    with open(path, "w") as fh:
        fh.write("\t".join(ELECTION_TSV_FIELDS) + "\n")
        for r in ordered:
            fh.write("\t".join(r.as_tsv_row()) + "\n")


def merge_election_shards(work_dir) -> Optional[str]:
    """Merge per-interval worker shards into ``junction_election.tsv``.

    Globally re-sorted (worker completion order cannot perturb the artifact);
    header-only when no contest was nominated anywhere. Shards are removed
    after a successful merge. Returns None when ``work_dir`` is missing.
    """
    import csv
    import shutil
    from pathlib import Path as _Path

    wd = _Path(work_dir)
    if not wd.is_dir():
        return None
    shard_dir = wd / "_jelection_shards"
    rows: List[dict] = []
    if shard_dir.is_dir():
        for shard in sorted(shard_dir.glob("*.tsv")):
            with shard.open() as fh:
                rows.extend(csv.DictReader(fh, delimiter="\t"))
    rows.sort(key=lambda r: (
        r["interval"], r["contest_id"], r["accused_junction"]))
    out = wd / "junction_election.tsv"
    tmp = out.with_suffix(".tsv.tmp")
    with tmp.open("w") as fh:
        fh.write("\t".join(ELECTION_TSV_FIELDS) + "\n")
        for r in rows:
            fh.write("\t".join(r.get(f, "") for f in ELECTION_TSV_FIELDS) + "\n")
    tmp.replace(out)
    if shard_dir.is_dir():
        shutil.rmtree(shard_dir)
    return str(out)
