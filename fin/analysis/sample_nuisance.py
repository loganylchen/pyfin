"""Sample-level nuisance models, estimated WITHIN the sample without truth.

The generative model of Part A conditions every channel likelihood on
sample-specific nuisance parameters ``eta_s``: how precisely this run places a
3' end, how often a polyA call is available, how the 5' end degrades, how
junctions err. None of these may be a constant carried over from a benchmark;
each is an instrument/sample property and is estimated from the sample's own
reads. This module holds the ones the existence test consumes.

Endpoint model (the 3'/polyA channel, A7 of the spec)
-----------------------------------------------------
For a polyA-PASS read ``r`` and candidate ``c`` with claimed TES ``t_c``, the
strand-aware signed distance ``d = o_c(e_r) - o_c(t_c)`` (``> 0`` = the read runs
PAST the claimed end) is modelled as a two-component mixture

    P(d) = (1 - eta) g(d) + eta b(d)

``g`` is the narrow endpoint-jitter kernel of a true TES; ``b`` is the broad
error kernel (chimeras, mapping truncation, read-through). Both are discrete
Laplace densities; ``eta`` is the broad weight. All three are fitted by EM, from
the distances to the INPUT annotation's claimed TES positions when at least two
of those exist (however far away they are -- there is no proximity filter), and
otherwise -- no annotation, or fewer than two reads on an annotated locus -- from
split-half distances: modes are built from one half of the sample's admitted 3'
ends and distances measured from the other half, so no read is measured against a
mode it helped define (which under-states the jitter and then penalises any TES
coordinate a few bp from the read cluster). Neither path sees a candidate, a
truth set or a condition name.

Measured before this was written (three grids, 2026-09-05): relative to the
nearest same-strand annotated TES, PASS reads have MAD 1 bp, q90 <= 8 bp, and a
> 25 bp read-through rate of 0.2-0.4 %. The kernel is the same on every grid;
that is why it is a nuisance parameter and not a tuning knob.

Abstention: by default a read that is not polyA-PASS has no 3' term (0 for every
candidate). Which verdicts carry a 3' term is what
``polya_end_accept_adapter`` / ``endpoint_qc_conditional`` change; see
:func:`fin.analysis.molecule_evidence.polya_end_evidential` for the two rules.
A sample where krill scored nothing has no endpoint model, and the channel is
reported unavailable rather than silently neutral.
"""
from __future__ import annotations

import bisect
import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

# A mode must be defined by more than the read being scored against it, so it
# needs >= 2 reads; two ends within this many bp are the same mode. This is the
# resolution at which two distinct polyA sites are distinguishable at all --
# below the measured jitter MAD there is no information -- not a tuning value.
MODE_MERGE_BP = 2
# A read's 5' start is evidence about a TRANSCRIPT's start only when other
# molecules end there too. The 3' channel gets its specificity from the polyA
# call -- a PASS says this molecule was polyadenylated at this position -- and
# the 5' side has no such per-read statement: a single read starting upstream of
# an annotated TSS is as likely to be read-through transcription, an unspliced
# precursor or a molecule from an upstream gene as it is to be a real
# alternative start. Measured on human dRNA (H9, chr1/19/21, GENCODE v44):
# 12.9% of reads start >50 bp upstream of the nearest annotated TSS, at a MEDIAN
# distance of 1232 bp, and 60% of those are alone at their start position while
# only 26.8% are shared by >= 3 reads. Requiring a start CLUSTER is the 5'
# analogue of requiring a polyA PASS.
TSS_CLUSTER_TOL_BP = 10

# A kernel is USABLE evidence only when it is a calibrated model of end
# jitter: enough ends to estimate three mixture parameters, and a narrow
# component that IS narrow. The narrow scale models the technical precision of
# a polyA-anchored 3' end / an at-TSS 5' end -- measured 0.8-2.5 bp on three
# datasets -- so a fit whose "narrow" component absorbs tens of bp has
# collapsed onto a bimodal or sparse sample (e.g. half the ends 100 bp off) and
# must abstain rather than call every 100-bp discrepancy jitter.
MIN_KERNEL_FIT = 50
NARROW_SCALE_MAX_BP = 25.0
BROAD_SCALE_MAX_BP = 500.0
# Cap on the broad weight. Fitted on ONE pooled population the cap is a
# regularizer -- it ASSUMES a sample is not more than half read-through, which is
# not something measured here, and it stops a bimodal pooled sample from fitting
# its way out of the channel entirely. Fitted PER
# polyA QC CLASS the cap instead FORCES a class to claim at least half its
# ends are jitter, which is false for a class that is mostly wrong: SIRV4
# NOREGION has only 34.5% of its ends within 20 bp of an annotated TES yet fits
# a narrow 1.8 bp component, so capped it claims 50% narrow mass and released it
# fits eta 0.695 (0.305 narrow weight per read). Releasing the cap therefore
# lets a weak class report a smaller within-read contrast. It is NOT what makes a
# class abstain: that is `available`, which reads n_fit and the narrow SCALE and
# never eta -- though the cap can move b_narrow, so the two are separable in
# practice rather than by construction. On the four grids measured, every PASS /
# ADAPTER / READ_FAILED class already fitted eta below 0.5 and is unchanged to
# three decimals; only the NOREGION classes move.
ETA_MAX = 0.5
ETA_MAX_PER_CLASS = 0.99


@dataclass(frozen=True)
class EndpointModel:
    """Fitted 3'-endpoint mixture: narrow Laplace(b_narrow), broad Laplace(b_broad), weight eta."""

    b_narrow: float
    b_broad: float
    eta: float
    n_fit: int
    version: str = "endpoint-laplace-mixture-v1"

    @staticmethod
    def _laplace(d: float, b: float) -> float:
        # discrete Laplace pmf on integers: (1-q)/(1+q) q^|d|, q = exp(-1/b)
        q = math.exp(-1.0 / max(b, 1e-6))
        return (1.0 - q) / (1.0 + q) * q ** abs(d)

    def loglik(self, d: float) -> float:
        p = ((1.0 - self.eta) * self._laplace(d, self.b_narrow)
             + self.eta * self._laplace(d, self.b_broad))
        return math.log(max(p, 1e-300))

    def narrow_posterior(self, d: float) -> float:
        """P(narrow | d): how much of this distance the true-TES kernel explains."""
        n = (1.0 - self.eta) * self._laplace(d, self.b_narrow)
        b = self.eta * self._laplace(d, self.b_broad)
        return n / max(n + b, 1e-300)

    @property
    def available(self) -> bool:
        return self.n_fit >= MIN_KERNEL_FIT and self.b_narrow <= NARROW_SCALE_MAX_BP


def _end_modes(ends: Sequence[int]) -> List[Tuple[int, int]]:
    """Cluster sorted end positions into (position, count) modes."""
    modes: List[Tuple[int, int]] = []
    for e in sorted(ends):
        if modes and e - modes[-1][0] <= MODE_MERGE_BP:
            pos, n = modes[-1]
            # running mean keeps the mode where the mass is
            modes[-1] = (int(round((pos * n + e) / (n + 1))), n + 1)
        else:
            modes.append((e, 1))
    return modes


def endpoint_distances(
    ends_by_locus: Mapping[Tuple[str, str], Sequence[int]],
    read_ids_by_locus: Optional[Mapping[Tuple[str, str], Sequence[str]]] = None,
) -> List[int]:
    """Split-half distances of PASS ends to the nearest mode of the OTHER half.

    Measuring a read's distance to a mode its own end helped define
    under-states the jitter (a mode of two reads sits between them), and the
    kernel then punishes any TES coordinate that is a few bp from the read
    cluster -- which is exactly where an annotated TES tends to sit relative to
    the observed cleavage. Reads are therefore split in two deterministic halves
    (by read-id hash, else alternating order): modes are built from one half,
    distances measured from the other, in both directions. Distances are signed
    in transcript orientation (positive = past the mode).
    """
    return [d for _rid, d in endpoint_distance_pairs(ends_by_locus, read_ids_by_locus)]


def endpoint_distance_pairs(
    ends_by_locus: Mapping[Tuple[str, str], Sequence[int]],
    read_ids_by_locus: Optional[Mapping[Tuple[str, str], Sequence[str]]] = None,
) -> List[Tuple[str, int]]:
    """:func:`endpoint_distances`, keeping each distance's read id so the fit can
    be conditioned on the read's polyA QC class. Order is unchanged."""
    import hashlib

    out: List[Tuple[str, int]] = []
    for key, ends in ends_by_locus.items():
        chrom, strand = key
        ids = list(read_ids_by_locus.get(key, ())) if read_ids_by_locus else []
        halves = ([], [])
        for i, e in enumerate(ends):
            if i < len(ids):
                h = int.from_bytes(hashlib.blake2b(ids[i].encode("utf-8"), digest_size=4).digest(), "big") & 1
                rid = ids[i]
            else:
                h = i & 1
                rid = ""
            halves[h].append((rid, e))
        for a, b in ((0, 1), (1, 0)):
            modes = [m for m in _end_modes([e for _r, e in halves[a]]) if m[1] >= 2]
            if not modes:
                continue
            centres = [m[0] for m in modes]
            for rid, e in halves[b]:
                k = bisect.bisect_left(centres, e)
                cands = [centres[j] for j in (k - 1, k) if 0 <= j < len(centres)]
                d = min((e - c for c in cands), key=abs)
                out.append((rid, d if strand == "+" else -d))
    return out


def fit_endpoint_model(distances: Sequence[int], *, max_iter: int = 200,
                       tol: float = 1e-8, eta_max: float = ETA_MAX) -> EndpointModel:
    """EM for the two-Laplace mixture on a list of signed endpoint distances.

    The caller decides what those are: distances to annotated TES positions,
    split-half distances to the sample's own read modes, or (through
    :func:`fit_start_model`) upstream 5' distances."""
    d = [int(x) for x in distances]
    n = len(d)
    if n < 2:
        return EndpointModel(b_narrow=1.0, b_broad=100.0, eta=0.0, n_fit=n)
    absd = sorted(abs(x) for x in d)
    med = absd[n // 2]
    b_n = max(float(med), 0.5)
    b_b = max(float(absd[int(0.95 * (n - 1))]), 10.0 * b_n, 20.0)
    eta = 0.05
    for _ in range(max_iter):
        w_sum = 0.0; s_n = 0.0; s_b = 0.0; n_n = 0.0; n_b = 0.0
        for x in d:
            gn = (1 - eta) * EndpointModel._laplace(x, b_n)
            gb = eta * EndpointModel._laplace(x, b_b)
            z = gn + gb
            if z <= 0:
                continue
            r = gb / z
            w_sum += r
            s_n += (1 - r) * abs(x); n_n += (1 - r)
            s_b += r * abs(x); n_b += r
        new_eta = min(max(w_sum / n, 1e-4), eta_max)
        new_bn = max(s_n / max(n_n, 1e-9), 0.5)
        # The broad component describes read-through / chimeric / truncation
        # errors, whose scale is that of a transcript (hundreds of bp), not the
        # distance to the next gene: an unbounded scale makes b(d) flat and the
        # channel uninformative at exactly the distances it exists to judge.
        new_bb = min(max(s_b / max(n_b, 1e-9), 2.0 * new_bn), BROAD_SCALE_MAX_BP)
        delta = max(abs(new_eta - eta), abs(new_bn - b_n), abs(new_bb - b_b))
        eta, b_n, b_b = new_eta, new_bn, new_bb
        if delta < tol:
            break
    return EndpointModel(b_narrow=b_n, b_broad=b_b, eta=eta, n_fit=n)


@dataclass(frozen=True)
class StartModel:
    """One-sided 5'-start model (spec A6, the part that needs no hazard).

    A read cannot begin UPSTREAM of the transcript it came from: 5' degradation
    only shortens. So for a read whose 5' end lies ``u > 0`` bp upstream of a
    candidate's claimed TSS, the candidate explains it only through the error
    kernel (mapping overhang, TSS mis-annotation, chimera). Reads starting at or
    downstream of the claimed TSS abstain here (their hazard term is the
    degradation channel proper, integrated separately). The upstream distances
    are modelled as a jitter/error Laplace mixture fitted within the sample
    against the input annotation's claimed TSS positions.
    """

    b_narrow: float
    b_broad: float
    eta: float
    n_fit: int
    version: str = "start-upstream-mixture-v1"

    def loglik_upstream(self, u: float) -> float:
        """log P(read starts u bp UPSTREAM of the claimed TSS | candidate), u > 0."""
        p = ((1.0 - self.eta) * EndpointModel._laplace(u, self.b_narrow)
             + self.eta * EndpointModel._laplace(u, self.b_broad))
        return math.log(max(p, 1e-300))

    def narrow_posterior(self, u: float) -> float:
        """P(jitter | u): how much of an upstream distance the at-TSS kernel explains."""
        n = (1.0 - self.eta) * EndpointModel._laplace(u, self.b_narrow)
        b = self.eta * EndpointModel._laplace(u, self.b_broad)
        return n / max(n + b, 1e-300)

    @property
    def available(self) -> bool:
        return self.n_fit >= MIN_KERNEL_FIT and self.b_narrow <= NARROW_SCALE_MAX_BP


def fit_start_model(upstream_distances: Sequence[int]) -> StartModel:
    m = fit_endpoint_model(upstream_distances)
    return StartModel(b_narrow=m.b_narrow, b_broad=m.b_broad, eta=m.eta, n_fit=m.n_fit)


@dataclass(frozen=True)
class SampleNuisanceModel:
    endpoint: EndpointModel
    polya_call_rate: float          # rho_s: evidential among scorable reads
    n_scorable: int
    start: Optional[StartModel] = None
    # One endpoint kernel per krill polyA QC class, when the evidence was
    # derived class-conditionally. A TUPLE, not a dict: this dataclass is frozen
    # and callers may use a fitted model as a dict key or a cache argument, which
    # a dict field would break. `qc_conditional` records WHICH mode produced this
    # fit, so an empty mapping (no class qualified) is not mistaken for "pooled
    # mode" -- otherwise this model and `apply_endpoint_likelihood` would
    # disagree about whether the pooled kernel may be used.
    endpoint_kernels: Tuple[Tuple[str, EndpointModel], ...] = ()
    qc_conditional: bool = False
    version: str = "sample-nuisance-v1"

    @property
    def endpoint_by_qc(self) -> Dict[str, EndpointModel]:
        return dict(self.endpoint_kernels)

    @property
    def polya_available(self) -> bool:
        return self.n_scorable > 0 and self.endpoint.available

    def endpoint_for(self, qc: Optional[str]) -> Optional[EndpointModel]:
        """The kernel that scores a read of QC class ``qc``, or None to abstain.

        Class-conditional: krill's verdict is a raw-signal reliability label, and
        each label's kernel is fitted from the sample's own ends, so how much a
        class is trusted is estimated per sample rather than fixed by a
        whitelist. A class abstains when its fit is not `available` -- fewer than
        MIN_KERNEL_FIT ends, or a "narrow" component wider than
        NARROW_SCALE_MAX_BP. That is a test of whether the class's ends CLUSTER
        tightly enough to be a calibrated jitter model, which is a necessary
        condition for carrying endpoint evidence and not a sufficient one; on the
        human and Sequin samples measured it is what removes NOREGION, on HEYA8
        NOREGION passes it (and its ends are in fact good there)."""
        if not self.qc_conditional:
            return self.endpoint if self.polya_available else None
        for name, m in self.endpoint_kernels:      # a handful of QC classes
            if name == qc:
                return m if m.available else None
        return None


@dataclass(frozen=True)
class NuisanceEvidence:
    """Per-read endpoint distances of ONE sample, derived once.

    Deriving the distances (nearest claimed TES / TSS per read) is sample work
    that does not depend on which interval is being judged; only the MIXTURE FIT
    does, through leave-one-interval-out exclusion. Separating the two keeps the
    per-interval cost at "filter a list and run the EM" instead of re-deriving
    every read's distance (and re-sorting the annotation) for every interval.
    """

    tes: Tuple[Tuple[str, int], ...]                 # (read_id, signed distance to claimed TES)
    tss: Tuple[Tuple[str, int], ...]                 # (read_id, upstream distance, >= 0)
    ends_by_locus: Mapping[Tuple[str, str], Tuple[int, ...]]     # evidential 3' ends (fallback)
    ids_by_locus: Mapping[Tuple[str, str], Tuple[str, ...]]
    n_scorable: int
    n_pass: int          # reads whose polyA call is 3'-end evidential
    has_claimed_tes: bool
    # read_id -> krill polyA QC verdict, present only when the evidence was
    # derived class-conditionally (every verdict admitted, the fit decides).
    qc_by_read: Mapping[str, str] = ()
    # which derivation produced this evidence. Recorded rather than inferred from
    # `qc_by_read`, which is empty both when the mode is off and when it is on
    # over a sample with no usable verdicts.
    qc_conditional: bool = False
    # In class-conditional mode, the evidence the NAME-based predicate alone
    # would have produced. The POOLED kernel is fitted from THIS, so it is
    # bit-for-bit the kernel it was before the widening -- not a re-fit on a
    # filtered list, whose float summation order would differ. Its consumer is
    # the discovery contradiction test, which asks an ABSOLUTE question
    # ("narrow_posterior(d) < 0.5") whose meaning is not invariant to eta: a
    # class fitted as LESS reliable has a lower narrow posterior everywhere, so
    # it would produce MORE contradictions and be protected by FEWER annotated
    # alternatives, which inverts the intent. The widening is an
    # ASSIGNMENT-channel change, and assignment asks a within-read contrast where
    # a lower-confidence class correctly moves the answer less.
    pooled: Optional["NuisanceEvidence"] = None
    # sorted 5' starts per (chrom, strand) of every primary read, so a candidate
    # upstream start can be asked how many molecules share it
    starts_by_locus: Mapping[Tuple[str, str], Tuple[int, ...]] = ()

    def start_support(self, chrom: str, strand: str, pos: int,
                      tol: int = TSS_CLUSTER_TOL_BP) -> int:
        """How many reads start within ``tol`` bp of ``pos`` on this locus."""
        arr = (self.starts_by_locus or {}).get((chrom, strand))
        if not arr:
            return 0
        return (bisect.bisect_right(arr, pos + tol)
                - bisect.bisect_left(arr, pos - tol))

    def start_share(self, chrom: str, strand: str, pos: int,
                    tol: int = TSS_CLUSTER_TOL_BP, window: int = 200) -> float:
        """Share of this locus's nearby starts that sit within ``tol`` of ``pos``.

        The per-LOCUS form of the question the sample-level gate asks: is this
        position a dominant start, or one shoulder of a smear? A count alone
        cannot tell those apart -- three reads is a cluster among five and noise
        among fifty -- and on 5'-degraded dRNA the second case is the common one.
        Returns 0.0 when nothing is within ``window``.
        """
        arr = (self.starts_by_locus or {}).get((chrom, strand))
        if not arr:
            return 0.0
        lo = bisect.bisect_left(arr, pos - window)
        hi = bisect.bisect_right(arr, pos + window)
        n_local = hi - lo
        if n_local <= 0:
            return 0.0
        return self.start_support(chrom, strand, pos, tol) / n_local


def annotated_tss_tracking(
    starts_by_locus: Mapping[Tuple[str, str], Sequence[int]],
    annotated_tss: Mapping[Tuple[str, str], Sequence[int]],
    *,
    min_reads: int = 3,
    window_bp: int = 200,
    mode_tol_bp: int = 10,
    track_bp: int = 25,
    min_mode_share: float = 0.5,
) -> Tuple[float, int]:
    """Fraction of annotated loci where this run RESOLVES a 5' end.

    Answers one question about a SAMPLE: does this run place 5' ends precisely
    enough for an outlier start cluster to be evidence? For every annotated TSS
    with at least ``min_reads`` read starts within ``window_bp``, find the
    densest start position (the mode of a ``mode_tol_bp`` window -- the mode,
    not the mean, so a degradation tail does not drag it). The locus counts as
    RESOLVED when EITHER

      * the mode is within ``track_bp`` of the annotated TSS, or
      * the mode holds at least ``min_mode_share`` of the local starts.

    The second clause is not slack, it removes a confound. A locus can fail to
    track for two opposite reasons: the run's 5' ends are imprecise (which is
    what this statistic is for), or the ANNOTATION's TSS is wrong -- and the
    second is exactly the case the 5' extent repair exists to fix. A tight
    cluster sitting away from the annotated TSS is evidence against the
    annotation, not against the reads. Measured: on a SIRV library the
    tracking-only fraction collapses 0.797 -> 0.438 when the annotation's 5'
    ends are truncated (`c_trunc5p50`) while this combined fraction is
    unchanged (0.957 -> 0.959), and a gate built on tracking alone therefore
    switched the repair off exactly where it was needed (-4.88 F1 on that
    family).

    Separation on the combined statistic, over 14 measurements: human dRNA
    0.824-0.871, SIRV4 0.918-0.971 (including two with a 5'-truncated
    annotation), Sequin 1.000. The gap is 0.871-0.918, i.e. 0.047 wide -- thin,
    and confounded with depth: the 0.871 human value is the lowest-depth BAM,
    and fewer local starts raise the share. These datasets were INSPECTED while
    designing this, so the separation is development-derived calibration, not a
    held-out result.

    The scoring unit is the annotated POSITION, deduplicated per locus, so a
    gene whose isoforms share a TSS votes once. A locus tracks when SOME densest
    window sits within ``track_bp`` of the annotated TSS: that is invariant under
    coordinate reflection, while a first-maximum tie-break makes the answer
    depend on which strand the locus lies on.

    What this measures is a LOCAL GENOMIC CONCENTRATION of read starts around a
    supplied annotated position -- not this transcript's start. Neighbouring
    promoters, overlapping genes and alternative isoforms contribute to the same
    window. The +-200 bp and +-10 bp scales are engineering choices, not
    established optima.

    Returns ``(fraction_resolved, n_loci_scored)``; ``(0.0, 0)`` when nothing
    could be scored, which callers must treat as "unknown", not as "unresolved".
    """
    if not starts_by_locus or not annotated_tss:
        return 0.0, 0
    tracked = scored = 0
    for key, positions in annotated_tss.items():
        arr = starts_by_locus.get(key)
        if not arr:
            continue
        # one vote per POSITION, not per transcript: annotated TSSs are shared
        # by the isoforms of a gene, so counting them per transcript would give
        # annotation-rich loci extra weight in the fraction
        for t in sorted(set(positions)):
            lo = bisect.bisect_left(arr, t - window_bp)
            hi = bisect.bisect_right(arr, t + window_bp)
            local = arr[lo:hi]
            if len(local) < min_reads:
                continue
            scored += 1
            counts = [(bisect.bisect_right(local, x + mode_tol_bp)
                       - bisect.bisect_left(local, x - mode_tol_bp))
                      for x in local]
            best_n = max(counts)
            # Ties must not be broken by genomic direction. Taking the FIRST
            # maximum makes the answer depend on which strand the locus lies on
            # -- the same geometry can score 0 on '+' and 1 on its mirror image
            # -- and taking the midpoint of all tied maxima lands between two
            # distant clusters, on nothing. The question is whether SOME densest
            # window sits on the annotated TSS, and that is invariant under
            # reflection because the positions and `t` reflect together.
            tied = [x for x, n in zip(local, counts) if n == best_n]
            if (any(abs(x - t) <= track_bp for x in tied)
                    or best_n / len(local) >= min_mode_share):
                tracked += 1
    return (tracked / scored if scored else 0.0), scored


def _nearest_signed(sorted_pos: Sequence[int], x: int) -> int:
    """Signed genomic distance from the nearest position in ``sorted_pos``."""
    k = bisect.bisect_left(sorted_pos, x)
    cands = [sorted_pos[j] for j in (k - 1, k) if 0 <= j < len(sorted_pos)]
    return min((x - t for t in cands), key=abs)


def sample_nuisance_evidence(
    polya_map: Mapping[str, Tuple[Optional[float], Optional[str]]],
    read_ends: Mapping[str, Tuple[str, str, int]],
    *,
    claimed_tes: Optional[Mapping[Tuple[str, str], Sequence[int]]] = None,
    read_starts: Optional[Mapping[str, Tuple[str, str, int]]] = None,
    claimed_tss: Optional[Mapping[Tuple[str, str], Sequence[int]]] = None,
    accept_adapter: bool = False,
    qc_conditional: bool = False,
) -> NuisanceEvidence:
    """Derive every read's endpoint distances once (annotation sorted once per locus).

    ``accept_adapter`` widens the 3' population to krill's ADAPTER calls;
    ``qc_conditional`` admits every krill VERDICT and records each read's class so
    the kernel can be fitted per class, and records which of them the NAME-based
    predicate would have admitted so the pooled kernel keeps its own population.
    The fit and the channel share one ADMISSION PREDICATE, so they cannot disagree
    about which polyA calls carry 3' evidence. They do not see identical read
    SETS: leave-one-interval-out exclusion, reads with no end coordinate, loci
    with no annotation, split-half mode eligibility and per-class availability all
    change who reaches the fit."""
    from fin.analysis.molecule_evidence import polya_end_evidential

    by_locus: Dict[Tuple[str, str], List[int]] = {}
    ids_by_locus: Dict[Tuple[str, str], List[str]] = {}
    qc_by_read: Dict[str, str] = {}
    n_scorable = n_pass = 0
    for rid in sorted(polya_map):
        n_scorable += 1
        if not polya_end_evidential(polya_map[rid], accept_adapter, qc_conditional):
            continue
        n_pass += 1
        if qc_conditional:
            qc_by_read[rid] = polya_map[rid][1]
        loc = read_ends.get(rid)
        if loc is None:
            continue
        chrom, strand, e = loc
        by_locus.setdefault((chrom, strand), []).append(int(e))
        ids_by_locus.setdefault((chrom, strand), []).append(rid)

    tes_pairs: List[Tuple[str, int]] = []
    if claimed_tes:
        sorted_tes = {k: sorted(v) for k, v in claimed_tes.items()}
        for key, ends in by_locus.items():
            pos = sorted_tes.get(key)
            if not pos:
                continue
            strand = key[1]
            for rid, e in zip(ids_by_locus[key], ends):
                d = _nearest_signed(pos, e)
                tes_pairs.append((rid, d if strand == "+" else -d))

    tss_pairs: List[Tuple[str, int]] = []
    if read_starts and claimed_tss:
        # signed 5' distance: u > 0 = the read starts UPSTREAM of the nearest
        # claimed TSS (transcript orientation). Fitted on ALL primary reads (no
        # polyA needed): the 5' end is what the read has, degraded or not.
        sorted_tss = {k: sorted(v) for k, v in claimed_tss.items()}
        for rid in sorted(read_starts):
            chrom, strand, s5 = read_starts[rid]
            pos = sorted_tss.get((chrom, strand))
            if not pos:
                continue
            d = _nearest_signed(pos, s5)
            u = -d if strand == "+" else d          # upstream is negative genomic on +
            # One-sided: a degraded read (downstream start) is not an upstream
            # observation; it contributes 0 (at/inside the model), so the narrow
            # component is the at-TSS mass and the broad one the upstream errors.
            tss_pairs.append((rid, max(u, 0)))

    starts_by_locus: Dict[Tuple[str, str], List[int]] = {}
    for rid in sorted(read_starts or {}):
        chrom, strand, s5 = read_starts[rid]
        starts_by_locus.setdefault((chrom, strand), []).append(int(s5))
    return NuisanceEvidence(
        tes=tuple(tes_pairs), tss=tuple(tss_pairs),
        ends_by_locus={k: tuple(v) for k, v in by_locus.items()},
        ids_by_locus={k: tuple(v) for k, v in ids_by_locus.items()},
        n_scorable=n_scorable, n_pass=n_pass, has_claimed_tes=bool(claimed_tes),
        qc_by_read=qc_by_read, qc_conditional=bool(qc_conditional),
        # the same derivation under the name-based predicate, so the pooled
        # kernel is computed by the identical code path on the identical inputs
        pooled=(sample_nuisance_evidence(polya_map, read_ends,
                                         claimed_tes=claimed_tes,
                                         accept_adapter=accept_adapter)
                if qc_conditional else None),
        starts_by_locus={k: tuple(sorted(v)) for k, v in starts_by_locus.items()},
    )


def _surviving_pairs(evidence: NuisanceEvidence,
                     excl: Iterable[str]) -> List[Tuple[str, int]]:
    """(read_id, endpoint distance) for the reads left after exclusion, from the
    annotation when there is one and from the split-half read modes otherwise."""
    excl = set(excl)
    pairs = [(rid, d) for rid, d in evidence.tes if rid not in excl]
    if evidence.has_claimed_tes and len(pairs) >= 2:
        return pairs
    # With fewer than two surviving read-to-annotation distances -- no
    # annotation, or an annotation no surviving read is measured against -- the
    # split-half read-mode distances are the only reference; they depend on the
    # surviving read set, so they are recomputed per fit.
    ends, ids = {}, {}
    for key, locus_ids in evidence.ids_by_locus.items():
        keep = [(rid, e) for rid, e in zip(locus_ids, evidence.ends_by_locus[key])
                if rid not in excl]
        if keep:
            ids[key] = [rid for rid, _ in keep]
            ends[key] = [e for _, e in keep]
    return endpoint_distance_pairs(ends, ids)


def _fit_pooled(evidence: NuisanceEvidence, excl: Iterable[str]) -> EndpointModel:
    return fit_endpoint_model([d for _rid, d in _surviving_pairs(evidence, excl)])


def fit_from_evidence(
    evidence: NuisanceEvidence, *, exclude_reads: Iterable[str] = (),
) -> SampleNuisanceModel:
    """Fit the mixtures on the evidence minus ``exclude_reads`` (leave-one-interval-out:
    the interval under test contributes nothing to the nuisance it is judged against)."""
    excl = set(exclude_reads)
    model = _fit_pooled(evidence.pooled or evidence, excl)
    by_qc: List[Tuple[str, EndpointModel]] = []
    if evidence.qc_by_read:
        # deriving the surviving distances is not free -- in the annotation-free
        # branch it rebuilds the split-half modes -- and only the per-class fits
        # need the WIDENED set, so it stays inside this branch
        pairs = _surviving_pairs(evidence, excl)
        # One kernel per QC class over the WIDENED distances (the pooled kernel
        # above stays on the name-admitted population), with the broad weight
        # free (ETA_MAX_PER_CLASS) so a class that is mostly read-through can say
        # so instead of being pinned at a cap that exists to protect a POOLED fit
        # from turning itself off.
        grouped: Dict[str, List[int]] = {}
        for rid, d in pairs:
            qc = evidence.qc_by_read.get(rid)
            if qc:
                grouped.setdefault(qc, []).append(d)
        for qc, ds in grouped.items():
            if len(ds) >= 2:
                by_qc.append((qc, fit_endpoint_model(ds, eta_max=ETA_MAX_PER_CLASS)))
    ups = [u for rid, u in evidence.tss if rid not in excl]
    start_model = fit_start_model(ups) if len(ups) >= 2 else None
    return SampleNuisanceModel(
        endpoint=model,
        endpoint_kernels=tuple(by_qc),
        qc_conditional=evidence.qc_conditional,
        polya_call_rate=(evidence.n_pass / evidence.n_scorable) if evidence.n_scorable else 0.0,
        n_scorable=evidence.n_scorable,
        start=start_model,
    )


def fit_sample_nuisance(
    polya_map: Mapping[str, Tuple[Optional[float], Optional[str]]],
    read_ends: Mapping[str, Tuple[str, str, int]],
    *,
    exclude_reads: Iterable[str] = (),
    claimed_tes: Optional[Mapping[Tuple[str, str], Sequence[int]]] = None,
    read_starts: Optional[Mapping[str, Tuple[str, str, int]]] = None,
    claimed_tss: Optional[Mapping[Tuple[str, str], Sequence[int]]] = None,
    accept_adapter: bool = False,
    qc_conditional: bool = False,
) -> SampleNuisanceModel:
    """Derive the sample's endpoint distances and fit the mixtures in one call.

    ``read_ends``: read_id -> (chrom, strand, three_prime_genomic).

    ``claimed_tes`` (chrom, strand) -> annotated TES positions of the INPUT
    annotation. The kernel scores read ends against CLAIMED TES coordinates, so
    it must be fitted to that kind of distance: an annotation's TES sits a few
    bp from the observed cleavage cluster, and a kernel fitted read-vs-read
    (0.7 bp) would penalise every annotated model for that convention offset.
    Corrupted entries in the input annotation land in the broad component.
    Without an annotation the split-half read-mode distances are used.

    Callers that fit repeatedly over the same sample (once per interval, with
    different exclusions) should derive :func:`sample_nuisance_evidence` ONCE and
    call :func:`fit_from_evidence` instead.
    """
    return fit_from_evidence(
        sample_nuisance_evidence(
            polya_map, read_ends, claimed_tes=claimed_tes,
            read_starts=read_starts, claimed_tss=claimed_tss,
            accept_adapter=accept_adapter, qc_conditional=qc_conditional),
        exclude_reads=exclude_reads,
    )


def tes_polya_loglik(molecule, candidate, endpoint: EndpointModel) -> Optional[float]:
    """The 3'/polyA channel term ``l^{3'/A}_{rc}`` for one (read, candidate).

    None = abstain: the read has no PASS call, or the model is unavailable.
    Otherwise ``log[(1-eta) g(d) + eta b(d)]`` with ``d`` the strand-aware
    signed distance from the candidate's claimed TES to the read's 3' end
    (positive = the read runs past the claimed end).
    """
    if not endpoint.available or not getattr(molecule, "polya_available", False):
        return None
    strand = getattr(candidate, "strand", "+")
    t_c = candidate.end if strand == "+" else candidate.start
    e_r = molecule.three_prime_pos
    d = (e_r - t_c) if strand == "+" else (t_c - e_r)
    return endpoint.loglik(d)


def tes_distance(molecule, candidate) -> int:
    strand = getattr(candidate, "strand", "+")
    t_c = candidate.end if strand == "+" else candidate.start
    e_r = molecule.three_prime_pos
    return (e_r - t_c) if strand == "+" else (t_c - e_r)
