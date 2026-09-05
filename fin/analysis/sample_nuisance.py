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
Laplace densities; ``eta`` is the broad weight. All three are fitted by EM on
split-half distances: modes are built from one half of the sample's
polyA-PASS 3' ends and distances measured from the other half -- no annotation,
no candidate, no condition name, and no read measured against a mode it helped
define (which under-states the jitter and then penalises any TES coordinate a
few bp from the read cluster).

Measured before this was written (three grids, 2026-09-05): relative to the
nearest same-strand annotated TES, PASS reads have MAD 1 bp, q90 <= 8 bp, and a
> 25 bp read-through rate of 0.2-0.4 %. The kernel is the same on every grid;
that is why it is a nuisance parameter and not a tuning knob.

Abstention: a read that is not PASS has no 3' term (0 for every candidate).
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
    import hashlib

    out: List[int] = []
    for key, ends in ends_by_locus.items():
        chrom, strand = key
        ids = list(read_ids_by_locus.get(key, ())) if read_ids_by_locus else []
        halves = ([], [])
        for i, e in enumerate(ends):
            if i < len(ids):
                h = int.from_bytes(hashlib.blake2b(ids[i].encode("utf-8"), digest_size=4).digest(), "big") & 1
            else:
                h = i & 1
            halves[h].append(e)
        for a, b in ((0, 1), (1, 0)):
            modes = [m for m in _end_modes(halves[a]) if m[1] >= 2]
            if not modes:
                continue
            centres = [m[0] for m in modes]
            for e in halves[b]:
                k = bisect.bisect_left(centres, e)
                cands = [centres[j] for j in (k - 1, k) if 0 <= j < len(centres)]
                d = min((e - c for c in cands), key=abs)
                out.append(d if strand == "+" else -d)
    return out


def fit_endpoint_model(distances: Sequence[int], *, max_iter: int = 200,
                       tol: float = 1e-8) -> EndpointModel:
    """EM for the two-Laplace mixture on split-half endpoint distances."""
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
        new_eta = min(max(w_sum / n, 1e-4), 0.5)
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
    polya_call_rate: float          # rho_s: PASS among scorable reads
    n_scorable: int
    start: Optional[StartModel] = None
    version: str = "sample-nuisance-v1"

    @property
    def polya_available(self) -> bool:
        return self.n_scorable > 0 and self.endpoint.available


def fit_sample_nuisance(
    polya_map: Mapping[str, Tuple[Optional[float], Optional[str]]],
    read_ends: Mapping[str, Tuple[str, str, int]],
    *,
    exclude_reads: Iterable[str] = (),
    claimed_tes: Optional[Mapping[Tuple[str, str], Sequence[int]]] = None,
    read_starts: Optional[Mapping[str, Tuple[str, str, int]]] = None,
    claimed_tss: Optional[Mapping[Tuple[str, str], Sequence[int]]] = None,
) -> SampleNuisanceModel:
    """Fit the sample model from krill polyA calls and strand-aware 3' ends.

    ``read_ends``: read_id -> (chrom, strand, three_prime_genomic). ``exclude_reads``
    implements leave-one-interval-out: the interval under test contributes
    nothing to the nuisance it is judged against.

    ``claimed_tes`` (chrom, strand) -> annotated TES positions of the INPUT
    annotation. The kernel scores read ends against CLAIMED TES coordinates, so
    it must be fitted to that kind of distance: an annotation's TES sits a few
    bp from the observed cleavage cluster, and a kernel fitted read-vs-read
    (0.7 bp) would penalise every annotated model for that convention offset.
    Corrupted entries in the input annotation land in the broad component.
    Without an annotation the split-half read-mode distances are used.
    """
    from fin.analysis.molecule_evidence import POLYA_PASS, polya_state

    excl = set(exclude_reads)
    by_locus: Dict[Tuple[str, str], List[int]] = {}
    ids_by_locus: Dict[Tuple[str, str], List[str]] = {}
    n_scorable = n_pass = 0
    for rid, call in polya_map.items():
        n_scorable += 1
        qc, _ = polya_state(call)
        if qc != POLYA_PASS:
            continue
        n_pass += 1
        if rid in excl:
            continue
        loc = read_ends.get(rid)
        if loc is None:
            continue
        chrom, strand, e = loc
        by_locus.setdefault((chrom, strand), []).append(int(e))
        ids_by_locus.setdefault((chrom, strand), []).append(rid)
    if claimed_tes:
        dists: List[int] = []
        for key, ends in by_locus.items():
            tes = sorted(claimed_tes.get(key, ()))
            if not tes:
                continue
            strand = key[1]
            for e in ends:
                k = bisect.bisect_left(tes, e)
                cands = [tes[j] for j in (k - 1, k) if 0 <= j < len(tes)]
                d = min((e - t for t in cands), key=abs)
                dists.append(d if strand == "+" else -d)
        model = fit_endpoint_model(dists) if len(dists) >= 2 else fit_endpoint_model(
            endpoint_distances(by_locus, ids_by_locus))
    else:
        model = fit_endpoint_model(endpoint_distances(by_locus, ids_by_locus))
    start_model = None
    if read_starts and claimed_tss:
        # signed 5' distance: u > 0 = the read starts UPSTREAM of the nearest
        # claimed TSS (transcript orientation). Fitted on ALL primary reads (no
        # polyA needed): the 5' end is what the read has, degraded or not.
        ups: List[int] = []
        for rid, (chrom, strand, s5) in read_starts.items():
            if rid in excl:
                continue
            tss = claimed_tss.get((chrom, strand))
            if not tss:
                continue
            tss = sorted(tss)
            k = bisect.bisect_left(tss, s5)
            cands = [tss[j] for j in (k - 1, k) if 0 <= j < len(tss)]
            d = min((s5 - t for t in cands), key=abs)
            u = -d if strand == "+" else d          # upstream is negative genomic on +
            # One-sided: a degraded read (downstream start) is not an upstream
            # observation; it contributes 0 (at/inside the model), so the narrow
            # component is the at-TSS mass and the broad one the upstream errors.
            ups.append(max(u, 0))
        if len(ups) >= 2:
            start_model = fit_start_model(ups)
    return SampleNuisanceModel(
        endpoint=model,
        polya_call_rate=(n_pass / n_scorable) if n_scorable else 0.0,
        n_scorable=n_scorable,
        start=start_model,
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
