"""Unified, evidence-based candidate EXISTENCE test (audit-only v1).

Today a candidate's existence is decided by a cascade of abundance proxies
(``min_abundance``, ``max_soft_mass_ratio``, ``min_fulllen_fraction``, the M2
support gate + its complete-chain rescue). Those constants are six-H9-sample
fits, and the second dataset showed exactly what that means: on Sequin the
SIRV working point halves recall (Sn 16.05 -> 8.85) while the quantification
quality of what it keeps is unchanged (matched Spearman 0.727 vs 0.718). It is
not an existence criterion; it is a cut calibrated to one depth distribution.

This module implements the replacement in the only order that can be believed:
**audit first, decide nothing**. The criterion is

    a candidate exists iff adding it produces a REPRODUCIBLE predictive
    improvement, on READS HELD OUT of the fit, for evidence the other
    candidates in its conflict component cannot explain.

Design + pre-registration: ``CANDIDATE_EXISTENCE_DESIGN.md``. Review that fixed
several errors in the first draft: ``.omc/codex_review_existence_design_2026-09-04.md``.

v1 scope (deliberately narrow, so a calibration failure is attributable):

* archetype ``gtf_multi_exon_m2`` ONLY -- externally pre-specified annotated
  multi-exon candidates. They were not proposed by these reads, so they are
  selection-honest by construction; novel candidates are NOT testable here
  because discovery/collapse used every read (cross-fitting the abundance
  alone does not repair that -- see the design doc, section 5.1).
* structural channel ONLY (M1 raw-AS tie mask + M2 junction NLL contrast,
  eligibility = the read actually covers the discrimination region). The 5'
  hazard and polyA channels are deliberately absent: with two channels at once
  a calibration failure cannot be attributed. (polyA is not even reachable
  here -- no per-read polyA call exists on CandidateSet/QuantOutput.)
* the evidence read is ``QuantOutput.structural_d_tx_raw`` -- the per-read
  junction NLLs over the raw AS tie set -- NOT the ``d_tx`` the EM consumes.
  With the diff-coverage gate on (the default) ``d_tx`` is a DECISION matrix:
  a distinguishing read is collapsed to one-hot and a fuzzy read carries the
  component-wide covered_vote ratio, so scoring it would both discard the
  strongest evidence and feed other reads' votes back in as if they were this
  read's likelihood.
* the test defines its OWN model: a shared-weight finite mixture over the
  structural kernel, fitted per fold. That is deliberately not the production
  per-read softmax (which has no fitted parameter at beta=0, hence nothing to
  cross-fit); H0 refits the remaining weights, because letting the alternatives
  re-absorb the removed candidate's mass is what makes the nested comparison
  honest. A single held-out read can therefore give a NEGATIVE delta -- that is
  the criterion working, not a bug.
* three-state outcome. ``p > alpha`` is NOT evidence of absence, so nothing is
  ever dropped: statuses are ``supported`` / ``unresolved`` (which includes
  every abstention) and ``contradicted``, which v1 never emits because it
  needs its own pre-registered negative-evidence test.

Nothing here mutates the candidate set, the responsibilities, or the outputs.
"""
from __future__ import annotations

import hashlib
import math
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

SCHEMA_VERSION = "candidate-existence-v1"
FOLD_SALT = "candidate-existence-v1"
N_FOLDS = 5
SHARD_DIR = "_existence_shards"
OUT_NAME = "candidate_existence.tsv"

STATUS_SUPPORTED = "supported"
STATUS_UNRESOLVED = "unresolved"
STATUS_CONTRADICTED = "contradicted"      # v1 never emits this

TSV_FIELDS = (
    "sample_id", "interval", "chrom", "strand", "candidate_id", "source",
    "family_id", "archetype", "component_id", "component_size",
    "dominant_alternative", "selection_honest", "identifiable",
    "identifiability_reason", "n_candidate_folds", "n_tested_folds",
    "n_discriminating_reads", "n_covered_reads", "n_finite_contrasts",
    "delta_predictive_ll", "fold_delta_ll", "structural_llr", "tss_llr",
    "tes_polya_llr", "null_method", "bootstrap_replicates",
    "bootstrap_exceedances", "p_value", "p_valid", "q_value",
    "fdr_family_size", "fdr_method", "status", "decision",
    "h1_converged", "h0_converged", "software_version", "fold_salt",
    "bootstrap_seed",
    # 3'/polyA channel audit (E3): per-candidate counts over held-out reads
    "polya_available", "n_polya_pass", "n_tes_near", "n_tes_readthrough",
    "endpoint_kernel_version", "negative_channel",
    # directional bootstrap + hierarchical control (E4)
    "p_minus", "component_p_plus", "component_p_minus",
    "component_q_plus", "component_q_minus", "holm_p_plus", "holm_p_minus",
    "component_fdr_method", "within_component_method", "directional_alpha",
    "n_component_reads",
)

# `d_tx` uses this sentinel for "this read cannot come from this candidate".
# Anything at or above it is treated as no support at all.
_MISSING_FLOOR = 1e5


@dataclass(frozen=True)
class ConflictComponent:
    """Candidates that compete for the same reads' discriminating evidence.

    Edges come from EVIDENCE, never from the fitted posterior: ``a`` and ``b``
    are linked when some read has both in its raw-AS tie set with a finite M2
    NLL for both. Using ``R[i, j] > eps`` instead would make the test family
    depend on the very fit being tested (and on eps).

    Coverage of the discrimination region is deliberately NOT an edge condition.
    The component delimits statistical DEPENDENCE -- which candidates compete for
    mass, hence the scope of the H0 refit -- and two candidates a read supports
    are coupled whether or not that read can tell them apart. Coverage governs
    whether a read may CONTRIBUTE EVIDENCE, and is applied per read at scoring.
    """

    component_id: str
    columns: Tuple[int, ...]              # indices into quant_output.cand_list
    read_rows: Tuple[int, ...]            # rows with finite support in this component


@dataclass
class CandidateExistenceRow:
    sample_id: str = ""
    interval: str = ""
    chrom: str = ""
    strand: str = ""
    candidate_id: str = ""
    source: str = ""
    family_id: str = ""
    archetype: str = ""
    component_id: str = ""
    component_size: int = 0
    dominant_alternative: str = ""
    selection_honest: bool = False
    identifiable: bool = False
    identifiability_reason: str = ""
    n_candidate_folds: int = 0
    n_tested_folds: int = 0
    n_discriminating_reads: int = 0
    n_covered_reads: int = 0
    n_finite_contrasts: int = 0
    delta_predictive_ll: Optional[float] = None
    fold_delta_ll: str = ""               # ";"-joined per-fold values
    structural_llr: Optional[float] = None
    tss_llr: Optional[float] = None       # v1: left EMPTY, never 0
    tes_polya_llr: Optional[float] = None
    null_method: str = ""
    bootstrap_replicates: int = 0
    bootstrap_exceedances: int = 0
    p_value: Optional[float] = None
    p_valid: bool = False
    q_value: Optional[float] = None
    fdr_family_size: int = 0
    fdr_method: str = ""
    status: str = STATUS_UNRESOLVED
    decision: str = "no_action"           # audit: always
    h1_converged: bool = True
    h0_converged: bool = True
    software_version: str = SCHEMA_VERSION
    fold_salt: str = FOLD_SALT
    bootstrap_seed: int = 0
    polya_available: bool = False
    n_polya_pass: int = 0
    n_tes_near: int = 0
    n_tes_readthrough: int = 0
    endpoint_kernel_version: str = ""
    negative_channel: str = ""
    p_minus: Optional[float] = None
    component_p_plus: Optional[float] = None
    component_p_minus: Optional[float] = None
    component_q_plus: Optional[float] = None
    component_q_minus: Optional[float] = None
    holm_p_plus: Optional[float] = None
    holm_p_minus: Optional[float] = None
    component_fdr_method: str = ""
    within_component_method: str = ""
    directional_alpha: float = 0.025
    n_component_reads: int = 0

    def as_tsv(self) -> str:
        out = []
        for key in TSV_FIELDS:
            v = getattr(self, key)
            if v is None:
                out.append("")
            elif isinstance(v, bool):
                out.append("true" if v else "false")
            elif isinstance(v, float):
                out.append(f"{v:.6g}")
            else:
                out.append(str(v))
        return "\t".join(out)


def read_fold(read_id: str, *, n_folds: int = N_FOLDS, salt: str = FOLD_SALT) -> int:
    """Deterministic fold of a read, stable across runs, machines and orderings."""
    digest = hashlib.blake2b(f"{salt}|{read_id}".encode("utf-8"), digest_size=8)
    return int.from_bytes(digest.digest(), "big") % n_folds


def component_folds(component_id: str, read_ids: Sequence[str], *,
                    max_folds: int = N_FOLDS, salt: str = FOLD_SALT):
    """Balanced, deterministic folds for ONE component's reads.

    A global hash-mod-K leaves small components with empty folds -- at Sequin's
    depth (a few reads per component) most folds have no test read at all, so
    the candidate is discarded for a purely technical reason. Here K adapts,
    ``K = min(max_folds, n_reads)``, and reads are dealt round-robin in a
    hash-sorted order, so every fold is non-empty and every read is scored
    exactly once while being absent from its own training set. At n <= max_folds
    this degenerates to leave-one-out, which is the correct limit.

    One fold map is shared by ALL candidates of the component: optimising the
    split per candidate would be a second selection effect.

    Returns ``(fold_by_read_id, K)``; K == 0 when a component has < 2 reads
    (no honest split exists, so the candidate must abstain).
    """
    ids = sorted(set(read_ids))
    if len(ids) < 2:
        return {}, 0
    k = min(max_folds, len(ids))
    ordered = sorted(ids, key=lambda r: hashlib.blake2b(
        f"{salt}|{component_id}|{r}".encode("utf-8"), digest_size=8).digest())
    return {r: i % k for i, r in enumerate(ordered)}, k


def _component_id(candidate_ids: Sequence[str]) -> str:
    payload = "|".join(sorted(candidate_ids)).encode("utf-8")
    return "cmp_" + hashlib.blake2b(payload, digest_size=6).hexdigest()


def build_conflict_components(
    quant_output,
    candidate_set,
    *,
    require_covered_contrast: bool = False,
) -> List[ConflictComponent]:
    """Group candidates into components joined by shared discriminating evidence.

    ``require_covered_contrast`` is retained only as an experiment switch; the
    criterion wants it False.
    """
    cand_list = list(getattr(quant_output, "cand_list", []) or [])
    n_c = len(cand_list)
    if n_c == 0:
        return []
    ties = getattr(quant_output, "ties_by_read", {}) or {}
    nlls = getattr(quant_output, "nlls_by_read", {}) or {}
    cover = getattr(quant_output, "cover_by_read", {}) or {}

    parent = list(range(n_c))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    rows_by_col: Dict[int, set] = {}       # every read supporting a candidate
    edge_cols: set = set()                 # candidates that have a scored contrast
    for i in sorted(ties):
        tie = ties[i]
        for j in tie:
            rows_by_col.setdefault(j, set()).add(i)
        if require_covered_contrast and cover.get(i) is not True:
            continue
        scored = [j for j in tie if j in (nlls.get(i) or {})]
        if len(scored) < 2:
            continue
        edge_cols.update(scored)
        first = scored[0]
        for j in scored[1:]:
            a, b = find(first), find(j)
            if a != b:
                parent[max(a, b)] = min(a, b)

    groups: Dict[int, List[int]] = {}
    for j in sorted(edge_cols):
        groups.setdefault(find(j), []).append(j)

    out: List[ConflictComponent] = []
    for cols in groups.values():
        if len(cols) < 2:
            continue                       # nothing to contrast against
        cols_sorted = tuple(sorted(
            cols, key=lambda j: getattr(cand_list[j], "candidate_id", str(j))))
        rows: set = set()
        for j in cols_sorted:
            rows |= rows_by_col.get(j, set())
        out.append(ConflictComponent(
            component_id=_component_id(
                [getattr(cand_list[j], "candidate_id", str(j)) for j in cols_sorted]),
            columns=cols_sorted,
            read_rows=tuple(sorted(rows)),
        ))
    out.sort(key=lambda c: c.component_id)
    return out


def _log_weights(d_row, cols) -> Dict[int, float]:
    """Unnormalised log-likelihood weights of one read over ``cols``.

    ``d_tx`` holds distances (0 = best) consumed by the EM as ``exp(-d)`` at
    sigma=1, so the log weight is simply ``-d``. Cells at the MISSING sentinel
    are absent, not zero-weighted: absence of evidence is abstention.
    """
    out = {}
    for j in cols:
        d = float(d_row[j])
        if d < _MISSING_FLOOR and math.isfinite(d):
            out[j] = -d
    return out


def _predictive_ll(weights: Dict[int, float], theta: Dict[int, float]) -> Optional[float]:
    """log sum_j theta_j * exp(logw_j) over the candidates the read can come from."""
    terms = [(theta[j], w) for j, w in weights.items() if theta.get(j, 0.0) > 0.0]
    if not terms:
        return None
    m = max(w for _, w in terms)
    s = sum(t * math.exp(w - m) for t, w in terms)
    if s <= 0.0:
        return None
    return m + math.log(s)


def _fit_theta(rows, d_tx, cols, *, max_iter: int = 200, tol: float = 1e-9):
    """EM for mixture weights over ``cols`` using ``rows`` (exact for the
    production default: beta=0 and no abundance feedback, so the model is a
    plain finite mixture over the fixed per-cell weights)."""
    theta = {j: 1.0 / len(cols) for j in cols}
    if not rows:
        return theta, True
    cached = [_log_weights(d_tx[i], cols) for i in rows]
    cached = [w for w in cached if w]
    if not cached:
        return theta, True
    converged = False
    for _ in range(max_iter):
        acc = {j: 0.0 for j in cols}
        n_used = 0
        for w in cached:
            terms = {j: theta[j] * math.exp(lw - max(w.values()))
                     for j, lw in w.items() if theta.get(j, 0.0) > 0.0}
            z = sum(terms.values())
            if z <= 0.0:
                continue
            n_used += 1
            for j, t in terms.items():
                acc[j] += t / z
        if n_used == 0:
            return theta, True
        new = {j: acc[j] / n_used for j in cols}
        delta = max(abs(new[j] - theta[j]) for j in cols)
        theta = new
        if delta < tol:
            converged = True
            break
    return theta, converged


def _fit_theta_w(rows, weights_of, cols, *, max_iter: int = 200, tol: float = 1e-9):
    """As :func:`_fit_theta`, but over an arbitrary per-read log-weight function
    (so channels can be composed before the mixture is fitted)."""
    theta = {j: 1.0 / len(cols) for j in cols}
    if not rows:
        return theta, True
    cached = [weights_of(i) for i in rows]
    cached = [w for w in cached if w]
    if not cached:
        return theta, True
    converged = False
    for _ in range(max_iter):
        acc = {j: 0.0 for j in cols}
        n_used = 0
        for w in cached:
            m = max(w.values())
            terms = {j: theta[j] * math.exp(lw - m)
                     for j, lw in w.items() if theta.get(j, 0.0) > 0.0}
            z = sum(terms.values())
            if z <= 0.0:
                continue
            n_used += 1
            for j, t in terms.items():
                acc[j] += t / z
        if n_used == 0:
            return theta, True
        new = {j: acc[j] / n_used for j in cols}
        delta = max(abs(new[j] - theta[j]) for j in cols)
        theta = new
        if delta < tol:
            converged = True
            break
    return theta, converged


def evaluate_candidate_existence(
    quant_output,
    candidate_set,
    components: Sequence[ConflictComponent],
    interval,
    *,
    sample_id: str = "",
    archetype: str = "gtf_multi_exon_m2",
    n_folds: int = N_FOLDS,
    fold_salt: str = FOLD_SALT,
    bootstrap_replicates: int = 4095,
    bootstrap_seed: int = 0,
    molecules=None,
    nuisance=None,
) -> List[CandidateExistenceRow]:
    """Cross-fitted held-out log-likelihood contrast for every testable candidate.

    Channels. The structural kernel (raw junction NLL over the raw AS tie) is
    always on. When ``molecules`` (per-read :class:`MoleculeEvidence`) and a
    ``nuisance`` model with an available endpoint fit are supplied, the 3'/polyA
    term of Part A is ADDED to each polyA-PASS read's per-candidate log weight:
    ``l_rc = l_struct + l_3'``. Reads without a PASS call contribute nothing on
    that channel (abstention). The statistic is computed twice, structural-only
    and full, so the artifact reports the channel's own contribution
    (``tes_polya_llr = T_full - T_structural``) and a calibration failure can be
    attributed to one channel.

    For each candidate C in each component: fit the mixture WITH and WITHOUT C on
    the reads of the other folds, then score only this fold's reads. Summing over
    folds gives ``T_C``; every read is held out exactly once. Candidates that
    cannot be distinguished from their alternatives on held-out reads are marked
    ``unidentifiable`` and RETAINED -- insufficient evidence is not evidence of
    absence.
    """
    rows_out: List[CandidateExistenceRow] = []
    cand_list = list(getattr(quant_output, "cand_list", []) or [])
    # EVIDENCE, not the EM's decision matrix (see the module docstring).
    kernel = getattr(quant_output, "structural_d_tx_raw", None)
    if kernel is None:
        kernel = getattr(quant_output, "d_tx", None)
    kept = list(getattr(quant_output, "kept_read_ids", []) or [])
    if kernel is None or not cand_list or not components:
        return rows_out
    ivl = (f"{interval.chrom}:{interval.start}-{interval.end}:{interval.strand}"
           if interval is not None else "")
    cover = getattr(quant_output, "cover_by_read", {}) or {}

    endpoint = getattr(nuisance, "endpoint", None) if nuisance is not None else None
    tes_on = bool(molecules) and endpoint is not None and endpoint.available
    if tes_on:
        from fin.analysis.sample_nuisance import tes_distance, tes_polya_loglik

    def weights_for(i: int, cols_, with_tes: bool) -> Dict[int, float]:
        """Per-read log weight over ``cols_`` = structural kernel (+ 3' term)."""
        w = _log_weights(kernel[i], cols_)
        if not with_tes or not w:
            return w
        mol = molecules.get(kept[i]) if i < len(kept) else None
        if mol is None or not mol.polya_available:
            return w                                   # abstain on this channel
        out = {}
        for j, lw in w.items():
            t = tes_polya_loglik(mol, cand_list[j], endpoint)
            out[j] = lw + (t if t is not None else 0.0)
        return out

    def crossfit(comp, cols_, j, others, folds, k_folds, with_tes):
        """Returns (per_fold_deltas, n_tested, n_disc, n_covered, ok1, ok0)."""
        per_fold, n_tested, n_disc, covered = [], 0, 0, 0
        ok1_all = ok0_all = True
        for f in range(k_folds):
            test_rows = [i for i in comp.read_rows if folds.get(i) == f]
            train_rows = [i for i in comp.read_rows if folds.get(i) not in (f, None)]
            if not test_rows or not train_rows:
                continue
            theta1, ok1 = _fit_theta_w(train_rows, lambda i: weights_for(i, cols_, with_tes), cols_)
            theta0, ok0 = _fit_theta_w(train_rows, lambda i: weights_for(i, others, with_tes), others)
            ok1_all = ok1_all and ok1
            ok0_all = ok0_all and ok0
            delta, used = 0.0, 0
            for i in test_rows:
                w_all = weights_for(i, cols_, with_tes)
                if j not in w_all or len(w_all) < 2:
                    continue          # cannot tell C from the alternatives
                # Coverage is a SCORING condition for the structural channel;
                # a polyA read can still contribute its 3' term uncovered.
                mol = molecules.get(kept[i]) if (with_tes and i < len(kept)) else None
                polya_read = bool(mol is not None and mol.polya_available)
                if cover and cover.get(i) is not True and not polya_read:
                    continue
                covered += 1
                l1 = _predictive_ll(w_all, theta1)
                l0 = _predictive_ll({k: v for k, v in w_all.items() if k != j}, theta0)
                if l1 is None or l0 is None:
                    continue
                delta += (l1 - l0)
                used += 1
            if used:
                per_fold.append(delta)
                n_tested += 1
                n_disc += used
        return per_fold, n_tested, n_disc, covered, ok1_all, ok0_all

    for comp in components:
        cols = list(comp.columns)
        # ONE fold map per component, shared by all its candidates: optimising
        # the split per candidate would be a second selection effect.
        comp_read_ids = [kept[i] for i in comp.read_rows if i < len(kept)]
        fold_of_id, k_folds = component_folds(
            comp.component_id, comp_read_ids, max_folds=n_folds, salt=fold_salt)
        folds = {i: fold_of_id[kept[i]] for i in comp.read_rows
                 if i < len(kept) and kept[i] in fold_of_id}
        for j in cols:
            cand = cand_list[j]
            row = CandidateExistenceRow(
                sample_id=sample_id, interval=ivl,
                chrom=getattr(cand, "chrom", ""), strand=getattr(cand, "strand", ""),
                candidate_id=getattr(cand, "candidate_id", str(j)),
                source=getattr(cand, "source", ""),
                family_id=getattr(cand, "family_id", "") or "",
                archetype=archetype, component_id=comp.component_id,
                component_size=len(cols), n_candidate_folds=k_folds,
                bootstrap_replicates=0, bootstrap_seed=bootstrap_seed,
                null_method="none_v1_audit", fold_salt=fold_salt,
                polya_available=tes_on,
                endpoint_kernel_version=(endpoint.version if tes_on else ""),
            )
            rows_out.append(row)

            row.selection_honest = getattr(cand, "source", "") == "gtf"
            if not row.selection_honest:
                row.identifiability_reason = "not_selection_honest_v1"
                continue
            if not getattr(cand, "intron_chain", None) or not cand.intron_chain.introns:
                row.identifiability_reason = "not_multi_exon"
                continue
            others = [k for k in cols if k != j]
            if not others:
                row.identifiability_reason = "no_alternative"
                continue
            if k_folds < 2:
                row.identifiability_reason = "no_honest_split"
                continue

            # structural-only statistic (the v1 baseline, kept for attribution)
            pf_s, nt_s, nd_s, cov_s, ok1, ok0 = crossfit(comp, cols, j, others, folds, k_folds, False)
            row.h1_converged, row.h0_converged = ok1, ok0
            t_struct = 2.0 * sum(pf_s) if pf_s else None
            if tes_on:
                pf_f, nt_f, nd_f, cov_f, ok1f, ok0f = crossfit(comp, cols, j, others, folds, k_folds, True)
                row.h1_converged = row.h1_converged and ok1f
                row.h0_converged = row.h0_converged and ok0f
                per_fold, n_tested, n_disc, covered = pf_f, nt_f, nd_f, cov_f
                # per-candidate 3' audit counts over the component's reads
                for i in comp.read_rows:
                    mol = molecules.get(kept[i]) if i < len(kept) else None
                    if mol is None or not mol.polya_available:
                        continue
                    if j not in _log_weights(kernel[i], cols):
                        continue
                    row.n_polya_pass += 1
                    d = tes_distance(mol, cand)
                    post = endpoint.narrow_posterior(d)
                    if post >= 0.5:
                        row.n_tes_near += 1
                    elif d > 0:
                        row.n_tes_readthrough += 1
                if row.n_tes_readthrough:
                    row.negative_channel = "tes_readthrough_polya"
            else:
                per_fold, n_tested, n_disc, covered = pf_s, nt_s, nd_s, cov_s
            row.n_tested_folds = n_tested
            row.n_discriminating_reads = n_disc
            row.n_covered_reads = covered
            row.n_finite_contrasts = n_disc
            if not per_fold:
                row.identifiability_reason = "no_held_out_discriminating_read"
                continue
            row.identifiable = True
            row.identifiability_reason = "ok"
            row.delta_predictive_ll = 2.0 * sum(per_fold)
            row.structural_llr = t_struct
            if tes_on:
                row.tes_polya_llr = (row.delta_predictive_ll - (t_struct or 0.0))
            row.fold_delta_ll = ";".join(f"{2.0 * v:.6g}" for v in per_fold)
            # No calibrated null is fitted yet, so there is deliberately NO
            # p-value: reporting one would be the exact overreach this module
            # exists to prevent.
            row.p_valid = False
            row.status = STATUS_UNRESOLVED
            theta_all, _ = _fit_theta_w(list(comp.read_rows),
                                        lambda i: weights_for(i, cols, tes_on), cols)
            alt = max(others, key=lambda k: (theta_all.get(k, 0.0),
                                             getattr(cand_list[k], "candidate_id", "")))
            row.dominant_alternative = getattr(cand_list[alt], "candidate_id", str(alt))

        # --- calibrated null: centred molecule-level bootstrap over the component
        if bootstrap_replicates > 0:
            _attach_bootstrap(comp, cols, rows_out, kept, cover, molecules, tes_on,
                              weights_for, n_folds, bootstrap_replicates,
                              bootstrap_seed, cand_list)
    return rows_out


def _attach_bootstrap(comp, cols, rows_out, kept, cover, molecules, tes_on,
                      weights_for, n_folds, replicates, seed, cand_list) -> None:
    """Build the component's (n x k) log-weight matrix and eligibility mask from
    the SAME per-read channel composition the observed statistic used, run the
    bootstrap for every testable column, and write p+/p- onto the rows.

    A row is testable iff it is identifiable (the audit has already decided
    that on the observed data). Rows that are not testable never enter the
    null -- they are not part of any family."""
    import numpy as np

    from fin.analysis.existence_bootstrap import NULL_METHOD, bootstrap_component

    rows_by_cid = {r.candidate_id: r for r in rows_out
                   if r.component_id == comp.component_id}
    testable = [j for j in cols
                if rows_by_cid.get(getattr(cand_list[j], "candidate_id", str(j)))
                and rows_by_cid[getattr(cand_list[j], "candidate_id", str(j))].identifiable]
    n = len(comp.read_rows)
    for r in rows_by_cid.values():
        r.n_component_reads = n
    if not testable or n < 2:
        return
    W = np.full((n, len(cols)), -np.inf)
    elig = np.zeros(n, dtype=bool)
    col_pos = {j: p for p, j in enumerate(cols)}
    for a, i in enumerate(comp.read_rows):
        w = weights_for(i, cols, tes_on)
        for j, lw in w.items():
            W[a, col_pos[j]] = lw
        mol = molecules.get(kept[i]) if (tes_on and molecules and i < len(kept)) else None
        polya_read = bool(mol is not None and mol.polya_available)
        elig[a] = (not cover) or cover.get(i) is True or polya_read
    res = bootstrap_component(
        W, elig, [col_pos[j] for j in testable], max_folds=n_folds,
        replicates=replicates, seed=seed, component_tag=comp.component_id)
    for j in testable:
        r = rows_by_cid[getattr(cand_list[j], "candidate_id", str(j))]
        pj = col_pos[j]
        r.null_method = NULL_METHOD
        r.bootstrap_replicates = res.replicates
        r.bootstrap_exceedances = res.exceed_plus.get(pj, 0)
        # the bootstrap re-derives T on its own (random balanced) folds; keep the
        # deterministic hash-fold value as the reported statistic and record the
        # bootstrap's own observed value for auditability via fold_delta_ll
        if res.degenerate or not res.replicates:
            r.p_valid = False
            r.identifiability_reason = "bootstrap_degenerate"
            continue
        r.p_value = res.p_plus[pj]
        r.p_minus = res.p_minus[pj]
        r.component_p_plus = res.component_p_plus
        r.component_p_minus = res.component_p_minus
        r.p_valid = True


def apply_sample_fdr(
    rows: Sequence[CandidateExistenceRow],
    *,
    alpha: float = 0.05,
    method: str = "by",
) -> List[CandidateExistenceRow]:
    """Sample-level, directional, hierarchical error control (spec B6).

    Family = rows with a VALID directional p (identifiable, selection-honest,
    bootstrap not degenerate). ``alpha`` is split over the two directions
    (0.025 each). Per direction:

      layer 1  component omnibus p (max statistic, same replicates) ->
               Benjamini-Yekutieli across the sample's components;
      layer 2  within a component whose q passes, Holm over its candidates.

    supported     : component q+ <= a_d AND Holm p+ <= a_d AND T > 0
    contradicted  : component q- <= a_d AND Holm p- <= a_d AND T < 0
                    AND a pre-registered negative channel is active
    unresolved    : everything else -- and it is RETAINED.

    ``fdr_family_size`` is the number of components actually corrected.
    Never per chromosome: chrIS alone would silently change the burden.
    """
    from fin.analysis.existence_bootstrap import by_adjust, holm_adjust

    a_d = alpha / 2.0
    valid = [r for r in rows if r.p_valid and r.identifiable and r.selection_honest
             and r.p_value is not None and r.p_minus is not None]
    comps: Dict[str, List[CandidateExistenceRow]] = {}
    for r in valid:
        comps.setdefault(r.component_id, []).append(r)
    for r in rows:
        r.fdr_method = method
        r.component_fdr_method = "BY"
        r.within_component_method = "Holm"
        r.directional_alpha = a_d
        r.fdr_family_size = len(comps)
        r.status = STATUS_UNRESOLVED
    if not comps:
        return list(rows)
    cids = sorted(comps)
    for direction in ("plus", "minus"):
        comp_p = [min((getattr(r, f"component_p_{direction}") or 1.0) for r in comps[c])
                  for c in cids]
        comp_q = by_adjust(comp_p)
        for c, q in zip(cids, comp_q):
            members = comps[c]
            pk = "p_value" if direction == "plus" else "p_minus"
            holm = holm_adjust([getattr(r, pk) for r in members])
            for r, hp in zip(members, holm):
                setattr(r, f"component_q_{direction}", q)
                setattr(r, f"holm_p_{direction}", hp)
    for r in valid:
        t = r.delta_predictive_ll or 0.0
        if (r.component_q_plus is not None and r.component_q_plus <= a_d
                and r.holm_p_plus is not None and r.holm_p_plus <= a_d and t > 0):
            r.status = STATUS_SUPPORTED
            r.q_value = r.component_q_plus
        elif (r.component_q_minus is not None and r.component_q_minus <= a_d
                and r.holm_p_minus is not None and r.holm_p_minus <= a_d and t < 0
                and r.negative_channel):
            r.status = STATUS_CONTRADICTED
            r.q_value = r.component_q_minus
        else:
            r.q_value = r.component_q_plus
    return list(rows)


def write_existence_shard(work_dir, interval, rows: Sequence[CandidateExistenceRow]) -> None:
    if not rows:
        return
    d = Path(work_dir) / SHARD_DIR
    d.mkdir(parents=True, exist_ok=True)
    name = (f"{interval.chrom}_{interval.start}_{interval.end}_"
            f"{interval.strand or 'u'}.tsv")
    tmp = d / (name + ".tmp")
    with tmp.open("w") as fh:
        for r in rows:
            fh.write(r.as_tsv() + "\n")
    tmp.replace(d / name)


def merge_existence_shards(work_dir, *, alpha: float = 0.05) -> Optional[str]:
    """Merge interval shards, apply the SAMPLE-level FDR, write the artifact."""
    wd = Path(work_dir)
    shard_dir = wd / SHARD_DIR
    raw: List[str] = []
    if shard_dir.is_dir():
        for shard in sorted(shard_dir.glob("*.tsv")):
            with shard.open() as fh:
                raw.extend(line.rstrip("\n") for line in fh if line.strip())
    if not raw:
        if shard_dir.is_dir():
            shutil.rmtree(shard_dir)
        return None
    rows = []
    for line in raw:
        parts = line.split("\t")
        row = CandidateExistenceRow()
        for key, val in zip(TSV_FIELDS, parts):
            cur = getattr(row, key)
            if isinstance(cur, bool):
                setattr(row, key, val == "true")
            elif isinstance(cur, int) and not isinstance(cur, bool):
                setattr(row, key, int(val) if val else 0)
            elif key in ("delta_predictive_ll", "structural_llr", "tss_llr",
                         "tes_polya_llr", "p_value", "q_value", "p_minus",
                         "component_p_plus", "component_p_minus",
                         "component_q_plus", "component_q_minus",
                         "holm_p_plus", "holm_p_minus"):
                setattr(row, key, float(val) if val else None)
            elif key == "directional_alpha":
                setattr(row, key, float(val) if val else 0.025)
            else:
                setattr(row, key, val)
        rows.append(row)
    rows = apply_sample_fdr(rows, alpha=alpha)
    rows.sort(key=lambda r: (r.chrom, r.interval, r.component_id, r.candidate_id))
    out = wd / OUT_NAME
    tmp = out.with_suffix(".tsv.tmp")
    with tmp.open("w") as fh:
        fh.write("\t".join(TSV_FIELDS) + "\n")
        for r in rows:
            fh.write(r.as_tsv() + "\n")
    tmp.replace(out)
    shutil.rmtree(shard_dir)
    return str(out)
