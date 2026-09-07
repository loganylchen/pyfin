"""The per-molecule evidence row ``x_r`` of the generative model.

Every decision the pipeline makes -- does a transcript model exist, how
abundant is it, which model did this read come from -- is a statement about
the same likelihood ``P(x_r | c, eta_s)``. Part A of the specification
(``.omc/codex_spec_full_algorithm_2026-09-05.md``) factorises it over four
observation channels of one read ``r``:

    x_r = (A_r, C_r, J_r, S_r, E_r, Q_r, L_r)

    A_r  read-to-candidate alignment score / raw best-AS tie mask
    C_r  primary CIGAR: junctions, local indels, span
    J_r  junction-signal NLL over contested windows (M2)
    S_r  strand-aware genomic 5' endpoint
    E_r  strand-aware genomic 3' endpoint
    Q_r  polyA call state: PASS / FAIL / MISSING
    L_r  polyA length (only to form a reliable PASS; never transcript evidence)

This module defines that row as an IMMUTABLE value and builds it from the
carriers the pipeline already produces. Two properties matter and are pinned
by tests:

* one row per molecule. Secondary and supplementary alignments are placements
  of the same query, not independent evidence, and never produce a second row.
* a channel that was not observed is ABSENT (``None``), never a default that
  looks like evidence. A missing polyA call, a read the signal stack did not
  score, a candidate outside the read's raw tie set -- each contributes 0 to
  every likelihood ratio, which is abstention. Filling in 0 NLL or a mean
  would manufacture support; filling in a penalty would manufacture refutation.

Nothing here scores anything. Channel likelihoods live next to their nuisance
models; this file is the data contract they share, and it is the resampling
unit of the existence test's bootstrap (a molecule is resampled as a whole,
never one channel of it).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Mapping, Optional, Sequence, Tuple

POLYA_PASS = "PASS"
POLYA_FAIL = "FAIL"
POLYA_MISSING = "MISSING"
# krill's own verdicts, before they are collapsed to the three states above.
# ADAPTER means a polyA region WAS located but the length estimate is
# contaminated by the adapter signal; the 3' channel does not use the length.
POLYA_QC_ADAPTER = "ADAPTER"
POLYA_QC_NOREGION = "NOREGION"


def polya_end_evidential(
    call: Optional[Tuple[Optional[float], Optional[str]]],
    accept_adapter: bool = False,
    qc_conditional: bool = False,
) -> bool:
    """Is this call admitted to the 3' endpoint channel?

    The 3'/polyA channel never reads the tail LENGTH -- it scores the read's
    genomic 3' end, which comes from the alignment -- so what it needs from
    krill is the existence statement: a polyA region was found on this molecule,
    therefore this end is a polyadenylation site and not a fragment. There are
    two ways to decide which verdicts make that statement, and this function
    implements both.

    BY VERDICT NAME (the default, and ``accept_adapter``).
    ``PASS`` says exactly that. ``ADAPTER`` says the region was found but the
    LENGTH is contaminated, which the channel does not care about, and the ends
    bear that out: reads within 20 bp of an annotated TES are 82.8% of PASS
    against 79.9% of ADAPTER on human dRNA (95.3 / 92.1 on HEYA8, 95.6 / 89.0
    on Sequin), while ``NOREGION`` -- no polyA region at all -- is 10.1% on
    human and lands a median 334 bp away (it is NOT bad everywhere; see below). ``READ_FAILED`` is krill declining to
    process the read: technical missingness, which abstains by the same rule
    that makes a missing call abstain rather than count against a model.

    Takes krill's ``(polya_length, polya_qc)`` so the PASS branch stays exactly
    :func:`polya_state`'s PASS -- a PASS verdict with no usable length is still
    not evidence -- and only the ADAPTER branch is new. An ADAPTER call has no
    trustworthy length by construction, so it is judged on the verdict alone.

    BY FITTED KERNEL (``qc_conditional``), which applies to the ASSIGNMENT
    channel only -- the discovery contradiction hook keeps the verdict-name rule
    above, because its test is an absolute posterior threshold whose meaning
    changes with the fitted eta (see `runner._read_extent_contradicted`). Every
    read krill returned a verdict for is admitted HERE, and the abstention
    decision moves to the kernel: each
    class is fitted separately and a class whose fit is not `available` scores
    nothing (see :meth:`SampleNuisanceModel.endpoint_for`). So under this mode
    the paragraph above is NOT the rule -- a NOREGION read is admitted by this
    function and abstains only if its class failed to fit. That is deliberate:
    the classes are not equally good on every run (NOREGION is 10.1% within
    20 bp of an annotated TES on the human sample and 74-75% on HEYA8), so which
    ones to trust is a per-sample estimate rather than a name on a list.
    """
    if qc_conditional:
        return call is not None and bool(call[1])
    if polya_state(call)[0] == POLYA_PASS:
        return True
    if not accept_adapter or call is None:
        return False
    return call[1] == POLYA_QC_ADAPTER


@dataclass(frozen=True)
class MoleculeEvidence:
    """One read's complete observation row. Frozen: consumers may not edit it."""

    read_id: str
    strand: str                                   # "+" / "-" as mapped
    # A_r: candidate_id -> raw mappy AS; raw_tie_ids = the best-AS tie set
    alignment_scores: Mapping[str, float]
    raw_tie_ids: Tuple[str, ...]
    # C_r
    cigar_junctions: Tuple[Tuple[int, int], ...]  # (donor, acceptor) genomic
    indel_spans: Tuple[Tuple[int, int], ...]      # D as [lo,hi), I as (p,p)
    reference_start: int
    reference_end: int
    # S_r / E_r: strand-aware. five_prime_pos is the read's transcriptional
    # 5' end in genomic coordinates (reference_start on +, reference_end on -).
    five_prime_pos: int
    three_prime_pos: int
    # J_r: candidate_id -> contested-window junction NLL (only scored ties)
    m2_nll: Mapping[str, float]
    # candidate_id -> covered discrimination-region tags; empty = no coverage
    covered_diff_regions: Mapping[str, Tuple[str, ...]]
    # Q_r / L_r
    polya_qc: str                                 # PASS / FAIL / MISSING
    polya_length: Optional[float]
    mapping_quality: int = 0

    # -- channel availability (abstention semantics) -------------------------
    def structural_available(self, cid: str) -> bool:
        return cid in self.raw_tie_ids

    def signal_available(self, cid: str) -> bool:
        return cid in self.m2_nll

    @property
    def polya_available(self) -> bool:
        """Only a PASS call carries 3' evidence. FAIL and MISSING both abstain:
        a failed call is a statement about the signal, not about the transcript."""
        return self.polya_qc == POLYA_PASS


def polya_state(call: Optional[Tuple[Optional[float], Optional[str]]]) -> Tuple[str, Optional[float]]:
    """Map krill's ``(polya_length, polya_qc)`` (or None) to (state, length)."""
    if call is None:
        return POLYA_MISSING, None
    length, qc = call
    if qc == POLYA_PASS and length is not None and length > 0:
        return POLYA_PASS, float(length)
    if qc is None:
        return POLYA_MISSING, None
    return POLYA_FAIL, (float(length) if length is not None else None)


def build_molecule_evidence(
    quant_output,
    candidate_set,
    *,
    polya_map: Optional[Mapping[str, Tuple[Optional[float], Optional[str]]]] = None,
    read_strands: Optional[Mapping[str, str]] = None,
) -> Dict[str, MoleculeEvidence]:
    """Assemble one frozen row per kept read from the pipeline's carriers.

    Sources (all already produced on the default path):
      * ``quant_output.cand_list / raw / ties_by_read / nlls_by_read /
        cover_by_read / kept_read_ids``
      * ``candidate_set.read_spans / indels_by_read``
      * per-read CIGAR junctions from ``candidate_set.read_chains`` when the
        discovery path recorded them, else empty (abstains structurally)
      * ``polya_map`` from :func:`fin.scoring.polya.compute_polya`; None means
        the signal stack did not run -> every row is MISSING.

    Reads absent from ``kept_read_ids`` (no candidate at AS>0) are not rows:
    they carry no candidate-comparable evidence.
    """
    cand_ids = [getattr(c, "candidate_id", str(j))
                for j, c in enumerate(getattr(quant_output, "cand_list", []) or [])]
    kept = list(getattr(quant_output, "kept_read_ids", []) or [])
    raw = getattr(quant_output, "raw", None)
    ties = getattr(quant_output, "ties_by_read", {}) or {}
    nlls = getattr(quant_output, "nlls_by_read", {}) or {}
    cover = getattr(quant_output, "cover_by_read", {}) or {}
    spans = getattr(candidate_set, "read_spans", None) or {}
    indels = getattr(candidate_set, "indels_by_read", None) or {}
    chains = getattr(candidate_set, "read_junctions", None) or {}
    strands = dict(getattr(candidate_set, "read_strands", None) or {})
    strands.update(read_strands or {})
    polya_map = polya_map or {}

    out: Dict[str, MoleculeEvidence] = {}
    for i, rid in enumerate(kept):
        if rid in out:
            continue                                # one row per molecule
        span = spans.get(rid)
        if span is None:
            continue
        rs, re = int(span[0]), int(span[1])
        strand = strands.get(rid) or getattr(candidate_set.interval, "strand", None) or "+"
        five, three = (rs, re) if strand == "+" else (re, rs)
        scores: Dict[str, float] = {}
        if raw is not None and i < raw.shape[0]:
            for j, cid in enumerate(cand_ids):
                v = float(raw[i, j])
                if v > float("-inf"):
                    scores[cid] = v
        tie_ids = tuple(cand_ids[j] for j in (ties.get(i) or []) if j < len(cand_ids))
        m2 = {cand_ids[j]: float(v) for j, v in (nlls.get(i) or {}).items()
              if j < len(cand_ids)}
        covered: Dict[str, Tuple[str, ...]] = {}
        if cover.get(i) is True:
            for cid in tie_ids:
                covered[cid] = ("diff",)
        qc, length = polya_state(polya_map.get(rid))
        out[rid] = MoleculeEvidence(
            read_id=rid, strand=strand,
            alignment_scores=dict(scores), raw_tie_ids=tie_ids,
            cigar_junctions=tuple(chains.get(rid, ())),
            indel_spans=tuple(indels.get(rid, ())),
            reference_start=rs, reference_end=re,
            five_prime_pos=five, three_prime_pos=three,
            m2_nll=m2, covered_diff_regions=covered,
            polya_qc=qc, polya_length=length,
        )
    return out


def resample_molecules(rows: Sequence[MoleculeEvidence], indices: Sequence[int],
                       replicate: int) -> Tuple[MoleculeEvidence, ...]:
    """Bootstrap copies of whole rows with stable replicate-local ids.

    The resampling unit is the molecule: its channels stay together. A duplicate
    keeps every field and gets the id ``"{id}#boot{replicate}#{copy}"`` so that
    downstream fold assignment (which hashes the id) is deterministic and two
    copies of one read can land in different folds, as independent draws would.
    """
    seen: Dict[str, int] = {}
    out = []
    for k in indices:
        base = rows[k]
        n = seen.get(base.read_id, 0)
        seen[base.read_id] = n + 1
        out.append(MoleculeEvidence(
            **{**base.__dict__, "read_id": f"{base.read_id}#boot{replicate}#{n}"}))
    return tuple(out)
