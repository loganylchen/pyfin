"""Evidence ranker over inference-time observable candidate features.

A frozen L2-logistic model scores every post-selection survivor from the
``candidate_evidence`` feature row with a RAW LOGIT score (ranking evidence,
NOT a calibrated probability - no reliability analysis has been performed and
the existing ``confidence`` field is never overwritten);
``ranking_mode="filter"`` then removes NOVEL candidates scoring below the
frozen operating point BEFORE junction snapping and the final-survivor
abundance refit, so released read mass is re-dealt by the ordinary refit
accounting.

Provenance and discipline:

* Trained only on the H9 replicate2_run2 tuning sample (frozen offline by
  ``experiments/prod_validation/fit_candidate_ranker.py``; grouped
  cross-validation by chromosome, mean AUC 0.809). GENCODE labels and
  NanoCount counts were used ONLY as training labels - every input feature
  here is observable at inference time.
* The operating point was chosen on the tuning frontier under the hard
  constraint "T1 honest F1 must not drop below the unfiltered profile" and
  frozen before any validation sample was scored.
* ``featurize`` is shared verbatim with the offline fit script, so training
  and inference features cannot drift.
* GTF and fusion candidates are always exempt; only ``source == "novel"``
  rows can be filtered.
"""
from __future__ import annotations

import logging
import math
from typing import Dict, Iterable, Mapping, Tuple

logger = logging.getLogger(__name__)


# Frozen model v1: exact full-precision constants exported from
# experiments/prod_validation/gencode/_goal_opt/ranker_v1/ranker_model.json
# (fit 2026-08-29 on H9 r2r2 precision evidence, 10,217 rows, grouped-CV AUC
# 0.809). Never hand-edit or round these values: the operating threshold is
# defined on this exact score scale.
RANKER_V1 = {
    "features": (
        'log_abundance',
        'log_num_reads',
        'soft_hard_ratio',
        'confidence',
        'family_share',
        'log_family_rank',
        'family_dominant_share',
        'is_subchain_of_sibling',
        'is_superchain_of_sibling',
        'log_weakest_junction',
        'n_junctions_below3',
        'canonical_fraction',
        'end5_support_frac',
        'end3_support_frac',
        'fulllen_frac',
        'ends_missing',
        'is_mono',
        'log_n_exons',
        'log_tx_length',
    ),
    "mean": (
        2.395851523874122,
        2.391541885703071,
        1.0717252422433197,
        0.935693403151609,
        0.7791465427229131,
        0.8141570124653513,
        0.8708796625232478,
        0.009102476265048448,
        0.01115787413134971,
        1.82641496161654,
        0.515121855730645,
        0.9237545267691103,
        0.3721283253401161,
        0.6931453655671934,
        0.2679803366937435,
        0.0001957521777429774,
        0.07624547323088969,
        2.0883552110037145,
        7.642853598039244,
    ),
    "std": (
        0.9632042216238322,
        1.0013137742302232,
        0.5611878668370605,
        0.11904122986137435,
        0.331118750565284,
        0.25451549214451435,
        0.1979792770946097,
        0.09497168625907294,
        0.10503987802838773,
        1.0413346852999372,
        1.6802753376353499,
        0.26539046901250996,
        0.27372517755212533,
        0.30882777952855855,
        0.26053260084400054,
        0.013989776939891454,
        0.26539046901250996,
        0.6967943483388723,
        0.6320766261647113,
    ),
    "weights": (
        0.29856789390542937,
        0.4716416797791935,
        -0.030291409169086174,
        0.040563226259139794,
        0.21812662021432216,
        -0.16851358802113245,
        0.11353462995315117,
        -0.08108621589363393,
        -0.020463962521821113,
        0.5630020961285614,
        -0.04348525298270244,
        0.3613196514845769,
        0.26666136509244204,
        0.05884456530481536,
        -0.0033399860648672852,
        0.023168477069488738,
        -0.37090723311826207,
        -0.13773801045041884,
        -0.03028946653501723,
    ),
    "bias": 1.6616456914194901,
    "score_threshold": 0.5682708140508252,
}


# Frozen training receipts: the identity of the data and script that produced
# the constants above. Mirrored bit-for-bit in
# experiments/prod_validation/models/candidate_ranker_v1.json and asserted
# unconditionally by the unit tests, so the provenance chain survives even if
# the local evidence artifacts are pruned.
RANKER_V1_RECEIPTS = {
    "evidence_sha256": (
        "21143c2cb9112e7029ac38a9c09bce75b2379a6b70d95dafef792e1e36a84725"
    ),
    "evidence_run_source_sha256": (
        "b8999bf9509fd280c5a939cdd039b35b959644c445dc92d6528943c71ac5c81e"
    ),
    "fit_script_sha256": (
        "964fc4122faa8c950d4815a1375de37897574b6df45f45a72b17d6f145a11fcf"
    ),
}


def featurize(row: Mapping) -> list[float]:
    """Evidence row -> model feature vector (shared by fit and inference).

    Accepts both the in-pipeline dict (numeric values) and the TSV row
    (string values); sentinel -1 marks a missing evidence source and maps to
    the neutral value with its indicator, never to a low observation.
    """
    mono = int(row["is_mono"]) == 1
    weakest = float(row["weakest_junction_support"])
    below3 = float(row["n_junctions_below3"])
    canon = float(row["canonical_fraction"])
    e5 = float(row["end5_support_frac"])
    e3 = float(row["end3_support_frac"])
    fl = float(row["fulllen_frac"])
    ends_missing = 1.0 if e5 < 0 else 0.0
    return [
        math.log1p(float(row["abundance"])),
        math.log1p(float(row["num_reads"])),
        min(max(float(row["soft_hard_ratio"]), 0.0), 10.0),
        float(row["confidence"]),
        max(float(row["family_share"]), 0.0),
        math.log1p(float(row["family_rank"])),
        max(float(row["family_dominant_share"]), 0.0),
        float(row["is_subchain_of_sibling"]),
        float(row["is_superchain_of_sibling"]),
        math.log1p(weakest) if weakest >= 0 else 0.0,
        below3 if below3 >= 0 else 0.0,
        canon if canon >= 0 else 0.0,
        max(e5, 0.0), max(e3, 0.0), max(fl, 0.0),
        ends_missing,
        1.0 if mono else 0.0,
        math.log1p(float(row["n_exons"])),
        math.log1p(float(row["tx_length"])),
    ]


def score_evidence_rows(rows: Iterable[Mapping]) -> Dict[str, float]:
    """Frozen-model raw logit per candidate_id (ranking evidence, not a
    calibrated probability)."""
    model = RANKER_V1
    mean, std, w = model["mean"], model["std"], model["weights"]
    scores: Dict[str, float] = {}
    for row in rows:
        x = featurize(row)
        z = model["bias"]
        for i in range(len(w)):
            z += w[i] * (x[i] - mean[i]) / std[i]
        scores[str(row["candidate_id"])] = z
    return scores


def ranking_filter(
    results: Dict[str, object],
    evidence_rows: Iterable[Mapping],
    *,
    threshold: float | None = None,
) -> Tuple[Dict[str, object], Dict[str, float], list[str]]:
    """Drop below-threshold NOVEL candidates; GTF/fusion always exempt.

    Returns ``(kept_results, scores, dropped_ids)``. A candidate without an
    evidence row (should not happen) is kept - missing evidence must never
    become a drop reason.
    """
    thr = RANKER_V1["score_threshold"] if threshold is None else float(threshold)
    scores = score_evidence_rows(evidence_rows)
    dropped = [
        cid for cid, qr in results.items()
        if getattr(qr, "source", "") == "novel"
        and cid in scores and scores[cid] < thr
    ]
    kept = {cid: qr for cid, qr in results.items() if cid not in set(dropped)}
    if dropped:
        logger.info(
            "Candidate ranking filter dropped %d/%d novel candidates below "
            "score %.4f", len(dropped),
            sum(1 for q in results.values() if getattr(q, "source", "") == "novel"),
            thr,
        )
    return kept, scores, sorted(dropped)
