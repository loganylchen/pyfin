"""krill whole-read polyA estimation (RNA002).

The 3' polyA tail is only in-frame in a WHOLE-read eventalign, so this runs krill
with ``polya=True`` over each read's full basecalled sequence (variant label
"self") and returns the estimated tail length + QC verdict per read. Used by the
polyA + 5'-proximity candidate-retention filter in the pipeline finalize stage.
"""
from __future__ import annotations

import logging
from typing import Dict, Optional, Tuple

from fin.scoring.krill_aligner import make_krill_aligner

logger = logging.getLogger(__name__)


def compute_polya(
    read_seqs: Dict[str, str],
    signal_path: str,
    pore: str = "rna002",
    use_gpu: bool = True,
    batch: int = 256,
) -> Dict[str, Tuple[Optional[float], Optional[str]]]:
    """Estimate per-read polyA length / QC via krill whole-read eventalign.

    Args:
        read_seqs: ``{read_id: full_read_sequence}``. Empty sequences are skipped.
        signal_path: SLOW5/BLOW5 (or POD5) signal file.
        pore: krill pore model (default ``rna002``).
        use_gpu: run krill on GPU when available.
        batch: reads per ``align_reads_variants`` call (bounds GPU memory).

    Returns:
        ``{read_id: (polya_length, polya_qc)}`` for reads krill scored. A read
        absent from the result was not scored (treat as no polyA support).

    Degrades gracefully: a missing/old krill, a GPU-unavailable host, an
    unsupported signal file, or a per-batch krill error never aborts the run —
    it returns whatever was scored (possibly empty), and the caller treats an
    empty map as "filter disabled" (no-op). With ``use_gpu=True`` on a host
    without a usable GPU, krill is retried once on CPU so the filter still runs.
    """
    try:
        import krill
    except Exception as exc:  # krill not installed / import failure
        logger.warning("polyA filter: krill import failed (%s); skipping", exc)
        return {}

    items = [(rid, seq) for rid, seq in read_seqs.items() if seq]
    out: Dict[str, Tuple[Optional[float], Optional[str]]] = {}
    if not items:
        return out

    aligner, eff_gpu = _make_polya_aligner(krill, pore, use_gpu)
    if aligner is None:
        return out

    for i in range(0, len(items), batch):
        chunk = items[i : i + batch]
        reads_variants = {rid: {"self": seq} for rid, seq in chunk}
        try:
            results = krill.align_reads_variants(
                signal_path, reads_variants, aligner=aligner, use_gpu=eff_gpu
            )
        except Exception as exc:  # bad signal file / per-batch krill error
            # A partial failure is unrecoverable: reads in the failed batch
            # would be indistinguishable from reads with no polyA support, so
            # the filter could silently DROP candidates (incl. annotated GTF
            # ones) on a krill error rather than on evidence. Disable the whole
            # filter (return {} -> caller no-ops) — recall-safe.
            logger.warning(
                "polyA filter: krill batch %d failed (%s); disabling filter (no-op)",
                i // batch, exc,
            )
            return {}
        for rid, res_list in results.items():
            if not res_list:
                continue
            r = res_list[0]
            out[rid] = (r.get("polya_length"), r.get("polya_qc"))
        logger.debug(
            "krill polyA: %d/%d reads", min(i + batch, len(items)), len(items)
        )
    return out


def _make_polya_aligner(krill, pore: str, use_gpu: bool):
    """Build a krill polyA aligner, retrying on CPU if a GPU build fails.

    Returns ``(aligner, effective_gpu)``. ``effective_gpu`` is the device the
    aligner actually built on, so the matching ``align_reads_variants`` call uses
    the same device. ``aligner`` is ``None`` (filter no-ops) only when even the
    CPU aligner cannot be constructed, so a GPU-less host still runs the filter
    rather than aborting.
    """
    return make_krill_aligner(
        krill, pore, use_gpu, hmm_confidence=True, polya=True
    )
