"""Shared krill aligner builder honoring ``use_gpu`` with a CPU fallback.

Every production krill eventalign step requests GPU when the pipeline is run with
``use_gpu=True`` (CLI ``--gpu``), but a GPU-less host (or a GPU build error) must
NOT abort the run -- it falls back to CPU. krill takes ``use_gpu`` on BOTH
``krill.Aligner(...)`` and the matching ``align_read_variants(...)`` call, so the
two must agree on the device: if the aligner GPU-init falls back to CPU, the
eventalign call has to use CPU too. This helper centralizes that logic and
returns the *effective* device so callers can keep the two in sync.
"""
from __future__ import annotations

import logging
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


def make_krill_aligner(
    krill,
    pore: str,
    use_gpu: bool,
    **kwargs,
) -> Tuple[Optional[object], bool]:
    """Build a krill aligner honoring ``use_gpu``, GPU->CPU fallback on init error.

    Args:
        krill: imported krill module.
        pore: krill pore model.
        use_gpu: request GPU; on init failure retry once on CPU.
        **kwargs: forwarded verbatim to ``krill.Aligner`` (e.g.
            ``hmm_confidence``, ``polya``, ``num_thread``).

    Returns:
        ``(aligner, effective_gpu)``. ``effective_gpu`` reflects the device the
        returned aligner actually built on, so the matching
        ``align_read_variants`` call can pass the same flag. ``(None, False)``
        only when even the CPU build fails -- caller should treat that as
        "krill unavailable" and degrade gracefully.
    """
    try:
        return krill.Aligner(pore=pore, use_gpu=use_gpu, **kwargs), use_gpu
    except Exception as exc:  # noqa: BLE001 - krill raises broad errors
        if use_gpu:
            logger.warning(
                "krill GPU aligner init failed (%s); retrying on CPU", exc
            )
            try:
                return krill.Aligner(pore=pore, use_gpu=False, **kwargs), False
            except Exception as exc2:  # noqa: BLE001
                logger.warning("krill CPU aligner init failed (%s); skipping", exc2)
                return None, False
        logger.warning("krill aligner init failed (%s); skipping", exc)
        return None, False
