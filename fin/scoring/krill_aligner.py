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
import os
from typing import Optional, Tuple

logger = logging.getLogger(__name__)


def krill_thread_count() -> int:
    """krill ``num_thread`` from the ``KRILL_THREADS`` env var (default 0 = auto).

    A malformed value must not change scoring behavior, so anything that does not
    parse as a non-negative int falls back to ``0`` with a warning rather than
    raising (which would silently disable a signal leg).
    """
    raw = os.environ.get("KRILL_THREADS", "0")
    try:
        n = int(raw)
    except (TypeError, ValueError):
        logger.warning("KRILL_THREADS=%r is not an int; using 0 (auto)", raw)
        return 0
    if n < 0:
        logger.warning("KRILL_THREADS=%r is negative; using 0 (auto)", raw)
        return 0
    return n


_WARNED_KRILL_CPU_ONLY = False


def _safe_gpu_count(krill):
    """``krill.gpu_device_count()`` guarded (returns the raw value or 'err')."""
    try:
        return krill.gpu_device_count()
    except Exception:  # noqa: BLE001
        return "err"


def krill_gpu_available(krill) -> bool:
    """True iff this krill BUILD can actually score on a GPU right now.

    A CPU-only krill build accepts ``use_gpu=True`` but silently scores on the CPU
    (only a ``UserWarning``), so the requested ``use_gpu`` reflects intent, NOT
    reality. The real backend is ``krill.built_with_cuda`` (a bool ATTRIBUTE, not a
    call -- calling it raises ``TypeError``) AND a live device
    (``gpu_device_count() > 0``; a CPU-only build returns ``-1``). Any error ->
    treat as CPU (fail-safe).
    """
    try:
        return bool(getattr(krill, "built_with_cuda", False)) and int(krill.gpu_device_count()) > 0
    except Exception:  # noqa: BLE001 - defensive; unknown krill build
        return False


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
    # HONEST backend: a CPU-only krill build would score on the CPU even when we
    # pass use_gpu=True (silent, only a UserWarning), and returning use_gpu=True
    # would make the caller mislabel the eventalign as GPU. Trust the build, not
    # the request: only claim GPU when krill can actually use one.
    want_gpu = bool(use_gpu)
    effective = want_gpu and krill_gpu_available(krill)
    if want_gpu and not effective:
        global _WARNED_KRILL_CPU_ONLY
        if not _WARNED_KRILL_CPU_ONLY:
            logger.warning(
                "krill --gpu requested but this krill build scores on CPU "
                "(built_with_cuda=%s, gpu_device_count=%s); signal eventalign runs "
                "on CPU. For GPU eventalign install krill==<ver>+cu122 from the "
                "cu122 index. (This warning is emitted once per process.)",
                getattr(krill, "built_with_cuda", False),
                _safe_gpu_count(krill),
            )
            _WARNED_KRILL_CPU_ONLY = True
    try:
        return krill.Aligner(pore=pore, use_gpu=effective, **kwargs), effective
    except Exception as exc:  # noqa: BLE001 - krill raises broad errors
        if effective:
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
