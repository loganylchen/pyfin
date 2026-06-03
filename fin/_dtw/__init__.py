"""
DTW (Dynamic Time Warping) module — krill backend.

This module previously wrapped a self-shipped CUDA DTW extension
(``._cuda_dtw``). It now delegates all DTW work to ``krill`` (same faithful
f5c-engine port used for eventalign), giving one signal dependency for the
whole stack.

Notes
-----
* ``use_open_start`` / ``use_open_end`` are accepted for backward
  compatibility but IGNORED — krill implements standard global DTW (no
  open-boundary parameter) and the project decided open boundaries are not
  needed.
* GPU is used only when krill was built with CUDA AND a device is visible;
  otherwise everything runs on CPU (krill's DTW ``use_gpu=True`` raises when
  no GPU is present, so we guard it).
"""

import logging
import subprocess
from typing import Union

import numpy as np

logger = logging.getLogger(__name__)

# Try to import krill (the DTW backend).
try:
    import krill as _krill

    _KRILL_OK = True
except ImportError as e:  # pragma: no cover - only when krill missing
    _KRILL_OK = False
    _import_error = str(e)


def _gpu_ok() -> bool:
    """True only if krill was built with CUDA and a device is visible."""
    if not _KRILL_OK:
        return False
    try:
        if not bool(getattr(_krill, "built_with_cuda", False)):
            return False
        return _krill.gpu_device_count() > 0
    except Exception:
        return False


_GPU = _gpu_ok()

# Backward-compat flag. True whenever krill DTW is usable (GPU or CPU); callers
# (signal_dtw) use this to route into the batched krill path instead of their
# pure-python CPU fallback.
CUDA_AVAILABLE = _KRILL_OK

if _KRILL_OK:
    logger.info("DTW: krill backend ENABLED (gpu=%s)", _GPU)
else:
    logger.info("DTW: krill not available — callers fall back to CPU")


def _require_krill() -> None:
    if not _KRILL_OK:
        raise RuntimeError(
            "krill DTW backend is not available.\n"
            f"Import error: {_import_error}\n\n"
            "Install with: pip install krill --no-deps "
            "--index-url https://loganylchen.github.io/krill-dist/simple/"
        )


def _varlen(segments) -> np.ndarray:
    """Pairwise DTW over variable-length 1D sequences via krill."""
    _require_krill()
    segs = [np.ascontiguousarray(s, dtype=np.float32) for s in segments]
    out = _krill.dtw_pairwise_varlen(segs, use_gpu=_GPU)
    return np.asarray(out, dtype=np.float64)


def dtw_distance(
    seq1: Union[np.ndarray, list],
    seq2: Union[np.ndarray, list],
    use_open_start: bool = False,
    use_open_end: bool = False,
) -> float:
    """DTW distance between two 1D sequences (krill backend).

    ``use_open_start`` / ``use_open_end`` are accepted but ignored.
    """
    _require_krill()
    a = np.ascontiguousarray(seq1, dtype=np.float32)
    b = np.ascontiguousarray(seq2, dtype=np.float32)
    if a.ndim != 1 or b.ndim != 1:
        raise ValueError("dtw_distance expects 1-dimensional sequences")
    if len(a) == 0 or len(b) == 0:
        raise ValueError("Sequences cannot be empty")
    return float(_krill.dtw_distance(a, b, use_gpu=_GPU))


def dtw_pairwise(
    sequences: Union[np.ndarray, list],
    use_open_start: bool = False,
    use_open_end: bool = False,
) -> np.ndarray:
    """Pairwise DTW for a batch of equal-length sequences (krill backend).

    ``use_open_start`` / ``use_open_end`` are accepted but ignored.
    """
    arr = np.asarray(sequences, dtype=np.float32)
    if arr.ndim != 2:
        raise ValueError(f"sequences must be 2D array, got shape {arr.shape}")
    if arr.shape[0] < 2:
        raise ValueError(f"Need at least 2 sequences, got {arr.shape[0]}")
    if arr.shape[1] == 0:
        raise ValueError("Sequence length cannot be 0")
    return _varlen([arr[i] for i in range(arr.shape[0])])


def dtw_pairwise_varlen(
    segments,
    use_open_start: bool = False,
    use_open_end: bool = False,
) -> np.ndarray:
    """Pairwise DTW for variable-length sequences (krill backend).

    ``use_open_start`` / ``use_open_end`` are accepted but ignored.
    """
    return _varlen(segments)


def cleanup():
    """No-op (krill manages its own resources)."""
    return None


def is_available() -> bool:
    """True if the krill DTW backend is importable."""
    return _KRILL_OK


def estimate_gpu_memory(num_sequences: int, max_length: int) -> int:
    """Estimate GPU bytes needed for pairwise varlen DTW (with 20% headroom)."""
    input_bytes = num_sequences * max_length * 4
    lengths_bytes = num_sequences * 8
    num_pairs = num_sequences * (num_sequences - 1) // 2
    pairs_bytes = num_pairs * 4
    max_parallel = num_sequences - 1
    cost_bytes = max_length * max_parallel * 4 * 2
    return int((input_bytes + lengths_bytes + pairs_bytes + cost_bytes) * 1.2)


def get_free_gpu_memory() -> int:
    """Free GPU memory in bytes (via nvidia-smi), or 0 if unavailable."""
    if not _GPU:
        return 0
    try:
        res = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.free",
             "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=5,
        )
        if res.returncode == 0:
            free_mb = int(res.stdout.strip().split("\n")[0])
            return free_mb * 1024 * 1024
    except (FileNotFoundError, subprocess.TimeoutExpired, ValueError):
        pass
    return 0


__all__ = [
    "dtw_distance",
    "dtw_pairwise",
    "dtw_pairwise_varlen",
    "estimate_gpu_memory",
    "get_free_gpu_memory",
    "cleanup",
    "is_available",
    "CUDA_AVAILABLE",
]
