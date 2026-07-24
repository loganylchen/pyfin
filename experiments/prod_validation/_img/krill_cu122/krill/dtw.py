"""Dynamic Time Warping (DTW) pairwise-distance primitive for krill.

This is an *additive downstream primitive*, separate from the eventalign /
resquiggle aligner (``Aligner``). It computes DTW **distances** between raw
signal segments (read-vs-read), e.g. for clustering per-reference-position
signals when detecting RNA modifications. It does NOT perform signal-to-
reference alignment and does not touch the ABEA aligner.

Backends
--------
* GPU: the vendored cuDTW++ warp-shuffle kernel (``krill._dtw``), built only
  when krill is compiled with CUDA. Signals are padded to a fixed length
  bucket in {127, 255, 511, 1023, 2047}; signals longer than 2047 are resampled
  (scipy) down to 2047 before the GPU call (lossy).
* CPU: a self-contained full DTW (``_krill.dtw_distance_cpu``), always
  available. Numerically consistent with the GPU kernel by construction
  (squared-euclidean local cost, 3-neighbour min, final sqrt). For consistency
  with the GPU, CPU inputs of length <= 2047 are padded to the same length
  bucket the GPU would use. Signals longer than 2047 are run at full length on
  the CPU (no resampling), so CPU and GPU results diverge in that regime.

The local-cost definition is identical for both backends; only the >2047
resample step differs.

Attribution
-----------
The GPU kernel is vendored from cuDTW++ (Apache-2.0); see
``src/dtw/cudtw/LICENSE`` in the source tree.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Backend detection
# ---------------------------------------------------------------------------

try:
    from . import _dtw as _gpu  # cuDTW++ warp-shuffle C-API module
    _GPU_EXTENSION_AVAILABLE = True
except ImportError:
    _gpu = None
    _GPU_EXTENSION_AVAILABLE = False

from ._krill import dtw_distance_cpu as _dtw_distance_cpu_ext
from ._krill import dtw_pairwise_cpu as _dtw_pairwise_cpu_ext
from ._krill import gpu_device_count as _gpu_device_count


_BUCKETS = (127, 255, 511, 1023, 2047)
_MAX_GPU_LEN = 2047


def _device_count():
    """Number of usable CUDA devices right now (0 on any error). Mirrors the
    runtime detection used by ``Aligner`` so DTW and the aligner agree."""
    try:
        n = _gpu_device_count() if callable(_gpu_device_count) else _gpu_device_count
        return int(n)
    except Exception:
        return 0


def _gpu_runtime_available():
    """True iff the GPU extension is importable AND a CUDA device is present.
    A CUDA-built wheel on a no-GPU host (login node, CI, Docker without
    ``--gpus``) reports False here, so auto-selection falls back to CPU."""
    return _GPU_EXTENSION_AVAILABLE and _device_count() > 0


def extension_available():
    """True if the GPU (cuDTW++) DTW extension was compiled into this build,
    regardless of whether a CUDA device is present at runtime."""
    return _GPU_EXTENSION_AVAILABLE


def backend():
    """Return the backend auto-selection will use right now: 'gpu' only if the
    cuDTW++ extension is built AND a CUDA device is available, else 'cpu'."""
    return "gpu" if _gpu_runtime_available() else "cpu"


def is_available():
    """True if GPU DTW is usable *now* (extension built and a device present),
    not merely compiled in. Matches ``Aligner`` runtime GPU detection."""
    return _gpu_runtime_available()


def cleanup():
    """Reset the CUDA device and free GPU DTW resources. No-op on CPU."""
    if _GPU_EXTENSION_AVAILABLE:
        _gpu.cleanup()


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _select_bucket(n):
    for b in _BUCKETS:
        if n <= b:
            return b
    return _BUCKETS[-1]  # caller resamples before reaching here


def _resample(sig, target_len):
    """Resample a 1-D signal to target_len (scipy, lazy import)."""
    from scipy.signal import resample
    return resample(sig.astype(np.float64), target_len).astype(np.float32)


def _as_1d_f32(sig, name="signal"):
    arr = np.ascontiguousarray(np.asarray(sig, dtype=np.float32))
    if arr.ndim != 1:
        raise ValueError("%s must be 1-dimensional, got shape %r"
                         % (name, arr.shape))
    if arr.size == 0:
        raise ValueError("%s cannot be empty" % name)
    return arr


def _pad_to(sig, width):
    """Front-load sig into a zero-filled float32 array of length width.
    Caller guarantees len(sig) <= width."""
    out = np.zeros(width, dtype=np.float32)
    out[:sig.shape[0]] = sig
    return out


def _resolve_use_gpu(use_gpu):
    """None -> auto (GPU iff a device is available now); True -> require a
    working GPU (raise otherwise); False -> CPU. Same contract as ``Aligner``."""
    if use_gpu is None:
        return _gpu_runtime_available()
    if use_gpu:
        if not _GPU_EXTENSION_AVAILABLE:
            raise RuntimeError(
                "use_gpu=True but the cuDTW++ GPU extension (krill._dtw) is "
                "not available; this krill build has no GPU DTW support. "
                "Use use_gpu=False for the CPU backend.")
        if _device_count() <= 0:
            raise RuntimeError(
                "use_gpu=True but no CUDA device is available "
                "(gpu_device_count() == 0). Use use_gpu=False or use_gpu=None "
                "(auto) for the CPU backend.")
        return True
    return False


def _cpu_width(global_max):
    """CPU target width: the GPU bucket for <=2047, else full length."""
    return _select_bucket(global_max) if global_max <= _MAX_GPU_LEN else global_max


def _cpu_pairwise_matrix(padded, num_threads=0):
    """Symmetric DTW matrix over a list of equal-length float32 arrays (CPU).
    Zero diagonal, sqrt(full-DTW-cost) off-diagonal — matching the GPU.
    Delegates to the multithreaded, GIL-released C++ pairwise kernel; the
    per-pair numerics are identical to looping ``dtw_distance_cpu``.
    ``num_threads`` <= 0 auto-selects the core count; > 0 pins an exact count."""
    n = len(padded)
    if n < 2:
        return np.zeros((n, n), dtype=np.float64)
    seqs = np.ascontiguousarray(np.stack(padded), dtype=np.float32)
    return np.asarray(_dtw_pairwise_cpu_ext(seqs, num_threads), dtype=np.float64)


# ---------------------------------------------------------------------------
# dtw_distance
# ---------------------------------------------------------------------------

def dtw_distance(seq1, seq2, use_gpu=None):
    """DTW distance between two 1-D signals.

    Parameters
    ----------
    seq1, seq2 : array-like
        1-D signals (converted to float32).
    use_gpu : bool or None
        None=auto (GPU iff available), True=require GPU, False=CPU.

    Returns
    -------
    float
    """
    seq1 = _as_1d_f32(seq1, "seq1")
    seq2 = _as_1d_f32(seq2, "seq2")
    on_gpu = _resolve_use_gpu(use_gpu)

    global_max = max(seq1.shape[0], seq2.shape[0])

    if on_gpu:
        if global_max > _MAX_GPU_LEN:
            if seq1.shape[0] > _MAX_GPU_LEN:
                seq1 = _resample(seq1, _MAX_GPU_LEN)
            if seq2.shape[0] > _MAX_GPU_LEN:
                seq2 = _resample(seq2, _MAX_GPU_LEN)
        return float(_gpu.dtw_distance(seq1, seq2))

    width = _cpu_width(global_max)
    return float(_dtw_distance_cpu_ext(_pad_to(seq1, width), _pad_to(seq2, width)))


# ---------------------------------------------------------------------------
# dtw_pairwise (equal-length batch)
# ---------------------------------------------------------------------------

def dtw_pairwise(sequences, use_gpu=None, num_threads=0):
    """Pairwise DTW distance matrix for a batch of equal-length signals.

    Parameters
    ----------
    sequences : array-like, shape (N, L)
        N signals, all of length L.
    use_gpu : bool or None
    num_threads : int
        CPU-only worker count (<=0 auto = core count); ignored on GPU.

    Returns
    -------
    np.ndarray, shape (N, N), float64
        Symmetric, zero diagonal.
    """
    seqs = np.ascontiguousarray(np.asarray(sequences, dtype=np.float32))
    if seqs.ndim != 2:
        raise ValueError("sequences must be 2-D (N, L), got shape %r" % (seqs.shape,))
    if seqs.shape[0] < 2:
        raise ValueError("need at least 2 sequences, got %d" % seqs.shape[0])
    if seqs.shape[1] == 0:
        raise ValueError("sequence length cannot be 0")

    on_gpu = _resolve_use_gpu(use_gpu)
    L = seqs.shape[1]

    if on_gpu:
        if L > _MAX_GPU_LEN:
            seqs = np.ascontiguousarray(
                np.stack([_resample(s, _MAX_GPU_LEN) for s in seqs]))
        return np.asarray(_gpu.dtw_pairwise(seqs), dtype=np.float64)

    width = _cpu_width(L)
    padded = [_pad_to(s, width) for s in seqs]
    return _cpu_pairwise_matrix(padded, num_threads)


# ---------------------------------------------------------------------------
# dtw_pairwise_varlen (variable-length)
# ---------------------------------------------------------------------------

def dtw_pairwise_varlen(signals, use_gpu=None, num_threads=0):
    """Pairwise DTW distance matrix for variable-length signals.

    Parameters
    ----------
    signals : list of 1-D array-like
        N signals of possibly differing lengths.
    use_gpu : bool or None
    num_threads : int
        CPU-only worker count (<=0 auto = core count); ignored on GPU.

    Returns
    -------
    np.ndarray, shape (N, N), float64
        Symmetric, zero diagonal.
    """
    prepped = [_as_1d_f32(s, "signals[%d]" % i) for i, s in enumerate(signals)]
    if len(prepped) < 2:
        raise ValueError("need at least 2 signals, got %d" % len(prepped))

    on_gpu = _resolve_use_gpu(use_gpu)

    if on_gpu:
        if max(s.shape[0] for s in prepped) > _MAX_GPU_LEN:
            prepped = [_resample(s, _MAX_GPU_LEN) if s.shape[0] > _MAX_GPU_LEN else s
                       for s in prepped]
        lengths = np.array([s.shape[0] for s in prepped], dtype=np.int64)
        max_len = int(lengths.max())
        padded = np.zeros((len(prepped), max_len), dtype=np.float32)
        for i, s in enumerate(prepped):
            padded[i, :s.shape[0]] = s
        return np.asarray(_gpu.dtw_pairwise_varlen(padded, lengths), dtype=np.float64)

    width = _cpu_width(max(s.shape[0] for s in prepped))
    padded = [_pad_to(s, width) for s in prepped]
    return _cpu_pairwise_matrix(padded, num_threads)


# ---------------------------------------------------------------------------
# dtw_multi_position_pairwise (batched per-position pairwise)
# ---------------------------------------------------------------------------

def dtw_multi_position_pairwise(position_signals, use_gpu=None,
                                num_streams=16, device_id=0, num_threads=0):
    """Pairwise DTW matrices for several positions in one batched call.

    Parameters
    ----------
    position_signals : list of list of 1-D array-like
        position_signals[p][r] is the signal for position p, read r.
    use_gpu : bool or None
    num_streams, device_id : int
        GPU-only tuning (ignored on CPU).
    num_threads : int
        CPU-only worker count (<=0 auto = core count); ignored on GPU.

    Returns
    -------
    list of np.ndarray
        One (n_p, n_p) float64 matrix per position.
    """
    if len(position_signals) < 1:
        raise ValueError("need at least 1 position, got 0")

    prepped = []
    counts = []
    for p, pos_sigs in enumerate(position_signals):
        ps = [_as_1d_f32(s, "position_signals[%d][%d]" % (p, i))
              for i, s in enumerate(pos_sigs)]
        prepped.append(ps)
        counts.append(len(ps))

    on_gpu = _resolve_use_gpu(use_gpu)

    if on_gpu:
        if any(s.shape[0] > _MAX_GPU_LEN for ps in prepped for s in ps):
            prepped = [[_resample(s, _MAX_GPU_LEN) if s.shape[0] > _MAX_GPU_LEN else s
                        for s in ps] for ps in prepped]
        global_max = max((s.shape[0] for ps in prepped for s in ps), default=0)
        total = sum(counts)
        padded = np.zeros((total, global_max), dtype=np.float32)
        lengths = np.empty(total, dtype=np.int64)
        idx = 0
        for ps in prepped:
            for s in ps:
                padded[idx, :s.shape[0]] = s
                lengths[idx] = s.shape[0]
                idx += 1
        counts_arr = np.array(counts, dtype=np.int64)
        flat = _gpu.dtw_multi_position_pairwise(
            padded, lengths, counts_arr,
            num_cuda_streams=num_streams, device_id=device_id)
        out = []
        off = 0
        for n in counts:
            out.append(np.asarray(flat[off:off + n * n],
                                  dtype=np.float64).reshape(n, n))
            off += n * n
        return out

    # CPU: one shared bucket across all positions (matches the GPU's global width)
    global_max = max((s.shape[0] for ps in prepped for s in ps), default=0)
    width = _cpu_width(global_max)
    out = []
    for ps in prepped:
        if len(ps) < 2:
            out.append(np.zeros((len(ps), len(ps)), dtype=np.float64))
            continue
        out.append(_cpu_pairwise_matrix([_pad_to(s, width) for s in ps],
                                        num_threads))
    return out


__all__ = [
    "dtw_distance",
    "dtw_pairwise",
    "dtw_pairwise_varlen",
    "dtw_multi_position_pairwise",
    "backend",
    "is_available",
    "cleanup",
]
