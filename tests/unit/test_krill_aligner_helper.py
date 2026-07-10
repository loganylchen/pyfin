"""Unit tests for the shared ``make_krill_aligner`` helper.

A fake krill module stub lets us exercise the GPU->CPU fallback and the HONEST
backend detection without a real krill install or a GPU. ``FakeKrill`` records
the ``use_gpu`` it was asked for, exposes the real krill backend-introspection
API (``built_with_cuda`` bool ATTRIBUTE + ``gpu_device_count()`` CALL) so the
helper only claims GPU when the build actually supports it, and can be configured
to raise on a GPU build (mimicking a GPU-less host).
"""
from fin import scoring  # noqa: F401 - ensure package import path
import fin.scoring.krill_aligner as ka
from fin.scoring.krill_aligner import make_krill_aligner


class _Aligner:
    def __init__(self, pore, use_gpu, **kwargs):
        self.pore = pore
        self.use_gpu = use_gpu
        self.kwargs = kwargs


class FakeKrill:
    def __init__(self, fail_gpu=False, fail_cpu=False, gpu_capable=True):
        self.fail_gpu = fail_gpu
        self.fail_cpu = fail_cpu
        self.calls = []  # list of use_gpu flags requested
        # Honest-backend introspection the helper now trusts instead of the
        # request: a CPU-only krill build reports built_with_cuda=False and
        # gpu_device_count()=-1 (mirrors the real krill API).
        self.built_with_cuda = gpu_capable
        self._ndev = 1 if gpu_capable else -1

    def gpu_device_count(self):
        return self._ndev

    def Aligner(self, pore, use_gpu, **kwargs):  # noqa: N802 - mimic krill API
        self.calls.append(use_gpu)
        if use_gpu and self.fail_gpu:
            raise RuntimeError("no GPU")
        if not use_gpu and self.fail_cpu:
            raise RuntimeError("CPU build failed")
        return _Aligner(pore, use_gpu, **kwargs)


def test_cpu_request_stays_cpu():
    krill = FakeKrill()
    aligner, eff = make_krill_aligner(krill, "rna002", use_gpu=False)
    assert aligner is not None
    assert eff is False
    assert aligner.use_gpu is False
    assert krill.calls == [False]  # no GPU attempt when CPU requested


def test_gpu_request_succeeds():
    krill = FakeKrill()
    aligner, eff = make_krill_aligner(krill, "rna002", use_gpu=True)
    assert aligner is not None
    assert eff is True
    assert aligner.use_gpu is True
    assert krill.calls == [True]


def test_gpu_init_falls_back_to_cpu():
    krill = FakeKrill(fail_gpu=True)
    aligner, eff = make_krill_aligner(krill, "rna002", use_gpu=True)
    assert aligner is not None
    assert eff is False  # effective device is CPU after fallback
    assert aligner.use_gpu is False
    assert krill.calls == [True, False]  # tried GPU, then CPU


def test_both_fail_returns_none():
    krill = FakeKrill(fail_gpu=True, fail_cpu=True)
    aligner, eff = make_krill_aligner(krill, "rna002", use_gpu=True)
    assert aligner is None
    assert eff is False


def test_cpu_only_build_never_attempts_gpu():
    # The bug this fixes: a CPU-only krill build (built_with_cuda=False) accepts
    # use_gpu=True but silently scores on the CPU. The helper must report the
    # HONEST effective device (False), NOT the request, and must not even ask the
    # aligner for GPU (so callers don't mislabel the eventalign as GPU).
    ka._WARNED_KRILL_CPU_ONLY = False  # deterministic: allow the one-time warning
    krill = FakeKrill(gpu_capable=False)
    aligner, eff = make_krill_aligner(krill, "rna002", use_gpu=True)
    assert aligner is not None
    assert eff is False               # honest: CPU-only build -> not GPU
    assert aligner.use_gpu is False
    assert krill.calls == [False]     # never attempted a GPU build


def test_gpu_available_helper():
    from fin.scoring.krill_aligner import krill_gpu_available
    assert krill_gpu_available(FakeKrill(gpu_capable=True)) is True
    assert krill_gpu_available(FakeKrill(gpu_capable=False)) is False
    assert krill_gpu_available(object()) is False  # unknown build -> CPU (fail-safe)


def test_kwargs_forwarded():
    krill = FakeKrill()
    aligner, _ = make_krill_aligner(
        krill, "rna002", use_gpu=False, hmm_confidence=True, polya=True
    )
    assert aligner.kwargs == {"hmm_confidence": True, "polya": True}
