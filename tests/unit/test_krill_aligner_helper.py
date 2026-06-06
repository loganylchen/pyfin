"""Unit tests for the shared ``make_krill_aligner`` helper.

A fake krill module stub lets us exercise the GPU->CPU fallback logic without a
real krill install or a GPU: ``FakeKrill`` records the ``use_gpu`` it was asked
for and can be configured to raise on a GPU build (mimicking a GPU-less host).
"""
from fin.scoring.krill_aligner import make_krill_aligner


class _Aligner:
    def __init__(self, pore, use_gpu, **kwargs):
        self.pore = pore
        self.use_gpu = use_gpu
        self.kwargs = kwargs


class FakeKrill:
    def __init__(self, fail_gpu=False, fail_cpu=False):
        self.fail_gpu = fail_gpu
        self.fail_cpu = fail_cpu
        self.calls = []  # list of use_gpu flags requested

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


def test_kwargs_forwarded():
    krill = FakeKrill()
    aligner, _ = make_krill_aligner(
        krill, "rna002", use_gpu=False, hmm_confidence=True, polya=True
    )
    assert aligner.kwargs == {"hmm_confidence": True, "polya": True}
