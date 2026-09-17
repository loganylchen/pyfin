"""The endpoint EM must produce bit-identical parameters after hoisting.

`fit_endpoint_model` runs once PER INTERVAL over the whole library's endpoint
distances, excluding only that interval's own reads -- about 127,646
observations on a real human library, up to 200 iterations. Timed with one span
around `_endpoint_context`, it cost ~7.0 s per interval, FLAT against the
interval's own read count (1 read cost the same as 111), and accounts for
essentially all of a genome-wide cell's 14.6 hours.

`_laplace` recomputed `math.exp(-1/b)` and `(1-q)/(1+q)` for every observation
even though both depend only on `b`, which is fixed within an iteration, and
`abs(x)` was recomputed on every one of up to 200 iterations.

Hoisting those is EXACT rather than approximate, and this file holds it to that:
a deterministic function of the same input returns the same float, and `q ** a`
memoised on the same pair returns the same float. What is NOT safe is
re-associating -- Codex verified on the live function that merely reversing a
200-distance sequence moves the fitted scales in the last bits
(0x1.6d789b108b1aep+2 against 0x1.6d789b108b1aap+2), so this EM depends on the
ORDER of its input, not just its multiset.
"""
import math
import random

import pytest

from fin.analysis.sample_nuisance import EndpointModel, fit_endpoint_model


def _hexes(model):
    return (float.hex(model.b_narrow), float.hex(model.b_broad),
            float.hex(model.eta), model.n_fit)


def _reference(distances, *, max_iter=200, tol=1e-8, eta_max=None):
    """The pre-hoist recurrence, verbatim, as the oracle."""
    from fin.analysis.sample_nuisance import BROAD_SCALE_MAX_BP, ETA_MAX

    eta_max = ETA_MAX if eta_max is None else eta_max
    d = [int(x) for x in distances]
    n = len(d)
    if n < 2:
        return EndpointModel(b_narrow=1.0, b_broad=100.0, eta=0.0, n_fit=n)
    absd = sorted(abs(x) for x in d)
    med = absd[n // 2]
    b_n = max(float(med), 0.5)
    b_b = max(float(absd[int(0.95 * (n - 1))]), 10.0 * b_n, 20.0)
    eta = 0.05
    for _ in range(max_iter):
        w_sum = 0.0; s_n = 0.0; s_b = 0.0; n_n = 0.0; n_b = 0.0
        for x in d:
            gn = (1 - eta) * EndpointModel._laplace(x, b_n)
            gb = eta * EndpointModel._laplace(x, b_b)
            z = gn + gb
            if z <= 0:
                continue
            r = gb / z
            w_sum += r
            s_n += (1 - r) * abs(x); n_n += (1 - r)
            s_b += r * abs(x); n_b += r
        new_eta = min(max(w_sum / n, 1e-4), eta_max)
        new_bn = max(s_n / max(n_n, 1e-9), 0.5)
        new_bb = min(max(s_b / max(n_b, 1e-9), 2.0 * new_bn), BROAD_SCALE_MAX_BP)
        delta = max(abs(new_eta - eta), abs(new_bn - b_n), abs(new_bb - b_b))
        eta, b_n, b_b = new_eta, new_bn, new_bb
        if delta < tol:
            break
    return EndpointModel(b_narrow=b_n, b_broad=b_b, eta=eta, n_fit=n)


def _sample(generator, n, spread):
    narrow = [int(generator.gauss(0, spread)) for _ in range(int(n * 0.9))]
    broad = [int(generator.gauss(0, spread * 40)) for _ in range(n - len(narrow))]
    values = narrow + broad
    generator.shuffle(values)
    return values


@pytest.mark.parametrize("seed,n,spread", [
    (1, 200, 3), (2, 1000, 8), (3, 5000, 25), (4, 300, 1), (5, 2000, 120),
])
def test_bit_identical_against_the_pre_hoist_recurrence(seed, n, spread):
    values = _sample(random.Random(seed), n, spread)
    assert _hexes(fit_endpoint_model(values)) == _hexes(_reference(values))


def test_the_inlined_kernel_reproduces_laplace_exactly():
    """The hot loop no longer calls `_laplace`; it reproduces its arithmetic.
    If `_laplace` is ever changed, this fails rather than silently drifting."""
    for b in (0.5, 1.0, 3.7, 40.0, 1e-7, 250.0):
        q = math.exp(-1.0 / max(b, 1e-6))
        c = (1.0 - q) / (1.0 + q)
        for x in (0, 1, -1, 7, -13, 250, -1000):
            assert float.hex(c * (q ** abs(x))) == float.hex(EndpointModel._laplace(x, b))


def test_the_grouping_is_load_bearing():
    """`(1 - eta) * (c * v)` is what the original computed, and the hoist keeps
    those parentheses. The algebraically equal `((1 - eta) * c) * v` is not
    reliably the same float: over 400,000 random draws in this shape the two
    groupings disagree in the last bit about 34% of the time. A single example
    would have been a coin flip -- my first version of this test asserted on one
    and it happened to agree -- so assert on the RATE."""
    generator = random.Random(7)
    differ = 0
    for _ in range(20000):
        eta, c, v = generator.random() * 0.2, generator.random(), generator.random()
        differ += (1 - eta) * (c * v) != ((1 - eta) * c) * v
    assert differ > 2000
    eta, c, v = 0.06476655296663247, 0.15084917392450192, 0.6509344730398537
    assert float.hex((1 - eta) * (c * v)) == "0x1.78263478a5acdp-4"
    assert float.hex(((1 - eta) * c) * v) == "0x1.78263478a5accp-4"


def test_order_changes_the_fit_so_the_multiset_is_not_the_key():
    """Recorded because it disqualifies a whole class of optimisation: any
    rewrite that reorders or aggregates observations is NOT byte-identical."""
    values = _sample(random.Random(11), 400, 6)
    assert _hexes(fit_endpoint_model(values)) != _hexes(fit_endpoint_model(values[::-1]))


@pytest.mark.parametrize("values", [[], [5], [0, 0], [3, 3, 3, 3]])
def test_degenerate_inputs_match(values):
    assert _hexes(fit_endpoint_model(values)) == _hexes(_reference(values))


def test_zero_distances_do_not_divide_by_zero():
    model = fit_endpoint_model([0] * 50)
    assert model.b_narrow >= 0.5 and model.n_fit == 50
