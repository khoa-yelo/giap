"""Hypothesis properties of the fgsea multilevel stage (issue M8-11; GIAP_SPEC.md §7.7.4, §6.6).

Interface as pinned in tests/unit/stats/test_fgsea_multilevel.py:
``fgsea_multilevel(prep, simple, sample_size, eps, rng) -> SimpleResult``.

For any stats and sets:

- only ``p`` and ``log2err`` of the sets in ``simple.needs_multilevel`` change; every other field
  and every other set is returned as it came;
- a refined p is finite with ``max(eps, tiny) <= p <= 1`` (never 0), and its ``log2err`` is
  ``multilevelError(p, sampleSize)``, NaN when ``isCpGeHalf`` is false, or NaN when p was floored
  at ``eps > 0``;
- the same seed gives the same result.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy.special import polygamma

from giap.rng import for_step
from giap.stats.fgsea import fgsea_simple, prepare

FAST = settings(max_examples=25, deadline=None)
SAMPLE_SIZE = 101
STREAM = "function.gsea"


def _ml_error(p: float) -> float:
    single = float(polygamma(1, (SAMPLE_SIZE + 1) / 2) - polygamma(1, SAMPLE_SIZE + 1))
    return math.sqrt(math.floor(-math.log2(p) + 1) * single) / math.log(2)


@st.composite
def problems(draw: st.DrawFn) -> tuple[np.ndarray, dict[str, list[int]], int, float]:
    n = draw(st.integers(80, 400))
    data_seed = draw(st.integers(0, 2**31 - 1))
    rng = np.random.default_rng(data_seed)
    x = rng.standard_normal(n)
    order = np.argsort(-x, kind="stable")
    n_sets = draw(st.integers(1, 5))
    sets: dict[str, list[int]] = {}
    for j in range(n_sets):
        k = draw(st.integers(5, min(40, n // 3)))
        kind = draw(st.sampled_from(["top", "bottom", "random", "mixed"]))
        if kind == "top":
            idx = order[:k]
        elif kind == "bottom":
            idx = order[::-1][:k]
        elif kind == "random":
            idx = rng.choice(n, size=k, replace=False)
        else:
            idx = np.concatenate([order[: k // 2], rng.choice(order[k:], k - k // 2, False)])
        sets[f"s{j}"] = [int(i) for i in idx]
    seed = draw(st.integers(0, 10**6))
    eps = draw(st.sampled_from([0.0, 0.0, 1e-8, 1e-4]))
    return x, sets, seed, eps


@FAST
@given(problems())
def test_multilevel_changes_only_the_refined_p_and_log2err(
    prob: tuple[np.ndarray, dict[str, list[int]], int, float],
) -> None:
    from giap.stats.fgsea import fgsea_multilevel

    x, sets, seed, eps = prob
    lab = [f"g{i}" for i in range(len(x))]
    prep = prepare(x, lab, {s: [lab[i] for i in v] for s, v in sets.items()}, 1, len(x))
    simple = fgsea_simple(prep, 1000, SAMPLE_SIZE, for_step(STREAM, seed))
    out = fgsea_multilevel(prep, simple, SAMPLE_SIZE, eps, for_step(STREAM, seed + 1))
    again = fgsea_multilevel(prep, simple, SAMPLE_SIZE, eps, for_step(STREAM, seed + 1))
    np.testing.assert_array_equal(out.p, again.p)
    np.testing.assert_array_equal(out.log2err, again.log2err)

    for f in ("es", "nes", "n_more_extreme", "mode_fraction", "needs_multilevel"):
        np.testing.assert_array_equal(getattr(out, f), getattr(simple, f), err_msg=f)
    for a, b in zip(out.leading_edge, simple.leading_edge, strict=True):
        np.testing.assert_array_equal(a, b)
    need = np.asarray(simple.needs_multilevel)
    np.testing.assert_array_equal(out.p[~need], simple.p[~need])
    np.testing.assert_array_equal(out.log2err[~need], simple.log2err[~need])
    for i in np.flatnonzero(need):
        p, err = float(out.p[i]), float(out.log2err[i])
        assert math.isfinite(p) and 0 < p <= 1, p
        assert p >= eps
        if eps > 0 and p == eps:
            assert math.isnan(err)
        elif not math.isnan(err):
            assert err == pytest.approx(_ml_error(p), rel=1e-9)
