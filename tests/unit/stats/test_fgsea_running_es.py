"""M9-17: the GSEA running sum behind V14(a) (docs/specs/M9-17-plot-support-analysis.md, example 1
and the definition "``running_es`` is fgsea's running sum for that set, computed by the same port
function the ``es`` came from; the max-deviation of ``running_es`` equals the row's ``es`` to
1e-9").

Interface assumed (test author): ``giap.stats.fgsea.running_es(stats, members) -> F64`` takes
what ``gsea_es`` takes (``Prepared.stats``: |stats| in rank order; ``members``: distinct 1-based
ranks) and returns the running sum after each of the n ranked positions (scoreType std,
gseaParam 1): a hit at rank i adds |stat_i| / N_R, a miss subtracts 1 / (n - m).
"""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

TOL = 1e-9


def _running(stats: np.ndarray, members: np.ndarray) -> np.ndarray:
    from giap.stats.fgsea import running_es

    out = np.asarray(running_es(np.asarray(stats, float), np.asarray(members, np.int64)))
    assert out.dtype == np.float64
    return out


def _oracle(stats: np.ndarray, members: np.ndarray) -> np.ndarray:
    n, hits = len(stats), set(int(m) for m in members)
    w = np.abs(stats)
    nr = sum(w[i - 1] for i in hits)
    nm = n - len(hits)
    cur, out = 0.0, []
    for i in range(1, n + 1):
        if i in hits:
            cur += (w[i - 1] / nr) if nr > 0 else 1.0 / len(hits)
        else:
            cur -= 1.0 / nm
        out.append(cur)
    return np.array(out)


def _max_dev(r: np.ndarray) -> float:
    lo, hi = float(min(0.0, r.min())), float(max(0.0, r.max()))
    return hi if hi > -lo else lo


def test_example_1_six_tus_set_tu1_tu3() -> None:
    """Centred scores 2, 1, 0.5, -0.5, -1, -2 (already in rank order); set {TU_1, TU_3} = ranks
    {1, 3}. N_R = 2.5, n - m = 4: 0.8, 0.55, 0.75, 0.5, 0.25, 0. ES = 0.8."""
    from giap.stats.fgsea import gsea_es, prepare

    scores = np.array([2.0, 1.0, 0.5, -0.5, -1.0, -2.0])
    labels = [f"TU_{i}" for i in range(1, 7)]
    prep = prepare(scores, labels, {"s": ["TU_1", "TU_3"]}, 1, 10)
    assert prep.members[0].tolist() == [1, 3]
    r = _running(prep.stats, prep.members[0])
    np.testing.assert_allclose(r, [0.8, 0.55, 0.75, 0.5, 0.25, 0.0], rtol=0, atol=TOL)
    es, _ = gsea_es(prep.stats, prep.members[0])
    assert es == pytest.approx(0.8, abs=TOL)
    assert _max_dev(r) == pytest.approx(es, abs=TOL)


def test_negative_set_runs_below_zero() -> None:
    """Set {TU_5, TU_6} at the bottom: the walk dips to -1 before the hits, ES = -1."""
    from giap.stats.fgsea import gsea_es

    stats = np.array([2.0, 1.0, 0.5, 0.5, 1.0, 2.0])  # |scores| in rank order
    members = np.array([5, 6])
    r = _running(stats, members)
    np.testing.assert_allclose(r, [-0.25, -0.5, -0.75, -1.0, -2 / 3, 0.0], rtol=0, atol=TOL)
    es, _ = gsea_es(stats, members)
    assert es == pytest.approx(-1.0, abs=TOL)
    assert _max_dev(r) == pytest.approx(es, abs=TOL)


def test_rejects_what_gsea_es_rejects() -> None:
    from giap.stats.fgsea import running_es

    stats = np.array([3.0, 2.0, 1.0])
    for bad in ([], [0], [4], [2, 2], [1, 2, 3]):
        with pytest.raises(ValueError):
            running_es(stats, np.array(bad, dtype=np.int64))


@st.composite
def _case(draw: st.DrawFn) -> tuple[np.ndarray, np.ndarray]:
    n = draw(st.integers(3, 60))
    floats = st.floats(-5, 5, allow_nan=False, allow_infinity=False, width=32)
    vals = draw(st.lists(floats, min_size=n, max_size=n))
    stats = np.sort(np.abs(np.array(vals, dtype=float)))[::-1].copy()
    m = draw(st.integers(1, n - 1))
    members = np.array(sorted(draw(st.sets(st.integers(1, n), min_size=m, max_size=m))))
    return stats, members


@settings(max_examples=200, deadline=None)
@given(_case())
def test_property_running_sum_matches_oracle_and_es(case: tuple[np.ndarray, np.ndarray]) -> None:
    """For any ranking and set: n values, the oracle's walk, ends at 0, and its max deviation is
    ``gsea_es``'s ES (fgsea's ``maxP == -minP`` tie gives ES 0, so it is skipped)."""
    from giap.stats.fgsea import gsea_es

    stats, members = case
    r = _running(stats, members)
    assert r.shape == (len(stats),)
    np.testing.assert_allclose(r, _oracle(stats, members), rtol=0, atol=1e-9)
    assert abs(r[-1]) <= 1e-9
    es, _ = gsea_es(stats, members)
    hi, lo = float(max(0.0, r.max())), float(min(0.0, r.min()))
    if abs(hi + lo) > 1e-9:
        assert _max_dev(r) == pytest.approx(es, abs=1e-9)


@settings(max_examples=50, deadline=None)
@given(_case())
def test_property_deterministic_and_pure(case: tuple[np.ndarray, np.ndarray]) -> None:
    stats, members = case
    s0, m0 = stats.copy(), members.copy()
    a = _running(stats, members)
    b = _running(stats, members)
    assert a.tobytes() == b.tobytes()
    assert np.array_equal(stats, s0) and np.array_equal(members, m0)
