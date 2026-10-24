"""``_ranks``, ``_rho`` and ``_auroc`` against their library oracles (PR #286 review, B1.2).

* ``_ranks(x)`` equals ``scipy.stats.rankdata(x, method="average")`` exactly (ties, ±0, n up to
  a few thousand);
* ``_rho(x, y)`` equals ``scipy.stats.spearmanr(x, y).statistic`` exactly when both inputs vary
  and n >= ``MIN_UNIVARIATE_ROWS``; NaN when either is constant or n is below that;
* ``_auroc(x, y)`` equals the Mann-Whitney pair count with ties ½, ``U / (n1 n0)``, exactly: as a
  correctly rounded ``Fraction`` (``_auroc`` makes one rounding, the final division; its rank sums
  are exact half-integers) and as ``scipy.stats.mannwhitneyu(...).statistic / (n1 n0)``. It also
  equals ``sklearn.metrics.roc_auc_score`` to 1e-12, not bit for bit: sklearn integrates the ROC
  curve with the trapezoid rule, a different float evaluation of the same rational (it differs in
  the last ulp on about a third of random inputs).
"""

from __future__ import annotations

import math
import warnings
from fractions import Fraction

import numpy as np
import numpy.typing as npt
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from scipy import stats
from sklearn.metrics import roc_auc_score

from giap.determinants.models import MIN_UNIVARIATE_ROWS, _auroc, _ranks, _rho

SETTINGS = settings(
    max_examples=80,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large],
)

F64 = npt.NDArray[np.float64]

# Few distinct values, including both zeros, so ties (and ±0 ties) are frequent.
TIED = st.sampled_from([-0.0, 0.0, 1.0, -1.0, 2.5, 1e-300, -1e300])
ANY = st.floats(allow_nan=False, allow_infinity=False, width=64)


@st.composite
def vectors(draw: st.DrawFn, min_size: int = 1) -> F64:
    """Either a Hypothesis list (adversarial values, heavy ties, ±0) or a seeded numpy draw of
    up to 3,000 values on a coarse grid (large n, many ties)."""
    if draw(st.booleans()):
        elems = draw(st.sampled_from([TIED, ANY, st.one_of(TIED, ANY)]))
        return np.asarray(draw(st.lists(elems, min_size=min_size, max_size=60)), dtype=np.float64)
    n = draw(st.integers(min_value=min_size, max_value=3000))
    r = np.random.default_rng(draw(st.integers(min_value=0, max_value=2**32 - 1)))
    levels = draw(st.integers(min_value=1, max_value=50))
    x = r.integers(0, levels, n).astype(np.float64) - levels // 2
    # Signed zeros mixed in: -0.0 == 0.0 must tie.
    return np.where((x == 0) & (r.random(n) < 0.5), -0.0, x)


@st.composite
def pairs(draw: st.DrawFn) -> tuple[F64, F64]:
    x = draw(vectors(min_size=MIN_UNIVARIATE_ROWS))
    r = np.random.default_rng(draw(st.integers(min_value=0, max_value=2**32 - 1)))
    kind = draw(st.sampled_from(["normal", "tied", "monotone", "const"]))
    if kind == "normal":
        y = r.normal(size=x.size)
    elif kind == "tied":
        y = r.integers(0, 3, x.size).astype(np.float64)
    elif kind == "monotone":
        y = np.where(np.isfinite(x * 2.0), x * 2.0, x)
    else:
        y = np.full(x.size, 1.5)
    return x, y


@st.composite
def labelled(draw: st.DrawFn) -> tuple[F64, npt.NDArray[np.bool_]]:
    """Scores with ties and both classes present, total size >= MIN_UNIVARIATE_ROWS."""
    n = draw(st.integers(min_value=MIN_UNIVARIATE_ROWS, max_value=400))
    r = np.random.default_rng(draw(st.integers(min_value=0, max_value=2**32 - 1)))
    levels = draw(st.integers(min_value=1, max_value=12))
    x = (r.integers(0, levels, n) / 4.0).astype(np.float64)
    if draw(st.booleans()):
        x = np.where(x == 0, -0.0, x)
    n1 = draw(st.integers(min_value=1, max_value=n - 1))
    y = np.zeros(n, dtype=bool)
    y[r.permutation(n)[:n1]] = True
    return x, y


@SETTINGS
@given(vectors())
def test_ranks_equal_rankdata_average(x: F64) -> None:
    got = _ranks(x)
    want = stats.rankdata(x, method="average")
    assert got.dtype == np.float64
    assert np.array_equal(got, want)


def test_ranks_examples() -> None:
    assert np.array_equal(_ranks(np.array([3.0, 1.0, 2.0])), [3.0, 1.0, 2.0])
    # Ties share the mean rank; -0.0 ties with 0.0.
    assert np.array_equal(
        _ranks(np.array([0.0, -0.0, 5.0, 5.0, 5.0, -1.0])), [2.5, 2.5, 5, 5, 5, 1]
    )
    assert np.array_equal(_ranks(np.full(4, 7.0)), [2.5, 2.5, 2.5, 2.5])
    assert np.array_equal(_ranks(np.array([1.0])), [1.0])


def _spearman(x: F64, y: F64) -> float:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # ConstantInputWarning
        return float(stats.spearmanr(x, y).statistic)


@SETTINGS
@given(pairs())
def test_rho_equals_spearmanr(xy: tuple[F64, F64]) -> None:
    x, y = xy
    with np.errstate(over="ignore"):  # ptp of +-1e300-scale values overflows to inf (non-zero)
        got = _rho(x, y)
        constant = np.ptp(x) == 0 or np.ptp(y) == 0
    if constant:
        assert math.isnan(got)
        assert math.isnan(_spearman(x, y))
    else:
        assert got == _spearman(x, y)


def test_rho_nan_for_constant_input_and_short_input() -> None:
    n = MIN_UNIVARIATE_ROWS
    v = np.arange(n, dtype=np.float64)
    assert math.isnan(_rho(np.full(n, 2.0), v))
    assert math.isnan(_rho(v, np.full(n, 2.0)))
    assert math.isnan(_rho(np.where(v % 2 == 0, 0.0, -0.0), v))  # ±0 only: constant
    assert math.isnan(_rho(v[: n - 1], v[: n - 1]))
    # Exact identity with scipy, including its last-ulp 0.9999999999999999 for rho = 1.
    assert _rho(v, v) == _spearman(v, v)
    assert _rho(v, -v) == _spearman(v, -v)
    assert abs(_rho(v, v) - 1.0) <= 1e-15


def _pair_count(x: F64, y: npt.NDArray[np.bool_]) -> Fraction:
    pos, neg = x[y][:, None], x[~y][None, :]
    twice_u = 2 * int((pos > neg).sum()) + int((pos == neg).sum())
    return Fraction(twice_u, 2 * pos.size * neg.size)


@SETTINGS
@given(labelled())
def test_auroc_equals_mann_whitney_exactly(xy: tuple[F64, npt.NDArray[np.bool_]]) -> None:
    x, y = xy
    got = _auroc(x, y)
    n1, n0 = int(y.sum()), int((~y).sum())
    assert got == float(_pair_count(x, y))
    assert got == float(stats.mannwhitneyu(x[y], x[~y]).statistic) / (n1 * n0)
    assert abs(got - roc_auc_score(y, x)) <= 1e-12


def test_auroc_examples_with_tied_scores() -> None:
    n = MIN_UNIVARIATE_ROWS
    y = np.arange(n) % 2 == 0
    assert _auroc(np.full(n, 3.0), y) == 0.5  # all tied
    assert _auroc(np.where(y, -0.0, 0.0), y) == 0.5  # ±0 tie
    assert _auroc(np.where(y, 1.0, 0.0), y) == 1.0
    assert _auroc(np.where(y, 0.0, 1.0), y) == 0.0
    # Three positives at 2.0 beat all five negatives; two at 1.0 tie with all five: U = 15 + 5.
    x = np.array([2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    yy = np.array([1, 0, 1, 0, 1, 0, 1, 0, 1, 0], dtype=bool)
    assert _auroc(x, yy) == float(Fraction(3 * 5 * 2 + 2 * 5, 2 * 25))
    assert math.isnan(_auroc(np.arange(n, dtype=float), np.zeros(n, dtype=bool)))
    assert math.isnan(_auroc(np.arange(n, dtype=float), np.ones(n, dtype=bool)))
