"""``_ElasticNetCV`` (inner paths with ``check_input=False``) equals ``ElasticNetCV`` exactly
(PR #286 review, B1.3 and N1).

The speed-up relies on sklearn's private contract that ``_path_residuals`` validates X, y, Xy and
the Gram matrix (Fortran X, C-ordered float64 Gram) before calling ``path``, so the checks it
skips are no-ops. Verified on scikit-learn 1.9.1 (uv.lock), scipy 1.18.1, numpy 2.5.4. A
sklearn release that breaks the contract fails here: either loudly (a Cython layout/dtype error)
or as a changed path.

Both branches run: n > p (``precompute="auto"`` builds the Gram matrix) and n <= p (no Gram).
The keyword arguments are those of ``gippa.determinants.models._fit_predict``.
"""

from __future__ import annotations

import numpy as np
import pytest
from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import GroupKFold

from gippa.determinants.models import _ElasticNetCV


def _kw(cv: list[tuple[np.ndarray, np.ndarray]]) -> dict[str, object]:
    return {
        "l1_ratio": [0.1, 0.5, 0.9],
        "cv": cv,
        "alphas": 100,
        "max_iter": 10000,
        "selection": "cyclic",
        "random_state": 0,
    }


@pytest.mark.parametrize(
    ("n", "p", "seed"),
    [
        (200, 12, 1),  # n > p: Gram branch
        (120, 30, 2),  # n > p: Gram branch, more columns
        (60, 60, 3),  # n == p: no Gram
        (60, 150, 4),  # n < p: no Gram
    ],
)
def test_prechecked_elastic_net_equals_sklearn(n: int, p: int, seed: int) -> None:
    r = np.random.default_rng(seed)
    X = r.normal(size=(n, p))
    X[:, 1] = np.round(X[:, 1])  # a tied column, as count-like features are
    y = X[:, 0] - 0.5 * X[:, 1] + r.normal(size=n)
    groups = r.integers(0, 20, n)
    cv = list(GroupKFold(3).split(X, y, groups))

    want = ElasticNetCV(**_kw(cv)).fit(X, y)
    got = _ElasticNetCV(**_kw(cv)).fit(X, y)

    assert np.array_equal(got.coef_, want.coef_)
    assert got.intercept_ == want.intercept_
    assert np.array_equal(got.mse_path_, want.mse_path_)
    assert np.array_equal(got.alphas_, want.alphas_)
    assert got.alpha_ == want.alpha_
    assert got.l1_ratio_ == want.l1_ratio_
    assert got.n_iter_ == want.n_iter_
    # The model is not degenerate (some coefficient is non-zero), so equality is informative.
    assert np.any(want.coef_ != 0)


def test_prechecked_path_is_the_one_used() -> None:
    """The override is the ``path`` staticmethod; the class is still an ``ElasticNetCV``."""
    assert issubclass(_ElasticNetCV, ElasticNetCV)
    assert _ElasticNetCV.path is not ElasticNetCV.path
