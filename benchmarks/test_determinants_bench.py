"""Runtime of the determinants univariate bootstrap and permutation null (GIAP_SPEC.md §7.8.2
steps 1 and 4; docs/specs/M8-08-determinants-models.md, M8-09-determinants-attribution-null.md).

Inputs are the spec example 1 world of ``tests/helpers/determinants_models_data.py`` (400 genes in
TUs of 1-5 genes, the six planted columns, signal 2), built once and not timed:

- ``univariate`` with the shipped ``stats.bootstrap`` = 1,000 (rho and AUROC per encoded column);
- ``null_p`` with 50 permutations of both models (shipped ``DeterminantsConfig``) on 4 threads.

The measured times belong under ``determinants_univariate_b1000_n400`` and
``determinants_null_b50_n400`` in ``benchmarks/baseline.json``. The budgets are about 3 x the
measured times.

Run with ``pytest benchmarks/test_determinants_bench.py --benchmark-only`` (needs ``tests/helpers``,
added below).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from ._unit_helpers import timed

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "tests" / "helpers"))
import determinants_models_data as d

SEED = 20261023
UNIVARIATE_BUDGET_S = 5.0
NULL_BUDGET_S = 45.0
NULL_B = 50
THREADS = 4


@pytest.fixture(scope="module")
def world() -> tuple[Any, ...]:
    rng = np.random.default_rng(SEED)
    X, y, tu, blocks = d.planted(rng)
    detected = (rng.random(len(y)) > 0.2).tolist()  # a fifth undetected: AUROC rows estimable
    feats = d.features({c: X[c].tolist() for c in X.columns}, tu=tu, abundance=y.tolist(),
                       detected=detected, tags=list(X.index))  # fmt: skip
    return X, y, tu, blocks, feats


def test_determinants_univariate_b1000_n400_budget(benchmark: Any, world: tuple[Any, ...]) -> None:
    from giap.determinants.models import univariate
    from giap.rng import for_step

    feats = world[4]

    def run() -> float:
        secs, tbl = timed(lambda: univariate(feats, 1000, for_step("determinants.bootstrap", SEED)))
        assert int(tbl["value"].notna().sum()) >= 12, tbl
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs <= UNIVARIATE_BUDGET_S, (
        f"univariate took {secs:.1f} s (budget {UNIVARIATE_BUDGET_S})"
    )


def test_determinants_null_b50_n400_budget(benchmark: Any, world: tuple[Any, ...]) -> None:
    from giap.config.model import DeterminantsConfig
    from giap.determinants.null import null_p

    X, y, tu, _, _ = world
    cfg = DeterminantsConfig()
    folds = d.folds_for(tu, cfg, SEED)
    observed = {"lightgbm": 0.5, "elastic_net": 0.5}

    def run() -> float:
        secs, res = timed(lambda: null_p(X, y, tu, folds[0], observed, cfg, NULL_B, SEED, THREADS))
        assert all(r.null_r2.shape == (NULL_B,) for r in res.values()), res
        return secs

    secs = benchmark.pedantic(run, rounds=1, iterations=1)
    benchmark.extra_info["secs"] = secs
    assert secs <= NULL_BUDGET_S, f"null_p took {secs:.1f} s (budget {NULL_BUDGET_S})"
