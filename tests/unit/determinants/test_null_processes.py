"""The elastic-net permutation null on worker processes (perf; ADR-0053 follow-up, ADR-0063).

With ``threads`` > 1 the null refits of both models run on one pool of spawned worker processes,
one BLAS/OpenMP thread each, collected in b order: the null R² arrays, p and b are bit-identical
to the serial path (``threads`` = 1) for any worker count.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_models_data as d

from giap.determinants import null as N
from giap.rng import for_step

SEED = 20261026


def _world(n: int = 90) -> tuple[pd.DataFrame, np.ndarray, list[str], np.ndarray]:
    rng = np.random.default_rng(5)
    tu = d.tu_ids(d.mixed_sizes(rng, n, 1, 3))
    X = pd.DataFrame(
        {"cai": rng.normal(size=n), "gravy": rng.normal(size=n)},
        index=pd.Index([f"G{i:05d}" for i in range(n)], name="locus_tag"),
    )
    y = X["cai"].to_numpy() + rng.normal(0.0, 1.0, n)
    folds = d.mod().tu_folds(tu, 5, 1, for_step("determinants.cv", SEED))
    return X, np.asarray(y, dtype=float), tu, np.asarray(folds[0], dtype=np.int64)


def _run(threads: int, models: tuple[str, ...] = ("elastic_net",), B: int = 11) -> dict[str, Any]:
    X, y, tu, f1 = _world()
    cfg = d.config(fast=True, n_estimators=10)
    return dict(N.null_p(X, y, tu, f1, {m: 0.1 for m in models}, cfg, B, SEED, threads))


def _key(res: dict[str, Any]) -> dict[str, tuple[bytes, float, int, int, str]]:
    return {m: (r.null_r2.tobytes(), r.p, r.b_count, r.n_perm, r.flags) for m, r in res.items()}


@pytest.mark.parametrize("threads", [2, 3])
def test_process_null_is_bit_identical_to_serial(threads: int) -> None:
    assert _key(_run(threads)) == _key(_run(1))


def test_both_models_bit_identical_across_threads() -> None:
    both = ("lightgbm", "elastic_net")
    assert _key(_run(3, both, B=7)) == _key(_run(1, both, B=7))


def test_threads_above_one_use_the_process_pool(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[int] = []
    real = N._process_columns

    def spy(*args: Any) -> list[list[float]]:
        calls.append(args[-1])
        return real(*args)

    monkeypatch.setattr(N, "_process_columns", spy)
    _run(2, B=5)
    assert calls == [2]
    calls.clear()
    _run(1, B=5)
    assert calls == []
