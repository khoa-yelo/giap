"""The switch between the worker-process null and the ADR-0053 thread/serial null (ADR-0063).

``null_p`` uses ``_process_columns`` exactly when ``threads`` > 1 and B >= ``_PROCESS_MIN_B``
(200 shipped). Pins the boundary (B == threshold uses processes, threshold - 1 does not), the
``threads`` > 1 condition and the shipped constant.
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

SEED = 20261027


def _world(n: int = 60) -> tuple[pd.DataFrame, np.ndarray, list[str], np.ndarray]:
    rng = np.random.default_rng(9)
    tu = d.tu_ids(d.mixed_sizes(rng, n, 1, 3))
    X = pd.DataFrame(
        {"cai": rng.normal(size=n), "gravy": rng.normal(size=n)},
        index=pd.Index([f"G{i:05d}" for i in range(n)], name="locus_tag"),
    )
    y = X["cai"].to_numpy() + rng.normal(0.0, 1.0, n)
    folds = d.mod().tu_folds(tu, 5, 1, for_step("determinants.cv", SEED))
    return X, np.asarray(y, dtype=float), tu, np.asarray(folds[0], dtype=np.int64)


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    """Records the ``threads`` of every ``_process_columns`` call (the real one still runs)."""
    seen: list[int] = []
    real = N._process_columns

    def spy(*args: Any) -> list[list[float]]:
        seen.append(args[-1])
        return real(*args)

    monkeypatch.setattr(N, "_process_columns", spy)
    return seen


def _run(threads: int, B: int) -> None:
    X, y, tu, f1 = _world()
    cfg = d.config(fast=True, n_estimators=5)
    N.null_p(X, y, tu, f1, {"elastic_net": 0.1}, cfg, B, SEED, threads)


def test_shipped_threshold_is_200() -> None:
    assert N._PROCESS_MIN_B == 200


def test_one_thread_never_uses_processes_at_or_above_the_threshold(
    calls: list[int], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(N, "_PROCESS_MIN_B", 3)
    _run(1, 3)
    _run(1, 4)
    assert calls == []


def test_b_equal_to_the_threshold_uses_processes(
    calls: list[int], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(N, "_PROCESS_MIN_B", 4)
    _run(2, 4)
    assert calls == [2]


def test_b_one_below_the_threshold_does_not(
    calls: list[int], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(N, "_PROCESS_MIN_B", 4)
    _run(2, 3)
    assert calls == []
