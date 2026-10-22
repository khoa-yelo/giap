"""M8-08 lead decision 4: ``cross_validate(..., tu=...)`` groups the elastic net's inner
``GroupKFold`` by TU, so no TU is split between an inner training and validation fold; without
``tu`` every row is its own group."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import determinants_models_data as d
from determinants_models_data import mod


def _data() -> tuple[pd.DataFrame, np.ndarray, list[str], np.ndarray]:
    rng = np.random.default_rng(11)
    n = 90
    tu = [f"TU_{i // 3:03d}" for i in range(n)]  # 30 TUs of 3 genes, rows in TU order
    X = pd.DataFrame(
        {"cai": rng.normal(size=n), "gravy": rng.normal(size=n)},
        index=pd.Index([f"G{i:05d}" for i in range(n)]),
    )
    y = X["cai"].to_numpy() + rng.normal(0, 0.5, n)
    fold = np.repeat(np.arange(5), n // 5)  # contiguous blocks of 18 rows: TUs never split
    return X, y, tu, np.vstack([fold, fold])


def _record(monkeypatch: pytest.MonkeyPatch) -> list[np.ndarray]:
    m = mod()
    seen: list[np.ndarray] = []
    real = m.GroupKFold

    class Spy(real):  # type: ignore[misc, valid-type]
        def split(self, X: Any, y: Any = None, groups: Any = None) -> Any:
            g = np.asarray(groups)
            seen.append(g)
            for tr, te in super().split(X, y, groups):
                assert not set(g[tr]) & set(g[te]), "a group is split across inner folds"
                yield tr, te

    monkeypatch.setattr(m, "GroupKFold", Spy)
    return seen


def test_tu_groups_the_inner_folds(monkeypatch: pytest.MonkeyPatch) -> None:
    X, y, tu, folds = _data()
    seen = _record(monkeypatch)
    mod().cross_validate(X, y, folds, {"cai": "codon", "gravy": "protein"}, d.config(), 3, 1, tu=tu)
    assert seen, "the elastic net never ran its inner GroupKFold"
    for g in seen:
        assert set(g) <= set(tu)
        assert len(set(g)) < len(g)  # TU labels, not one group per row


def test_without_tu_every_row_is_its_own_group(monkeypatch: pytest.MonkeyPatch) -> None:
    X, y, _, folds = _data()
    seen = _record(monkeypatch)
    mod().cross_validate(X, y, folds, {"cai": "codon", "gravy": "protein"}, d.config(), 3, 1)
    assert seen
    for g in seen:
        assert len(set(g.tolist())) == len(g)
