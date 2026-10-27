"""M8-05 lead decision 7 (statistician pre-code review): degenerate redundancy re-tests.

- the re-test calls ``prepare`` with min_size 1, so a 1-TU remainder is still tested;
- a 0-TU remainder (the parent is exactly its significant descendants) is not tested:
  ``redundant = True`` and ``size_filtered``;
- a NaN re-test p gives ``redundant`` NA; the main-run p, q and NES are unchanged."""

# ruff: noqa: F811
from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd
from helpers.gsea_data import gene_sets
from test_gsea import ALPHA, FakeEngine, _family, _flags, _go_frame, _row, _run, fake  # noqa: F401

from gippa.stats.multitest import bh


def test_zero_tu_retest_is_redundant_and_size_filtered(fake: Callable[..., FakeEngine]) -> None:
    genes, tags = _go_frame()
    sets = {"GO:C": list(tags[:15]), "GO:P": list(tags[:15]), "GO:N": list(tags[-40:])}
    parents = {"GO:C": {"GO:P"}, "GO:P": set(), "GO:N": set()}
    eng = fake({15: 0.001, 40: 0.5})
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)])
    p = _row(out, "GO:P")
    assert p["q_value"] <= ALPHA
    assert p["redundant"] == True  # noqa: E712
    assert "size_filtered" in _flags(p["flags"])
    assert all(s != [0] and s for s, _ in eng.calls), eng.calls  # no fgsea call on 0 TUs


def test_nan_retest_leaves_redundant_na_and_main_stats(fake: Callable[..., FakeEngine]) -> None:
    genes, tags = _go_frame()
    sets, parents = _family(tags, 15, 12, [40, 41, 42])
    fake({15: 0.001, 27: 0.002, 12: float("nan"), 40: 0.5, 41: 0.6, 42: 0.7})
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)])
    p = _row(out, "GO:P")
    assert pd.isna(p["redundant"]), p["redundant"]
    assert "unestimated" not in _flags(p["flags"])
    assert p["p_value"] == 0.002 and np.isfinite(p["nes"])
    np.testing.assert_allclose(p["q_value"], bh([0.001, 0.002, 0.5, 0.6, 0.7])[1], atol=1e-15)


def test_one_tu_retest_is_still_run(fake: Callable[..., FakeEngine]) -> None:
    genes, tags = _go_frame()
    sets, parents = _family(tags, 15, 1, [40, 41, 42])
    eng = fake({15: 0.001, 16: 0.002, 1: 0.3, 40: 0.5, 41: 0.6, 42: 0.7})
    out = _run(genes, [gene_sets("go_bp", sets, parents=parents)])
    assert eng.calls[-1][0] == [1]
    assert _row(out, "GO:P")["redundant"] == True  # noqa: E712
