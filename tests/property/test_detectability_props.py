"""Hypothesis properties of the M4-13 detectability model (docs/specs/M4-13-detectability.md).

- P1 ``tu_folds``: a partition of the distinct TU IDs into k folds whose sizes differ by ≤ 1,
  equal to the oracle, and independent of input order / duplicates.
- P2 ``calibration_diagnostics``: equal-count bins (sizes differ by ≤ 1, sum n), bin ``mean_p``
  non-decreasing, Brier = mean (p - d)^2 in [0, 1], counts right; with distinct p the overall
  metrics do not depend on the input order.
- P3 (acceptance criterion): permuting the gene rows of all inputs leaves ``p_detect`` (per
  locus_tag), the folds and the chosen C unchanged (folds depend on TU IDs only).
- P4 eligible ``p_detect`` lies in [0, 1] (open interval when fitted).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import detectability_data as dd
import detectability_oracle as O

from giap.config.model import DetectabilityConfig
from giap.rng import for_step


def _mod() -> Any:
    from giap.genetable import detectability

    return detectability


tu_lists = st.lists(st.integers(1, 60).map(lambda i: f"TU_{i}"), min_size=1, max_size=80)


@settings(max_examples=100, deadline=None)
@given(
    tus=tu_lists, k=st.integers(1, 7), seed=st.integers(0, 2**32 - 1), perm_seed=st.integers(0, 99)
)
def test_tu_folds_partition(tus: list[str], k: int, seed: int, perm_seed: int) -> None:
    f = _mod().tu_folds(tus, k, np.random.default_rng(seed))
    assert set(f) == set(tus)
    assert all(0 <= v < k for v in f.values())
    sizes = np.bincount(list(f.values()), minlength=k)
    assert sizes.max() - sizes.min() <= 1
    assert f == O.tu_folds(tus, k, np.random.default_rng(seed))
    shuffled = list(np.random.default_rng(perm_seed).permutation(tus)) + tus[:3]
    assert _mod().tu_folds(shuffled, k, np.random.default_rng(seed)) == f


def _values(diag: pd.DataFrame, group_type: str, metric: str) -> pd.Series:
    sel = diag[(diag["group_type"] == group_type) & (diag["metric"] == metric)]
    return sel.set_index(sel["group"].astype(str))["value"].astype(float)


@settings(max_examples=100, deadline=None)
@example(n=20, n_bins=1, seed=81781)  # completely separated: no recalibration MLE (decision 17)
@given(
    n=st.integers(20, 300),
    n_bins=st.integers(1, 10),
    seed=st.integers(0, 2**32 - 1),
)
def test_calibration_diagnostics_bins(n: int, n_bins: int, seed: int) -> None:
    r = np.random.default_rng(seed)
    p = r.random(n)
    d = (r.random(n) < p).astype(int)
    diag = _mod().calibration_diagnostics(d, p, n_bins=n_bins)
    sizes = _values(diag, "bin", "n")
    assert sorted(int(g) for g in sizes.index) == list(range(1, n_bins + 1))
    sizes = sizes.loc[[str(i) for i in range(1, n_bins + 1)]].to_numpy()
    assert sizes.sum() == n and sizes.max() - sizes.min() <= 1
    mp = _values(diag, "bin", "mean_p").loc[[str(i) for i in range(1, n_bins + 1)]].to_numpy()
    assert np.all(np.diff(mp) >= -1e-15)
    brier = _values(diag, "overall", "brier").iloc[0]
    assert abs(brier - np.mean((p - d) ** 2)) < 1e-12
    assert 0.0 <= brier <= 1.0
    assert _values(diag, "overall", "n_eligible").iloc[0] == n
    assert _values(diag, "overall", "n_detected").iloc[0] == d.sum()
    perm = r.permutation(n)
    diag2 = _mod().calibration_diagnostics(d[perm], p[perm], n_bins=n_bins)
    metrics = ["brier"]
    lp = np.log(p / (1 - p))
    separated = 0 < d.sum() < n and (
        lp[d == 1].min() >= lp[d == 0].max() or lp[d == 1].max() <= lp[d == 0].min()
    )
    if 0 < d.sum() < n:  # the recalibration MLE needs both classes and no separation (decision 17)
        metrics += ["calibration_intercept", "calibration_slope"]
    for m in metrics:
        a, b = _values(diag, "overall", m).iloc[0], _values(diag2, "overall", m).iloc[0]
        if separated and m != "brier":
            assert np.isnan(a) and np.isnan(b), m
            continue
        assert abs(a - b) <= 1e-6 * max(1.0, abs(a)), m
    assert _values(diag, "overall", "recalibration_separated").iloc[0] == int(separated)


@settings(max_examples=6, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(data_seed=st.integers(0, 10_000), perm_seed=st.integers(0, 10_000), seed=st.integers(0, 50))
def test_row_order_invariance(data_seed: int, perm_seed: int, seed: int) -> None:
    gs = dd.synthetic(120, seed=data_seed)
    cfg = DetectabilityConfig()
    a = _mod().fit_detectability(
        gs.detected, gs.cov, gs.tu_of_gene, cfg, for_step("detectability", seed)
    )
    perm = np.random.default_rng(perm_seed).permutation(len(gs.detected))
    tags = gs.detected.index[perm]
    b = _mod().fit_detectability(
        gs.detected.loc[tags],
        gs.cov.loc[np.random.default_rng(perm_seed + 1).permutation(gs.cov.index)],
        gs.tu_of_gene.loc[tags[::-1]],
        cfg,
        for_step("detectability", seed),
    )
    assert a.fold_of_tu == b.fold_of_tu
    assert a.chosen_C == b.chosen_C
    assert a.constant == b.constant
    pa = a.p_detect.to_numpy(float)
    pb = b.p_detect.loc[a.p_detect.index].to_numpy(float)
    np.testing.assert_allclose(pb, pa, rtol=1e-7, atol=1e-9)
    assert np.all((pa >= 0) & (pa <= 1))
    if not a.constant:
        assert np.all((pa > 0) & (pa < 1))


def test_reverse_separation_has_no_recalibration_fit() -> None:
    """Decision 17 covers both directions: every detected logit below every undetected one."""
    p = np.array([0.05, 0.1, 0.15, 0.2, 0.6, 0.7, 0.8, 0.9])
    d = np.array([1, 1, 1, 1, 0, 0, 0, 0])
    diag = _mod().calibration_diagnostics(d, p, n_bins=2)
    for m in ("calibration_intercept", "calibration_slope"):
        assert np.isnan(_values(diag, "overall", m).iloc[0]), m
    assert _values(diag, "overall", "recalibration_separated").iloc[0] == 1
