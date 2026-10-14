"""TU-block bootstrap, percentile CI and circular block bootstrap (issue M4-04; GIAP_SPEC.md
§7.3.5(a)-(b), §8.2 "resample TUs, not genes", §8.3; decision 14 sets B = 1000 by default).

Conventions pinned here: randomness only from the ``np.random.Generator`` passed in (callers get
it from ``giap.rng.for_step``, as in ``giap.stats.permutation``); ``estimate = stat_fn(df)``;
each replicate draws, per stratum, as many TU IDs as the stratum has, with replacement, and
concatenates their rows (a TU drawn twice contributes its rows twice, in its internal order); NaN
replicate statistics are dropped and fewer than B/2 valid gives CI (NaN, NaN). Circular blocks
(Politis-Romano): ceil(n / block_len) start positions are drawn with ``rng.integers`` uniform on
[0, n); each block wraps modulo n; concatenated and truncated to n. For n = 0 the result is shape
(B, 0) whatever ``block_len`` is (resolved ambiguity: the n = 0 rule wins over ``block_len > n``).
"""

from __future__ import annotations

import dataclasses
import math
from collections import Counter
from collections.abc import Callable
from typing import Any

import numpy as np
import pandas as pd
import pytest

from giap.rng import for_step
from giap.stats.bootstrap import (
    BootCI,
    circular_block_bootstrap,
    percentile_ci,
    tu_block_bootstrap,
)

NAN = float("nan")


def rng(seed: int = 0) -> np.random.Generator:
    return for_step("test.bootstrap", seed)


def frame(tus: dict[str, list[float]], **extra: dict[str, Any]) -> pd.DataFrame:
    """Gene rows grouped by TU with a within-TU position column ``pos``."""
    rows = []
    for tu, vals in tus.items():
        for i, v in enumerate(vals):
            row = {"tu": tu, "pos": i, "v": float(v)}
            for col, mapping in extra.items():
                row[col] = mapping[tu]
            rows.append(row)
    return pd.DataFrame(rows)


class Recorder:
    """stat_fn wrapper that keeps a copy of every frame it is called on."""

    def __init__(self, fn: Callable[[pd.DataFrame], float]) -> None:
        self.fn = fn
        self.frames: list[pd.DataFrame] = []

    def __call__(self, df: pd.DataFrame) -> float:
        self.frames.append(df.reset_index(drop=True).copy())
        return self.fn(df)


def split_draws(resampled: pd.DataFrame, original: pd.DataFrame) -> list[str]:
    """Split a resampled frame into whole-TU draws; fail if any TU is split or reordered."""
    groups = {tu: g.reset_index(drop=True) for tu, g in original.groupby("tu", sort=False)}
    draws: list[str] = []
    i = 0
    n = len(resampled)
    while i < n:
        tu = resampled.at[i, "tu"]
        assert tu in groups, f"unknown TU {tu!r}"
        g = groups[tu]
        block = resampled.iloc[i : i + len(g)].reset_index(drop=True)
        assert len(block) == len(g), f"TU {tu!r} truncated"
        pd.testing.assert_frame_equal(block[g.columns], g, check_dtype=False)
        draws.append(str(tu))
        i += len(g)
    return draws


# ---------------------------------------------------------------------------- BootCI, percentile_ci


def test_bootci_is_frozen_slotted_dataclass_with_pinned_fields() -> None:
    assert dataclasses.is_dataclass(BootCI)
    assert [f.name for f in dataclasses.fields(BootCI)] == [
        "estimate",
        "ci_low",
        "ci_high",
        "n_valid",
        "B",
    ]
    ci = BootCI(1.0, 0.0, 2.0, 10, 10)
    with pytest.raises(dataclasses.FrozenInstanceError):
        ci.estimate = 3.0  # type: ignore[misc]
    assert not hasattr(ci, "__dict__")


def test_example_4_percentile_ci() -> None:
    assert percentile_ci(np.arange(101)) == pytest.approx((2.5, 97.5), abs=1e-12)
    lo, hi = percentile_ci(np.arange(101), level=0.90)
    assert (lo, hi) == pytest.approx((5.0, 95.0), abs=1e-12)


def test_percentile_ci_uses_linear_quantiles() -> None:
    boot = np.array([4.0, 1.0, 3.0, 2.0])
    lo, hi = percentile_ci(boot, level=0.5)
    expected = np.quantile(boot, [0.25, 0.75], method="linear")
    assert (lo, hi) == (float(expected[0]), float(expected[1]))
    assert isinstance(lo, float) and isinstance(hi, float)


@pytest.mark.parametrize("level", [0.0, 1.0, -0.1, 1.5])
def test_percentile_ci_level_outside_open_unit_interval_raises(level: float) -> None:
    with pytest.raises(ValueError):
        percentile_ci(np.arange(10), level=level)


# ---------------------------------------------------------------------------- tu_block_bootstrap

EX5 = {"A": [1.0, 1.0], "B": [0.0], "C": [1.0]}


def test_example_5_sum_over_three_tus() -> None:
    df = frame(EX5)
    rec = Recorder(lambda d: float(d["v"].sum()))
    ci = tu_block_bootstrap(df, "tu", rec, B=2000, rng=rng(5))
    assert isinstance(ci, BootCI)
    assert ci.estimate == 3.0
    assert ci.B == 2000
    assert ci.n_valid == 2000
    totals = {"A": 2.0, "B": 0.0, "C": 1.0}
    reps = [f for f in rec.frames]
    assert len(reps) >= 2000
    sums = set()
    for f in reps:
        draws = split_draws(f, df)
        assert len(draws) == 3
        s = float(f["v"].sum())
        assert s == sum(totals[t] for t in draws)
        sums.add(s)
    assert sums <= {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0}
    assert len(sums) == 7  # P(any value missing in 2000 draws) is negligible
    # sum of 3 iid draws from {2, 0, 1}: symmetric about 3, tails 1/27 each
    assert 0.0 <= ci.ci_low <= 1.0
    assert 5.0 <= ci.ci_high <= 6.0
    assert ci.ci_low <= ci.estimate <= ci.ci_high


def test_example_6_single_tu_gives_degenerate_ci() -> None:
    df = frame({"A": [1.0, 1.0]})
    rec = Recorder(lambda d: float(d["v"].sum()))
    ci = tu_block_bootstrap(df, "tu", rec, B=200, rng=rng(6))
    assert ci.estimate == ci.ci_low == ci.ci_high == 2.0
    assert ci.n_valid == 200
    for f in rec.frames:
        assert split_draws(f, df) == ["A"]


def test_example_7_strata_keep_group_sizes() -> None:
    tus = {f"i{k}": [float(k)] * (k % 3 + 1) for k in range(2)}
    tus.update({f"o{k}": [float(-k)] * (k % 2 + 1) for k in range(5)})
    inside = {t: t.startswith("i") for t in tus}
    df = frame(tus, inside=inside)
    rec = Recorder(
        lambda d: float(d.loc[d["inside"], "v"].mean() - d.loc[~d["inside"], "v"].mean())
    )
    ci = tu_block_bootstrap(df, "tu", rec, B=300, rng=rng(7), strata_col="inside")
    assert ci.n_valid == 300
    for f in rec.frames:
        draws = split_draws(f, df)
        c = Counter(inside[t] for t in draws)
        assert c[True] == 2 and c[False] == 5


def test_tu_drawn_twice_contributes_rows_twice() -> None:
    df = frame({"A": [1.0, 2.0, 3.0], "B": [10.0]})
    rec = Recorder(lambda d: float(len(d)))
    tu_block_bootstrap(df, "tu", rec, B=400, rng=rng(8))
    lengths = {len(f) for f in rec.frames}
    # draws AA -> 6 rows, AB/BA -> 4 rows, BB -> 2 rows
    assert lengths == {2, 4, 6}
    for f in rec.frames:
        assert len(split_draws(f, df)) == 2


def test_tu_ids_can_be_integers_and_rows_unsorted() -> None:
    df = pd.DataFrame({"tu": [7, 3, 7, 3, 9], "pos": [0, 0, 1, 1, 0], "v": [1.0, 2, 3, 4, 5]})
    rec = Recorder(lambda d: float(d["v"].sum()))
    ci = tu_block_bootstrap(df, "tu", rec, B=200, rng=rng(9))
    assert ci.estimate == 15.0
    totals = {7: 4.0, 3: 6.0, 9: 5.0}
    for f in rec.frames:
        per_tu = f.groupby("tu")["v"].sum()
        counts = f.groupby("tu").size()
        sizes = df.groupby("tu").size()
        for tu, cnt in counts.items():
            assert cnt % sizes[tu] == 0  # whole TUs only
            assert per_tu[tu] == totals[tu] * (cnt // sizes[tu])
        assert sum(counts[t] // sizes[t] for t in counts.index) == 3


def test_nan_replicates_are_dropped_and_counted() -> None:
    # stat is NaN whenever TU "B" is missing from the resample (prob (2/3)^3 = 0.30)
    df = frame({"A": [1.0], "B": [2.0], "C": [4.0]})

    def stat(d: pd.DataFrame) -> float:
        return float(d["v"].sum()) if (d["tu"] == "B").any() else NAN

    ci = tu_block_bootstrap(df, "tu", stat, B=1000, rng=rng(10))
    assert ci.B == 1000
    assert 600 <= ci.n_valid <= 780
    assert math.isfinite(ci.ci_low) and math.isfinite(ci.ci_high)
    assert ci.estimate == 7.0


def test_fewer_than_half_valid_gives_nan_ci() -> None:
    # NaN unless all three TUs present (prob 6/27 = 0.22 < 1/2)
    df = frame({"A": [1.0], "B": [2.0], "C": [4.0]})

    def stat(d: pd.DataFrame) -> float:
        return float(d["v"].sum()) if d["tu"].nunique() == 3 else NAN

    ci = tu_block_bootstrap(df, "tu", stat, B=1000, rng=rng(11))
    assert ci.estimate == 7.0
    assert ci.n_valid < 500
    assert math.isnan(ci.ci_low) and math.isnan(ci.ci_high)


def test_ci_matches_percentile_of_recorded_replicates() -> None:
    df = frame({f"t{k}": [float(k), float(k * k % 5)] for k in range(8)})
    rec = Recorder(lambda d: float(d["v"].mean()))
    ci = tu_block_bootstrap(df, "tu", rec, B=500, rng=rng(12), level=0.9)
    reps = np.array([float(f["v"].mean()) for f in rec.frames])
    # the estimate call is the frame equal to df; drop exactly one such call
    assert len(reps) == 501
    est_idx = [i for i, f in enumerate(rec.frames) if f.equals(df.reset_index(drop=True))]
    assert est_idx
    boot = np.delete(reps, est_idx[0])
    lo, hi = percentile_ci(boot, level=0.9)
    assert (ci.ci_low, ci.ci_high) == pytest.approx((lo, hi), abs=1e-12)


def test_determinism_same_seed_same_result_different_seed_differs() -> None:
    df = frame({f"t{k}": [float(k), float(-k) / 3] for k in range(12)})

    def stat(d: pd.DataFrame) -> float:
        return float(d["v"].sum())

    a = tu_block_bootstrap(df, "tu", stat, B=300, rng=for_step("regions", 42))
    b = tu_block_bootstrap(df, "tu", stat, B=300, rng=for_step("regions", 42))
    c = tu_block_bootstrap(df, "tu", stat, B=300, rng=for_step("regions", 43))
    assert a == b
    assert (a.ci_low, a.ci_high) != (c.ci_low, c.ci_high)


def test_does_not_mutate_input_frame() -> None:
    df = frame(EX5)
    before = df.copy()
    tu_block_bootstrap(df, "tu", lambda d: float(d["v"].sum()), B=50, rng=rng(13))
    pd.testing.assert_frame_equal(df, before)


def test_na_tu_id_raises() -> None:
    df = pd.DataFrame({"tu": ["A", None, "B"], "v": [1.0, 2.0, 3.0]})
    with pytest.raises(ValueError):
        tu_block_bootstrap(df, "tu", lambda d: float(d["v"].sum()), B=10, rng=rng(14))
    df2 = pd.DataFrame({"tu": pd.array(["A", pd.NA], dtype="string"), "v": [1.0, 2.0]})
    with pytest.raises(ValueError):
        tu_block_bootstrap(df2, "tu", lambda d: float(d["v"].sum()), B=10, rng=rng(14))


def test_stat_fn_errors_propagate() -> None:
    calls = {"n": 0}

    def stat(d: pd.DataFrame) -> float:
        calls["n"] += 1
        if calls["n"] > 3:
            raise RuntimeError("boom")
        return float(d["v"].sum())

    with pytest.raises(RuntimeError, match="boom"):
        tu_block_bootstrap(frame(EX5), "tu", stat, B=20, rng=rng(15))


@pytest.mark.parametrize("level", [0.0, 1.0, 1.2])
def test_tu_block_level_validation(level: float) -> None:
    with pytest.raises(ValueError):
        tu_block_bootstrap(
            frame(EX5), "tu", lambda d: float(d["v"].sum()), B=10, rng=rng(16), level=level
        )


# ---------------------------------------------------------------------------- circular blocks


class ForcedStarts:
    """Stand-in Generator whose ``integers`` returns the given start positions in order."""

    def __init__(self, starts: list[int]) -> None:
        self.starts = list(starts)
        self.pos = 0

    def integers(self, low: int, high: int | None = None, size: Any = None, **_: Any) -> Any:
        hi = low if high is None else high
        lo = 0 if high is None else low
        if size is None:
            v = self.starts[self.pos]
            self.pos += 1
            assert lo <= v < hi
            return np.int64(v)
        k = int(np.prod(size))
        vals = np.asarray(self.starts[self.pos : self.pos + k], dtype=np.int64)
        assert vals.size == k, "unexpected number of start draws"
        self.pos += k
        assert ((vals >= lo) & (vals < hi)).all()
        return vals.reshape(size)


def test_example_8_forced_starts_wrap_and_truncate() -> None:
    fake = ForcedStarts([8, 0, 4, 9])
    out = circular_block_bootstrap(np.arange(10), 3, B=1, rng=fake)  # type: ignore[arg-type]
    assert out.shape == (1, 10)
    assert out.dtype == np.float64
    np.testing.assert_array_equal(out[0], [8, 9, 0, 0, 1, 2, 4, 5, 6, 9])
    assert fake.pos == 4  # ceil(10 / 3) starts per replicate


def test_forced_starts_two_replicates_exact_division() -> None:
    fake = ForcedStarts([3, 0, 2, 1])
    out = circular_block_bootstrap(np.array([10.0, 11, 12, 13]), 2, B=2, rng=fake)  # type: ignore[arg-type]
    np.testing.assert_array_equal(out, [[13, 10, 10, 11], [12, 13, 11, 12]])


def test_example_9_block_len_one_is_iid_resampling() -> None:
    x = np.array([0.5, 1.5, 2.5, 3.5, 4.5])
    out = circular_block_bootstrap(x, 1, B=4000, rng=rng(9))
    assert out.shape == (4000, 5)
    assert np.isin(out, x).all()
    # iid: each value has marginal probability 1/5; some rows repeat a value
    freq = np.array([(out == v).mean() for v in x])
    np.testing.assert_allclose(freq, 0.2, atol=0.01)
    assert (np.array([len(set(r)) for r in out]) < 5).mean() > 0.9  # 1 - 5!/5^5 = 0.96


@pytest.mark.parametrize(("n", "L"), [(10, 3), (12, 4), (7, 7), (9, 2), (5, 5)])
def test_rows_are_runs_of_consecutive_values_mod_n(n: int, L: int) -> None:
    out = circular_block_bootstrap(np.arange(n), L, B=300, rng=rng(n * 100 + L))
    assert out.shape == (300, n)
    idx = out.astype(int)
    for row in idx:
        for s in range(0, n, L):
            block = row[s : s + L]
            assert ((np.diff(block) % n) == 1).all(), f"block {block} not consecutive mod {n}"


def test_marginals_are_uniform_under_circular_scheme() -> None:
    n = 11
    out = circular_block_bootstrap(np.arange(n), 4, B=5000, rng=rng(21))
    freq = np.bincount(out.astype(int).ravel(), minlength=n) / out.size
    np.testing.assert_allclose(freq, 1 / n, atol=0.012)


def test_tiny_replicons() -> None:
    one = circular_block_bootstrap(np.array([2.0]), 1, B=5, rng=rng(1))
    np.testing.assert_array_equal(one, np.full((5, 1), 2.0))
    two = circular_block_bootstrap(np.array([1.0, 2.0]), 2, B=200, rng=rng(2))
    rows = {tuple(r) for r in two}
    assert rows == {(1.0, 2.0), (2.0, 1.0)}


def test_empty_input_gives_b_by_zero() -> None:
    out = circular_block_bootstrap(np.array([], dtype=float), 1, B=7, rng=rng(3))
    assert out.shape == (7, 0)
    assert out.dtype == np.float64


@pytest.mark.parametrize("L", [0, -1, 6])
def test_block_len_out_of_range_raises(L: int) -> None:
    with pytest.raises(ValueError):
        circular_block_bootstrap(np.arange(5), L, B=3, rng=rng(4))


def test_circular_determinism() -> None:
    x = np.linspace(-1.0, 1.0, 50)
    a = circular_block_bootstrap(x, 5, B=100, rng=for_step("islands", 1))
    b = circular_block_bootstrap(x, 5, B=100, rng=for_step("islands", 1))
    c = circular_block_bootstrap(x, 5, B=100, rng=for_step("islands", 2))
    np.testing.assert_array_equal(a, b)
    assert not np.array_equal(a, c)


def test_circular_does_not_mutate_input() -> None:
    x = np.arange(6, dtype=float)
    circular_block_bootstrap(x, 2, B=10, rng=rng(5))
    np.testing.assert_array_equal(x, np.arange(6))
