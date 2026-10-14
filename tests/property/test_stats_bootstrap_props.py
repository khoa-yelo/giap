"""Hypothesis properties of the bootstrap and Hodges-Lehmann functions (issue M4-04; §8.2, §8.3).

- ``hodges_lehmann``: equals the direct median of pairwise differences; shift-equivariant
  (HL(x + c, y) = HL(x, y) + c, HL(x, y - c) = HL(x, y) + c); antisymmetric
  (HL(x, y) = -HL(y, x)); invariant to input order; scale-equivariant for c > 0.
- ``tu_block_bootstrap``: replicates contain only whole TUs (each drawn TU's rows intact and in
  their internal order) and exactly as many TU draws as TUs; ci_low <= ci_high; for a sum
  statistic the 95 % CI contains the estimate; n_valid == B for a finite statistic.
- ``circular_block_bootstrap``: shape (B, n), values drawn from x, each block of ``block_len``
  (the last possibly shorter) is a run of consecutive positions modulo n.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from hypothesis import assume, given, settings
from hypothesis import strategies as st

from giap.rng import for_step
from giap.stats.bootstrap import circular_block_bootstrap, tu_block_bootstrap
from giap.stats.robust import hodges_lehmann

FAST = settings(max_examples=100, deadline=None)

INTS = st.lists(st.integers(-50, 50).map(float), min_size=1, max_size=25)
FLOATS = st.lists(
    st.floats(-1e3, 1e3, allow_nan=False, allow_infinity=False), min_size=1, max_size=25
)


def direct(x: list[float], y: list[float]) -> float:
    return float(np.median(np.subtract.outer(np.asarray(x), np.asarray(y))))


# ---------------------------------------------------------------------------- hodges_lehmann


@FAST
@given(FLOATS, FLOATS)
def test_hl_equals_direct_median(x: list[float], y: list[float]) -> None:
    assert hodges_lehmann(x, y) == direct(x, y)


@FAST
@given(INTS, INTS, st.integers(-1000, 1000).map(float))
def test_hl_shift_equivariant_exact_on_integers(x: list[float], y: list[float], c: float) -> None:
    base = hodges_lehmann(x, y)
    assert hodges_lehmann([v + c for v in x], y) == base + c
    assert hodges_lehmann(x, [v - c for v in y]) == base + c


@FAST
@given(FLOATS, FLOATS, st.floats(-1e3, 1e3, allow_nan=False))
def test_hl_shift_equivariant_floats(x: list[float], y: list[float], c: float) -> None:
    base = hodges_lehmann(x, y)
    shifted = hodges_lehmann([v + c for v in x], y)
    assert abs(shifted - (base + c)) <= 1e-9 * max(
        1.0, abs(base), abs(c), *map(abs, x), *map(abs, y)
    )


@FAST
@given(FLOATS, FLOATS)
def test_hl_antisymmetric(x: list[float], y: list[float]) -> None:
    assert hodges_lehmann(x, y) == -hodges_lehmann(y, x)


@FAST
@given(INTS, INTS, st.randoms(use_true_random=False))
def test_hl_order_invariant(x: list[float], y: list[float], r: object) -> None:
    xs, ys = list(x), list(y)
    r.shuffle(xs)  # type: ignore[attr-defined]
    r.shuffle(ys)  # type: ignore[attr-defined]
    assert hodges_lehmann(xs, ys) == hodges_lehmann(x, y)


@FAST
@given(INTS, INTS, st.sampled_from([0.5, 2.0, 4.0]))
def test_hl_scale_equivariant(x: list[float], y: list[float], k: float) -> None:
    assert hodges_lehmann([k * v for v in x], [k * v for v in y]) == k * hodges_lehmann(x, y)


@FAST
@given(INTS, INTS)
def test_hl_between_min_and_max_difference(x: list[float], y: list[float]) -> None:
    d = np.subtract.outer(np.asarray(x), np.asarray(y))
    assert d.min() <= hodges_lehmann(x, y) <= d.max()


# ---------------------------------------------------------------------------- tu_block_bootstrap

TU_SIZES = st.lists(st.integers(1, 4), min_size=1, max_size=8)


@st.composite
def tu_frames(draw: st.DrawFn) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    sizes = draw(TU_SIZES)
    rows = []
    for t, k in enumerate(sizes):
        for i in range(k):
            v = draw(st.integers(-20, 20))
            rows.append({"tu": f"T{t}", "pos": i, "v": float(v)})
    # interleave the TUs' rows in the input; resampled TUs keep their internal order anyway
    perm = draw(st.permutations(range(len(rows))))
    df = pd.DataFrame([rows[j] for j in perm]).reset_index(drop=True)
    df["pos"] = df.groupby("tu").cumcount()  # internal order = order of appearance in the input
    groups = {
        tu: g.sort_values("pos").reset_index(drop=True) for tu, g in df.groupby("tu", sort=False)
    }
    return df, groups


def _draws(resampled: pd.DataFrame, groups: dict[str, pd.DataFrame]) -> int:
    """Number of whole-TU draws; asserts each TU's rows appear as whole copies of the TU."""
    r = resampled.reset_index(drop=True)
    total = 0
    for tu, cnt in r.groupby("tu").size().items():
        g = groups[str(tu)]
        assert cnt % len(g) == 0, f"TU {tu} split"
        copies = cnt // len(g)
        sub = r[r["tu"] == tu].sort_values("pos", kind="stable").reset_index(drop=True)
        expected = pd.concat([g] * copies).sort_values("pos", kind="stable").reset_index(drop=True)
        pd.testing.assert_frame_equal(sub[g.columns], expected, check_dtype=False)
        total += copies
    return total


@FAST
@given(tu_frames(), st.integers(0, 2**31 - 1))
def test_tu_bootstrap_replicates_are_whole_tus(
    data: tuple[pd.DataFrame, dict[str, pd.DataFrame]], seed: int
) -> None:
    df, groups = data
    seen: list[pd.DataFrame] = []

    def stat(d: pd.DataFrame) -> float:
        seen.append(d.copy())
        return float(d["v"].sum())

    ci = tu_block_bootstrap(df, "tu", stat, B=30, rng=for_step("prop.tu", seed))
    assert ci.n_valid == 30 and ci.B == 30
    assert ci.estimate == float(df["v"].sum())
    assert ci.ci_low <= ci.ci_high
    for f in seen:
        assert _draws(f, groups) == len(groups)
        # within each TU copy the internal order is preserved: pos runs 0..k-1 in sequence
        for tu, g in f.reset_index(drop=True).groupby("tu", sort=False):
            k = len(groups[str(tu)])
            pos = g["pos"].to_numpy()
            assert (pos == np.tile(np.arange(k), len(pos) // k)).all()


@FAST
@given(
    st.lists(st.floats(-100, 100, allow_nan=False), min_size=2, max_size=15),
    st.integers(0, 2**31 - 1),
)
def test_tu_bootstrap_sum_ci_contains_estimate(vals: list[float], seed: int) -> None:
    assume(len(set(vals)) > 1)
    df = pd.DataFrame({"tu": [f"T{i}" for i in range(len(vals))], "v": vals})
    ci = tu_block_bootstrap(
        df, "tu", lambda d: float(d["v"].sum()), B=400, rng=for_step("prop.tu.ci", seed)
    )
    assert ci.ci_low <= ci.estimate + 1e-9 and ci.estimate - 1e-9 <= ci.ci_high


# ---------------------------------------------------------------------------- circular blocks


@FAST
@given(st.integers(1, 40), st.data(), st.integers(1, 20), st.integers(0, 2**31 - 1))
def test_circular_rows_are_consecutive_runs(n: int, data: st.DataObject, B: int, seed: int) -> None:
    L = data.draw(st.integers(1, n))
    x = np.arange(n, dtype=float) * 1.5 - 3.0
    out = circular_block_bootstrap(x, L, B, for_step("prop.cbb", seed))
    assert out.shape == (B, n)
    assert np.isin(out, x).all()
    idx = np.rint((out + 3.0) / 1.5).astype(int)
    for row in idx:
        for s in range(0, n, L):
            assert ((np.diff(row[s : s + L]) % n) == 1).all()
