"""The ``columns=`` path of ``tu_block_bootstrap`` equals its DataFrame path (PR #286 review, B1.1).

``tu_block_bootstrap(df, tu, f, B, rng, columns=cols)`` calls ``f(*arrays)`` on the same rows,
in the same order, as the DataFrame path calls ``g(frame)`` with ``g(d) = f(*(d[c] for c in
cols))``. The statistics here depend on row order (a position-weighted sum, the first row), so
a replicate that drops, adds, reorders or re-draws rows on either path changes the ``BootCI``.
The generator's final state is compared too, which pins the number of draws.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import numpy.typing as npt
import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from giap.stats.bootstrap import BootCI, tu_block_bootstrap

SETTINGS = settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])

F64 = npt.NDArray[np.float64]


def _weighted(x: F64) -> float:
    """Order-sensitive: ``sum(x * arange(n))``."""
    return float(np.sum(x * np.arange(x.size)))


def _weighted_pair(x: F64, y: F64) -> float:
    """Order-sensitive in both columns, and in their pairing."""
    n = np.arange(x.size)
    return float(np.sum(x * n) - 3.0 * np.sum(y * n * n) + x[0] - y[-1])


@st.composite
def frames(draw: st.DrawFn) -> pd.DataFrame:
    n = draw(st.integers(min_value=1, max_value=40))
    n_tu = draw(st.integers(min_value=1, max_value=n))
    seed = draw(st.integers(min_value=0, max_value=2**32 - 1))
    r = np.random.default_rng(seed)
    tu = r.integers(0, n_tu, n)
    return pd.DataFrame(
        {
            "tu": [f"tu{t}" for t in tu],
            "x": np.round(r.normal(size=n), 2),
            "y": r.integers(-3, 4, n).astype(float),
            # A TU's stratum is that of its first row; rows may disagree.
            "st": r.integers(0, draw(st.integers(min_value=1, max_value=4)), n),
        }
    )


def _both(
    df: pd.DataFrame,
    fn: Callable[..., float],
    cols: tuple[str, ...],
    B: int,
    seed: int,
    strata: str | None,
) -> tuple[BootCI, BootCI, dict[str, object], dict[str, object]]:
    r_df = np.random.default_rng(seed)
    r_col = np.random.default_rng(seed)

    def on_frame(d: pd.DataFrame) -> float:
        return fn(*(d[c].to_numpy() for c in cols))

    a = tu_block_bootstrap(df, "tu", on_frame, B, r_df, strata_col=strata)
    b = tu_block_bootstrap(df, "tu", fn, B, r_col, strata_col=strata, columns=cols)
    return a, b, r_df.bit_generator.state, r_col.bit_generator.state


@SETTINGS
@given(
    df=frames(),
    B=st.integers(min_value=1, max_value=30),
    seed=st.integers(min_value=0, max_value=2**32 - 1),
    strata=st.sampled_from([None, "st"]),
)
def test_columns_path_equals_frame_path_one_column(
    df: pd.DataFrame, B: int, seed: int, strata: str | None
) -> None:
    a, b, sa, sb = _both(df, _weighted, ("x",), B, seed, strata)
    assert a == b
    assert sa == sb


@SETTINGS
@given(
    df=frames(),
    B=st.integers(min_value=1, max_value=30),
    seed=st.integers(min_value=0, max_value=2**32 - 1),
    strata=st.sampled_from([None, "st"]),
)
def test_columns_path_equals_frame_path_two_columns(
    df: pd.DataFrame, B: int, seed: int, strata: str | None
) -> None:
    a, b, sa, sb = _both(df, _weighted_pair, ("x", "y"), B, seed, strata)
    assert a == b
    assert sa == sb


@SETTINGS
@given(
    df=frames(),
    B=st.integers(min_value=1, max_value=30),
    seed=st.integers(min_value=0, max_value=2**32 - 1),
    strata=st.sampled_from([None, "st"]),
)
def test_columns_are_passed_in_the_order_given(
    df: pd.DataFrame, B: int, seed: int, strata: str | None
) -> None:
    """``columns=("y", "x")`` passes y first: equals the frame path with the arguments swapped."""
    a, b, sa, sb = _both(df, _weighted_pair, ("y", "x"), B, seed, strata)
    assert a == b
    assert sa == sb


def test_columns_path_example_strata_on_and_off() -> None:
    """A fixed example (several rows per TU, two strata) so the property has a concrete anchor."""
    df = pd.DataFrame(
        {
            "tu": ["a", "a", "b", "c", "c", "c", "d", "e", "e", "f"],
            "x": [0.5, -1.0, 2.0, 0.25, 3.0, -0.75, 1.5, -2.0, 0.0, 4.0],
            "y": [1.0, 0.0, 2.0, -1.0, 0.0, 3.0, 1.0, 2.0, -2.0, 0.0],
            "st": [0, 0, 1, 1, 1, 1, 0, 1, 1, 0],
        }
    )
    for strata in (None, "st"):
        a, b, sa, sb = _both(df, _weighted_pair, ("x", "y"), 200, 7, strata)
        assert a == b
        assert sa == sb
        assert a.n_valid == 200
