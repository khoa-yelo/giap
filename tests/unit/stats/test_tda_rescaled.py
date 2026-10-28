"""``stats.tda.picked_fdr_rescaled``: ADR-0009 picked FDR with decoy wins rescaled by 1/f
(code review #237 B3; M7-11 lead decisions 10 and 13; M4-12 ``_picked_q`` uses it).

FDR(t) = (D(t)/f + 1)/T(t) over the picked list (ties go to the decoy, both within a gene and
across the list), made monotone and clipped to 1; f = 1 gives ``picked_fdr``; f must lie in
(0, 1], so a paired share of 0 raises (decision 10).
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.stats.tda import picked_fdr, picked_fdr_rescaled

NAN = math.nan


def close(actual: object, expected: object, atol: float = 1e-12) -> None:
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float), expected, rtol=0, atol=atol, equal_nan=True
    )


# ---------------------------------------------------------------- f outside (0, 1]


@pytest.mark.parametrize("f", [0.0, -0.1, 1.1, NAN, -1.0, math.inf])
def test_paired_share_outside_unit_interval_raises(f: float) -> None:
    with pytest.raises(ValueError, match="paired_share"):
        picked_fdr_rescaled([3.0, 2.0], [1.0, 2.5], f)


@pytest.mark.parametrize("f", [0.0, -0.1, 1.1, NAN])
def test_paired_share_is_checked_even_on_empty_input(f: float) -> None:
    with pytest.raises(ValueError, match="paired_share"):
        picked_fdr_rescaled([], [], f)


def test_paired_share_one_is_accepted() -> None:
    picked_fdr_rescaled([3.0], [1.0], 1.0)


# ---------------------------------------------------------------- f = 1 is picked_fdr

T_EX = [20.0, 19.0, 18.0, 17.0, 16.0, 9.0, 8.0, 6.0, NAN, 4.0]
D_EX = [NAN, NAN, NAN, NAN, 5.0, 9.0, NAN, NAN, 6.0, NAN]


def _same(a, b) -> None:  # type: ignore[no-untyped-def]
    np.testing.assert_array_equal(a.decoy_wins, b.decoy_wins)
    close(a.q, b.q)
    assert (a.n_target_wins, a.n_decoy_wins) == (b.n_target_wins, b.n_decoy_wins)


def test_f_one_equals_picked_fdr_on_the_hand_example() -> None:
    _same(picked_fdr_rescaled(T_EX, D_EX, 1.0), picked_fdr(T_EX, D_EX))


_score = st.one_of(st.integers(0, 6).map(float), st.just(NAN), st.just(-math.inf))


@settings(max_examples=300, deadline=None)
@given(st.lists(st.tuples(_score, _score), max_size=30))
def test_f_one_equals_picked_fdr_property(pairs: list[tuple[float, float]]) -> None:
    t = [a for a, _ in pairs]
    d = [b for _, b in pairs]
    _same(picked_fdr_rescaled(t, d, 1.0), picked_fdr(t, d))


@settings(max_examples=300, deadline=None)
@given(
    st.lists(st.tuples(_score, _score), max_size=30),
    st.floats(0.05, 1.0),
    st.floats(0.05, 1.0),
)
def test_rescaled_q_is_monotone_in_f_and_keeps_the_picks(
    pairs: list[tuple[float, float]], f1: float, f2: float
) -> None:
    """Smaller f inflates decoy wins: q never decreases; picks and counts never change; q is a
    monotone function of the target score and lies in (0, 1]."""
    lo, hi = sorted((f1, f2))
    t = [a for a, _ in pairs]
    d = [b for _, b in pairs]
    a, b = picked_fdr_rescaled(t, d, lo), picked_fdr_rescaled(t, d, hi)
    np.testing.assert_array_equal(a.decoy_wins, b.decoy_wins)
    assert (a.n_target_wins, a.n_decoy_wins) == (b.n_target_wins, b.n_decoy_wins)
    np.testing.assert_array_equal(np.isnan(a.q), np.isnan(b.q))
    ok = ~np.isnan(a.q)
    assert np.all(a.q[ok] >= b.q[ok] - 1e-12)
    assert np.all((a.q[ok] > 0) & (a.q[ok] <= 1))
    ts = np.asarray(t, dtype=float)[ok]
    qs = a.q[ok]
    o = np.argsort(-ts, kind="stable")
    assert np.all(np.diff(qs[o]) >= -1e-12)  # better score, never a larger q


# ---------------------------------------------------------------- hand-computed, f = 0.8


def test_hand_computed_f_08_with_a_tie_between_a_target_and_a_decoy_win() -> None:
    """Picked list: target wins 20, 19, 18, 17, 16 (its decoy 5 loses), 8, 6, 4; decoy wins 9
    (a within-gene tie 9 = 9 goes to the decoy) and 6 (a decoy-only gene, tied with the target
    win at 6, so it counts at t = 6). With f = 0.8 each decoy win counts 1.25:

    t:   20   19   18    17    16    8               6                 4
    T:    1    2    3     4     5    6               7                 8
    D:    0    0    0     0     0    1               2                 2
    FDR:  1   1/2  1/3   1/4   1/5  (1.25+1)/6=.375 (2.5+1)/7=.5     3.5/8=.4375

    Monotone from the loose end: q = .2 for 20..16, .375 at 8, .4375 at 6 and 4.
    """
    res = picked_fdr_rescaled(T_EX, D_EX, 0.8)
    close(res.q, [0.2, 0.2, 0.2, 0.2, 0.2, NAN, 0.375, 0.4375, NAN, 0.4375])
    np.testing.assert_array_equal(
        res.decoy_wins, [False, False, False, False, False, True, False, False, True, False]
    )
    assert (res.n_target_wins, res.n_decoy_wins) == (8, 2)


def test_hand_computed_f_08_before_monotone_clips_at_one() -> None:
    """Target 3 over its decoy 1, decoy-only win 5, target 2: FDR(3) = (1.25 + 1)/1 = 2.25 and
    FDR(2) = (1.25 + 1)/2 = 1.125, both clipped to 1."""
    res = picked_fdr_rescaled([3.0, NAN, 2.0], [1.0, 5.0, NAN], 0.8)
    close(res.q, [1.0, NAN, 1.0])
    assert (res.n_target_wins, res.n_decoy_wins) == (2, 1)


def test_hand_computed_f_05_no_decoy_win_above() -> None:
    """Targets 4 and 3 above a decoy win at 2, target 1: f = 0.5 gives FDR 1, 1/2, then
    (2 + 1)/3 = 1 at t = 1 (clipped), so q = [.5, .5, 1]."""
    res = picked_fdr_rescaled([4.0, 3.0, NAN, 1.0], [NAN, NAN, 2.0, NAN], 0.5)
    close(res.q, [0.5, 0.5, NAN, 1.0])


# ---------------------------------------------------------------- degenerate


def test_empty_input() -> None:
    res = picked_fdr_rescaled([], [], 0.8)
    assert res.q.shape == (0,) and res.decoy_wins.shape == (0,)
    assert (res.n_target_wins, res.n_decoy_wins) == (0, 0)


def test_no_target_win() -> None:
    """Every gene is won by the decoy (a higher decoy, a tie, a decoy-only gene) or is empty."""
    res = picked_fdr_rescaled([1.0, 2.0, NAN, NAN], [3.0, 2.0, 0.5, NAN], 0.8)
    close(res.q, [NAN, NAN, NAN, NAN])
    np.testing.assert_array_equal(res.decoy_wins, [True, True, True, False])
    assert (res.n_target_wins, res.n_decoy_wins) == (0, 3)


def test_length_mismatch_raises_as_picked_fdr() -> None:
    with pytest.raises(ValueError):
        picked_fdr_rescaled([1.0, 2.0], [1.0], 0.8)
