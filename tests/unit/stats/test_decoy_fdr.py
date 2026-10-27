"""Distance-reweighted decoy FDR, monotone q and the local fdr (issue M6-01; §7.4.4, §8.2, §8.3).

Spec examples 4-10 of docs/specs/M6-01-stats-tdist-decoy-fdr.md, hand-computed, plus every row of
its "Degenerate cases" table that concerns the decoy FDR. Scores are ``s = -|Delta|`` (higher is
better, as everywhere in ``gippa.stats``), bins come from ``bin_index`` and

    FDR(t) = pi0 * sum_b N_NO,b * (#{decoy_b : s >= t} + 1) / (N_decoy,b + 1) / #{target : s >= t}

with ``N_NO,b`` the count of **all** targets in bin b (spec question 1, proposed default). Readings
this file pins, all three the spec's own "proposed default" (questions 1-3):

- ``q`` is the running minimum of the element-wise maximum FDR **from the loose end**, i.e.
  ``q_i = min over t <= s_i of max_d FDR_d(t)``, clipped to [0, 1]: the only direction that makes q
  an upper bound on the FDP of the list it selects (``gippa.stats.tda.monotone_q`` already does this
  for the separate-search q, and the two must agree);
- the maximum is taken **first** and monotonised after (step 5), so q is >= the element-wise maximum
  of the per-decoy monotone curves; both orders are computed here so the order is pinned;
- ``fdr`` is the raw step function, unclipped (example 7 reports 10/3 while ``q`` is 1.0).

The statistical claims (FDP control, the lfdr under a shared null) are in
tests/calibration/test_operon_fdr_calibration.py.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

import numpy as np
import numpy.typing as npt
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy import stats as sps

from gippa.stats.decoy_fdr import (
    DecoyFdr,
    _lfdr_densities,
    _reflected_density,
    bin_index,
    local_fdr,
    max_decoy_fdr,
    reweighted_fdr,
)
from gippa.stats.tda import monotone_q

FAST = settings(max_examples=40, deadline=None)
NAN = float("nan")
INF = float("inf")
#: ``operons.distance_bins_bp`` (spec step 3): five half-open bins
EDGES: tuple[float, ...] = (-INF, 0.0, 50.0, 150.0, 300.0, INF)


def points(
    spec: Mapping[int, Sequence[tuple[float, int]]],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]:
    """(scores, bins) from ``{bin: [(score, count), ...]}``."""
    s: list[float] = []
    b: list[int] = []
    for bin_, items in spec.items():
        for score, count in items:
            s.extend([score] * count)
            b.extend([bin_] * count)
    return np.array(s, dtype=np.float64), np.array(b, dtype=np.int64)


# Examples 4-6: 200 targets in bin 0 and 100 in bin 1; three score levels 1.0 > 0.0 > -1.0.
# At t = 1.0: 40 targets, 5 + 2 convergent decoys, 1 + 1 sandwich decoys.
# At t = 0.0: 100 targets, 20 + 15 convergent decoys, 1 + 1 sandwich decoys.
TARGET_S, TARGET_B = points({0: [(1.0, 30), (0.0, 40), (-1.0, 130)],
                             1: [(1.0, 10), (0.0, 20), (-1.0, 70)]})  # fmt: skip
CV_S, CV_B = points({0: [(1.0, 5), (0.0, 15), (-1.0, 80)],
                     1: [(1.0, 2), (0.0, 13), (-1.0, 85)]})  # fmt: skip
SW_S, SW_B = points({0: [(1.0, 1), (-1.0, 49)],
                     1: [(1.0, 1), (-1.0, 19)]})  # fmt: skip
EX4 = 0.3712871287128713  # (200*6/101 + 100*3/101)/40
EX5 = 0.5742574257425743  # (200*21/101 + 100*16/101)/100
EX6 = 0.4341736694677871  # (200*2/51 + 100*2/21)/40
DECOYS: Mapping[str, tuple[npt.NDArray[np.float64], npt.NDArray[np.int64]]] = {
    "convergent": (CV_S, CV_B),
    "sandwich": (SW_S, SW_B),
}


def at(values: npt.NDArray[np.float64], score: float) -> npt.NDArray[np.float64]:
    """The ``values`` of every target whose score is ``score`` (they must all be equal)."""
    sel = values[np.equal(TARGET_S, score)]
    assert sel.size > 0
    return sel


# ------------------------------------------------------------------------ bin_index


def test_example8_bin_index() -> None:
    got = bin_index([-5.0, 0.0, 49.0, 50.0, 1000.0, NAN], EDGES)
    assert got.tolist() == [0, 1, 1, 2, 4, -1]
    assert got.dtype == np.int64


def test_bins_are_half_open_on_the_left() -> None:
    got = bin_index([0.0, 49.999, 50.0, 149.0, 150.0, 299.0, 300.0, 1e9], EDGES)
    assert got.tolist() == [1, 1, 2, 2, 3, 3, 4, 4]


def test_bin_index_on_an_empty_input() -> None:
    got = bin_index([], EDGES)
    assert got.shape == (0,)
    assert got.dtype == np.int64


@FAST
@given(d=st.lists(st.floats(-1000.0, 100_000.0) | st.just(NAN), max_size=30))
def test_property_bin_index_is_monotone_and_in_range(d: list[float]) -> None:
    got = bin_index(d, EDGES)
    assert got.shape == (len(d),)
    assert np.all((got >= -1) & (got <= len(EDGES) - 2))
    finite = [x for x in d if not math.isnan(x)]
    idx = bin_index(sorted(finite), EDGES)
    assert np.all(np.diff(idx) >= 0)
    for x, i in zip(finite, bin_index(finite, EDGES), strict=True):
        assert EDGES[i] <= x < EDGES[i + 1]


# -------------------------------------------------------------------- reweighted_fdr


def test_example4_and_5_convergent_decoys() -> None:
    fdr = reweighted_fdr(TARGET_S, TARGET_B, CV_S, CV_B, n_bins=2)
    assert fdr.shape == TARGET_S.shape
    assert at(fdr, 1.0) == pytest.approx(EX4, rel=1e-12)
    assert at(fdr, 0.0) == pytest.approx(EX5, rel=1e-12)
    # every decoy is above the loosest threshold, so the +1 terms cancel exactly: FDR = pi0
    assert at(fdr, -1.0) == pytest.approx(1.0, rel=1e-12)


def test_example6_sandwich_decoys_give_a_larger_fdr() -> None:
    fdr = reweighted_fdr(TARGET_S, TARGET_B, SW_S, SW_B, n_bins=2)
    assert at(fdr, 1.0) == pytest.approx(EX6, rel=1e-12)
    assert EX6 > EX4


def test_example7_a_bin_with_targets_and_no_decoys_contributes_its_target_count() -> None:
    """10 targets in bin 0, 3 of them above t, no decoy in that bin: FDR = 10/3, unclipped."""
    t_s, t_b = points({0: [(1.0, 3), (-1.0, 7)]})
    d_s, d_b = points({1: [(1.0, 1)]})  # the decoy set is not empty, but bin 0 has no decoy
    fdr = reweighted_fdr(t_s, t_b, d_s, d_b, n_bins=2)
    assert float(fdr[t_s == 1.0][0]) == pytest.approx(10.0 / 3.0, rel=1e-12)
    assert float(fdr[t_s == -1.0][0]) == pytest.approx(1.0, rel=1e-12)


def test_an_empty_decoy_set_leaves_the_maximally_conservative_value() -> None:
    t_s, t_b = points({0: [(1.0, 3), (-1.0, 7)]})
    fdr = reweighted_fdr(t_s, t_b, [], [], n_bins=1)
    assert float(fdr[t_s == 1.0][0]) == pytest.approx(10.0 / 3.0, rel=1e-12)


def test_pi0_multiplies_the_fdr_exactly() -> None:
    one = reweighted_fdr(TARGET_S, TARGET_B, CV_S, CV_B, n_bins=2)
    half = reweighted_fdr(TARGET_S, TARGET_B, CV_S, CV_B, n_bins=2, pi0=0.5)
    assert half == pytest.approx(0.5 * one, rel=1e-12)


def test_ties_in_the_target_score_share_one_fdr() -> None:
    fdr = reweighted_fdr(TARGET_S, TARGET_B, CV_S, CV_B, n_bins=2)
    for score in (1.0, 0.0, -1.0):
        assert len(set(at(fdr, score).tolist())) == 1


def test_targets_without_a_bin_or_a_score_are_excluded() -> None:
    t_s = np.concatenate([TARGET_S, [2.0, NAN]])
    t_b = np.concatenate([TARGET_B, [-1, 0]])
    fdr = reweighted_fdr(t_s, t_b, CV_S, CV_B, n_bins=2)
    assert math.isnan(float(fdr[-1])), "a NaN score is excluded"
    assert math.isnan(float(fdr[-2])), "bin -1 (NaN distance) is excluded everywhere"
    # the excluded pair changes nothing for the rest
    assert fdr[:-2] == pytest.approx(
        reweighted_fdr(TARGET_S, TARGET_B, CV_S, CV_B, n_bins=2), rel=1e-12, nan_ok=True
    )


def test_decoys_without_a_bin_are_excluded() -> None:
    d_s = np.concatenate([CV_S, [5.0, 5.0, 5.0]])
    d_b = np.concatenate([CV_B, [-1, -1, -1]])
    assert reweighted_fdr(TARGET_S, TARGET_B, d_s, d_b, n_bins=2) == pytest.approx(
        reweighted_fdr(TARGET_S, TARGET_B, CV_S, CV_B, n_bins=2), rel=1e-12
    )


def test_no_targets_gives_an_empty_result() -> None:
    fdr = reweighted_fdr([], [], CV_S, CV_B, n_bins=2)
    assert fdr.shape == (0,)
    assert fdr.dtype == np.float64


@pytest.mark.parametrize("pi0", [0.0, -0.5, 1.5, NAN])
def test_pi0_outside_the_unit_interval_is_rejected(pi0: float) -> None:
    with pytest.raises(ValueError):
        reweighted_fdr(TARGET_S, TARGET_B, CV_S, CV_B, n_bins=2, pi0=pi0)
    with pytest.raises(ValueError):
        max_decoy_fdr(TARGET_S, TARGET_B, DECOYS, n_bins=2, pi0=pi0)
    with pytest.raises(ValueError):
        local_fdr(TARGET_S, CV_S, np.ones_like(CV_S), pi0=pi0)


def test_mismatched_lengths_are_rejected() -> None:
    with pytest.raises(ValueError):
        reweighted_fdr(TARGET_S, TARGET_B[:-1], CV_S, CV_B, n_bins=2)
    with pytest.raises(ValueError):
        reweighted_fdr(TARGET_S, TARGET_B, CV_S, CV_B[:-1], n_bins=2)


# --------------------------------------------------------------------- max_decoy_fdr


def test_example6_the_larger_decoy_fdr_wins_and_is_recorded() -> None:
    res = max_decoy_fdr(TARGET_S, TARGET_B, DECOYS, n_bins=2)
    assert isinstance(res, DecoyFdr)
    assert at(res.fdr, 1.0) == pytest.approx(EX6, rel=1e-12)
    assert set(at(res.which, 1.0).tolist()) == {"sandwich"}
    assert at(res.fdr, 0.0) == pytest.approx(EX5, rel=1e-12)
    assert set(at(res.which, 0.0).tolist()) == {"convergent"}


def test_q_is_the_monotone_maximum_clipped_to_one() -> None:
    res = max_decoy_fdr(TARGET_S, TARGET_B, DECOYS, n_bins=2)
    assert at(res.q, 1.0) == pytest.approx(EX6, rel=1e-12)
    assert at(res.q, 0.0) == pytest.approx(EX5, rel=1e-12)
    assert at(res.q, -1.0) == pytest.approx(1.0, rel=1e-12)
    assert np.all((res.q >= 0.0) & (res.q <= 1.0))


def test_example7_q_saturates_at_one() -> None:
    t_s, t_b = points({0: [(1.0, 3), (-1.0, 7)]})
    d_s, d_b = points({1: [(1.0, 1)]})
    res = max_decoy_fdr(t_s, t_b, {"convergent": (d_s, d_b)}, n_bins=2)
    assert float(res.fdr[t_s == 1.0][0]) == pytest.approx(10.0 / 3.0, rel=1e-12)
    assert float(res.q[t_s == 1.0][0]) == 1.0


def test_the_lfdr_field_has_the_target_shape_and_is_a_probability_or_nan() -> None:
    """The spec gives ``DecoyFdr.lfdr`` but does not say which decoy set feeds it, so
    only its shape and range are pinned here."""
    res = max_decoy_fdr(TARGET_S, TARGET_B, DECOYS, n_bins=2)
    assert res.lfdr.shape == TARGET_S.shape
    assert np.all(np.isnan(res.lfdr) | ((res.lfdr >= 0.0) & (res.lfdr <= 1.0)))


def test_ties_share_fdr_q_and_which_and_the_tie_break_is_deterministic() -> None:
    first = max_decoy_fdr(TARGET_S, TARGET_B, DECOYS, n_bins=2)
    second = max_decoy_fdr(TARGET_S, TARGET_B, DECOYS, n_bins=2)
    for score in (1.0, 0.0, -1.0):
        assert len(set(at(first.fdr, score).tolist())) == 1
        assert len(set(at(first.q, score).tolist())) == 1
        assert len(set(at(first.which, score).tolist())) == 1
    assert first.which.tolist() == second.which.tolist()
    assert set(first.which.tolist()) <= set(DECOYS)


def test_a_decoy_set_without_any_decoy_is_skipped() -> None:
    empty: tuple[list[float], list[int]] = ([], [])
    res = max_decoy_fdr(TARGET_S, TARGET_B, {"convergent": (CV_S, CV_B), "sandwich": empty}, 2)
    only_cv = max_decoy_fdr(TARGET_S, TARGET_B, {"convergent": (CV_S, CV_B)}, n_bins=2)
    assert res.fdr == pytest.approx(only_cv.fdr, rel=1e-12)
    assert set(res.which.tolist()) == {"convergent"}


def test_no_decoy_set_survives_gives_all_nan() -> None:
    empty: tuple[list[float], list[int]] = ([], [])
    res = max_decoy_fdr(TARGET_S, TARGET_B, {"convergent": empty, "sandwich": empty}, n_bins=2)
    assert np.isnan(res.q).all()
    assert np.isnan(res.fdr).all()
    assert res.which.shape == TARGET_S.shape


def test_no_targets_gives_empty_arrays() -> None:
    res = max_decoy_fdr([], [], DECOYS, n_bins=2)
    assert res.q.shape == (0,)
    assert res.fdr.shape == (0,)
    assert res.lfdr.shape == (0,)
    assert res.which.shape == (0,)


def test_one_decoy_set_reproduces_reweighted_fdr() -> None:
    res = max_decoy_fdr(TARGET_S, TARGET_B, {"convergent": (CV_S, CV_B)}, n_bins=2)
    assert res.fdr == pytest.approx(
        reweighted_fdr(TARGET_S, TARGET_B, CV_S, CV_B, n_bins=2), rel=1e-12
    )


# --------------------------------------------------- properties of the decoy FDR


def random_scenario(
    rng: np.random.Generator, n_bins: int, n_target: int, n_decoy: int
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.int64], npt.NDArray[np.float64],
           npt.NDArray[np.int64]]:  # fmt: skip
    t_s = -np.abs(rng.normal(size=n_target))
    t_b = rng.integers(0, n_bins, size=n_target)
    d_s = -np.abs(rng.normal(size=n_decoy))
    d_b = rng.integers(0, n_bins, size=n_decoy)
    return t_s, t_b.astype(np.int64), d_s, d_b.astype(np.int64)


@FAST
@given(seed=st.integers(0, 2**32 - 1), n_bins=st.integers(1, 5), n_t=st.integers(1, 60))
def test_property_q_is_non_increasing_in_the_score(seed: int, n_bins: int, n_t: int) -> None:
    rng = np.random.default_rng(seed)
    t_s, t_b, d_s, d_b = random_scenario(rng, n_bins, n_t, 40)
    res = max_decoy_fdr(t_s, t_b, {"convergent": (d_s, d_b)}, n_bins=n_bins)
    order = np.argsort(-t_s, kind="stable")
    q = res.q[order]
    assert np.all(np.diff(q) >= -1e-12), "q must not fall as the score falls"
    assert np.all((q >= 0.0) & (q <= 1.0))


@FAST
@given(seed=st.integers(0, 2**32 - 1), n_bins=st.integers(1, 5))
def test_property_maximise_then_monotonise(seed: int, n_bins: int) -> None:
    """Acceptance: q is the monotone curve of the element-wise maximum, which dominates the
    element-wise maximum of the per-decoy monotone curves (the order of operations is pinned)."""
    rng = np.random.default_rng(seed)
    t_s, t_b, cv_s, cv_b = random_scenario(rng, n_bins, 50, 40)
    sw_s, sw_b = random_scenario(rng, n_bins, 1, 25)[2:]
    decoys = {"convergent": (cv_s, cv_b), "sandwich": (sw_s, sw_b)}
    res = max_decoy_fdr(t_s, t_b, decoys, n_bins=n_bins)
    raw = {k: reweighted_fdr(t_s, t_b, s, b, n_bins=n_bins) for k, (s, b) in decoys.items()}
    stacked = np.vstack(list(raw.values()))
    assert res.fdr == pytest.approx(stacked.max(axis=0), rel=1e-12)

    order = np.argsort(-t_s, kind="stable")
    expected = np.empty_like(res.q)
    expected[order] = monotone_q(res.fdr[order])
    assert res.q == pytest.approx(expected, rel=1e-12)

    per_decoy = []
    for value in raw.values():
        one = np.empty_like(value)
        one[order] = monotone_q(value[order])
        per_decoy.append(one)
    assert np.all(res.q >= np.vstack(per_decoy).max(axis=0) - 1e-12)


@FAST
@given(seed=st.integers(0, 2**32 - 1), n_bins=st.integers(1, 4), pi0=st.floats(0.05, 1.0))
def test_property_fdr_at_the_loosest_threshold_is_pi0(seed: int, n_bins: int, pi0: float) -> None:
    """Every decoy is above the lowest target score, so the +1 terms cancel and FDR = pi0."""
    rng = np.random.default_rng(seed)
    t_s, t_b, d_s, d_b = random_scenario(rng, n_bins, 40, 30)
    d_s = d_s + 10.0  # every decoy is above every target, so every bin ratio is exactly 1
    fdr = reweighted_fdr(t_s, t_b, d_s, d_b, n_bins=n_bins, pi0=pi0)
    assert float(fdr[np.argmin(t_s)]) == pytest.approx(pi0, rel=1e-12)


@FAST
@given(seed=st.integers(0, 2**32 - 1), n_bins=st.integers(1, 4))
def test_property_scaling_every_count_converges_to_the_plus_one_free_fdr(
    seed: int, n_bins: int
) -> None:
    """Acceptance: multiplying every count by k leaves the FDR unchanged up to the +1 terms."""
    rng = np.random.default_rng(seed)
    t_s, t_b, d_s, d_b = random_scenario(rng, n_bins, 25, 20)
    mid = int(np.argsort(-t_s, kind="stable")[t_s.size // 2])  # a threshold with real counts
    values = []
    for k in (1, 10, 200):
        fdr = reweighted_fdr(
            np.tile(t_s, k), np.tile(t_b, k), np.tile(d_s, k), np.tile(d_b, k), n_bins=n_bins
        )
        values.append(float(fdr[: t_s.size][mid]))
    limit = _plus_one_free(t_s, t_b, d_s, d_b, n_bins, float(t_s[mid]))
    assert values[0] >= values[1] >= values[2] >= limit - 1e-9
    assert values[2] == pytest.approx(limit, abs=0.02 * max(1.0, limit))


def _plus_one_free(
    t_s: npt.NDArray[np.float64],
    t_b: npt.NDArray[np.int64],
    d_s: npt.NDArray[np.float64],
    d_b: npt.NDArray[np.int64],
    n_bins: int,
    t: float,
) -> float:
    """The same sum without the +1 corrections (0/0 bins contribute nothing)."""
    total = 0.0
    for b in range(n_bins):
        n_no = float((t_b == b).sum())
        big = float(((d_b == b) & (d_s >= t)).sum())
        all_d = float((d_b == b).sum())
        total += n_no * (big / all_d if all_d else 1.0)
    return total / float((t_s >= t).sum())


@FAST
@given(seed=st.integers(0, 2**32 - 1), n_bins=st.integers(1, 4))
def test_property_results_follow_the_input_order(seed: int, n_bins: int) -> None:
    rng = np.random.default_rng(seed)
    t_s, t_b, d_s, d_b = random_scenario(rng, n_bins, 30, 25)
    perm = rng.permutation(t_s.size)
    plain = max_decoy_fdr(t_s, t_b, {"convergent": (d_s, d_b)}, n_bins=n_bins)
    shuffled = max_decoy_fdr(t_s[perm], t_b[perm], {"convergent": (d_s, d_b)}, n_bins=n_bins)
    assert shuffled.q == pytest.approx(plain.q[perm], rel=1e-12, nan_ok=True)
    assert shuffled.fdr == pytest.approx(plain.fdr[perm], rel=1e-12, nan_ok=True)
    assert shuffled.which.tolist() == plain.which[perm].tolist()


# ------------------------------------------------------------------------ local_fdr


def test_example9_same_distribution_gives_lfdr_near_one() -> None:
    rng = np.random.default_rng(609)
    target = -np.abs(rng.normal(size=2000))
    decoy = -np.abs(rng.normal(size=2000))
    lfdr = local_fdr(target, decoy, np.ones(2000))
    assert lfdr.shape == target.shape
    assert np.all(lfdr <= 1.0)
    assert 0.8 <= float(np.nanmedian(lfdr)) <= 1.2


def test_example10_lfdr_is_small_at_small_abs_delta_and_rises() -> None:
    """Example 10, read as a trend: a KDE ratio is a noisy estimate, so "rising monotonically over
    the central 90 %" is checked as a rank correlation plus strictly increasing decile means."""
    rng = np.random.default_rng(610)
    target = -np.abs(rng.normal(scale=0.05, size=2000))
    decoy = -np.abs(rng.normal(scale=1.0, size=2000))
    lfdr = local_fdr(target, decoy, np.ones(2000))
    order = np.argsort(-target)  # by increasing |Delta|
    assert float(lfdr[order[0]]) < 0.2, "the smallest |Delta| is clearly enriched"
    lo, hi = int(0.05 * 2000), int(0.95 * 2000)
    central_x, central_y = (-target)[order][lo:hi], lfdr[order][lo:hi]
    assert float(sps.spearmanr(central_x, central_y).statistic) > 0.95
    deciles = [float(np.mean(chunk)) for chunk in np.array_split(central_y, 10)]
    assert deciles == sorted(deciles)
    assert deciles[-1] > 3.0 * deciles[0]


@pytest.mark.parametrize("n", [0, 1, 4])
def test_fewer_than_five_values_leaves_the_lfdr_undefined(n: int) -> None:
    rng = np.random.default_rng(611)
    many = -np.abs(rng.normal(size=50))
    few = -np.abs(rng.normal(size=n))
    assert np.isnan(local_fdr(many, few, np.ones(n))).all()
    out = local_fdr(few, many, np.ones(50))
    assert out.shape == (n,)
    assert np.isnan(out).all()


def test_local_fdr_is_invariant_to_the_weight_normalisation() -> None:
    rng = np.random.default_rng(612)
    target = -np.abs(rng.normal(size=300))
    decoy = -np.abs(rng.normal(size=400))
    w = rng.uniform(0.1, 3.0, size=400)
    assert local_fdr(target, decoy, w) == pytest.approx(
        local_fdr(target, decoy, 7.5 * w), rel=1e-10, nan_ok=True
    )


def test_local_fdr_weights_move_the_null_density() -> None:
    """A weighting that upweights the decoys nearest 0 raises the lfdr of the best targets."""
    rng = np.random.default_rng(613)
    target = -np.abs(rng.normal(scale=0.3, size=500))
    decoy = -np.abs(rng.normal(scale=1.0, size=800))
    near = np.where(-decoy < 0.3, 10.0, 1.0)
    best = np.argsort(-target)[:50]
    flat = local_fdr(target, decoy, np.ones(800))
    tilted = local_fdr(target, decoy, near)
    assert float(np.nanmean(tilted[best])) > float(np.nanmean(flat[best]))


def test_local_fdr_scales_with_pi0() -> None:
    rng = np.random.default_rng(614)
    target = -np.abs(rng.normal(size=400))
    decoy = -np.abs(rng.normal(size=400))
    one = local_fdr(target, decoy, np.ones(400))
    half = local_fdr(target, decoy, np.ones(400), pi0=0.5)
    free = one < 1.0 - 1e-9
    assert free.sum() > 0
    assert half[free] == pytest.approx(0.5 * one[free], rel=1e-10)


@FAST
@given(seed=st.integers(0, 2**32 - 1), n=st.integers(5, 80))
def test_property_local_fdr_is_a_probability(seed: int, n: int) -> None:
    rng = np.random.default_rng(seed)
    target = -np.abs(rng.normal(size=n))
    decoy = -np.abs(rng.normal(size=n + 3))
    lfdr = local_fdr(target, decoy, rng.uniform(0.5, 2.0, size=n + 3))
    assert lfdr.shape == (n,)
    assert np.all(np.isnan(lfdr) | ((lfdr >= 0.0) & (lfdr <= 1.0)))


@FAST
@given(seed=st.integers(0, 2**32 - 1))
def test_property_local_fdr_is_scale_equivariant(seed: int) -> None:
    """|Delta| has no natural unit: rescaling target and decoy together leaves the ratio alone
    (both densities pick up the same Jacobian and Silverman's bandwidth rescales with them)."""
    rng = np.random.default_rng(seed)
    target = -np.abs(rng.normal(size=300))
    decoy = -np.abs(rng.normal(size=300))
    plain = local_fdr(target, decoy, np.ones(300))
    scaled = local_fdr(3.7 * target, 3.7 * decoy, np.ones(300))
    assert scaled == pytest.approx(plain, rel=1e-8, nan_ok=True)


def test_reflected_density_integrates_to_one_on_the_half_line() -> None:
    """Lead decision 8: the reflection at 0 keeps all the mass on [0, inf) (a plain KDE would
    leak about half a kernel below 0)."""
    sample = np.abs(np.random.default_rng(3).normal(size=500))
    w = np.ones_like(sample)
    x = np.linspace(0.0, 8.0, 4001)
    f = _reflected_density(sample, w, x)
    assert np.trapezoid(f, x) == pytest.approx(1.0, abs=1e-3)


def test_lfdr_densities_are_finite_and_non_negative_at_every_target() -> None:
    rng = np.random.default_rng(4)
    t, d = rng.normal(size=300), rng.normal(size=300)
    f0, f_no = _lfdr_densities(t, d, np.ones_like(d))
    assert f0.shape == f_no.shape == t.shape
    assert np.all(np.isfinite(f0)) and np.all(f0 >= 0.0)
    assert np.all(np.isfinite(f_no)) and np.all(f_no > 0.0)
