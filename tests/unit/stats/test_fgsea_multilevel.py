"""fgsea port II: the multilevel splitting stage (issue M8-11; GIAP_SPEC.md §7.7.4, §8.2, §8.3,
§6.6; docs/specs/M8-11-stats-fgsea-multilevel.md examples 1-6 and every degenerate row).

Interface pinned here (spec "Public interface"):

- ``giap.stats.fgsea.fgsea_multilevel(prep, simple, sample_size, eps, rng) -> SimpleResult``: the
  simple result with ``p`` and ``log2err`` upgraded for exactly the sets in
  ``simple.needs_multilevel``. Every other field (``es``, ``nes``, ``leading_edge``,
  ``n_more_extreme``, ``mode_fraction``, ``needs_multilevel``) is returned unchanged, and so are
  ``p`` / ``log2err`` of the sets that did not need the stage;
- ``giap.stats.fgsea.MultilevelResult`` is a frozen dataclass with fields ``p``, ``log2err``,
  ``is_cp_ge_half`` (the per-set kernel output; its construction is not pinned);
- ``sample_size < 3`` raises ``ValueError``; an even ``sample_size`` behaves as ``sample_size + 1``
  (same draws, same output);
- p of a multilevel set: ``0 < p <= 1``; with ``eps > 0`` a p below eps is reported as ``eps`` with
  ``log2err`` NaN; otherwise ``log2err`` is ``multilevelError(p, sampleSize)`` or NaN (the
  ``isCpGeHalf`` false case, which keeps p);
- one random stream per size group, spawned in ascending size order: adding a group of a *larger*
  size never changes the result of a smaller group (the i-th child of ``spawn`` does not depend on
  how many children are spawned);
- one ruler per (size, sign): within a size group and sign, ``p * denomProb`` (fgsea's
  ``cppMPval``, ``denomProb = (modeFraction + 1)/(nPerm + 1)``) is non-increasing in |ES|.

Test data come from fixed ``numpy.random.default_rng(seed)`` draws (as spec example 1 does); every
generator handed to GIAP comes from ``giap.rng.for_step("function.gsea", seed)``.
Golden R values are M8-13's job; nothing here is compared with R.
"""

from __future__ import annotations

import dataclasses
import itertools
import math

import numpy as np
import pytest
from scipy.special import polygamma

from giap.rng import for_step
from giap.stats.fgsea import Prepared, SimpleResult, fgsea_simple, prepare

STREAM = "function.gsea"
SAMPLE_SIZE = 101
N_PERM = 1000


def _ml():  # type: ignore[no-untyped-def]
    from giap.stats.fgsea import fgsea_multilevel

    return fgsea_multilevel


def multilevel_error(p: float, sample_size: int = SAMPLE_SIZE) -> float:
    """Spec item 3: sqrt(floor(-log2 p + 1) (trigamma((s + 1)/2) - trigamma(s + 1))) / ln 2."""
    single = float(polygamma(1, (sample_size + 1) / 2) - polygamma(1, sample_size + 1))
    return math.sqrt(math.floor(-math.log2(p) + 1) * single) / math.log(2)


def labels(n: int) -> list[str]:
    return [f"g{i:05d}" for i in range(n)]


def normal_stats(n: int, seed: int = 1) -> np.ndarray:
    return np.random.default_rng(seed).standard_normal(n)


def top(x: np.ndarray, k: int) -> list[int]:
    return [int(i) for i in np.argsort(-x, kind="stable")[:k]]


def bottom(x: np.ndarray, k: int) -> list[int]:
    return [int(i) for i in np.argsort(x, kind="stable")[:k]]


def prep_of(x: np.ndarray, sets: dict[str, list[int]], lo: int = 1, hi: int = 10**6) -> Prepared:
    lab = labels(len(x))
    return prepare(x, lab, {s: [lab[i] for i in idx] for s, idx in sets.items()}, lo, hi)


def simple_of(prep: Prepared, seed: int = 1, n_perm: int = N_PERM) -> SimpleResult:
    return fgsea_simple(prep, n_perm, SAMPLE_SIZE, for_step(STREAM, seed))


def ml(
    prep: Prepared,
    simple: SimpleResult,
    seed: int = 2,
    eps: float = 0.0,
    sample_size: int = SAMPLE_SIZE,
) -> SimpleResult:
    return _ml()(prep, simple, sample_size, eps, for_step(STREAM, seed))


def assert_unchanged_fields(simple: SimpleResult, out: SimpleResult) -> None:
    np.testing.assert_array_equal(out.es, simple.es)
    np.testing.assert_array_equal(out.nes, simple.nes)
    np.testing.assert_array_equal(out.n_more_extreme, simple.n_more_extreme)
    np.testing.assert_array_equal(out.mode_fraction, simple.mode_fraction)
    np.testing.assert_array_equal(out.needs_multilevel, simple.needs_multilevel)
    assert len(out.leading_edge) == len(simple.leading_edge)
    for a, b in zip(out.leading_edge, simple.leading_edge, strict=True):
        np.testing.assert_array_equal(a, b)
    keep = ~np.asarray(simple.needs_multilevel)
    np.testing.assert_array_equal(out.p[keep], simple.p[keep])
    np.testing.assert_array_equal(out.log2err[keep], simple.log2err[keep])


def assert_multilevel_values(simple: SimpleResult, out: SimpleResult, eps: float = 0.0) -> None:
    for i in np.flatnonzero(simple.needs_multilevel):
        p, err = float(out.p[i]), float(out.log2err[i])
        assert math.isfinite(p) and 0 < p <= 1, (i, p)
        if eps > 0 and p == eps:
            assert math.isnan(err), (i, err)
            continue
        assert p >= eps, (i, p)
        if not math.isnan(err):  # NaN only in the isCpGeHalf-false case
            assert err == pytest.approx(multilevel_error(p), rel=1e-9), (i, p, err)


def assert_same(a: SimpleResult, b: SimpleResult) -> None:
    np.testing.assert_array_equal(a.p, b.p)
    np.testing.assert_array_equal(a.log2err, b.log2err)


# --------------------------------------------------------------------------- interface


def test_multilevel_result_is_a_frozen_dataclass_with_the_spec_fields() -> None:
    from giap.stats.fgsea import MultilevelResult

    assert dataclasses.is_dataclass(MultilevelResult)
    names = [f.name for f in dataclasses.fields(MultilevelResult)]
    assert names == ["p", "log2err", "is_cp_ge_half"]
    r = MultilevelResult(np.array([0.1]), np.array([0.5]), np.array([True]))
    with pytest.raises(dataclasses.FrozenInstanceError):
        r.p = np.array([0.2])  # type: ignore[misc]


def test_fgsea_multilevel_is_exported() -> None:
    import giap.stats.fgsea as fg

    assert "fgsea_multilevel" in fg.__all__
    assert "MultilevelResult" in fg.__all__


@pytest.mark.parametrize("sample_size", [0, 1, 2, -5])
def test_sample_size_below_three_raises(sample_size: int) -> None:
    x = normal_stats(300)
    prep = prep_of(x, {"top": top(x, 15)})
    simple = simple_of(prep)
    with pytest.raises(ValueError):
        ml(prep, simple, sample_size=sample_size)


def test_multilevel_error_reference_values() -> None:
    """The spec's quoted values of multilevelError(p, 101): 0.64 at 1e-6, 1.18 at 1e-20."""
    assert multilevel_error(1e-6) == pytest.approx(0.64, abs=0.005)
    assert multilevel_error(1e-20) == pytest.approx(1.18, abs=0.005)


# --------------------------------------------------------------------------- examples


@pytest.fixture(scope="module")
def example_1() -> tuple[Prepared, SimpleResult, SimpleResult]:
    """The top 50 of 2,000 N(0, 1) stats (seed 1) as a set; nPerm 1000."""
    x = normal_stats(2000, seed=1)
    prep = prep_of(x, {"top50": top(x, 50)})
    simple = simple_of(prep, seed=1)
    return prep, simple, ml(prep, simple, seed=2)


def test_example_1_simple_p_at_the_floor_and_needs_multilevel(
    example_1: tuple[Prepared, SimpleResult, SimpleResult],
) -> None:
    _, simple, _ = example_1
    assert simple.n_more_extreme[0] == 0
    assert simple.p[0] == pytest.approx(1 / (simple.mode_fraction[0] + 1), rel=1e-12)
    assert bool(simple.needs_multilevel[0])


def test_example_1_multilevel_p_below_1e_minus_10(
    example_1: tuple[Prepared, SimpleResult, SimpleResult],
) -> None:
    _, simple, out = example_1
    p = float(out.p[0])
    assert 0 < p < 1e-10, p
    assert out.log2err[0] == pytest.approx(multilevel_error(p), rel=1e-9)
    assert out.log2err[0] > multilevel_error(1e-10)
    assert_unchanged_fields(simple, out)


def test_example_2_two_seeds_agree_within_the_m8_13_criterion(
    example_1: tuple[Prepared, SimpleResult, SimpleResult],
) -> None:
    prep, simple, a = example_1
    b = ml(prep, simple, seed=3)
    pa, pb = float(a.p[0]), float(b.p[0])
    assert pa != pb  # two streams, two estimates
    err = max(float(a.log2err[0]), float(b.log2err[0]))
    assert abs(math.log2(pa) - math.log2(pb)) < 3 * math.sqrt(2) * err, (pa, pb, err)


def test_example_3_a_set_with_simple_p_0_3_is_not_multilevel_and_unchanged() -> None:
    x = normal_stats(2000, seed=1)
    pick = np.random.default_rng(3)  # test data: random 50-sets until one has simple p near 0.3
    for _ in range(200):
        idx = [int(i) for i in pick.choice(2000, size=50, replace=False)]
        prep = prep_of(x, {"random50": idx})
        simple = simple_of(prep, seed=1)
        if 0.2 <= simple.p[0] <= 0.4:
            break
    else:  # pragma: no cover - data guard
        pytest.fail("no random set with simple p near 0.3")
    assert not bool(simple.needs_multilevel[0])
    out = ml(prep, simple, seed=2)
    assert_unchanged_fields(simple, out)
    assert out.p[0] == simple.p[0]
    assert out.log2err[0] == simple.log2err[0]


def test_example_4_bottom_50_negative_es_reversed_ruler() -> None:
    x = normal_stats(2000, seed=1)
    prep = prep_of(x, {"bottom50": bottom(x, 50)})
    simple = simple_of(prep, seed=1)
    assert simple.es[0] < 0
    assert bool(simple.needs_multilevel[0])
    out = ml(prep, simple, seed=2)
    assert 0 < out.p[0] < 1e-10
    assert out.log2err[0] == pytest.approx(multilevel_error(float(out.p[0])), rel=1e-9)
    assert_unchanged_fields(simple, out)


def test_example_4_mirrored_stats_give_a_symmetric_p(
    example_1: tuple[Prepared, SimpleResult, SimpleResult],
) -> None:
    """Negating every stat turns the top-50 set into the bottom-50 set with the same |stats|
    profile: the exact p is identical, so the two estimates agree within the tolerance."""
    _, simple_pos, pos = example_1
    x = -normal_stats(2000, seed=1)
    prep = prep_of(x, {"mirror": bottom(x, 50)})
    simple = simple_of(prep, seed=1)
    assert simple.es[0] == pytest.approx(-simple_pos.es[0], abs=1e-12)
    neg = ml(prep, simple, seed=2)
    pp, pn = float(pos.p[0]), float(neg.p[0])
    err = max(float(pos.log2err[0]), float(neg.log2err[0]))
    assert abs(math.log2(pp) - math.log2(pn)) < 3 * math.sqrt(2) * err, (pp, pn, err)


def test_example_5_identical_across_numba_thread_counts() -> None:
    numba = pytest.importorskip("numba")
    max_threads = int(numba.config.NUMBA_NUM_THREADS)
    if max_threads < 2:
        pytest.skip("needs NUMBA_NUM_THREADS >= 2")
    x = normal_stats(1500, seed=4)
    o = np.argsort(-x, kind="stable")
    sets = {
        "a40": [int(i) for i in o[:40]],
        "b25": [int(i) for i in o[:20]] + [int(i) for i in o[700:705]],
        "c25": bottom(x, 25),
        "d60": [int(i) for i in o[5:65]],
        "e30": [int(i) for i in o[::50]],
    }
    prep = prep_of(x, sets)
    simple = simple_of(prep, seed=5)
    assert int(np.sum(simple.needs_multilevel)) >= 3
    before = numba.get_num_threads()
    try:
        numba.set_num_threads(1)
        one = ml(prep, simple, seed=6)
        numba.set_num_threads(min(8, max_threads))
        many = ml(prep, simple, seed=6)
    finally:
        numba.set_num_threads(before)
    assert_same(one, many)


def test_example_6_eps_floors_p_and_drops_log2err(
    example_1: tuple[Prepared, SimpleResult, SimpleResult],
) -> None:
    prep, simple, _ = example_1
    out = ml(prep, simple, seed=2, eps=1e-10)
    assert out.p[0] == 1e-10
    assert math.isnan(out.log2err[0])
    assert_unchanged_fields(simple, out)


# --------------------------------------------------------------------------- degenerate rows


def test_no_set_needs_multilevel_returns_the_simple_result() -> None:
    x = normal_stats(800, seed=7)
    o = np.argsort(-x, kind="stable")
    prep = prep_of(x, {f"s{j}": [int(i) for i in o[j::40]] for j in range(5)})
    simple = simple_of(prep, seed=1)
    assert not np.any(simple.needs_multilevel)
    out = ml(prep, simple, seed=2)
    assert_unchanged_fields(simple, out)
    np.testing.assert_array_equal(out.p, simple.p)
    np.testing.assert_array_equal(out.log2err, simple.log2err)


def test_no_sets_at_all() -> None:
    x = normal_stats(100, seed=7)
    prep = prep_of(x, {"tiny": [0]}, lo=10, hi=50)  # dropped by the size filter
    assert prep.names == ()
    simple = simple_of(prep, seed=1)
    out = ml(prep, simple, seed=2)
    assert out.p.shape == (0,) and out.log2err.shape == (0,)


def test_unestimated_sets_are_never_refined() -> None:
    """modeFraction < 10 sets (p NaN) are excluded even when the mask says otherwise."""
    x = normal_stats(600, seed=8)
    prep = prep_of(x, {"top": top(x, 20), "mid": [int(i) for i in np.argsort(-x)[::30]]})
    simple = simple_of(prep, seed=1)
    forced = dataclasses.replace(
        simple,
        p=np.array([np.nan, simple.p[1]]),
        nes=np.array([np.nan, simple.nes[1]]),
        log2err=np.array([np.nan, simple.log2err[1]]),
        mode_fraction=np.array([3, simple.mode_fraction[1]], dtype=np.int64),
        needs_multilevel=np.array([True, bool(simple.needs_multilevel[1])]),
    )
    out = ml(prep, forced, seed=2)
    assert math.isnan(out.p[0]) and math.isnan(out.log2err[0]) and math.isnan(out.nes[0])


def test_extreme_es_top_30_of_3000_reaches_a_positive_p() -> None:
    x = normal_stats(3000, seed=9)
    prep = prep_of(x, {"top30": top(x, 30)})
    simple = simple_of(prep, seed=1)
    assert bool(simple.needs_multilevel[0])
    out = ml(prep, simple, seed=2)
    p, err = float(out.p[0]), float(out.log2err[0])
    assert math.isfinite(p) and 0 < p < 1e-10, p  # never 0 by construction of the levels
    assert err == pytest.approx(multilevel_error(p), rel=1e-9)
    assert err > multilevel_error(1e-10)  # log2err grows with depth


def test_one_size_group_with_both_signs_runs_two_rulers() -> None:
    x = normal_stats(2000, seed=1)
    prep = prep_of(x, {"up": top(x, 50), "down": bottom(x, 50)})
    assert list(prep.sizes) == [50, 50]
    simple = simple_of(prep, seed=1)
    assert np.all(simple.needs_multilevel)
    out = ml(prep, simple, seed=2)
    i_up, i_down = prep.names.index("up"), prep.names.index("down")
    assert out.es[i_up] > 0 > out.es[i_down]
    for i in (i_up, i_down):
        assert 0 < out.p[i] < 1e-10
    assert_unchanged_fields(simple, out)
    assert_multilevel_values(simple, out)


def test_even_sample_size_is_rounded_up() -> None:
    x = normal_stats(1000, seed=10)
    prep = prep_of(x, {"top": top(x, 25), "down": bottom(x, 30)})
    simple = simple_of(prep, seed=1)
    assert np.all(simple.needs_multilevel)
    even = ml(prep, simple, seed=2, sample_size=100)
    odd = ml(prep, simple, seed=2, sample_size=101)
    assert_same(even, odd)


# --------------------------------------------------------------------------- determinism, streams


@pytest.fixture(scope="module")
def mixed() -> tuple[Prepared, SimpleResult]:
    """Two size groups (25 and 40), both signs, plus sets that do not need the stage."""
    x = normal_stats(1200, seed=11)
    o = np.argsort(-x, kind="stable")
    sets = {
        "u25": [int(i) for i in o[:25]],
        "d25": bottom(x, 25),
        "u40": [int(i) for i in o[3:43]],
        "n25": [int(i) for i in o[::48]][:25],
        "n40": [int(i) for i in o[1::30]][:40],
    }
    prep = prep_of(x, sets)
    simple = simple_of(prep, seed=1)
    return prep, simple


def test_same_seed_same_result(mixed: tuple[Prepared, SimpleResult]) -> None:
    prep, simple = mixed
    assert_same(ml(prep, simple, seed=12), ml(prep, simple, seed=12))


def test_different_seed_different_estimates(mixed: tuple[Prepared, SimpleResult]) -> None:
    prep, simple = mixed
    a, b = ml(prep, simple, seed=12), ml(prep, simple, seed=13)
    need = np.asarray(simple.needs_multilevel)
    assert need.any()
    assert not np.array_equal(a.p[need], b.p[need])


def test_global_numpy_state_is_not_used(mixed: tuple[Prepared, SimpleResult]) -> None:
    prep, simple = mixed
    np.random.seed(1)  # noqa: NPY002 - the legacy global state must not matter
    a = ml(prep, simple, seed=12)
    np.random.seed(999)  # noqa: NPY002
    b = ml(prep, simple, seed=12)
    assert_same(a, b)


def test_inputs_are_not_mutated(mixed: tuple[Prepared, SimpleResult]) -> None:
    prep, simple = mixed
    p, err = simple.p.copy(), simple.log2err.copy()
    ml(prep, simple, seed=12)
    np.testing.assert_array_equal(simple.p, p)
    np.testing.assert_array_equal(simple.log2err, err)


def test_mixed_invariants(mixed: tuple[Prepared, SimpleResult]) -> None:
    prep, simple = mixed
    out = ml(prep, simple, seed=12)
    assert_unchanged_fields(simple, out)
    assert_multilevel_values(simple, out)


def _subset(prep: Prepared, simple: SimpleResult, keep: list[int]) -> tuple[Prepared, SimpleResult]:
    sp = dataclasses.replace(
        prep,
        names=tuple(prep.names[i] for i in keep),
        members=tuple(prep.members[i] for i in keep),
        sizes=prep.sizes[keep],
    )
    ss = SimpleResult(
        *(np.asarray(getattr(simple, f))[keep] for f in (
            "es", "nes", "p", "log2err", "n_more_extreme", "mode_fraction", "needs_multilevel"
        )),
        tuple(simple.leading_edge[i] for i in keep),
    )  # fmt: skip
    return sp, ss


def test_a_larger_size_group_does_not_change_a_smaller_one() -> None:
    """Groups get spawned child streams in ascending size order (spec item 5)."""
    x = normal_stats(1500, seed=14)
    o = np.argsort(-x, kind="stable")
    prep = prep_of(x, {"small": [int(i) for i in o[:30]], "large": [int(i) for i in o[2:62]]})
    simple = simple_of(prep, seed=1)
    assert np.all(simple.needs_multilevel)
    i_small = prep.names.index("small")
    both = ml(prep, simple, seed=15)
    alone_prep, alone_simple = _subset(prep, simple, [i_small])
    alone = ml(alone_prep, alone_simple, seed=15)
    assert both.p[i_small] == alone.p[0]
    assert both.log2err[i_small] == alone.log2err[0] or (
        math.isnan(both.log2err[i_small]) and math.isnan(alone.log2err[0])
    )


def test_one_ruler_per_size_and_sign_orders_the_cpp_p_by_es() -> None:
    """Sets of one size and sign share a ruler, so p * denomProb is monotone in |ES|."""
    x = normal_stats(2000, seed=16)
    o = np.argsort(-x, kind="stable")
    sets = {
        "a": [int(i) for i in o[:40]],
        "b": [int(i) for i in o[:30]] + [int(i) for i in o[100:110]],
        "c": [int(i) for i in o[:20]] + [int(i) for i in o[200:220]],
        "d": [int(i) for i in o[10:50]],
    }
    prep = prep_of(x, sets)
    assert set(prep.sizes.tolist()) == {40}
    simple = simple_of(prep, seed=1)
    out = ml(prep, simple, seed=17)
    idx = [i for i in range(len(prep.names)) if simple.needs_multilevel[i] and out.p[i] < 1]
    assert len(idx) >= 3
    denom = (simple.mode_fraction + 1) / (N_PERM + 1)
    cpp = out.p * denom
    by_es = sorted(idx, key=lambda i: float(simple.es[i]))
    for lo_i, hi_i in itertools.pairwise(by_es):
        assert cpp[hi_i] <= cpp[lo_i] * (1 + 1e-12), (prep.names[lo_i], prep.names[hi_i])


@pytest.mark.parametrize("eps", [1e-3, 1e-6])
def test_eps_is_a_floor_for_every_refined_set(
    mixed: tuple[Prepared, SimpleResult], eps: float
) -> None:
    prep, simple = mixed
    out = ml(prep, simple, seed=12, eps=eps)
    assert_unchanged_fields(simple, out)
    assert_multilevel_values(simple, out, eps=eps)
    need = np.asarray(simple.needs_multilevel)
    assert np.all(out.p[need] >= eps)
    at = need & (out.p == eps)
    assert at.any()  # u25 / d25 are far below 1e-3
    assert np.all(np.isnan(out.log2err[at]))
