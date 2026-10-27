"""fgsea port I: the simple permutation stage, NES, the simple p, ``log2err`` and the multilevel
decision (issue M8-02; GIPPA_SPEC.md §7.7.4, §2.4 G15, §6.6).

Oracles from fgsea **1.36.2** (Bioconductor 3.22 tarball): ``R/fgsea.R`` ``fgseaSimpleImpl`` and
``R/fgseaMultilevel.R`` (simple-stage part), ``src/fastGSEA.cpp`` ``calcGseaStatCumulativeBatch``.

- ``leZero`` / ``geZero`` count random ES ``<= 0`` / ``>= 0``: a random ES of exactly 0 counts on
  **both** sides (pinned on N = 3, where a single hit at rank 2 is an ES = 0 tie);
- NES = ES / geZeroMean if ES > 0, ES / |leZeroMean| if ES <= 0;
- p = min((1 + nLeEs)/(1 + nLeZero), (1 + nGeEs)/(1 + nGeZero)). For ES != 0 this is the
  sign-conditional p (1 + nMoreExtreme)/(1 + modeFraction) exactly; for ES = 0 it is 1;
- nMoreExtreme = nGeEs if ES > 0 else nLeEs; modeFraction = nGeZero if ES >= 0 else nLeZero;
- modeFraction < 10: p, NES, log2err NaN and the set never goes to the multilevel stage;
- log2err = sqrt(trigamma(nMoreExtreme + 1) - trigamma(nPerm + 1)) / ln 2;
- needs_multilevel = multilevelError(p0, sampleSize) < 0.5 max(crude - left, right - crude).

R's random stream is not reproduced (spec item 6); these tests use exact count identities and
hand-computed limits with Monte-Carlo tolerances of >= 5 standard errors.
"""

from __future__ import annotations

import numpy as np
import pytest
from helpers.fgsea_oracle import simple_and_multilevel_errors, simple_log2err

from gippa.rng import for_step
from gippa.stats.fgsea import Prepared, SimpleResult, fgsea_simple, gsea_es, prepare

STREAM = "function.gsea"
SAMPLE_SIZE = 101


def labels(n: int) -> list[str]:
    return [f"g{i:05d}" for i in range(n)]


def run(prep: Prepared, n_perm: int, seed: int = 20261021) -> SimpleResult:
    return fgsea_simple(prep, n_perm, SAMPLE_SIZE, for_step(STREAM, seed))


def check_internal_consistency(prep: Prepared, res: SimpleResult, n_perm: int) -> None:
    """Every identity that holds whatever the random draws."""
    k = len(prep.names)
    for field in ("es", "nes", "p", "log2err", "n_more_extreme", "mode_fraction"):
        assert np.asarray(getattr(res, field)).shape == (k,), field
    assert np.asarray(res.needs_multilevel).shape == (k,)
    assert np.asarray(res.needs_multilevel).dtype == np.bool_
    assert len(res.leading_edge) == k
    order = np.asarray(prep.order)
    for i in range(k):
        es, le = gsea_es(prep.stats, prep.members[i])
        assert res.es[i] == es
        assert [int(x) for x in res.leading_edge[i]] == [int(order[r - 1]) for r in le]
        nme, mf = int(res.n_more_extreme[i]), int(res.mode_fraction[i])
        assert 0 <= mf <= n_perm
        if mf < 10:
            assert np.isnan(res.p[i]) and np.isnan(res.nes[i]) and np.isnan(res.log2err[i])
            assert not res.needs_multilevel[i]
            continue
        if np.isnan(res.nes[i]):  # NES mean 0: p NaN as well
            assert np.isnan(res.p[i])
            continue
        assert 1 / (n_perm + 1) <= res.p[i] <= 1
        if es != 0:
            assert 0 <= nme <= mf
            assert res.p[i] == pytest.approx((1 + nme) / (1 + mf), rel=0, abs=1e-15)
            assert np.sign(res.nes[i]) == np.sign(es)
        else:
            assert res.p[i] == 1.0
            assert res.nes[i] == 0.0
        np.testing.assert_allclose(res.log2err[i], simple_log2err(nme, n_perm), rtol=0, atol=1e-12)
        simple_err, mult_err = simple_and_multilevel_errors(nme, n_perm, SAMPLE_SIZE)
        if abs(mult_err - simple_err) > 1e-9:
            assert bool(res.needs_multilevel[i]) == (mult_err < simple_err), i


def example_6() -> Prepared:
    rng = np.random.default_rng(1)
    stats = rng.standard_normal(1000)
    lab = labels(1000)
    sets = {
        f"set{j:02d}": [lab[i] for i in rng.choice(1000, size=15, replace=False)] for j in range(20)
    }
    return prepare(stats, lab, sets, 15, 500)


def assert_same(a: SimpleResult, b: SimpleResult) -> None:
    for field in ("es", "nes", "p", "log2err", "n_more_extreme", "mode_fraction"):
        np.testing.assert_array_equal(getattr(a, field), getattr(b, field), err_msg=field)
    np.testing.assert_array_equal(a.needs_multilevel, b.needs_multilevel)
    assert len(a.leading_edge) == len(b.leading_edge)
    for x, y in zip(a.leading_edge, b.leading_edge, strict=True):
        np.testing.assert_array_equal(x, y)


# ---------------------------------------------------------------------------------------------
# Example 6: 1,000 normal stats, 20 random sets of size 15, n_perm 1000
# ---------------------------------------------------------------------------------------------


def test_example_6_consistency_floor_and_signs() -> None:
    prep = example_6()
    assert len(prep.names) == 20
    res = run(prep, 1000)
    assert isinstance(res, SimpleResult)
    check_internal_consistency(prep, res, 1000)
    defined = ~np.isnan(res.p)
    assert defined.sum() == 20  # balanced normal stats: every set estimated
    assert np.all(res.p[defined] >= 1 / 1001)
    assert np.all(np.sign(res.nes[defined]) == np.sign(res.es[defined]))
    # |NES| of a null random set is O(1)
    assert np.all(np.abs(res.nes[defined]) < 4)


def test_example_6_same_seed_identical() -> None:
    prep = example_6()
    assert_same(run(prep, 1000, seed=7), run(prep, 1000, seed=7))


def test_example_6_different_seed_differs() -> None:
    prep = example_6()
    a, b = run(prep, 1000, seed=7), run(prep, 1000, seed=8)
    np.testing.assert_array_equal(a.es, b.es)  # the observed ES never depends on the seed
    assert not np.array_equal(a.mode_fraction, b.mode_fraction) or not np.array_equal(
        a.n_more_extreme, b.n_more_extreme
    )


def test_example_6_identical_across_numba_thread_counts() -> None:
    numba = pytest.importorskip("numba")
    max_threads = int(numba.config.NUMBA_NUM_THREADS)
    if max_threads < 2:
        pytest.skip("needs NUMBA_NUM_THREADS >= 2")
    prep = example_6()
    before = numba.get_num_threads()
    try:
        numba.set_num_threads(1)
        one = run(prep, 1000, seed=11)
        numba.set_num_threads(min(8, max_threads))
        many = run(prep, 1000, seed=11)
    finally:
        numba.set_num_threads(before)
    assert_same(one, many)


def test_n_perm_not_a_multiple_of_the_chunk() -> None:
    prep = example_6()
    res = run(prep, 250)
    check_internal_consistency(prep, res, 250)
    assert np.all(np.asarray(res.mode_fraction) <= 250)
    assert_same(res, run(prep, 250))


def test_does_not_mutate_prepared() -> None:
    prep = example_6()
    stats0, order0 = np.array(prep.stats, copy=True), np.array(prep.order, copy=True)
    members0 = [np.array(m, copy=True) for m in prep.members]
    sizes0 = np.array(prep.sizes, copy=True)
    run(prep, 200)
    np.testing.assert_array_equal(prep.stats, stats0)
    np.testing.assert_array_equal(prep.order, order0)
    np.testing.assert_array_equal(prep.sizes, sizes0)
    for m, m0 in zip(prep.members, members0, strict=True):
        np.testing.assert_array_equal(m, m0)


# ---------------------------------------------------------------------------------------------
# Example 7: the top 50 of 2,000 normal stats -> nMoreExtreme 0, needs multilevel
# ---------------------------------------------------------------------------------------------


def test_example_7_top_set_needs_multilevel() -> None:
    rng = np.random.default_rng(2)
    stats = rng.standard_normal(2000)
    lab = labels(2000)
    top = [lab[i] for i in np.argsort(-stats, kind="stable")[:50]]
    prep = prepare(stats, lab, {"top": top}, 15, 500)
    res = run(prep, 1000)
    check_internal_consistency(prep, res, 1000)
    assert res.es[0] > 0.9
    assert int(res.n_more_extreme[0]) == 0
    mf = int(res.mode_fraction[0])
    assert res.p[0] == pytest.approx(1 / (1 + mf), rel=0, abs=1e-15)
    assert res.p[0] >= 1 / 1001
    assert bool(res.needs_multilevel[0]) is True
    np.testing.assert_allclose(res.log2err[0], simple_log2err(0, 1000), rtol=0, atol=1e-12)
    # leading edge: input indices, top hit first
    le = [int(x) for x in res.leading_edge[0]]
    assert le[0] == int(np.argmax(stats))


# ---------------------------------------------------------------------------------------------
# Example 8: random ES of almost one sign -> modeFraction < 10 -> p, NES, log2err NaN
# ---------------------------------------------------------------------------------------------


def test_example_8_mode_fraction_below_ten_is_nan() -> None:
    # stats 2^(-i/4): within any random 15-set the top member dominates the weight, so the random
    # ES is >= 0 unless the first hit falls below about rank N/2 (probability ~ 2^-15).
    n = 1000
    stats = 2.0 ** (-np.arange(n) / 4)
    lab = labels(n)
    sets = {"bottom": lab[-15:], "top": lab[:15]}
    prep = prepare(stats, lab, sets, 15, 500)
    res = run(prep, 1000)
    check_internal_consistency(prep, res, 1000)
    b, t = prep.names.index("bottom"), prep.names.index("top")
    assert res.es[b] == pytest.approx(-1.0, abs=1e-12)
    assert int(res.mode_fraction[b]) < 10
    assert np.isnan(res.p[b]) and np.isnan(res.nes[b]) and np.isnan(res.log2err[b])
    assert not res.needs_multilevel[b]
    # the ES and leading edge are still reported (fgsea keeps them)
    assert [int(x) for x in res.leading_edge[b]] == list(range(n - 1, n - 16, -1))
    # the top set is estimated
    assert int(res.mode_fraction[t]) >= 10
    assert np.isfinite(res.p[t]) and np.isfinite(res.nes[t]) and res.nes[t] > 0


# ---------------------------------------------------------------------------------------------
# Example 9 and argument validation
# ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("n_perm", [0, -1])
def test_example_9_n_perm_below_one_raises(n_perm: int) -> None:
    prep = prepare([5, 4, 3, 2, 1], ["a", "b", "c", "d", "e"], {"s": ["a", "b"]}, 1, 4)
    with pytest.raises(ValueError):
        fgsea_simple(prep, n_perm, SAMPLE_SIZE, for_step(STREAM, 1))


@pytest.mark.parametrize("sample_size", [2, 0])
def test_sample_size_below_three_raises(sample_size: int) -> None:
    prep = prepare([5, 4, 3, 2, 1], ["a", "b", "c", "d", "e"], {"s": ["a", "b"]}, 1, 4)
    with pytest.raises(ValueError):
        fgsea_simple(prep, 100, sample_size, for_step(STREAM, 1))


def test_no_surviving_set_gives_zero_length_result() -> None:
    prep = prepare([5, 4, 3, 2, 1], ["a", "b", "c", "d", "e"], {"none": ["x"]}, 1, 4)
    res = fgsea_simple(prep, 100, SAMPLE_SIZE, for_step(STREAM, 1))
    for field in ("es", "nes", "p", "log2err", "n_more_extreme", "mode_fraction"):
        assert np.asarray(getattr(res, field)).shape == (0,), field
    assert np.asarray(res.needs_multilevel).shape == (0,)
    assert res.leading_edge == ()


# ---------------------------------------------------------------------------------------------
# Hand-computed permutation limits on tiny lists (single-hit random ES are exact fractions)
# ---------------------------------------------------------------------------------------------
# A single hit at rank r of N: tops = 1 - (r-1)/(N-1), bottoms = -(r-1)/(N-1). With K = 1 the
# random rank is uniform on 1..N.


def test_n3_random_es_zero_counts_on_both_sides() -> None:
    # N 3: random ES is 1 (r 1), 0 (r 2, the tie), -1 (r 3), each with probability 1/3.
    # Set {a} (ES 1): nMoreExtreme = n1, modeFraction = nGeZero = n1 + n2,
    # geZeroMean = n1 / (n1 + n2) -> NES = modeFraction / nMoreExtreme -> 2; p -> 1/2.
    # Excluding the zeros from geZero would give NES 1 instead.
    n_perm = 20000
    prep = prepare([3, 2, 1], ["a", "b", "c"], {"top": ["a"], "bottom": ["c"]}, 1, 2)
    res = run(prep, n_perm)
    check_internal_consistency(prep, res, n_perm)
    top, bot = prep.names.index("top"), prep.names.index("bottom")
    assert res.es[top] == 1.0 and res.es[bot] == -1.0
    for i, sign in ((top, 1.0), (bot, -1.0)):
        nme, mf = int(res.n_more_extreme[i]), int(res.mode_fraction[i])
        np.testing.assert_allclose(res.nes[i], sign * mf / nme, rtol=1e-12, atol=0)
        assert res.nes[i] == pytest.approx(2.0 * sign, abs=0.1)  # sd ~ 0.017
        assert res.p[i] == pytest.approx(0.5, abs=0.03)  # sd ~ 0.004
        assert mf == pytest.approx(2 * n_perm / 3, abs=500)  # sd ~ 67
        assert nme == pytest.approx(n_perm / 3, abs=500)


def test_n3_observed_es_zero_tie() -> None:
    # Example 4's set {y} (ES 0 tie): NES 0 (the ES <= 0 branch), p = min(1, 1) = 1,
    # nMoreExtreme = nLeEs = nLeZero, modeFraction = nGeZero (ES >= 0), empty leading edge.
    n_perm = 20000
    prep = prepare([2, 2, 1], ["x", "y", "z"], {"s": ["y"]}, 1, 2)
    res = run(prep, n_perm)
    check_internal_consistency(prep, res, n_perm)
    assert res.es[0] == 0.0
    assert res.nes[0] == 0.0
    assert res.p[0] == 1.0
    assert len(res.leading_edge[0]) == 0
    assert int(res.n_more_extreme[0]) == pytest.approx(2 * n_perm / 3, abs=500)
    assert int(res.mode_fraction[0]) == pytest.approx(2 * n_perm / 3, abs=500)


def test_n4_nes_normalisation() -> None:
    # N 4: random single-hit ES is 1, 2/3, -2/3, -1 (r = 1..4), each 1/4.
    # Set {a}: geZeroMean = (n1 + 2/3 n2)/(n1 + n2) -> 5/6, NES -> 6/5; p -> n1/(n1 + n2) -> 1/2.
    # Exactly: nMoreExtreme = n1, modeFraction = n1 + n2.
    # A second set of size 2 makes K = 2: the size-1 random ES is the first draw's, still uniform.
    n_perm = 20000
    sets = {"a": ["a"], "d": ["d"], "ab": ["a", "b"]}
    prep = prepare([4, 3, 2, 1], ["a", "b", "c", "d"], sets, 1, 3)
    res = run(prep, n_perm)
    check_internal_consistency(prep, res, n_perm)
    for name, sign in (("a", 1.0), ("d", -1.0)):
        i = prep.names.index(name)
        n1, mf = int(res.n_more_extreme[i]), int(res.mode_fraction[i])
        n2 = mf - n1
        np.testing.assert_allclose(res.nes[i], sign * mf / (n1 + 2 / 3 * n2), rtol=1e-12, atol=0)
        assert res.nes[i] == pytest.approx(1.2 * sign, abs=0.03)  # sd ~ 0.003
        assert res.p[i] == pytest.approx(0.5, abs=0.03)
    # ES of {a, b}: S = [1, 2], NR = 7, tops [4/7, 1] -> 1
    assert res.es[prep.names.index("ab")] == pytest.approx(1.0, abs=1e-12)


def test_leading_edge_in_input_indices_and_fgsea_order() -> None:
    # Examples 2 and 3 through fgsea_simple: input indices, the negative edge reversed.
    prep = prepare(
        [5, 4, 3, 2, 1], ["a", "b", "c", "d", "e"], {"up": ["a", "b"], "down": ["d", "e"]}, 1, 4
    )
    res = run(prep, 200)
    assert [int(x) for x in res.leading_edge[prep.names.index("up")]] == [0, 1]
    assert [int(x) for x in res.leading_edge[prep.names.index("down")]] == [4, 3]
    prep = prepare([3, -1, 2, -2], ["a", "b", "c", "d"], {"s": ["a", "d"]}, 1, 3)
    res = run(prep, 200)
    assert [int(x) for x in res.leading_edge[0]] == [0]
    assert res.es[0] == pytest.approx(0.6, abs=1e-12)
