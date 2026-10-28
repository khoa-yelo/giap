"""fgsea port I: stat preparation, the observed ES and the leading edge (issue M8-02; GIPPA_SPEC.md
§7.7.4, §2.4 G15, §8.3 row ``fgsea_multilevel``).

Oracles are hand-computed from fgsea **1.36.2** (Bioconductor 3.22 source tarball
``fgsea_1.36.2.tar.gz``, sha256 d6996ff5...a525f5a): ``R/fgsea.R`` ``calcGseaStat`` and
``preparePathwaysAndStats``, ``R/util.R`` ``preparePathways``. Pinned behaviours (M8-02 lead
decision 4, checked against 1.36.2, not master):

- stats sorted decreasing and **stable** (R's radix order keeps ties in input order), then ``abs``;
- the ES = 0 tie (``maxP == -minP``) gives ES 0 and an **empty** leading edge;
- the negative-ES leading edge is ``rev(S[which.min(bottoms) .. m])`` (bottom hit first);
- the positive-ES leading edge stops at the **first** maximum (``which.max``);
- the NR = 0 branch uses equal hit weights 1/m;
- ``minSize`` is raised to 1 and ``maxSize`` lowered to N - 1 (``preparePathways``).

Integer scaling: fgsea 1.36.2 has **no** ``2^30 / sum`` scaling and no ``round`` in
``preparePathwaysAndStats`` (that step appears only after 1.36.2; GitHub master 1.39.4 lines 64-72).
``test_no_integer_scaling_in_1_36_2`` pins the 1.36.2 arithmetic on a case where rounding would move
the ES by 2.4e-10, well outside the 1e-12 tolerance.
"""

from __future__ import annotations

import copy

import numpy as np
import pytest
from helpers.fgsea_oracle import is_exact_tie, reference_gsea_stat

from gippa.stats.fgsea import Prepared, gsea_es, prepare

ATOL = 1e-12
ABC = ["a", "b", "c", "d", "e"]


def close(actual: object, expected: object, atol: float = ATOL) -> None:
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float), expected, rtol=0, atol=atol, equal_nan=True
    )


def ints(a: object) -> list[int]:
    arr = np.asarray(a)
    assert arr.ndim == 1
    assert arr.size == 0 or np.issubdtype(arr.dtype, np.integer), arr.dtype
    return [int(x) for x in arr]


def es_of(stats: list[float], labels: list[str], members: list[str]) -> tuple[float, list[int]]:
    """ES and leading edge (as input indices) of one set via ``prepare`` + ``gsea_es``."""
    prep = prepare(stats, labels, {"s": members}, 1, len(stats) - 1)
    assert prep.names == ("s",)
    es, le = gsea_es(prep.stats, prep.members[0])
    order = np.asarray(prep.order)
    return float(es), [int(order[r - 1]) for r in ints(le)]


# ---------------------------------------------------------------------------------------------
# Spec examples 1-5 (ES to 1e-12, leading edges and their order)
# ---------------------------------------------------------------------------------------------


def test_example_1_top_set_es_one() -> None:
    # tops [5/9, 1], bottoms [0, 5/9]
    es, le = es_of([5, 4, 3, 2, 1], ABC, ["a", "b"])
    close(es, 1.0)
    assert le == [0, 1]


def test_example_2_bottom_set_reversed_leading_edge() -> None:
    # tops [-1/3, 0], bottoms [-1, -1/3]; leading edge rev(S[1..2]) = [e, d]
    es, le = es_of([5, 4, 3, 2, 1], ABC, ["d", "e"])
    close(es, -1.0)
    assert le == [4, 3]


def test_example_2_gsea_es_returns_reversed_ranks() -> None:
    es, le = gsea_es(np.array([5.0, 4.0, 3.0, 2.0, 1.0]), np.array([4, 5], dtype=np.int64))
    close(es, -1.0)
    assert ints(le) == [5, 4]


def test_example_3_mixed_signs_sorted_then_abs() -> None:
    # sorted a(3), c(2), b(-1), d(-2) -> |.| = [3, 2, 1, 2]; S = [1, 4], NR = 5
    # tops [0.6, 0], bottoms [0, -0.4]
    prep = prepare([3, -1, 2, -2], ["a", "b", "c", "d"], {"s": ["a", "d"]}, 1, 3)
    assert ints(prep.order) == [0, 2, 1, 3]
    close(prep.stats, [3.0, 2.0, 1.0, 2.0], atol=0)
    assert ints(prep.members[0]) == [1, 4]
    es, le = gsea_es(prep.stats, prep.members[0])
    close(es, 0.6)
    assert ints(le) == [1]
    assert es_of([3, -1, 2, -2], ["a", "b", "c", "d"], ["a", "d"]) == (pytest.approx(0.6), [0])


def test_example_4_stable_tie_order_gives_es_zero_tie() -> None:
    # x, y, z kept in input order: S = [2], tops [0.5], bottoms [-0.5] -> maxP == -minP -> ES 0
    prep = prepare([2, 2, 1], ["x", "y", "z"], {"s": ["y"]}, 1, 2)
    assert ints(prep.order) == [0, 1, 2]
    assert ints(prep.members[0]) == [2]
    es, le = gsea_es(prep.stats, prep.members[0])
    assert es == 0.0
    assert ints(le) == []


def test_example_4_swapped_labels_puts_y_first() -> None:
    # labels [y, x, z]: y is now input 0 and keeps rank 1 -> tops [1], bottoms [0] -> ES 1
    es, le = es_of([2, 2, 1], ["y", "x", "z"], ["y"])
    close(es, 1.0)
    assert le == [0]


def test_example_5_equal_weights_first_maximum() -> None:
    # S = [1, 3]; tops [0.5, 0.5], bottoms [0, 0]: ES 0.5, which.max -> first -> [a]
    es, le = es_of([1, 1, 1, 1], ["a", "b", "c", "d"], ["a", "c"])
    close(es, 0.5)
    assert le == [0]


def test_negative_es_leading_edge_starts_at_first_minimum() -> None:
    # stats [1, 1, 1, 1], S = [2, 4]: rCumSum [1/2, 1], tops [1/2 - 1/2, 1 - 2/2] = [0, 0]
    # bottoms [-1/2, -1/2]: ES -1/2, which.min(bottoms) = 1 (first) -> rev(S[1..2]) = [4, 2]
    es, le = gsea_es(np.ones(4), np.array([2, 4], dtype=np.int64))
    close(es, -0.5)
    assert ints(le) == [4, 2]


def test_negative_es_partial_leading_edge() -> None:
    # stats [4, 3, 2, 1, 1, 1] (N 6), S = [1, 5, 6], NR = 6, N - m = 3
    # rCumSum [4/6, 5/6, 1]; tops [4/6, 5/6 - 1, 1 - 1] = [2/3, -1/6, 0]
    # bottoms [0, -1/3, -1/6]; maxP 2/3 > 1/3 -> ES 2/3, leading edge [1]
    es, le = gsea_es(np.array([4.0, 3, 2, 1, 1, 1]), np.array([1, 5, 6], dtype=np.int64))
    close(es, 2 / 3)
    assert ints(le) == [1]
    # S = [2, 3, 6], NR = 6: rCumSum [1/2, 5/6, 1]; tops [1/2 - 1/3, 5/6 - 1/3, 1 - 1]
    # = [1/6, 1/2, 0]; bottoms [-1/3, 1/6, -1/6]; maxP 1/2 > 1/3 -> ES 1/2, edge [2, 3]
    es, le = gsea_es(np.array([4.0, 3, 2, 1, 1, 1]), np.array([2, 3, 6], dtype=np.int64))
    close(es, 0.5)
    assert ints(le) == [2, 3]
    # S = [4, 5, 6] (weights 1 each): rCumSum [1/3, 2/3, 1]
    # tops [1/3 - 1, 2/3 - 1, 1 - 1] = [-2/3, -1/3, 0]; bottoms [-1, -2/3, -1/3]
    # ES -1, which.min = 1 -> rev(S) = [6, 5, 4]
    es, le = gsea_es(np.array([4.0, 3, 2, 1, 1, 1]), np.array([4, 5, 6], dtype=np.int64))
    close(es, -1.0)
    assert ints(le) == [6, 5, 4]


def test_negative_es_edge_excludes_hits_above_the_minimum() -> None:
    # stats [1]*6, S = [1, 5, 6]: rCumSum [1/3, 2/3, 1]; N - m = 3
    # tops [1/3, 2/3 - 1, 1 - 1] = [1/3, -1/3, 0]; bottoms [0, -2/3, -1/3]
    # maxP 1/3 < 2/3 -> ES -2/3; which.min(bottoms) = 2 -> rev(S[2..3]) = [6, 5]
    es, le = gsea_es(np.ones(6), np.array([1, 5, 6], dtype=np.int64))
    close(es, -2 / 3)
    assert ints(le) == [6, 5]


# ---------------------------------------------------------------------------------------------
# Degenerate cases
# ---------------------------------------------------------------------------------------------


def test_nr_zero_branch_equal_weights() -> None:
    # stats [1, 0, 0, 0], S = [2, 4]: NR = 0 -> rCumSum [1/2, 1]
    # tops [1/2 - 1/2, 1 - 2/2] = [0, 0]; bottoms = tops - 1/2 = [-1/2, -1/2]
    # ES -1/2; which.min first -> rev(S[1..2]) = [4, 2]
    es, le = gsea_es(np.array([1.0, 0.0, 0.0, 0.0]), np.array([2, 4], dtype=np.int64))
    close(es, -0.5)
    assert ints(le) == [4, 2]


def test_nr_zero_branch_through_prepare() -> None:
    # a(1), then b, c, d tied at 0 in input order; set {b, d} -> ranks [2, 4]
    es, le = es_of([1, 0, 0, 0], ["a", "b", "c", "d"], ["d", "b"])
    close(es, -0.5)
    assert le == [3, 1]


def test_all_selected_raises() -> None:
    with pytest.raises(ValueError):
        gsea_es(np.array([3.0, 2.0, 1.0]), np.array([1, 2, 3], dtype=np.int64))


def test_es_zero_tie_symmetric_single_hit() -> None:
    # N 5, S = [3]: tops [1 - 2/4] = [1/2], bottoms [-1/2] -> tie -> 0, empty edge
    es, le = gsea_es(np.array([5.0, 4, 3, 2, 1]), np.array([3], dtype=np.int64))
    assert es == 0.0
    assert ints(le) == []


# ---------------------------------------------------------------------------------------------
# Integer scaling: absent in fgsea 1.36.2
# ---------------------------------------------------------------------------------------------


def test_no_integer_scaling_in_1_36_2() -> None:
    """1.36.2 computes ES on the float |stats|.

    stats [0.9, 0.7, 0.31, 0.1], set {a, c}: S = [1, 3], NR = 1.21, tops [0.9/1.21, 1 - 1/2];
    ES = 90/121 = 0.743801652892562. Master's ``round(s * floor(2^30 / sum))`` would give
    0.7438016526496332 (2.4e-10 lower).
    """
    prep = prepare([0.9, 0.7, 0.31, 0.1], ["a", "b", "c", "d"], {"s": ["a", "c"]}, 1, 3)
    close(prep.stats, [0.9, 0.7, 0.31, 0.1], atol=0)
    es, le = gsea_es(prep.stats, prep.members[0])
    close(es, 90 / 121, atol=1e-13)
    assert abs(es - 0.7438016526496332) > 1e-10
    assert ints(le) == [1]


def test_prepared_stats_are_abs_of_sorted_floats() -> None:
    raw = [0.3, -2.5, 1.0e-7, 4.25, -0.3]
    prep = prepare(raw, ["a", "b", "c", "d", "e"], {"s": ["a", "b"]}, 1, 4)
    # decreasing, stable: d 4.25, a 0.3, c 1e-7, e -0.3, b -2.5
    assert ints(prep.order) == [3, 0, 2, 4, 1]
    close(prep.stats, [4.25, 0.3, 1.0e-7, 0.3, 2.5], atol=0)
    assert np.all(np.isfinite(np.asarray(prep.stats, dtype=float)))


# ---------------------------------------------------------------------------------------------
# prepare: mapping, filters, validation
# ---------------------------------------------------------------------------------------------


def test_prepare_drops_absent_members_and_duplicates() -> None:
    prep = prepare([3, -1, 2, -2], ["a", "b", "c", "d"], {"s": ["d", "a", "a", "zz"]}, 1, 3)
    assert ints(prep.members[0]) == [1, 4]
    assert ints(prep.sizes) == [2]
    assert prep.dropped == ()


def test_prepare_members_are_sorted_ranks() -> None:
    prep = prepare([5, 4, 3, 2, 1], ABC, {"s": ["e", "c", "a"]}, 1, 4)
    assert ints(prep.members[0]) == [1, 3, 5]


def test_prepare_size_filter_inclusive_and_after_dropping_absent() -> None:
    sets = {
        "one": ["a"],
        "two": ["a", "b"],
        "two_after_drop": ["c", "d", "x", "y", "z"],
        "three": ["a", "b", "c"],
        "four": ["a", "b", "c", "d"],
    }
    prep = prepare([5, 4, 3, 2, 1], ABC, sets, 2, 3)
    assert prep.names == ("two", "two_after_drop", "three")
    assert ints(prep.sizes) == [2, 2, 3]
    assert set(prep.dropped) == {"one", "four"}
    assert len(prep.members) == 3


def test_prepare_empty_overlap_is_dropped_even_with_min_size_zero() -> None:
    # preparePathways: minSize <- max(minSize, 1)
    prep = prepare([5, 4, 3, 2, 1], ABC, {"none": ["x", "y"], "s": ["a"]}, 0, 4)
    assert prep.names == ("s",)
    assert prep.dropped == ("none",)


def test_prepare_full_set_is_dropped_even_with_max_size_n() -> None:
    # preparePathways: maxSize <- min(maxSize, length(universe) - 1)
    prep = prepare([5, 4, 3, 2, 1], ABC, {"all": ABC, "s": ["a"]}, 1, 5)
    assert prep.names == ("s",)
    assert prep.dropped == ("all",)


def test_prepare_no_surviving_set() -> None:
    prep = prepare([5, 4, 3, 2, 1], ABC, {"none": ["x"]}, 1, 4)
    assert prep.names == ()
    assert prep.members == ()
    assert len(prep.sizes) == 0
    assert prep.dropped == ("none",)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_prepare_non_finite_raises(bad: float) -> None:
    with pytest.raises(ValueError):
        prepare([1.0, bad, 0.5], ["a", "b", "c"], {"s": ["a"]}, 1, 2)


def test_prepare_duplicate_labels_raise() -> None:
    with pytest.raises(ValueError):
        prepare([1.0, 2.0, 3.0], ["a", "b", "a"], {"s": ["a"]}, 1, 2)


def test_prepare_returns_prepared_and_does_not_mutate_inputs() -> None:
    stats = np.array([3.0, -1.0, 2.0, -2.0])
    labels = ["a", "b", "c", "d"]
    sets = {"s": ["d", "a", "zz"], "t": ["b", "c"]}
    stats0, labels0, sets0 = stats.copy(), list(labels), copy.deepcopy(sets)
    prep = prepare(stats, labels, sets, 1, 3)
    assert isinstance(prep, Prepared)
    np.testing.assert_array_equal(stats, stats0)
    assert labels == labels0
    assert sets == sets0


def test_gsea_es_does_not_mutate_inputs() -> None:
    stats = np.array([4.0, 3, 2, 1, 1, 1])
    members = np.array([2, 3, 6], dtype=np.int64)
    s0, m0 = stats.copy(), members.copy()
    gsea_es(stats, members)
    np.testing.assert_array_equal(stats, s0)
    np.testing.assert_array_equal(members, m0)


# ---------------------------------------------------------------------------------------------
# Agreement with a slow pure-Python transcription of calcGseaStat on 200 random inputs
# ---------------------------------------------------------------------------------------------


def _random_case(rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    n = int(rng.integers(2, 60))
    m = int(rng.integers(1, n))
    kind = int(rng.integers(0, 4))
    if kind == 0:
        x = rng.normal(size=n)
    elif kind == 1:  # ties
        x = np.round(rng.normal(size=n), 1)
    elif kind == 2:  # many zeros
        x = rng.normal(size=n) * (rng.random(n) < 0.3)
    else:  # heavy tail
        x = rng.standard_cauchy(size=n)
    r = np.abs(np.sort(x)[::-1])
    sel = np.sort(rng.choice(n, size=m, replace=False) + 1).astype(np.int64)
    return r, sel


def test_matches_reference_calc_gsea_stat_200_inputs() -> None:
    """200 random inputs (ties, zeros, heavy tails). Inputs whose ES is an exact rational tie
    (maxP == -minP before rounding) are left to the hand examples: float rounding decides them."""
    rng = np.random.default_rng(20261021)
    compared = 0
    for _ in range(200):
        r, sel = _random_case(rng)
        if is_exact_tie(r, sel):
            continue
        want_es, want_le = reference_gsea_stat(r, sel)
        es, le = gsea_es(r, sel)
        close(es, want_es)
        assert ints(le) == want_le, (r, sel)
        compared += 1
    assert compared >= 150
