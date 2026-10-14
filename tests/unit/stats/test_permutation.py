"""Permutation p-values, circular shifts and Besag-Clifford (issue M4-03; GIAP_SPEC.md §7.3.5(b),
§8.2, §8.3, decision 14, Appendix B.4).

Conventions pinned here: the +1 rule p = (1 + b)/(1 + B); NaN nulls are dropped; two-sided p is
centred on the pooled mean (sum(null) + obs)/(B + 1) (lead decision 1, PR #23); ties (within
1e-12 * max(1, |x|)) count as exceedances; linear contigs keep only non-wrapping translations;
Besag-Clifford stops at the h-th null >= obs with p = h / L.
"""

from __future__ import annotations

import inspect
import re
from collections.abc import Callable, Iterable
from pathlib import Path

import numpy as np
import pytest

from giap.rng import for_step
from giap.stats.permutation import (
    apply_shift,
    besag_clifford,
    choose_shifts,
    circular_shifts,
    perm_p,
)

NAN = float("nan")
STATS_DIR = Path(__file__).resolve().parents[3] / "src" / "giap" / "stats"


def close(actual: object, expected: object, atol: float = 1e-12) -> None:
    np.testing.assert_allclose(
        np.asarray(actual, dtype=float), expected, rtol=0, atol=atol, equal_nan=True
    )


def ints(a: object) -> list[int]:
    arr = np.asarray(a)
    assert arr.dtype == np.int64
    assert arr.ndim == 1
    return [int(x) for x in arr]


def sequence_fn(
    values: Iterable[float],
) -> tuple[Callable[[np.random.Generator], float], list[int]]:
    """A stat_fn yielding ``values`` in order (then 0.0 forever); calls[0] counts the calls."""
    it = iter(values)
    calls = [0]

    def fn(rng: np.random.Generator) -> float:
        assert isinstance(rng, np.random.Generator)
        calls[0] += 1
        return float(next(it, 0.0))

    return fn, calls


# ---------------------------------------------------------------------------- interface (pinned)


def test_signatures_are_pinned() -> None:
    assert list(inspect.signature(perm_p).parameters) == ["obs", "null", "alternative"]
    assert inspect.signature(perm_p).parameters["alternative"].default == "two-sided"
    assert list(inspect.signature(circular_shifts).parameters) == ["n", "idx", "circular"]
    assert list(inspect.signature(apply_shift).parameters) == ["idx", "k", "n"]
    assert list(inspect.signature(choose_shifts).parameters) == [
        "valid",
        "n_max_exhaustive",
        "n",
        "B",
        "rng",
    ]
    assert list(inspect.signature(besag_clifford).parameters) == [
        "stat_fn",
        "obs",
        "h",
        "max_B",
        "rng",
    ]


# ---------------------------------------------------------------------------- perm_p examples


def test_example1_two_sided_centred_on_pooled_mean() -> None:
    # centre = (15 + 5.5)/6 = 3.41667; |obs dev| = 2.0833; null |dev| 2.4167 is the only one >= it
    p = perm_p(5.5, [1, 2, 3, 4, 5])
    assert isinstance(p, float)
    close(p, 2 / 6)


def test_example2_two_sided_counts_equal_deviations() -> None:
    close(perm_p(1.0, [1, 2, 3, 4, 5]), 3 / 6)


def test_example3_greater_and_less_count_ties() -> None:
    close(perm_p(4, [1, 2, 3, 4, 5], "greater"), 0.5)
    close(perm_p(4, [1, 2, 3, 4, 5], "less"), 5 / 6)


def test_example4_appendix_b4_value() -> None:
    null = np.linspace(-1.0, 1.0, 99)  # mean 0, every |dev| <= 1 < 5
    close(perm_p(5.0, null), 0.01)
    close(perm_p(5.0, null, "greater"), 0.01)
    close(perm_p(-5.0, null, "less"), 0.01)


def test_example4_appendix_b4_end_to_end_on_100_gene_replicon() -> None:
    n = 100
    x = np.zeros(n)
    x[[10, 11, 12]] = [3.0, 4.0, 5.0]  # region holds all the signal
    idx = np.array([10, 11, 12])
    shifts = circular_shifts(n, idx, circular=True)
    assert len(shifts) == 99  # K = 99
    shifts2, exhaustive = choose_shifts(shifts, 20_000, n, 1_000, for_step("context", 1))
    assert bool(exhaustive) is True
    assert ints(shifts2) == list(range(1, 100))
    obs = float(x[idx].sum())
    null = np.array([x[apply_shift(idx, int(k), n)].sum() for k in shifts2])
    centre = (null.sum() + obs) / (len(null) + 1)
    assert (np.abs(null - centre) < abs(obs - centre)).all()
    close(perm_p(obs, null), 0.01)


# ---------------------------------------------------------------------------- perm_p edges


def test_obs_beyond_everything_hits_floor_not_zero() -> None:
    p = perm_p(1e9, [0.0, 1.0, 2.0], "greater")
    close(p, 1 / 4)
    assert p > 0


def test_all_ties_give_one() -> None:
    for alt in ("two-sided", "greater", "less"):
        close(perm_p(3.0, [3.0, 3.0, 3.0], alt), 1.0)


def test_float_noise_counts_as_tie() -> None:
    a = 0.1 + 0.2 + 0.3  # 0.6000000000000001
    b = 0.3 + 0.2 + 0.1  # 0.6
    assert a != b
    close(perm_p(a, [b], "greater"), 1.0)
    close(perm_p(b, [a], "less"), 1.0)
    # two-sided: obs deviation equal to a null deviation up to float noise
    null = np.array([0.0, b, -a])  # centre = b/4 = 0.15; |dev| of b equals obs's up to noise
    close(perm_p(a, null), 3 / 4)


def test_tie_tolerance_is_relative_to_magnitude() -> None:
    big = 1e6
    # 1e-9 apart at magnitude 1e6 is within 1e-12 * 1e6 = 1e-6 -> tie
    close(perm_p(big + 1e-9, [big], "greater"), 1.0)
    # a real difference is not a tie
    close(perm_p(big + 1.0, [big], "greater"), 0.5)


def test_nan_null_entries_are_dropped() -> None:
    close(perm_p(4, [1, 2, NAN, 3, 4, 5, NAN], "greater"), 0.5)
    close(perm_p(5.5, [NAN, 1, 2, 3, 4, 5]), 2 / 6)


@pytest.mark.parametrize("alt", ["two-sided", "greater", "less"])
def test_nan_obs_or_empty_null_gives_nan(alt: str) -> None:
    assert np.isnan(perm_p(NAN, [1.0, 2.0], alt))  # type: ignore[arg-type]
    assert np.isnan(perm_p(1.0, [], alt))  # type: ignore[arg-type]
    assert np.isnan(perm_p(1.0, [NAN, NAN], alt))  # type: ignore[arg-type]


def test_single_null_value() -> None:
    close(perm_p(2.0, [1.0], "greater"), 0.5)
    close(perm_p(0.0, [1.0], "greater"), 1.0)
    # two-sided with one null: centre = (1 + 7)/2 = 4, both |dev| = 3 -> tie -> p = 1
    close(perm_p(7.0, [1.0]), 1.0)


def test_accepts_numpy_input() -> None:
    close(perm_p(np.float64(4), np.array([1, 2, 3, 4, 5]), "greater"), 0.5)


# ---------------------------------------------------------------------------- circular_shifts


def test_example5_circular_all_nonzero_shifts() -> None:
    assert ints(circular_shifts(10, [2, 3, 4], circular=True)) == list(range(1, 10))


def test_example6_linear_keeps_non_wrapping_translations() -> None:
    assert ints(circular_shifts(10, [2, 3, 4], circular=False)) == [1, 2, 3, 4, 5, 8, 9]


def test_example7_linear_spanning_region_has_no_shift() -> None:
    assert ints(circular_shifts(5, [0, 4], circular=False)) == []


@pytest.mark.parametrize("circular", [True, False])
def test_single_gene_replicon_has_no_shift(circular: bool) -> None:
    assert ints(circular_shifts(1, [0], circular=circular)) == []


def test_two_gene_replicon() -> None:
    assert ints(circular_shifts(2, [0], circular=True)) == [1]
    assert ints(circular_shifts(2, [0], circular=False)) == [1]  # 0 -> 1, no wrap
    assert ints(circular_shifts(2, [1], circular=False)) == [1]  # 1 -> 0, a left translation
    assert ints(circular_shifts(2, [0, 1], circular=False)) == []
    assert ints(circular_shifts(2, [0, 1], circular=True)) == [1]


def test_non_contiguous_set_translates_whole() -> None:
    # max 5 + k < 10 -> k in 1..4; min 1 + k >= 10 -> k = 9
    assert ints(circular_shifts(10, [1, 5], circular=False)) == [1, 2, 3, 4, 9]
    assert ints(circular_shifts(10, [5, 1], circular=False)) == [1, 2, 3, 4, 9]  # order-free
    assert ints(circular_shifts(10, [1, 5], circular=True)) == list(range(1, 10))


def test_region_at_end_of_linear_contig_shifts_left_only() -> None:
    assert ints(circular_shifts(6, [4, 5], circular=False)) == [2, 3, 4, 5]


def test_region_at_start_of_linear_contig_shifts_right_only() -> None:
    assert ints(circular_shifts(6, [0, 1], circular=False)) == [1, 2, 3, 4]


@pytest.mark.parametrize("circular", [True, False])
@pytest.mark.parametrize(
    ("n", "idx"),
    [(10, []), (10, [10]), (10, [-1]), (10, [3, 3]), (10, [0, 2, 2]), (0, [0])],
)
def test_invalid_idx_raises(n: int, idx: list[int], circular: bool) -> None:
    with pytest.raises(ValueError):
        circular_shifts(n, idx, circular=circular)


# ---------------------------------------------------------------------------- apply_shift


def test_example8_wraps_and_preserves_size() -> None:
    out = apply_shift([98, 99], 3, 100)
    assert ints(out) == [1, 2]


def test_apply_shift_is_elementwise_mod_n() -> None:
    assert ints(apply_shift([8, 9], 1, 10)) == [9, 0]
    assert ints(apply_shift([2, 3, 4], 8, 10)) == [0, 1, 2]
    assert ints(apply_shift(np.array([1, 5], dtype=np.int32), 4, 10)) == [5, 9]


# ---------------------------------------------------------------------------- choose_shifts


def test_exhaustive_when_n_at_most_threshold() -> None:
    valid = np.arange(1, 100)
    shifts, exhaustive = choose_shifts(valid, 100, 100, 10, for_step("context", 7))
    assert bool(exhaustive) is True
    assert ints(shifts) == list(range(1, 100))


def test_random_without_replacement_when_n_above_threshold() -> None:
    valid = np.arange(1, 1000)
    shifts, exhaustive = choose_shifts(valid, 500, 1000, 50, for_step("context", 7))
    assert bool(exhaustive) is False
    s = ints(shifts)
    assert len(s) == 50
    assert len(set(s)) == 50
    assert s == sorted(s)
    assert set(s) <= set(range(1, 1000))


def test_random_takes_all_when_B_exceeds_valid() -> None:
    valid = [1, 2, 3, 4, 5, 8, 9]
    shifts, _ = choose_shifts(valid, 5, 10, 1_000, for_step("context", 7))
    assert ints(shifts) == valid


def test_unsorted_valid_comes_back_sorted() -> None:
    shifts, _ = choose_shifts([9, 1, 5], 100, 10, 1_000, for_step("context", 7))
    assert ints(shifts) == [1, 5, 9]


@pytest.mark.parametrize(("n_max", "n"), [(100, 10), (5, 10)])
def test_empty_valid_gives_empty(n_max: int, n: int) -> None:
    shifts, _ = choose_shifts([], n_max, n, 100, for_step("context", 7))
    assert ints(shifts) == []


def test_choose_shifts_deterministic_via_for_step() -> None:
    valid = np.arange(1, 50_000)
    a, _ = choose_shifts(valid, 20_000, 50_000, 1_000, for_step("context.regions", 42))
    b, _ = choose_shifts(valid, 20_000, 50_000, 1_000, for_step("context.regions", 42))
    c, _ = choose_shifts(valid, 20_000, 50_000, 1_000, for_step("context.regions", 43))
    assert ints(a) == ints(b)
    assert ints(a) != ints(c)


# ---------------------------------------------------------------------------- besag_clifford


def test_example9_stops_at_h_th_exceedance() -> None:
    fn, calls = sequence_fn([1, 2, 6, 0, 0, 0, 7, 9, 9, 9])
    p, b_used = besag_clifford(fn, 5, 2, 100, for_step("bc", 1))
    close(p, 2 / 7)
    assert b_used == 7
    assert calls[0] == 7


def test_example10_never_exceeds_hits_max_B() -> None:
    fn, calls = sequence_fn([])
    p, b_used = besag_clifford(fn, 5, 2, 100, for_step("bc", 1))
    close(p, 1 / 101)
    assert b_used == 100
    assert calls[0] == 100


def test_fewer_than_h_exceedances_uses_plus_one_rule() -> None:
    fn, _ = sequence_fn([1, 2, 6, 0, 0])
    p, b_used = besag_clifford(fn, 5, 2, 5, for_step("bc", 1))
    close(p, (1 + 1) / (5 + 1))
    assert b_used == 5


def test_h_th_exceedance_on_last_allowed_draw_stops_sequentially() -> None:
    fn, _ = sequence_fn([1, 2, 6, 0, 0, 0, 7])
    p, b_used = besag_clifford(fn, 5, 2, 7, for_step("bc", 1))
    close(p, 2 / 7)
    assert b_used == 7


def test_null_equal_to_obs_counts_as_exceedance() -> None:
    fn, _ = sequence_fn([5, 5])
    p, b_used = besag_clifford(fn, 5, 2, 100, for_step("bc", 1))
    close(p, 1.0)
    assert b_used == 2


def test_h_one_stops_at_first_exceedance() -> None:
    fn, _ = sequence_fn([0, 0, 0, 10])
    p, b_used = besag_clifford(fn, 5, 1, 100, for_step("bc", 1))
    close(p, 1 / 4)
    assert b_used == 4


@pytest.mark.parametrize(("h", "max_B"), [(0, 10), (-1, 10), (11, 10)])
def test_invalid_h_or_max_B_raises(h: int, max_B: int) -> None:
    fn, _ = sequence_fn([])
    with pytest.raises(ValueError):
        besag_clifford(fn, 5, h, max_B, for_step("bc", 1))


def test_besag_clifford_passes_rng_and_is_deterministic() -> None:
    def fn(rng: np.random.Generator) -> float:
        return float(rng.normal())

    a = besag_clifford(fn, 1.5, 10, 1_000, for_step("context.bc", 9))
    b = besag_clifford(fn, 1.5, 10, 1_000, for_step("context.bc", 9))
    assert a == b
    p, b_used = a
    assert 0 < p <= 1
    assert 10 <= b_used <= 1_000


# ---------------------------------------------------------------------------- no global RNG


def test_permutation_module_exists_and_stats_use_no_global_rng() -> None:
    assert (STATS_DIR / "permutation.py").is_file()
    allowed = re.compile(r"np\.random\.Generator\b")
    for path in sorted(STATS_DIR.glob("*.py")):
        for i, line in enumerate(path.read_text().splitlines(), 1):
            stripped = allowed.sub("", line)
            assert "np.random." not in stripped, f"{path.name}:{i}: {line.strip()}"
            assert "numpy.random" not in stripped, f"{path.name}:{i}: {line.strip()}"
