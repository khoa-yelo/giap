"""M9-22 ``giap.compare.abundance.abundance_pairs``: pairwise Spearman rho with a TU-block
bootstrap CI and Bland-Altman bias and limits (docs/specs/M9-22-compare-tables.md "Tables",
"Degenerate cases", examples 2 and 7, acceptance "ties in Spearman; ddof"; GIAP_SPEC.md §7.9.3,
§8.2, §8.3).

File named ``test_compare_abundance.py``: ``tests/unit/test_abundance.py`` exists and tests/ has
no ``__init__.py``, so basenames must be unique.

Interface assumed (test-author reading of the spec):

- ``abundance_pairs(runs, tus, B, rng)``: ``runs`` maps label -> genes frame in label order;
  ``tus`` is ``tu_id`` indexed by ``locus_tag`` (the TUs of the first run, see
  ``helpers.compare_data.tus_of``); returns one row per unordered pair in label order with the
  ``compare_abundance`` columns. Abundance is ``abundance_log2``.
- Shared genes = ``detected`` in both and ``abundance_imputed`` false in both.
- The CI is ``giap.stats.bootstrap.tu_block_bootstrap`` over the shared genes (genes-table row
  order, ``tus`` as the TU column) with the Spearman statistic, drawing from ``rng`` (one pair:
  the oracle test below reproduces it exactly).
- ``flags`` is a ``;``-joined string (NA or empty when none); ``n_boot`` = ``B``.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from helpers.compare_data import ABUNDANCE_COLUMNS, flag_set, gene, genes_frame, tus_of
from hypothesis import assume, given, settings
from hypothesis import strategies as st
from scipy import stats as sps

from giap.rng import for_step

FAST = settings(max_examples=40, deadline=None)
B_SMALL = 50


def _pair(
    a: list[float | None],
    b: list[float | None],
    *,
    tus: list[str] | None = None,
    type_a: str = "maxlfq",
    type_b: str = "maxlfq",
    labels: tuple[str, str] = ("A", "B"),
) -> dict[str, pd.DataFrame]:
    tus = tus or [f"TU_{i}" for i in range(len(a))]

    def frame(vals: list[float | None], kind: str) -> pd.DataFrame:
        return genes_frame([gene(f"g{i:03d}", start=100 * i, log2=v, tu=tus[i],
                                 abundance_type=kind) for i, v in enumerate(vals)])  # fmt: skip

    return {labels[0]: frame(a, type_a), labels[1]: frame(b, type_b)}


def _run(runs: dict[str, pd.DataFrame], B: int = B_SMALL, seed: int = 7) -> pd.DataFrame:
    from giap.compare.abundance import abundance_pairs

    first = next(iter(runs.values()))
    return abundance_pairs(runs, tus_of(first), B, for_step("compare", seed))


def _row(runs: dict[str, pd.DataFrame], **kw: int) -> pd.Series:
    out = _run(runs, **kw)
    assert len(out) == 1
    return out.iloc[0]


def _na(x: object) -> bool:
    return x is None or x is pd.NA or (isinstance(x, float) and math.isnan(x))


# ------------------------------------------------------------------- example 2


def test_example2_rho_bias_and_limits() -> None:
    row = _row(_pair([1, 2, 3, 4, 5], [2, 1, 4, 3, 5]))
    assert set(ABUNDANCE_COLUMNS) <= set(row.index)
    assert row["run_a"] == "A" and row["run_b"] == "B"
    assert int(row["n_shared"]) == 5
    assert row["spearman_rho"] == pytest.approx(0.8)
    assert row["ba_bias_log2"] == pytest.approx(0.0, abs=1e-12)
    assert row["ba_loa_low"] == pytest.approx(-1.96)  # SD 1 needs ddof = 1
    assert row["ba_loa_high"] == pytest.approx(1.96)
    assert row["method"] == "spearman_tu_bootstrap"
    assert int(row["n_boot"]) == B_SMALL


def test_example2_ci_is_ordered_and_in_range() -> None:
    row = _row(_pair([1, 2, 3, 4, 5], [2, 1, 4, 3, 5]))
    assert -1.0 <= row["ci_low"] <= row["ci_high"] <= 1.0


def test_bias_sign_is_a_minus_b() -> None:
    row = _row(_pair([3, 4, 5, 6], [1, 2, 3, 4.5]))
    d = np.array([2, 2, 2, 1.5])
    assert row["ba_bias_log2"] == pytest.approx(d.mean())
    sd = d.std(ddof=1)
    assert row["ba_loa_low"] == pytest.approx(d.mean() - 1.96 * sd)
    assert row["ba_loa_high"] == pytest.approx(d.mean() + 1.96 * sd)
    assert row["ba_loa_high"] != pytest.approx(d.mean() + 1.96 * d.std(ddof=0))


def test_spearman_ties_use_average_ranks() -> None:
    a = [1.0, 2.0, 2.0, 3.0, 4.0, 5.0, 5.0]
    b = [1.0, 3.0, 2.0, 2.0, 5.0, 4.0, 4.0]
    row = _row(_pair(a, b))
    expected = sps.spearmanr(a, b).statistic
    assert row["spearman_rho"] == pytest.approx(expected)
    # ordinal (first-occurrence) ranks would give a different value
    ordinal = sps.pearsonr(sps.rankdata(a, "ordinal"), sps.rankdata(b, "ordinal")).statistic
    assert row["spearman_rho"] != pytest.approx(ordinal)


# ------------------------------------------------------------- shared genes


def test_imputed_and_undetected_genes_are_not_shared() -> None:
    a = genes_frame([gene(f"g{i}", start=100 * i, log2=v) for i, v in enumerate([1, 2, 3, 4, 5])]
                    + [gene("g5", start=500, log2=9.0, imputed=True),
                       gene("g6", start=600, log2=9.0),
                       gene("g7", start=700, log2=None, detected=False)])  # fmt: skip
    b = genes_frame([gene(f"g{i}", start=100 * i, log2=v) for i, v in enumerate([2, 1, 4, 3, 5])]
                    + [gene("g5", start=500, log2=-9.0),
                       gene("g6", start=600, log2=-9.0, detected=False),
                       gene("g7", start=700, log2=-9.0)])  # fmt: skip
    row = _row({"A": a, "B": b})
    assert int(row["n_shared"]) == 5
    assert row["spearman_rho"] == pytest.approx(0.8)


def test_two_shared_genes_is_too_few() -> None:
    row = _row(_pair([1, 2], [2, 1]))
    assert int(row["n_shared"]) == 2
    for col in ("spearman_rho", "ci_low", "ci_high", "ba_bias_log2", "ba_loa_low", "ba_loa_high"):
        assert _na(row[col]), col
    assert "too_few_shared" in flag_set(row["flags"])


def test_three_shared_genes_is_enough() -> None:
    row = _row(_pair([1, 2, 3], [1, 3, 2]))
    assert row["spearman_rho"] == pytest.approx(0.5)
    assert "too_few_shared" not in flag_set(row["flags"])


def test_one_tu_gives_ci_few_units() -> None:
    row = _row(_pair([1, 2, 3, 4], [1, 3, 2, 4], tus=["TU_1"] * 4))
    assert row["spearman_rho"] == pytest.approx(sps.spearmanr([1, 2, 3, 4], [1, 3, 2, 4])[0])
    assert _na(row["ci_low"]) and _na(row["ci_high"])
    assert "ci_few_units" in flag_set(row["flags"])
    assert not _na(row["ba_bias_log2"])


def test_two_tus_give_a_ci_flagged_few_units() -> None:
    """Review B1 on #232 (decision 11): 2-9 units keep the CI but flag ``ci_few_units``."""
    row = _row(_pair([1, 2, 3, 4, 5, 6], [1, 3, 2, 4, 6, 5], tus=["T1"] * 3 + ["T2"] * 3))
    assert "ci_few_units" in flag_set(row["flags"])
    assert row["n_units"] == 2
    assert math.isfinite(row["ci_low"]) and math.isfinite(row["ci_high"])


def test_ten_tus_are_not_few_units() -> None:
    a = list(range(1, 21))
    b = [x + (0.5 if x % 2 else -0.5) for x in a]
    row = _row(_pair(a, b, tus=[f"T{i // 2}" for i in range(20)]))
    assert row["n_units"] == 10
    assert "ci_few_units" not in flag_set(row["flags"])


# --------------------------------------------------------- scales and missing


def test_example7_scale_differs_keeps_ba_columns() -> None:
    row = _row(_pair([1, 2, 3, 4, 5], [2, 1, 4, 3, 5], type_a="maxlfq", type_b="ibaq"))
    assert "scale_differs" in flag_set(row["flags"])
    assert row["ba_bias_log2"] == pytest.approx(0.0, abs=1e-12)
    assert row["ba_loa_high"] == pytest.approx(1.96)
    assert row["spearman_rho"] == pytest.approx(0.8)


def test_same_type_different_scale_is_flagged() -> None:
    runs = _pair([1, 2, 3, 4], [1, 3, 2, 4])
    runs["B"]["abundance_scale"] = "molar"
    assert "scale_differs" in flag_set(_row(runs)["flags"])


def test_same_type_and_scale_is_not_flagged() -> None:
    assert "scale_differs" not in flag_set(_row(_pair([1, 2, 3, 4], [1, 3, 2, 4]))["flags"])


def test_run_without_abundance_gives_na_row_missing_in() -> None:
    runs = _pair([1, 2, 3, 4, 5], [None] * 5, type_b="none")
    row = _row(runs)
    for col in ("spearman_rho", "ci_low", "ci_high", "ba_bias_log2", "ba_loa_low", "ba_loa_high"):
        assert _na(row[col]), col
    assert "missing_in:B" in flag_set(row["flags"])


# ------------------------------------------------------------- pairs, seeds


def test_three_runs_give_three_pairs_in_label_order() -> None:
    vals = [1.0, 2.0, 3.0, 4.0, 5.0]
    runs = {lab: genes_frame([gene(f"g{i}", start=100 * i, log2=v + k * (i % 2))
                              for i, v in enumerate(vals)])
            for k, lab in enumerate(["Z", "A", "M"])}  # fmt: skip
    out = _run(runs)
    assert list(zip(out["run_a"], out["run_b"], strict=True)) == [("Z", "A"), ("Z", "M"),
                                                                   ("A", "M")]  # fmt: skip


def test_example6_same_seed_same_rows() -> None:
    runs = _pair(list(np.linspace(0, 10, 20)), list(np.linspace(0, 10, 20)[::-1] ** 0.5),
                 tus=[f"T{i // 2}" for i in range(20)])  # fmt: skip
    pd.testing.assert_frame_equal(_run(runs, seed=11), _run(runs, seed=11))


def test_ci_is_the_tu_block_bootstrap_of_spearman() -> None:
    """Genes are carried with their TU: the CI equals ``tu_block_bootstrap`` on the shared genes
    with the first run's TUs and the same generator."""
    from giap.stats.bootstrap import tu_block_bootstrap

    gen = np.random.default_rng(3)
    n = 30
    a = list(gen.normal(size=n))
    b = list(np.asarray(a) + gen.normal(scale=0.7, size=n))
    tus = [f"TU_{i // 3}" for i in range(n)]
    row = _row(_pair(a, b, tus=tus), B=200, seed=5)
    df = pd.DataFrame({"tu": tus, "a": a, "b": b})

    def rho(d: pd.DataFrame) -> float:
        return float(sps.spearmanr(d["a"], d["b"]).statistic)

    ci = tu_block_bootstrap(df, "tu", rho, B=200, rng=for_step("compare", 5))
    assert row["spearman_rho"] == pytest.approx(ci.estimate)
    assert row["ci_low"] == pytest.approx(ci.ci_low)
    assert row["ci_high"] == pytest.approx(ci.ci_high)


# ------------------------------------------------------------------ properties

_vals = st.floats(min_value=-20, max_value=40, allow_nan=False, allow_infinity=False)


@st.composite
def _pairs(draw: st.DrawFn) -> tuple[list[float], list[float]]:
    n = draw(st.integers(min_value=3, max_value=25))
    a = draw(st.lists(_vals, min_size=n, max_size=n))
    b = draw(st.lists(_vals, min_size=n, max_size=n))
    return a, b


@FAST
@given(_pairs())
def test_property_matches_scipy_and_bland_altman_formula(ab: tuple[list[float], ...]) -> None:
    a, b = ab
    assume(len(set(a)) > 1 and len(set(b)) > 1)
    row = _row(_pair(list(a), list(b)), B=20)
    assert row["spearman_rho"] == pytest.approx(sps.spearmanr(a, b).statistic, abs=1e-9)
    d = np.asarray(a) - np.asarray(b)
    assert row["ba_bias_log2"] == pytest.approx(d.mean(), abs=1e-9)
    half = 1.96 * d.std(ddof=1)
    assert row["ba_loa_low"] == pytest.approx(d.mean() - half, abs=1e-9)
    assert row["ba_loa_high"] == pytest.approx(d.mean() + half, abs=1e-9)
    assert -1.0 - 1e-12 <= row["spearman_rho"] <= 1.0 + 1e-12


@FAST
@given(_pairs())
def test_property_swapping_runs_negates_bias_keeps_rho(ab: tuple[list[float], ...]) -> None:
    a, b = ab
    assume(len(set(a)) > 1 and len(set(b)) > 1)
    fwd = _row(_pair(list(a), list(b)), B=20)
    rev = _row(_pair(list(b), list(a)), B=20)
    assert rev["spearman_rho"] == pytest.approx(fwd["spearman_rho"], abs=1e-9)
    assert rev["ba_bias_log2"] == pytest.approx(-fwd["ba_bias_log2"], abs=1e-9)
    assert rev["ba_loa_low"] == pytest.approx(-fwd["ba_loa_high"], abs=1e-9)


@FAST
@given(_pairs(), st.floats(min_value=0.1, max_value=5), st.floats(min_value=-5, max_value=5))
def test_property_rho_invariant_to_monotone_rescaling(
    ab: tuple[list[float], ...], scale: float, shift: float
) -> None:
    """A run on another log scale (``scale * x + shift``) keeps rho; the bias moves by the
    shift only when ``scale == 1`` (checked with the pure shift below)."""
    a, b = ab
    assume(len(set(a)) > 1 and len(set(b)) > 1)
    base = _row(_pair(list(a), list(b)), B=20)
    moved = _row(_pair([scale * x + shift for x in a], list(b)), B=20)
    assume(len({scale * x + shift for x in a}) == len(set(a)))  # no float collisions
    assert moved["spearman_rho"] == pytest.approx(base["spearman_rho"], abs=1e-9)
    shifted = _row(_pair([x + shift for x in a], list(b)), B=20)
    assert shifted["ba_bias_log2"] == pytest.approx(base["ba_bias_log2"] + shift, abs=1e-6)
    width = base["ba_loa_high"] - base["ba_loa_low"]
    assert shifted["ba_loa_high"] - shifted["ba_loa_low"] == pytest.approx(width, abs=1e-6)


@FAST
@given(_pairs(), st.integers(min_value=0, max_value=2**31 - 1))
def test_property_ci_ordered_and_bounded(ab: tuple[list[float], ...], seed: int) -> None:
    a, b = ab
    assume(len(set(a)) > 1 and len(set(b)) > 1)
    n = len(a)
    row = _row(_pair(list(a), list(b), tus=[f"T{i // 2}" for i in range(n)]), B=20, seed=seed)
    if not _na(row["ci_low"]):
        assert -1.0 - 1e-12 <= row["ci_low"] <= row["ci_high"] <= 1.0 + 1e-12


def test_loa_contain_about_95_percent_of_normal_differences() -> None:
    """Bland-Altman claim: under normal differences the limits hold ~95 % of them
    (n = 4000; binomial SD 0.0034, band +- 0.015)."""
    gen = for_step("test.compare.loa", 20261022)
    n = 4000
    a = gen.normal(10.0, 2.0, n)
    b = a - gen.normal(0.3, 0.5, n)
    row = _row(_pair(list(a), list(b)), B=20)
    d = a - b
    inside = np.mean((d >= row["ba_loa_low"]) & (d <= row["ba_loa_high"]))
    assert abs(inside - 0.95) < 0.015
    assert row["ba_bias_log2"] == pytest.approx(0.3, abs=0.05)
