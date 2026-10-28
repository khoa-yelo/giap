"""De novo precision per score bin, bootstrap CI and the de novo π̂ (issue M3-11).

GIPPA_SPEC.md §7.2.5 ("De novo precision, estimated from the genome"; step 4, m = 0 de novo),
§8.1, §9.3 (`mapping_fdr.tsv` de novo rows), §9.4 (`c_assumed`);
docs/specs/M3-11-denovo-precision.md
examples 1-9, edge cases and properties; M3-index lead decision 7 (drop unreportable classes within
a bin and renormalise); ADR-0010 (π̂ handed to M3-10 must be the conservative bound, 1 - ci_low).

Conventions fixed by these tests where the spec leaves them open:
- length-class keys are the M3-10 labels: "7", "8", "9", "10-11", "12-14", "15-19", "20-40";
- ``c_by_class`` reads ``peptide``, ``q_value``, ``is_decoy`` and ``target_exact`` (bool) of the
  reference peptides (``length_class`` is supplied too); the class comes from the peptide length;
- ``score_bins`` labels peptides without a score with a negative bin (no bin; π̂ = 1).
"""

from __future__ import annotations

import dataclasses
from collections.abc import Sequence

import numpy as np
import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.mapping.denovo_precision import (
    BinPrecision,
    c_by_class,
    denovo_precision,
    pi_wrong_denovo,
    score_bins,
)
from gippa.rng import for_step

CLASSES = [[7, 7], [8, 8], [9, 9], [10, 11], [12, 14], [15, 19], [20, 40]]
LABELS = ["7", "8", "9", "10-11", "12-14", "15-19", "20-40"]
SEED = 20261014


def _rng(seed: int = SEED) -> np.random.Generator:
    return for_step("placement", seed)


def _block(
    n: int, n_target: int, n_decoy: int, cls: str, score: float = 5.0, pep_len: int = 8
) -> pd.DataFrame:
    """n de novo peptides of one class: the first n_target match the target exactly, the last
    n_decoy match the shuffled decoy exactly (disjoint when n_target + n_decoy <= n)."""
    target = np.zeros(n, dtype=bool)
    target[:n_target] = True
    decoy = np.zeros(n, dtype=bool)
    if n_decoy:
        decoy[-n_decoy:] = True
    return pd.DataFrame(
        {
            "peptide": ["A" * pep_len] * n,
            "score": np.full(n, score, dtype=float),
            "length_class": [cls] * n,
            "target_exact": target,
            "decoy_exact": decoy,
        }
    )


def _one(rows: list[BinPrecision]) -> BinPrecision:
    assert len(rows) == 1, rows
    return rows[0]


# ---------------------------------------------------------------- interface


def test_bin_precision_is_frozen_dataclass() -> None:
    b = BinPrecision(
        score_bin=0, n=10, precision=0.5, ci_low=0.4, ci_high=0.6, c_value=0.99, flags=()
    )
    assert dataclasses.is_dataclass(b)
    with pytest.raises(dataclasses.FrozenInstanceError):
        b.precision = 0.7  # type: ignore[misc]
    assert [f.name for f in dataclasses.fields(b)] == [
        "score_bin",
        "n",
        "precision",
        "ci_low",
        "ci_high",
        "c_value",
        "flags",
    ]


# ---------------------------------------------------------------- examples 1-5 (estimator)


def test_example1_one_class() -> None:
    df = _block(100, 60, 5, "8")
    row = _one(denovo_precision(df, {"8": 0.99}, False, _rng()))
    assert row.score_bin == 0
    assert row.n == 100
    assert row.precision == pytest.approx(0.55 / 0.94, abs=1e-6)
    assert row.precision == pytest.approx(0.585106, abs=1e-6)
    assert row.c_value == pytest.approx(0.99)
    assert row.flags == ()
    assert row.ci_low is not None and row.ci_high is not None
    assert 0.0 <= row.ci_low <= row.precision <= row.ci_high <= 1.0


def test_example2_two_classes_weighted() -> None:
    df = pd.concat(
        [_block(140, 84, 7, "8"), _block(60, 24, 6, "10-11", pep_len=10)], ignore_index=True
    )
    row = _one(denovo_precision(df, {"8": 0.99, "10-11": 0.90}, False, _rng()))
    assert row.n == 200
    assert row.precision == pytest.approx(0.7 * 0.55 / 0.94 + 0.3 * 0.375, abs=1e-6)
    assert row.precision == pytest.approx(0.522074, abs=1e-6)
    # c_value = Σ_l w_bl c_l
    assert row.c_value == pytest.approx(0.7 * 0.99 + 0.3 * 0.90)


def test_example2_row_order_invariant() -> None:
    df = pd.concat(
        [_block(140, 84, 7, "8"), _block(60, 24, 6, "10-11", pep_len=10)], ignore_index=True
    )
    shuffled = df.sample(frac=1.0, random_state=3).reset_index(drop=True)
    c = {"8": 0.99, "10-11": 0.90}
    a = _one(denovo_precision(df, c, False, _rng()))
    b = _one(denovo_precision(shuffled, c, False, _rng()))
    assert a.precision == pytest.approx(b.precision, abs=1e-12)


def test_example3_negative_clipped_to_zero() -> None:
    row = _one(denovo_precision(_block(100, 2, 5, "8"), {"8": 0.99}, False, _rng()))
    assert row.precision == 0.0
    assert row.ci_low is not None and row.ci_low >= 0.0


def test_example4_above_one_clipped_to_one() -> None:
    row = _one(denovo_precision(_block(100, 99, 0, "8"), {"8": 0.95}, False, _rng()))
    assert row.precision == 1.0
    assert row.ci_high is not None and row.ci_high <= 1.0


def test_example5_unreportable_bin() -> None:
    df = _block(100, 30, 5, "8")
    row = _one(denovo_precision(df, {"8": 0.08}, False, _rng()))
    assert row.precision is None
    assert row.ci_low is None and row.ci_high is None
    assert "unreported" in row.flags
    assert row.n == 100
    pi = pi_wrong_denovo([row], np.zeros(100, dtype=int))
    np.testing.assert_array_equal(pi, np.ones(100))


def test_unreportable_class_dropped_and_weights_renormalised() -> None:
    # lead decision 7 (M3-11 Q1): class "10-11" has c - r = 0.08 - 0.05 < 0.05 → dropped; the bin
    # precision is the class-"8" estimate alone (weights renormalised to 1), not 0.7 times it.
    df = pd.concat(
        [_block(140, 84, 7, "8"), _block(60, 20, 3, "10-11", pep_len=10)], ignore_index=True
    )
    row = _one(denovo_precision(df, {"8": 0.99, "10-11": 0.08}, False, _rng()))
    assert row.precision == pytest.approx(0.55 / 0.94, abs=1e-6)
    assert "unreported" not in row.flags


def test_reportability_threshold_boundary() -> None:
    # c - r = 0.05 exactly is reportable ("< 0.05" is the exclusion rule): r = 0.05, c = 0.10
    df = _block(100, 8, 5, "8")
    row = _one(denovo_precision(df, {"8": 0.10}, False, _rng()))
    assert row.precision is not None
    assert row.precision == pytest.approx((0.08 - 0.05) / 0.05, abs=1e-6)
    # c - r = 0.04 → unreported
    row = _one(denovo_precision(df, {"8": 0.09}, False, _rng()))
    assert row.precision is None and "unreported" in row.flags


# ---------------------------------------------------------------- examples 6-7 (c_l)


def _ref(n: int, n_match: int, pep_len: int, q: float = 0.001, decoy: bool = False) -> pd.DataFrame:
    match = np.zeros(n, dtype=bool)
    match[:n_match] = True
    label = next(lab for (lo, hi), lab in zip(CLASSES, LABELS, strict=True) if lo <= pep_len <= hi)
    return pd.DataFrame(
        {
            "peptide": ["K" * pep_len] * n,
            "q_value": np.full(n, q),
            "is_decoy": np.full(n, decoy),
            "target_exact": match,
            "length_class": [label] * n,
        }
    )


def test_example6_no_reference_peptides() -> None:
    c, assumed = c_by_class(None, CLASSES)
    assert assumed is True
    assert set(c) == set(LABELS)
    assert all(v == pytest.approx(0.99) for v in c.values())
    df = pd.concat(
        [_block(50, 30, 2, "8", score=1.0), _block(50, 20, 2, "9", score=2.0, pep_len=9)],
        ignore_index=True,
    )
    rows = denovo_precision(df, c, assumed, _rng())
    assert rows
    assert all("c_assumed" in r.flags for r in rows)
    assert all(r.c_value == pytest.approx(0.99) for r in rows)


def test_example7_reference_class_8() -> None:
    c, assumed = c_by_class(_ref(200, 198, 8), CLASSES)
    assert assumed is False
    assert c["8"] == pytest.approx(0.99)
    rows = denovo_precision(_block(100, 60, 5, "8"), c, assumed, _rng())
    assert all("c_assumed" not in r.flags for r in rows)


def test_reference_excludes_high_q_and_decoys() -> None:
    ref = pd.concat(
        [
            _ref(200, 198, 8),
            _ref(100, 0, 8, q=0.05),  # q > 0.01 → excluded
            _ref(100, 0, 8, decoy=True),  # decoy → excluded
            _ref(100, 100, 9),
        ],
        ignore_index=True,
    )
    c, assumed = c_by_class(ref, CLASSES)
    assert assumed is False
    assert c["8"] == pytest.approx(0.99)
    assert c["9"] == pytest.approx(1.0)


def test_reference_q_boundary_included() -> None:
    ref = pd.concat([_ref(100, 100, 8), _ref(100, 0, 8, q=0.01)], ignore_index=True)
    c, _ = c_by_class(ref, CLASSES)
    assert c["8"] == pytest.approx(0.5)  # q ≤ 0.01 is inclusive


def test_reference_class_by_length_range() -> None:
    ref = pd.concat([_ref(10, 9, 10), _ref(10, 7, 11)], ignore_index=True)
    c, _ = c_by_class(ref, CLASSES)
    assert c["10-11"] == pytest.approx(16 / 20)


# ---------------------------------------------------------------- example 8 (π̂, ADR-0010)


def test_example8_pi_wrong_is_one_minus_lower_bound() -> None:
    b = BinPrecision(
        score_bin=0, n=50, precision=0.9, ci_low=0.85, ci_high=0.94, c_value=0.99, flags=()
    )
    pi = pi_wrong_denovo([b], np.array([0, 0, 0]))
    np.testing.assert_allclose(pi, [0.15, 0.15, 0.15])


def test_pi_wrong_uses_each_peptides_bin_and_missing_is_one() -> None:
    bins = [
        BinPrecision(0, 30, 0.5, 0.4, 0.6, 0.99, ()),
        BinPrecision(1, 30, None, None, None, 0.99, ("unreported",)),
        BinPrecision(2, 30, 0.9, 0.8, 0.95, 0.99, ()),
    ]
    # bin 3 has no row, -1 = no score → π̂ = 1
    pi = pi_wrong_denovo(bins, np.array([2, 0, 1, 3, -1, 2]))
    np.testing.assert_allclose(pi, [0.2, 0.6, 1.0, 1.0, 1.0, 0.2])


def test_pi_wrong_is_conservative_end_of_ci() -> None:
    # ADR-0010: π̂ = 0.5π breaks placement FDR, so the π̂ handed to M3-10 is 1 - ci_low, which is
    # never below 1 - precision (nor 1 - ci_high).
    rng = np.random.default_rng(1)
    df = pd.concat(
        [
            _block(n, int(n * m), int(n * 0.05), "8", score=float(s))
            for s, (n, m) in enumerate([(40, 0.3), (60, 0.5), (80, 0.7), (100, 0.9), (30, 0.95)])
        ],
        ignore_index=True,
    )
    df = df.sample(frac=1.0, random_state=rng.integers(1 << 31)).reset_index(drop=True)
    rows = denovo_precision(df, {"8": 0.99}, False, _rng())
    bin_of = score_bins(df["score"].to_numpy())
    pi = pi_wrong_denovo(rows, bin_of)
    by_bin = {r.score_bin: r for r in rows}
    for i, b in enumerate(bin_of):
        r = by_bin[int(b)]
        assert r.precision is not None and r.ci_low is not None
        assert pi[i] == pytest.approx(1.0 - r.ci_low)
        assert pi[i] >= 1.0 - r.precision - 1e-12
    assert np.all((pi >= 0.0) & (pi <= 1.0))


# ---------------------------------------------------------------- example 9 and bins


def test_example9_deciles_of_1_to_100() -> None:
    bins = score_bins(np.arange(1, 101, dtype=float))
    assert bins.shape == (100,)
    np.testing.assert_array_equal(bins, np.repeat(np.arange(10), 10))
    assert np.all(bins[:10] == 0)


def test_bins_right_closed_on_quantile_edges() -> None:
    scores = np.arange(1, 101, dtype=float)
    edges = np.quantile(scores, np.arange(1, 10) / 10)
    bins = score_bins(scores)
    for k, e in enumerate(edges):
        assert bins[scores <= e].max() <= k
        assert bins[scores > e].min() >= k + 1


def test_all_scores_equal_one_bin() -> None:
    bins = score_bins(np.full(37, 4.2))
    np.testing.assert_array_equal(bins, np.zeros(37, dtype=int))
    df = _block(37, 20, 1, "8", score=4.2)
    rows = denovo_precision(df, {"8": 0.99}, False, _rng())
    assert len(rows) == 1 and rows[0].n == 37


def test_tied_scores_share_a_bin_and_fewer_bins() -> None:
    scores = np.r_[np.full(50, 1.0), np.arange(2, 52, dtype=float)]
    bins = score_bins(scores)
    assert len(set(bins[:50].tolist())) == 1
    assert len(np.unique(bins)) < 10


def test_missing_score_forms_no_bin() -> None:
    scores = np.r_[np.arange(1, 21, dtype=float), [np.nan, np.nan]]
    bins = score_bins(scores)
    assert np.all(bins[-2:] < 0)
    assert np.all(bins[:-2] >= 0)
    df = pd.concat([_block(20, 12, 1, "8"), _block(2, 2, 0, "8")], ignore_index=True).assign(
        score=scores
    )
    rows = denovo_precision(df, {"8": 0.99}, False, _rng())
    assert sum(r.n for r in rows) == 20
    pi = pi_wrong_denovo(rows, bins)
    np.testing.assert_array_equal(pi[-2:], [1.0, 1.0])


def test_bins_reported_in_order_and_n_sums() -> None:
    df = pd.concat(
        [_block(30, 15, 1, "8", score=float(s)) for s in range(1, 11)], ignore_index=True
    )
    rows = denovo_precision(df, {"8": 0.99}, False, _rng())
    assert [r.score_bin for r in rows] == list(range(10))
    assert all(r.n == 30 for r in rows)


# ---------------------------------------------------------------- edge cases


def test_empty_input() -> None:
    df = _block(0, 0, 0, "8")
    assert denovo_precision(df, {"8": 0.99}, False, _rng()) == []
    assert score_bins(np.array([], dtype=float)).shape == (0,)
    assert pi_wrong_denovo([], np.array([], dtype=int)).shape == (0,)


def test_all_correct() -> None:
    row = _one(denovo_precision(_block(50, 50, 0, "8"), {"8": 0.99}, False, _rng()))
    assert row.precision == 1.0 and row.ci_low == 1.0 and row.ci_high == 1.0
    np.testing.assert_allclose(pi_wrong_denovo([row], np.zeros(3, dtype=int)), 0.0, atol=1e-12)


def test_all_wrong() -> None:
    row = _one(denovo_precision(_block(50, 0, 0, "8"), {"8": 0.99}, False, _rng()))
    assert row.precision == 0.0 and row.ci_low == 0.0 and row.ci_high == 0.0
    np.testing.assert_allclose(pi_wrong_denovo([row], np.zeros(3, dtype=int)), 1.0)


def test_single_peptide_bin_reported_degenerate_ci() -> None:
    row = _one(denovo_precision(_block(1, 1, 0, "8"), {"8": 0.99}, False, _rng()))
    assert row.n == 1
    assert row.precision == 1.0
    assert row.ci_low == row.ci_high == 1.0
    assert "unreported" not in row.flags


def test_exact_columns_only() -> None:
    # tolerant-hit columns, if present, are ignored: m and r use target_exact/decoy_exact only
    df = _block(100, 60, 5, "8").assign(target_tolerant=True, decoy_tolerant=True, n_diff=1)
    row = _one(denovo_precision(df, {"8": 0.99}, False, _rng()))
    assert row.precision == pytest.approx(0.585106, abs=1e-6)


# ---------------------------------------------------------------- bootstrap and determinism


def test_bootstrap_deterministic_with_for_step() -> None:
    df = pd.concat([_block(40, 20, 2, "8", score=float(s)) for s in range(5)], ignore_index=True)
    a = denovo_precision(df, {"8": 0.99}, False, for_step("placement", 11))
    b = denovo_precision(df, {"8": 0.99}, False, for_step("placement", 11))
    assert a == b


def test_bootstrap_depends_on_seed() -> None:
    df = _block(80, 40, 4, "8")
    a = _one(denovo_precision(df, {"8": 0.99}, False, for_step("placement", 1)))
    b = _one(denovo_precision(df, {"8": 0.99}, False, for_step("placement", 2)))
    assert a.precision == b.precision
    assert (a.ci_low, a.ci_high) != (b.ci_low, b.ci_high)


def test_bootstrap_ci_width_matches_binomial_se() -> None:
    # one class, r = 0: prec = m / c; percentile 95 % CI ≈ ± 1.96 sqrt(m(1-m)/n) / c
    n, m, c = 400, 0.5, 0.99
    row = _one(denovo_precision(_block(n, int(n * m), 0, "8"), {"8": c}, False, _rng()))
    assert row.ci_low is not None and row.ci_high is not None
    half = 1.96 * np.sqrt(m * (1 - m) / n) / c
    assert (row.ci_high - row.ci_low) / 2 == pytest.approx(half, rel=0.2)


def test_n_boot_respected_ci_from_resamples() -> None:
    # with n_boot = 1 the "interval" is one resample: lo == hi
    row = _one(denovo_precision(_block(60, 30, 3, "8"), {"8": 0.99}, False, _rng(), n_boot=1))
    assert row.ci_low == row.ci_high


def test_default_n_boot_is_1000() -> None:
    import inspect

    assert inspect.signature(denovo_precision).parameters["n_boot"].default == 1000
    assert inspect.signature(score_bins).parameters["n_bins"].default == 10


# ---------------------------------------------------------------- Hypothesis properties


@st.composite
def _frames(draw: st.DrawFn) -> tuple[pd.DataFrame, dict[str, float]]:
    n = draw(st.integers(1, 120))
    cls = draw(st.lists(st.sampled_from(["8", "10-11"]), min_size=n, max_size=n))
    scores = draw(st.lists(st.sampled_from([1.0, 2.0, 3.0, 4.5, 7.0, 9.0]), min_size=n, max_size=n))
    tgt = draw(st.lists(st.booleans(), min_size=n, max_size=n))
    dec = draw(st.lists(st.booleans(), min_size=n, max_size=n))
    c = {
        "8": draw(st.floats(0.05, 1.0)),
        "10-11": draw(st.floats(0.05, 1.0)),
    }
    df = pd.DataFrame(
        {
            "peptide": ["PEPTIDEK"] * n,
            "score": scores,
            "length_class": cls,
            "target_exact": tgt,
            "decoy_exact": dec,
        }
    )
    return df, c


def _check_rows(rows: Sequence[BinPrecision], n_total: int) -> None:
    assert sum(r.n for r in rows) == n_total
    for r in rows:
        if r.precision is None:
            assert "unreported" in r.flags and r.ci_low is None and r.ci_high is None
            continue
        assert 0.0 <= r.precision <= 1.0
        assert r.ci_low is not None and r.ci_high is not None
        assert 0.0 <= r.ci_low <= r.ci_high <= 1.0


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(_frames(), st.integers(0, 2**31 - 1))
def test_property_precision_in_unit_interval(
    data: tuple[pd.DataFrame, dict[str, float]], seed: int
) -> None:
    df, c = data
    rows = denovo_precision(df, c, False, for_step("placement", seed), n_boot=50)
    _check_rows(rows, len(df))
    pi = pi_wrong_denovo(rows, score_bins(df["score"].to_numpy()))
    assert pi.shape == (len(df),)
    assert np.all((pi >= 0.0) & (pi <= 1.0))


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(_frames(), st.integers(0, 2**31 - 1))
def test_property_row_order_invariant(
    data: tuple[pd.DataFrame, dict[str, float]], seed: int
) -> None:
    df, c = data
    perm = np.random.default_rng(seed).permutation(len(df))
    a = denovo_precision(df, c, False, for_step("placement", 0), n_boot=20)
    b = denovo_precision(
        df.iloc[perm].reset_index(drop=True), c, False, for_step("placement", 0), n_boot=20
    )
    assert [(r.score_bin, r.n, r.flags) for r in a] == [(r.score_bin, r.n, r.flags) for r in b]
    for x, y in zip(a, b, strict=True):
        if x.precision is None:
            assert y.precision is None
        else:
            assert y.precision == pytest.approx(x.precision, abs=1e-12)
            assert y.c_value == pytest.approx(x.c_value, abs=1e-12)


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    st.integers(20, 200),
    st.floats(0.0, 1.0),
    st.floats(0.0, 0.3),
    st.floats(0.5, 1.0),
    st.integers(0, 2**31 - 1),
)
def test_property_ci_contains_estimate_n_ge_20(
    n: int, m: float, r: float, c: float, seed: int
) -> None:
    n_t = round(n * m)
    n_d = min(round(n * r), n - n_t)
    df = _block(n, n_t, n_d, "8")
    row = _one(denovo_precision(df, {"8": c}, False, for_step("placement", seed), n_boot=200))
    if row.precision is None:
        return
    assert row.ci_low is not None and row.ci_high is not None
    assert row.ci_low - 1e-12 <= row.precision <= row.ci_high + 1e-12


@settings(max_examples=60, deadline=None)
@given(
    st.lists(
        st.one_of(st.floats(-1e6, 1e6, allow_nan=False), st.sampled_from([0.0, 1.0, 2.0])),
        min_size=1,
        max_size=300,
    )
)
def test_property_score_bins_monotone_and_bounded(scores: list[float]) -> None:
    s = np.asarray(scores, dtype=float)
    bins = score_bins(s)
    assert bins.shape == s.shape
    assert bins.min() == 0 and bins.max() <= 9
    order = np.argsort(s, kind="stable")
    assert np.all(np.diff(bins[order]) >= 0)
    # equal scores share a bin
    for v in np.unique(s):
        assert len(set(bins[s == v].tolist())) == 1
