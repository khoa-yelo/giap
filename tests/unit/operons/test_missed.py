"""Missed operons: NO-pair scores, CV and SW decoys, q-values and local fdr (issue M6-06; §7.4.4,
§8.1 row "Missed operon", §9.3 ``pairs.tsv``).

Spec examples 1-9 and the "Degenerate cases" table of docs/specs/M6-06-missed-operons.md on
hand-built pair frames (tests/helpers/missed_data.py), with lead decisions 1-5 (2026-10-21):
``p_value`` NA, imputed abundances stay among the targets, the minimum is 20 scored NO pairs, and
(M6-01 decision 12) the family is near-powerless at realistic effects, so nothing here asserts
power at 0.25-1 SD: example 6 plants |Delta| <= 1e-3 SD.

Readings this file pins (the spec leaves them open):

- ``rows`` holds NO pairs only (decoys are never tested); a scored target has
  ``test = "missed_operon_decoy"``, and the degenerate paths leave NO rows with ``test = "none"``
  carrying the new flag, as M6-05 does for ``too_few_op_pairs``.
- ``rows`` has a ``which`` column (the decoy set that decided the row's FDR) since examples 2 and 8
  assert it, besides the §9.2 test columns of spec step 4.
- ``bin_counts`` has one row per distance bin, in bin order, with columns ``N_NO``, ``N_cv``,
  ``N_sw``.
- Eligibility (``partial``), imputation and the relative scale live in the gene table, which the
  spec's ``missed_operons(pairs, cfg)`` does not take: it is passed as the keyword ``genes=``
  (absent = every gene eligible, molar scale).
- The performance target (20,000 + 20,000 pairs < 2 s) belongs to benchmarks/ (M5-02 decision 12).
"""

from __future__ import annotations

import logging
import math

import numpy as np
import pandas as pd
import pytest
from helpers.missed_data import (
    CV,
    EX1,
    EX2,
    SW,
    TARGETS,
    block,
    cfg,
    example_frame,
    genes_for,
    levels,
    null_frame,
    stack,
)
from helpers.pairs_data import BASE, frame, genes_of, ring
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.operons.missed import MIN_NO_PAIRS, MissedOperons, missed_operons

FAST = settings(max_examples=30, deadline=None)
TEST = "missed_operon_decoy"
ROW_COLS = (
    "pair_id", "test", "stat", "effect", "ci_low", "ci_high", "p_value", "q_value", "method",
    "null_model", "n_perm", "lfdr", "which", "flags",
)  # fmt: skip
CV_ONLY = {"decoys": ["convergent"]}


def scored(res: MissedOperons) -> pd.DataFrame:
    return res.rows[res.rows["test"] == TEST]


def flagged(res: MissedOperons, flag: str) -> bool:
    """Every row carries ``flag`` (and there is at least one row to carry it)."""
    return len(res.rows) > 0 and all(flag in list(f) for f in res.rows["flags"])


def by_id(res: MissedOperons) -> pd.DataFrame:
    return scored(res).set_index("pair_id").sort_index()


def q_at(res: MissedOperons, pairs: pd.DataFrame, abs_delta: float) -> np.ndarray:
    ids = pairs.loc[(pairs["class"] == "NO") & np.isclose(pairs["delta_log2"].abs(), abs_delta)]
    return by_id(res).loc[ids["pair_id"], "q_value"].to_numpy()


def noise(seed: int, n_no: int = 40, n_d: int = 40) -> pd.DataFrame:
    """A random frame: NO, CV and SW pairs with Delta ~ N(0, 1) over bins 0-4."""
    r = np.random.default_rng(seed)
    parts = [
        block(cls, r.normal(size=k), int(b))
        for cls, n in (("NO", n_no), ("CV", n_d), ("SW", n_d))
        for b, k in zip(*np.unique(r.integers(0, 5, size=n), return_counts=True), strict=True)
    ]
    return stack(*parts)


# ------------------------------------------------------------------ examples 1-3, 8


def test_example1_cv_decoys_give_the_m6_01_example4_fdr() -> None:
    pairs = example_frame(sandwich=False)
    res = missed_operons(pairs, cfg())
    assert isinstance(res, MissedOperons)
    assert res.n_targets == 300
    assert res.n_decoys["convergent"] == 200
    assert q_at(res, pairs, 0.1) == pytest.approx(EX1, rel=1e-9)
    assert set(by_id(res)["which"]) == {"convergent"}


def test_example2_sandwich_decoys_win_where_their_fdr_is_larger() -> None:
    pairs = example_frame()
    res = missed_operons(pairs, cfg())
    assert res.n_decoys == {"convergent": 200, "sandwich": 70}
    assert q_at(res, pairs, 0.1) == pytest.approx(EX2, rel=1e-9)
    top = pairs.loc[(pairs["class"] == "NO") & np.isclose(pairs["delta_log2"].abs(), 0.1)]
    assert set(by_id(res).loc[top["pair_id"], "which"]) == {"sandwich"}


def test_the_rows_carry_the_spec_step4_test_columns() -> None:
    pairs = example_frame()
    res = missed_operons(pairs, cfg())
    assert set(ROW_COLS) <= set(res.rows.columns)
    s = scored(res)
    assert len(s) == 300
    assert set(s["pair_id"]) == set(pairs.loc[pairs["class"] == "NO", "pair_id"])
    d = pairs.set_index("pair_id").loc[s["pair_id"], "delta_log2"].to_numpy()
    np.testing.assert_allclose(s["stat"].to_numpy(dtype=float), -np.abs(d), rtol=0, atol=1e-12)
    np.testing.assert_allclose(s["effect"].to_numpy(dtype=float), np.abs(d), rtol=0, atol=1e-12)
    # lead decision 1: a decoy FDR has no p; no CI for a single pair
    for col in ("p_value", "ci_low", "ci_high", "n_perm"):
        assert s[col].isna().all(), col
    assert set(s["method"]) == {"decoy_fdr"}
    assert set(s["null_model"]) == {"convergent+sandwich_reweighted"}
    q = s["q_value"].to_numpy(dtype=float)
    assert np.all((q >= 0.0) & (q <= 1.0))


def test_example3_the_most_co_abundant_target_has_the_smallest_q() -> None:
    r = np.random.default_rng(3)
    pairs = stack(
        block("NO", r.uniform(0.5, 2.0, 30), 1), block("NO", [0.02], 1),
        block("CV", r.normal(size=40), 1), block("SW", r.normal(size=40), 1),
    )  # fmt: skip
    s = by_id(missed_operons(pairs, cfg()))
    best = pairs.loc[np.isclose(pairs["delta_log2"], 0.02), "pair_id"].item()
    assert s.loc[best, "q_value"] == s["q_value"].min()
    assert s.loc[best, "stat"] == s["stat"].max()


def test_example8_one_decoy_set_from_config() -> None:
    pairs = example_frame()
    res = missed_operons(pairs, cfg(**CV_ONLY))
    assert set(by_id(res)["which"]) == {"convergent"}
    assert q_at(res, pairs, 0.1) == pytest.approx(EX1, rel=1e-9)
    assert "sandwich" not in res.n_decoys or res.n_decoys["sandwich"] == 0


def test_pi0_from_config_is_honoured() -> None:
    pairs = example_frame(sandwich=False)
    res = missed_operons(pairs, cfg(pi0=0.5))
    assert q_at(res, pairs, 0.1) == pytest.approx(0.5 * EX1, rel=1e-9)


def test_an_empty_decoy_list_is_a_config_error() -> None:
    with pytest.raises(ValueError):
        missed_operons(example_frame(), cfg(decoys=[]))


# ---------------------------------------------------------------------- targets, decoys


def test_only_both_detected_finite_no_pairs_are_targets() -> None:
    base = levels(TARGETS, "NO") + levels(CV, "CV") + levels(SW, "SW")
    extra = stack(
        *base,
        block("NO", [0.0], 1, detected=(True, False)),
        block("NO", [0.0], 1, detected=(False, False)),
        block("OP", [0.0] * 5, 1), block("DV", [0.0] * 5, 1),
        block("CV", [0.0], 1, detected=(False, True)),
    )  # fmt: skip
    nan_ab = block("NO", [0.0, 0.0], 1)
    nan_ab.loc[0, "abundance_log2_a"] = np.nan
    nan_ab.loc[1, "abundance_log2_a"] = -np.inf  # log2(0): lead decision 5
    nan_ab["delta_log2"] = nan_ab["abundance_log2_b"] - nan_ab["abundance_log2_a"]
    pairs = stack(extra, nan_ab)
    res = missed_operons(pairs, cfg())
    assert res.n_targets == 300
    assert res.n_decoys == {"convergent": 200, "sandwich": 70}
    assert len(scored(res)) == 300
    assert set(res.rows["pair_id"]) <= set(pairs.loc[pairs["class"] == "NO", "pair_id"])
    assert q_at(res, pairs, 0.1) == pytest.approx(EX2, rel=1e-9)


def test_rna_between_and_imputed_abundances_stay_among_the_targets() -> None:
    pairs = example_frame()
    no = pairs.index[pairs["class"] == "NO"]
    for i in no[:10]:
        pairs.at[i, "flags"] = ["rna_between"]
    genes = genes_for(pairs, abundance_imputed=dict.fromkeys(pairs.loc[no[10:20], "gene_a"], True))
    res = missed_operons(pairs, cfg(), genes=genes)
    assert res.n_targets == 300
    assert set(pairs.loc[no[:20], "pair_id"]) <= set(scored(res)["pair_id"])


def test_partial_genes_are_never_eligible() -> None:
    pairs = example_frame()
    extra = stack(pairs, block("NO", [0.0] * 30, 0))
    partial = extra.loc[extra.index[-30:], "gene_b"]
    genes = genes_for(extra, partial=dict.fromkeys(partial, True))
    res = missed_operons(extra, cfg(), genes=genes)
    assert res.n_targets == 300
    assert not set(extra.loc[extra.index[-30:], "pair_id"]) & set(scored(res)["pair_id"])
    assert q_at(res, extra[extra.index < len(pairs)], 0.1) == pytest.approx(EX2, rel=1e-9)


def test_bin_counts_sum_to_the_scored_targets_and_decoys() -> None:
    res = missed_operons(example_frame(), cfg())
    bc = res.bin_counts
    assert len(bc) == len(cfg().operons.distance_bins_bp) - 1
    assert bc["N_NO"].tolist() == [200, 100, 0, 0, 0]
    assert bc["N_cv"].tolist() == [100, 100, 0, 0, 0]
    assert bc["N_sw"].tolist() == [50, 20, 0, 0, 0]
    assert bc["N_NO"].sum() == res.n_targets
    assert bc["N_cv"].sum() + bc["N_sw"].sum() == sum(res.n_decoys.values())


# ------------------------------------------------------------------- degenerate cases


@pytest.mark.parametrize("n", [0, 1, MIN_NO_PAIRS - 1])
def test_example4_too_few_no_pairs_gives_no_rows(n: int) -> None:
    assert MIN_NO_PAIRS == 20
    pairs = stack(block("NO", np.linspace(0, 1, n), 1), *levels(CV, "CV"), *levels(SW, "SW"))
    res = missed_operons(pairs, cfg())
    assert scored(res).empty
    if n:
        assert flagged(res, "too_few_no_pairs")
    assert len(res.bin_counts) == 5  # still written
    assert res.bin_counts["N_NO"].sum() == n


def test_undetected_pairs_do_not_count_towards_the_minimum() -> None:
    pairs = stack(
        block("NO", np.linspace(0, 1, 19), 1), block("NO", [0.1] * 6, 1, detected=(True, False)),
        *levels(CV, "CV"),
    )  # fmt: skip
    res = missed_operons(pairs, cfg())
    assert scored(res).empty
    assert flagged(res, "too_few_no_pairs")


def test_exactly_the_minimum_is_scored() -> None:
    pairs = stack(block("NO", np.linspace(0, 1, MIN_NO_PAIRS), 1), *levels(CV, "CV"))
    assert len(scored(missed_operons(pairs, cfg()))) == MIN_NO_PAIRS


def test_one_empty_decoy_set_is_dropped_with_a_log_line(caplog: pytest.LogCaptureFixture) -> None:
    pairs = example_frame(sandwich=False)
    with caplog.at_level(logging.INFO):
        res = missed_operons(pairs, cfg())  # both sets configured, no SW pair in the frame
    assert set(by_id(res)["which"]) == {"convergent"}
    assert q_at(res, pairs, 0.1) == pytest.approx(EX1, rel=1e-9)
    assert any("sandwich" in r.getMessage().lower() for r in caplog.records)


def test_both_decoy_sets_empty_gives_no_rows() -> None:
    res = missed_operons(stack(*levels(TARGETS, "NO")), cfg())
    assert scored(res).empty
    assert flagged(res, "too_few_no_pairs")


def test_example7_targets_in_a_bin_without_decoys(caplog: pytest.LogCaptureFixture) -> None:
    r = np.random.default_rng(7)
    pairs = stack(
        block("NO", r.normal(size=50), 4), block("CV", r.normal(size=50), 0),
        block("SW", r.normal(size=50), 0),
    )  # fmt: skip
    with caplog.at_level(logging.INFO):
        res = missed_operons(pairs, cfg())
    np.testing.assert_array_equal(scored(res)["q_value"].to_numpy(dtype=float), 1.0)
    bc = res.bin_counts
    assert bc["N_NO"].tolist() == [0, 0, 0, 0, 50]
    assert bc["N_cv"].tolist() == [50, 0, 0, 0, 0]
    assert bc["N_sw"].tolist() == [50, 0, 0, 0, 0]
    msgs = [r.getMessage().lower() for r in caplog.records if r.levelno >= logging.WARNING]
    assert any("bin" in m and "decoy" in m for m in msgs), caplog.text


def test_all_delta_identical_gives_one_tied_q() -> None:
    pairs = stack(
        block("NO", [0.3] * 25, 1), block("CV", [0.3] * 25, 1), block("SW", [-0.3] * 9, 2)
    )
    q = scored(missed_operons(pairs, cfg()))["q_value"].to_numpy(dtype=float)
    assert q.size == 25
    assert np.unique(q).size == 1
    assert 0.0 <= q[0] <= 1.0


def test_fewer_than_five_decoys_leaves_lfdr_na_and_q_unaffected() -> None:
    pairs = stack(block("NO", np.linspace(0.01, 1, 30), 1), block("CV", [0.05, 0.4, 0.9], 1))
    s = scored(missed_operons(pairs, cfg()))
    assert s["lfdr"].isna().all()
    q = s["q_value"].to_numpy(dtype=float)
    assert np.all(np.isfinite(q) & (q >= 0.0) & (q <= 1.0))


def test_lfdr_is_a_probability_when_defined() -> None:
    lf = scored(missed_operons(noise(11, 60, 60), cfg()))["lfdr"].to_numpy(dtype=float)
    assert np.isfinite(lf).all()
    assert np.all((lf >= 0.0) & (lf <= 1.0))


def test_relative_scale_gives_no_rows() -> None:
    pairs = example_frame()
    tags = [*pairs["gene_a"], *pairs["gene_b"]]
    genes = genes_for(
        pairs, abundance_type=dict.fromkeys(tags, "maxlfq"),
        abundance_scale=dict.fromkeys(tags, "relative"),
    )  # fmt: skip
    res = missed_operons(pairs, cfg(), genes=genes)
    assert scored(res).empty
    assert flagged(res, "relative_scale")


def test_the_m6_04_pair_frame_is_accepted() -> None:
    """The real ``pair_frame`` output (5 genes, one NO pair at most) runs and is too small."""
    genes = genes_of(BASE)
    pairs = frame(genes, ring(), {("G1", "G2"): 0.1, ("G3", "G4"): 0.9, ("G5", "G1"): 0.1})
    res = missed_operons(pairs, cfg())
    assert scored(res).empty
    assert res.n_targets < MIN_NO_PAIRS


# -------------------------------------------------------------- determinism, properties


def test_example9_same_frame_twice_gives_identical_rows() -> None:
    pairs = noise(9, 80, 80)
    a, b = missed_operons(pairs, cfg()), missed_operons(pairs.copy(), cfg())
    pd.testing.assert_frame_equal(a.rows, b.rows)
    pd.testing.assert_frame_equal(a.bin_counts, b.bin_counts)


@FAST
@given(st.integers(0, 2**31 - 1), st.integers(20, 80), st.integers(5, 60))
def test_property_q_is_non_increasing_in_stat_and_in_unit_interval(
    seed: int, n_no: int, n_d: int
) -> None:
    s = scored(missed_operons(noise(seed, n_no, n_d), cfg()))
    stat, q = s["stat"].to_numpy(dtype=float), s["q_value"].to_numpy(dtype=float)
    assert np.all((q >= 0.0) & (q <= 1.0))
    order = np.argsort(stat, kind="stable")
    assert np.all(np.diff(q[order]) <= 1e-12)


@FAST
@given(st.integers(0, 2**31 - 1))
def test_property_dropping_a_decoy_set_only_lowers_or_keeps_q(seed: int) -> None:
    pairs = noise(seed, 50, 30)
    both = by_id(missed_operons(pairs, cfg()))["q_value"]
    for one in (["convergent"], ["sandwich"]):
        single = by_id(missed_operons(pairs, cfg(decoys=one)))["q_value"]
        assert (single <= both + 1e-12).all(), one


@FAST
@given(st.integers(0, 2**31 - 1))
def test_property_invariant_to_row_order_and_pair_id_relabelling(seed: int) -> None:
    pairs = noise(seed, 40, 30)
    perm = np.random.default_rng(seed + 1).permutation(len(pairs))
    shuffled = pairs.iloc[perm].reset_index(drop=True)
    new_ids = {p: f"X{k}" for k, p in enumerate(np.random.default_rng(seed + 2).permutation(
        pairs["pair_id"].to_numpy()))}  # fmt: skip
    shuffled["pair_id"] = shuffled["pair_id"].map(new_ids)
    a = by_id(missed_operons(pairs, cfg()))
    b = by_id(missed_operons(shuffled, cfg()))
    b.index = b.index.map({v: k for k, v in new_ids.items()})
    b = b.sort_index()
    for col in ("stat", "q_value", "lfdr"):
        np.testing.assert_allclose(
            a[col].to_numpy(dtype=float), b[col].to_numpy(dtype=float), rtol=1e-9, equal_nan=True
        )
    assert (a["which"] == b["which"]).all()


# -------------------------------------------------- examples 5 and 6 (slow, fixed seeds)

SEEDS_EX5 = range(20)
#: one fixed seed (the first); a single frame's FDP is random - over seeds 0-9 it is 0.03-0.15
#: with mean 0.08 - so the FDR claim proper is M6-09's gate, not this test's
SEED_EX6 = 0


@pytest.mark.slow
def test_example5_shared_null_almost_never_calls() -> None:
    """No target at q <= 0.1 in >= 90 % of seeds; every call is false, so FDP is 1 or 0."""
    fdp = []
    for seed in SEEDS_EX5:
        s = scored(missed_operons(null_frame(np.random.default_rng(seed)), cfg()))
        fdp.append(float((s["q_value"] <= 0.1).any()))
    assert 1.0 - np.mean(fdp) >= 0.9
    assert np.mean(fdp) <= 0.1


@pytest.mark.slow
def test_example6_planted_tiny_effects_are_called_with_fdp_control() -> None:
    """100 planted pairs at |Delta| <= 1e-3 SD (M6-01 decision 12: the only effect this family
    detects); power is asserted only here, never at 0.25-1 SD."""
    pairs = null_frame(np.random.default_rng(SEED_EX6), planted=100)
    planted = pairs.attrs["planted"]
    assert len(planted) == 100
    s = scored(missed_operons(pairs, cfg()))
    called = set(s.loc[s["q_value"] <= 0.1, "pair_id"])
    assert len(called & planted) >= 50
    assert len(called - planted) / max(len(called), 1) <= 0.1
    assert not math.isnan(s["q_value"].min())
